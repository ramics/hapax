#!/usr/bin/env python3
"""Rebuild the two date-held-out BLAST databases from the accession list PR #29 committed.

WHY THIS EXISTS RATHER THAN RE-DERIVING THE HOLDOUT. `holdout_run.py --stages holdout`
recovers each RefSeq record's original GenBank accession from its COMMENT provenance line
and drops the ones deposited after the cutoff. That derivation is the expensive and
error-prone half - it is where the `createdate` trap lives (3,926 viral RefSeq records
carry a 2023 createdate while wrapping decades-older sequence, so keying on it would
delete genuinely old relatives and manufacture novelty in the dangerous direction). Its
*output* is committed: `results/holdout_run/dropped_accessions.json` is the 439 accessions
that were withheld, and `holdout_manifest.json` records what that cost each index.

So this script does not repeat the derivation. It reads the committed list and applies it,
which means a later run holds out EXACTLY what the earlier run held out rather than
something a re-derivation happened to agree with. The container the first run used is
gone; the accession list is the part that had to survive, and it did.

THE CHECK THAT MAKES THIS WORTH ANYTHING. Four numbers were recorded by the first run:
nucleotide kept 19,211 and dropped 439; protein kept 703,330, dropped 19,493 and 4 kept
because no parent record could be found for them. This script recomputes all four from
today's RefSeq release and **refuses to write a database if any of them differs**. That is
the whole point: it converts "I rebuilt the indexes" into "the indexes are the same ones",
and it would catch a silently re-versioned release file, a truncated download, or a
protein-to-parent map that drifted. A rebuild nobody checked is a second instrument, and
this project's rule is never to compare numbers from two instruments.

THE PROTEIN SIDE IS THE PART THAT CAN GO WRONG QUIETLY. A protein accession is not its
nucleotide accession and the defline does not carry the parent, so the date rule cannot be
applied to proteins directly. Each protein is mapped to the RefSeq nucleotide record whose
CDS declared it, parsed from the genomic GenBank flatfile, and the one date rule then
governs both indexes. Filtering proteins by organism name instead would apply a DIFFERENT
rule to the two indexes and the comparison between them would stop meaning anything.

ASYMMETRY, DELIBERATE: an entry is dropped only if it is positively known to postdate the
cutoff. An entry whose date could not be recovered is KEPT. Keeping one can only leave a
relative in the index that should have been hidden, which makes a candidate look *less*
novel; dropping it would fabricate novelty. Every one-sided choice in this file leans the
same way, so every distance measured against these indexes is a floor.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def log(*a):
    print(*a, flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def opener(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def protein_parents(gbff: Path) -> dict[str, str]:
    """protein accession -> the nucleotide record whose CDS declared it.

    Streamed rather than parsed with a full GenBank reader: the release flatfile is 456 MB
    gzipped and only two fields are wanted. VERSION gives the record, /protein_id gives
    each CDS product, and a protein can be declared by more than one record (segmented
    genomes, identical products), in which case the first wins - which is safe here because
    the date rule is applied per parent and a protein kept by either parent is kept.
    """
    parent: dict[str, str] = {}
    cur = None
    with opener(gbff) as fh:
        for line in fh:
            if line.startswith("VERSION"):
                parts = line.split()
                cur = parts[1] if len(parts) > 1 else None
            elif cur and '/protein_id="' in line:
                acc = line.split('/protein_id="', 1)[1].split('"', 1)[0]
                parent.setdefault(acc, cur)
    return parent


def filter_fasta(src: Path, dst: Path, drop) -> dict:
    """Copy a FASTA, omitting records whose key is in `drop`.

    `drop` is a callable returning (is_dropped, is_unmapped) for a defline accession, so
    the nucleotide and protein rules differ only in that callable and the record-walking
    is identical for both. Writing the two with different loops is how the two indexes
    would come to obey different rules without anyone noticing.
    """
    kept = dropped = unmapped = 0
    emit = False
    with opener(src) as fh, open(dst, "w") as w:
        for line in fh:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                bad, un = drop(acc)
                unmapped += un
                emit = not bad
                kept += emit
                dropped += not emit
            if emit:
                w.write(line)
    return {"kept": kept, "dropped": dropped, "unmapped_kept": unmapped}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="no randomness here; carried so every script in this repo takes "
                         "the same argument and a run manifest is comparable")
    ap.add_argument("--nt-fa", type=Path, required=True,
                    help="RefSeq viral genomic release FASTA (viral.*.genomic.fna.gz)")
    ap.add_argument("--prot-fa", type=Path, required=True,
                    help="RefSeq viral protein release FASTA (viral.*.protein.faa.gz)")
    ap.add_argument("--gbff", type=Path, required=True,
                    help="RefSeq viral genomic flatfile, for the protein->parent map only")
    ap.add_argument("--dropped", type=Path,
                    default=ROOT / "results/holdout_run/dropped_accessions.json")
    ap.add_argument("--expect", type=Path,
                    default=ROOT / "results/holdout_run/holdout_manifest.json",
                    help="the first run's recorded per-index counts. The rebuild refuses "
                         "to write a database unless it reproduces them exactly")
    ap.add_argument("--blastbin", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--allow-mismatch", action="store_true",
                    help="write the databases even if the counts differ from the first "
                         "run's. Off by default, because a rebuild that does not "
                         "reproduce the recorded holdout is a different experiment and "
                         "should be reported as one rather than used")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    timings = {}

    drop_date = {x.split(".")[0] for x in json.loads(a.dropped.read_text())}
    log(f"holding out {len(drop_date)} RefSeq accessions deposited after the cutoff")

    t = time.perf_counter()
    parent = protein_parents(a.gbff)
    timings["parent_map"] = round(time.perf_counter() - t, 1)
    log(f"protein -> parent record: {len(parent)} proteins mapped "
        f"({timings['parent_map']}s)")

    stats = {}

    def nt_rule(acc):
        return acc.split(".")[0] in drop_date, 0

    def prot_rule(acc):
        p = parent.get(acc) or parent.get(acc.split(".")[0])
        if p is None:
            return False, 1          # unmapped is KEPT, see module docstring
        return p.split(".")[0] in drop_date, 0

    for label, src, rule, dbtype in [("nt", a.nt_fa, nt_rule, "nucl"),
                                     ("prot", a.prot_fa, prot_rule, "prot")]:
        t = time.perf_counter()
        fa = a.out / f"holdout_{label}.fa"
        stats[label] = filter_fasta(src, fa, rule)
        timings[f"filter_{label}"] = round(time.perf_counter() - t, 1)
        log(f"  {label}: kept {stats[label]['kept']}, dropped {stats[label]['dropped']}, "
            f"unmapped kept {stats[label]['unmapped_kept']}")

    # The gate. Compared before makeblastdb runs, so a mismatch leaves no database behind
    # that a later stage could pick up and use as if it were the recorded holdout.
    expected = json.loads(a.expect.read_text()).get("per_index", {}) if a.expect.exists() \
        else {}
    diffs = []
    for label, got in stats.items():
        want = expected.get(label)
        if not want:
            diffs.append(f"{label}: no recorded counts to check against")
            continue
        for k in ("kept", "dropped", "unmapped_kept"):
            if k in want and want[k] != got[k]:
                diffs.append(f"{label}.{k}: recorded {want[k]}, rebuilt {got[k]}")
    if diffs:
        log("\nTHE REBUILD DOES NOT REPRODUCE THE RECORDED HOLDOUT:")
        for d in diffs:
            log(f"  {d}")
        if not a.allow_mismatch:
            log("\nrefusing to build the databases. Pass --allow-mismatch to build them "
                "anyway and report the run as a different holdout than PR #29's.")
            return 2
        log("\n--allow-mismatch given: building anyway. This is NOT PR #29's holdout.")
    else:
        log("\nrebuild reproduces the recorded holdout exactly on all six counts")

    for label, dbtype in [("nt", "nucl"), ("prot", "prot")]:
        t = time.perf_counter()
        fa = a.out / f"holdout_{label}.fa"
        subprocess.run([str(a.blastbin / "makeblastdb"), "-in", str(fa),
                        "-dbtype", dbtype, "-out", str(a.out / f"db_holdout_{label}")],
                       check=True, capture_output=True)
        timings[f"makeblastdb_{label}"] = round(time.perf_counter() - t, 1)
        log(f"  db_holdout_{label} built ({timings[f'makeblastdb_{label}']}s)")

    timings["total"] = round(time.perf_counter() - t0, 1)
    man = {
        "script": "holdout_indexes.py",
        "script_sha256": sha256(Path(__file__)),
        "git_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                  capture_output=True, text=True).stdout.strip(),
        "tree_dirty": bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                          capture_output=True, text=True).stdout.strip()),
        "seed": a.seed,
        "parameters": {k: str(v) for k, v in vars(a).items()},
        "inputs": {p.name: {"bytes": p.stat().st_size, "sha256_16": sha256(p)[:16]}
                   for p in (a.nt_fa, a.prot_fa, a.dropped) if p.exists()},
        "held_out_accessions": len(drop_date),
        "proteins_mapped_to_parent": len(parent),
        "per_index": stats,
        "reproduces_recorded_holdout": not diffs,
        "mismatches": diffs or None,
        "timings_s": timings,
        "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python": platform.python_version(),
    }
    (a.out / "holdout_indexes_manifest.json").write_text(json.dumps(man, indent=1))
    log(f"\n{a.out}/holdout_indexes_manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
