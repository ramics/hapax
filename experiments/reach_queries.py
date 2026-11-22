"""Build the query panel for the reach comparison, and its two controls.

THE QUESTION THIS PANEL EXISTS TO ANSWER. For a contig that sits far from anything
published, which search method recovers the most distant real relative, and at what cost?
Every arm of `experiments/reach_arms.py` is scored on exactly these sequences, in one pass,
so a difference between arms is a difference between instruments rather than between query
sets. That is the whole reason the panel is built once, here, and uploaded rather than
re-derived on the instance.

WHAT GOES IN, AND WHY EACH PILE IS NEEDED.

1. **The test arm**: the seven near-complete divergent genomes of
   `results/divergent_contigs/REPORT.md` §3 — 3.0 to 4.5 kb at 116-467x mapped depth, even
   coverage end to end, polymerase identity 44.9-84.8%. These are the material the panel is
   about: complete *Fiersviridae* gene sets whose nearest named relative is at family-level
   distance.

2. **The ground-truth arm**: the withheld-genome contig of `results/holdout_cascade/` —
   3,990 nt, assembled from a library carrying the genome at ~452x, whose true nucleotide
   divergence from its nearest *available* relative is 30.28% because the genome itself was
   deposited after the index's cutoff and withheld. **It is the only query here whose
   distance nobody chose**, so it is what calibrates "most distant relative recovered"
   against a known answer. Its already-referenced companion from the same library
   (`leak_already_referenced`, 6,383 nt at 1,188x) ships with it as the zero-distance end of
   the same instrument.

3. **The negative arm**: seven near-reference contigs from the same assembly, same family
   (two or more of `RNA_replicase_B`, `Levi_coat`, `Phage_mat-A`), same length band, at
   95.2-99.2% identity to their subject. **This is the control that decides every row.** A
   method that reaches further on the test arm and no further here has extended reach; a
   method that reaches further on both is reading its own threshold, and this project has
   been caught by exactly that shape three times in one day (CLAUDE.md, "the tell in all
   three was the control moving with the test"). Matching the family matters as much as
   matching the length: a plant virus at 98% identity has a different architecture and would
   make the control unreadable, so the plant-virus controls in the same length band are
   deliberately not used.

4. **The shuffled null**: every record above, dinucleotide-shuffled, at five seeds. The
   convention here is five seeds on any headline number, and this panel is the one place it
   bites: the real contigs are fixed sequence with nothing to re-draw, so the only
   stochastic element in the whole comparison is the shuffle. Five of them is what turns
   "the null read zero" from one draw into a range. The implementation is
   imported from `experiments/holdout_run.py` — the fixed Altschul-Erikson walk that
   `results/fail_loudly_audit/` measures at a 0.00% return-the-input rate — and NOT
   reimplemented. Two hand-rolled walks in this repo answered a dead end with `return seq`,
   put real sequence into the null pile, and manufactured two kills in two days. The rate is
   a property of the input as much as of the code, so this script measures the fallback rate
   on *these* sequences and refuses to write a panel if anything falls through.

WHAT THE NULL CAN AND CANNOT DO HERE. A dinucleotide shuffle preserves stop-codon density,
so the shuffled records yield few and short ORFs. That is the honest behaviour and it is
also a limit: the per-ORF arms have a nearly empty null, so the null calibrates *existence*
of a significant alignment and not any field conditional on one. `results/divergent_known/`
§2 reached the same wall from the other side and said so; this script records the ORF counts
so the report can state it rather than imply it.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The seven of `results/divergent_contigs/REPORT.md` §3, by name rather than by filter.
#: A filter would be re-derivable but it would also drift: the panel those seven numbers
#: were published against is this list, and a re-run that silently selected a different
#: seven would not be a re-run. Depth and polymerase identity are carried here only so the
#: manifest can be read without opening the workup file.
DIVERGENT = [
    ("WW_211013_20-22hmix_S5_L001_NODE_625_length_3484_cov_127.079032",
     "Kahshuvirus borborenecus", 467.2, 44.9),
    ("WW_211013_20-22hmix_S5_L001_NODE_645_length_3452_cov_103.860171",
     "Fiersviridae", 402.4, 48.6),
    ("WW_211013_20-22hmix_S5_L001_NODE_829_length_3058_cov_22.692974",
     "Koteshevirus lutadaptatum", 296.8, 46.6),
    ("WW_211013_20-22hmix_S5_L001_NODE_416_length_3725_cov_82.670027",
     "Sholavirus sp. 'caenicola'", 281.7, 83.8),
    ("WW_211013_20-22hmix_S5_L001_NODE_246_length_3930_cov_37.719226",
     "Rodtovirus pelovivens", 170.7, 45.6),
    ("WW_211013_20-22hmix_S5_L001_NODE_223_length_3996_cov_43.609997",
     "Dugnivirus solivicinum", 142.4, 47.8),
    ("WW_211013_20-22hmix_S5_L001_NODE_112_length_4526_cov_26.779468",
     "Chimpavirus luticola", 116.2, 84.8),
]

#: Pfam domains that say a contig is the same architecture as the test arm. Two of three is
#: the bar: `Levi_coat` is the domain most often missed on a real genome because the coat
#: protein is short, so requiring all three would select for annotation luck.
FIERSVIRIDAE_DOMAINS = {"RNA_replicase_B", "Levi_coat", "Phage_mat-A"}

#: Minimum whole-contig protein identity for the negative arm. 95% is inside the band where
#: `results/divergent_contigs/` §2 measures a contig's own genes agreeing to within 2.35
#: points, which is what makes "the control did not move" a statement about the method
#: rather than about gene choice.
CONTROL_MIN_PIDENT = 95.0


def read_fasta(path: Path) -> dict[str, str]:
    """Name to sequence, name being the first whitespace-delimited token of the header."""
    out: dict[str, str] = {}
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def load_shuffle():
    """The fixed dinucleotide walk, imported from the module that owns it.

    Imported by file path rather than by package, because `experiments/` is not a package
    and the alternative - copying forty lines of Eulerian walk into this file - is the exact
    mistake `results/fail_loudly_audit/` exists to document. Four implementations of this
    function in one tree is how two of them came to be wrong.
    """
    spec = importlib.util.spec_from_file_location(
        "holdout_run_shuffle", ROOT / "experiments" / "holdout_run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.dinuc_shuffle


def pick_controls(workup: dict, lengths: list[int], want: int) -> list[dict]:
    """Near-reference contigs of the same family, matched to the test arm's lengths.

    Greedy nearest-length assignment rather than a top-n by length: the test arm spans
    3,058 to 4,526 nt and a top-n would pile the control at one end of that, which would
    confound length with arm in a comparison whose whole point is that length drives
    alignment statistics. Deterministic, and it takes each control at most once.
    """
    pool = []
    for row in workup["rows"]:
        if row["arm"] != "control_near_reference":
            continue
        if not row.get("polymerase_pident"):
            continue
        if (row.get("whole_contig_pident") or 0) < CONTROL_MIN_PIDENT:
            continue
        domains = {p["id"] for orf in row.get("orf_rows", [])
                   for p in (orf.get("pfam") or [])}
        if len(domains & FIERSVIRIDAE_DOMAINS) < 2:
            continue
        pool.append(row)
    chosen: list[dict] = []
    for target in lengths:
        remaining = [r for r in pool if r not in chosen]
        if not remaining:
            break
        remaining.sort(key=lambda r: (abs(r["contig_nt"] - target), r["contig"]))
        chosen.append(remaining[0])
    if len(chosen) < want:
        raise SystemExit(
            f"only {len(chosen)} near-reference controls of the same family in the length "
            f"band; the negative arm would be smaller than the test arm and the table "
            f"would not be readable")
    return chosen


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1,
                    help="the panel's own seed, recorded in the manifest. The real arms "
                         "are fixed sequence, so this seeds the shuffles and nothing else")
    ap.add_argument("--null-seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5],
                    help="one null pile per seed. Five, because a null that reads zero "
                         "once is a draw and a null that reads zero five times is a range")
    ap.add_argument("--identified", type=Path, default=ROOT / "data/reach/identified.fa",
                    help="the assembly's identified contigs; the test and control arms are "
                         "both drawn from it. Staged from "
                         "s3://.../hapax-bigindex/shared/contig_sets/")
    ap.add_argument("--holdout", type=Path,
                    default=Path("/mnt/project-files/holdout_cascade/holdout_contigs.fa"))
    ap.add_argument("--workup", type=Path,
                    default=ROOT / "results/divergent_contigs/workup_s1.json")
    ap.add_argument("--out", type=Path, default=Path("/tmp/claude-0/wk/reach_queries"))
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    dinuc_shuffle = load_shuffle()
    identified = read_fasta(args.identified)
    holdout = read_fasta(args.holdout)
    workup = json.loads(args.workup.read_text())

    records: list[dict] = []
    for i, (name, taxon, depth, pol) in enumerate(DIVERGENT, start=1):
        if name not in identified:
            raise SystemExit(f"{name} is not in {args.identified}; the test arm would be "
                             f"short and the table would silently describe six genomes")
        records.append(dict(qid=f"div{i}", arm="divergent", contig=name,
                            taxon=taxon, mapped_depth=depth, polymerase_pident=pol,
                            seq=identified[name]))

    controls = pick_controls(workup, [len(r["seq"]) for r in records], len(records))
    for i, row in enumerate(controls, start=1):
        name = row["contig"]
        if name not in identified:
            raise SystemExit(f"control {name} is not in {args.identified}")
        records.append(dict(qid=f"ctl{i}", arm="near_reference", contig=name,
                            taxon=row["taxon"],
                            mapped_depth=row["read_support"]["mapped_mean_depth"],
                            polymerase_pident=row["polymerase_pident"],
                            seq=identified[name]))

    for qid, header, arm, taxon in (
            ("hold", "candidate_withheld", "holdout", "withheld genome, 30.28% diverged"),
            ("leak", "leak_already_referenced", "holdout_reference",
             "already-referenced companion from the same library")):
        match = [k for k in holdout if k == header] or [k for k in holdout if header in k]
        if not match:
            raise SystemExit(f"{header} is not in {args.holdout}; the only query with "
                             f"ground truth would be missing")
        records.append(dict(qid=qid, arm=arm, contig=match[0], taxon=taxon,
                            mapped_depth=None, polymerase_pident=None,
                            seq=holdout[match[0]]))

    # The null, and the number that says whether it is one. `fallback` counts records the
    # walk could not shuffle; anything above zero means real sequence would be entering the
    # control pile, and this refuses rather than writing it. Measured per seed, because the
    # walk's dead end depends on the draw as well as on the sequence.
    shuffled: dict[int, list[dict]] = {}
    fallback: dict[int, list[str]] = {}
    for null_seed in args.null_seeds:
        seed_rng = random.Random((args.seed, null_seed).__hash__())
        pile, missed = [], []
        for rec in records:
            seq = dinuc_shuffle(rec["seq"].upper(), seed_rng)
            if seq is None:
                missed.append(rec["qid"])
                continue
            pile.append(dict(qid=f"sh{null_seed}_{rec['qid']}", arm="shuffled",
                             null_seed=null_seed, contig=rec["contig"],
                             taxon=f"dinucleotide shuffle of {rec['qid']}, seed {null_seed}",
                             mapped_depth=None, polymerase_pident=None, seq=seq))
        shuffled[null_seed], fallback[null_seed] = pile, missed
    if any(fallback.values()):
        raise SystemExit(
            f"the dinucleotide walk returned nothing for {dict(fallback)}: those records "
            f"would be absent from the null rather than shuffled into it. Raise `tries` or "
            f"exclude them explicitly - do not let the null be quietly smaller than the "
            f"test, which is how two kills were manufactured in two days")

    args.out.mkdir(parents=True, exist_ok=True)
    piles = [("panel", records)] + [(f"null_s{s}", shuffled[s]) for s in args.null_seeds]
    for stem, pile in piles:
        text = "".join(
            f">{r['qid']} {r['arm']} {r['contig']} | {r['taxon']}\n"
            + "\n".join(r["seq"][i:i + 70] for i in range(0, len(r["seq"]), 70)) + "\n"
            for r in pile)
        (args.out / f"{stem}.fa").write_text(text)

    try:
        sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"],
                                    capture_output=True, text=True,
                                    check=True).stdout.strip())
    except Exception:
        sha, dirty = "unknown", True

    manifest = dict(
        unit="one contig",
        question="which search method recovers the most distant real relative, at what cost",
        seed=args.seed,
        null_seeds=args.null_seeds,
        arms={arm: sum(1 for r in records if r["arm"] == arm)
              for arm in ("divergent", "near_reference", "holdout", "holdout_reference")}
        | {"shuffled": sum(len(p) for p in shuffled.values())},
        shuffle=dict(
            implementation="holdout_run.dinuc_shuffle (fixed Altschul-Erikson walk)",
            measured_at="results/fail_loudly_audit/ reports 0.00% fallback on read windows "
                        "and whole contigs; the rate is a property of the input, so it is "
                        "measured again here",
            fallback_records={s: len(v) for s, v in fallback.items()},
            fallback_rate=0.0,
            records_shuffled_per_seed=len(records),
            note="a dinucleotide shuffle preserves stop-codon density, so the null yields "
                 "few and short ORFs. It calibrates existence of an alignment, not any "
                 "field conditional on one",
        ),
        rows=[{k: v for k, v in r.items() if k != "seq"} | {"nt": len(r["seq"])}
              for r in records + [x for s in args.null_seeds for x in shuffled[s]]],
        git_sha=sha,
        tree_dirty=dirty,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (args.out / "panel.json").write_text(json.dumps(manifest, indent=1))

    for arm, n in manifest["arms"].items():
        print(f"  {arm:20s} {n}")
    print(f"shuffle fallback per seed, out of {len(records)} records: "
          f"{ {s: len(v) for s, v in fallback.items()} }")
    print(f"written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
