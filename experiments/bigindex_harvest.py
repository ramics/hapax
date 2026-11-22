"""Summarise what a 319 GB nucleotide index does to the divergent contig panel.

WHY THIS EXISTS. `results/divergent_contigs/` ranks contigs by how far their proteins sit
from anything published, using the indexes that fit in this container: `viral.dmnd` and
`sprot.dmnd`, a few hundred megabytes between them. The obvious objection is that the
distance is a property of the index, not of the contig, and that a full-size index would
collapse the panel. So the same 284 contigs went to NCBI core_nt (319 GB, 2026-10-22) on
EC2, divergent arm and near-reference arm together, and this script reads the harvest.

The near-reference arm is the calibration: it is contigs that already match a reference
closely, so core_nt must recover them. An arm that does not come back at ~100% identity
over most of its length means the search, not the panel, is what moved.

TWO THRESHOLDS, both stated on every row this prints. A hit counts only at
E < 1e-5 (`SIGNIFICANT_EVALUE`), because identity on a 40 nt chance alignment is
meaningless and this repo has had to retract a number quoted without one. And identity is
read beside alignment fraction, never alone: 97% over 6% of a contig is a conserved
fragment, not an identification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SIGNIFICANT_EVALUE = 1e-5
# "Collapsed" means the big index turned the contig into an identification: strain-level
# identity over at least half the contig. "Still distant" means it did not get close on
# either axis. The band between the two is deliberately named rather than split, because
# nothing here makes one cut in it defensible.
COLLAPSE_PIDENT, COLLAPSE_FRAC = 90.0, 0.5
DISTANT_PIDENT = 80.0


def sha256_of(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def quant(xs: list[float]) -> dict:
    if not xs:
        return {"n": 0}
    s = sorted(xs)
    return {"n": len(s), "min": round(s[0], 3), "p25": round(s[len(s) // 4], 3),
            "median": round(statistics.median(s), 3),
            "p75": round(s[(3 * len(s)) // 4], 3), "max": round(s[-1], 3)}


def arm(rows: list[dict]) -> dict:
    sig = [r for r in rows if r.get("hit") and r["evalue"] < SIGNIFICANT_EVALUE]
    collapsed = [r for r in sig
                 if r["pident"] >= COLLAPSE_PIDENT and r["aln_frac"] >= COLLAPSE_FRAC]
    distant = [r for r in sig
               if r["pident"] < DISTANT_PIDENT or r["aln_frac"] < COLLAPSE_FRAC]
    # An explicit census, because an empty table and a search that never ran look the
    # same. Every contig queried is in exactly one of these counts.
    return {
        "contigs_queried": len(rows),
        "with_any_alignment": sum(1 for r in rows if r.get("hit")),
        "significant": len(sig),
        "no_significant_hit": len(rows) - len(sig),
        "identity_pct": quant([r["pident"] for r in sig]),
        "alignment_fraction_of_contig": quant([r["aln_frac"] for r in sig]),
        "alignment_length_nt": quant([float(r["aln_len"]) for r in sig]),
        "collapsed_to_identification": {
            "n": len(collapsed),
            "share_of_arm": round(len(collapsed) / len(rows), 4) if rows else None,
            "definition": f">= {COLLAPSE_PIDENT}% identity over >= "
                          f"{COLLAPSE_FRAC:.0%} of the contig",
        },
        "still_distant": {
            "n": len(distant),
            "share_of_arm": round(len(distant) / len(rows), 4) if rows else None,
            "definition": f"< {DISTANT_PIDENT}% identity, or covered over < "
                          f"{COLLAPSE_FRAC:.0%} of the contig",
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, required=True,
                    help="rows_*.jsonl harvested from the EC2 run")
    ap.add_argument("--timing", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=1, help="unused; the harvest is a read")
    a = ap.parse_args()
    t0 = time.time()

    rows = [json.loads(ln) for ln in a.rows.read_text().splitlines() if ln.strip()]
    timing = json.loads(a.timing.read_text())
    classes = sorted({r["class"] for r in rows})
    out = {
        "seed": a.seed,
        "rows_file": str(a.rows),
        "significant_evalue": SIGNIFICANT_EVALUE,
        "search": {k: timing.get(k) for k in
                   ("label", "program", "db", "db_version", "queries", "threads",
                    "blastn_task", "seconds", "exit_code")},
        "search_command": timing.get("command"),
        "arms": {c: arm([r for r in rows if r["class"] == c]) for c in classes},
        "all": arm(rows),
    }
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                             capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                    check=True, capture_output=True,
                                    text=True).stdout.strip())
    except Exception:
        sha, dirty = "unknown", True
    out["git_sha"] = sha + ("-dirty" if dirty else "")
    out["tree_dirty"] = dirty
    out["script_sha256"] = sha256_of(Path(__file__).resolve())
    out["timings_seconds"] = {"total": round(time.time() - t0, 3)}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"wrote {a.out}", file=sys.stderr)
    for c, d in out["arms"].items():
        print(f"{c}: {d['significant']}/{d['contigs_queried']} significant, "
              f"median {d['identity_pct'].get('median')}% over "
              f"{d['alignment_length_nt'].get('median')} nt "
              f"({d['alignment_fraction_of_contig'].get('median')} of contig), "
              f"collapsed {d['collapsed_to_identification']['n']}, "
              f"still distant {d['still_distant']['n']}", file=sys.stderr)


if __name__ == "__main__":
    main()
