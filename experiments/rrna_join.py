"""Where the ribosomal reads sit in the peptide-absence residue.

Two tables, measured independently and on the same reads:

- `results/rrna_content/rrna_s{seed}.tsv` (this branch) calls a read ribosomal by positive
  exact 31-mer containment against SILVA, which does not depend on a gene caller.
- `results/absence_residual/reads_s{seed}.tsv` (branch `claude/project-thread-clcuow`,
  commit 9e63679) carries, per read, whether every peptide 10-mer was absent from the
  protein index, whether the read maps to the sample's own assembly, that contig's
  coverage, and whether the read falls inside a Prodigal-predicted CDS.

They join row for row on read_id and mate: same read blocks, same mate-2 trim, names
untouched.

**The tension this is run to settle.** 34.4% of reads are ribosomal, and only 23.0% have
no peptide 10-mer at all. A ribosomal read should have no true peptide match, and a random
10-mer's chance of appearing in a 2.2e8-peptide index is around 2e-5 across a 1e13 space,
so chance matching cannot account for the gap. At least eleven points of the library is
therefore ribosomal AND matching some peptide in the index, which should not be possible
against a database of real proteins.

**Advance prediction, written before the join was run.** The index arm that produced 23.0%
includes proteins predicted by Prodigal from the sample's own assembly. Ribosomal RNA is
the most abundant thing in a total-RNA library, so it assembles into high-coverage contigs,
and a prokaryotic gene caller run over those contigs will call spurious genes on ribosomal
sequence. Those spurious peptides then match the ribosomal reads and remove them from the
absent set. So:

1. Ribosomal reads will be far less often absent than the arithmetic of a clean index
   predicts; nothing like all of them will be absent.
2. The ribosomal reads that are NOT absent will be concentrated on mapped reads inside a
   predicted CDS, because that is where the spurious peptide came from.
3. in_cds among mapped ribosomal reads will be substantial rather than near zero, which
   would be a direct measurement of the Prodigal miscall rather than an inference from it.

If instead the non-absent ribosomal reads are mostly unmapped, or mostly outside any CDS,
the prediction is wrong and the gap is more likely to be false positives in the SILVA call,
which the held-out control bounds separately.

Both the flag and the raw hit count are carried through, because this project has
repeatedly found the gradient more informative than the threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import statistics as st
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path, key_cols=("read_id", "mate")) -> dict:
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["read_id"], row["mate"])] = row
    return out


def rate(rows, pred) -> float:
    sel = [r for r in rows if pred(r)]
    return len(sel) / len(rows) if rows else 0.0


def cond(rows, where, pred) -> tuple[float, int]:
    sel = [r for r in rows if where(r)]
    if not sel:
        return 0.0, 0
    return sum(1 for r in sel if pred(r)) / len(sel), len(sel)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0,
                    help="0 joins all five seeds; the inputs are fixed files, so the flag "
                         "selects which rather than randomising anything")
    ap.add_argument("--rrna", type=Path, default=ROOT / "results" / "rrna_content")
    ap.add_argument("--residual", type=Path,
                    default=ROOT / "results" / "absence_residual",
                    help="directory holding reads_s1..s5.tsv; these were produced on branch "
                         "claude/project-thread-clcuow at commit 9e63679 and have since "
                         "merged to main, so the default is now the in-repo path")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    seeds = [args.seed] if args.seed else [1, 2, 3, 4, 5]
    per_seed = []

    for s in seeds:
        mine = load(args.rrna / f"rrna_s{s}.tsv")
        theirs = load(args.residual / f"reads_s{s}.tsv")
        shared = set(mine) & set(theirs)
        rows = []
        for k in shared:
            a, b = mine[k], theirs[k]
            rows.append({
                "rrna": a["rrna"] == "1",
                "hits": int(a["kmer_hits"]),
                "absent": b["absent"] == "1",
                "mapped": b["mapped"] == "1",
                "in_cds": b["in_cds"] == "1",
                "cov": float(b["contig_cov"]) if b["contig_cov"] not in ("", "NA") else None,
            })

        d = {
            "seed": s,
            "rows_joined": len(rows),
            "rows_mine": len(mine),
            "rows_theirs": len(theirs),
            "rrna_fraction": rate(rows, lambda r: r["rrna"]),
            "absent_fraction": rate(rows, lambda r: r["absent"]),
        }

        # 1. Are ribosomal reads absent? The arithmetic says they cannot all be.
        d["absent_given_rrna"], d["n_rrna"] = cond(
            rows, lambda r: r["rrna"], lambda r: r["absent"])
        d["absent_given_not_rrna"], d["n_not_rrna"] = cond(
            rows, lambda r: not r["rrna"], lambda r: r["absent"])

        # 2. Where do the non-absent ribosomal reads sit?
        d["mapped_given_rrna_not_absent"], d["n_rrna_not_absent"] = cond(
            rows, lambda r: r["rrna"] and not r["absent"], lambda r: r["mapped"])
        d["in_cds_given_rrna_not_absent_mapped"], d["n_rrna_not_absent_mapped"] = cond(
            rows, lambda r: r["rrna"] and not r["absent"] and r["mapped"],
            lambda r: r["in_cds"])
        d["in_cds_given_rrna_absent_mapped"], d["n_rrna_absent_mapped"] = cond(
            rows, lambda r: r["rrna"] and r["absent"] and r["mapped"],
            lambda r: r["in_cds"])

        # 3. The cell the other branch asked for: rRNA rate in each absent x mapped bin.
        for name, where in (
            ("absent_mapped", lambda r: r["absent"] and r["mapped"]),
            ("absent_unmapped", lambda r: r["absent"] and not r["mapped"]),
            ("present_mapped", lambda r: not r["absent"] and r["mapped"]),
            ("present_unmapped", lambda r: not r["absent"] and not r["mapped"]),
        ):
            d[f"rrna_given_{name}"], d[f"n_{name}"] = cond(
                rows, where, lambda r: r["rrna"])

        # 4. The gradient rather than the threshold.
        for name, where in (
            ("absent_mapped", lambda r: r["absent"] and r["mapped"]),
            ("absent_unmapped", lambda r: r["absent"] and not r["mapped"]),
        ):
            sel = [r["hits"] for r in rows if where(r)]
            d[f"median_hits_{name}"] = st.median(sel) if sel else 0.0
            d[f"mean_hits_{name}"] = st.mean(sel) if sel else 0.0

        # 5. Coverage of the contigs ribosomal reads map to, to check "abundant".
        cov = [r["cov"] for r in rows if r["rrna"] and r["mapped"] and r["cov"]]
        d["median_cov_rrna_mapped"] = st.median(cov) if cov else None
        cov2 = [r["cov"] for r in rows if not r["rrna"] and r["mapped"] and r["cov"]]
        d["median_cov_not_rrna_mapped"] = st.median(cov2) if cov2 else None

        per_seed.append(d)
        print(f"seed {s}: joined {len(rows):,}  rRNA {d['rrna_fraction']:.4f}  "
              f"absent {d['absent_fraction']:.4f}  "
              f"P(absent|rRNA) {d['absent_given_rrna']:.4f}  "
              f"P(rRNA|absent,mapped) {d['rrna_given_absent_mapped']:.4f}  "
              f"P(rRNA|absent,unmapped) {d['rrna_given_absent_unmapped']:.4f}")

    keys = [k for k, v in per_seed[0].items() if isinstance(v, float)]
    summary = {k: {"mean": st.mean([d[k] for d in per_seed]),
                   "sd": st.stdev([d[k] for d in per_seed]) if len(per_seed) > 1 else 0.0,
                   "min": min(d[k] for d in per_seed),
                   "max": max(d[k] for d in per_seed)} for k in keys}

    out = {"seeds": seeds, "per_seed": per_seed, "across_seeds": summary,
           "residual_source": "claude/project-thread-clcuow commit 9e63679"}
    (args.out / "rrna_join.json").write_text(json.dumps(out, indent=2))
    (args.out / "rrna_join_manifest.json").write_text(json.dumps({
        "script": "experiments/rrna_join.py",
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }, indent=2))

    print("\nacross seeds, mean (sd):")
    for k in ("absent_given_rrna", "absent_given_not_rrna",
              "mapped_given_rrna_not_absent", "in_cds_given_rrna_not_absent_mapped",
              "in_cds_given_rrna_absent_mapped",
              "rrna_given_absent_mapped", "rrna_given_absent_unmapped",
              "rrna_given_present_mapped", "rrna_given_present_unmapped",
              "median_hits_absent_mapped", "median_hits_absent_unmapped",
              "median_cov_rrna_mapped", "median_cov_not_rrna_mapped"):
        if k in summary:
            print(f"  {k:<40} {summary[k]['mean']:.4f} ({summary[k]['sd']:.4f})")


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
