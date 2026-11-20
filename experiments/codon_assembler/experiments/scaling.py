"""Wall clock versus read count: is candidate generation near-linear?

The design bet is that indexing, not language, is the lever (CLAUDE.md). If
the minimizer index works, candidates per read stay bounded as depth grows
(they are capped by `--max-partners`), so total work should scale roughly
linearly with read count rather than quadratically. Coverage is raised on a
fixed world so that the per-read neighbourhood gets denser, which is the
worst case for an overlap assembler.

Usage: python experiments/scaling.py --seed 1 --out results/scaling
"""

from __future__ import annotations

import argparse
import csv
import shutil
import time
from pathlib import Path

from common import ensure_release_build, run_hapax, write_manifest

from python import swarm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--coverages", type=float, nargs="+", default=[25, 50, 100, 200, 400])
    ap.add_argument("--diversity", type=float, default=0.10)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    ensure_release_build()
    a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.perf_counter()
    for s in range(a.seed, a.seed + a.n_seeds):
        for cov in a.coverages:
            wdir = a.out / "work" / f"s{s}_c{cov:g}"
            w = swarm.build_world(s * 1000 + 7, a.diversity, coverage=cov)
            swarm.write_world(w, wdir, {"seed": s, "coverage": cov})
            for mode in ("codon", "nucleotide"):
                st = run_hapax(wdir / "reads.fa", wdir / mode, ["--mode", mode, "-t", str(a.threads)])
                rows.append({"seed": s, "coverage": cov, "mode": mode, "reads": st["reads"], "candidates": st["candidates"], "tier2": st["tier2"], "hmm_fraction": st["hmm_fraction"], **{f"t_{k}": v for k, v in st["seconds"].items()}, "wall": st["wall"]})
            shutil.rmtree(wdir)
            print(f"seed {s} cov {cov} done", flush=True)
    with open(a.out / "runs.tsv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    write_manifest(a.out, {"first": a.seed, "n": a.n_seeds}, {k: str(v) for k, v in vars(a).items()}, {"total": time.perf_counter() - t0})
    shutil.rmtree(a.out / "work", ignore_errors=True)


if __name__ == "__main__":
    main()
