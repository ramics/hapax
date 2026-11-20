"""Commit a subsample of the joint space itself, so the figure has provenance.

The repo's rule is that figures are generated from committed result files, never from a
working directory, because this project has had to withdraw claims whose prose drifted from
its outputs. The coverage and composition tables are 24 MB of intermediate and are
gitignored with the rest of `data/`, so this writes the points the figure plots - a
seeded subsample of each population, with every axis - into `results/`.
"""

from __future__ import annotations

import argparse
import json
import random
import hashlib
import subprocess
from pathlib import Path

from residue_structure_analyse import AXES, load

ROOT = Path(__file__).resolve().parents[1]


def _provenance() -> dict:
    """Git SHA, whether the tree was dirty, and a content hash of this script.

    The SHA alone is not provenance when the tree is dirty: it names HEAD, not the code that
    ran, so two runs can carry the same SHA and different code. That is not hypothetical -
    it is how a codon-table fix landing midway through a multi-seed run produced a results
    directory nobody could sort afterwards. The content hash pins the code whatever the tree
    looked like, and unlike a file timestamp it survives the container.
    """
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                         text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip())
    return {"git_sha": sha + ("-dirty" if dirty else ""),
            "tree_dirty": dirty,
            "script_sha256": hashlib.sha256(
                Path(__file__).read_bytes()).hexdigest()[:16]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--block", type=int, default=1, help="which seed block to draw")
    ap.add_argument("--work", type=Path, default=ROOT / "data/residue_structure")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_structure")
    ap.add_argument("--per-population", type=int, default=3000)
    ap.add_argument("--all-pairs", action="store_true",
                    help="write the coverage axes for EVERY pair of the block, not only "
                         "the residue ones, gzipped. The subtrahend-ladder check needs "
                         "residue defined at several rungs, and each rung is a different "
                         "subset of the same block")
    ap.add_argument("--join-table", action="store_true",
                    help="also write every residue pair of this block with its coverage "
                         "axes, keyed by read name, so another experiment can join on it")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    rows, dropped = load(args.work)
    names = set(open(ROOT / f"data/residue_structure/residue_names_s{args.block}.txt")
                .read().split())
    pops: dict[str, list[dict]] = {"residue": [], "subtracted": [], "null": []}
    for r in rows:
        if r["block"] == args.block:
            pops["residue" if r["name"] in names else "subtracted"].append(r)
        elif r["block"] == 9:
            pops["null"].append(r)
        elif r["block"] == 0:
            pops.setdefault(r["name"].split("|")[0], []).append(r)

    out = {"axes": AXES, "block": args.block, "pairs_dropped_no_kmers": dropped,
           "populations": {}}
    for k, v in pops.items():
        keep = v if len(v) <= args.per_population else rng.sample(v, args.per_population)
        out["populations"][k] = {
            "n_total": len(v),
            "n_plotted": len(keep),
            "points": [[round(r[a], 4) for a in AXES] for r in keep],
        }
    (args.out).mkdir(parents=True, exist_ok=True)
    (args.out / f"joint_points_s{args.block}.json").write_text(json.dumps(out))
    if args.join_table:
        # Per-read coverage axes for the whole residue population, not a subsample: the
        # cascade thread holds per-read contig depth from its mapping step, and the
        # prediction to test is that a high cov_step read has a mate on a contig of very
        # different depth - uneven coverage along a molecule rather than a junction.
        cols = ["log_cov", "cov_spread", "zero_frac", "mate_disagree", "cov_step"]
        with open(args.out / f"coverage_axes_s{args.block}.tsv", "w") as h:
            h.write("name\t" + "\t".join(cols) + "\n")
            for r in rows:
                if r["block"] == args.block and r["name"] in names:
                    h.write(r["name"] + "\t"
                            + "\t".join(f"{r[c]:.4f}" for c in cols) + "\n")
    if args.all_pairs:
        import gzip as _gz
        cols = ["log_cov", "cov_spread", "zero_frac", "mate_disagree", "cov_step"]
        with _gz.open(args.out / f"coverage_axes_allpairs_s{args.block}.tsv.gz", "wt") as h:
            h.write("name\tresidue\t" + "\t".join(cols) + "\n")
            for r in rows:
                if r["block"] == args.block:
                    h.write(f"{r['name']}\t{int(r['name'] in names)}\t"
                            + "\t".join(f"{r[c]:.4f}" for c in cols) + "\n")
    (args.out / f"points_manifest_s{args.block}.json").write_text(json.dumps({
        "script": "experiments/residue_structure_points.py",
        **_provenance(),
        "parameters": {k: str(v) for k, v in vars(args).items()},
    }, indent=2))
    print({k: (v["n_total"], v["n_plotted"]) for k, v in out["populations"].items()})


if __name__ == "__main__":
    main()
