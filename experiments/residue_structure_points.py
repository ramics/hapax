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
import subprocess
from pathlib import Path

from residue_structure_analyse import AXES, load

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--block", type=int, default=1, help="which seed block to draw")
    ap.add_argument("--work", type=Path, default=ROOT / "data/residue_structure")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_structure")
    ap.add_argument("--per-population", type=int, default=3000)
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
    (args.out / f"points_manifest_s{args.block}.json").write_text(json.dumps({
        "script": "experiments/residue_structure_points.py",
        "git_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                  capture_output=True, text=True).stdout.strip(),
        "parameters": {k: str(v) for k, v in vars(args).items()},
    }, indent=2))
    print({k: (v["n_total"], v["n_plotted"]) for k, v in out["populations"].items()})


if __name__ == "__main__":
    main()
