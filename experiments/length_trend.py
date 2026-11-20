"""Does divergence affect contig length, as a continuous outcome?

The binary readout - did any contig cross 1,500 nt - collapses the run to one
bit, and at 15 seeds an arm it can only detect a difference of about 0.44 in the
crossing rate. The quantities that were thrown away are on disk: longest contig
and mass over 1,500 nt per run. Regressing log length on divergence uses them,
and the single-lineage controls are the null anchor, since divergence between two
lineages cannot affect a run that contains only one of them. A slope there is
evidence of a confound - different genomes at different rungs - rather than of an
effect.

Reported with an OLS interval and an exact-style permutation test on the slope,
because n is small and the residuals are not obviously normal.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from statistics import mean


def load(paths: list[Path]) -> list[dict]:
    out = []
    for p in paths:
        with p.open() as fh:
            out += [json.loads(line) for line in fh if line.strip()]
    return out


def ols(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    n = len(xs)
    mx, my = mean(xs), mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    resid = [y - (a + b * x) for x, y in zip(xs, ys)]
    se = math.sqrt((sum(r * r for r in resid) / (n - 2)) / sxx)
    return a, b, se


def perm_p(xs: list[float], ys: list[float], b_obs: float, draws: int, seed: int) -> float:
    rng = random.Random(seed)
    ys = list(ys)
    hits = 0
    for _ in range(draws):
        rng.shuffle(ys)
        _, b, _ = ols(xs, ys)
        hits += abs(b) >= abs(b_obs)
    return (hits + 1) / (draws + 1)


def arm(tag: str, data: list[dict], key: str, log1p: bool, draws: int, seed: int) -> None:
    xs = [r["divergence"] for r in data]
    ys = [math.log1p(r[key]) if log1p else math.log(max(r[key], 1)) for r in data]
    a, b, se = ols(xs, ys)
    t = 2.09 if len(xs) < 25 else 2.01
    print(f"{tag:32} n={len(xs):3}  slope per 1% divergence = {b / 100:+.4f} "
          f"[{(b - t * se) / 100:+.4f}, {(b + t * se) / 100:+.4f}]  "
          f"t={b / se:+.2f}  perm P={perm_p(xs, ys, b, draws, seed):.4f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, nargs="+", required=True)
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rows = load(args.rows)
    mix = [r for r in rows if r["mode"] == "both"]
    ctl = [r for r in rows if r["mode"] != "both"]
    print(f"{len(mix)} mixture runs, {len(ctl)} single-lineage control runs")
    print("\nlog(longest contig) on divergence:")
    arm("mixtures", mix, "longest", False, args.draws, args.seed)
    if ctl:
        arm("controls (null anchor)", ctl, "longest", False, args.draws, args.seed)
    print("\nlog1p(mass >= 1500 nt) on divergence:")
    arm("mixtures", mix, "mass_over_1500", True, args.draws, args.seed)


if __name__ == "__main__":
    main()
