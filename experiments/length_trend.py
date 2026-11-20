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


def rung_mean_fit(tag: str, data: list[dict], key: str, log1p: bool,
                  weighted: bool = False) -> None:
    """Fit on the rung means rather than on pooled runs.

    Note the algebra before reading anything into the weighted version:
    weighting each rung mean by its inverse variance gives weights proportional
    to n_i, which reproduces the pooled OLS slope exactly. It is the same
    estimator, so it cannot remove the unbalance it was meant to remove - it only
    confirms the arithmetic. Giving every rung equal weight does balance the
    design, at the cost of efficiency, and is the fit to quote when the concern is
    that one contrast with many runs dominates.
    """
    from collections import defaultdict
    g: dict[float, list[float]] = defaultdict(list)
    for r in data:
        g[r["divergence"]].append(math.log1p(r[key]) if log1p else math.log(max(r[key], 1)))
    rungs = sorted(g)
    within = [y - mean(g[d]) for d in rungs for y in g[d]]
    dof = len(within) - len(rungs)
    s2 = sum(v * v for v in within) / dof if dof > 0 else 0.0
    xs = rungs
    ys = [mean(g[d]) for d in rungs]
    ws = ([len(g[d]) / s2 if s2 else 1.0 for d in rungs] if weighted
          else [1.0 for _ in rungs])
    sw = sum(ws)
    mx = sum(w * x for w, x in zip(ws, xs)) / sw
    my = sum(w * y for w, y in zip(ws, ys)) / sw
    sxx = sum(w * (x - mx) ** 2 for w, x in zip(ws, xs))
    b = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys)) / sxx
    if weighted:
        se = math.sqrt(1.0 / sxx)
    else:  # residual-based se, since equal weights carry no variance information
        a = my - b * mx
        resid = [y - (a + b * x) for x, y in zip(xs, ys)]
        dofr = len(xs) - 2
        se = math.sqrt((sum(r * r for r in resid) / dofr) / sxx) if dofr > 0 else float("nan")
    t = 2.57 if len(xs) <= 7 else 2.0  # t(5) at 95% for seven points
    print(f"{tag:32} rungs={len(xs)}  slope per 1% divergence = {b / 100:+.4f} "
          f"[{(b - t * se) / 100:+.4f}, {(b + t * se) / 100:+.4f}]  t={b / se:+.2f}")


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
    rung_mean_fit("mixtures, rung means (equal)", mix, "longest", False)
    rung_mean_fit("mixtures, rung means (inv-var)", mix, "longest", False, weighted=True)
    arm("mixtures, pooled runs", mix, "longest", False, args.draws, args.seed)
    if ctl:
        arm("controls (null anchor)", ctl, "longest", False, args.draws, args.seed)
    print("\nlog1p(mass >= 1500 nt) on divergence:")
    rung_mean_fit("mixtures, rung means (equal)", mix, "mass_over_1500", True)
    arm("mixtures, pooled runs", mix, "mass_over_1500", True, args.draws, args.seed)


if __name__ == "__main__":
    main()
