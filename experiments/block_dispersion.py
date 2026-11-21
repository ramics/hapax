"""How much wider than binomial the spread between file blocks really is, and why 5 blocks cannot say.

WHY THIS EXISTS. Every cascade-derived number in this repo is five contiguous blocks of one
flow cell, not a random sample (`subtract_known.py` skips (seed-1) x 1.5M pairs), so a
binomial interval on one block understates the spread if the library has position structure.
The companion thread measured that on the cascade's classes and found the abundant ones
inflated ~1.4-1.6x and the residue ones near 1. This re-derives it from the whole-library
label table, for two reasons.

**First, to check it.** A second instrument on the same question is worth more than agreement
in prose.

**Second, because the obvious estimator gets it exactly backwards.** The between-block SD of
FIVE blocks is a variance estimate from five samples: its own sampling error is enormous, and
on this library it reads 0.44-0.65x binomial for the abundant classes - i.e. "no inflation,
if anything less" - where all 87 non-overlapping windows of the same width give 1.44x and
1.61x. Five numbers cannot measure a variance. The same five-sample spread is what this
repo's "min-max across seeds" reports, so those ranges are not wrong but they are narrow,
and the factor to widen them by has to come from the full set of windows.

The per-window series also gives the lag-1 autocorrelation, which is the direct evidence that
the structure is positional: adjacent windows correlate on the abundant classes and not on the
residue classes.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance  # noqa: E402

CLASSES = {"none": 0, "rrna": 1, "identified": 2, "unident": 3, "other_ctg": 4}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="recorded; the pass is deterministic")
    ap.add_argument("--classes", type=Path, default=ROOT / "data/work/library_label_classes.npz")
    ap.add_argument("--window", type=int, default=100000, help="the cascade's block size")
    ap.add_argument("--block-stride", type=int, default=1500000, help="the cascade's seed stride")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    d = np.load(args.classes)
    spot, c1, c2 = d["spot"], d["class1"], d["class2"]
    # Pair-level class, in the cascade's own precedence: a pair is ribosomal if either mate
    # is, then identified, then unidentifiable-contig, then outside-split, else nothing.
    pair = np.zeros(len(spot), dtype=np.int8)
    for code in (4, 3, 2, 1):
        pair[(c1 == code) | (c2 == code)] = code

    n = args.window
    out: dict = {"seed": args.seed, "window": n, "library_pairs": int(len(spot)), "classes": {}}
    for name, code in CLASSES.items():
        whole = float((pair == code).mean())
        binom = float(np.sqrt(whole * (1 - whole) / n))
        win = []
        for lo in range(0, int(spot.max()) - n, n):
            m = (spot > lo) & (spot <= lo + n)
            if m.sum() >= n * 0.9:
                win.append(float((pair[m] == code).mean()))
        w = np.array(win)
        wc = w - w.mean()
        five = np.array([float((pair[(spot > s * args.block_stride)
                                    & (spot <= s * args.block_stride + n)] == code).mean())
                         for s in range(5)])
        out["classes"][name] = {
            "whole_library": whole, "binomial_sd": binom, "windows": int(w.size),
            "between_window_sd": float(w.std(ddof=1)),
            "inflation_all_windows": float(w.std(ddof=1) / binom),
            "inflation_from_5_blocks_only": float(five.std(ddof=1) / binom),
            "lag1_autocorrelation": float((wc[:-1] * wc[1:]).sum() / (wc * wc).sum()),
            "five_block_mean": float(five.mean()),
            "five_block_minus_whole_pct_points": float((five.mean() - whole) * 100),
        }
    out["reading"] = ("Widen a binomial interval on rrna or identified by its "
                      "inflation_all_windows factor; the residue classes (none, unident) are "
                      "near 1 and barely move. Never take the factor from 5 blocks: that "
                      "column is a 5-sample variance and reads below 1 where the truth is 1.6.")
    out["manifest"] = {"script": "experiments/block_dispersion.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()}}
    (args.out / "block_dispersion.json").write_text(json.dumps(out, indent=1))
    print(f"{'class':11s} {'whole':>8s} {'nwin':>4s} {'inflation':>9s} {'from 5 blk':>10s} {'lag1':>6s} {'5blk-whole':>10s}")
    for name, v in out["classes"].items():
        print(f"{name:11s} {v['whole_library']*100:7.3f}% {v['windows']:4d} "
              f"{v['inflation_all_windows']:9.2f} {v['inflation_from_5_blocks_only']:10.2f} "
              f"{v['lag1_autocorrelation']:+6.2f} {v['five_block_minus_whole_pct_points']:+10.3f}")


if __name__ == "__main__":
    main()
