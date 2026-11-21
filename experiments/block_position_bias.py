"""Is a contiguous block of a FASTQ a random sample of the library?

WHY THIS EXISTS. Every percentage this project reports comes from a CONTIGUOUS block of
read pairs: `subtract_known.py` skips `(seed - 1) * 1,500,000` pairs and then reads the
next `--pairs` in file order, so "5 seeds" means five regions of one flow cell rather
than five random draws. Everything downstream inherits that. The per-read report samples
at random, but only WITHIN the block; the web view samples at random, but only within the
report. The randomness is nested inside a non-random selection, and nothing in the repo
has ever checked the outer layer.

Two different things could be wrong and they need separating, because only the second is
fatal:

  LOCATION   the block sits somewhere unrepresentative, so a reported fraction is biased.
             Block 1 is the worst suspect: it is the start of the file, which for an
             Illumina run is the first tile, and tile 1 of lane 1 is where edge effects,
             cluster-density artefacts and the optical duplicates of a flow cell corner
             live.
  DISPERSION the blocks differ from each other by more than sampling noise, so the
             Wilson intervals this project prints - which assume the block is a simple
             random sample of the library - are too narrow, and a five-seed spread
             understates the real uncertainty.

HOW IT IS MEASURED WITHOUT RE-RUNNING ANYTHING. The residue thread labelled all 8,731,224
pairs of SRR23973193 against SILVA and all 82,947 of the sample's own contigs, in file
order (`library_label_classes.npz`). That is the content of every pair in the library at
the cascade's own 40-base operating point, so the question "does content vary with
position in the file" can be answered exactly, over the whole library, with no alignment
and no sampling of its own. The labels are a coarser instrument than the full cascade -
no Kraken2, no DIAMOND - but they carry the two classes that dominate every headline here,
ribosomal and identified, and the unplaced class that is the residue's upper bound.

THE CONTROL IS THE POINT, as always. The same statistics are computed on uniform random
samples of the same size drawn from the whole library. If a contiguous block reads like a
random sample, the block design is sound and this script says so with a number. If the
between-block variance exceeds the binomial expectation, the overdispersion factor it
reports is how much every interval in the repo should be widened, which is a correction
that can be applied without re-running anything.

WHAT IT CANNOT SAY. It is blind to any bias that does not change these five classes: a
block could hold a different mix of viral families at identical ribosomal content. And it
describes one library. A second library would need its own labels.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LABELS = Path("/mnt/project-files/shared/resasm/library_label_classes.npz")


def wilson(k: int, n: int, z: float = 1.96) -> dict:
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {"point": p, "lo": max(0.0, c - h), "hi": min(1.0, c + h)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    ap.add_argument("--pairs", type=int, default=100_000,
                    help="block size, matching the cascade's --pairs")
    ap.add_argument("--block", type=int, default=1_500_000,
                    help="the cascade's per-seed skip")
    ap.add_argument("--seeds", type=int, default=5,
                    help="how many cascade seeds to locate")
    ap.add_argument("--draws", type=int, default=200,
                    help="uniform random samples of the same size, as the control")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=ROOT / "results/block_position_bias")
    args = ap.parse_args()

    if not args.labels.exists():
        raise SystemExit(f"missing labels: {args.labels}")
    t0 = time.time()
    z = np.load(args.labels)
    codes = [str(c) for c in z["class_codes"]]
    # A pair is scored by mate 1 alone. The two mates are not interchangeable here either
    # - mate 2 is 93 nt against mate 1's 107 and maps at a different rate - and pooling
    # them would average two different measurements, which is the standing rule.
    cls = z["class1"]
    n_lib = cls.size
    rng = np.random.default_rng(args.seed)

    def fracs(idx_or_slice) -> np.ndarray:
        c = np.bincount(cls[idx_or_slice], minlength=len(codes))
        return c / c.sum()

    library = fracs(slice(None))

    # ---- every block of the library, end to end
    n_blocks = n_lib // args.pairs
    tiled = cls[:n_blocks * args.pairs].reshape(n_blocks, args.pairs)
    counts = np.stack([np.bincount(row, minlength=len(codes)) for row in tiled])
    block_f = counts / args.pairs

    # ---- the blocks the cascade actually uses
    used = []
    for s in range(1, args.seeds + 1):
        off = (s - 1) * args.block
        if off + args.pairs > n_lib:
            break
        used.append({"seed": s, "offset_pairs": off,
                     "fractions": dict(zip(codes, fracs(
                         slice(off, off + args.pairs)).round(6).tolist()))})

    # ---- the control: uniform random samples of the same size
    draws = np.stack([fracs(rng.choice(n_lib, size=args.pairs, replace=False))
                      for _ in range(args.draws)])

    # ---- dispersion. Binomial is what a Wilson interval assumes; the ratio of the
    # observed between-block variance to it is the factor every interval in the repo
    # would need widening by, and it is reported per class because the classes differ.
    out_classes = {}
    for i, c in enumerate(codes):
        p = float(library[i])
        binom_sd = float(np.sqrt(p * (1 - p) / args.pairs))
        obs_sd = float(block_f[:, i].std(ddof=1))
        ctl_sd = float(draws[:, i].std(ddof=1))
        first = float(block_f[0, i])
        used_mean = float(np.mean([u["fractions"][c] for u in used])) if used else None
        out_classes[c] = {
            "library_fraction": p,
            "block_sd_observed": obs_sd,
            "binomial_sd_expected": binom_sd,
            "overdispersion": obs_sd / binom_sd if binom_sd else None,
            "random_draw_sd": ctl_sd,
            "random_draw_sd_over_binomial": ctl_sd / binom_sd if binom_sd else None,
            "block_min": float(block_f[:, i].min()),
            "block_max": float(block_f[:, i].max()),
            "block_1_fraction": first,
            "block_1_z_against_blocks": (first - float(block_f[:, i].mean())) / obs_sd
            if obs_sd else None,
            "block_1_z_against_binomial": (first - p) / binom_sd if binom_sd else None,
            "cascade_seeds_mean": used_mean,
            "cascade_seeds_minus_library_points": (100 * (used_mean - p)
                                                   if used_mean is not None else None),
            "lag1_autocorrelation": float(np.corrcoef(block_f[:-1, i],
                                                      block_f[1:, i])[0, 1]),
        }

    out = {
        "parameters": {"labels": str(args.labels), "pairs_per_block": args.pairs,
                       "cascade_block_skip": args.block, "blocks_in_library": n_blocks,
                       "library_pairs": n_lib, "random_draws": args.draws,
                       "seed": args.seed,
                       "scored_by": "mate 1 only; the mates are not interchangeable and "
                                    "a pooled label would average two measurements",
                       "label_operating_point":
                           "SILVA 138.2 SSU+LSU and all 82,947 of this sample's own "
                           "contigs at >=40 matching bases, the cascade's own threshold; "
                           "no Kraken2 and no DIAMOND, so 'none' is an upper bound on the "
                           "residue rather than the residue"},
        "class_codes": codes,
        "library_fractions": dict(zip(codes, library.round(6).tolist())),
        "cascade_blocks": used,
        "classes": out_classes,
        "readout": "overdispersion is the factor by which the observed between-block "
                   "spread exceeds what a simple random sample of the same size would "
                   "give. 1.0 means a contiguous block IS a random sample for that class; "
                   "above 1.0 is how much every Wilson interval in this repo is too "
                   "narrow. random_draw_sd_over_binomial is the same statistic computed "
                   "on real random samples and must sit at 1.0, or the instrument is "
                   "wrong rather than the blocks",
        "wall_seconds": time.time() - t0,
        **provenance(__file__),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / f"block_position_bias_s{args.seed}.json"
    dest.write_text(json.dumps(out, indent=1))
    np.savez_compressed(args.out / f"block_fractions_s{args.seed}.npz",
                        block_fractions=block_f, class_codes=np.array(codes))
    for c, v in out_classes.items():
        print(f"{c:12s} library {100*v['library_fraction']:6.3f}%  "
              f"blocks {100*v['block_min']:6.3f}-{100*v['block_max']:6.3f}%  "
              f"overdispersion {v['overdispersion']:6.2f}x  "
              f"(control {v['random_draw_sd_over_binomial']:.2f}x)  "
              f"block1 z {v['block_1_z_against_binomial']:+7.1f}  "
              f"5 seeds off by {v['cascade_seeds_minus_library_points']:+.3f} pts")
    print(f"written {dest}")


if __name__ == "__main__":
    main()
