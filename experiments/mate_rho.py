"""The mate correlation of an absence-based novelty score, as a function of how
far apart the mates are, on real wastewater reads.

This is the measurement the paired-end case reduces to. From
results/readlevel_features/REPORT.md §3.1: every pair-level rule's merit is set
by the mate correlation rho, and the merit of requiring both mates to pass swings
four orders of magnitude between rho_fp = 0 and rho_fp = 0.1. The earlier
read-level arm discarded mate identity by scoring reads as an unpaired pool, so
rho has never been measured here. Insert length is the one knob that moves it,
which is why this sits in the wide-insert thread: if rho_fp is already near zero
at the library's natural insert, width buys nothing; if it falls with separation
inside the range the material can deliver, width is a real lever.

Design. Every read in SRR23973193 is background (no spike), so a flag on any of
these reads is a false flag and the correlation measured between mates is
directly rho_fp. The score is six-frame peptide k-mer absence
(experiments/peptide_absence.py), the one single-read feature with novelty
polarity. Three separations, two of them controls that bracket the real one:

- **gap 0**: the two halves of a single mate. Perfectly adjacent sequence, so
  this is rho at zero separation and the ceiling the real pairs are compared to.
- **the real gap**: mate 1 against mate 2, binned by the pair's measured insert
  length, which gives rho as a function of separation over the range the library
  actually delivers (median 194 nt, p99 485 - results/insert_profile/).
- **gap infinity**: mate 1 of one pair against mate 2 of another. Independent
  molecules from different organisms, so this is the rho = 0 floor, and it also
  catches any correlation that is an artefact of the scoring rather than of the
  biology.

Both correlations are reported: Pearson r on the continuous score, and the phi
coefficient of the thresholded flag, which is what the model in the feature map
calls rho. Thresholds are percentiles of the score's own distribution, because an
absence score against a viral-protein index has no absolute calibration - and
that is a limitation to state rather than hide: with a complete known-protein
index the absolute absence rates would be lower, so the levels here are not a
detector's false-positive rate. The quantity this script is for is how the
correlation depends on separation, which the index choice shifts far less than it
shifts the levels.
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

from insert_profile import read_pairs, revcomp  # same trimming and streaming
from peptide_absence import Scorer, build_index


def phi(a: np.ndarray, b: np.ndarray) -> float | None:
    """Correlation of two binary flags: the model's rho."""
    pa, pb = a.mean(), b.mean()
    den = (pa * (1 - pa) * pb * (1 - pb)) ** 0.5
    if den == 0:
        return None
    return float(((a * b).mean() - pa * pb) / den)


def pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    if a.size < 3 or a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def flag_stats(s1: np.ndarray, s2: np.ndarray, qs: list[float]) -> dict:
    """rho and the pair rates at several operating points."""
    out = {}
    pool = np.concatenate([s1, s2])
    for q in qs:
        thr = float(np.quantile(pool, q))
        f1 = (s1 >= thr).astype(float)
        f2 = (s2 >= thr).astype(float)
        p = float((f1.mean() + f2.mean()) / 2)
        both = float((f1 * f2).mean())
        either = float(((f1 + f2) > 0).mean())
        out[f"q{q:g}"] = {
            "threshold": round(thr, 4),
            "marginal_flag_rate": round(p, 5),
            "rho_flags": round(phi(f1, f2), 4) if phi(f1, f2) is not None else None,
            "p_both": round(both, 6),
            "p_both_if_independent": round(p * p, 6),
            "both_excess_x": round(both / (p * p), 2) if p > 0 else None,
            "p_either": round(either, 5),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    ap.add_argument("--ref", type=Path, required=True,
                    help="contigs to measure each pair's insert against")
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--block", type=int, default=2_000_000)
    ap.add_argument("--trim-r2", type=int, default=14)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--faa", type=Path, nargs="+",
                    default=[Path("data/prot/viral.1.protein.faa.gz")])
    ap.add_argument("--gpff", type=Path, default=Path("data/prot/viral.1.protein.gpff.gz"))
    ap.add_argument("--cache", type=Path, default=Path("data/pepidx"))
    ap.add_argument("--quantiles", type=float, nargs="+", default=[0.9, 0.99, 0.999])
    ap.add_argument("--out", type=Path, default=Path("results/mate_rho"))
    args = ap.parse_args()

    rng = random.Random(args.seed)
    np.random.seed(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    idx = build_index(args.faa, args.k, args.cache, None, args.gpff)
    scorer = Scorer(idx, args.k)

    import mappy
    aln = mappy.Aligner(str(args.ref), preset="sr", n_threads=4)
    if not aln:
        raise SystemExit(f"could not index {args.ref}")

    skip = (args.seed - 1) * args.block
    rows = []           # one per pair: insert (or None), score1, score2
    half_rows = []      # gap-0 control: the two halves of mate 1
    disp = Counter()
    for name, r1, r2 in read_pairs(args.r1, args.r2, args.pairs, args.trim_r2, skip):
        disp["pairs"] += 1
        s1, n1, _ = scorer.score_nt_best_frame(r1)
        s2, n2, _ = scorer.score_nt_best_frame(r2)
        if not (n1 and n2):
            disp["unscorable"] += 1
            continue
        h = len(r1) // 2
        ha, na, _ = scorer.score_nt_best_frame(r1[:h])
        hb, nb, _ = scorer.score_nt_best_frame(r1[h:])
        if na and nb:
            half_rows.append((ha, hb))
        ins = None
        a = [x for x in aln.map(r1) if x.is_primary]
        b = [x for x in aln.map(r2) if x.is_primary]
        if a and b and a[0].ctg == b[0].ctg and a[0].strand != b[0].strand:
            span = max(a[0].r_en, b[0].r_en) - min(a[0].r_st, b[0].r_st)
            span += a[0].q_st + (len(r1) - a[0].q_en) + b[0].q_st + (len(r2) - b[0].q_en)
            if 0 < span < 2000:
                ins = span
                disp["insert_measured"] += 1
        rows.append((ins, s1, s2))

    S1 = np.array([r[1] for r in rows])
    S2 = np.array([r[2] for r in rows])
    res = {
        "seed": args.seed,
        "k": args.k,
        "index_kmers": int(np.asarray(idx).size),
        "pairs": dict(disp),
        "score_distribution": {
            "mate1": {"mean": round(float(S1.mean()), 4), "median": round(float(np.median(S1)), 4),
                      "p90": round(float(np.quantile(S1, 0.9)), 4),
                      "p99": round(float(np.quantile(S1, 0.99)), 4),
                      "fraction_1_0": round(float((S1 >= 0.999).mean()), 4)},
            "mate2": {"mean": round(float(S2.mean()), 4), "median": round(float(np.median(S2)), 4),
                      "p90": round(float(np.quantile(S2, 0.9)), 4),
                      "p99": round(float(np.quantile(S2, 0.99)), 4),
                      "fraction_1_0": round(float((S2 >= 0.999).mean()), 4)},
        },
        "separations": {},
    }

    # gap 0: the halves of one mate
    if half_rows:
        HA = np.array([h[0] for h in half_rows])
        HB = np.array([h[1] for h in half_rows])
        res["separations"]["gap_0_halves_of_mate1"] = {
            "n": int(HA.size),
            "note": "two 53 nt halves of the same read: adjacent sequence, rho ceiling",
            "r_continuous": round(pearson(HA, HB), 4) if pearson(HA, HB) is not None else None,
            "flags": flag_stats(HA, HB, args.quantiles),
        }

    # the real gap, binned by measured insert
    bins = [(107, 150), (150, 200), (200, 250), (250, 300), (300, 400), (400, 600)]
    by_bin = defaultdict(list)
    for ins, s1, s2 in rows:
        if ins is None:
            continue
        for lo, hi in bins:
            if lo <= ins < hi:
                by_bin[(lo, hi)].append((ins, s1, s2))
                break
    for (lo, hi), vals in sorted(by_bin.items()):
        a = np.array([v[1] for v in vals])
        b = np.array([v[2] for v in vals])
        if a.size < 200:
            continue
        gaps = [v[0] - 107 - 93 for v in vals]
        res["separations"][f"insert_{lo}_{hi}"] = {
            "n": int(a.size),
            "mean_insert": round(float(np.mean([v[0] for v in vals])), 1),
            "mean_gap_between_mates": round(float(np.mean(gaps)), 1),
            "r_continuous": round(pearson(a, b), 4) if pearson(a, b) is not None else None,
            "flags": flag_stats(a, b, args.quantiles),
        }

    # all measured pairs pooled, and all pairs including unmapped
    for label, sel in (("all_insert_measured", [r for r in rows if r[0] is not None]),
                       ("all_pairs", rows)):
        a = np.array([r[1] for r in sel])
        b = np.array([r[2] for r in sel])
        res["separations"][label] = {
            "n": int(a.size),
            "r_continuous": round(pearson(a, b), 4) if pearson(a, b) is not None else None,
            "flags": flag_stats(a, b, args.quantiles),
        }

    # gap infinity: mate 1 of one pair against mate 2 of another
    perm = list(range(len(rows)))
    rng.shuffle(perm)
    a = np.array([rows[i][1] for i in range(len(rows))])
    b = np.array([rows[perm[i]][2] for i in range(len(rows))])
    res["separations"]["gap_infinity_shuffled_pairing"] = {
        "n": int(a.size),
        "note": "different molecules, different organisms: the rho = 0 floor",
        "r_continuous": round(pearson(a, b), 4) if pearson(a, b) is not None else None,
        "flags": flag_stats(a, b, args.quantiles),
    }

    (args.out / f"mate_rho_s{args.seed}.json").write_text(json.dumps(res, indent=2))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/mate_rho.py", "git_sha": sha, "seed": args.seed,
        "params": {k: str(v) for k, v in vars(args).items()},
        "wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2))
    for k, v in res["separations"].items():
        print(k, v["n"], "r=", v["r_continuous"],
              {q: (d["rho_flags"], d["both_excess_x"]) for q, d in v["flags"].items()})


if __name__ == "__main__":
    main()
