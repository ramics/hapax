"""What the §9 pair threshold actually does on real Berlin background.

`results/peptide_absence_noise/REPORT.md` §9 derived a pair threshold of 16 of 26 peptide
10-mers absent from the sequencing-error term alone, and §10 recorded that the organism
term would erode it. "Erode" turned out to be far too gentle a word, and this script is
the check that establishes by how much, so the number in the report is reproducible rather
than quoted from another thread's message.

It reads the per-pair scores measured by the mate-separation thread
(`experiments/mate_rho.py`, PR #10) against the viral + Swiss-Prot index of 221M peptide
10-mers: 40,000 real SRR23973193 pairs per seed, 5 seeds, mate identity kept, score =
fraction of a read's peptide 10-mers absent, minimum over the six frames. Every read in
that run is background, so every flag is a false flag.

Those .npz files live on branch `claude/project-thread-clcuow`. Extract them with
`git show <branch>:results/mate_rho/viral_sprot/pair_scores_s<N>.npz` and point --scores at
the directory; the paths are recorded in the manifest so the provenance survives the
branch being merged or renamed.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics as st
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="unused: the inputs are fixed files, but every experiment takes "
                         "--seed as the repo requires; the 5 input seeds are read in full")
    ap.add_argument("--scores", type=Path, required=True,
                    help="directory holding pair_scores_s1..s5.npz from PR #10")
    ap.add_argument("--label", default="unlabelled",
                    help="which index arm these scores came from, recorded in the output")
    ap.add_argument("--k-aa", type=int, default=10)
    ap.add_argument("--kmers-per-read", type=int, default=26,
                    help="26 peptide 10-mers in a 107 nt read at frame offset 1 or 2")
    ap.add_argument("--thresholds", default="16,21,26",
                    help="16 is the §9 error-only pair threshold; 21 and 26 are the §5 "
                         "single-read thresholds for 100/M and 1/M")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "peptide_absence_noise")
    args = ap.parse_args()

    t0 = time.perf_counter()
    files = sorted(args.scores.glob("*.npz"))
    if not files:
        raise SystemExit(f"no .npz under {args.scores}")

    rows = []
    sat_single, sat_pair, corr = [], [], []
    for f in files:
        d = np.load(f)
        a, b = d["score1"], d["score2"]
        sat_single.append(float(np.mean(np.concatenate([a, b]) >= 1.0)))
        sat_pair.append(float(np.mean((a >= 1.0) & (b >= 1.0))))
        corr.append(float(np.corrcoef(a, b)[0, 1]))

    out = {
        "source": str(args.scores),
        "source_branch": "claude/project-thread-clcuow (PR #10)",
        "source_files": [f.name for f in files],
        "index_arm": args.label,
        "pairs_per_seed": int(np.load(files[0])["score1"].shape[0]),
        "seeds": len(files),
        "background_fully_absent_single_read": _sp(sat_single),
        "background_fully_absent_both_mates": _sp(sat_pair),
        "mate_correlation_unconditional": _sp(corr),
        "thresholds": [],
    }

    for T in [int(x) for x in args.thresholds.split(",")]:
        frac = T / args.kmers_per_read
        single, both = [], []
        for f in files:
            d = np.load(f)
            a, b = d["score1"], d["score2"]
            single.append(float(np.mean(a >= frac)))
            both.append(float(np.mean((a >= frac) & (b >= frac))))
        s_mean, b_mean = st.mean(single), st.mean(both)
        rows.append(T)
        out["thresholds"].append({
            "threshold_kmers": T,
            "score_fraction": frac,
            "single_read_rate_per_million": s_mean * 1e6,
            "both_mates_rate_per_million": b_mean * 1e6,
            "independence_would_give_per_million": (s_mean ** 2) * 1e6,
            # NOT a gain: >1 means the observed both-mates rate EXCEEDS what independent
            # mates would give, i.e. correlation working against the conjunction.
            "excess_over_independence": (b_mean / (s_mean ** 2)) if s_mean else None,
            "over_budget_factor_vs_1_per_million": b_mean * 1e6,
        })

    (args.out / f"absence_threshold_check_{args.label}.json").write_text(
        json.dumps(out, indent=2))
    (args.out / f"absence_threshold_check_{args.label}_manifest.json").write_text(json.dumps({
        "script": "experiments/absence_threshold_check.py",
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))

    print(f"background fully absent, single read: {out['background_fully_absent_single_read']['mean']:.4f}")
    print(f"background fully absent, both mates : {out['background_fully_absent_both_mates']['mean']:.4f}")
    print(f"mate correlation, unconditional     : {out['mate_correlation_unconditional']['mean']:.4f}")
    for t in out["thresholds"]:
        print(f"T>={t['threshold_kmers']:>2}: single {t['single_read_rate_per_million']:>12,.0f}/M  "
              f"both mates {t['both_mates_rate_per_million']:>12,.0f}/M  "
              f"(independence {t['independence_would_give_per_million']:>12,.0f}/M, "
              f"excess {t['excess_over_independence']:.2f}x)")


def _sp(xs: list[float]) -> dict:
    return {"mean": st.mean(xs), "min": min(xs), "max": max(xs),
            "sd": st.stdev(xs) if len(xs) > 1 else 0.0}


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
