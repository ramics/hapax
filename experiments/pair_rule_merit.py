"""What a pair rule can do at one false positive per million, given the measured
correlations.

Imo fixed the operating point on 2026-11-20: **1 false flag per million reads**,
reached directly rather than by loosening the threshold and leaning on grouping.
That makes the pair rule's arithmetic load-bearing, and it is now arithmetic with
measured inputs rather than assumed ones:

- `rho_fp`, the mate correlation of FALSE flags, measured on real Berlin pairs
  (`results/mate_rho/`), by insert bin;
- `rho_true`, the mate correlation of NOVELTY, measured per insert length on
  four family holdouts (`results/novelty_blocks/`);
- `p_novel_either`, the chance a pair has at least one mate in novel sequence,
  from the same runs;
- the per-read operating point, which no measurement here fixes, so it is swept.

The two rules, for a marginal per-read rate p and correlation rho:

    P(both pass)   = p^2 + rho * p * (1 - p)
    P(either pass) = 2p - p^2 - rho * p * (1 - p)

The conjunction suppresses false flags only to the extent the mates are
independent; at rho_fp = 0.9 it barely helps, which is the measured regime. The
disjunction gains sensitivity as rho_true falls, which is what a wide insert
delivers. Both effects are read off the same measured curves here, so the output
is a trade-off table rather than a claim.

This script computes; it measures nothing. Its inputs are the committed JSONs of
the runs above, so every row traces to a measurement or to the swept per-read
rate.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from collections import defaultdict
from pathlib import Path


def both(p: float, rho: float) -> float:
    return p * p + rho * p * (1 - p)


def either(p: float, rho: float) -> float:
    return 2 * p - p * p - rho * p * (1 - p)


def load_novelty(dirs: Path) -> dict[str, dict[int, dict[str, float]]]:
    out: dict[str, dict[int, dict[str, float]]] = {}
    for f in sorted(glob.glob(str(dirs / "novelty_blocks_*_s?.json"))):
        d = json.load(open(f))
        acc = d["accession"]
        agg = out.setdefault(acc, {})
        for I, v in d["sweep"].items():
            row = agg.setdefault(int(I), defaultdict(list))
            row["rho_true"].append(v["rho_novelty"])
            row["p_either_novel"].append(v["p_either_novel"])
    return {a: {I: {k: round(st.fmean(v), 4) for k, v in row.items()}
                for I, row in sorted(rows.items())} for a, rows in out.items()}


def load_rho_fp(path: Path) -> dict[str, float]:
    d = json.load(open(path))
    return {k: v["r_continuous"] for k, v in d["separations"].items()
            if v.get("r_continuous") is not None}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--novelty-dir", type=Path, default=Path("results/novelty_blocks"))
    ap.add_argument("--rho-fp", type=Path, default=None,
                    help="a mate_rho_s*.json; without it, rho_fp is swept")
    ap.add_argument("--target-fp-per-million", type=float, default=1.0)
    ap.add_argument("--per-read-fp", type=float, nargs="+",
                    default=[1e-6, 1e-5, 1e-4, 1e-3])
    ap.add_argument("--per-read-sens", type=float, default=0.5,
                    help="per-read sensitivity on a read that IS in novel sequence")
    ap.add_argument("--inserts", type=int, nargs="+", default=[107, 200, 300, 500, 700])
    ap.add_argument("--out", type=Path, default=Path("results/mate_width/pair_rule_merit.json"))
    args = ap.parse_args()

    nov = load_novelty(args.novelty_dir)
    rho_fp_measured = load_rho_fp(args.rho_fp) if args.rho_fp else {}
    # map an insert length to the measured bin it falls in
    bin_for = {107: "insert_107_150", 200: "insert_200_250", 300: "insert_300_400",
               500: "insert_400_600", 700: "insert_400_600"}

    res = {"per_read_sens_on_novel_read": args.per_read_sens,
           "rho_fp_source": str(args.rho_fp) if args.rho_fp else "swept",
           "rho_fp_measured": rho_fp_measured,
           "rows": []}

    for acc, rows in nov.items():
        for I in args.inserts:
            if I not in rows:
                continue
            r_true = rows[I]["rho_true"]
            p_nov = rows[I]["p_either_novel"]
            rfp = rho_fp_measured.get(bin_for.get(I, ""), None)
            for p_f in args.per_read_fp:
                for rho_fp in ([rfp] if rfp is not None else [0.0, 0.25, 0.5, 0.9]):
                    if rho_fp is None:
                        continue
                    # sensitivity: the pair must contain a novel mate AND that mate
                    # must pass. Under the conjunction both mates must pass, which
                    # needs both to be novel - the quantity measured as p_both.
                    s_or = p_nov * either(args.per_read_sens, r_true)
                    s_and = p_nov * both(args.per_read_sens, r_true)
                    fp_or = either(p_f, rho_fp)
                    fp_and = both(p_f, rho_fp)
                    res["rows"].append({
                        "accession": acc, "insert": I,
                        "rho_true_novelty": r_true, "p_either_novel": p_nov,
                        "per_read_fp": p_f, "rho_fp": round(rho_fp, 4),
                        "either_rule": {
                            "pair_sens": round(s_or, 4),
                            "pair_fp_per_million": round(fp_or * 1e6, 3),
                            "meets_target": fp_or * 1e6 <= args.target_fp_per_million,
                        },
                        "both_rule": {
                            "pair_sens": round(s_and, 4),
                            "pair_fp_per_million": round(fp_and * 1e6, 3),
                            "meets_target": fp_and * 1e6 <= args.target_fp_per_million,
                            "fp_reduction_vs_single_read": round(p_f / fp_and, 2)
                            if fp_and > 0 else None,
                            "merit": round((p_f / fp_and) / (args.per_read_sens / (s_and / p_nov)), 2)
                            if fp_and > 0 and s_and > 0 else None,
                        },
                    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(res, indent=2))
    hdr = f"{'genome':14s}{'ins':>5s}{'rho_t':>7s}{'rho_fp':>8s}{'readFP':>9s}" \
          f"{'OR sens':>9s}{'OR FP/M':>9s}{'AND sens':>9s}{'AND FP/M':>9s}{'merit':>7s}"
    print(hdr)
    for r in res["rows"]:
        print(f"{r['accession']:14s}{r['insert']:5d}{r['rho_true_novelty']:7.3f}"
              f"{r['rho_fp']:8.3f}{r['per_read_fp']:9.0e}"
              f"{r['either_rule']['pair_sens']:9.3f}{r['either_rule']['pair_fp_per_million']:9.2f}"
              f"{r['both_rule']['pair_sens']:9.3f}{r['both_rule']['pair_fp_per_million']:9.2f}"
              f"{str(r['both_rule']['merit']):>7s}")


if __name__ == "__main__":
    main()
