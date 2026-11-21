"""Aggregate the subtraction cascade across seeds into the tables the report quotes.

Separate from the experiment so the report's numbers come from the committed JSON rather
than from prose, which is the failure mode the project's process note of 2026-11-20 names
first: check report prose against committed outputs. Every table this prints carries the
spread across seeds, because the repo rule is five seeds and a mean with no spread hides
whether the five agree.
"""

from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sp(xs):
    xs = list(xs)
    if not xs:
        return None
    return {"mean": st.mean(xs), "min": min(xs), "max": max(xs),
            "sd": st.stdev(xs) if len(xs) > 1 else 0.0, "n": len(xs)}


def fmt(s, pct=True, dp=2):
    if s is None:
        return "-"
    m = s["mean"] * (100 if pct else 1)
    lo = s["min"] * (100 if pct else 1)
    hi = s["max"] * (100 if pct else 1)
    return f"{m:.{dp}f} ({lo:.{dp}f}-{hi:.{dp}f})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1, help="unused; inputs are fixed files")
    ap.add_argument("--results", type=Path, default=ROOT / "results/subtract_known")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    runs = [json.loads(p.read_text())
            for p in sorted(args.results.glob("subtract_s*.json"))]
    if not runs:
        raise SystemExit("no runs")
    tiers = [c["tier"] for c in runs[0]["curve"]]

    lines = [f"Seeds: {len(runs)}; pairs per seed: {runs[0]['pairs']:,}",
             "",
             "## Cascade, cheapest-first (mean % of pairs, min-max across seeds)",
             "",
             "| tier | index GB | alone | marginal | cumulative | residue | wall s |",
             "|---|---|---|---|---|---|---|"]
    for i, tier in enumerate(tiers):
        rows = [r["curve"][i] for r in runs]
        lines.append(
            f"| {tier} | {rows[0]['index_bytes'] / 1e9:.3f} | "
            f"{fmt(sp(x['tier_alone_removed'] for x in rows))} | "
            f"{fmt(sp(x['marginal_removed'] for x in rows))} | "
            f"{fmt(sp(x['cumulative_removed'] for x in rows))} | "
            f"{fmt(sp(x['residue'] for x in rows))} | "
            f"{st.mean([x['wall_seconds'] or 0 for x in rows]):.1f} |")
    lines += ["",
              f"Union removed: {fmt(sp(r['union_removed'] for r in runs))} %",
              f"Residue pairs per seed: "
              f"{fmt(sp(r['residue_pairs'] / r['pairs'] for r in runs))} %",
              ""]

    lines += ["## Threshold sweeps (mean % of pairs removed, min-max)", ""]
    for tier in tiers:
        metas = [r["tiers"][tier] for r in runs if tier in r["tiers"]]
        for key in ("confidence_sweep", "mlen_sweep", "bitscore_sweep_of_library",
                    "entropy_sweep"):
            if key in metas[0] and metas[0][key]:
                keys = list(metas[0][key])
                cells = " | ".join(fmt(sp(m[key][k] for m in metas)) for k in keys)
                lines.append(f"**{tier}** ({key}): " + " | ".join(keys))
                lines.append("  " + cells)
                lines.append("")

    lines += ["## Removal per gigabyte, in cascade order", "",
              "| tier | cumulative GB | cumulative removed % | marginal % per GB |",
              "|---|---|---|---|"]
    for i, tier in enumerate(tiers):
        rows = [r["curve"][i] for r in runs]
        gb = rows[0]["cumulative_index_bytes"] / 1e9
        marg = st.mean([x["marginal_removed"] for x in rows]) * 100
        own_gb = rows[0]["index_bytes"] / 1e9
        per = marg / own_gb if own_gb > 0.0005 else float("inf")
        lines.append(f"| {tier} | {gb:.3f} | "
                     f"{st.mean([x['cumulative_removed'] for x in rows]) * 100:.2f} | "
                     f"{per:,.1f} |")

    text = "\n".join(lines)
    print(text)
    if args.out:
        args.out.write_text(text + "\n")


if __name__ == "__main__":
    main()
