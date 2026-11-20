"""Aggregate the holdout spike's seeds into the table the report quotes.

Kept separate from the run so the tables can be regenerated without re-running anything,
and so the spread across seeds is computed in one place rather than in prose.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fmt(xs) -> str:
    xs = list(xs)
    if not xs:
        return "-"
    m = st.mean(xs)
    return f"{m * 100:.1f} ({min(xs) * 100:.1f}-{max(xs) * 100:.1f})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results/spike_holdout")
    ap.add_argument("--out", type=Path, default=ROOT / "results/spike_holdout/TABLES.md")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.results / "spike_s*.json")))
    if not files:
        raise SystemExit(f"no spike results in {args.results}")
    ds = [json.load(open(f)) for f in files]
    genomes = list(ds[0]["rungs"][0]["genomes"])
    lines = [f"Seeds: {len(ds)}; {ds[0]['spike_pairs_per_rung']:,} simulated pairs per "
             f"genome per rung; genomes: {', '.join(genomes)}.", "",
             "## Spike survival through tier 1, by divergence",
             "", "| divergence | " + " | ".join(genomes) + " | all |",
             "|---|" + "---|" * (len(genomes) + 1)]
    for i, r0 in enumerate(ds[0]["rungs"]):
        cells, pooled = [], []
        for g in genomes:
            vals = [d["rungs"][i]["genomes"][g]["survival"]["survives_tier1"] for d in ds]
            pooled += vals
            cells.append(fmt(vals))
        lines.append(f"| {r0['divergence'] * 100:.0f}% | " + " | ".join(cells)
                     + f" | {fmt(pooled)} |")
    lines += ["", "Mean % of simulated pairs surviving every index the cascade subtracts "
                  "(min-max across seeds).", ""]

    lines += ["## Which tier removes the spike", "",
              "| divergence | " + " | ".join(ds[0]["rungs"][0]["genomes"][genomes[0]]
                                             ["survival"]["removed_by_tier"]) + " |",
              "|---|" + "---|" * len(ds[0]["rungs"][0]["genomes"][genomes[0]]["survival"]
                                     ["removed_by_tier"])]
    tiers = list(ds[0]["rungs"][0]["genomes"][genomes[0]]["survival"]["removed_by_tier"])
    for i, r0 in enumerate(ds[0]["rungs"]):
        cells = []
        for t in tiers:
            vals = [d["rungs"][i]["genomes"][g]["survival"]["removed_by_tier"].get(t, 0.0)
                    for d in ds for g in genomes]
            cells.append(f"{st.mean(vals) * 100:.1f}")
        lines.append(f"| {r0['divergence'] * 100:.0f}% | " + " | ".join(cells) + " |")

    a = ds[0]["one_in_a_million_arithmetic"]
    lines += ["", "## At the real abundance", "",
              f"- A virus at 1 in a million contributes about "
              f"{a['expected_spike_pairs_in_run']:.1f} read pairs in this "
              f"{a['run_pairs']:,}-pair run.",
              f"- Per-read survival above is measured at "
              f"{ds[0]['spike_pairs_per_rung']:,} pairs per rung, because a simulation at "
              f"the real rate cannot resolve the quantity; the sentence above is "
              f"arithmetic on top of it.",
              f"- Library recurrence measured with exactly the number of spike pairs "
              f"1-in-a-million implies: "
              f"{ds[0]['abundance_at_real_rate']['spike_pairs_in_sketch']} pairs in "
              f"{ds[0]['abundance_at_real_rate']['library_pairs_in_sketch']:,}, median "
              f"31-mer count "
              f"{(ds[0]['abundance_at_real_rate']['spike_median_31mer_count'] or {}).get('mean')}"
              f" against a background whose median is reported in the fields run."]

    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
