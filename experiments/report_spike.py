"""Aggregate the holdout spike's seeds into the table the report quotes.

Kept separate from the run so the tables can be regenerated without re-running anything,
and so the spread across seeds is computed in one place rather than in prose.

WHY THE PER-TIER TABLE IS NOW SPLIT BY GENOME. Two of the three spiked genomes are already
in this sample - pepper mild mottle virus assembles at 285x and human astrovirus as a
6,823 nt contig at 11.7x - and only human rhinovirus 1 is absent from it. `spike_holdout.py`
recorded that in a comment and measured all three separately, but this script pooled the
per-tier column across genomes, and the pooled number is the one the repo quoted.

It was wrong in a specific and predictable way. The sample's own assembly removed 100% and
96.4% of the two resident genomes at 0% divergence and **0.0% of the clean holdout at every
rung**, because a virus absent from the sample cannot be in the sample's assembly. So every
point the pooled table attributed to `selfasm_identified` or `cross_sample` was an artefact
of spiking material that was already there - including the 20.8% at 10% divergence that has
been quoted as evidence the sample's own data catches novel viruses. It does not, by
construction.

What survives the correction is the headline: the clean holdout's survival is
0.0 / 1.0 / 12.5 / 79.8 / 99.6 against the pooled 0.0 / 0.6 / 12.7 / 82.2 / 99.5, inside the
seed ranges at every rung. The own-assembly tier was redundant with the public indexes for
the resident genomes, so removing its contribution changes which tier gets the credit and
not how much is removed. For a genuinely novel virus the protein tier does essentially all
of the work above 5% divergence (91.3 / 84.8 / 20.0 / 0.4 on the holdout), which is the
"keep a protein tier" conclusion in a stronger form than the pooled table could state it.

A guard below refuses to write a table for a genome whose residency is not declared, because
the whole error was a status that lived in a comment rather than in the output.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Whether each spiked genome is already present in SRR23973193, which decides whether the
#: sample's own assembly can remove it. Declared here rather than in a comment, because the
#: pooled per-tier table was wrong for exactly as long as this lived in prose.
RESIDENCY = {
    "NC_003630.1": ("pepper mild mottle virus", "resident, 285x", False),
    "NC_001943.1": ("human astrovirus 1", "resident, 11.7x as a 6,823 nt contig", False),
    "NC_038311.1": ("human rhinovirus 1", "absent from this sample's assembly", True),
}


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
    undeclared = [g for g in genomes if g not in RESIDENCY]
    if undeclared:
        raise SystemExit(
            f"no residency declared for {undeclared}. Whether a spiked genome is already in "
            f"this sample decides whether the sample's own assembly can remove it, and a "
            f"pooled per-tier number is wrong by exactly that amount. Add it to RESIDENCY.")
    holdouts = [g for g in genomes if RESIDENCY[g][2]]
    lines = [f"Seeds: {len(ds)}; {ds[0]['spike_pairs_per_rung']:,} simulated pairs per "
             f"genome per rung.", "",
             "**Read the per-genome tier table, not the pooled one.** Of the spiked "
             "genomes, " + "; ".join(f"`{g}` {RESIDENCY[g][0]} is {RESIDENCY[g][1]}"
                                     for g in genomes) + ". A genome already in the sample "
             "is removed by the sample's own assembly for a reason that has nothing to do "
             "with detecting a novel virus, so only the holdout column "
             f"({', '.join(holdouts) or 'none'}) answers the question the spike is for.",
             "",
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

    tiers0 = list(ds[0]["rungs"][0]["genomes"][genomes[0]]["survival"]["removed_by_tier"])
    lines += ["## Which tier removes the spike, per genome", "",
              "The column that answers the spike's question is the holdout's. "
              "`selfasm_identified` and `cross_sample` read 0.0 on it at every rung, "
              "because a virus absent from the sample cannot be in the sample's assembly.",
              "", "| divergence | genome | " + " | ".join(tiers0) + " |",
              "|---|---|" + "---|" * len(tiers0)]
    for i, r0 in enumerate(ds[0]["rungs"]):
        for g in genomes:
            mark = " **(holdout)**" if RESIDENCY[g][2] else ""
            cells = [f"{st.mean([d['rungs'][i]['genomes'][g]['survival']['removed_by_tier'].get(t, 0.0) for d in ds]) * 100:.1f}"
                     for t in tiers0]
            lines.append(f"| {r0['divergence'] * 100:.0f}% | {g}{mark} | "
                         + " | ".join(cells) + " |")
    lines += ["", "## Which tier removes the spike, pooled over genomes "
                  "(**inflated — kept for comparison only**)", "",
              "Two of the three genomes are resident, so this table credits the sample's "
              "own assembly with removals that a novel virus would not have given it. It is "
              "the table this repo quoted before the split above existed.", ""]
    lines += ["",
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
              f"- Library recurrence is measured separately, with exactly the number of "
              f"spike pairs 1-in-a-million implies, in results/spike_abundance/."]

    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
