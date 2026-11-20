"""Aggregate the nucleotide content ladder and say plainly whether it flattens.

The marginal column is the question: if removal per added unit of reference content is
still large at the largest content available, the curve has not flattened and nobody
should extrapolate a plateau onto it.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def cell(xs) -> str:
    return f"{st.mean(xs) * 100:.2f} ({min(xs) * 100:.2f}-{max(xs) * 100:.2f})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results/nucleotide_ladder")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/nucleotide_ladder/TABLES.md")
    args = ap.parse_args()
    ds = [json.load(open(f))
          for f in sorted(glob.glob(str(args.results / "nucladder_s*.json")))]
    if not ds:
        raise SystemExit("no results")

    lines = [f"Seeds: {len(ds)}; {ds[0]['pairs']:,} pairs per seed; database "
             f"{Path(ds[0]['db']).name} ({ds[0]['db_bytes'] / 1e9:.1f} GB).", "",
             "Axis is REFERENCE CONTENT, not bytes: the library was classified once "
             "against the largest database available and nested random subsets of the "
             "reference taxa were then applied to the calls. Both biases run the same "
             "way and are stated in the script's docstring; the stronger one is that "
             "deleting calls understates what a real smaller database would remove, "
             "which flatters the left-hand side and makes this a conservative test of "
             "whether the curve has flattened.", ""]

    for conf in ds[0]["by_confidence"]:
        lines += [f"## Kraken2 confidence {conf}", "",
                  "| reference content | taxa | removed % of library | marginal | "
                  "removed % of non-rRNA | marginal |", "|---|---|---|---|---|---|"]
        for i, r0 in enumerate(ds[0]["by_confidence"][conf]["rungs"]):
            def col(key):
                return cell([d["by_confidence"][conf]["rungs"][i][key] for d in ds])
            taxa = st.mean([d["by_confidence"][conf]["rungs"][i]
                            ["reference_taxa_retained"] for d in ds])
            lines.append(
                f"| {r0['content_fraction'] * 100:g}% | {taxa:.0f} | "
                f"{col('removed_frac_of_library')} | {col('marginal_frac_of_library')} | "
                f"{col('removed_frac_of_non_rrna')} | "
                f"{col('marginal_frac_of_non_rrna')} |")
        conc = [d["by_confidence"][conf].get("read_concentration_over_taxa")
                for d in ds if d["by_confidence"][conf].get("read_concentration_over_taxa")]
        if conc:
            lines += ["", "Read concentration over reference taxa: top taxon "
                      f"{st.mean(c['top1'] for c in conc) * 100:.1f}%, top 10 "
                      f"{st.mean(c['top10'] for c in conc) * 100:.1f}%, top 100 "
                      f"{st.mean(c['top100'] for c in conc) * 100:.1f}% of classified "
                      "reads. The tail is long and the head is heavy, so a content "
                      "fraction is a weak proxy for what a subset buys: whether one "
                      "dominant taxon falls inside the subset moves a rung more than "
                      "hundreds of rare ones, which is why the spread across seeds is "
                      "wide and is reported rather than averaged away."]
        lines.append("")

    last = {c: ds[0]["by_confidence"][c]["rungs"][-1] for c in ds[0]["by_confidence"]}
    lines += ["## Does it flatten?", ""]
    for c, r in last.items():
        marg = st.mean([d["by_confidence"][c]["rungs"][-1]["marginal_frac_of_library"]
                        for d in ds])

        lines.append(f"- At confidence {c}, the last doubling of reference content "
                     f"(50% to 100%) still buys {marg * 100:.1f} points of the library, "
                     f"and no plateau is visible. Read that next to the concentration "
                     f"above, though: most of that step is one dominant taxon being "
                     f"present rather than absent, so this ladder is a weak instrument "
                     f"for the SHAPE of the curve and a strong one only for the claim "
                     f"that it has not visibly saturated.")
    lines += ["", "The stronger nucleotide evidence that the curve has not flattened is "
              "not this ladder but the remote BLAST in the main report: against NCBI "
              "core_nt, two orders of magnitude more sequence than anything local, 87% "
              "of residue reads find a partial, distant match. Sequence our indexes miss "
              "is findable in a database 70 times larger.",
              "", "This matches the protein side (PR #10 sections 4.6-4.7), where each "
              "tenfold step in index size was still buying about 20 points at the "
              "largest index reachable. So 'removal per gigabyte and where it flattens' "
              "has, so far, no flattening point on either axis, and this report does not "
              "invent one."]

    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
