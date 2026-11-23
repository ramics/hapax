#!/usr/bin/env python3
"""Draw the README overview figure from a table of measured quantities.

Why this is a generator and not a hand-drawn SVG: CLAUDE.md's rule is that a number
typed into a figure is a number without provenance, and the first cut of this figure
had nine of them typed straight into the drawing. Here every quantity appears exactly
once, in STAGES or LIBRARY_SPLIT below, beside the results directory it was read from,
and the geometry of the one to-scale element (panel b) is computed from the percentages
rather than measured off by eye -- so a figure that disagrees with the record is a
one-line diff in a table rather than an arithmetic error hidden in path coordinates.

What this script does NOT do, stated rather than implied: it does not parse the result
files. The values are transcribed from the committed REPORT.md of each directory, which
is why each carries its `source` and why `--check-sources` (on by default) refuses to
draw when a cited directory is absent from the tree. Parsing them out of the manifests
is the obvious next step and is not done: the five directories cited here write five
different shapes, and a parser per shape is a bigger change than the figure.

Usage:
    python3 experiments/make_overview_figure.py --out docs/hapax-overview.svg
    python3 experiments/make_overview_figure.py --check   # verify committed == generated
"""

from __future__ import annotations

import argparse
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------------------
# The record. One row per stage of `hapax run`, each with one input, one output and one
# measured quantity. The stage names are the driver's own (SPEC 13.1), grouped where one
# line of the figure covers several of its stages -- `score the distance` covers
# `divergent`, `shortlist` and `reads`; `report the residue` covers `census` and `layers`.
# --------------------------------------------------------------------------------------
STAGES = [
    {
        "name": "resolve",
        "io": "library accession → trim, block size, which panels apply",
        "number": "14 nt",
        "unit": "of non-template prefix on mate 2; the mode in 5 of 5 seeds",
        "source": "results/nontemplate_prefix",
        "accent": "open",
    },
    {
        "name": "merge",
        "io": "read pairs → one read per overlapping molecule",
        "number": "35.1%",
        "unit": "of pairs merge, which lifts recognition at 20% divergence from 10.5% to 38.4%",
        "source": "results/profile_reach",
        "accent": "open",
    },
    {
        "name": "assemble",
        "io": "the whole library → contigs, checked for recovery in the same pass",
        "number": "20 of 25",
        "unit": "depositor genomes recovered at 90% breadth or better; 53 zero-read controls stay empty",
        "source": "results/path_settings",
        "accent": "open",
    },
    {
        "name": "identify",
        "io": "reads and contigs against indexes, cheapest first → tier table and the evidence for each removal",
        "number": "77.10%",
        "unit": "of read pairs removed as already known; 77.02 to 77.28 over 5 seeds",
        "source": "results/subtract_recut",
        "accent": "known",
    },
    {
        "name": "score the distance",
        "io": "every row, removed or not → distance, closest relative, query coverage, E-value",
        "number": "0.00%",
        "unit": "on shuffled sequence, reported next to every field it belongs to",
        "source": "results/divergent_contigs",
        "accent": "known",
    },
    {
        "name": "report the residue",
        "io": "what nothing named → one row per pair, no threshold applied, each field beside its control",
        "number": "depth 1",
        "unit": "median depth for pairs that never assembled; set by read length and database reach",
        "source": "results/residue_census",
        "accent": "residue",
    },
    {
        "name": "subtract other plants",
        "io": "512 libraries from other plants → witnessed in the panel or not, paired by plant",
        "number": "87.5%",
        "unit": "of a held-out library removed, and still falling at 512 libraries",
        "source": "results/crossplant",
        "accent": "open",
    },
]

# Panel b, to scale. `points` are percentage points of read pairs; they are the marginal
# column of the cascade table, so they sum to the union and the bar needs no fudging.
# The residue's two tiers are quoted to one decimal in the committed table and do not sum
# to the residue exactly, so they are written as labels under the block and not as
# divisions of it -- a division would assert a precision the record does not have.
LIBRARY_SPLIT = {
    "source": "results/subtract_recut",
    "identified": [
        ("ribosomal, 37.93 points", 37.93),
        ("the sample's own contigs, 30.33", 30.33),
        ("everything else, 8.84", 8.84),
    ],
    "residue": 22.90,
    "residue_note": "≈8.4 assembled but unnamed, ≈11.1 never assembled",
    "cost_note": (
        "0.015 GB of the sample's own contigs removes 30.33 points of the library; 16 GB of "
        "Kraken2 removes 0.42. Indexes are therefore consulted in order of cost."
    ),
}

# Panel c's rug is one mark per contig in a schematic, not a dataset: positions come from
# a fixed arithmetic sequence so the drawing is byte-identical on every run. It carries no
# number and makes no claim about where contigs actually fall -- that is the web view's job.
RUG = (22, 38, 57, 93, 110, 148, 166, 193, 224, 246)

INK = "#14181d"
KNOWN = "#2f6fd0"
RESIDUE = "#e2a33c"
RESIDUE_TEXT = "#b07000"

W = 1120
BAR_W = 1040.0  # px for 100% of read pairs
BAR_X = 40
BAR_Y = 672
ROW_TOP = 166
ROW_STEP = 62


def _accent(kind: str) -> tuple[str, str]:
    """Fill and stroke for a stage marker. Filled markers are the two stages that produce
    the headline split; the residue stage carries the residue colour so panel a and panel
    b can be read against each other without a legend."""
    if kind == "known":
        return KNOWN, KNOWN
    if kind == "residue":
        return "#c77700", "#c77700"
    return "#ffffff", KNOWN


def check_sources() -> list[str]:
    cited = {s["source"] for s in STAGES} | {LIBRARY_SPLIT["source"]}
    return sorted(c for c in cited if not (REPO / c).is_dir())


def render() -> str:
    out: list[str] = []
    add = out.append

    height = 918
    add(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {height}" '
        f'width="{W}" height="{height}" role="img" aria-labelledby="figtitle figdesc">'
    )
    add('  <title id="figtitle">hapax: the stages of one run, where a library ends up, '
        'and what one row says</title>')
    desc = (
        "Panel a lists the seven stages of a run, each with its input, its output and one "
        "measured quantity: "
        + "; ".join(f"{s['name']}, {s['number']} {s['unit']}" for s in STAGES)
        + ". Panel b is a to-scale bar of 100% of read pairs: "
        + f"{100 - LIBRARY_SPLIT['residue']:.2f}% identified, made up of "
        + ", ".join(label for label, _ in LIBRARY_SPLIT["identified"])
        + f", and {LIBRARY_SPLIT['residue']:.2f}% residue, of which {LIBRARY_SPLIT['residue_note']}"
        + ". Panel c is the one distance axis every row sits on, running from far from "
        "anything published to matched, each row carrying what it is closest to."
    )
    add(f'  <desc id="figdesc">{desc}</desc>')
    add("")
    add("  <defs>")
    add("    <style>")
    add('      text { font-family: "Helvetica Neue", Helvetica, Arial, ui-sans-serif, '
        "system-ui, sans-serif; fill: " + INK + "; }")
    add("      .figh { font-size: 21px; font-weight: 700; }")
    add("      .figs { font-size: 13px; fill: #5d6671; }")
    add("      .pan  { font-size: 13px; font-weight: 700; }")
    add("      .pant { font-size: 13px; font-weight: 700; fill: " + INK + "; }")
    add("      .st   { font-size: 14.5px; font-weight: 700; letter-spacing: 0.2px; }")
    add("      .io   { font-size: 12px; fill: #5d6671; }")
    add("      .num  { font-size: 14.5px; font-weight: 700; font-variant-numeric: tabular-nums; }")
    add("      .unit { font-size: 11.5px; fill: #5d6671; }")
    add("      .ax   { font-size: 11.5px; fill: #5d6671; }")
    add("      .note { font-size: 11.5px; fill: #717a85; }")
    add("      .seg  { font-size: 11.5px; }")
    add("      .idx  { font-size: 11px; font-weight: 700; fill: #97a0aa; "
        "font-variant-numeric: tabular-nums; }")
    add("    </style>")
    add('    <marker id="arr" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" '
        'markerHeight="6" orient="auto">')
    add('      <path d="M0 1 L7 4 L0 7 z" fill="#4a525c"/>')
    add("    </marker>")
    add("  </defs>")
    add("")
    add(f'  <rect x="0.5" y="0.5" width="{W - 1}" height="{height - 1}" fill="#ffffff" '
        'stroke="#dfe4ea"/>')
    add("")
    add('  <text class="figh" x="40" y="48">hapax</text>')
    add('  <text class="figs" x="40" y="72">Which sequence in a sample is already known, '
        "which assembled without a name, and how far the remainder sits from published "
        "sequence, with the cost of every step.</text>")
    add('  <line x1="40" y1="92" x2="1080" y2="92" stroke="#dfe4ea"/>')
    add("")

    # ---- panel a -------------------------------------------------------------------
    add('  <text class="pan" x="40" y="122">a<tspan class="pant" dx="12">One run, stage by '
        "stage</tspan></text>")
    add('  <text class="ax" x="1080" y="122" text-anchor="end">input, output and one '
        "measured quantity per stage</text>")
    spine_end = ROW_TOP + ROW_STEP * (len(STAGES) - 1) + 30
    add(f'  <line x1="54" y1="140" x2="54" y2="{spine_end}" stroke="#e7ebf0"/>')
    add("")
    add("  <g>")
    for i, stage in enumerate(STAGES):
        cy = ROW_TOP + i * ROW_STEP
        fill, stroke = _accent(stage["accent"])
        add(f'    <circle cx="54" cy="{cy}" r="4.5" fill="{fill}" stroke="{stroke}" '
            'stroke-width="2"/>')
        add(f'    <text class="idx" x="40" y="{cy + 4}" text-anchor="end">{i + 1}</text>')
        add(f'    <text class="st" x="72" y="{cy - 3}">{stage["name"]}</text>')
        add(f'    <text class="io" x="72" y="{cy + 15}">{stage["io"]}</text>')
        add(f'    <text class="num" x="1080" y="{cy - 3}" text-anchor="end">{stage["number"]}</text>')
        add(f'    <text class="unit" x="1080" y="{cy + 15}" text-anchor="end">{stage["unit"]}</text>')
        if i < len(STAGES) - 1:
            add(f'    <line x1="72" y1="{cy + 30}" x2="1080" y2="{cy + 30}" stroke="#f0f3f6"/>')
        add("")
    add("  </g>")
    add('  <line x1="40" y1="580" x2="1080" y2="580" stroke="#e7ebf0"/>')
    add('  <text class="note" x="40" y="598">A stage missing an input is named as skipped '
        "and the run is marked incomplete. Nothing falls back to a cheaper index, and no "
        "sample&#39;s parameters are guessed.</text>")
    add("")

    # ---- panel b -------------------------------------------------------------------
    identified = sum(points for _, points in LIBRARY_SPLIT["identified"])
    residue = LIBRARY_SPLIT["residue"]
    if abs(identified + residue - 100.0) > 0.005:
        raise SystemExit(
            f"panel b does not close: {identified:.2f} identified + {residue:.2f} residue"
        )
    px = BAR_W / 100.0

    add('  <text class="pan" x="40" y="636">b<tspan class="pant" dx="12">One library '
        'by removal route</tspan><tspan class="ax" dx="14" font-weight="400">100% of read pairs, '
        "to scale</tspan></text>")
    add("")
    add(f'  <g transform="translate({BAR_X},{BAR_Y})">')
    add(f'    <rect x="0" y="0" width="{identified * px:.1f}" height="24" fill="{KNOWN}"/>')
    add(f'    <rect x="{identified * px:.1f}" y="0" width="{residue * px:.1f}" height="24" '
        f'fill="{RESIDUE}"/>')
    run = 0.0
    for label, points in LIBRARY_SPLIT["identified"][:-1]:
        run += points
        add(f'    <line x1="{run * px:.1f}" y1="0" x2="{run * px:.1f}" y2="24" stroke="#ffffff"/>')
    add(f'    <rect x="0" y="0" width="{BAR_W:.0f}" height="24" fill="none" stroke="#c9d0d8"/>')
    add("")
    add(f'    <text class="ax" x="0" y="-8" font-weight="700" fill="{INK}">{identified:.2f}% '
        "identified</text>")
    add(f'    <text class="ax" x="{BAR_W:.0f}" y="-8" text-anchor="end" font-weight="700" '
        f'fill="{RESIDUE_TEXT}">{residue:.2f}% residue</text>')
    add("")
    run = 0.0
    for label, points in LIBRARY_SPLIT["identified"]:
        add(f'    <text class="seg" x="{run * px + 4:.0f}" y="42" fill="{KNOWN}">{label}</text>')
        run += points
    add(f'    <text class="seg" x="{BAR_W:.0f}" y="62" text-anchor="end" fill="{RESIDUE_TEXT}">'
        f'{LIBRARY_SPLIT["residue_note"]}</text>')
    add(f'    <text class="note" x="0" y="86">{LIBRARY_SPLIT["cost_note"]}</text>')
    add("  </g>")
    add("")

    # ---- panel c -------------------------------------------------------------------
    add('  <text class="pan" x="40" y="788">c<tspan class="pant" dx="12">What a row '
        "reports</tspan></text>")
    add('  <g transform="translate(40,828)">')
    add('    <line x1="0" y1="0" x2="300" y2="0" stroke="#4a525c" marker-end="url(#arr)"/>')
    add('    <g stroke="#4a525c">')
    for a, b in zip(RUG[0::2], RUG[1::2]):
        add(f'      <line x1="{a}" y1="-8" x2="{a}" y2="0"/>'
            f'<line x1="{b}" y1="-8" x2="{b}" y2="0"/>')
    add("    </g>")
    add('    <text class="ax" x="0" y="16">far from published sequence</text>')
    add('    <text class="ax" x="300" y="16" text-anchor="end">matched</text>')
    add('    <text class="note" x="360" y="-14">One mark per contig, no threshold applied. A '
        "row reports the distance to its closest</text>")
    add('    <text class="note" x="360" y="2">published relative and the identity of that '
        "relative; where one of the two is absent,</text>")
    add('    <text class="note" x="360" y="18">the report states which. No novelty score, no '
        "threat score, no per-read call. Search</text>")
    add('    <text class="note" x="360" y="34">establishes prior knowledge of a sequence; '
        "novelty does not follow from a missing hit.</text>")
    add("  </g>")
    add("")
    add('  <line x1="40" y1="880" x2="1080" y2="880" stroke="#eef1f4"/>')
    add('  <text class="note" x="40" y="898">Berlin wastewater SRR23973193, 5 seeds. '
        "Percentages are of read pairs. Every number comes from a committed directory under "
        '<tspan font-weight="700">results/</tspan>, with the script and seed that produced '
        "it.</text>")
    add("</svg>")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="docs/hapax-overview.svg")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the committed file differs from what this script draws")
    ap.add_argument("--allow-missing-sources", action="store_true",
                    help="draw even when a cited results/ directory is absent")
    args = ap.parse_args()

    missing = check_sources()
    if missing and not args.allow_missing_sources:
        print("cited result directories are missing, refusing to draw:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        return 2

    svg = render()
    target = REPO / args.out
    if args.check:
        if not target.exists():
            print(f"{args.out} does not exist", file=sys.stderr)
            return 1
        if target.read_text() != svg:
            print(f"{args.out} differs from this script's output; re-run without --check",
                  file=sys.stderr)
            return 1
        print(f"{args.out} matches")
        return 0

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(svg)
    print(f"wrote {args.out} ({len(svg)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
