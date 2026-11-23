#!/usr/bin/env python3
"""Draw the README overview figure: what a run does and what it reports, on any library.

The figure carries no sample's results. An earlier version drew the Berlin cascade to scale
and put seven measured quantities in the margin, which made the front door a result from one
library rather than a description of the tool (Imo, 2026-11-23: "these are SPECIFIC SAMPLE
INFO. please make it GENERIC. THIS IS THE FRONT DOOR."). The split between identified,
unnameable and unplaced is a property of a library and is measured per run; a reader arriving
at the repo cannot tell from a Berlin bar what their own water would do, and a page shaped
like one sample cannot be compared with another (SPEC 13.4).

So the content here is the pipeline, the three output classes, and the refusals attached to
each. All of it is prose about the tool, which is why this script is a generator rather than
a drawing: the text lives in STAGES and CLASSES below, where it can be read against SPEC.md
in one screen, and `tests/test_overview_figure.py` fails if the committed SVG drifts from it.
Nothing in it is a number read from a results directory, and the test asserts that too.

Usage:
    python3 experiments/make_overview_figure.py --out docs/hapax-overview.svg
    python3 experiments/make_overview_figure.py --check   # verify committed == generated
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------------------
# The stages of `hapax run`, in order, with what each one is given and what it writes.
# The names are the driver's own (SPEC 13.1). Two lines cover several of its stages:
# `score the distance` covers divergent, shortlist and reads; `report the residue` covers
# census and layers. `spec` is the section a reader should check the line against.
# --------------------------------------------------------------------------------------
STAGES = [
    {
        "name": "resolve",
        "io": "library accession → its trim, block size and panels",
        "writes": "every value with how it was established, or a refusal to run",
        "spec": "4.1",
    },
    {
        "name": "merge",
        "io": "read pairs → one read per overlapping molecule",
        "writes": "a merged read where a pair overlaps, per-mate rows where it does not",
        "spec": "4.2",
    },
    {
        "name": "assemble",
        "io": "the whole library → contigs",
        "writes": "contigs with mapped depth and evenness, and recovery against a truth panel",
        "spec": "4.3",
    },
    {
        "name": "identify",
        "io": "reads and contigs against indexes, cheapest first → tiers",
        "writes": "the tier table and the evidence behind every removal, per route",
        "spec": "4.4",
    },
    {
        "name": "score the distance",
        "io": "every contig, removed or kept → one row on one scale",
        "writes": "identity and coverage per gene, nucleotide and protein, with the index used",
        "spec": "4.5, 5",
    },
    {
        "name": "report the residue",
        "io": "the reads left after identification → one row per pair",
        "writes": "every row kept, each field beside its shuffled and known-positive controls",
        "spec": "4.6, 7",
    },
    {
        "name": "subtract other plants",
        "io": "a panel of libraries from elsewhere → witnessed or not",
        "writes": "the witnessed fraction at both operating points, paired by plant",
        "spec": "7.1",
    },
]

# The three output classes of SPEC 5, 6 and 7, in increasing order of distance
# from published sequence, which is the order the scale in panel c runs in
# (Imo, 2026-11-23: "swap A and B, so everything is increasing order of
# divergence"). `note` is drawn only where a class carries one; B has one
# because the gene-level ranking is where this project's novelty work sits. `qualifier` is the sentence that travels with
# the class and keeps a reader from reading more into it than was measured. Each one is
# written as what the row states, because the first version wrote them as refusals ("no
# species call") and a page of negatives says less about the tool than the same facts
# stated forwards (Imo, 2026-11-23: "please remove ALL negatives").
CLASSES = [
    {
        "key": "A",
        "name": "novel strains of known viruses",
        "reports": "the near end of the scale, treated as a first-class case",
        "qualifier": "lineages inside the assembler's fusion limit are reported as fused",
        "note": "",
        "accent": "known",
    },
    {
        "key": "B",
        "name": "divergent but related contigs",
        "reports": "distance from published sequence, ranked on the conserved gene, with "
                   "the closest relative named",
        "qualifier": "the per-gene spread travels with the rank, and a shortlist is cut on "
                     "coverage",
        "note": "One contig's genes sit at different distances, and the gene that is "
                "conserved across a family is the one a novel lineage shows up on.",
        "accent": "known",
    },
    {
        "key": "C",
        "name": "the unplaceable residue",
        "reports": "distance measured, identity still open, with the depth it was seen at",
        "qualifier": "carries the index it was measured against and the depth it was seen at",
        "note": "",
        "accent": "residue",
    },
]


# The pill row. It used to list what the tool never emits; stated forwards, the same four
# facts say what every row does carry, which is the more useful sentence on a front page.
# SPEC 1 and 7.2 keep the refusals themselves, in full.
CARRIES = [
    "a distance",
    "the closest relative, with coverage and E-value",
    "the operating point behind it",
    "its controls",
]

# Panel c's rug is one mark per contig in a schematic. Positions are fixed so the drawing is
# byte-identical on every run, and they make no claim about where contigs fall on any library.
RUG = (22, 38, 57, 93, 110, 148, 166, 193, 224, 246)

INK = "#14181d"
KNOWN = "#2f6fd0"
RESIDUE = "#c77700"

W = 1120
# Layout, derived rather than typed: the canvas has to grow when a stage or a class is added
# to the tables above, and the first version of this figure was clipped at the bottom because
# its height was a literal.
ROW_TOP, ROW_STEP = 162, 58
NOTE_Y = ROW_TOP + ROW_STEP * (len(STAGES) - 1) + 50
B_TOP = NOTE_Y + 58
CLASS_STEP = 62
NOTE_EXTRA = 20  # a class that carries a note needs a line under its row
C_TOP = (B_TOP + 36 + CLASS_STEP * len(CLASSES)
         + NOTE_EXTRA * sum(1 for c in CLASSES if c["note"]) + 24)
REFUSAL_Y = C_TOP + 100
H = REFUSAL_Y + 92


def _accent(kind: str) -> str:
    return KNOWN if kind == "known" else RESIDUE


def render() -> str:
    out: list[str] = []
    add = out.append

    add(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" height="{H}" role="img" aria-labelledby="figtitle figdesc">'
    )
    add('  <title id="figtitle">hapax: the stages of a run, the three output classes, and '
        "the scale they share</title>")
    desc = (
        "Panel a lists the seven stages of a run, with what each is given and what it writes: "
        + "; ".join(f"{s['name']}, {s['io']}, writes {s['writes']}" for s in STAGES)
        + ". Panel b lists the three output classes, what each reports and the qualifier that "
        "travels with it: "
        + "; ".join(f"{c['key']}, {c['name']}, reports {c['reports']}, {c['qualifier']}"
                    for c in CLASSES)
        + ". Panel c is the distance scale all three sit on, running from matched to far "
        "from published sequence, and the row beneath it lists what every output carries: "
        + ", ".join(CARRIES) + "."
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
    add("      .wr   { font-size: 12px; fill: #3c444e; }")
    add("      .ax   { font-size: 11.5px; fill: #5d6671; }")
    add("      .note { font-size: 11.5px; fill: #717a85; }")
    add("      .idx  { font-size: 11px; font-weight: 700; fill: #97a0aa; }")
    add("      .cls  { font-size: 14px; font-weight: 700; }")
    add("      .ref  { font-size: 11.5px; fill: #8a5a12; }")
    add("    </style>")
    add('    <marker id="arr" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" '
        'markerHeight="6" orient="auto">')
    add('      <path d="M0 1 L7 4 L0 7 z" fill="#4a525c"/>')
    add("    </marker>")
    add("  </defs>")
    add("")
    add(f'  <rect x="0.5" y="0.5" width="{W - 1}" height="{H - 1}" fill="#ffffff" '
        'stroke="#dfe4ea"/>')
    add("")
    add('  <text class="figh" x="40" y="48">hapax</text>')
    add('  <text class="figs" x="40" y="72">Known sequence is identified and removed. The '
        "remainder is measured against published sequence.</text>")
    add('  <line x1="40" y1="92" x2="1080" y2="92" stroke="#dfe4ea"/>')
    add("")

    # ---- panel a: the stages -------------------------------------------------------
    add('  <text class="pan" x="40" y="122">a<tspan class="pant" dx="12">A run, stage by '
        "stage</tspan></text>")
    add('  <text class="ax" x="1080" y="122" text-anchor="end">given → produced, and what '
        "the run writes</text>")
    top, step = ROW_TOP, ROW_STEP
    add(f'  <line x1="54" y1="138" x2="54" y2="{top + step * (len(STAGES) - 1) + 26}" '
        'stroke="#e7ebf0"/>')
    add("")
    add("  <g>")
    for i, stage in enumerate(STAGES):
        cy = top + i * step
        add(f'    <circle cx="54" cy="{cy}" r="4.5" fill="#ffffff" stroke="{KNOWN}" '
            'stroke-width="2"/>')
        add(f'    <text class="idx" x="40" y="{cy + 4}" text-anchor="end">{i + 1}</text>')
        add(f'    <text class="st" x="72" y="{cy - 2}">{stage["name"]}</text>')
        add(f'    <text class="io" x="72" y="{cy + 16}">{stage["io"]}</text>')
        add(f'    <text class="wr" x="1080" y="{cy - 2}" text-anchor="end">{stage["writes"]}'
            "</text>")
        add(f'    <text class="ax" x="1080" y="{cy + 16}" text-anchor="end">SPEC &#167;'
            f'{stage["spec"]}</text>')
        if i < len(STAGES) - 1:
            add(f'    <line x1="72" y1="{cy + 28}" x2="1080" y2="{cy + 28}" stroke="#f0f3f6"/>')
        add("")
    add("  </g>")
    note_y = NOTE_Y
    add(f'  <line x1="40" y1="{note_y - 18}" x2="1080" y2="{note_y - 18}" stroke="#e7ebf0"/>')
    add(f'  <text class="note" x="40" y="{note_y}">Removal applies to read pairs, so every '
        "contig is measured, named or unnamed. How much a tier removes is a property of the "
        "library and is measured per run.</text>")
    add("")

    # ---- panel b: the output classes ------------------------------------------------
    b_top = B_TOP
    add(f'  <text class="pan" x="40" y="{b_top}">b<tspan class="pant" dx="12">Three output '
        "classes, on one scale</tspan></text>")
    add(f'  <text class="ax" x="1080" y="{b_top}" text-anchor="end">what each reports, and '
        "the qualifier it carries</text>")
    add("")
    add("  <g>")
    cy = b_top + 36
    for cls in CLASSES:
        colour = _accent(cls["accent"])
        add(f'    <rect x="40" y="{cy - 14}" width="4" height="32" fill="{colour}"/>')
        add(f'    <text class="cls" x="60" y="{cy - 1}" fill="{colour}">{cls["key"]}'
            f'<tspan class="cls" dx="10" fill="{INK}">{cls["name"]}</tspan></text>')
        add(f'    <text class="io" x="60" y="{cy + 17}">{cls["reports"]}</text>')
        add(f'    <text class="ref" x="1080" y="{cy - 1}" text-anchor="end">'
            f'{cls["qualifier"]}</text>')
        if cls["note"]:
            add(f'    <text class="note" x="60" y="{cy + 34}">{cls["note"]}</text>')
        cy += CLASS_STEP + (NOTE_EXTRA if cls["note"] else 0)
    add("  </g>")
    add("")

    # ---- panel c: the scale ---------------------------------------------------------
    c_top = C_TOP
    add(f'  <text class="pan" x="40" y="{c_top}">c<tspan class="pant" dx="12">What a row '
        "reports</tspan></text>")
    add(f'  <g transform="translate(40,{c_top + 40})">')
    add('    <line x1="0" y1="0" x2="300" y2="0" stroke="#4a525c" marker-end="url(#arr)"/>')
    add('    <g stroke="#4a525c">')
    for a, b in zip(RUG[0::2], RUG[1::2]):
        add(f'      <line x1="{a}" y1="-8" x2="{a}" y2="0"/>'
            f'<line x1="{b}" y1="-8" x2="{b}" y2="0"/>')
    add("    </g>")
    add('    <text class="ax" x="0" y="16">matched</text>')
    add('    <text class="ax" x="300" y="16" text-anchor="end">far from published '
        "sequence</text>")
    add('    <text class="note" x="360" y="-14">One mark per contig, every contig kept. A row '
        "gives the distance to its closest</text>")
    add('    <text class="note" x="360" y="2">published relative and the identity of that '
        "relative, and states which of the</text>")
    add('    <text class="note" x="360" y="18">two it has. Search establishes prior knowledge, '
        "and score ordering tracks</text>")
    add('    <text class="note" x="360" y="34">proximity to the database, so the report stops '
        "at the distance it measured.</text>")
    add("  </g>")
    add("")

    # ---- what every output carries ---------------------------------------------------
    r_y = REFUSAL_Y
    add(f'  <text class="ax" x="40" y="{r_y + 19}">Every row carries</text>')
    x = 160
    for carried in CARRIES:
        width = 16 + int(len(carried) * 6.6)
        add(f'    <rect x="{x}" y="{r_y}" width="{width}" height="28" rx="14" fill="#fffaf2" '
            'stroke="#e6d3b3"/>')
        add(f'    <text class="ref" x="{x + width / 2:.0f}" y="{r_y + 19}" '
            f'text-anchor="middle">{carried}</text>')
        x += width + 14
    add("")
    add(f'  <line x1="40" y1="{H - 38}" x2="1080" y2="{H - 38}" stroke="#eef1f4"/>')
    add(f'  <text class="note" x="40" y="{H - 20}">Stage by stage in <tspan '
        'font-weight="700">SPEC.md</tspan>, with the measurement behind each one in <tspan '
        'font-weight="700">results/</tspan>.</text>')
    add("</svg>")
    return "\n".join(out) + "\n"


def numbers_in(svg: str) -> list[str]:
    """Percentages and counts in the drawing's text, which the figure must not carry.

    Section references (SPEC 4.1) and the geometry are not results; a percentage or a
    quantity with units is. This is the check that keeps the front door generic as the
    figure is edited, since the easy mistake is to answer a reviewer with a number.
    """
    text = " ".join(re.findall(r">([^<>]+)<", svg))
    return re.findall(r"\b\d[\d,.]*\s*(?:%|GB|nt|points?|reads?|pairs?|libraries)\b", text)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="docs/hapax-overview.svg")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the committed file differs from what this script draws")
    args = ap.parse_args()

    svg = render()
    stray = numbers_in(svg)
    if stray:
        print("the figure carries sample numbers, refusing to draw:", stray, file=sys.stderr)
        return 2

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
