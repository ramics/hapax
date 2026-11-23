#!/usr/bin/env python3
"""Draw the README figure of what a run writes, by unit.

This picture replaced a paragraph. The README described the output as running prose over
four units at once -- library fractions, per-read fields, per-subject counts and the
divergence ladder -- and a reader could not tell which field belonged to which unit
(Imo, 2026-11-23: "so confusing. replace with a picture?"). The units are the natural
rows of a table, so the table is the figure.

Same rules as the overview figure beside it: the text lives in OUTPUTS, where it can be
read against SPEC.md, the committed SVG is generated and never hand-edited, and the script
refuses to draw a percentage or a count, because a fraction is a property of whichever
library produced it and this is the front door (SPEC 13.4).

Usage:
    python3 experiments/make_outputs_figure.py --out docs/hapax-outputs.svg
    python3 experiments/make_outputs_figure.py --check
"""

from __future__ import annotations

import argparse
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

from make_overview_figure import INK, KNOWN, RESIDUE, numbers_in  # noqa: E402

# One row per unit a run writes at. `row` is what one row of that table is, `carries` the
# fields on it, `spec` the section to check the line against. Order is coarsest unit first,
# which is also the order a reader opens them in.
OUTPUTS = [
    {
        "unit": "library",
        "row": "one breakdown per run: identified, assembled but unnameable, unplaced",
        "carries": "each fraction with its interval and the operating point that produced it",
        "spec": "4.4, 7.1",
        "accent": "known",
    },
    {
        "unit": "contig",
        "row": "one row per contig, whether or not an index named it",
        "carries": "distance to the closest published relative, identity and coverage per "
                   "gene, nucleotide and protein, with depth and evenness",
        "spec": "5.1",
        "accent": "known",
    },
    {
        "unit": "subject",
        "row": "one row per published relative the residue points at",
        "carries": "how many reads point at it and from how far, which is where a novel "
                   "strain of a known virus appears",
        "spec": "6",
        "accent": "known",
    },
    {
        "unit": "read pair",
        "row": "one row per surviving pair, un-thresholded",
        "carries": "distance with coverage and E-value, kind of sequence, recurrence in the "
                   "library, and what it sits next to on the molecule",
        "spec": "7.1",
        "accent": "residue",
    },
    {
        "unit": "divergence ladder",
        "row": "one rung per divergence, where a library carries a genome of known distance",
        "carries": "what the cascade recovered at each rung, on an ordinal scale",
        "spec": "8",
        "accent": "residue",
    },
]

W = 1120
ROW_TOP, ROW_STEP = 148, 72
NOTE_Y = ROW_TOP + ROW_STEP * (len(OUTPUTS) - 1) + 56
H = NOTE_Y + 46


def _wrap(text: str, width: int) -> list[str]:
    """Greedy wrap in characters, because an SVG text element does not wrap itself.

    The first version of this table right-aligned two columns against each other and the
    long cells overlapped in the middle of the figure. Columns are left-aligned and wrapped
    on a fixed width so the drawing is byte-identical on every run.
    """
    lines: list[str] = []
    line = ""
    for word in text.split():
        candidate = f"{line} {word}".strip()
        if len(candidate) > width and line:
            lines.append(line)
            line = word
        else:
            line = candidate
    if line:
        lines.append(line)
    return lines


def render() -> str:
    out: list[str] = []
    add = out.append

    add(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" height="{H}" role="img" aria-labelledby="outtitle outdesc">'
    )
    add('  <title id="outtitle">What a hapax run writes, by unit</title>')
    desc = (
        "A table of what a run writes, one row per unit: "
        + "; ".join(f'{o["unit"]}, {o["row"]}, carrying {o["carries"]} (SPEC {o["spec"]})'
                    for o in OUTPUTS)
        + ". Every field is reported beside its shuffled null and a known-positive "
        "comparison, and a verdict is made at contig and subject level."
    )
    add(f'  <desc id="outdesc">{desc}</desc>')
    add("")
    add("  <defs>")
    add("    <style>")
    add('      text { font-family: "Helvetica Neue", Helvetica, Arial, ui-sans-serif, '
        "system-ui, sans-serif; fill: " + INK + "; }")
    add("      .figh { font-size: 21px; font-weight: 700; }")
    add("      .figs { font-size: 13px; fill: #5d6671; }")
    add("      .un   { font-size: 14.5px; font-weight: 700; }")
    add("      .io   { font-size: 12px; fill: #5d6671; }")
    add("      .wr   { font-size: 12px; fill: #3c444e; }")
    add("      .ax   { font-size: 11.5px; fill: #5d6671; }")
    add("      .note { font-size: 11.5px; fill: #717a85; }")
    add("    </style>")
    add("  </defs>")
    add("")
    add(f'  <rect x="0.5" y="0.5" width="{W - 1}" height="{H - 1}" fill="#ffffff" '
        'stroke="#dfe4ea"/>')
    add("")
    add('  <text class="figh" x="40" y="48">What a run writes</text>')
    add('  <text class="figs" x="40" y="72">Every field is reported beside its shuffled '
        "null and a known-positive comparison.</text>")
    add('  <line x1="40" y1="92" x2="1080" y2="92" stroke="#dfe4ea"/>')
    add('  <text class="ax" x="40" y="118">unit</text>')
    add('  <text class="ax" x="300" y="118">one row is</text>')
    add('  <text class="ax" x="660" y="118">what the row carries</text>')
    add("")

    add("  <g>")
    for i, o in enumerate(OUTPUTS):
        cy = ROW_TOP + i * ROW_STEP
        colour = KNOWN if o["accent"] == "known" else RESIDUE
        add(f'    <rect x="40" y="{cy - 14}" width="4" height="46" fill="{colour}"/>')
        add(f'    <text class="un" x="60" y="{cy + 2}">{o["unit"]}</text>')
        add(f'    <text class="ax" x="60" y="{cy + 22}">SPEC &#167;{o["spec"]}</text>')
        for j, line in enumerate(_wrap(o["row"], 52)):
            add(f'    <text class="io" x="300" y="{cy + 2 + 18 * j}">{line}</text>')
        for j, line in enumerate(_wrap(o["carries"], 66)):
            add(f'    <text class="wr" x="660" y="{cy + 2 + 18 * j}">{line}</text>')
        if i < len(OUTPUTS) - 1:
            add(f'    <line x1="40" y1="{cy + 44}" x2="1080" y2="{cy + 44}" '
                'stroke="#f0f3f6"/>')
        add("")
    add("  </g>")
    add(f'  <line x1="40" y1="{NOTE_Y - 20}" x2="1080" y2="{NOTE_Y - 20}" stroke="#e7ebf0"/>')
    add(f'  <text class="note" x="40" y="{NOTE_Y}">A verdict is made at contig and subject '
        "level, where there is enough sequence to support one. Per-read fields are reported "
        "and left un-thresholded.</text>")
    add("</svg>")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="docs/hapax-outputs.svg")
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
