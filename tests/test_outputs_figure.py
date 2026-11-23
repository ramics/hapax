"""Tests for the figure of what a run writes, which is generated and not hand-drawn.

It replaced a paragraph that mixed four units of output into one block of prose, so the
thing worth asserting is that every row still names its unit, what one row of that unit is,
the fields it carries and the section to check it against. The front-door rules from the
overview figure apply here too: the committed SVG is the script's output, and no sample's
numbers appear on it.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import make_outputs_figure as fig  # noqa: E402
import make_overview_figure as overview  # noqa: E402


def test_committed_svg_is_what_the_script_draws():
    committed = (ROOT / "docs" / "hapax-outputs.svg").read_text()
    assert committed == fig.render(), (
        "docs/hapax-outputs.svg was hand-edited; re-run "
        "python3 experiments/make_outputs_figure.py"
    )


def test_the_figure_carries_no_sample_numbers():
    assert overview.numbers_in(fig.render()) == []


def test_every_unit_states_its_row_fields_and_section():
    units = [o["unit"] for o in fig.OUTPUTS]
    assert units == ["library", "contig", "subject", "read pair", "divergence ladder"]
    for o in fig.OUTPUTS:
        assert o["row"] and o["carries"] and o["spec"]


def test_columns_do_not_overlap():
    """The first version right-aligned two long columns against each other and they
    collided mid-figure. Wrapping is in characters, so the guard is on line length."""
    for o in fig.OUTPUTS:
        assert all(len(line) <= 52 for line in fig._wrap(o["row"], 52))
        assert all(len(line) <= 66 for line in fig._wrap(o["carries"], 66))
