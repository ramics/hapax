"""Tests for the README overview figure, which is generated and not hand-drawn.

Two failure modes are under test, and both have happened in this repo to other figures.
The first is prose drift: the committed SVG is edited by hand, the generator's table is
not, and the picture quietly stops agreeing with the record it cites. The second is a
figure that does not close -- panel b is drawn to scale from percentage points, so a
value edited in the table without its neighbours produces a bar that no longer sums to
the library and nothing in the drawing says so.

The third assertion is the cheap one that keeps the citations honest: every quantity in
the table names a results/ directory, and that directory has to exist in the tree.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import make_overview_figure as fig  # noqa: E402


def test_committed_svg_is_what_the_script_draws():
    committed = (ROOT / "docs" / "hapax-overview.svg").read_text()
    assert committed == fig.render(), (
        "docs/hapax-overview.svg was hand-edited; re-run "
        "python3 experiments/make_overview_figure.py"
    )


def test_panel_b_closes_to_the_whole_library():
    identified = sum(points for _, points in fig.LIBRARY_SPLIT["identified"])
    assert abs(identified + fig.LIBRARY_SPLIT["residue"] - 100.0) < 0.005


def test_every_quantity_cites_a_directory_that_exists():
    assert fig.check_sources() == []
    for stage in fig.STAGES:
        assert stage["source"].startswith("results/")
        assert stage["number"] and stage["unit"]


def test_the_figure_makes_no_verdict_claim():
    """The refusals are part of the deliverable (SPEC 1, 7.2), so a figure that drops
    them is a regression even though nothing else would notice."""
    svg = fig.render()
    for refusal in ("No novelty score.", "No threat score.", "No per-read call."):
        assert refusal in svg
    assert "It does not establish novelty." in svg
    for forbidden in ("threat level", "novelty score of", "likely novel", "probably new"):
        assert forbidden not in svg
