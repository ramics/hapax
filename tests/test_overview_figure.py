"""Tests for the README overview figure, which is generated and not hand-drawn.

Three failure modes, all of which have happened here. The first is prose drift: the committed
SVG is edited by hand, the generator's tables are not, and the picture stops agreeing with
the spec it cites. The second is a sample leaking onto the front door -- the first version of
this figure drew one library's cascade to scale, which made the repo's first page a result
rather than a description of the tool. The third is a refusal quietly falling off the figure
during a wording pass, which nothing else would notice.
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


def test_the_figure_carries_no_sample_numbers():
    """The front door describes the report. A percentage or a count on it would be a
    property of whichever library produced it, and a reader cannot tell from one library's
    bar what their own water would do (SPEC 13.4)."""
    assert fig.numbers_in(fig.render()) == []


def test_every_stage_states_input_output_and_what_it_writes():
    for stage in fig.STAGES:
        assert "→" in stage["io"], stage["name"]
        assert stage["writes"]
        assert stage["spec"]


def test_classes_run_from_the_near_end_outwards():
    """Increasing distance from published sequence, the same direction panel c's scale
    runs in (Imo, 2026-11-23). The near end is first and carries the fusion limit."""
    assert [c["name"] for c in fig.CLASSES][0] == "novel strains of known viruses"
    assert [c["name"] for c in fig.CLASSES][-1] == "the unplaceable residue"
    svg = fig.render()
    assert '<text class="ax" x="0" y="16">matched</text>' in svg


def test_the_gene_level_ranking_is_on_the_divergent_class():
    """What class B ranks on is the project's own novelty result (SPEC 2.2), so the figure
    states it rather than leaving the class as a distance alone."""
    b = next(c for c in fig.CLASSES if c["key"] == "B")
    assert "conserved gene" in b["reports"]
    assert b["note"] and "different distances" in b["note"]


def test_every_output_class_carries_its_qualifier():
    """What keeps a class honest travels with it. The qualifiers are written forwards, as
    what the row states, so the test asserts they are present rather than looking for a
    "no ..." string: the front door says what the report contains (SPEC 1, 7.2 keep the
    refusals in full)."""
    assert [c["key"] for c in fig.CLASSES] == ["A", "B", "C"]
    for cls in fig.CLASSES:
        assert cls["reports"] and cls["qualifier"]


def test_the_figure_claims_no_verdict_and_states_the_limit():
    svg = fig.render()
    assert "score ordering tracks" in svg
    assert "the report stops at the distance it measured" in svg
    for word in ("novelty score", "threat score", "risk label"):
        assert word not in svg, f"{word} belongs in SPEC, not on the figure"
