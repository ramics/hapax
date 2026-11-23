"""Tests for the paired-search statistics. Each one is a way the result could be faked.

Run with `python3 -m pytest tests/test_pair_search.py`. Nothing here touches an index: what
is under test is the arithmetic of the joint modes and the geometry that places a
concatenated alignment, because those are this experiment's own code rather than a tool's.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from pair_search import (BANDS, conditional_ladder, e_to_p, fisher4,  # noqa: E402
                         pair_modes)


def row(q, s, e, qs=1, qe=100, bits=None):
    return {"qseqid": q, "sseqid": s, "evalue": e, "bitscore": bits if bits else -e,
            "qstart": qs, "qend": qe, "program": "blastn"}


GEOM = {"len1": 100, "len2": 100, "cat_junction": 100, "catn_junction": 100,
        "catn_spacer": 10}


def test_e_to_p_stays_a_probability():
    """An E-value of 3 is not a p-value of 3, and Fisher breaks silently if it is used so."""
    assert e_to_p(1e-9) == pytest.approx(1e-9, rel=1e-6)
    assert 0 < e_to_p(3.0) < 1.0
    assert e_to_p(1000.0) == pytest.approx(1.0)


def test_fisher_is_monotone_and_beats_neither_mate_alone_when_one_is_noise():
    """Two weak mates combine; one strong mate plus noise comes out WORSE than the mate.

    The second half is the trade the conversion table exists to show: a joint statistic
    dilutes a good single mate with a useless partner, so a joint mode is only safe to
    deploy as `min(best_mate, joint)` and the report must say so.
    """
    both_weak = fisher4(e_to_p(1e-4), e_to_p(1e-4))
    assert both_weak < 1e-4
    assert fisher4(e_to_p(1e-5), e_to_p(1e-5)) < both_weak
    diluted = fisher4(e_to_p(1e-9), e_to_p(5.0))
    assert diluted > 1e-9


def test_agree_is_a_subset_of_best_mate_so_it_cannot_add_sensitivity():
    """`agree` reports the WORSE mate, so its statistic is never below `best_mate`'s.

    Its value, if any, is precision. A report that quotes it as a sensitivity gain has
    misread it.
    """
    rows = [row("p/1", "S1", 1e-20), row("p/2", "S1", 1e-3)]
    m = pair_modes(rows, [], [], GEOM)
    assert m["best_mate"][0] == 1e-20
    assert m["agree"][0] == 1e-3
    assert m["agree"][0] >= m["best_mate"][0]


def test_agree_refuses_when_the_mates_name_different_subjects():
    rows = [row("p/1", "S1", 1e-20), row("p/2", "S2", 1e-20)]
    m = pair_modes(rows, [], [], GEOM)
    assert m["agree"][0] == float("inf")


def test_fisher_needs_the_same_subject_from_both_mates():
    """One mate hitting a subject is not a joint observation, however good it is."""
    rows = [row("p/1", "S1", 1e-30), row("p/2", "S2", 1e-30)]
    m = pair_modes(rows, [], [], GEOM)
    assert m["fisher"][0] == float("inf")
    assert m["_census"]["subjects_both_mates"] == 0


def test_junction_census_places_every_hsp():
    """An HSP inside mate 1, inside mate 2, and one spanning the join, each counted once."""
    cat = [row("p/cat", "S1", 1e-9, 1, 90), row("p/cat", "S1", 1e-9, 110, 190),
           row("p/cat", "S1", 1e-9, 60, 150)]
    m = pair_modes([], cat, [], GEOM)
    j = m["_census"]["junction"]
    assert j["concat:mate1"] == 1
    assert j["concat:mate2"] == 1
    assert j["concat:crosses"] == 1


def test_band_is_past_the_overlap_limit_for_every_read_length_in_the_library():
    """`near` must not need a per-pair nudge: a nudge makes the band depend on the read."""
    lo, _hi = BANDS["near"]
    assert lo > 2 * 150


def test_conditional_ladder_spreads_inside_the_band_and_names_its_mode():
    """The first pass put every pair at one insert length because the ladder had one knot."""
    ladder, note = conditional_ladder(*BANDS["near"])
    lens = {v for _u, v in ladder}
    assert len(lens) > 5
    assert min(lens) >= BANDS["near"][0] and max(lens) <= BANDS["near"][1]
    assert note["mode"] == "conditioned"
    assert 0.05 < note["band_mass_of_library"] < 0.30

    _l, wide = conditional_ladder(*BANDS["wide"])
    assert wide["mode"] == "counterfactual_uniform"
    assert wide["band_mass_of_library"] == 0.0
