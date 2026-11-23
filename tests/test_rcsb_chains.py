"""Tripwires on the Foldseek-target parser in `rcsb_chains`.

Every case here is one the first run actually produced and the old entry-level code got
wrong without erroring, which is why they are pinned rather than assumed. The network half
of the module is not exercised: these are the pure-function traps.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from rcsb_chains import split_target  # noqa: E402


def test_the_ordinary_shape_keeps_its_chain():
    assert split_target("8GH6_A") == ("8GH6", "A")


def test_a_lower_case_or_multi_character_auth_chain_survives():
    # 10PX_1i and 7Z4O_EEE are real subjects from the first run. Truncating or upper-casing
    # the chain would resolve each to a different molecule of the same entry.
    assert split_target("10PX_1i") == ("10PX", "1i")
    assert split_target("7Z4O_EEE") == ("7Z4O", "EEE")


def test_an_nmr_ensemble_puts_the_chain_last():
    # Splitting on the first underscore gave 40 subjects in the first run a chain of
    # "MODEL", which resolves to nothing at all.
    assert split_target("2MTV_MODEL_10_A") == ("2MTV", "A")
    assert split_target("2MTV_MODEL_1_B") == ("2MTV", "B")


def test_a_chainless_target_reports_no_chain_rather_than_guessing():
    # 78 subjects in the first run carried no chain suffix. The caller resolves these only
    # where the entry holds exactly one protein entity; it must not be handed a fake chain.
    assert split_target("13XK") == ("13XK", "")


def test_a_cif_suffix_is_not_mistaken_for_a_chain():
    assert split_target("8GH6.cif.gz_A") == ("8GH6", "A")
