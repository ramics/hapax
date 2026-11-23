"""Tests for the per-library facts, one per way this has gone wrong or could.

Run with `python3 -m pytest tests/test_library_facts.py`. Nothing here touches a FASTQ:
what is under test is whether a value is resolved, kept or refused, because the failure
mode this file exists for is silent - a 14 nt trim on a library with no non-template
prefix deletes real bases from every mate, in the direction of less recognition, and
nothing in the output says so.

The two assertions that matter most are the negatives. A library with no entry must
raise, and a read path naming a DIFFERENT run must raise even when --library is left at
its default, because that is the case the ~20 scripts carrying `default=14` actually hit.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import library_facts as L  # noqa: E402


def parser(**extra):
    ap = argparse.ArgumentParser()
    L.add_argument(ap)
    ap.add_argument("--r1", type=Path, default=None)
    ap.add_argument("--r2", type=Path, default=None)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--block", type=int, default=None)
    for flag, kw in extra.items():
        ap.add_argument(flag, **kw)
    return ap


BERLIN = ["--r1", "data/berlin/SRR23973193_1.fastq.gz",
          "--r2", "data/berlin/SRR23973193_2.fastq.gz"]


# ------------------------------------------------------------------ the committed values

def test_the_file_parses_and_every_value_names_what_established_it():
    d = json.loads((ROOT / "experiments/libraries.json").read_text())
    for acc, lib in d["libraries"].items():
        for key in ("trim_r2", "block"):
            if key in lib:
                assert f"{key}_basis" in lib, f"{acc}: {key} with no basis"
                assert len(lib[f"{key}_basis"]) > 40, f"{acc}: {key}_basis is not a basis"


def test_berlin_resolves_to_the_values_every_committed_run_used():
    """The point of the whole exercise is that nothing moved. 14 and 1.5 M, as before."""
    args = parser().parse_args(BERLIN)
    L.resolve_args(args, ["trim_r2", "block"])
    assert (args.trim_r2, args.block) == (14, 1_500_000)


def test_the_holdout_library_is_not_trimmed():
    """It is single-end. A 14 nt trim here would delete real sequence from every read."""
    args = parser().parse_args(["--library", "SRR23421544",
                                "--r1", "reads/SRR23421544.fastq.gz"])
    L.resolve_args(args, ["trim_r2", "block"])
    assert (args.trim_r2, args.block) == (0, 1_000_000)


def test_the_map_back_trim_is_its_own_operating_point():
    """divergent_known trimmed 16 where the measured prefix is 14. Pinned, not reconciled."""
    args = parser().parse_args(BERLIN)
    L.resolve_args(args, {"trim_r2": "trim_r2_for_mapback"})
    assert args.trim_r2 == 16


# ---------------------------------------------------------------------- the refusals

def test_an_unknown_library_raises_rather_than_inheriting_berlins():
    args = parser().parse_args(["--library", "SRR99999999"])
    with pytest.raises(L.UnknownLibrary) as e:
        L.resolve_args(args, ["trim_r2"])
    assert "libraries.json" in str(e.value)


def test_a_read_path_naming_another_run_is_refused_at_the_default_library():
    """The hole the mechanical fix would otherwise leave: --library defaults to Berlin, so
    pointing --r1 at another library's FASTQ and saying nothing used to be silent."""
    args = parser().parse_args(["--r1", "reads/SRR23421544.fastq.gz"])
    with pytest.raises(L.UnknownLibrary) as e:
        L.resolve_args(args, ["trim_r2"])
    assert "SRR23421544" in str(e.value)


def test_an_underscore_does_not_hide_the_accession():
    """\\b would not match `ww_SRR23421544`, which is exactly the shape a staged copy has."""
    assert not L.paths_name_library("SRR23973193", [Path("data/ww_SRR23421544/r1.fq.gz")])


def test_a_path_naming_no_run_cannot_disagree():
    """`data/berlin` and a work directory are the library's reads under another name;
    refusing those would refuse every committed command in the repo."""
    assert L.paths_name_library("SRR23973193", [Path("data/berlin"), Path("/mnt/work")])


def test_a_missing_key_raises_by_name():
    with pytest.raises(L.UnknownLibrary) as e:
        L.fact("SRR23421544", "r1")
    assert "SRR23421544" in str(e.value) and "r1" in str(e.value)


# ------------------------------------------------------------- stating it still overrides

def test_an_explicit_value_is_kept_and_never_checked():
    args = parser().parse_args(["--r1", "reads/SRR23421544.fastq.gz", "--trim-r2", "7"])
    L.resolve_args(args, ["trim_r2"])
    assert args.trim_r2 == 7


def test_zero_is_a_value_and_not_an_absence():
    """The holdout's trim is 0. A truth test here would read it as 'not given' and go
    back to the library, which is the same bug in a different coat."""
    args = parser().parse_args(["--library", "SRR23421544", "--trim-r2", "0",
                                "--r1", "reads/SRR23421544.fastq.gz"])
    L.resolve_args(args, ["trim_r2"])
    assert args.trim_r2 == 0


def test_resolve_refuses_when_nothing_names_a_library():
    with pytest.raises(L.UnknownLibrary):
        L.resolve(None, "trim_r2", None)


# ------------------------------------------------- a read path is a fact, not a validator
#
# It was a validator first, which is why three scripts still carried their own defaults,
# all three naming `data/ww_berlin/` where this file and data/MANIFEST.md say
# `data/berlin/` (results/PROVENANCE.md §5). A path that names one library's reads is the
# same shape of defect as a default that names its trim.

def test_a_read_path_resolves_from_the_library():
    args = parser().parse_args(["--library", "SRR23973193"])
    L.resolve_args(args, ["r1", "r2", "trim_r2"])
    assert args.r1 == ROOT / "data/berlin/SRR23973193_1.fastq.gz"
    assert args.r2 == ROOT / "data/berlin/SRR23973193_2.fastq.gz"
    assert args.trim_r2 == 14


def test_a_resolved_read_path_is_absolute_so_the_cwd_cannot_change_it():
    """Staged copies live under scratch directories and EC2 mounts; a relative default
    read from the wrong cwd is a missing file at best and another library's reads at
    worst."""
    args = parser().parse_args(["--library", "SRR23973193"])
    L.resolve_args(args, ["r1"])
    assert args.r1.is_absolute()


def test_a_stated_read_path_still_wins():
    args = parser().parse_args(["--library", "SRR23973193",
                                "--r1", "scratch/SRR23973193_1.fastq.gz"])
    L.resolve_args(args, ["r1", "r2"])
    assert args.r1 == Path("scratch/SRR23973193_1.fastq.gz")
    assert args.r2 == ROOT / "data/berlin/SRR23973193_2.fastq.gz"


def test_one_mate_naming_another_run_refuses_the_other_mate():
    """Agreement is decided before any path is filled. Giving mate 1 of one run and
    letting mate 2 default must not hand over the named library's mate 2."""
    args = parser().parse_args(["--library", "SRR23973193",
                                "--r1", "data/holdout/SRR23421544_1.fastq.gz"])
    with pytest.raises(L.UnknownLibrary) as e:
        L.resolve_args(args, ["r1", "r2"])
    assert "SRR23421544" in str(e.value)


def test_a_library_with_no_read_path_refuses_rather_than_borrowing_one():
    """The holdout library's entry has no r1: nothing in this repo has staged it under a
    recorded path, and borrowing Berlin's is the failure this file exists to stop."""
    args = parser().parse_args(["--library", "SRR23421544"])
    with pytest.raises(L.UnknownLibrary) as e:
        L.resolve_args(args, ["r1"])
    assert "r1" in str(e.value)


def test_the_three_scripts_carry_no_read_path_default_any_more():
    """Name-based and deliberately crude: the defect was a literal, so the test reads for
    literals.

    SCOPED TO THE READ PATHS, and that scope is a finding rather than convenience. The
    reads are documented under `data/berlin/` - 62 committed manifests against 18 for
    `data/ww_berlin/` - but the ASSEMBLY subsets are documented the other way: every
    committed manifest that names one says `data/ww_berlin/spades_sample*.fa`, 11 to 13
    occurrences and none under `data/berlin/`. So the directory was not simply renamed,
    and `absence_residual.py --ref` and `rrna_contigs.py --contigs` are NOT stale: they
    name the spelling their own results were produced under. Correcting them to match the
    reads would point them at a path nothing documents, which is the opposite of the fix.
    They are still sample-specific defaults and belong in libraries.json as an `assembly`
    fact; that needs somebody who can see `data/` on disk, since the one thing this
    container cannot check is which spelling exists (results/PROVENANCE.md §5).
    """
    for name in ("absence_curve.py", "absence_residual.py", "rrna_match_source.py"):
        text = (ROOT / "experiments" / name).read_text()
        for line in text.splitlines():
            if '"--r1"' in line or '"--r2"' in line:
                assert "default=None" in line, f"{name}: {line.strip()}"
                assert "ww_berlin" not in line, f"{name}: {line.strip()}"
