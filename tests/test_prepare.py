"""Measuring, assembling and splitting as stages of the run. PLAN.md step 5.

WHAT THIS STEP REMOVED. `hapax run` started after the hardest three steps: a caller had to
supply `--contigs`, `--unidentified-contigs`, `--identified-contigs`, `--rrna-contigs`,
`--assembly`, `--cross-sample-ref` and a measured `--trim-r2`, and the only thing that knew
how to produce them was 250 lines of shell plus EC2 UserData. The tool therefore could not
be run on a new input by reading the repo (PLAN.md 1.2).

The tests here are about the one value in that list that can do silent damage. A trim is a
property of a sequencing run; written as `default=14` it deletes fourteen real bases from
every mate of a library that has no non-template prefix, in the direction of less
recognition. So the resolution order is tested, and so is the refusal at the end of it:
there is no fourth source.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from hapax.config import ConfigError, load
from hapax.prepare import Prepared, resolve_trim

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

BASE = """
[run]
seed = 1
seeds = 1

[inputs]
label = "{label}"
"""


def _cfg(tmp_path, label="fixture", extra=""):
    p = tmp_path / "run.toml"
    p.write_text(BASE.format(label=label) + extra)
    return load(p)


def test_a_stated_trim_wins_and_carries_where_it_came_from(tmp_path):
    cfg = _cfg(tmp_path, extra='\n[library]\ntrim_r2 = 7\n'
                                'trim_r2_basis = "measured on this kit"\n')
    value, source = resolve_trim(cfg, {"layout": "paired"})
    assert value == 7
    assert "measured on this kit" in source, (
        "the run has to record WHERE the trim came from, not only what it was: a stated "
        "value and a looked-up one are different claims")


def test_a_known_library_is_read_from_the_registry(tmp_path):
    cfg = _cfg(tmp_path, label="SRR23973193")
    value, source = resolve_trim(cfg, {"layout": "paired"})
    assert value == 14
    assert "libraries.json" in source


def test_a_single_end_library_needs_no_measurement(tmp_path):
    """The one case with an answer that is not a measurement: there is no mate 2."""
    cfg = _cfg(tmp_path, extra="")
    value, source = resolve_trim(cfg, {"layout": "single"})
    assert value == 0
    assert "no mate 2" in source


def test_an_unknown_library_with_no_stated_trim_is_refused(tmp_path):
    """The fourth source does not exist, and the refusal says what to do instead."""
    cfg = _cfg(tmp_path)
    with pytest.raises(ConfigError) as e:
        resolve_trim(cfg, {"layout": "paired"})
    message = str(e.value)
    assert "no right one to assume" in message
    assert "trim_r2_basis" in message
    assert "nontemplate_prefix_map.py" in message, (
        "the refusal has to name the instrument that would measure it, or it is a dead "
        "end rather than a next step")


def test_the_registry_is_never_consulted_for_a_library_the_reads_contradict(tmp_path):
    """A stated label is the override and is not checked; an absent one is not guessed.

    This is the hole `library_facts.resolve_args` exists for, one layer out: a caller who
    points the reads at one library and names another gets the named library's trim unless
    something checks. Here the check is that there is no default at all, so an unknown
    label cannot inherit anything.
    """
    cfg = _cfg(tmp_path, label="SRR00000000")
    with pytest.raises(ConfigError):
        resolve_trim(cfg, {"layout": "paired"})


def test_one_seed_has_no_stride_and_that_is_recorded_not_refused():
    """A one-seed instrument must be measurable; a stride only separates two or more.

    `stride_for` raises for a single seed, correctly, because there is nothing for a
    stride to keep apart. Refusing the whole measurement on that basis would mean a
    fixture or a pilot could not be measured at all, so the absence is recorded with its
    reason. A headline result still needs five seeds.
    """
    import measure_library

    with pytest.raises(ValueError, match="two or more seeds"):
        measure_library.stride_for(10_000, 1_000, 1)
    # And five blocks that do not fit are a refusal rather than overlapping blocks, which
    # is the case that produced a directory looking like a five-seed result.
    with pytest.raises(ValueError, match="DISJOINT"):
        measure_library.stride_for(9_987, 9_000, 5)


def test_the_prepared_manifest_says_the_run_produced_its_own_inputs(tmp_path):
    """What a reader of the run directory needs: these paths were not handed in."""
    p = Prepared(trim_r2=0, trim_r2_source="single-end", assembly=tmp_path / "a.fa",
                 identified_norrna=tmp_path / "i.fa", rrna=tmp_path / "r.fa",
                 unidentified=tmp_path / "u.fa", removal_evidence=None,
                 measured={"pairs": 10})
    m = p.as_manifest()
    assert m["trim_r2_source"] == "single-end"
    assert "identification and never by assembly" in m["note"]
    assert m["measured_off_the_reads"] == {"pairs": 10}
    assert json.dumps(m)                       # a manifest that cannot be written is none


def test_the_driver_no_longer_needs_the_contig_flags_from_a_caller():
    """The flags are still on the driver; what changed is who fills them.

    They are not removed from `hapax.py` in this step, because the driver is archived in
    step 7 and editing it changes the hash that pins five committed seed manifests. What
    this asserts is that the package supplies every one of them, so no caller has to.
    """
    from hapax.cli import legacy_argv

    cfg_path = Path(__file__).with_name("_tmp_prepare.toml")
    cfg_path.write_text(BASE.format(label="fixture"))
    try:
        cfg = load(cfg_path)
        prepared = Prepared(trim_r2=3, trim_r2_source="x", assembly=Path("/a.fa"),
                            identified_norrna=Path("/i.fa"), rrna=Path("/r.fa"),
                            unidentified=Path("/u.fa"), removal_evidence=None,
                            measured={"block": 1234})
        argv = legacy_argv(cfg, Path("/r1.fq"), Path("/r2.fq"), Path("/out"),
                           ["cascade"], prepared)
    finally:
        cfg_path.unlink()
    for flag in ("--trim-r2", "--assembly", "--contigs", "--identified-contigs",
                 "--rrna-contigs", "--unidentified-contigs", "--cross-sample-ref"):
        assert flag in argv, f"{flag} is not supplied by the run, so a caller still must"
    assert argv[argv.index("--trim-r2") + 1] == "3"
    # The cross-sample reference is the SAMPLE'S OWN identified contigs. Nothing from
    # another library reaches that flag: the sample's own assembly is the most efficient
    # subtrahend measured here and also removes 0.0% of a genuinely novel virus, so what
    # it is matters more than that it is set.
    assert argv[argv.index("--cross-sample-ref") + 1] == "/i.fa"
    assert argv[argv.index("--block") + 1] == "1234", (
        "a stride measured off the file has to reach the stages, or the measurement was "
        "decoration")
