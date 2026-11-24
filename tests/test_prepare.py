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


def _cfg(tmp_path, label="fixture", extra="", base=BASE):
    p = tmp_path / "run.toml"
    p.write_text(base.format(label=label) + extra)
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
                           ["cascade"], prepared,
                           {"classify": Path("/k2/classify"),
                            "diamond": Path("/b/diamond"),
                            "minimap2": Path("/b/minimap2")})
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


# ------------------------------------------------------- where the external binaries are


def test_a_stated_tool_path_is_used_as_it_is(tmp_path):
    """Stating it is the override, so it is not searched for or checked against PATH."""
    from hapax import tools

    binary = tmp_path / "classify"
    binary.write_text("#!/bin/sh\n")
    cfg = _cfg(tmp_path, extra=f'\n[tools]\nclassify = "{binary}"\n')
    assert tools.resolve(cfg, "classify") == binary


def test_a_stated_tool_path_that_is_not_there_is_refused(tmp_path):
    from hapax import tools

    cfg = _cfg(tmp_path, extra='\n[tools]\nclassify = "/nowhere/classify"\n')
    with pytest.raises(ConfigError, match="nothing there"):
        tools.resolve(cfg, "classify")


def test_a_missing_tool_is_refused_and_names_how_to_supply_it(tmp_path, monkeypatch):
    """The fault this replaces: nine scripts default --classify-bin to a path inside one
    container, so a machine without that directory got a FileNotFoundError out of a
    subprocess rather than a statement about its own installation."""
    from hapax import tools

    monkeypatch.setattr(tools.shutil, "which", lambda name: None)
    cfg = _cfg(tmp_path)
    with pytest.raises(ConfigError, match=r"\[tools\]"):
        tools.resolve(cfg, "classify")
    with pytest.raises(ConfigError, match="minimap2"):
        tools.resolve(cfg, "minimap2")


def test_every_absent_tool_is_named_in_one_refusal(tmp_path, monkeypatch):
    from hapax import tools

    monkeypatch.setattr(tools.shutil, "which", lambda name: None)
    cfg = _cfg(tmp_path)
    with pytest.raises(ConfigError) as e:
        tools.resolve_all(cfg, ("classify", "diamond", "minimap2"))
    for name in ("classify", "diamond", "minimap2"):
        assert name in str(e.value), (
            "one refusal per missing tool is the same fix discovered three times")


def test_classify_is_found_through_the_kraken2_wrapper(tmp_path, monkeypatch):
    """`classify` is on no distribution's PATH; the wrapper records where it is.

    Written as a test rather than trusted because the two layouts this repo has met put
    it in different places: bioconda in <prefix>/libexec, the staged tarball in
    <prefix>/share/kraken2-<version>/libexec.
    """
    from hapax import tools

    prefix = tmp_path / "env"
    (prefix / "bin").mkdir(parents=True)
    libexec = prefix / "share" / "kraken2-2.17.2" / "libexec"
    libexec.mkdir(parents=True)
    (libexec / "classify").write_text("binary")
    wrapper = prefix / "bin" / "kraken2"
    wrapper.write_text('#!/usr/bin/env perl\nmy $KRAKEN2_DIR = "%s";\n' % libexec)
    monkeypatch.setattr(tools.shutil, "which",
                        lambda name: str(wrapper) if name == "kraken2" else None)
    cfg = _cfg(tmp_path)
    assert tools.resolve(cfg, "classify") == libexec / "classify"


def test_only_the_drivers_own_stages_reach_its_stages_flag(tmp_path):
    """The driver's `--stages` is an argparse choice list, so a stage it does not have is
    an exit code 2 before anything runs rather than a stage that is skipped.

    Three stages were adopted by the package in PLAN.md step 6 and `report` was never the
    driver's, and the translation excluded `report` by name. A by-name exclusion has to be
    kept in step with every later adoption, and was not.
    """
    from legacy import driver

    from hapax.cli import legacy_argv
    from hapax.stages import R

    cfg = _cfg(tmp_path)
    argv = legacy_argv(cfg, Path("/r1.fq"), Path("/r2.fq"), tmp_path, R.production())
    passed = argv[argv.index("--stages") + 1:]
    assert passed, "no stage reached the driver"
    unknown = sorted(set(passed) - set(driver().STAGE_ORDER))
    assert not unknown, f"the driver has no such stage: {unknown}"
    assert set(passed) == {n for n in R.production() if R.stages[n].in_driver}


def test_the_two_protein_tiers_reach_the_driver_in_order(tmp_path):
    """Viral then Swiss-Prot. The distance stage reports which tier held the nearest
    relative and refuses a list of any other length, so the order is part of the answer
    and cannot come from iterating a table."""
    from hapax.cli import legacy_argv

    cfg = _cfg(tmp_path, extra='\n[databases]\nsprot = "/b/sprot.dmnd"\n'
                               'protein = "/a/viral.dmnd"\n')
    argv = legacy_argv(cfg, Path("/r1.fq"), None, tmp_path, ["cascade"])
    at = argv.index("--dmnd")
    assert argv[at + 1:at + 3] == ["/a/viral.dmnd", "/b/sprot.dmnd"]


def test_a_field_left_out_on_purpose_is_stated_and_passed_on(tmp_path):
    """A field whose reference is absent stops the run, because a dropped field reads as
    a zero. Naming it makes the omission an operating point rather than a hole."""
    from hapax.cli import legacy_argv

    cfg = _cfg(tmp_path, extra='\n[databases]\n',
               base='[run]\nseed = 1\nseeds = 1\nreport_without = ["crispr"]\n\n'
                    '[inputs]\nlabel = "{label}"\n')
    argv = legacy_argv(cfg, Path("/r1.fq"), None, tmp_path, ["reads"])
    assert argv[argv.index("--report-without") + 1] == "crispr"
    assert "--report-without" not in legacy_argv(
        _cfg(tmp_path), Path("/r1.fq"), None, tmp_path, ["reads"])


def test_the_manifest_says_whether_the_run_assembled_its_own_contigs(tmp_path):
    """`assembled_at_run_time`, in the file a reader opens.

    The convention already existed and what it is worth is visible in
    results/holdout_cascade/, whose assembly recorded it as false and pointed at a file
    that does not exist: only the contig names survive, so the command line behind a
    committed result cannot be recovered. A run that builds its own contigs and does not
    say so is indistinguishable from one handed them by a caller.
    """
    from hapax.prepare import Prepared

    p = Prepared(trim_r2=0, trim_r2_source="x", assembly=Path("/a.fa"),
                 identified_norrna=Path("/i.fa"), rrna=Path("/r.fa"),
                 unidentified=Path("/u.fa"), removal_evidence=None,
                 assembly_record={"assembled_at_run_time": True, "contigs": 5})
    assert p.as_manifest()["assembly_record"]["assembled_at_run_time"] is True
    assert Prepared(trim_r2=0, trim_r2_source="x", assembly=Path("/a.fa"),
                    identified_norrna=Path("/i.fa"), rrna=Path("/r.fa"),
                    unidentified=Path("/u.fa"),
                    removal_evidence=None).as_manifest()["assembly_record"] == {}, (
        "an empty record must not read as an assembly this run built")
