"""The report refuses rather than printing an absence. PLAN.md step 4.

THE DEFECT THESE PIN. `load()` returns None for an absent file and the panels print what
is missing, so a run that produced almost nothing still built a page that rendered. The
first production pass over an untouched library finished with zero failed stages and
built exactly that page. On a fixture it is a test passing on a broken tool; on real water
it tells a reader something is missing from THEIR SAMPLE when what is missing is a stage
of ours.

A render check cannot catch it, and this is the measured reason the contract is at file
level. Run against that same committed directory, the generator's render check reports all
five sections drawn, because the panels DO reach the DOM: they reach it carrying a census
of an absent file. The file-level contract reports three of the five as not drawable. Both
are true, and only the second is the question a user has.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from hapax.report import ReportError, build, inspect, missing_sections
from hapax.schemas import SCHEMAS, SchemaError, validate
from hapax.stages import R, SECTIONS

ROOT = Path(__file__).resolve().parents[1]
#: A real run directory, produced before three of its stages existed. The point of using
#: it rather than a synthetic one is that its gaps are the gaps production actually had.
PROD = ROOT / "results/prod_SRR23702504/run_s1"


def test_every_declared_output_has_a_schema():
    """An output with no schema is an output nothing checks, which is how a page came to
    render a column of dashes while looking finished."""
    for name, stage in R.stages.items():
        for out in stage.outputs:
            assert out.schema in SCHEMAS, (
                f"{name} declares {out.filename} with schema {out.schema!r}, which is "
                f"not in hapax.schemas.SCHEMAS")


def test_an_unknown_schema_is_refused_rather_than_skipped():
    with pytest.raises(SchemaError, match="no schema named"):
        validate("not_a_schema", {})


def test_a_present_but_empty_file_fails_its_contract():
    """Existence is not the contract. A `counts` of {} is a section with no content."""
    assert validate("layers", {"counts": {}, "totals": {}, "rule": "x",
                               "against_known_positive": {}})
    assert not validate("layers", {"counts": {"a": 1}, "totals": {"n": 1}, "rule": "x",
                                   "against_known_positive": {"a": 1}})


def test_a_measured_absence_is_allowed_where_the_stage_says_so():
    """An explicit census of nothing found is a result; an empty table with no state is
    not. The two look identical in a file, which is why the schema requires the field
    that distinguishes them rather than requiring content."""
    searched_nothing_found = {"state": "searched, no significant subject",
                              "searched": True, "subjects": [],
                              "significance_evalue": 1e-5}
    assert not validate("strains", searched_nothing_found)
    # The same file with the state removed is the case that must fail: an empty table and
    # a search that never ran are then indistinguishable.
    del searched_nothing_found["state"]
    assert validate("strains", searched_nothing_found)


def test_a_tier_table_without_its_operating_point_does_not_render():
    """Every number ships with its operating point or is unreadable, including this one."""
    assert validate("tier_table", {"pairs": 10, "operating_point": {}, "breakdown": {},
                                   "curve": {}, "interval_method": "wilson"})
    assert not validate("tier_table", {"pairs": 10, "operating_point": {"k": 31},
                                       "breakdown": {}, "curve": {},
                                       "interval_method": "wilson"})


@pytest.mark.skipif(not (PROD / "manifest.json").exists(),
                    reason="the committed production run directory is not in this tree")
def test_the_committed_production_run_is_reported_as_incomplete():
    """The run that finished with zero failed stages is three sections short.

    This is the number PLAN.md 1.5 is about, stated the way a user would ask for it: not
    how many stages exited zero, but how many sections of the report could be filled.
    """
    states = inspect(PROD, seed=1)
    absent = missing_sections(states)
    assert set(absent) == {
        "Watchlist screen", "What this water is",
        "How far the assembled viruses sit from anything published"}, sorted(absent)
    # Named by the stage that would fill them, so the message says what to run and not
    # only what is missing.
    blamed = {s.stage for bad in absent.values() for s in bad}
    assert blamed == {"profile", "taxa"}, blamed
    assert len(absent) < len(SECTIONS), "a run with nothing in it would prove nothing"


@pytest.mark.skipif(not (PROD / "manifest.json").exists(),
                    reason="the committed production run directory is not in this tree")
def test_the_build_refuses_an_incomplete_run_and_writes_nothing(tmp_path):
    out = tmp_path / "report.html"
    with pytest.raises(ReportError) as e:
        build(PROD, out=out, seed=1)
    assert "3 of 5 sections cannot be drawn" in str(e.value)
    assert "--allow-partial" in str(e.value)
    assert not out.exists(), (
        "the page was written before the refusal. A file that renders is a file somebody "
        "reads, so the check has to come first.")


def test_a_run_directory_that_is_not_one_is_refused(tmp_path):
    with pytest.raises(ReportError, match="not a hapax run directory"):
        build(tmp_path)


def test_a_manifest_that_does_not_say_which_library_ran_is_refused(tmp_path):
    """A page built from it could not say whose water it is."""
    from hapax.report import descriptor

    (tmp_path / "manifest.json").write_text(json.dumps({"tool": "hapax run"}))
    with pytest.raises(ReportError, match="which library"):
        descriptor(tmp_path, tmp_path / "page.json")


def test_the_descriptor_states_ignorance_rather_than_supplying_facts(tmp_path):
    """An unstated fact is left out of the sentence; it is never another library's.

    The first version of this passed `{}` whenever the run recorded no facts, which is a
    positive claim that nothing is known - and it silenced the registry for a library the
    registry does know, titling the Yeongdok page by its accession and dropping "Yeongdok
    plant effluent, South Korea" from a run whose facts are committed.
    """
    from hapax.report import descriptor

    (tmp_path / "manifest.json").write_text(json.dumps(
        {"tool": "hapax run", "parameters": {"library": "SRR23702504"}}))
    d = json.loads(descriptor(tmp_path, tmp_path / "page.json").read_text())
    assert "sample_facts" not in d, (
        "a run that recorded nothing about its water must leave the key out, so the "
        "registry answers for a library it knows")
    assert d["sample"]["accession"] == "SRR23702504"
    assert set(d["inputs"]) == {"spectrum", "spectrum_big_index"}
    for v in d["inputs"].values():
        assert "not_run" in v, "a declared absence and an undeclared one are not the same"


def test_no_report_section_depends_on_a_panel_a_library_may_not_have():
    """A section fed by a per-library panel could not be drawn for most water.

    Two stages resolve a reference out of `experiments/libraries.json` and skip with a
    committed reason when this library has none: `panel` (the depositors' truth genomes)
    and `crossplant` (the 530-library background). Neither output fills a report section
    today, which is why a run on water in no panel still builds a complete page — and the
    only reason this is not load-bearing is that nobody has attached a section to one.

    If one is attached, the report's own contract makes the absence a refusal, and the
    page for every production sample becomes incomplete for a property of the water
    rather than a hole in the tool. The two absences are not interchangeable and the
    report cannot currently tell them apart, so this asserts the condition instead.
    """
    per_library_panels = {"panel", "crossplant"}
    attached = {name: [o.filename for o in R.stages[name].outputs if o.section]
                for name in per_library_panels if name in R.stages}
    offenders = {k: v for k, v in attached.items() if v}
    assert not offenders, (
        f"{offenders} now fills a report section from a panel a library may not have. "
        f"Either the report must distinguish 'this library has no such panel, and "
        f"libraries.json says so' from 'a stage did not run', or the output stays off "
        f"the page.")


def test_an_absent_enrichment_carries_what_the_run_recorded(tmp_path):
    """The cause of an absence is read from the run, and a missing record says so.

    The registry's `absence_means` explains what the page shows without the file, which
    is always true. It also carried a guess at the cause - "which not every machine has
    staged" - and on SRR34217566 both indexes were staged and the stage failed on a stale
    read-path default. The two findings call for different fixes, so the recorded status
    travels beside the note.
    """
    from hapax import report

    (tmp_path / "enrichments.json").write_text(json.dumps({"stages": [
        {"stage": "workup", "status": "failed", "reason": "FileNotFoundError: berlin"},
        {"stage": "architecture", "status": "skipped", "reason": "it reads the workup"},
    ]}))
    failed = report._why_absent(tmp_path, "workup")
    assert failed["recorded"] == "failed" and "berlin" in failed["recorded_reason"]
    assert report._why_absent(tmp_path, "names")["recorded"].startswith("nothing:")
    # No record at all is its own answer, not the same as a stage nobody ran.
    assert report._why_absent(tmp_path / "elsewhere", "workup")["recorded"].startswith(
        "nothing: this run wrote no")


def test_the_descriptor_declares_the_spectrum_the_run_actually_produced(tmp_path):
    """A hardcoded absence is a claim about this water that nothing measured.

    `inputs.spectrum` read "the divergent-contig workup is not yet a stage of hapax run"
    long after the divergent stage was one. A complete five-stage run of SRR34217566 then
    drew its headline distance spectrum as an empty axis, and told three of that figure's
    controls the stage had not run, beside a run directory holding 1.7 MB of that stage's
    output. It is the harder of the two directions to notice, because a declared absence
    renders as scrupulousness rather than as a gap.
    """
    from hapax.report import descriptor

    (tmp_path / "manifest.json").write_text(json.dumps(
        {"tool": "hapax run", "parameters": {"library": "SRR23702504"}}))
    d = json.loads(descriptor(tmp_path, tmp_path / "page.json", seed=2).read_text())
    assert "not_run" in d["inputs"]["spectrum"], (
        "no divergent output in the run directory, so the absence is the honest reading")
    assert "PLAN.md" not in d["inputs"]["spectrum"]["not_run"], (
        "the reason must describe this run, not the state of the tool when it was typed")

    (tmp_path / "divergent_known_s2.json").write_text(json.dumps({"rows": []}))
    d = json.loads(descriptor(tmp_path, tmp_path / "page.json", seed=2).read_text())
    assert d["inputs"]["spectrum"] == {
        "path": str(tmp_path / "divergent_known_s2.json")}, (
        "the stage ran and its output is sitting there; the page has to be told")


def test_the_report_build_passes_its_seed_to_the_descriptor(tmp_path):
    """Seed 2's page must not be told about seed 1's contigs, or the other way about."""
    import inspect as _inspect

    from hapax import report

    src = _inspect.getsource(report.build)
    assert 'descriptor(rundir, rundir / "page_descriptor.json", seed)' in src, (
        "the descriptor resolves per-seed filenames, so the seed has to reach it")
