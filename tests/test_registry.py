"""The stage registry and the run configuration. PLAN.md step 3.

What these tests are defending. Three of the faults in PLAN.md 1 are not bugs in a
calculation; they are things nothing was in a position to check. A stage order typed in
three places could disagree with itself and did. A report section could name a file no
stage writes, and three did. A per-library fact could reach a run as a default rather than
as a measurement, and one did for twenty scripts. Each of those is now a declaration, and
a declaration is checkable.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from hapax import __version__
from hapax.config import ConfigError, load
from hapax.registry import Output, Registry, RegistryError, Stage
from hapax.stages import ORPHAN_REPORT_INPUTS, R, REPORT_INPUTS, SECTIONS

ROOT = Path(__file__).resolve().parents[1]
MAKE_WEB = ROOT / "experiments" / "make_web.py"


def test_the_order_is_derived_and_respects_every_declaration():
    order = R.production()
    assert order, "no production stages"
    where = {name: i for i, name in enumerate(order)}
    for name in order:
        for need in R.stages[name].needs:
            assert where[need] < where[name], (
                f"{name} runs before {need}, which it needs")


def test_the_order_is_the_same_every_time():
    """A run whose stage order depends on dictionary iteration has a manifest that cannot
    be compared with another run's."""
    assert R.production() == R.production() == R.production()


def test_a_cycle_is_refused_rather_than_ordered_somehow():
    reg = Registry()
    reg.add(Stage(name="a", why="x", needs=("b",)))
    reg.add(Stage(name="b", why="x", needs=("a",)))
    with pytest.raises(RegistryError, match="cycle"):
        reg.order()


def test_a_stage_declared_twice_is_refused():
    reg = Registry()
    reg.add(Stage(name="a", why="x"))
    with pytest.raises(RegistryError, match="twice"):
        reg.add(Stage(name="a", why="x"))


def test_a_subset_that_drops_a_need_is_refused_not_quietly_completed():
    """`--stages cascade` must not silently run four other stages.

    A subset run that is not the subset the caller named is how the old `--stages` came to
    hide the reason a stage was missing: four stages were recorded as "not selected" on a
    library that forbids two of them anyway, and one reason hid the other.
    """
    with pytest.raises(RegistryError, match="needs"):
        R.order(["ranking"])          # needs reads and strains
    assert R.order(["cascade"]) == ["cascade"]


def test_the_validation_stages_are_not_in_the_default_set():
    """SPEC 13.4. A stage that scores a run against sequence whose distance was fixed
    before it has no answer on water nobody has deposited from, which is every sample
    production is for."""
    validation = {n for n, s in R.stages.items() if s.validation}
    assert validation == {"calibration", "panel"}
    assert not validation & set(R.production())
    # Still selectable by name, which is how a fixture library runs them.
    assert "panel" in R.order(["cascade", "panel"])


def test_the_registry_agrees_with_the_legacy_driver_about_production():
    """One list, not two. The old one was typed in three places and the launcher's copy
    is the one that ran, so a validation stage ran in a production pass."""
    from legacy import driver

    hapax = driver()
    # `report` is a stage here and was not one in the driver: the page was built by hand
    # afterwards, which is why a run could finish with every stage green and no page.
    mine = {n for n in R.production()
            if R.stages[n].in_driver and n != "report"}
    assert mine == set(hapax.PRODUCTION_STAGES), (
        f"registry: {sorted(mine)}\n"
        f"driver:   {sorted(hapax.PRODUCTION_STAGES)}")


def test_a_stage_the_package_runs_is_marked_as_the_packages():
    """A stage marked as the driver's that the driver does not have would run nowhere.

    That is not hypothetical: the stage list was typed in three places and the launcher's
    copy is the one that ran, so a stage left out of one copy ran nowhere and nothing said
    so.
    """
    from legacy import driver

    declared_by_driver = set(driver().STAGE_ORDER)
    for name, stage in R.stages.items():
        if name == "report":
            continue
        if stage.in_driver:
            assert name in declared_by_driver, (
                f"{name} is marked as the driver's stage and the driver does not declare "
                f"it, so nothing runs it")
        else:
            assert name not in declared_by_driver, (
                f"{name} is marked as the package's stage and the driver also declares "
                f"it, so it would run twice")


def test_every_report_section_is_filled_by_some_stage():
    """A section with no stage behind it is a hole in the tool and not in the water."""
    filled = R.sections()
    missing = [s for s in SECTIONS if s not in filled]
    assert not missing, f"no stage fills: {missing}"


def test_the_section_names_are_the_generators_own():
    literal = MAKE_WEB.read_text()
    start = literal.index("const SUMMARY = [")
    block = literal[start:literal.index("];", start)]
    for section in SECTIONS:
        assert section in block, f"{section!r} is not a section of the report"
    assert block.count("'") // 2 == len(SECTIONS)


def _what_the_report_reads() -> set[str]:
    """Every run-directory file the generator reads, taken from its source.

    Read off the source rather than listed by hand, so that a read added to the renderer
    without a producer fails this suite instead of rendering as an absence on a page
    about somebody's water.
    """
    src = MAKE_WEB.read_text()
    found = set(re.findall(r'rundir / "([^"]+)"', src))
    found |= set(re.findall(r'_one\(rundir,\s*\n?\s*"([^"]+)"', src))
    found |= set(re.findall(r'run_rows\.jsonl\(rundir, "([^"]+)"', src))
    # A directory the page is written into, not a file it reads.
    return {f for f in found if "." in f}


def test_every_file_the_report_reads_is_declared():
    declared = set(REPORT_INPUTS) | set(ORPHAN_REPORT_INPUTS)
    undeclared = sorted(_what_the_report_reads() - declared)
    assert not undeclared, (
        f"the report reads these out of a run directory and nothing declares them: "
        f"{undeclared}. Declare the stage that writes each, or add it to "
        f"ORPHAN_REPORT_INPUTS with the script that would.")


def test_every_declared_report_input_has_a_producing_stage():
    """The census PLAN.md step 6 emptied.

    It held six. Three were sections rendering a census of an absent file, and are stages
    now. Three more, found by writing the map down, were enrichments whose sections render
    without them: stages now too, with optional outputs that state what their absence
    costs. The two kinds have to look different, because treating an enrichment as a hole
    would make every run without a Pfam database refuse, and treating a hole as an
    enrichment is what the page was doing when a run with three empty sections built
    cleanly.
    """
    produced = {o.filename for s in R.stages.values() for o in s.outputs}
    produced |= {"manifest.json", "blind_spots.json"}        # written by the run itself
    for pattern, who in REPORT_INPUTS.items():
        canonical = pattern.replace("_s*.json", "_s{seed}.json")
        assert canonical in produced or who == "the run itself", (
            f"{pattern} is declared as coming from {who}, and no stage output matches it")
    assert ORPHAN_REPORT_INPUTS == {}, (
        f"a report input has no producing stage: {sorted(ORPHAN_REPORT_INPUTS)}. That is "
        f"a section that will render an absence on a page about somebody's water.")


def test_an_optional_output_has_to_say_what_its_absence_costs():
    """An absence with no stated consequence is the thing a reader cannot weigh."""
    from hapax.registry import Output

    with pytest.raises(ValueError, match="absence_means"):
        Output("x.json", "ranking", optional=True)
    for name, stage in R.stages.items():
        for out in stage.outputs:
            if out.optional:
                assert len(out.absence_means or "") > 40, (
                    f"{name}/{out.filename} is optional and does not say what the page "
                    f"shows instead")


def test_the_enrichments_are_optional_and_the_measurements_are_not():
    """Which of the six were holes and which were enrichments, asserted rather than
    remembered. Moving one across this line changes whether a run refuses."""
    optional = {out.filename for s in R.stages.values() for out in s.outputs
                if out.optional}
    assert optional == {"divergent_contig_workup_s{seed}.json", "architecture_s{seed}.json",
                        "subject_names.json"}, sorted(optional)
    required_sections = {out.filename for s in R.stages.values() for out in s.outputs
                         if out.section and not out.optional}
    for f in ("public_health_profile_s{seed}.json",
              "known_strain_divergence_s{seed}.json", "ranking.json"):
        assert f in required_sections, (
            f"{f} filled a section that used to render empty; if it becomes optional a "
            f"run with three empty sections builds cleanly again")


def test_every_stage_says_why_it_exists():
    """CLAUDE.md: the reasoning is the research. A stage with no reason is a stage nobody
    can tell is still needed."""
    for name, stage in R.stages.items():
        assert len(stage.why) > 60, f"{name} does not say why it exists"


# ----------------------------------------------------------------------- the config file

GOOD = """
[run]
seed = 2
threads = 8

[inputs]
label = "fixture"

[databases]
protein = "/tmp/prot.dmnd"
"""


def test_a_config_states_what_a_run_needs(tmp_path):
    p = tmp_path / "run.toml"
    p.write_text(GOOD)
    cfg = load(p)
    assert cfg.seed == 2 and cfg.threads == 8
    assert cfg.label == "fixture"
    assert cfg.database("protein").name == "prot.dmnd"


def test_a_database_that_was_not_named_is_refused_not_guessed(tmp_path):
    p = tmp_path / "run.toml"
    p.write_text(GOOD)
    with pytest.raises(ConfigError, match="operating point"):
        load(p).database("classifier")


def test_an_unknown_key_is_refused_rather_than_ignored(tmp_path):
    """An ignored setting is one the caller believes is in force and is not."""
    p = tmp_path / "run.toml"
    p.write_text(GOOD + '\n[run]\n' if False else GOOD.replace("threads = 8",
                                                               "thread = 8"))
    with pytest.raises(ConfigError, match="does not know"):
        load(p)
    p.write_text(GOOD + "\n[nonsense]\nx = 1\n")
    with pytest.raises(ConfigError, match="tables this version does not know"):
        load(p)


def test_a_stated_library_fact_without_its_basis_is_refused(tmp_path):
    """The registry's rule, moved rather than dropped.

    A config copied from another run's config is a default with extra steps, and a trim
    of 14 on a library with no non-template prefix deletes fourteen real bases from every
    mate, silently, in the direction of less recognition.
    """
    p = tmp_path / "run.toml"
    p.write_text(GOOD + '\n[library]\ntrim_r2 = 14\n')
    with pytest.raises(ConfigError, match="trim_r2_basis"):
        load(p)
    p.write_text(GOOD + '\n[library]\ntrim_r2 = 14\ntrim_r2_basis = "measured in '
                        'results/nontemplate_prefix_map/"\n')
    assert load(p).fact("trim_r2") == 14


def test_a_basis_for_a_value_nobody_stated_is_refused(tmp_path):
    p = tmp_path / "run.toml"
    p.write_text(GOOD + '\n[library]\ntrim_r2_basis = "measured somewhere"\n')
    with pytest.raises(ConfigError, match="describes nothing"):
        load(p)


def test_an_unstated_fact_is_absent_and_never_another_librarys(tmp_path):
    p = tmp_path / "run.toml"
    p.write_text(GOOD)
    assert load(p).fact("trim_r2") is None, (
        "an unstated trim must be absent, so that a stage measures it or refuses; the "
        "one thing it must never be is the discriminating dataset's 14")


def test_the_version_is_available_for_a_manifest_to_record():
    assert re.match(r"^\d+\.\d+\.\d+", __version__), __version__
