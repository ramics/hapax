"""Reads to report, on a fixture, in one command. PLAN.md step 2.

WHAT THIS TEST IS FOR, AND WHY IT FAILS TODAY ON PURPOSE.

Nothing in this repo has ever run from reads to a finished report. The existing suite
stops at a unit or a mocked stage because CI installs no external tool and no index, and
the consequence is measured: the first production pass over a library nobody had touched
finished with **zero failed stages** and built a page that was mostly a census of absent
files. A stage exiting zero is not a filled report section, and until something checks
the second thing, "it ran" means nothing (PLAN.md 1.1, 1.3, 1.5).

So this test asserts the definition of done from PLAN.md 4, against the fixture:

    hapax run --config run.toml R1 R2 --out DIR

one command, from reads, producing a report whose every section is filled from a file
the run itself wrote, with nothing supplied or edited by hand.

That command does not exist yet. The driver today starts AFTER assembly and demands
`--contigs`, `--unidentified-contigs` and `--cross-sample-ref` from its caller, and three
report sections have no stage behind them at all. This test is therefore RED until
PLAN.md step 6, and that is its job: it is the thing that says how much of the rebuild is
left, in the terms a user would use, rather than a plan in a document agreeing with
itself.

HOW IT REPORTS. The work lands in several PRs, so a single pass/fail would say nothing
between them. The test walks four phases, each one a plan step, and fails at the first
one that does not hold with a message naming the step that lands it. A later phase is
not attempted once an earlier one fails, because its failure would be a consequence and
not a finding.

WHERE IT RUNS. The phases need MEGAHIT, DIAMOND, Kraken2 and minimap2, which the normal
`pytest` job does not have (see .github/workflows/tests.yml, and TOOLS.md for the pins).
It therefore SKIPS unless `HAPAX_E2E=1`, which the end-to-end CI job sets after installing
them. The skip is so that a developer with no bioconda environment still gets a green
suite; it is not a way of not knowing, because the CI job runs on every pull request and
is red until step 6.

WHAT IT MUST NEVER DO. The fixture knows which arm every read came from, and that is
ground truth. It is used here to say which output class a missing section lost, and it
never reaches a report: a known distance validates this tool and is not output of it
(SPEC 13.4). Nothing in this file may be relaxed into a report assertion about real water.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures"
TOOLS = ["megahit", "diamond", "kraken2", "kraken2-build", "minimap2"]

#: The report's section list, which is the contract this test checks against.
#:
#: It is one literal in the generator (experiments/make_web.py, `const SUMMARY`) and one
#: here, and `test_the_section_list_matches_the_generator` below asserts they agree, so a
#: section added to the report without a stage to fill it fails this test rather than
#: rendering as an absence. The list is identical on every sample by design: a page that
#: grows and drops sections by sample cannot be compared with another one, and a missing
#: section and a measured absence must not look the same (SPEC 13.4).
SECTIONS = ["Watchlist screen", "What this water is",
            "How far the assembled viruses sit from anything published",
            "The unnamed remainder", "Methods"]

#: Which fixture arm each section needs something from, for the failure message. The
#: fixture carries all three output classes precisely so that a section cannot pass by
#: having nothing to say.
SECTION_NEEDS = {
    "Watchlist screen": "the identifiable arms (MS2, E. coli) and what survives them",
    "What this water is": "the classifier and SILVA arms over the E. coli background",
    "How far the assembled viruses sit from anything published":
        "the divergent arm: mutated PMMoV against the real PMMoV in the protein index",
    "The unnamed remainder": "the unplaceable arm, at 3x, and the depth-1 pile",
    "Methods": "every other section's operating points",
}

requires_e2e = pytest.mark.skipif(
    os.environ.get("HAPAX_E2E") != "1",
    reason="the end-to-end phases need MEGAHIT, DIAMOND, Kraken2 and minimap2; the "
           "end-to-end CI job installs them and sets HAPAX_E2E=1")


def test_the_fixture_is_present_and_pinned():
    """The committed fixture matches its own manifest, byte for byte.

    This one runs everywhere, with no external tool, because a fixture that has drifted
    from its manifest makes every phase below unreadable: a failure would be a question
    about the fixture rather than about the tool. `make_fixture.py` writes the manifest,
    so this is the cheap half of the provenance rule (CLAUDE.md): the content hash pins
    the code whatever the tree looked like, and the file hashes pin the data.
    """
    import hashlib

    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    script = (FIXTURE / "make_fixture.py").read_bytes()
    assert manifest["script_sha256"] == hashlib.sha256(script).hexdigest(), (
        "make_fixture.py has changed since the fixture was built; rebuild it with "
        "--seed {} so the manifest and the data agree".format(manifest["seed"]))
    for rel, rec in manifest["files"].items():
        p = FIXTURE / rel
        assert p.exists(), f"{rel} is in the manifest and not in the tree"
        got = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        assert got == rec["sha256_16"], f"{rel} does not match the manifest"
        assert p.stat().st_size == rec["bytes"], f"{rel} is a different size"


def test_the_fixture_carries_all_three_output_classes():
    """All three classes are in the library, so no section can pass by having no input.

    The measured protein identity of the divergent arm is read rather than asserted
    against a number typed here: the mutation procedure decides it and the manifest
    records it. What is asserted is the BAND - that the arm sits where a protein search
    still finds a relative and a nucleotide search mostly does not, which is the wall
    this project has measured (SPEC 4.2) and where the holdout positive control sits at
    66.5%. An arm at 99% would make the divergent section a tier-1 subtraction, and one
    at 30% would make it unplaceable, and either way the section under test would be
    measuring something else.
    """
    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    identities = [m["protein_identity"]
                  for m in manifest["divergent_arm"]["measured_per_orf"]]
    assert identities, "the divergent arm records no measured identity"
    for ident in identities:
        assert 0.55 <= ident <= 0.85, (
            f"the divergent arm is at {ident:.3f} protein identity, outside the band this "
            f"fixture exists to exercise; change DIVERGENT_NT_RATE and rebuild")
    assert manifest["coverage"]["unplaceable"] <= 5, (
        "the unplaceable arm has to be shallow enough to produce a depth-1 pile, which is "
        "the measured boundary the residue report describes")
    for rel in ("db/nt_known.fa", "db/prot_known.faa", "db/silva_toy.fa"):
        assert (FIXTURE / rel).exists(), rel
    divergent = (FIXTURE / "genomes/divergent_pmmov.fa").read_text()
    assert "NC_003630.1" not in (FIXTURE / "db/nt_known.fa").read_text(), (
        "the real PMMoV is in the nucleotide index, so the divergent arm would be "
        "subtracted as known and the section under test would have no input")
    assert divergent.count(">") == 1


def test_the_section_list_matches_the_generator():
    """One section list, in the generator, and this test reads it rather than restating it.

    A report section with no stage behind it is a hole in the tool and not in the water,
    and the way that hole has twice been introduced is a section added to the renderer
    with nothing to fill it. So the contract is checked in both directions: the list here
    is the list there, and the phases below require every entry of it to be filled.
    """
    src = (ROOT / "experiments" / "make_web.py").read_text()
    start = src.index("const SUMMARY = [")
    literal = src[start:src.index("];", start)]
    for section in SECTIONS:
        assert section in literal, (
            f"the generator's section list does not contain {section!r}; this test and "
            f"the generator disagree about what the report is")
    quoted = literal.count("'") // 2
    assert quoted == len(SECTIONS), (
        f"the generator lists {quoted} sections and this test knows {len(SECTIONS)}; a "
        f"section was added to the report without a stage to fill it, or this list is "
        f"stale")


@requires_e2e
def test_the_external_tools_are_installed():
    """The phases below are meaningless without them, so their absence is a failure here.

    In the end-to-end job this cannot be a skip. A test that skips when its tools are
    missing and runs in a job whose whole purpose is to install them reports a broken
    install as a pass, which is the same shape as a stage exiting zero while writing
    nothing.
    """
    missing = [t for t in TOOLS if shutil.which(t) is None]
    assert not missing, f"not installed: {missing}; see TOOLS.md and the end-to-end job"


@requires_e2e
def test_reads_to_report_in_one_command(tmp_path):
    """The definition of done, PLAN.md 4, against the fixture. Red until step 6.

    Four phases, each a plan step. The first that does not hold fails with the step that
    lands it, and nothing after it is attempted: once the command cannot be invoked, a
    report assertion would fail for a reason that is not about the report.
    """
    dbs = tmp_path / "dbs"
    build = subprocess.run([str(FIXTURE / "build_dbs.sh"), str(dbs)],
                           capture_output=True, text=True, timeout=1800)
    assert build.returncode == 0, (
        f"the toy indexes did not build, so nothing below is readable:\n{build.stderr[-4000:]}")

    out = tmp_path / "run"
    config = tmp_path / "run.toml"
    config.write_text(CONFIG.format(dbs=dbs, fixture=FIXTURE, out=out))

    # PHASE 1, PLAN.md step 3: one installed command that takes a config.
    entry = shutil.which("hapax")
    cmd = ([entry] if entry else [sys.executable, str(ROOT / "hapax.py")])
    run = subprocess.run(
        [*cmd, "run", "--config", str(config),
         str(FIXTURE / "reads/fixture_R1.fastq.gz"),
         str(FIXTURE / "reads/fixture_R2.fastq.gz"), "--out", str(out),
         "--report"],
        capture_output=True, text=True, timeout=5400, cwd=ROOT)
    assert entry is not None, (
        "there is no installed `hapax` entry point; `pip install .` is the whole install "
        "in PLAN.md 4 and the package skeleton is step 3")
    assert run.returncode == 0, (
        "PLAN.md step 3/5: `hapax run --config` did not complete from reads alone. Today "
        "the driver starts after assembly and demands --contigs, --unidentified-contigs "
        f"and --cross-sample-ref from its caller (PLAN.md 1.2).\n{run.stderr[-4000:]}")

    # PHASE 2, PLAN.md step 5: the run assembled its own input, and no stage failed.
    manifest = json.loads((out / "manifest.json").read_text())
    failed = [s for s in manifest.get("stages", []) if s.get("status") == "failed"]
    assert not failed, "stages failed: " + json.dumps(
        [{k: s.get(k) for k in ("stage", "reason")} for s in failed], indent=1)
    empty = [s["stage"] for s in manifest.get("stages", [])
             if s.get("status") not in ("failed", "skipped") and not s.get("outputs")]
    assert not empty, (
        f"these stages exited cleanly and wrote nothing: {empty}. A clean exit code is "
        f"not an output, which is PLAN.md 1.5.")
    assembly = manifest.get("assembly") or {}
    assert assembly.get("assembled_at_run_time") is True, (
        "PLAN.md step 5: the contigs were supplied rather than assembled by the run, so "
        "this is not one command from reads")

    # PHASE 3, PLAN.md step 4: the report builds from the run directory alone, and
    # refuses rather than printing an absence.
    page = out / "report.html"
    report = subprocess.run([*cmd, "report", str(out)], capture_output=True, text=True,
                            timeout=1800, cwd=ROOT)
    assert report.returncode == 0, (
        "PLAN.md step 4/6: the report refused this run. That refusal is correct "
        f"behaviour and the run is what is incomplete.\n{report.stdout[-4000:]}"
        f"\n{report.stderr[-4000:]}")
    assert page.exists(), (
        "`hapax report` returned 0 and wrote no page, which is the one outcome that must "
        "not be possible")

    # PHASE 4, PLAN.md step 6: every section filled from a file the run wrote.
    check = json.loads((out / "render_check.json").read_text())
    drawn = check.get("top") or []
    missing = [s for s in SECTIONS if s not in drawn]
    assert not missing, "\n".join(
        [f"PLAN.md step 6: {len(missing)} of {len(SECTIONS)} report sections did not "
         f"render from a file this run wrote:"]
        + [f"  {s}: needs {SECTION_NEEDS[s]}" for s in missing])
    assert not check.get("sections_not_drawable"), (
        f"these sections could not be drawn from a file this run wrote: "
        f"{check['sections_not_drawable']}. On real water that tells a reader something "
        f"is missing from their sample when what is missing is a stage (SPEC 13.4).")
    # An enrichment this run could not produce is allowed and must be RECORDED with what
    # the page shows instead. The fixture stages no Pfam library, so the per-gene workup
    # and the architecture call are expected here; the point is that each one says what
    # its absence costs rather than leaving a reader to wonder.
    for record in check.get("enrichments_absent") or []:
        assert record["means"], (
            f"{record['file']} is absent with no stated consequence, which is the thing "
            f"a reader cannot weigh")


#: The run configuration, which is the whole of what a new input needs (PLAN.md 2).
#:
#: It is written here rather than committed as tests/fixtures/run.toml because the
#: database paths are build outputs under the test's own temporary directory, and a
#: committed config carrying a path from someone's machine is how a sample fact became a
#: default in roughly twenty scripts (CLAUDE.md).
CONFIG = """\
# The fixture, as one run configuration: everything a new input needs, in one file.
#
# `trim_r2` is STATED here with its basis rather than left to be measured, and the reason
# is the honest one. Measuring a non-template 5' prefix needs reads anchored to a
# reference, and the instrument that does it refuses below 2,000 anchored mates, because
# at lower depth the histogram describes the reference's reach into the water rather than
# the library's chemistry. This fixture's rRNA arm is a 16 kb window, so the measurement
# would refuse, which is correct. What the run does measure off the reads is the pair
# count, the modal read length per mate, the layout and the seed stride.
[run]
seed = 1
# One seed, stated. The seeds read disjoint blocks of the file, so five of them need a
# library five blocks long: 9,987 pairs cannot hold five blocks of 9,000, and the
# measurement refuses a stride that would make them overlap rather than silently
# producing a directory that looks like a five-seed result. A headline result needs five
# seeds; a fixture is not a result.
seeds = 1
threads = 2
pairs = 9000
report_pairs = 4000
out = "{out}"

[inputs]
label = "fixture"
sample_facts = {{ matrix = "simulated reads", place = "the CI fixture", molecule = "RNA", basis = "tests/fixtures/manifest.json: three real genomes, 9,987 simulated pairs, built by make_fixture.py --seed 1. Not water, and no number from it belongs in results/." }}

[library]
trim_r2 = 0
trim_r2_basis = "the fixture's reads are simulated with no non-template prefix (tests/fixtures/make_fixture.py), so any trim would delete real template. The prefix instrument cannot measure it here: the rRNA arm is a 16 kb window and anchoring falls below its 2,000-mate floor."

[databases]
# The FASTA the protein index was built from, so the name map reads the same source the
# search used rather than a second one that could disagree with it.
protein_fasta = "{fixture}/db/prot_known.faa"
nucleotide = "{dbs}/nt_known.mmi"
protein = "{dbs}/prot_known.dmnd"
rrna = "{dbs}/silva_toy.mmi"
classifier = "{dbs}/k2_toy"
"""
