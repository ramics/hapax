"""The three promoted stages, as `hapax run` sees them.

WHY THE DRIVER IS TESTED SEPARATELY FROM THE SCRIPTS. Every defect this driver has had was
in the joint, not in a script: three indexes validated and never passed, a block stride that
was never exposed, a stage that reported `ok` for work it had not done. So what is asserted
here is the joint - that a stage is registered everywhere it has to be, that an unknown
library stops the run before anything runs, that a library with no panel is `skipped` with
the committed reason rather than `ok`, and that the shortlist arrives in the run directory
carrying the published cell.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HAPAX = ROOT / "hapax.py"
ROWS = ROOT / "results" / "holdout_cascade" / "contig_rows_s1.tsv"

# Loaded from its path, not by name: `hapax` is the package now (see tests/legacy.py).
from legacy import driver  # noqa: E402

hapax = driver()

NEW = ("shortlist", "panel", "crossplant")


def test_the_new_stages_are_registered_everywhere_they_have_to_be():
    for name in NEW:
        assert name in hapax.STAGE_ORDER
        assert name in hapax.SCRIPT_STAGES, (
            f"{name} runs a script, so it must be preflighted by calling that script")
        assert hasattr(hapax.Run, name)
    for script in ("removal_shortlist.py", "truth_panel.py", "crossplant_stage.py"):
        assert script in hapax.PREFLIGHT_SCRIPTS
        assert script in hapax.SAMPLE_SPECIFIC, (
            "each of these resolves sample-specific reference data, so the flag that "
            "decides which sample must be one hapax always passes")


def evidence_from_committed_rows(path: Path) -> None:
    cols = ["contig", "contig_nt", "route", "index", "subject_or_taxon", "pident",
            "aln_aa", "evalue", "bitscore", "query_coverage", "matching_bases",
            "confidence", "tier"]
    with open(ROWS) as fh, open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["arm"] != "background" or r["tier"] != "identified":
                continue
            if r["holdout_prot_significant"] != "True":
                continue
            if not r["holdout_prot_identity"] or not r["holdout_prot_coverage"]:
                continue
            row = {"contig": r["contig"], "contig_nt": r["length"], "route": "protein",
                   "index": "viral", "subject_or_taxon": "s",
                   "pident": r["holdout_prot_identity"],
                   "evalue": r["holdout_prot_evalue"],
                   "bitscore": r["holdout_prot_bitscore"],
                   "query_coverage": r["holdout_prot_coverage"],
                   "tier": "identified_norrna"}
            out.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


def run_hapax(tmp_path: Path, *extra: str) -> tuple[int, dict, str]:
    reads = tmp_path / "r.fq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    out = tmp_path / "run"
    proc = subprocess.run(
        [sys.executable, "-I", str(HAPAX), "run", str(reads), str(reads),
         "--seed", "1", "--trim-r2", "14", "--out", str(out), *extra],
        cwd=ROOT, capture_output=True, text=True)
    manifest = {}
    if (out / "manifest.json").exists():
        manifest = json.loads((out / "manifest.json").read_text())
    return proc.returncode, manifest, proc.stdout + proc.stderr


def stage(manifest: dict, name: str) -> dict:
    return next(s for s in manifest["stages"] if s["stage"] == name)


def test_library_is_required():
    proc = subprocess.run([sys.executable, "-I", str(HAPAX), "run", "a", "b", "--seed",
                           "1", "--trim-r2", "14", "--out", "x"],
                          cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode != 0
    assert "--library" in proc.stderr


def test_an_unknown_library_refuses_the_whole_run(tmp_path):
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR00000001",
                                    "--stages", "panel", "crossplant")
    assert code == 2, log
    assert manifest["status"] == "refused to start"
    assert any("is not in experiments/libraries.json" in p
               for p in manifest["problems"])


def test_a_library_with_no_truth_panel_is_skipped_not_ok(tmp_path):
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR23973193",
                                    "--stages", "panel")
    assert code == 0, log
    s = stage(manifest, "panel")
    assert s["status"] == "skipped"
    assert "divergence ladder" in s["reason"]
    assert manifest["complete"] is False, (
        "a run whose ground-truth check did not happen is incomplete, and the manifest "
        "has to say so rather than reading as a clean run")


@pytest.mark.skipif(not ROWS.exists(), reason="the holdout rows are not in this tree")
def test_the_shortlist_arrives_in_the_run_directory_with_the_published_cell(tmp_path):
    ev = tmp_path / "removal_evidence.tsv"
    evidence_from_committed_rows(ev)
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR23973193",
                                    "--removal-evidence", str(ev),
                                    "--stages", "shortlist")
    assert code == 0, log
    s = stage(manifest, "shortlist")
    assert s["status"] == "ok", s.get("reason")
    assert "removal_shortlist.json" in s["outputs"]
    payload = json.loads((tmp_path / "run" / "removal_shortlist.json").read_text())
    assert payload["pile"]["contigs_scored"] == 7601
    assert payload["band"]["contigs"] == 69
    # The null was dropped by name, and the run directory says so rather than omitting it.
    assert payload["null"]["state"].startswith("not supplied")
    assert s["provenance"]["script_sha256"]


def test_the_shortlist_is_skipped_without_its_input(tmp_path):
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR23973193",
                                    "--stages", "shortlist")
    assert code == 0, log
    s = stage(manifest, "shortlist")
    assert s["status"] == "skipped"
    assert "removal-evidence" in s["reason"]


def test_an_assembly_without_minimap2_refuses_before_anything_runs(tmp_path):
    """An assembly to score and no aligner to score it with is caught at preflight.

    Not mid-run: the panel stage is late in the order, and this project has already paid
    for discovering a missing tool after the expensive stages had run.
    """
    fa = tmp_path / "asm.fa"
    fa.write_text(">c\nACGT\n")
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR23421544",
                                    "--assembly", str(fa), "--stages", "panel")
    assert code == 2, log
    assert manifest["status"] == "refused to start"
    assert any("minimap2" in p for p in manifest["problems"]), manifest["problems"]


def test_the_source_digest_is_read_once_not_once_per_contig(tmp_path, monkeypatch):
    """The tier-2 stage hashes its input file once, whatever the contig count.

    THE DEFECT THIS PINS, measured on the first production library. Every row of
    `contigs.jsonl` carries the digest of the fasta it was read from, and the digest was
    taken inside the per-contig writer. So the stage re-read its whole input once per
    contig: 194,703 unidentifiable contigs against a 104 MB file is ~20 TB of SHA-256 for
    a stage whose own arithmetic is a length and a GC fraction. It ran 1h45 and was about
    two thirds done, and it is a PER-SEED cost, so five seeds would have spent most of a
    day on it.

    The output is byte-identical either way, which is why nothing caught it: only a clock
    tells, and a correct-looking run directory is exactly the failure mode this repo keeps
    being caught by. So the assertion is on the number of digests, not on the timing - a
    count is deterministic and a timing is not.
    """
    fa = tmp_path / "unid.fa"
    fa.write_text("".join(f">c{i}\nACGTACGTGG\n" for i in range(50)))
    calls = []
    real = hapax.sha256
    monkeypatch.setattr(hapax, "sha256", lambda p: calls.append(str(p)) or real(p))

    run = hapax.Run.__new__(hapax.Run)
    run.args = argparse.Namespace(unidentified_contigs=[fa])
    run.dir = tmp_path / "run"
    run.dir.mkdir()
    run.preflighting = False
    run.stages = []
    entry = hapax.Run.contigs(run)

    assert entry["status"] == "ok" and entry["rows"] == 50
    assert calls.count(str(fa)) == 1, (
        f"the input fasta was hashed {calls.count(str(fa))} times for 50 contigs; it is "
        f"one file read once, and this is how the stage came to cost hours")
    rows = [json.loads(l) for l in (run.dir / "contigs.jsonl").read_text().splitlines()]
    assert len({r["operating_point"]["source_sha256_16"] for r in rows}) == 1


def test_every_library_aware_script_is_told_which_library():
    """The driver requires --library and must hand it to every script that takes one.

    THE DEFECT, measured on the first production library. `hapax run --library X` passed
    `--library` to two of the thirteen stages. The other scripts take it through
    `library_facts.add_argument`, whose default is the discriminating dataset so that every
    committed command reproduces unchanged - so a stage that was not told asked for
    BERLIN's facts while running another library's reads. The `linkage` stage needed a
    block stride, was not given one, and refused because the read path named SRR23702504
    and `--library` said SRR23973193: the path check catching it one layer in. Under a
    neutral read path there is no disagreement to catch and it would have used Berlin's
    stride in silence.

    The set is checked against the scripts' own source rather than retyped, so a script
    that becomes library-aware without being added fails here.
    """
    stage_scripts = set(hapax.SAMPLE_SPECIFIC) | set(hapax.PREFLIGHT_SCRIPTS)
    takes_library = set()
    for script in sorted(stage_scripts):
        src = (ROOT / "experiments" / script).read_text()
        if "library_facts.add_argument" in src or 'add_argument("--library"' in src:
            takes_library.add(script)
    assert takes_library == hapax.LIBRARY_AWARE, (
        f"LIBRARY_AWARE disagrees with the scripts: missing {takes_library - hapax.LIBRARY_AWARE}, "
        f"stale {hapax.LIBRARY_AWARE - takes_library}")


def test_the_library_reaches_a_stage_that_was_not_passed_one(tmp_path):
    """`Run.call` appends --library, so no stage falls back to the default library."""
    run = hapax.Run.__new__(hapax.Run)
    run.args = argparse.Namespace(library="SRR23702504", report_without=[])
    run.preflighting = True
    run.problems = []
    seen = {}
    run._preflight_script = lambda name, script, argv: seen.setdefault(script, argv)
    for script in sorted(hapax.LIBRARY_AWARE):
        # every sample-specific flag the script requires, so the refusal above does not
        # short-circuit the call and leave the assertion vacuous
        argv = ["--seed", "1"]
        for flag, _field in hapax.SAMPLE_SPECIFIC.get(script, []):
            if flag != "--library":
                argv += [flag, "x"]
        hapax.Run.call(run, "s", script, argv, {}, None)
    assert set(seen) == hapax.LIBRARY_AWARE, f"only reached {sorted(seen)}"
    for script, argv in seen.items():
        assert "--library" in argv, f"{script} ran without --library"
        assert argv[argv.index("--library") + 1] == "SRR23702504"


def test_the_driver_has_no_private_git_helper_that_can_return_nothing():
    """A manifest must not be able to record an empty commit or an unchecked clean tree.

    The driver carried its own four-line `git()` that returned "" on any exception. On an
    instance, which has a payload tarball and no git, that wrote `git_sha: ""` into every
    seed manifest of `results/prod_SRR23702504/` and, worse, `tree_dirty: false`, because
    `bool("")` is False - a positive claim about a tree nobody looked at. `provenance.py`
    exists for that machine and raises rather than guessing. This is the shuffle tripwire
    applied to provenance: a new private copy fails here.
    """
    src = (ROOT / "hapax.py").read_text()
    assert "def git(" not in src, "hapax.py has grown a private git helper again"
    assert "from provenance import git_sha, tree_dirty" in src
    import provenance
    assert hapax.git_sha is provenance.git_sha
    assert hapax.tree_dirty is provenance.tree_dirty


def test_provenance_refuses_rather_than_reporting_a_clean_tree_it_cannot_see(monkeypatch):
    """With no git and no commit passed in, both calls raise instead of answering."""
    import provenance
    monkeypatch.setattr(provenance, "_git", lambda *a: None)
    monkeypatch.delenv(provenance.SHA_ENV, raising=False)
    for fn in (provenance.git_sha, provenance.tree_dirty):
        with pytest.raises(provenance.ProvenanceError):
            fn()
    # and with the launcher's commit present, the tree is clean on that commit's authority
    monkeypatch.setenv(provenance.SHA_ENV, "deadbee")
    assert provenance.git_sha() == "deadbee"
    assert provenance.tree_dirty() is False

def test_a_deselected_stage_records_whose_choice_it_was():
    """`--stages` defaults to every stage, so a missing stage was asked for.

    The holdout fixture's run directory carries four stages recorded as "not selected
    with --stages", and the page built from it then reports an incomplete run without
    saying that the subset was typed by hand rather than forced by the library. Two of
    those four (linkage, census) this library forbids anyway - both are mate-linkage
    instruments and it is single-end - so one reason was hiding the other. A deselected
    stage now names the default it departs from, and the library's own constraint where
    there is one.
    """
    src = HAPAX.read_text()
    assert 'reason="not selected with --stages"' not in src
    assert "deselected on the command line" in src
    assert "could not have run on this library either" in src


def test_every_report_section_has_a_stage_that_fills_it():
    """A report section with no stage behind it is a hole in the tool, not in the water.

    The first production pass over an untouched library built a page whose every summary
    section was a census of an absent file, because the three scripts that fill them grew
    as experiments over one sample and were never adopted into `hapax run`: the watchlist
    screen and the composition read public_health_profile_s*.json, the distance line reads
    known_strain_divergence_s*.json, and the unnamed remainder reads ranking.json. Each
    now has a stage, and this test fails if a later section arrives without one.
    """
    src = HAPAX.read_text()
    for script, stage in (("public_health_profile.py", "profile"),
                          ("known_strain_divergence.py", "taxa"),
                          ("concern_ranking.py", "ranking"),
                          ("divergent_known.py", "divergent")):
        assert script in src, script
        assert f'"{stage}"' in src, stage
        assert f"def {stage}(self)" in src, stage
    # And they are in the default run, not only available behind a flag.
    import re
    order = re.search(r"STAGE_ORDER = \[(.*?)\]", src, re.S).group(1)
    for stage in ("profile", "taxa", "ranking"):
        assert f'"{stage}"' in order, stage


#: Every flag that names an index, and the config key that supplies it. Which index was
#: consulted is an operating point: a removal fraction means nothing without it, and one
#: Kraken2 database removes 42.8%, 21.7% or 3.5% of the same library depending only on
#: its confidence setting. Each of these defaulted to a path under this checkout's own
#: `data/`, which is the `default=14` fault applied to the quantity that matters most.
INDEX_FLAGS = {"--k2-dbs": "classifier", "--silva": "rrna", "--dmnd": "protein",
               "--pfam": "pfam"}


def test_no_index_flag_carries_a_default():
    """A default index is a guessed operating point, and this repo has withdrawn a figure
    whose index silently differed from the one named."""
    parser_actions = {}
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()):
        with pytest.raises(SystemExit):
            sys.argv = ["hapax", "run", "--help"]
            hapax.main()
    # Read the declarations rather than the help text: the help is prose and the default
    # is the thing that fires.
    source = HAPAX.read_text()
    for flag in INDEX_FLAGS:
        start = source.index(f'r.add_argument("{flag}"')
        decl = source[start:source.index("r.add_argument", start + 10)]
        assert "ROOT /" not in decl, (
            f"{flag} defaults to a path inside this checkout. A caller who forgets it "
            f"would measure against whatever that path held and the manifest would name "
            f"what was asked for")
        assert "default=[]" in decl or "default=None" in decl, flag
        parser_actions[flag] = decl
    assert len(parser_actions) == len(INDEX_FLAGS)


@pytest.mark.parametrize("flag,key", sorted(INDEX_FLAGS.items()))
def test_a_stage_refuses_the_index_it_was_not_given_naming_the_config_key(flag, key,
                                                                         tmp_path):
    """Leaving the flag off has the same end as passing the wrong index.

    The stage script below carries its own default, so a tier with no index given still
    runs, still removes some fraction of the library, and reports it against an index
    nobody chose. The refusal has to name the setting, because a reader of a message that
    quotes `data/pfam/Pfam-A.hmm` has to work out for themselves what that file is.

    The refusal arrives from the preflight, before anything runs, which is where a run
    that would not have done what was asked is supposed to stop.
    """
    stage_for = {"--k2-dbs": "cascade", "--silva": "cascade", "--dmnd": "cascade",
                 "--pfam": "divergent"}
    # Every index except the one under test, so the message is about that one. The paths
    # only have to exist: what is being tested is the absence, not the content.
    others = []
    for other in INDEX_FLAGS:
        if other == flag:
            continue
        # `--dmnd` is an ORDERED PAIR, viral then Swiss-Prot, and the distance stage
        # refuses a list of any other length: it reports which of the two produced the
        # nearest relative, which a single index cannot answer.
        names = ["viral", "sprot"] if other == "--dmnd" else ["db"]
        paths = []
        for name in names:
            stand_in = tmp_path / f"{other.strip('-')}_{name}.db"
            stand_in.write_text("x")
            paths.append(str(stand_in))
        others += [other, *paths]
    # The `divergent` stage skips before it looks at a profile library when there is no
    # identified pile to measure the distance over, and that order is right, so this
    # gives it one.
    pile = tmp_path / "identified.fa"
    pile.write_text(">c\nACGT\n")
    code, manifest, log = run_hapax(tmp_path, "--library", "SRR23973193",
                                    "--stages", stage_for[flag], "--contigs",
                                    "--identified-contigs", str(pile), *others)
    assert code == 2, log
    assert manifest["status"] == "refused to start"
    problems = " ".join(manifest["problems"])
    assert flag in problems, problems
    assert key in problems, (
        f"the refusal names {flag} and not the config key that supplies it; {problems}")
