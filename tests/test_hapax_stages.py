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

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
HAPAX = ROOT / "hapax.py"
ROWS = ROOT / "results" / "holdout_cascade" / "contig_rows_s1.tsv"

import hapax  # noqa: E402

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
    assert any("will not guess" in p for p in manifest["problems"])


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
