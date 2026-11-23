"""The cross-plant stage's rules, asserted as refusals rather than trusted as documentation.

Four of the five rules can be checked without a pass: a query library inside the panel it is
being subtracted against, a window width that disagrees with the panel's query length, a
report missing the paired per-plant strata or an operating point, and a library with no panel
at all. The fifth - read down the depth column - is a sentence, and the test asserts it
travels in the output rather than only in the report.

The committed `results/crossplant/plant_summary.tsv` is used as the positive case for the
strata check, so this stage's idea of a readable report is the published one and not a
stricter rule invented here.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
SCRIPT = ROOT / "experiments" / "crossplant_stage.py"
SUMMARY = ROOT / "results" / "crossplant" / "plant_summary.tsv"
PANEL = ROOT / "results" / "crossplant" / "panel.tsv"

import crossplant_stage  # noqa: E402


def run(tmp_path: Path, *extra: str) -> tuple[int, dict, str]:
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, "-I", str(SCRIPT), "--work",
                           str(tmp_path / "work"), "--out", str(out), *extra],
                          cwd=ROOT, capture_output=True, text=True)
    payload = {}
    if (out / "crossplant_stage.json").exists():
        payload = json.loads((out / "crossplant_stage.json").read_text())
    return proc.returncode, payload, proc.stdout + proc.stderr


@pytest.mark.skipif(not SUMMARY.exists(), reason="the committed report is not in this tree")
def test_the_published_report_satisfies_the_strata_check():
    with open(SUMMARY) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert crossplant_stage.missing_paired_rows(rows) == []


def test_a_report_without_the_paired_strata_is_rejected():
    rows = [{"op": "most1", "stratum": "within_plant_other_water"},
            {"op": "most1", "stratum": "other_plant_other_project"}]
    missing = crossplant_stage.missing_paired_rows(rows)
    assert any("paired_within_minus_other_project" in m for m in missing)
    assert any(m.startswith("any1/") for m in missing), (
        "both operating points are required: they disagree by a factor of four")


@pytest.mark.skipif(not PANEL.exists(), reason="the panel table is not in this tree")
def test_a_query_library_inside_the_panel_refuses():
    """Rule 1. Own-library reads are abundant by construction, so this is a refusal.

    Asserted against the real panel's own first member rather than a fixture library: the
    rule has to fire on the libraries that are actually in there.
    """
    with open(PANEL) as fh:
        members = {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}
    member = next(iter(members))
    why = crossplant_stage.own_library_refusal(member, members)
    assert why and "its own library" in why and "54%" in why
    # And the query sample of the committed pass is not a member, which is the other half
    # of the same rule: it was held out by construction.
    assert crossplant_stage.own_library_refusal("SRR23973193", members) is None


@pytest.mark.skipif(not PANEL.exists(), reason="the panel table is not in this tree")
def test_a_window_width_that_disagrees_with_the_panel_refuses(tmp_path):
    """Rule 3. One width across both arms, or the arms are two instruments."""
    code, _, log = run(tmp_path, "--library", "SRR23973193", "--stages", "queries",
                       "--contig-dir", str(tmp_path), "--window", "150")
    assert code == 2, log
    assert "two instruments" in log


def test_a_library_with_no_panel_is_skipped_with_its_reason(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR23421544", "--stages", "queries",
                         "--contig-dir", str(tmp_path))
    assert code == 0, log
    assert out["status"] == "skipped"
    assert "pooled reading" in out["reason"]


def test_an_unknown_library_refuses(tmp_path):
    code, _, log = run(tmp_path, "--library", "SRR12121212", "--stages", "queries",
                       "--contig-dir", str(tmp_path))
    assert code == 2, log
    assert "is not in experiments/libraries.json" in log


@pytest.mark.skipif(not PANEL.exists(), reason="the panel table is not in this tree")
def test_the_queries_arm_cuts_panel_width_windows_and_records_the_panel(tmp_path):
    """Rule 3 the other way round: the arm that runs must use the panel's own width."""
    cdir = tmp_path / "contigs"
    cdir.mkdir()
    (cdir / "sample_contigs_unidentified.fa").write_text(
        "".join(f">u{i}\n{'ACGT' * 200}\n" for i in range(5)))
    (cdir / "sample_contigs_identified.fa").write_text(
        "".join(f">i{i}\n{'ACGT' * 200}\n" for i in range(5)))
    code, out, log = run(tmp_path, "--library", "SRR23973193", "--stages", "queries",
                         "--contig-dir", str(cdir))
    assert code == 0, log
    q = out["stages"]["queries"]
    assert q["state"] == "built" and q["window_nt"] == 90 and q["rows"] > 0
    rows = [l.split("\t") for l in
            open(Path(q["file"])).read().splitlines()]
    assert rows and all(len(r[1]) == 90 and len(r[2]) == 90 for r in rows)
    assert out["panel"]["query_library_is_not_a_member"] is True
    assert out["panel"]["subject_depth_pairs"] == 2_000_000
    staged = out["panel"]["staged_for_reuse"]
    assert staged["results_prefix"].startswith("s3://")
    assert "re-downloads" in staged["what_is_NOT_staged"], (
        "what is not staged is a cost and has to be stated, or a silent re-download "
        "reads as reuse")
    assert "depth" in " ".join(out["rules_enforced"]).lower()
    assert out["stages"]["report"]["state"] == "not selected"


@pytest.mark.skipif(not PANEL.exists(), reason="the panel table is not in this tree")
def test_the_scannable_subject_count_matches_the_committed_pass():
    """530 of 536 rows, and the six left out are left out rather than scanned shallower."""
    with open(PANEL) as fh:
        members = {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert len(members) == 536
    assert len(crossplant_stage.eligible_subjects(members, 2_000_000)) == 530
