"""The truth-panel stage must compose the committed panel exactly, and refuse to guess one.

The composition is the part that can be checked for free: the present and absent arms come
out of `results/holdout_candidates/depth.tsv` by code, so a run of this stage on
SRR23421544 must return the published 25 present and 53 zero-read controls at the published
depth and divergence range. Scoring needs minimap2 and an assembly and is skipped where
neither is present; what is asserted here instead is that the stage refuses rather than
scoring a library whose panel it does not know, and that the sentences about what the panel
cannot test travel in its output rather than only in a report.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "truth_panel.py"
DEPTH = ROOT / "results" / "holdout_candidates" / "depth.tsv"


def run(tmp_path: Path, *extra: str) -> tuple[int, dict, str]:
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, "-I", str(SCRIPT), "--work",
                           str(tmp_path / "work"), "--out", str(out), *extra],
                          cwd=ROOT, capture_output=True, text=True)
    payload = {}
    if (out / "truth_panel.json").exists():
        payload = json.loads((out / "truth_panel.json").read_text())
    return proc.returncode, payload, proc.stdout + proc.stderr


@pytest.mark.skipif(not DEPTH.exists(), reason="the depth table is not in this tree")
def test_composes_the_committed_panel(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR23421544")
    assert code == 0, log
    c = out["composition"]
    assert (c["present"], c["absent_control"]) == (25, 53)
    assert c["mean_depth_range"] == [0.445, 700.665]
    assert c["divergence_percent_range"] == [0.07, 30.38]
    assert out["scoring"]["state"] == "not scored"


@pytest.mark.skipif(not DEPTH.exists(), reason="the depth table is not in this tree")
def test_what_the_panel_cannot_test_travels_with_the_number(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR23421544")
    assert code == 0, log
    cannot = " ".join(out["what_this_panel_cannot_test"])
    assert "3.8%" in cannot, "the sub-3.8% fusion limit must be stated where it is reported"
    assert "depth 1" in cannot
    assert isinstance(out["residue_fraction_from_the_same_pass"], str), (
        "with no residue fraction passed, the field must say so rather than be absent")
    assert "8.79" in out["residue_fraction_from_the_same_pass"]


@pytest.mark.skipif(not DEPTH.exists(), reason="the depth table is not in this tree")
def test_residue_fraction_travels_when_given(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR23421544",
                         "--residue-fraction", "0.3542")
    assert code == 0, log
    assert out["residue_fraction_from_the_same_pass"] == pytest.approx(0.3542)


def test_a_library_with_no_panel_is_skipped_with_its_reason(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR23973193")
    assert code == 0, log
    assert out["status"] == "skipped"
    assert "divergence ladder" in out["reason"]
    assert "found nothing" in out["not_the_same_as"]


def test_an_unknown_library_refuses_to_start(tmp_path):
    code, out, log = run(tmp_path, "--library", "SRR99999999")
    assert code == 2, log
    assert "is not in experiments/libraries.json" in log
    assert out == {}


@pytest.mark.skipif(not DEPTH.exists(), reason="the depth table is not in this tree")
def test_scoring_without_minimap2_refuses_rather_than_scoring_nothing(tmp_path):
    fa = tmp_path / "contigs.fa"
    fa.write_text(">c1\nACGTACGTAC\n")
    code, _, log = run(tmp_path, "--library", "SRR23421544", "--contigs", str(fa),
                       "--minimap2", str(tmp_path / "no_such_minimap2"))
    assert code == 2, log
    assert "minimap2" in log
