"""Tests for the one-generator rule: every sample's page comes from the same code.

Run with `python3 -m pytest tests/test_make_web.py`. Each test is a way a page has
already said, or could say, something about the wrong sample: a spectrum read from repo
paths whatever run was passed, Berlin's library row printed on the holdout page, and a
ground truth that a sample could acquire because a file of the right name existed.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import make_web  # noqa: E402

HOLDOUT = ROOT / "results/web/pages/holdout_SRR23421544.json"
BERLIN = ROOT / "results/web/pages/berlin_SRR23973193.json"


def _write(tmp_path, d) -> Path:
    p = tmp_path / "page.json"
    p.write_text(json.dumps(d))
    return p


def _holdout() -> dict:
    return json.loads(HOLDOUT.read_text())


def test_every_committed_descriptor_declares_every_input():
    for f in sorted((ROOT / "results/web/pages").glob("*.json")):
        d = json.loads(f.read_text())
        assert set(d["inputs"]) == set(make_web.PAGE_INPUTS), f


def test_an_undeclared_input_is_a_refusal_not_a_default(tmp_path):
    d = _holdout()
    del d["inputs"]["truth_gene_map"]
    with pytest.raises(SystemExit, match="truth_gene_map is not declared"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_declared_path_that_is_absent_is_not_read_as_not_run(tmp_path):
    d = _holdout()
    d["inputs"]["spectrum"] = {"path": "results/nowhere/divergent_known_s1.json"}
    with pytest.raises(SystemExit, match="does not exist"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_reference_point_needs_its_true_distance(tmp_path):
    d = json.loads(BERLIN.read_text())
    d["rundir"] = "results/holdout_cascade/run_s1"
    d["sample"]["accession"] = "SRR23421544"
    d["inputs"]["spectrum_reference_truth"] = {"not_run": "x"}
    with pytest.raises(SystemExit, match="come together"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_descriptor_cannot_dress_another_samples_run(tmp_path):
    d = _holdout()
    d["sample"]["accession"] = "SRR23973193"
    with pytest.raises(SystemExit, match="not in the run's own input path"):
        make_web.load_page(_write(tmp_path, d), None)


def test_no_ground_truth_unless_named_by_path(tmp_path):
    # The run directory has a gene map sitting next to it; declaring not_run must win.
    d = _holdout()
    d["inputs"]["truth_gene_map"] = {"not_run": "pretend this sample has no answer"}
    page = make_web.load_page(_write(tmp_path, d), None)
    assert page["_inputs"]["truth_gene_map"] is None


def test_a_row_measured_on_one_library_is_not_printed_on_another():
    berlin = {m["id"]: dict(m["rows"]) for m in make_web.resolve_methods("SRR23973193")}
    hold = {m["id"]: dict(m["rows"]) for m in make_web.resolve_methods("SRR23421544")}
    other = {m["id"]: dict(m["rows"]) for m in make_web.resolve_methods("SRR00000000")}
    assert "Berlin" in berlin["library"]["Input"]
    assert "SRR23973193" not in json.dumps(hold)
    assert hold["cascade"]["Bias"] == make_web.NOT_MEASURED
    # A procedure row is the tool's, and reads the same everywhere.
    assert berlin["cascade"]["Procedure"] == hold["cascade"]["Procedure"]
    # A sample nobody has written about still gets every entry, numbered the same.
    assert list(other) == list(berlin)
    assert other["library"]["Input"] == make_web.NOT_MEASURED


def test_every_entry_has_the_six_fixed_rows():
    for m in make_web.resolve_methods("SRR23973193"):
        assert len(m["rows"]) == 6, m["id"]


def test_the_holdout_page_builds_and_draws_no_berlin_spectrum(tmp_path):
    out = tmp_path / "holdout.html"
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/make_web.py"),
                        str(HOLDOUT), "--out", str(out)],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    t = out.read_text()
    a = t.index('type="application/json">') + len('type="application/json">')
    b = json.loads(t[a:t.index("</script>", a)])
    assert b["spectrum"] is None
    assert b["page"]["sample"]["accession"] == "SRR23421544"
    assert "Berlin" not in json.dumps(b["methods"])
    assert "Sewage from one treatment plant in Berlin" not in t


def test_a_page_that_throws_is_not_written(tmp_path, monkeypatch):
    # The {} absent-stage bug: one throw silently dropped every later panel. The build
    # now runs the page and refuses it.
    monkeypatch.setattr(make_web, "PAGE", make_web.PAGE.replace(
        "const METHODS = D.methods || [];", "const METHODS = D.methods || []; null.x;", 1))
    out = tmp_path / "throw.html"
    monkeypatch.setattr(sys, "argv", ["make_web.py", str(HOLDOUT), "--out", str(out)])
    if subprocess.run(["node", "-e", "require('playwright')"],
                      capture_output=True).returncode:
        pytest.skip("no playwright here, so the render check cannot run")
    with pytest.raises(SystemExit, match="threw while rendering"):
        make_web.main()
    assert not out.exists()


def _summary_titles() -> list[str]:
    """The one fixed section list, read out of the generator rather than retyped here."""
    src = (ROOT / "experiments/make_web.py").read_text()
    body = src[src.index("const SUMMARY = ["):]
    return [m for m in __import__("re").findall(r"'([^']+)'",
                                                body[:body.index("];")])]


def test_one_report_format_the_same_sections_in_the_same_order(tmp_path):
    """SPEC 13.4: a page never grows or drops a section because of what a sample has.

    The failure this rules out is a reader learning where a figure lives on one library
    and not finding it on the next - and, worse, a missing section and a measured
    absence looking identical. The holdout run is the thin one in the repository: it has
    no spectrum inputs at all, so if any section could vanish for want of data it would
    vanish here. The build prints the sections it drew, in order, and this asserts that
    line against the generator's single list.
    """
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/make_web.py"),
                        str(HOLDOUT), "--out", str(tmp_path / "p.html")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stdout + r.stderr
    line = [ln for ln in r.stdout.splitlines() if ln.startswith("summary:")]
    assert line, r.stdout
    drawn = [t.strip() for t in line[0][len("summary:"):].split("|")]
    assert drawn == _summary_titles()
