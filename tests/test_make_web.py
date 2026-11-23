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
    del d["inputs"]["spectrum"]
    with pytest.raises(SystemExit, match="spectrum is not declared"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_declared_path_that_is_absent_is_not_read_as_not_run(tmp_path):
    d = _holdout()
    d["inputs"]["spectrum"] = {"path": "results/nowhere/divergent_known_s1.json"}
    with pytest.raises(SystemExit, match="does not exist"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_descriptor_cannot_dress_another_samples_run(tmp_path):
    d = _holdout()
    d["sample"]["accession"] = "SRR23973193"
    with pytest.raises(SystemExit, match="not in the run's own input path"):
        make_web.load_page(_write(tmp_path, d), None)


def test_the_report_has_no_ground_truth_section_on_any_sample():
    """Imo, 2026-11-23: "BY DEFAULT WE KNOW NOTHING. WE ALWAYS KNOW NOTHING."

    On unknown wastewater there is never a genome whose distance was fixed before the
    run, so the production report has no section for one - not a filled one on the
    fixture library, and not an empty one everywhere else. Known-distance material
    validates the tool and lives in tests/ and results/holdout_run/. This asserts the
    generator carries no input and no panel that only a fixture could fill, because the
    page went both wrong ways once already.
    """
    src = (ROOT / "experiments/make_web.py").read_text()
    for gone in ("truth_gene_map", "holdout_gene_map", "spectrum_reference_point",
                 "spectrum_reference_truth", "Genomes with a distance fixed"):
        assert gone not in src.replace("NO GROUND-TRUTH SECTION", ""), gone
    assert not [k for k in make_web.PAGE_INPUTS if "truth" in k]


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


def test_a_sample_nobody_has_written_about_gets_its_sentences_derived(tmp_path):
    """SPEC 13.4's operational test, one layer in: no hand-written English per sample.

    The descriptor used to carry the page title and the one line saying what the water
    was. That made a sentence a prerequisite for a page, and it was the only
    unattributable English on a page of numbers that all carry provenance. Both now come
    from the library's own facts.
    """
    import page_descriptor
    d = _holdout()
    d.pop("title", None)
    d["sample"].pop("lede", None)
    page = make_web.load_page(_write(tmp_path, d), None)
    want = page_descriptor.derive("SRR23421544")
    assert page["title"] == want["title"]
    assert page["sample"]["lede"] == want["lede"]
    assert page["overrides"] == []
    # And the facts behind the sentence are the library's, with their own basis recorded.
    assert want["facts"]["collected"] == "2022-05-10"
    assert "PRJNA917124" in want["facts"]["basis"]


def test_a_fact_the_library_does_not_carry_is_left_out_not_guessed(tmp_path):
    import page_descriptor
    # Berlin has no collection date in libraries.json, so no clause about one appears.
    assert "collected" not in page_descriptor.lede("SRR23973193")
    # A library with no sample block at all says so rather than describing nothing.
    assert "no sample block" in page_descriptor.lede("SRR23973265") or \
        page_descriptor.sample_facts("SRR23973265")


def test_a_hand_written_sentence_is_recorded_as_an_override(tmp_path):
    d = _holdout()
    d["sample"]["lede"] = "Something somebody typed."
    page = make_web.load_page(_write(tmp_path, d), None)
    assert page["sample"]["lede"] == "Something somebody typed."
    assert page["overrides"] == ["sample.lede"]
    assert page["derived"]["lede"].startswith("Treated wastewater")


def test_a_descriptor_for_a_library_nobody_has_described_refuses(tmp_path):
    d = _holdout()
    d["sample"]["accession"] = "SRR99999999"
    with pytest.raises(SystemExit, match="not in experiments/libraries.json"):
        make_web.load_page(_write(tmp_path, d), None)
