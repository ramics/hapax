"""Tests for the playground's controls: SPEC 13.3, one test per rule it enforces.

Run with `python3 -m pytest tests/test_web_controls.py`. The run directory used is the
committed holdout run, because it is the one in the repository: a test that needed Berlin's
run directory would be a test nobody can execute.

Each test is a way a control could manufacture a reading. A range wider than the
measurement; a shipped marker on a point the run never used; a smooth slider over a
threshold that was sampled at four values; a control that quietly disappears on a run that
did not measure it. The last is the shape this project has paid for twice, so absent and
unmeasured are asserted to be different states from refused.
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
import web_controls  # noqa: E402

HOLDOUT_PAGE = ROOT / "results/web/pages/holdout_SRR23421544.json"
RUN = ROOT / "results/holdout_cascade/run_s1"


def _bundle():
    page = make_web.load_page(HOLDOUT_PAGE, None)
    return make_web.bundle(page["_rundir"], 40, 40, 1, page)


def _declared():
    return json.loads(web_controls.CONTROLS_FILE.read_text())


def _resolve_with(tmp_path, doc, monkeypatch):
    f = tmp_path / "controls.json"
    f.write_text(json.dumps(doc))
    monkeypatch.setattr(web_controls, "CONTROLS_FILE", f)
    page = make_web.load_page(HOLDOUT_PAGE, None)
    return make_web.bundle(page["_rundir"], 40, 40, 1, page)


def test_every_control_resolves_to_a_range_or_a_stated_absence():
    # Neither is not an option: a control with no points and no bounds and no reason is a
    # widget whose range the page would have to invent.
    for c in _bundle()["controls"]["controls"]:
        assert (c.get("absent") or c.get("points") or c.get("bounds")
                or c.get("points_from_control")), c["id"]


def test_a_discrete_controls_points_are_the_runs_own_measured_points():
    b = _bundle()
    grid = {c["id"]: c for c in b["controls"]["controls"]}["sweep_tier"]["grid"]
    tt = json.loads((RUN / "tier_table.json").read_text())
    for tier in (c for c in tt["curve"] if c.get("threshold_sweep")):
        name, sweep = next(iter(tier["threshold_sweep"].items()))
        got = grid[tier["tier"]]
        assert [p["value"] for p in got["points"]] == list(sweep)
        assert [p["removed"] for p in got["points"]] == list(sweep.values())
        # And the shipped marker is the point the run actually ran at, not the first one.
        op = tt["operating_point"][web_controls.SWEEP_OPERATING_FIELD[name]]
        assert abs(float(got["shipped"]) - float(op)) < 1e-9


def test_kraken_confidence_is_four_points_and_not_a_smooth_axis():
    # The rule's own case: a factor of twelve in removal lives between the measured
    # points, so the control snaps to them. If a run ever sweeps confidence more finely
    # this test should be updated from the run, not relaxed.
    b = _bundle()
    grid = {c["id"]: c for c in b["controls"]["controls"]}["sweep_tier"]["grid"]
    k2 = [g for t, g in grid.items() if g["quantity"] == "confidence"]
    assert k2, "this run consulted no Kraken2 database"
    for g in k2:
        assert len(g["points"]) <= 5
        assert g["shipped"] in [p["value"] for p in g["points"]]


def test_a_continuous_range_does_not_run_past_the_measurement():
    for c in _bundle()["controls"]["controls"]:
        b = c.get("bounds")
        if not b:
            continue
        assert b["min"] <= b["max"]
        assert b["n_values"] > 0
        # Widening happens only to reach the shipped default, and says so when it does.
        if not b["widened_to_include_shipped"] and c.get("shipped") is not None:
            assert b["min"] <= c["shipped"] <= b["max"], c["id"]


def test_a_control_this_run_cannot_measure_is_absent_and_names_what_it_needs():
    got = {c["id"]: c for c in _bundle()["controls"]["controls"]}
    # Per-pair tier membership is not emitted by any run yet, so the union control states
    # what it would take. Absent is a census entry; it is not a refusal and it is not zero.
    assert "cascade_membership.jsonl" in got["union_tiers"]["absent"]
    assert "order" in got["union_tiers"]["absent"]
    # And the holdout library has no cross-sample occupancy measured at all.
    assert got["min_occupancy"].get("absent")
    assert "bounds" not in got["min_occupancy"]


def test_a_default_off_the_measured_grid_is_a_refusal(tmp_path, monkeypatch):
    # The failure this rules out: the page marks a setting as shipped that the run never
    # used, and a reader cites it.
    monkeypatch.setattr(web_controls, "SWEEP_OPERATING_FIELD",
                        dict(web_controls.SWEEP_OPERATING_FIELD,
                             confidence_sweep="diamond_bitscore"))
    with pytest.raises(SystemExit, match="not one of them"):
        _bundle()


def test_a_sweep_with_no_operating_point_field_is_a_refusal(monkeypatch):
    monkeypatch.setattr(web_controls, "SWEEP_OPERATING_FIELD",
                        {k: v for k, v in web_controls.SWEEP_OPERATING_FIELD.items()
                         if k != "mlen_sweep"})
    with pytest.raises(SystemExit, match="cannot tie to an operating point"):
        _bundle()


def test_a_control_declared_twice_is_a_refusal(tmp_path, monkeypatch):
    doc = _declared()
    doc["controls"].append(dict(doc["controls"][0]))
    with pytest.raises(SystemExit, match="declared twice"):
        _resolve_with(tmp_path, doc, monkeypatch)


def test_a_control_with_no_stated_arms_is_a_refusal(tmp_path, monkeypatch):
    doc = _declared()
    doc["controls"][0]["arms"] = {"control": "the identified reads"}
    with pytest.raises(SystemExit, match="arms.null is not stated"):
        _resolve_with(tmp_path, doc, monkeypatch)


def test_a_control_no_panel_reads_is_a_refusal(tmp_path, monkeypatch):
    doc = _declared()
    doc["controls"][0]["panels"] = []
    with pytest.raises(SystemExit, match="no panel reads it"):
        _resolve_with(tmp_path, doc, monkeypatch)


def test_a_control_with_no_resolver_is_a_refusal_not_a_blank_slider(tmp_path,
                                                                    monkeypatch):
    doc = _declared()
    doc["controls"].append({
        "id": "invented_floor", "label": "Invented floor", "kind": "continuous",
        "unit": "bits", "placement": "detail", "panels": ["significance"],
        "arms": {"control": "none", "null": "none"}})
    with pytest.raises(SystemExit, match="no resolver"):
        _resolve_with(tmp_path, doc, monkeypatch)


def test_the_cannot_move_table_is_not_allowed_to_be_empty(tmp_path, monkeypatch):
    doc = _declared()
    doc["fixed"] = []
    with pytest.raises(SystemExit, match="cannot-move table is empty"):
        _resolve_with(tmp_path, doc, monkeypatch)


def test_the_built_page_renders_every_control_it_resolved(tmp_path):
    out = tmp_path / "holdout.html"
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/make_web.py"),
                        str(HOLDOUT_PAGE), "--out", str(out)],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    t = out.read_text()
    b = json.loads(t[t.index('type="application/json">')
                     + len('type="application/json">'):t.index("</script>")])
    for c in b["controls"]["controls"]:
        assert f"data-ctl=\"{c['id']}\"" in t or f"'{c['id']}'" in t, c["id"]
    # The render check is the one that knows whether they reached the DOM, so its own
    # count is asserted here rather than re-derived - where it ran. Without a browser the
    # build says so and does not fake it, and this test asserts that rather than the
    # count, because a test that passed on a skipped check would be a test of nothing.
    if "render not checked" in r.stdout:
        assert "node or playwright is not available" in r.stdout
    else:
        assert "controls over" in r.stdout


def test_the_shipped_default_of_every_control_is_the_runs_operating_point():
    b = _bundle()
    got = {c["id"]: c for c in b["controls"]["controls"]}
    assert got["evalue_floor"]["shipped"] == b["significant_evalue"]


def test_a_run_with_fields_dropped_by_name_still_builds(tmp_path):
    """A field dropped on purpose carries its reason, and the page reads around it.

    `--report-without` is how a run says a field was not searched, which is a different
    state from searched-and-found-nothing and has to stay different. Every per-read
    reader in the generator used to assume a mapping, so one run with a field dropped by
    name would have taken the whole build down - and the reason belongs in the census the
    page prints, not in a traceback.
    """
    import shutil
    run = tmp_path / "run"
    shutil.copytree(RUN, run)
    rows = []
    for line in (RUN / "reads.jsonl").read_text().splitlines():
        d = json.loads(line)
        for f in ("cross_sample", "library_abundance", "motif", "chimaera"):
            d[f] = f"not searched: --report-without {f}"
        rows.append(json.dumps(d))
    (run / "reads.jsonl").write_text("\n".join(rows) + "\n")
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/make_web.py"),
                        str(HOLDOUT_PAGE), "--rundir", str(run),
                        "--out", str(tmp_path / "p.html")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stdout + r.stderr
    # And the control over a field that was not searched states that in the run's own
    # words, rather than offering a cut over numbers that do not exist. hapax writes the
    # reason as a plain string into the field where the value would be.
    page = make_web.load_page(HOLDOUT_PAGE, run)
    b = make_web.bundle(run, 20, 20, 1, page)
    assert b["not_searched"]["cross_sample"].startswith("not searched")
    got = {c["id"]: c for c in b["controls"]["controls"]}
    assert "--report-without cross_sample" in got["min_occupancy"]["absent"]
