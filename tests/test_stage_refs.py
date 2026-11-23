"""The per-library resolver's three answers, and that it raises where a default would fit.

WHY THESE ARE ASSERTED RATHER THAN ASSUMED. The bug family this repo keeps finding is a
default where a refusal belongs, and the three answers here are easy to collapse into two:
an unknown library and a library with no panel are different facts, and a stage that treated
them alike would either score the wrong sample's genomes or report a measured absence as a
configuration error. The committed reference file is also asserted to be internally
consistent, because a null with no reason beside it is indistinguishable from a forgotten
entry.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import library_facts  # noqa: E402
import stage_refs  # noqa: E402


def test_unknown_library_raises_and_names_the_known_ones():
    """And it raises library_facts' exception, not a second one of this module's own."""
    with pytest.raises(stage_refs.UnknownLibrary) as e:
        stage_refs.resolve("SRR00000000", "truth_panel")
    assert "SRR23421544" in str(e.value)
    assert stage_refs.UnknownLibrary is library_facts.UnknownLibrary, (
        "one exception type for an unknown library, or a caller has to catch two")


def test_a_resolved_panel_carries_its_files_and_its_record():
    r = stage_refs.resolve("SRR23421544", "truth_panel")
    assert r["state"] == "resolved" and r["panel"] == "depositor_25"
    assert set(r["files"]) == {"genomes", "depth_table", "divergence_table"}
    assert stage_refs.missing_files(r) == [], "a named reference file is absent from the tree"
    assert r["record"]["what_this_panel_cannot_test"]


def test_a_null_panel_is_a_measured_absence_with_a_reason():
    r = stage_refs.resolve("SRR23973193", "truth_panel")
    assert r["state"] == "absent"
    assert "divergence ladder" in r["reason"]


def test_every_null_in_the_committed_file_has_a_reason():
    """The pointers live in libraries.json now; the panel definitions stay here."""
    refs = stage_refs.load()
    assert "libraries" not in refs, (
        "the per-library pointers moved into libraries.json; two homes for one fact is "
        "how they drift apart")
    libs = json.loads((ROOT / "experiments" / "libraries.json").read_text())["libraries"]
    for name, entry in libs.items():
        if "truth_panel" not in entry and "crossplant_panel" not in entry:
            continue
        for kind in ("truth_panel", "crossplant_panel"):
            if entry.get(kind) is None:
                assert entry.get(f"{kind}_absent_because"), (
                    f"{name}'s {kind} is null with no reason, which cannot be told apart "
                    f"from a forgotten entry")
            else:
                table = {"truth_panel": "truth_panels",
                         "crossplant_panel": "crossplant_panels"}[kind]
                assert entry[kind] in refs[table], f"{name} names an undefined {kind}"


def test_a_null_without_a_reason_raises(tmp_path, monkeypatch):
    libs = json.loads((ROOT / "experiments" / "libraries.json").read_text())
    libs["libraries"]["SRR23973193"].pop("truth_panel_absent_because")
    p = tmp_path / "libraries.json"
    p.write_text(json.dumps(libs))
    monkeypatch.setattr(library_facts, "FACTS", p)
    with pytest.raises(stage_refs.MissingReference):
        stage_refs.resolve("SRR23973193", "truth_panel")


def test_a_library_that_says_nothing_about_a_panel_is_not_an_absence(tmp_path, monkeypatch):
    """An incomplete entry and a measured absence are different facts.

    Answering "absent" for a library nobody has considered would report an unmeasured gap
    as a measured one, and the stage would then write a skip reason nobody wrote.
    """
    libs = json.loads((ROOT / "experiments" / "libraries.json").read_text())
    entry = libs["libraries"]["SRR23973193"]
    entry.pop("crossplant_panel", None)
    p = tmp_path / "libraries.json"
    p.write_text(json.dumps(libs))
    monkeypatch.setattr(library_facts, "FACTS", p)
    with pytest.raises(stage_refs.MissingReference) as e:
        stage_refs.resolve("SRR23973193", "crossplant_panel")
    assert "no crossplant_panel key" in str(e.value)


def test_the_band_cannot_be_taken_without_its_qualification():
    b = stage_refs.band()
    assert b["coverage"] == [0.25, 0.5] and b["identity_percent"] == [60, 70]
    assert b["chosen_after_seeing_the_positive"] is True
    assert b["measured_on"] == "SRR23421544"
    assert "93.46th percentile" in b["what_must_not_be_done_instead"]


def test_a_missing_reference_file_is_a_refusal_not_an_empty_resolution(tmp_path):
    refs = json.loads((ROOT / "experiments" / "stage_references.json").read_text())
    refs["truth_panels"]["depositor_25"]["genomes"] = "results/nowhere/absent.fa"
    p = tmp_path / "refs.json"
    p.write_text(json.dumps(refs))
    r = stage_refs.resolve("SRR23421544", "truth_panel", stage_refs.load(p))
    assert any("absent.fa" in m for m in stage_refs.missing_files(r))
