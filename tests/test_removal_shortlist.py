"""The shortlist stage must reproduce the committed cell, and refuse where a default would fit.

WHY THE FIXTURE IS THE COMMITTED ROWS AND NOT A SYNTHETIC PILE. The only thing worth
asserting about this stage is that promoting the measurement into code did not move the
measurement: `results/holdout_cascade/contig_rows_s1.tsv` is the pile the published cell was
counted over, so the test converts those rows into the `removal_evidence.tsv` shape the
stage reads and asserts the published figures come back - 7,601 scored, 69 in the cell,
0.91% of the pile, and zero of the leak arm's 34 known contigs inside it. If a refactor
moves any of them, the published number and the stage have diverged and this fails.

The refusal tests are the other half, because the bug family here is a default where a
refusal belongs: an absent null, a pile with no protein route, an empty table.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "removal_shortlist.py"
ROWS = ROOT / "results" / "holdout_cascade" / "contig_rows_s1.tsv"

EV_COLS = ["contig", "contig_nt", "route", "index", "subject_or_taxon", "pident",
           "aln_aa", "evalue", "bitscore", "query_coverage", "matching_bases",
           "confidence", "tier"]


def write_evidence(path: Path, rows: list[dict]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(EV_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in EV_COLS) + "\n")


def fixture_rows(arms=("background",), tier="identified") -> list[dict]:
    """The committed holdout rows, in the evidence shape, protein level.

    `arm` and `tier` select exactly the pile the report's section 4b counted: background
    contigs the cascade identified, with a significant protein hit against the held-out
    index. The leak arm is the same conversion with arm="leak".
    """
    out = []
    with open(ROWS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["arm"] not in arms:
                continue
            if tier is not None and r["tier"] != tier:
                continue
            if r["holdout_prot_significant"] != "True":
                continue
            if not r["holdout_prot_identity"] or not r["holdout_prot_coverage"]:
                continue
            out.append({"contig": r["contig"], "contig_nt": r["length"],
                        "route": "protein", "index": "holdout_viral_prot",
                        "subject_or_taxon": "withheld_index_subject",
                        "pident": r["holdout_prot_identity"],
                        "aln_aa": "", "evalue": r["holdout_prot_evalue"],
                        "bitscore": r["holdout_prot_bitscore"],
                        "query_coverage": r["holdout_prot_coverage"],
                        "tier": "identified_norrna"})
    return out


def run(tmp_path: Path, evidence: Path, *extra: str) -> tuple[int, dict, str]:
    out = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, "-I", str(SCRIPT), "--evidence", str(evidence),
         "--out", str(out), *extra],
        cwd=ROOT, capture_output=True, text=True)
    payload = {}
    if (out / "removal_shortlist.json").exists():
        payload = json.loads((out / "removal_shortlist.json").read_text())
    return proc.returncode, payload, proc.stderr + proc.stdout


@pytest.mark.skipif(not ROWS.exists(), reason="the holdout rows are not in this tree")
def test_reproduces_the_published_cell(tmp_path):
    ev = tmp_path / "removal_evidence.tsv"
    rows = fixture_rows()
    write_evidence(ev, rows)
    code, out, log = run(tmp_path, ev, "--without-null")
    assert code == 0, log
    # The published pile and cell: section 4b of results/holdout_cascade/REPORT.md.
    assert out["pile"]["contigs_scored"] == 7601
    assert out["band"]["contigs"] == 69
    assert out["band"]["percent_of_scored_pile"] == pytest.approx(0.91, abs=0.005)
    assert out["pile"]["identity_median_percent"] == pytest.approx(40.0, abs=0.05)
    # The band is the shipped one and says so, with the sentence that qualifies it.
    assert out["operating_point"]["band_is_shipped_default"] is True
    assert out["operating_point"]["band_provenance"][
        "chosen_after_seeing_the_positive"] is True


@pytest.mark.skipif(not ROWS.exists(), reason="the holdout rows are not in this tree")
def test_the_candidate_lands_in_the_band_and_is_not_its_top_row(tmp_path):
    """The whole point: the withheld genome is ON the list, and not at the top of it."""
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, fixture_rows(arms=("background", "candidate"), tier=None))
    out_dir = tmp_path / "out"
    code, out, log = run(tmp_path, ev, "--without-null")
    assert code == 0, log
    listed = list(csv.DictReader(open(out_dir / "removal_shortlist.tsv"), delimiter="\t"))
    cand = [r for r in listed if "NODE_18516" in r["contig"]]
    assert len(cand) == 1, "the withheld genome is not on the shortlist"
    assert float(cand[0]["explains_fraction"]) == pytest.approx(0.4624, abs=1e-4)
    assert float(cand[0]["identity_percent"]) == pytest.approx(66.5, abs=1e-6)
    by_bits = sorted(listed, key=lambda r: -float(r["bitscore"]))
    assert "NODE_18516" not in by_bits[0]["contig"], (
        "the candidate ranks fourth by bitscore inside its own cell at nucleotide level "
        "and second at protein; a stage whose top row is the candidate has changed the "
        "measurement, and the report's point is that the cell is defensible, not the rank")
    assert out["band"]["ordering_is_validated"] is False


@pytest.mark.skipif(not ROWS.exists(), reason="the holdout rows are not in this tree")
def test_known_material_from_the_same_water_does_not_enter_the_band(tmp_path):
    """The direction check. The leak arm is known sequence from this library."""
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, fixture_rows(arms=("background", "leak"), tier=None))
    code, out, log = run(tmp_path, ev, "--without-null", "--near-reference-identity", "80")
    assert code == 0, log
    leak = fixture_rows(arms=("leak",), tier=None)
    assert len(leak) == 34
    assert all(float(r["pident"]) >= 80 for r in leak), (
        "no leak contig sits below 80% identity at protein level; if one does, the arm "
        "has changed and the band's direction is no longer pinned")
    assert out["control_arm"]["contigs_inside_the_band"] == 0


@pytest.mark.skipif(not ROWS.exists(), reason="the holdout rows are not in this tree")
def test_moving_the_band_moves_the_count_and_is_recorded_as_moved(tmp_path):
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, fixture_rows())
    code, out, log = run(tmp_path, ev, "--without-null",
                         "--band-coverage", "0.0", "1.01",
                         "--band-identity", "0", "100.01")
    assert code == 0, log
    assert out["band"]["contigs"] == 7601
    assert out["operating_point"]["band_is_shipped_default"] is False


def test_refuses_without_a_null(tmp_path):
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, [{"contig": "c1", "contig_nt": 1000, "route": "protein",
                         "index": "viral", "subject_or_taxon": "s", "pident": 66.5,
                         "evalue": 1e-99, "bitscore": 800, "query_coverage": 0.46,
                         "tier": "identified_norrna"}])
    code, _, log = run(tmp_path, ev)
    assert code == 2, log
    assert "without-null" in log


def test_refuses_a_pile_with_no_protein_route(tmp_path):
    """A nucleotide-only pile is a refusal, not an empty shortlist."""
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, [{"contig": "c1", "contig_nt": 900, "route": "kraken2",
                         "index": "k2_pluspf16", "subject_or_taxon": "10239",
                         "confidence": 0.1, "tier": "identified_norrna"}])
    code, _, log = run(tmp_path, ev, "--without-null")
    assert code != 0
    assert "no protein-route rows" in log


def test_refuses_an_empty_table(tmp_path):
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, [])
    code, _, log = run(tmp_path, ev, "--without-null")
    assert code != 0
    assert "no rows" in log or "empty" in log


def test_a_protein_row_with_no_coverage_is_counted_not_floored(tmp_path):
    """Blank coverage must not read as 0.0, which would put the row in the lowest cell."""
    ev = tmp_path / "removal_evidence.tsv"
    write_evidence(ev, [
        {"contig": "c1", "contig_nt": 1000, "route": "protein", "index": "viral",
         "subject_or_taxon": "s", "pident": 66.5, "evalue": 1e-99, "bitscore": 800,
         "query_coverage": 0.46, "tier": "identified_norrna"},
        {"contig": "c2", "contig_nt": 1000, "route": "protein", "index": "viral",
         "subject_or_taxon": "s", "pident": 65.0, "evalue": 1e-99, "bitscore": 700,
         "query_coverage": "", "tier": "identified_norrna"}])
    code, out, log = run(tmp_path, ev, "--without-null")
    assert code == 0, log
    assert out["census"]["protein_rows_without_coverage"] == 1
    assert out["pile"]["contigs_scored"] == 1
    assert all(c["contigs"] == 0 for c in out["joint_table"]
               if c["coverage_hi"] <= 0.25)
