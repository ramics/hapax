"""The three new steps a brand-new library needs, and what each of them refuses.

Running production over a library nobody has measured adds exactly three things to the
ordering `holdout_cascade_run.sh` already had: a trim measured rather than inherited, a
trimmed copy for the assembler, and a header relabelling. Each has a refusal in it, and a
refusal that is not tested is a refusal that will be discovered by a launch.

What is NOT tested here is the numbers: nothing in this file runs MEGAHIT, Kraken2 or a
search. It tests which inputs make a step stop, in the shape `tests/test_library_facts.py`
uses - the negatives included, because the negatives are the design.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import library_facts  # noqa: E402
import stage_refs  # noqa: E402
import trim_from_prefix  # noqa: E402

PROD_LIBRARY = "SRR23702504"


def seed_row(mode=14, share=0.51, anchored=30_000, r1_mode=0):
    return {"libraries": {PROD_LIBRARY: {
        "run": PROD_LIBRARY, "note": "test",
        "anchored": {"r1": anchored, "r2": anchored},
        "r1_five_prime_clip": {"mode": r1_mode, "mode_share": 0.40, "n": anchored},
        "r2_five_prime_clip": {"mode": mode, "mode_share": share, "n": anchored},
    }}}


# ------------------------------------------------------------------ the trim decision
def test_agreeing_seeds_name_the_mode():
    d = trim_from_prefix.decide([seed_row() for _ in range(5)], PROD_LIBRARY, 0.25, 2_000)
    assert d["problems"] == []
    assert d["trim_r2"] == 14


def test_a_mode_that_moves_between_seeds_is_refused():
    rows = [seed_row(), seed_row(), seed_row(mode=9), seed_row(), seed_row()]
    d = trim_from_prefix.decide(rows, PROD_LIBRARY, 0.25, 2_000)
    assert d["trim_r2"] is None
    assert any("not the same in every seed" in p for p in d["problems"])


def test_a_flat_histogram_is_not_a_prefix():
    rows = [seed_row(share=0.04) for _ in range(5)]
    d = trim_from_prefix.decide(rows, PROD_LIBRARY, 0.25, 2_000)
    assert d["trim_r2"] is None
    assert any("flat histogram" in p for p in d["problems"])


def test_too_little_anchored_describes_the_reference_not_the_library():
    rows = [seed_row(anchored=200) for _ in range(5)]
    d = trim_from_prefix.decide(rows, PROD_LIBRARY, 0.25, 2_000)
    assert d["trim_r2"] is None
    assert any("anchored" in p for p in d["problems"])


def test_both_mates_spiking_at_the_same_length_is_an_instrument_property():
    rows = [seed_row(r1_mode=14) for _ in range(5)]
    d = trim_from_prefix.decide(rows, PROD_LIBRARY, 0.25, 2_000)
    assert d["trim_r2"] is None
    assert any("instrument property" in p for p in d["problems"])


def test_a_library_missing_from_the_seed_files_is_refused():
    d = trim_from_prefix.decide([seed_row() for _ in range(5)], "SRR00000000", 0.25, 2_000)
    assert d["trim_r2"] is None
    assert any("no row for" in p for p in d["problems"])


def test_the_cli_prints_the_integer_and_exits_two_on_a_refusal(tmp_path):
    good = [tmp_path / f"nontemplate_prefix_s{i}.json" for i in (1, 2)]
    for p in good:
        p.write_text(json.dumps(seed_row()))
    cmd = [sys.executable, "-I", str(ROOT / "experiments/trim_from_prefix.py"),
           "--library", PROD_LIBRARY, "--prefix-json", *[str(p) for p in good],
           "--out", str(tmp_path / "decision.json")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "14"               # a shell reads this and nothing else

    (tmp_path / "nontemplate_prefix_s2.json").write_text(json.dumps(seed_row(mode=7)))
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 2
    assert "REFUSING" in r.stderr
    assert json.loads((tmp_path / "decision.json").read_text())["trim_r2"] is None


def test_one_seed_cannot_show_agreement(tmp_path):
    p = tmp_path / "nontemplate_prefix_s1.json"
    p.write_text(json.dumps(seed_row()))
    r = subprocess.run([sys.executable, "-I",
                        str(ROOT / "experiments/trim_from_prefix.py"),
                        "--library", PROD_LIBRARY, "--prefix-json", str(p)],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "every seed" in r.stderr


# ------------------------------------------------------------------ the trimmed copy
def fastq(path: Path, seqs: list[str]) -> None:
    with open(path, "w") as h:
        for i, s in enumerate(seqs):
            h.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


def run_trim(tmp_path, trim, r2_seqs=None):
    r1, r2 = tmp_path / "x_1.fq", tmp_path / "x_2.fq"
    fastq(r1, ["A" * 20, "C" * 20])
    fastq(r2, r2_seqs or ["G" * 20, "T" * 20])
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/trim_mate2.py"),
                        "--r1", str(r1), "--r2", str(r2), "--trim-r2", str(trim),
                        "--out-dir", str(tmp_path / "out"), "--prefix", "x"],
                       capture_output=True, text=True)
    return r, tmp_path / "out/x_1.fq", tmp_path / "out/x_2.fq"


def test_only_mate_two_is_cut(tmp_path):
    r, o1, o2 = run_trim(tmp_path, 5)
    assert r.returncode == 0, r.stderr
    assert [len(l.strip()) for l in o1.read_text().splitlines()[1::4]] == [20, 20]
    assert [len(l.strip()) for l in o2.read_text().splitlines()[1::4]] == [15, 15]
    assert o2.read_text().splitlines()[1] == "G" * 15


def test_zero_is_an_instruction_not_an_absence(tmp_path):
    r, o1, o2 = run_trim(tmp_path, 0)
    assert r.returncode == 0, r.stderr
    assert o2.read_text().splitlines()[1] == "G" * 20
    man = json.loads((tmp_path / "out/x_trim_manifest.json").read_text())
    assert man["trim_r2"] == 0


def test_a_mate_shorter_than_the_prefix_keeps_the_record_count(tmp_path):
    r, o1, o2 = run_trim(tmp_path, 15, r2_seqs=["G" * 20, "T" * 10])
    assert r.returncode == 0, r.stderr
    assert len(o1.read_text().splitlines()) == len(o2.read_text().splitlines())
    man = json.loads((tmp_path / "out/x_trim_manifest.json").read_text())
    assert man["mate2_shorter_than_the_prefix"] == 1


def test_an_unknown_library_is_refused_rather_than_given_berlins_trim(tmp_path):
    r1, r2 = tmp_path / "SRR00000000_1.fq", tmp_path / "SRR00000000_2.fq"
    fastq(r1, ["A" * 20])
    fastq(r2, ["G" * 20])
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/trim_mate2.py"),
                        "--r1", str(r1), "--r2", str(r2), "--library", "SRR00000000",
                        "--out-dir", str(tmp_path / "out"), "--prefix", "x"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "SRR00000000" in r.stderr and "libraries.json" in r.stderr


# ------------------------------------------------------------------ the relabelling
def test_megahit_headers_become_parsable_and_a_foreign_header_is_refused(tmp_path):
    src = tmp_path / "final.contigs.fa"
    src.write_text(">k141_3 flag=1 multi=7.2500 len=486\nACGT\n")
    dst = tmp_path / "renamed.fa"
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/megahit_rename.py"),
                        "--contigs", str(src), "--out", str(dst), "--tag", "SRRTEST"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    head = dst.read_text().splitlines()[0]
    assert head.startswith(">SRRTEST_NODE_1_length_486_cov_7.2500")
    assert "k141_3" in head                      # traceable to the assembly it came from

    src.write_text(">NODE_1_length_486_cov_7.25\nACGT\n")
    r = subprocess.run([sys.executable, "-I", str(ROOT / "experiments/megahit_rename.py"),
                        "--contigs", str(src), "--out", str(dst), "--tag", "SRRTEST"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "not a MEGAHIT contig name" in r.stderr


# ------------------------------------------------- the new library's committed entries
def test_the_prod_library_has_its_facts_and_its_trim_is_flagged_unverified():
    """0 is a stated value on a library the ribosomal instrument could not measure."""
    e = library_facts.library(PROD_LIBRARY)
    assert e["layout"] == "paired"
    assert e["pairs"] == 32_125_614
    assert e["trim_r2"] == 0
    assert e["trim_r2_unverified"] is True
    assert "0.13-0.14% of pairs" in e["trim_r2_basis"]
    block = library_facts.fact(PROD_LIBRARY, "block")
    assert 4 * block + 200_000 <= e["pairs"]      # five seeds fit, which is the whole point


def test_provenance_can_name_the_commit_on_a_machine_with_no_checkout(monkeypatch):
    """An instance has a tarball and no git, and a manifest must still say which code ran."""
    import provenance as pv
    monkeypatch.setattr(pv, "_git", lambda *a: None)
    monkeypatch.delenv(pv.SHA_ENV, raising=False)
    with pytest.raises(pv.ProvenanceError):
        pv.git_sha()
    monkeypatch.setenv(pv.SHA_ENV, "abc123")
    assert pv.git_sha() == "abc123"
    assert pv.tree_dirty() is False
    # A machine WITH a checkout ignores it, so it cannot dress a dirty tree as clean.
    monkeypatch.undo()
    monkeypatch.setenv(pv.SHA_ENV, "deadbeef")
    assert pv.git_sha() != "deadbeef"


def test_the_prefix_measurement_resolves_an_accession_through_the_committed_entry():
    import nontemplate_prefix_map as nt
    spec = nt.spec_for(PROD_LIBRARY)
    assert spec["run"] == PROD_LIBRARY
    assert spec["r1"].endswith(f"{PROD_LIBRARY}_1.fastq.gz")
    assert nt.spec_for("rna_2021")["run"] == "SRR23973193"        # the nicknames still work
    with pytest.raises(library_facts.UnknownLibrary):
        nt.spec_for("SRR00000000")


def test_truth_is_absent_with_a_reason_and_the_crossplant_panel_is_legal():
    truth = stage_refs.resolve(PROD_LIBRARY, "truth_panel")
    assert truth["state"] == "absent" and "No depositor-assembled" in truth["reason"]
    cp = stage_refs.resolve(PROD_LIBRARY, "crossplant_panel")
    assert cp["state"] == "resolved" and cp["panel"] == "cp3"
    assert stage_refs.missing_files(cp) == []


def test_the_query_library_is_not_a_member_of_the_panel_it_subtracts_against():
    """The one refusal that would invent recurrence: 54% false growth against 6% swapped."""
    panel = ROOT / "results/crossplant/panel.tsv"
    runs = {line.split("\t")[0] for line in panel.read_text().splitlines()[1:]}
    assert PROD_LIBRARY not in runs
    plants = {line.split("\t")[10].strip() for line in panel.read_text().splitlines()[1:]}
    # and the plant IS in it, or there is no within-plant stratum and the paired
    # per-plant comparison the measurement requires cannot be drawn.
    assert stage_refs.resolve_library(PROD_LIBRARY)["plant"] in plants
