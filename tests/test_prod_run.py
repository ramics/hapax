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


# ------------------------------------------------- the facts that come off the file
def write_fastq(path, n, length):
    with open(path, "w") as h:
        for i in range(n):
            h.write(f"@r{i}\n{'A' * length}\n+\n{'I' * length}\n")


def test_the_file_answers_count_length_layout_and_stride(tmp_path):
    """Three of the four facts a new library needed by hand are in the file itself."""
    import measure_library as ml
    r1, r2 = tmp_path / "SRRX_1.fq", tmp_path / "SRRX_2.fq"
    write_fastq(r1, 1_000, 151)
    write_fastq(r2, 1_000, 151)
    m = ml.measure(r1, r2, pairs_per_seed=100, seeds=5)
    assert m["pairs"] == 1_000
    assert m["read_length"] == [151, 151]
    assert m["layout"] == "paired"
    assert 4 * m["block"] + 100 <= 1_000
    # and what the file cannot answer is named rather than left out
    assert "trim_r2" in m["not_measurable_from_the_file"]


def test_unequal_mates_are_not_pairs(tmp_path):
    import measure_library as ml
    r1, r2 = tmp_path / "SRRX_1.fq", tmp_path / "SRRX_2.fq"
    write_fastq(r1, 10, 100)
    write_fastq(r2, 9, 100)
    with pytest.raises(ValueError, match="not pairs"):
        ml.measure(r1, r2, pairs_per_seed=2, seeds=5)


def test_a_stride_below_the_block_would_overlap_the_seeds_and_is_refused():
    """The arithmetic passes happily on overlapping blocks, so it is checked."""
    import measure_library as ml
    assert ml.stride_for(32_125_614, 200_000, 5) >= 200_000
    assert ml.stride_for(4_290_429, 200_000, 5) >= 200_000
    with pytest.raises(ValueError, match="DISJOINT"):
        ml.stride_for(400, 100, 5)


def test_a_single_end_library_is_recorded_as_single(tmp_path):
    import measure_library as ml
    r1 = tmp_path / "SRRX_1.fq"
    write_fastq(r1, 500, 145)
    m = ml.measure(r1, None, pairs_per_seed=50, seeds=5)
    assert m["layout"] == "single" and m["read_length"] == [145]


# ------------------------------------------------- validation is not production
def test_the_depositor_panel_and_the_spiked_ladder_are_not_production_stages():
    """Unknown water has no depositor, so a production run does not ask (SPEC 13.4)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("hapax_mod", ROOT / "hapax.py")
    hp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hp)
    assert set(hp.VALIDATION_STAGES) == {"panel", "calibration"}
    assert "panel" not in hp.PRODUCTION_STAGES
    assert "calibration" not in hp.PRODUCTION_STAGES
    # every other stage is still there, in the same order
    assert hp.PRODUCTION_STAGES == [s for s in hp.STAGE_ORDER
                                    if s not in hp.VALIDATION_STAGES]
    # and the production ordering a run actually uses does not select them
    run_sh = (ROOT / "experiments/prod_run.sh").read_text()
    stages_line = [l for l in run_sh.splitlines() if l.startswith("STAGES=")][0]
    assert "panel" not in stages_line.replace("crossplant", "")
    assert "calibration" not in stages_line


def test_the_stage_list_is_typed_once(tmp_path):
    """One list, in hapax.py. The first production pass ran `panel` because of a copy.

    `--stages` was typed in three places: hapax.py's default, prod_run.sh's `STAGES`, and
    the launcher's `--stages`. Two were changed when the validation stages came out of
    production and the third - the launcher's, the one that builds the command that
    actually runs - was not, so launch 5 ran `panel` in a production pass. It refused,
    correctly, and said why; the defect is that the disagreement was possible at all.
    """
    launcher = (ROOT / "experiments/ec2_prod.py").read_text()
    assert "PRODUCTION_STAGES" in launcher, "the launcher must import the list, not type it"
    for stage in ("contigs", "linkage", "shortlist"):
        assert f'"{stage}' not in launcher.replace("PRODUCTION_STAGES", ""), (
            f"{stage} is named literally in the launcher: a second stage list")
    # the shell passes --stages only when a caller overrode it
    run_sh = (ROOT / "experiments/prod_run.sh").read_text()
    assert "${STAGES:+--stages $STAGES}" in run_sh
    # and what the launcher would actually pass IS that list, read from the launcher
    # itself rather than from a copy of the expectation
    import importlib.util
    loaded = []
    for name, path in (("hapax_mod2", ROOT / "hapax.py"),
                       ("ec2_prod_mod", ROOT / "experiments/ec2_prod.py")):
        s = importlib.util.spec_from_file_location(name, path)
        m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        loaded.append(m)
    hp_stages, launcher = loaded[0].PRODUCTION_STAGES, loaded[1]
    assert launcher.PRODUCTION_STAGES == hp_stages, "the launcher's list is not hapax.py's"
    # what it actually ships is that list minus the stages the instance has no reference
    # data for, each of which has to say why rather than simply be absent
    omitted = set(hp_stages) - set(launcher.LAUNCHER_STAGES)
    assert omitted == set(launcher.NO_REFERENCE_ON_THIS_INSTANCE)
    for stage, reason in launcher.NO_REFERENCE_ON_THIS_INSTANCE.items():
        assert stage in hp_stages and len(reason) > 40, stage
    assert "panel" not in launcher.LAUNCHER_STAGES
    assert "calibration" not in launcher.LAUNCHER_STAGES


def test_a_library_with_no_entry_is_not_a_refusal_for_a_launcher():
    """resolve() must refuse for a stage and must not for a launcher reporting what will run."""
    with pytest.raises(stage_refs.UnknownLibrary):
        stage_refs.resolve("SRR00000000", "crossplant_panel")
    r = stage_refs.resolve_optional("SRR00000000", "crossplant_panel")
    assert r["state"] == "unregistered" and "SRR00000000" in r["reason"]


# ------------------------------------------------- the payload is the instance's whole repo
def payload_tree(tmp_path, library=PROD_LIBRARY):
    """A directory holding exactly what the launcher ships, and nothing else."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("ec2_prod", ROOT / "experiments/ec2_prod.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    tree = tmp_path / "payload"
    for rel in mod.payload_files(ROOT, library):
        src = ROOT / rel
        assert src.exists(), f"{rel} is in the payload list and not on disk"
        dst = tree / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(src.read_bytes())
    return tree


def test_every_preflighting_script_imports_from_the_payload_alone(tmp_path):
    """The one test that would have caught the first production run's refusal.

    A payload that is missing a module looks fine until a stage starts, which on this run
    was two hours and one assembly in. `residue_report.py` fails at IMPORT time there: it
    reads a constant out of `readlevel_features.py`, which imports `common` from the
    archived assembler tree, which the payload excluded. --help exercises exactly the
    module-level work that failed and nothing expensive.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("hapax_mod", ROOT / "hapax.py")
    hp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hp)
    tree = payload_tree(tmp_path)
    for script in sorted(hp.PREFLIGHT_SCRIPTS):
        r = subprocess.run([sys.executable, "-I", f"experiments/{script}", "--help"],
                           cwd=tree, capture_output=True, text=True)
        assert r.returncode == 0, f"{script} cannot start from the payload:\n{r.stderr}"


def test_the_reference_tables_a_stage_reads_are_in_the_payload(tmp_path):
    """Committed results read as reference data, which nothing under results/ shipped."""
    tree = payload_tree(tmp_path)
    cp = stage_refs.resolve(PROD_LIBRARY, "crossplant_panel")
    for path in cp["files"].values():
        rel = Path(path).resolve().relative_to(ROOT)
        assert (tree / rel).exists(), f"{rel} is read by the crossplant stage and not shipped"


def test_the_generated_tables_carry_numbers_and_not_a_column_of_dashes(tmp_path):
    """The table generator is checked for CELLS, not for exiting zero.

    THE DEFECT THIS PINS. `report_prod.py` read the cascade curve's four fields out of
    `tiers["tiers"][name]`, which holds each tier's own internal sweep, while the per-tier
    removal lives one row per tier in `tiers["curve"]`. So the one table that says what each
    index bought rendered as `| artefact | - | - | - | - |` for every tier, on both libraries
    it had ever been run on, and the page looked finished. A generated file is this repo's
    defence against prose drifting from its outputs, which makes a silently empty generated
    file worse than no file: an empty cell and a measurement that never happened are
    indistinguishable, which is the failure `record no significant hits as an explicit
    census` already exists for one layer in.

    So the assertion is that the tier table holds numbers. It runs against whatever
    production run directories are committed, and skips if none is.
    """
    import re
    runs = sorted(ROOT.glob("results/prod_*/run_s*/tiers.json"))
    if not runs:
        pytest.skip("no production run directory is committed")
    res = runs[0].parent.parent
    library = res.name.replace("prod_", "")
    out = tmp_path / "TABLES.md"
    r = subprocess.run([sys.executable, "-I", "experiments/report_prod.py",
                        "--library", library, "--results", str(res), "--out", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = out.read_text()
    tier_rows = [l for l in text.splitlines()
                 if l.startswith("| ") and re.search(r"\| (rrna|prot_sprot) \|", l)]
    assert tier_rows, "the cascade tier table has no rows at all"
    for row in tier_rows:
        cells = [c.strip() for c in row.strip("|").split("|")[1:4]]
        assert all(re.fullmatch(r"[\d.]+( \([\d.]+-[\d.]+\))?", c) for c in cells), (
            f"a cascade tier row carries no number: {row}")
    # the sections whose absence is invisible, each of which was missing on the first run
    for heading in ("## The assembly", "## The cascade", "## Measured off the file",
                    "## The mate-2 non-template prefix"):
        assert heading in text, f"{heading} did not render from {res}"


def test_the_host_summary_never_invites_a_pairing_of_two_maxima():
    """The closest row's identity must ship with that row's OWN alignment length.

    `host_confidence_summary` carried `best_identity_percent` and `best_aligned_aa` as
    maxima taken separately over the rows, and they are different subjects on every seed of
    the first production run. Printed as one match they read "100.0% identity over 53 amino
    acids", where the 100% row is 21 aa and the 53 aa row is something else - a subject that
    does not exist, flattering hardest where identity is least trustworthy. A web page did
    pair them. The invariant is checked on committed output rather than a fixture, because
    the defect was in what the field shape invited, not in any arithmetic.
    """
    rankings = sorted(ROOT.glob("results/prod_*/run_s*/ranking.json"))
    if not rankings:
        pytest.skip("no ranking is committed")
    for path in rankings:
        r = json.loads(path.read_text())
        h = r["host_confidence_summary"]
        rows = r["subjects"]
        if not rows:
            continue
        assert "closest_row_aligned_aa" in h, f"{path} has no closest_row_aligned_aa"
        assert "best_fields_are_independent" in h, f"{path} does not say the maxima are independent"
        closest = max(rows, key=lambda x: x["median_identity_percent"] or 0)
        assert h["best_identity_percent"] == closest["median_identity_percent"]
        assert h["closest_row_aligned_aa"] == closest["median_aligned_aa"], (
            f"{path}: closest_row_aligned_aa is not the closest row's own length")
        assert h["longest_alignment_aa"] == max((x["median_aligned_aa"] or 0) for x in rows)


def test_a_seed_directory_with_no_manifest_is_not_counted_as_a_seed(tmp_path):
    """A run in flight leaves a directory, and the census counted it as a finished seed.

    Fetching a five-seed run while seed 5 was still going produced four manifests and a
    fifth directory holding one file. The stage census read "0 stages failed across 5
    seeds" and "0 of 5 complete" off that, while the cascade table two sections above it
    correctly said four: the same run described over two different denominators, and the
    larger one looked like the authoritative figure. The partial seed is named instead,
    so a count is always over the seeds that produced a manifest.
    """
    import shutil
    runs = sorted(ROOT.glob("results/prod_*/run_s*/manifest.json"))
    if len(runs) < 2:
        pytest.skip("no production run directory with two finished seeds is committed")
    src = runs[0].parent.parent
    res = tmp_path / src.name
    shutil.copytree(src, res)
    for d in sorted(res.glob("run_s*")):
        if not (d / "manifest.json").exists():
            shutil.rmtree(d)
    kept = sorted(d.name for d in res.glob("run_s*"))
    (res / "run_s99").mkdir()
    (res / "run_s99" / "invocation.json").write_text("{}")
    out = tmp_path / "TABLES.md"
    r = subprocess.run([sys.executable, "-I", "experiments/report_prod.py",
                        "--library", src.name.replace("prod_", ""),
                        "--results", str(res), "--out", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = out.read_text()
    assert f"across {len(kept)} seeds" in text, text[-900:]
    assert f"complete`: **0 of {len(kept)}**" in text or \
        f"of {len(kept)}**" in text, text[-900:]
    assert "s99" in text, "the directory with no manifest is not named anywhere"


def test_a_results_directory_of_unfinished_seeds_is_refused(tmp_path):
    """Nothing to average over is a refusal, not a table of empty cells."""
    res = tmp_path / "prod_SRR23702504"
    (res / "run_s1").mkdir(parents=True)
    r = subprocess.run([sys.executable, "-I", "experiments/report_prod.py",
                        "--library", "SRR23702504", "--results", str(res),
                        "--out", str(tmp_path / "t.md")],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode != 0
    assert "not one manifest.json" in (r.stderr + r.stdout)
