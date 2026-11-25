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
    # The number is Berlin's and does not travel. What does travel is the design: a
    # shuffled control arm and the direction of a bias are properties of the procedure.
    # Printing neither is what 28 of this page's 90 rows used to do on a library nobody
    # had written about, nine of them for rows that state no number in the first place.
    assert "77.10%" in berlin["cascade"]["Bias"]
    assert "77.10%" not in hold["cascade"]["Bias"]
    assert hold["cascade"]["Bias"].endswith(make_web.NO_VALUE)
    assert "identifiable fraction exceeds the nameable" in hold["cascade"]["Bias"]
    # A procedure row is the tool's, and reads the same everywhere.
    assert berlin["cascade"]["Procedure"] == hold["cascade"]["Procedure"]
    # A sample nobody has written about still gets every entry, numbered the same.
    assert list(other) == list(berlin)
    assert other["library"]["Input"].endswith(make_web.NO_VALUE)
    bare = [(i, lab) for i, e in other.items() for lab, txt in e.items()
            if txt == make_web.NOT_MEASURED]
    assert not bare, (
        f"{len(bare)} rows say only that nothing was measured, where the row's design "
        f"holds for any water: {bare}")


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


def test_the_section_list_is_one_literal_with_no_per_sample_branch():
    """The list is in the page code once, and nothing in the render path keys on a sample.

    A second list, or a filter on the accession, would be two report formats wearing one
    generator - which is what the page did before SPEC 13.4, and what Imo asked to be
    able to check. The accession reaches the page only as data to print and as the key
    the Methods rows are looked up under, so no run identifier may appear in the code.
    """
    src = (ROOT / "experiments/make_web.py").read_text()
    assert src.count("const SUMMARY = [") == 1
    page = make_web.PAGE
    assert "SRR" not in page, "a run accession is written into the page code"
    # Comments may name a library to record a bug that printed its numbers elsewhere;
    # code may not, so the check is on the code with the comments stripped.
    code = "\n".join(ln for ln in page.split("\n")
                     if not ln.lstrip().startswith(("//", "*", "/*")))
    assert "Berlin" not in code and "berlin" not in code


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
    if "render not checked" in r.stdout:
        # The sections are assembled in the page at render time, so which ones were drawn
        # is only observable with a browser. Without one this is not asserted rather than
        # asserted weakly, and the line below still holds the generator to one list.
        pytest.skip("no playwright here, so the drawn sections cannot be read back")
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


def test_a_descriptor_whose_accession_is_not_the_runs_own_refuses(tmp_path):
    """A descriptor must not dress one run's directory in another sample's words."""
    d = _holdout()
    d["sample"]["accession"] = "SRR99999999"
    with pytest.raises(SystemExit, match="is not in the run's own input path"):
        make_web.load_page(_write(tmp_path, d), None)


def test_a_library_nobody_has_described_gets_a_page_that_says_so(tmp_path):
    """The rule changed here, on purpose, and this is the reason.

    This used to refuse: a library with no entry in `experiments/libraries.json` could
    not have a page. That is right for a registry lookup and wrong for this tool. The
    registry raising is correct, because the alternative is another library's facts, and
    a PRODUCTION sample is by definition one the registry does not know - so the old
    behaviour meant the only samples that could have a page were the ones somebody had
    already described, which is a fixture-shaped tool (SPEC 13.4).

    So an unknown library now gets a page that states it knows nothing about the water,
    and the statement is RECORDED as an override rather than made quietly, so nobody can
    read the sentence as a derived fact. Nothing is guessed in either direction.
    """
    d = _holdout()
    d["sample"] = {"accession": "SRR99999999"}
    # A run directory of its own, whose reads name the same thing, so the agreement check
    # above does not fire first and hide what this test is about.
    run = tmp_path / "run"
    run.mkdir()
    (run / "manifest.json").write_text(json.dumps(
        {"inputs": {"r1": "/x/SRR99999999_1.fastq.gz"}}))
    page = make_web.load_page(_write(tmp_path, d), run)
    assert page["unknown_library"] is True
    assert "unknown_library" in page["overrides"]
    assert "no sample facts" in page["derived"]["lede"]
    assert page["derived"]["facts"] == {}


def test_a_committed_fact_with_no_clause_refuses_rather_than_vanishing(tmp_path,
                                                                       monkeypatch):
    """Caught by the thread building the first untouched sample's page.

    `country` sat in four sample blocks with its basis and had no clause, so the sentence
    dropped it silently - and the library that the genericness test turns on is the one
    outside Europe. A fact somebody wrote down that the renderer ignores is the same
    asymmetry page_descriptor.py was written to close, so it raises now.
    """
    import page_descriptor
    assert "Germany" in page_descriptor.lede("SRR23973193")
    monkeypatch.setattr(page_descriptor, "CLAUSES",
                        tuple(c for c in page_descriptor.CLAUSES if c[0] != "country"))
    with pytest.raises(Exception, match="no clause for"):
        page_descriptor.lede("SRR23973193")


def test_a_gzipped_rows_file_is_read_and_two_spellings_refuse(tmp_path):
    """The first real production directory built a page with zero per-read rows.

    It writes reads.jsonl.gz; every earlier run here writes reads.jsonl, and the reader
    looked for the plain name only. So the residue table, the scatter and both search
    arms drew empty and the page said nothing was there, for a run over 32 M pairs. A
    directory carrying both spellings refuses, because nothing says which one the run
    wrote.
    """
    import gzip
    import shutil
    run = tmp_path / "run"
    shutil.copytree(RUN if "RUN" in globals() else
                    ROOT / "results/holdout_cascade/run_s1", run)
    plain = run / "reads.jsonl"
    with open(plain, "rb") as src, gzip.open(run / "reads.jsonl.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    page = make_web.load_page(HOLDOUT, run)
    with pytest.raises(SystemExit, match="nothing says which"):
        make_web.bundle(run, 20, 20, 1, page)
    plain.unlink()
    b = make_web.bundle(run, 20, 20, 1, page)
    assert b["scatter"]["residue"]["n_total"] > 0


def test_the_row_reader_is_imported_and_not_re_typed():
    """One copy of gzipped-or-plain, as for the shuffle and the search parser.

    Four scripts read a run directory's row dumps by name. A second copy of this logic
    is how the shuffle bug shipped twice in two days, so the reader lives in
    experiments/run_rows.py and anything reading those dumps imports it.
    """
    import run_rows
    assert "run_rows" in (ROOT / "experiments/make_web.py").read_text()
    for f in sorted((ROOT / "experiments").glob("*.py")):
        if f.name == "run_rows.py":
            continue
        src = f.read_text()
        for name in ("reads.jsonl", "layers.jsonl", "contigs.jsonl"):
            if f'{name}.gz' in src:
                assert "run_rows" in src, f"{f.name} opens {name}.gz on its own"
    assert run_rows.jsonl(ROOT / "results/holdout_cascade/run_s1", "reads.jsonl")


def test_the_chimaera_rate_is_over_the_contigs_the_flag_could_fire_on():
    """The flag needs two named ORFs, so the arm's contig count is the wrong denominator.

    Measured on SRR34217566: 6 of 15 divergent contigs and 11 of 30 near-reference
    controls carry the flag, which reads 40.0% against 36.7% over the arms and 100%
    against 84.6% over the contigs that could be flagged. The first invites "the flag
    separates the arms a little"; only the second is the quantity the Methods entry
    describes. The page derives it rather than trusting the stage, because runs produced
    before the stage carried the field would otherwise fall back to the wider
    denominator, which reads lower in the direction of a separation.
    """
    rows = ([{"arm": "divergent", "orfs": 2, "orfs_with_a_subject": 2,
              "subject_disagreement_flag": True}] * 6
            + [{"arm": "divergent", "orfs": 1, "orfs_with_a_subject": 1,
                "subject_disagreement_flag": False}] * 9
            + [{"arm": "control_near_reference", "orfs": 3, "orfs_with_a_subject": 3,
                "subject_disagreement_flag": True}] * 11
            + [{"arm": "control_near_reference", "orfs": 3, "orfs_with_a_subject": 3,
                "subject_disagreement_flag": False}] * 2
            + [{"arm": "control_near_reference", "orfs": 1, "orfs_with_a_subject": 0,
                "subject_disagreement_flag": False}] * 17)
    out = make_web._flaggable(
        rows, {"divergent": {"contigs": 15},
               "control_near_reference": {"contigs": 30}})
    assert out["divergent"]["contigs_that_could_be_flagged"] == 6
    assert out["divergent"]["fraction_flagged_of_the_contigs_that_could_be"] == 1.0
    assert out["control_near_reference"]["contigs_that_could_be_flagged"] == 13
    assert (out["control_near_reference"]
            ["fraction_flagged_of_the_contigs_that_could_be"] == 0.8462)


def test_an_arm_with_nothing_testable_is_not_reported_as_a_zero_rate():
    """No contig could be flagged and no contig was flagged render the same as 0.0."""
    out = make_web._flaggable(
        [{"arm": "divergent", "orfs": 1, "orfs_with_a_subject": 1,
          "subject_disagreement_flag": False}],
        {"divergent": {"contigs": 1}})
    assert out["divergent"]["contigs_that_could_be_flagged"] == 0
    assert out["divergent"]["fraction_flagged_of_the_contigs_that_could_be"] is None


def test_a_flag_count_above_its_own_denominator_is_refused():
    """The impossible case, asserted where the quantity is computed."""
    import pytest
    with pytest.raises(ValueError, match="two or more named ORFs"):
        make_web._flaggable(
            [{"arm": "divergent", "orfs": 1, "orfs_with_a_subject": 1,
              "subject_disagreement_flag": True}],
            {"divergent": {"contigs": 1}})


def _page_function(name: str) -> str:
    """One function lifted out of the page code, so node can be asked what it returns."""
    src = make_web.PAGE
    start = src.index(f"function {name}(")
    depth, i = 0, src.index("{", start)
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start:j + 1]
    raise AssertionError(f"{name} does not close")


def test_a_dated_kraken2_directory_still_counts_as_a_tier_that_names():
    """The tier name carries the database release, and a literal map stopped matching.

    `hapax fetch-db` names each Kraken2 directory after the release it downloaded, so the
    cascade's tiers arrive as `k2_viral_20261022.d` rather than `k2_viral`. Two panels
    held a literal map of the undated names: every Kraken2 read moved out of the named
    band into recognised-but-unnamed, the named band read 2.24% against the cascade's
    4.07%, and the reconciliation panel asserted that its two named rows agreed across a
    fifteen-point gap. A tier is classified by what it is, and an unrecognised one is
    reported as unknown so it cannot land in a band by default.
    """
    if subprocess.run(["node", "-e", "0"], capture_output=True).returncode:
        pytest.skip("no node here")
    probe = _page_function("tierKind") + """
const cases = {
  'k2_viral_20261022.d': 'names', 'k2_pluspf_16_GB_20261022.d': 'names',
  'k2_viral': 'names', 'prot_viral': 'names', 'prot_sprot': 'names',
  'rrna': 'recognises', 'artefact': 'recognises',
  'selfasm_sample_contigs_rrna': 'recognises',
  'selfasm_sample_contigs_identified_norrna': 'recognises',
  'observe_sample_contigs_unidentified': 'observes',
  'something_new': 'unknown', '': 'unknown'};
const bad = Object.entries(cases).filter(([k, v]) => tierKind(k) !== v)
  .map(([k, v]) => k + ' read ' + tierKind(k) + ', wanted ' + v);
if (bad.length) { console.log(bad.join('; ')); process.exit(1); }
console.log('ok');
"""
    r = subprocess.run(["node", "-e", probe], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_the_page_never_pairs_the_two_independent_maxima():
    """The producer ships a field saying not to, and the page did it anyway.

    `host_confidence_summary` carries `best_identity_percent` and `best_aligned_aa` as
    maxima over the rows taken separately, usually of different subjects, and carries
    `best_fields_are_independent` saying in words that they must not be printed as one
    match. The warning panel honoured that and took the closest row's own length; the
    headline card above it printed the two maxima together and read "96.3%, and over only
    60 amino acids" where that 96.3% row is a 31 aa alignment. A pairing of two maxima
    describes a subject that does not exist, so the page code may not read the maximum
    length at all: `closest_row_aligned_aa`, or the chosen row's own field, is the only
    length that may sit beside an identity.
    """
    page = make_web.PAGE
    code = "\n".join(ln for ln in page.split("\n")
                      if not ln.lstrip().startswith(("//", "*", "/*")))
    assert "best_aligned_aa" not in code, (
        "the page reads the maximum aligned length, which belongs to a different "
        "subject than the maximum identity")
    assert "longest_alignment_aa" not in code
