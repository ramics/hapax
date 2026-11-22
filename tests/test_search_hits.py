"""search_hits.py: the parser every new script reads BLAST and DIAMOND output through.

Run with `python3 -m pytest tests/test_search_hits.py`. Pure Python, no search tools: the
tabular files are synthesised with the awkward cases the copies disagree on: tied
bitscores, several HSPs per pair, gapped alignments, both mates of a pair, and HSPs either
side of the significance floor.

What is checked, and why:
- Best-hit SELECTION is identical to the grandfathered copies, so switching a script to
  this module picks the same rows.
- Union coverage is identical to the copies on the same HSPs, and any-subject coverage is
  never below per-subject coverage when both see the same floor. A negative gap is
  impossible under one definition, and PR #33's table carries hundreds of them because its
  two unions used different floors.
- Coverage is in the query's own units, so a translated search cannot cap at 1/3.
- Refusals: a missing file, a row with the wrong field count, an unknown program.
"""

from __future__ import annotations

import ast
import collections
import random
import sys
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments"
sys.path.insert(0, str(EXP))

import search_hits as sh  # noqa: E402

FLOOR = 1e-5


def _extract(fname: str, func: str, extra: dict | None = None):
    path = EXP / fname
    if not path.exists():
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == func:
            ns = {"Path": Path, "collections": collections, "SIGNIFICANCE_EVALUE": FLOOR}
            ns.update(extra or {})
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), ns)
            return ns[func]
    return None


def _synth(seed: int = 1, n_queries: int = 60, mates: bool = False):
    """Rows in DEFAULT_COLUMNS order for a DIAMOND blastx-like search of nt queries."""
    r = random.Random(seed)
    rows = []
    for qi in range(n_queries):
        qlen = r.choice([93, 107, 300, 1500, 3990])
        names = [f"q{qi}/1", f"q{qi}/2"] if mates else [f"q{qi}"]
        for q in names:
            for _ in range(r.randint(1, 6)):
                s = f"s{r.randint(0, 5)}"
                a = r.randint(1, max(1, qlen - 30))
                b = min(qlen, a + r.randint(30, qlen))
                if r.random() < 0.3:
                    a, b = b, a                                   # reverse-frame HSP
                aa = (abs(b - a) + 1) // 3 + r.randint(0, 4)        # residues, gaps add
                bits = r.choice([30.0, 45.5, 50.0, 50.0, 88.2, 120.0])   # ties on purpose
                ev = r.choice([1e-30, 1e-6, 1e-5, 2e-5, 0.1, 5.0])      # either side of floor
                slen = r.randint(aa, aa + 400)
                ss = r.randint(1, slen - aa + 1)
                rows.append([q, s, f"{r.uniform(25, 100):.1f}", aa, r.randint(0, aa // 2),
                             r.randint(0, 3), a, b, ss, ss + aa - 1, f"{ev:g}", f"{bits}",
                             qlen, slen])
    return rows


def _write(tmp_path: Path, rows, name="hits.tsv") -> Path:
    p = tmp_path / name
    p.write_text("".join("\t".join(map(str, x)) + "\n" for x in rows))
    return p


# ------------------------------------------------------------------------- selection
def test_best_selection_matches_holdout_cascade(tmp_path):
    old = _extract("holdout_cascade.py", "parse_best")
    if old is None:
        pytest.skip("holdout_cascade.parse_best is gone")
    p = _write(tmp_path, _synth())
    new = sh.best_per_query(sh.read_hits(p, program="diamond_blastx"))
    ref = old(p, aa=True) if "aa" in old.__code__.co_varnames else old(p)
    assert set(new) == set(ref)
    for q in new:
        assert (new[q]["sseqid"], new[q]["bitscore"], new[q]["evalue"], new[q]["qstart"]) == \
               (ref[q]["subject"], ref[q]["bitscore"], ref[q]["evalue"], ref[q]["qstart"])


def test_best_selection_matches_profile_reach(tmp_path):
    old = _extract("profile_reach.py", "parse_diamond")
    if old is None:
        pytest.skip("profile_reach.parse_diamond is gone")
    six = [[x[0], x[1], x[2], x[3], x[10], x[11]] for x in _synth(seed=2)]
    p = _write(tmp_path, six)
    cols = ("qseqid", "sseqid", "pident", "length", "evalue", "bitscore")
    new = sh.best_per_query(sh.read_hits(p, cols, program="diamond_blastx"))
    ref = old(p)
    assert {q: (v["sseqid"], v["length"]) for q, v in new.items()} == \
           {q: (v["subject"], v["length"]) for q, v in ref.items()}


def test_pair_key_takes_the_better_mate_first_seen_on_tie(tmp_path):
    """What subtract_known, residue_report and full_residue_search do inline:
    `name = q.rsplit("/", 1)[0]`, strict `>` on bitscore."""
    rows = _synth(seed=3, mates=True)
    p = _write(tmp_path, rows)
    new = sh.best_per_query(sh.read_hits(p, program="diamond_blastx"), key=sh.pair_key)
    ref: dict[str, tuple] = {}
    for x in rows:                                  # the inline copy, verbatim in shape
        name, bits = x[0].rsplit("/", 1)[0], float(x[11])
        if bits > ref.get(name, (-1.0,))[0]:
            ref[name] = (bits, x[1], x[0])
    assert {k: (v["bitscore"], v["sseqid"], v["qseqid"]) for k, v in new.items()} == ref


# -------------------------------------------------------------------------- coverage
def test_union_any_matches_holdout_cascade(tmp_path):
    old = _extract("holdout_cascade.py", "union_coverage_any")
    if old is None:
        pytest.skip("holdout_cascade.union_coverage_any is gone")
    p = _write(tmp_path, _synth(seed=4))
    new = sh.union_query_coverage(sh.read_hits(p, program="diamond_blastx"),
                                  per_subject=False, evalue_below=FLOOR)
    ref = old(p, FLOOR)
    assert set(new) == set(ref)
    for q in new:
        assert round(new[q]["covered_fraction"], 4) == ref[q]["covered_fraction"]
        assert new[q]["significant_subjects"] == ref[q]["significant_subjects"]


def test_union_per_subject_matches_holdout_cascade(tmp_path):
    """Since PR #33's fix, holdout_cascade.union_coverage takes the same floor as its
    any-subject union. Matches at the floor and with no floor at all."""
    old = _extract("holdout_cascade.py", "union_coverage")
    if old is None:
        pytest.skip("holdout_cascade.union_coverage is gone")
    p = _write(tmp_path, _synth(seed=5))
    rows = sh.read_hits(p, program="diamond_blastx")
    for floor in (FLOOR, float("inf")):
        new = sh.union_query_coverage(rows, per_subject=True,
                                      evalue_below=None if floor == float("inf") else floor)
        ref = old(p, floor)
        assert {q: round(v["covered_fraction"], 4) for q, v in new.items()} == \
               {q: v["covered_fraction"] for q, v in ref.items()}


def test_gap_matches_coverage_gap_null_with_matched_floors(tmp_path):
    """coverage_gap_null.py floors both unions; with the same floor here, same numbers."""
    old = _extract("coverage_gap_null.py", "gap_stats") or _extract(
        "coverage_gap_null.py", "gaps")
    if old is None:
        pytest.skip("coverage_gap_null's union function is not extractable by name")
    ten = [[x[0], x[1], x[2], x[3], x[6], x[7], x[10], x[11], x[12], x[13]]
           for x in _synth(seed=6)]
    p = _write(tmp_path, ten)
    cols = ("qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore",
            "qlen", "slen")
    rows = sh.read_hits(p, cols, program="diamond_blastx")
    per = sh.union_query_coverage(rows, per_subject=True, evalue_below=FLOOR)
    anyu = sh.union_query_coverage(rows, per_subject=False, evalue_below=FLOOR)
    ref = old(p, 0) if old.__code__.co_argcount >= 2 else old(p)
    for q, v in ref.items():
        assert round(per[q]["covered_fraction"], 4) == v["per_subject"]
        assert round(anyu[q]["covered_fraction"], 4) == v["any_subject"]


@pytest.mark.parametrize("floor", [None, FLOOR, 0.1])
def test_gap_is_never_negative_under_one_floor(tmp_path, floor):
    p = _write(tmp_path, _synth(seed=7, n_queries=200))
    rows = sh.read_hits(p, program="diamond_blastx")
    per = sh.union_query_coverage(rows, per_subject=True, evalue_below=floor)
    anyu = sh.union_query_coverage(rows, per_subject=False, evalue_below=floor)
    assert set(per) == set(anyu)
    for q in per:
        assert anyu[q]["covered_fraction"] >= per[q]["covered_fraction"]


def test_mismatched_floors_can_go_negative(tmp_path):
    """The shape of PR #33's table, reproduced: per-subject unfloored, any-subject floored."""
    p = _write(tmp_path, _synth(seed=7, n_queries=200))
    rows = sh.read_hits(p, program="diamond_blastx")
    per = sh.union_query_coverage(rows, per_subject=True)
    anyu = sh.union_query_coverage(rows, per_subject=False, evalue_below=FLOOR)
    assert any(anyu[q]["covered_fraction"] < per[q]["covered_fraction"] for q in anyu)


def test_translated_coverage_is_in_query_units(tmp_path):
    """A full-length blastx HSP covers the whole query: 1.0, not 0.333."""
    row = ["c1", "p1", "66.5", "1330", "0", "0", "1", "3990", "1", "1330", "1e-200",
           "840", "3990", "1330"]
    p = _write(tmp_path, [row])
    r = sh.read_hits(p, program="diamond_blastx")[0]
    assert sh.query_coverage(r) == 1.0
    assert sh.subject_coverage(r) == 1.0
    assert sh.UNITS["diamond_blastx"] == ("nt", "aa", "aa")


def test_query_coverage_never_exceeds_one(tmp_path):
    p = _write(tmp_path, _synth(seed=8, n_queries=200))
    for r in sh.read_hits(p, program="diamond_blastx"):
        assert 0 < sh.query_coverage(r) <= 1.0


# -------------------------------------------------------------------------- refusals
def test_missing_file_is_not_zero_hits(tmp_path):
    with pytest.raises(sh.SearchOutputError):
        sh.read_hits(tmp_path / "never_written.tsv", program="blastn")


def test_empty_file_is_zero_hits(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text("")
    assert sh.read_hits(p, program="blastn") == []


def test_wrong_field_count_raises(tmp_path):
    p = _write(tmp_path, [["q", "s", "99.0", "100", "1e-50", "200"]])
    with pytest.raises(sh.SearchOutputError):
        sh.read_hits(p, program="blastn")              # DEFAULT_COLUMNS has 14


def test_unknown_program_raises(tmp_path):
    p = _write(tmp_path, _synth()[:1])
    with pytest.raises(sh.SearchOutputError):
        sh.read_hits(p, program="blast")


def test_columns_must_include_evalue(tmp_path):
    p = _write(tmp_path, [["q", "s", "200"]])
    with pytest.raises(sh.SearchOutputError):
        sh.read_hits(p, ("qseqid", "sseqid", "bitscore"), program="blastn")


def test_outfmt_args_round_trip():
    assert sh.outfmt_args(diamond=True)[:2] == ["--outfmt", "6"]
    assert sh.outfmt_args(diamond=False) == ["-outfmt", "6 " + " ".join(sh.DEFAULT_COLUMNS)]


# ------------------------------------------------------------------ the gap, asserted
def test_coverage_gap_matches_the_two_unions(tmp_path):
    p = _write(tmp_path, _synth(seed=9, n_queries=200))
    rows = sh.read_hits(p, program="diamond_blastx")
    g = sh.coverage_gap(rows, evalue_below=FLOOR)
    per = sh.union_query_coverage(rows, per_subject=True, evalue_below=FLOOR)
    anyu = sh.union_query_coverage(rows, per_subject=False, evalue_below=FLOOR)
    assert set(g) == set(per)
    for q, v in g.items():
        assert v["gap"] >= 0
        assert v["gap"] == anyu[q]["covered_fraction"] - per[q]["covered_fraction"]


def test_floor_is_strict_at_exactly_the_floor():
    """The > versus >= bug: a row AT the floor is not significant, on every path."""
    at = {"evalue": FLOOR}
    assert not sh.passes_floor(at, FLOOR)
    assert sh.passes_floor({"evalue": FLOOR * 0.999}, FLOOR)
    assert sh.passes_floor(at, None)


def test_coverage_gap_stops_the_run_on_a_violation(monkeypatch, tmp_path):
    """Simulate the drift the function exists to prevent: a per-subject union that saw a
    row the any-subject union did not. It must raise, not write a negative gap."""
    p = _write(tmp_path, _synth(seed=10, n_queries=50))
    rows = sh.read_hits(p, program="diamond_blastx")
    real = sh.union_query_coverage

    def drifted(rs, *, per_subject, evalue_below=None):
        if not per_subject:                         # any-subject side loses one HSP
            rs = [r for r in rs if r is not rows[0]]
        return real(rs, per_subject=per_subject, evalue_below=evalue_below)

    monkeypatch.setattr(sh, "union_query_coverage", drifted)
    # Make query 0's widest HSP visible to the per-subject side only.
    q0 = [r for r in rows if r["qseqid"] == rows[0]["qseqid"]]
    widest = max(q0, key=lambda r: abs(r["qend"] - r["qstart"]))
    rows.remove(widest)
    rows.insert(0, widest)
    with pytest.raises(sh.InvariantViolation):
        sh.coverage_gap(rows, evalue_below=None)
