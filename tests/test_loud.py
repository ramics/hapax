"""Tests for the fail-loudly pass: every check here is one this project was bitten by.

Run with `python3 -m pytest tests/test_loud.py`. Needs mappy (`pip install -r requirements.txt`). External
tools are replaced by stubs written into a temp directory, so nothing here touches a real
index; what is under test is whether a stage refuses, and what it records, not the numbers.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT))

import loud  # noqa: E402
import hapax  # noqa: E402

PY = sys.executable


# ------------------------------------------------------------------------------ fixtures

def _fastq_pair(d: Path, n: int = 50) -> tuple[Path, Path]:
    import random
    rng = random.Random(1)
    r1, r2 = d / "r_1.fastq.gz", d / "r_2.fastq.gz"
    with gzip.open(r1, "wt") as a, gzip.open(r2, "wt") as b:
        for i in range(n):
            s1 = "".join(rng.choice("ACGT") for _ in range(107))
            s2 = "".join(rng.choice("ACGT") for _ in range(107))
            a.write(f"@p{i}/1\n{s1}\n+\n{'I' * 107}\n")
            b.write(f"@p{i}/2\n{s2}\n+\n{'I' * 107}\n")
    return r1, r2


def _exe(path: Path, body: str) -> Path:
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def _stub_kraken(d: Path) -> Path:
    """Classifies nothing, one line per pair, like Kraken2's per-read output."""
    return _exe(d / "classify", """
import sys
a = sys.argv
out, fq1 = a[a.index("-O") + 1], a[-2]
with open(fq1) as f, open(out, "w") as o:
    for i, line in enumerate(f):
        if i % 4 == 0:
            o.write(f"U\\t{line[1:].split()[0]}\\t0\\t107|107\\t0:73\\n")
""")


def _stub_diamond(d: Path) -> Path:
    """Searches and finds nothing: writes an empty table, which is a real answer."""
    return _exe(d / "diamond", """
import sys
a = sys.argv
if len(a) > 1 and a[1] == "version":
    print("diamond stub"); sys.exit(0)
open(a[a.index("--out") + 1], "w").close()
""")


def _k2db(d: Path, name: str) -> Path:
    db = d / name
    db.mkdir()
    for f in ("hash.k2d", "taxo.k2d", "opts.k2d"):
        (db / f).write_bytes(b"x")
    return db


def _fasta(path: Path, n: int = 3) -> Path:
    import random
    rng = random.Random(7)
    with open(path, "w") as h:
        for i in range(n):
            h.write(f">c{i}\n{''.join(rng.choice('ACGT') for _ in range(800))}\n")
    return path


@pytest.fixture()
def kit(tmp_path):
    r1, r2 = _fastq_pair(tmp_path)
    return Namespace(
        d=tmp_path, r1=r1, r2=r2, classify=_stub_kraken(tmp_path),
        diamond=_stub_diamond(tmp_path), k2=_k2db(tmp_path, "k2_tiny"),
        silva=_fasta(tmp_path / "silva.fa"), asm=_fasta(tmp_path / "asm.fa"),
        dmnd=tmp_path / "tiny.dmnd")


def _needs_mappy():
    """Fail, by name, when the cascade cannot get past its own preflight in this container.

    Not a skip: mappy is pip-installable and declared in requirements.txt, so a container
    without it is misconfigured, not a place where the test cannot run. Without this the
    failure read `assert 2 == 0`, with the cascade's correct refusal buried in its stderr.
    """
    import importlib.util
    if importlib.util.find_spec("mappy") is None:
        pytest.fail("mappy is not installed in this interpreter, so subtract_known.py "
                    "refuses to start (as it should); pip install -r requirements.txt",
                    pytrace=False)


def _sk(kit, *extra, work=None, out=None):
    kit.dmnd.write_bytes(b"x")
    return [PY, "-I", str(ROOT / "experiments/subtract_known.py"),
            "--seed", "1", "--pairs", "50", "--r1", str(kit.r1), "--r2", str(kit.r2),
            "--k2-dbs", str(kit.k2), "--classify-bin", str(kit.classify),
            "--silva", str(kit.silva), "--selfasm", str(kit.asm),
            "--dmnd", str(kit.dmnd), "--diamond-bin", str(kit.diamond),
            "--threads", "1", "--rows-out",
            "--work", str(work or kit.d / "work"), "--out", str(out or kit.d / "out"),
            *extra]


# ------------------------------------------------------------------------------ loud.py

def test_preflight_lists_every_problem_at_once(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    notexe = tmp_path / "tool"
    notexe.write_text("")
    pf = loud.Preflight("t")
    pf.file(tmp_path / "nope.fa", "a")
    pf.file(empty, "b")
    pf.file(tmp_path, "c")
    pf.exe(notexe, "d")
    pf.kraken2_db(tmp_path / "nodb")
    assert len(pf.problems) == 5
    with pytest.raises(SystemExit) as e:
        pf.done()
    assert e.value.code == loud.EXIT_PREFLIGHT


def test_cache_goes_stale_when_its_source_changes(tmp_path):
    src = tmp_path / "ref.fa"
    src.write_text(">a\nACGT\n")
    cache = tmp_path / "ref.mmi"
    cache.write_text("idx")
    assert not loud.cache_is_current(cache, [src]), "an unmarked cache is stale"
    loud.mark_cache(cache, [src])
    assert loud.cache_is_current(cache, [src])
    src.write_text(">a\nACGTACGT\n")
    assert not loud.cache_is_current(cache, [src])


def test_consulted_refuses_impossible_counts(tmp_path):
    c = loud.Consulted()
    with pytest.raises(SystemExit):
        c.searched("x", tmp_path, queries=3, hits=4)


def test_run_tool_shows_stderr(capsys, tmp_path):
    bad = _exe(tmp_path / "bad", "import sys; sys.stderr.write('index is corrupt\\n'); "
                                 "sys.exit(3)")
    with pytest.raises(SystemExit):
        loud.run_tool([bad], "bad tool")
    assert "index is corrupt" in capsys.readouterr().err


# ------------------------------------------------------------------------ subtract_known

def test_cascade_refuses_an_absent_index_before_reading_anything(kit):
    cmd = _sk(kit)
    cmd[cmd.index("--k2-dbs") + 1] = str(kit.d / "k2_missing")
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 2
    assert "k2_missing" in p.stderr
    assert not (kit.d / "out").exists(), "refused runs must not write results"


def test_cascade_refuses_operating_confidence_outside_sweep(kit):
    p = subprocess.run(_sk(kit, "--k2-confidence", "0.0", "0.2"),
                       capture_output=True, text=True)
    assert p.returncode == 2 and "operating-confidence" in p.stderr


def test_cascade_records_what_it_consulted(kit):
    _needs_mappy()
    p = subprocess.run(_sk(kit), capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    m = json.loads((kit.d / "out" / "manifest_s1.json").read_text())
    assert m["script_sha256"] and len(m["script_sha256"]) == 16
    got = {e["path"]: e for e in m["consulted"]}
    for path in (kit.k2, kit.silva, kit.asm, kit.dmnd):
        assert str(path) in got, f"{path} requested but not in consulted"
        assert got[str(path)]["state"] == "searched"
    # A protein search that returned nothing is a count of zero beside its query count.
    assert got[str(kit.dmnd)]["hits"] == 0 and got[str(kit.dmnd)]["queries"] > 0


def test_stale_index_cache_is_rebuilt_not_served(kit):
    """Same work dir, same file stem, different reference: the old cache must not be used."""
    _needs_mappy()
    work = kit.d / "work"
    assert subprocess.run(_sk(kit, work=work), capture_output=True).returncode == 0
    side = work / "asm.mmi.source.json"
    before = json.loads(side.read_text())
    _fasta(kit.asm, n=5)
    assert subprocess.run(_sk(kit, work=work), capture_output=True).returncode == 0
    assert json.loads(side.read_text()) != before


# ------------------------------------------------------------------------------- hapax.py

def _run(tmp_path, **kw) -> hapax.Run:
    args = Namespace(out=tmp_path / "run", work=None, seed=1, stages=hapax.STAGE_ORDER,
                     dmnd=[tmp_path / "v.dmnd"], r1=tmp_path / "r1", r2=tmp_path / "r2")
    for k, v in kw.items():
        setattr(args, k, v)
    return hapax.Run(args)


def test_stage_fails_when_a_requested_index_is_not_in_its_argv(tmp_path):
    run = _run(tmp_path)
    idx = tmp_path / "k2_viral"
    argv = ["--seed", "1", "--r1", "a", "--r2", "b", "--trim-r2", "0", "--selfasm"]
    e = run.call("cascade", "subtract_known.py", argv, {}, None, requested=[idx])
    assert e["status"] == "failed" and str(idx) in e["reason"]


def test_stage_fails_when_its_manifest_does_not_name_a_requested_index(tmp_path):
    run = _run(tmp_path)
    idx = tmp_path / "v.dmnd"
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "out.json").write_text("{}")
    (stage / "manifest.json").write_text(json.dumps({
        "script_sha256": "x", "consulted": [{"path": "/somewhere/else.dmnd",
                                             "state": "searched"}]}))
    # A script that exits 0 and does nothing: what was run is not under test, only what
    # its manifest claims.
    script = ROOT / "experiments" / "_test_noop.py"
    script.write_text("")
    try:
        e = run.call("reads", script.name, ["--dmnd", str(idx)],
                     {"x.json": stage / "out.json"}, stage / "manifest.json",
                     requested=[idx])
    finally:
        script.unlink()
    assert e["status"] == "failed" and "does not include" in e["reason"]


def test_stage_fails_without_a_script_hash(tmp_path):
    run = _run(tmp_path)
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "manifest.json").write_text(json.dumps({"git_sha": "abc"}))
    script = ROOT / "experiments" / "_test_noop.py"
    script.write_text("")
    try:
        e = run.call("linkage", script.name, [], {}, stage / "manifest.json")
    finally:
        script.unlink()
    assert e["status"] == "failed" and "script_sha256" in e["reason"]


def _reads(run, rows):
    with open(run.dir / "reads.jsonl", "w") as h:
        for r in rows:
            h.write(json.dumps(r) + "\n")


def test_strains_never_writes_an_empty_negative_when_nothing_was_searched(tmp_path):
    run = _run(tmp_path)
    _reads(run, [{"name": "a", "set": "residue"}, {"name": "b", "set": "known"}])
    (run.dir / "fields.json").write_text(json.dumps({"consulted": []}))
    e = run.strains()
    assert e["status"] == "failed"
    s = json.loads((run.dir / "strains.json").read_text())
    assert s["state"] == "not searched" and s["subjects"] is None


def test_strains_counts_a_searched_control_with_no_hits_as_searched(tmp_path):
    run = _run(tmp_path)
    hit = {"subject": "S1", "pident": 70.0, "aln_len": 35, "evalue": 1e-9,
           "bitscore": 60.0}
    _reads(run, [{"name": "a", "set": "residue", "closest_v": hit},
                 {"name": "b", "set": "known"}])
    (run.dir / "fields.json").write_text(json.dumps({"consulted": [
        {"role": "closest_v residue", "path": "v.dmnd", "state": "searched",
         "queries": 1, "hits": 1},
        {"role": "closest_v known", "path": "v.dmnd", "state": "searched",
         "queries": 1, "hits": 0}]}))
    assert run.strains()["status"] == "ok"
    s = json.loads((run.dir / "strains.json").read_text())
    assert s["known_pile_searched_reads"] == 1
    assert s["subjects"][0]["reading"].startswith("only the residue points here")


def test_finish_counts_known_rows_by_field_not_by_substring(tmp_path):
    run = _run(tmp_path, trim_r2=14)
    # A known row whose subject name contains the word "residue" used to be counted as
    # residue; a residue row with no "set" key used to be counted as known.
    _reads(run, [{"name": "a", "set": "known", "subject": "residue_like_protein"},
                 {"name": "b", "set": "residue"}, {"name": "c", "set": "known"}])
    run.finish()
    m = json.loads((run.dir / "manifest.json").read_text())
    assert m["known_positive_rows"] == 2


def test_hapax_refuses_to_start_on_a_missing_input(tmp_path):
    p = subprocess.run([PY, str(ROOT / "hapax.py"), "run", str(tmp_path / "no_1.fq.gz"),
                        str(tmp_path / "no_2.fq.gz"), "--seed", "1",
                        "--out", str(tmp_path / "run"), "--stages", "linkage",
                        # --library is required for the same reason --trim-r2 is: the
                        # panels a run is checked against are properties of the library.
                        "--library", "SRR23973193",
                        "--trim-r2", "0"],
                       capture_output=True, text=True)
    assert p.returncode == 2, p.stderr
    m = json.loads((tmp_path / "run" / "manifest.json").read_text())
    assert m["status"] == "refused to start" and m["problems"]
    assert not (tmp_path / "run" / "stages" / "linkage.log").exists()


# ------------------------------------------------------------- sample-specific defaults

def _hapax(tmp_path, kit, *extra):
    return subprocess.run(
        [PY, str(ROOT / "hapax.py"), "run", str(kit.r1), str(kit.r2), "--seed", "1",
         "--out", str(tmp_path / "run"), "--trim-r2", "0", "--pairs", "50",
         "--library", "SRR23973193",
         "--k2-dbs", str(kit.k2), "--silva", str(kit.silva), "--dmnd", str(kit.dmnd),
         "--classify-bin", str(kit.classify), "--diamond-bin", str(kit.diamond), *extra],
        capture_output=True, text=True)


def test_cascade_refuses_without_this_samples_contigs(tmp_path, kit):
    kit.dmnd.write_bytes(b"x")
    p = _hapax(tmp_path, kit, "--stages", "cascade")
    assert p.returncode == 2 and "--contigs not given" in p.stderr


def test_report_refuses_another_samples_assembly_by_default(tmp_path, kit):
    """The chimaera test used to map against Berlin's assembly for any library."""
    kit.dmnd.write_bytes(b"x")
    p = _hapax(tmp_path, kit, "--contigs", str(kit.asm), "--stages", "cascade", "reads",
               "--report-without", "rfam", "pfam", "crispr", "cross_sample")
    assert p.returncode == 2 and "--assembly" in p.stderr


def test_census_refuses_without_this_samples_label_contigs(tmp_path, kit):
    kit.dmnd.write_bytes(b"x")
    p = _hapax(tmp_path, kit, "--contigs", str(kit.asm), "--stages", "cascade", "census")
    assert p.returncode == 2 and "--rrna-contigs" in p.stderr


def test_call_refuses_a_stage_missing_a_sample_specific_flag(tmp_path):
    run = _run(tmp_path, report_without=[])
    e = run.call("census", "linkage_census.py", ["--r1", "x", "--r2", "y"], {}, None)
    assert e["status"] == "failed" and "--trim-r2" in e["reason"]


def test_trim_is_required(tmp_path):
    p = subprocess.run([PY, str(ROOT / "hapax.py"), "run", "a", "b", "--seed", "1",
                        "--out", str(tmp_path / "r")], capture_output=True, text=True)
    assert p.returncode == 2 and "--trim-r2" in p.stderr
