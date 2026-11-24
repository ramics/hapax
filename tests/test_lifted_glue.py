"""The package's lifted glue gives exactly what the committed scripts give.

`hapax.assembly_input` and `hapax.measure` are copies of `experiments/trim_mate2.py`,
`experiments/megahit_rename.py` and `experiments/measure_library.py`, kept because those
scripts are the committed producers of results whose `script_sha256` must keep matching.
A copy is safe only while it is the same instrument, and CLAUDE.md records what happens
otherwise: two shuffles that each looked correct and differed. So each test runs the
script and the function on one input and compares the outputs byte for byte, the refusals
included. When one of these fails, change both copies or neither.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from hapax.assembly_input import GlueError, rename_megahit, trim_mate2
from hapax.measure import measure_to

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments"
READS = ROOT / "tests/fixtures/reads"
R1, R2 = READS / "fixture_R1.fastq.gz", READS / "fixture_R2.fastq.gz"


def script(name, *args):
    return subprocess.run([sys.executable, "-I", str(EXP / name), *map(str, args)],
                          capture_output=True, text=True)


def fastq(path, seqs):
    path.write_text("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs)))


@pytest.mark.parametrize("trim", [0, 14, 200])
def test_trim_is_byte_identical_to_the_script(tmp_path, trim):
    r = script("trim_mate2.py", "--r1", R1, "--r2", R2, "--trim-r2", trim,
               "--out-dir", tmp_path / "s", "--prefix", "x")
    assert r.returncode == 0, r.stderr
    man = trim_mate2(R1, R2, trim, tmp_path / "p", "x")
    for mate in ("1", "2"):
        assert (tmp_path / f"s/x_{mate}.fq").read_bytes() == \
            (tmp_path / f"p/x_{mate}.fq").read_bytes()
    smap = json.loads((tmp_path / "s/x_trim_manifest.json").read_text())
    for key in ("pairs", "trim_r2", "mate2_shorter_than_the_prefix"):
        assert smap[key] == man[key], key
    assert man["pairs"] > 0


def test_trim_refuses_unequal_mates_like_the_script(tmp_path):
    fastq(tmp_path / "a.fq", ["A" * 20, "A" * 20])
    fastq(tmp_path / "b.fq", ["G" * 20])
    r = script("trim_mate2.py", "--r1", tmp_path / "a.fq", "--r2", tmp_path / "b.fq",
               "--trim-r2", 3, "--out-dir", tmp_path / "s", "--prefix", "x")
    assert r.returncode != 0 and "more records than" in r.stderr
    with pytest.raises(GlueError, match="more records than"):
        trim_mate2(tmp_path / "a.fq", tmp_path / "b.fq", 3, tmp_path / "p", "x")


def _megahit(path, headers):
    path.write_text("".join(f">{h}\nACGTACGT\nACG\n" for h in headers))


def test_rename_is_byte_identical_to_the_script(tmp_path):
    src = tmp_path / "final.contigs.fa"
    _megahit(src, ["k141_12 flag=1 multi=3.7261 len=486", "k141_3 flag=0 multi=12.0 len=11"])
    r = script("megahit_rename.py", "--contigs", src, "--out", tmp_path / "s.fa",
               "--tag", "LIB")
    assert r.returncode == 0, r.stderr
    man = rename_megahit(src, tmp_path / "p.fa", "LIB")
    assert (tmp_path / "s.fa").read_bytes() == (tmp_path / "p.fa").read_bytes()
    assert man["contigs"] == 2


@pytest.mark.parametrize("headers", [["k141_1 flag=1 len=486"], []])
def test_rename_refuses_what_the_script_refuses(tmp_path, headers):
    src = tmp_path / "final.contigs.fa"
    _megahit(src, headers)
    r = script("megahit_rename.py", "--contigs", src, "--out", tmp_path / "s.fa",
               "--tag", "LIB")
    assert r.returncode != 0
    with pytest.raises(GlueError):
        rename_megahit(src, tmp_path / "p.fa", "LIB")


VOLATILE = {"provenance", "wall_seconds"}


@pytest.mark.parametrize("seeds,per_seed", [(1, 2000), (2, 1000)])
def test_measure_matches_the_script(tmp_path, seeds, per_seed):
    r = script("measure_library.py", "--r1", R1, "--r2", R2, "--seeds", seeds,
               "--pairs-per-seed", per_seed, "--accession", "fixture",
               "--out", tmp_path / "s.json")
    assert r.returncode == 0, r.stderr
    got = measure_to(R1, R2, per_seed, seeds, tmp_path / "p.json", accession="fixture")
    want = json.loads((tmp_path / "s.json").read_text())
    got = json.loads(json.dumps(got))
    assert {k: v for k, v in got.items() if k not in VOLATILE} == \
        {k: v for k, v in want.items() if k not in VOLATILE}
    assert got["pairs"] > 0
    assert got["provenance"]["producer"] == "hapax.measure.measure_to"


def test_measure_refuses_overlapping_seeds_like_the_script(tmp_path):
    r = script("measure_library.py", "--r1", R1, "--r2", R2, "--seeds", 5,
               "--pairs-per-seed", 10_000_000)
    assert r.returncode != 0 and "DISJOINT" in r.stderr
    with pytest.raises(ValueError, match="DISJOINT"):
        measure_to(R1, R2, 10_000_000, 5, tmp_path / "p.json")


def test_provenance_reads_the_checkout_not_the_install_location(tmp_path):
    """The end-to-end job installs the package into site-packages, which is no checkout.
    The first version of `hapax.provenance` asked git about its own file's directory and
    refused there; the commit has to come from the checkout `repo_root` finds."""
    import shutil

    shutil.copytree(ROOT / "src/hapax", tmp_path / "site/hapax")
    env = {"HAPAX_REPO": str(ROOT), "PATH": __import__("os").environ["PATH"],
           "PYTHONPATH": str(tmp_path / "site")}
    r = subprocess.run([sys.executable, "-c",
                        "import json, hapax.provenance as p; "
                        "print(json.dumps(p.record(p.__file__, 'x')))"],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()
    assert json.loads(r.stdout)["git_sha"].startswith(head)
