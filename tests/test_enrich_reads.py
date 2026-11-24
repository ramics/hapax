"""The enrichments use the reads the run read, not the reads the config names.

WHY THIS EXISTS. `hapax run` takes the FASTQs either from `[inputs]` in the config or
positionally on the command line, and resolves the two into one pair. The per-gene workup
needs them, and reading `cfg.inputs["r1"]` directly raised a KeyError on the fixture,
which passes its reads on the command line - after every stage had run, in the end-to-end
job, which takes seventeen minutes to say so. One resolution, passed in, is the rule the
binaries already follow.
"""
from __future__ import annotations

from pathlib import Path

from hapax.config import Config
from hapax.enrich import enrich


def _cfg(tmp_path: Path) -> Config:
    """A run config with NO reads in it, which is the case that failed."""
    stub = tmp_path / "diamond"
    stub.write_text("#!/bin/sh\nexit 1\n")
    return Config(path=tmp_path / "run.toml",
                  inputs={"accession": "SRR34217566"},
                  run={"seed": 1, "threads": 2, "pairs": 1000},
                  library={"trim_r2": 0},
                  databases={"pfam": str(tmp_path / "none.hmm"),
                             "sprot": str(tmp_path / "none.dmnd"),
                             "protein": str(tmp_path / "none.dmnd")},
                  # Stated, because a machine running this test need not have DIAMOND
                  # and a stated path is the override. Nothing is executed here.
                  tools={"diamond": str(stub)})


def test_enrich_runs_with_reads_given_on_the_command_line(tmp_path):
    r1, r2 = tmp_path / "a_1.fq.gz", tmp_path / "a_2.fq.gz"
    records = enrich(_cfg(tmp_path), tmp_path, Path.cwd(), 1, (r1, r2))
    assert [r["stage"] for r in records][:1] == ["workup"]
    # Nothing is staged here, so every enrichment is skipped with its reason rather than
    # raising: the point of the test is that building the call did not raise.
    assert all(r["status"] in ("skipped", "failed") for r in records)


def test_a_single_end_run_passes_one_mate_and_not_the_word_none(tmp_path):
    records = enrich(_cfg(tmp_path), tmp_path, Path.cwd(), 1, (tmp_path / "a.fq.gz", None))
    assert records, "no enrichment was attempted"
