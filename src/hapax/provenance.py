"""What a manifest written by package code records about the code that wrote it.

WHY THE PACKAGE HAS ITS OWN. `experiments/provenance.py` hashes the calling SCRIPT, because
in `experiments/` a script is the unit that produced a number (CLAUDE.md). Package code is
versioned instead (PLAN.md section 2): a manifest names the package version, the commit, and
a content hash of the module that wrote it, so two runs from a dirty tree can still be told
apart by inspection. It raises rather than writing "unknown", for the same reason the
script version does: a manifest that cannot say which code ran should not be written.

The commit is read from git when there is a checkout, and from `HAPAX_GIT_SHA` only when
git cannot answer, which is the instance case the script version documents.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

from . import __version__

ROOT = Path(__file__).resolve().parents[2]
SHA_ENV = "HAPAX_GIT_SHA"


class ProvenanceError(RuntimeError):
    """Raised instead of writing a manifest that cannot identify its code."""


def _git(*args: str):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    except FileNotFoundError:
        return None


def commit() -> tuple[str, bool]:
    """HEAD and whether tracked files differ from it."""
    head = _git("rev-parse", "HEAD")
    if head is None or head.returncode != 0:
        env = (os.environ.get(SHA_ENV) or "").strip()
        if env:
            return env, False
        raise ProvenanceError(
            f"git cannot name the commit in {ROOT} and {SHA_ENV} is not set; a manifest "
            f"that cannot say which code ran should not be written")
    diff = _git("diff", "--quiet", "HEAD")
    if diff is None or diff.returncode not in (0, 1):
        raise ProvenanceError(f"git diff failed in {ROOT}")
    return head.stdout.strip(), diff.returncode == 1


def record(module_file: str | Path, producer: str) -> dict:
    """The provenance block for a manifest written by `producer` in `module_file`."""
    sha, dirty = commit()
    return {"producer": producer, "hapax_version": __version__,
            "git_sha": sha + ("-dirty" if dirty else ""), "tree_dirty": dirty,
            "module_sha256": hashlib.sha256(Path(module_file).read_bytes()).hexdigest()}
