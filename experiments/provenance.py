"""What a run manifest needs to identify the code that produced it, imports included.

WHY THIS FILE EXISTS. `subtract_known.provenance(script)` records the git SHA, a dirty
flag and `script_sha256`, which is a content hash of the calling script only. 47 scripts
import their mapping, Kraken2, filtering and statistics code from `subtract_known.py`. When
that shared code changes, their `script_sha256` does not change, and the git SHA is the
only thing left to tell two runs apart. A `-dirty` SHA cannot do that. So the hash this
repo relies on to pin code "whatever the tree looked like" has never pinned the code that
does most of the work. Consolidating more code into shared modules makes that hole bigger,
which is why this module comes with the consolidation and not after it.

This keeps the three existing keys with the same meaning, so readers of manifests are
unaffected. It adds two keys:

- `modules_sha256`: a content hash of every module loaded from this repo at the moment
  the manifest is written, keyed by path relative to the repo root. That is the code the
  run could actually have called.
- `code_sha256`: one hash over the script and all of those modules. Two runs with equal
  `code_sha256` ran the same Python from this repo. It says nothing about external tools
  or their indexes; those belong in `loud.Consulted`.

It is called when the manifest is written, like the original. Editing a module while a
run is in flight therefore stamps the new hash onto the old code's output, as CLAUDE.md
already warns for scripts. And unlike the original, it raises instead of writing
"unknown": a manifest that cannot say which code ran should not be written.
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ProvenanceError(RuntimeError):
    """Raised instead of writing a manifest that cannot identify its code."""


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)


def git_sha() -> str:
    """HEAD, with "-dirty" appended when tracked files differ from it."""
    r = _git("rev-parse", "HEAD")
    if r.returncode != 0:
        raise ProvenanceError(f"git rev-parse failed in {ROOT}: {r.stderr.strip()}")
    return r.stdout.strip() + ("-dirty" if tree_dirty() else "")


def tree_dirty() -> bool:
    r = _git("diff", "--quiet", "HEAD")
    if r.returncode not in (0, 1):
        raise ProvenanceError(f"git diff failed in {ROOT}: {r.stderr.strip()}")
    return r.returncode == 1


def local_modules(exclude: Path | None = None) -> dict[str, str]:
    """sha256 of every module currently imported from inside this repo, by relative path.

    Read from `sys.modules`, so it is what this process has loaded, not what it might
    import later. Call it when the manifest is written, after the work is done.
    """
    out = {}
    for mod in list(sys.modules.values()):
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        p = Path(f).resolve()
        if p.suffix != ".py" or ROOT not in p.parents:
            continue
        if exclude is not None and p == Path(exclude).resolve():
            continue
        out[str(p.relative_to(ROOT))] = _sha(p)
    return dict(sorted(out.items()))


def provenance(script, git_sha_override: str | None = None) -> dict:
    """git_sha, tree_dirty and script_sha256 as before, plus modules_sha256 and code_sha256.

    `script_sha256` keeps its 16-character truncation so existing manifests and the
    aggregators that compare it across seeds read the new ones unchanged.

    WHY `git_sha_override` EXISTS, added 2026-11-23. An EC2 instance has no checkout and no
    git binary - shipping a code payload rather than a working tree is the point - so
    `git_sha()` there does not fail to find a commit, it fails to find *git*, and this
    module raises rather than writing "unknown", which is right. The only true answer on an
    instance is the commit its launcher built the payload from, and a launcher knows it. So
    it is passed in, and `git_sha_source` records which of the two answered, because "the
    launcher said so" and "the repo said so" are different strengths of claim and a field
    that silently means either is worse than no field. The content hashes are unaffected:
    they pin the code whatever the commit says, which is why a launcher-supplied commit is
    safe here and `-dirty` is still refused below.
    """
    script = Path(script).resolve()
    if not script.exists():
        raise ProvenanceError(f"script {script} does not exist")
    mods = local_modules(exclude=script)
    whole = hashlib.sha256()
    whole.update(_sha(script).encode())
    for rel, h in mods.items():
        whole.update(f"\0{rel}\0{h}".encode())
    if git_sha_override is None:
        sha, dirty, source = git_sha(), tree_dirty(), "repo"
    else:
        sha = str(git_sha_override).strip()
        if not sha or sha.endswith("-dirty") or sha == "unknown":
            raise ProvenanceError(
                f"git_sha_override is {git_sha_override!r}, which is not a commit. A "
                f"manifest that cannot identify its code is worse than no manifest, and "
                f"'-dirty' names HEAD rather than the code that ran.")
        dirty, source = False, "launcher"
    return {"git_sha": sha, "tree_dirty": dirty, "git_sha_source": source,
            "script_sha256": _sha(script)[:16],
            "modules_sha256": {k: v[:16] for k, v in mods.items()},
            "code_sha256": whole.hexdigest()[:16]}
