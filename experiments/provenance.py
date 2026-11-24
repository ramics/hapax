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
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ProvenanceError(RuntimeError):
    """Raised instead of writing a manifest that cannot identify its code."""


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(*args: str) -> subprocess.CompletedProcess | None:
    """None when there is no git at all, which is an instance rather than a checkout."""
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    except FileNotFoundError:
        return None


#: For a machine with NO CHECKOUT. An EC2 instance gets a tarball of the code and no git,
#: so every script importing this module raised there - measured, on the first production
#: launch, where it crashed the step that decides the mate-2 trim and made the shell report
#: a refusal that had not happened. "unknown" is what this module exists to refuse, and
#: there is a better answer available: the launcher keys the payload by the commit and
#: refuses to launch when the shipped files are dirty, so the commit IS known and the tree
#: IS clean, from outside. The launcher passes it in this variable.
#:
#: It is read only when git cannot answer. A machine with a checkout ignores it, so it
#: cannot be used to dress a dirty tree up as a clean one locally.
SHA_ENV = "HAPAX_GIT_SHA"


def _env_sha() -> str | None:
    v = (os.environ.get(SHA_ENV) or "").strip()
    return v or None


def git_sha() -> str:
    """HEAD, with "-dirty" appended when tracked files differ from it."""
    r = _git("rev-parse", "HEAD")
    if r is None or r.returncode != 0:
        env = _env_sha()
        if env:
            return env
        detail = "git is not installed" if r is None else r.stderr.strip()
        raise ProvenanceError(
            f"git rev-parse failed in {ROOT}: {detail}. On a machine with no checkout, set "
            f"{SHA_ENV} to the commit the code was shipped from; a manifest that cannot say "
            f"which code ran should not be written.")
    return r.stdout.strip() + ("-dirty" if tree_dirty() else "")


def tree_dirty() -> bool:
    r = _git("diff", "--quiet", "HEAD")
    if r is None or r.returncode not in (0, 1):
        if _env_sha():
            # No checkout to compare against. The launcher's refusal to ship a dirty tree
            # is what stands behind this, and the SHA above says which commit.
            return False
        detail = "git is not installed" if r is None else r.stderr.strip()
        raise ProvenanceError(f"git diff failed in {ROOT}: {detail}")
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
        # A virtualenv created INSIDE the repo root puts third-party code under it, and
        # "every module loaded from this repo" then quietly includes site-packages. The
        # first production run recorded `venv/lib/python3.12/site-packages/
        # _distutils_hack/__init__.py` in its manifests beside this repo's own modules,
        # which makes `code_sha256` change when a dependency changes - the dependency set
        # is `pip_freeze.txt` and `loud.Consulted`'s job, and mixing the two means neither
        # hash answers its own question.
        if "site-packages" in p.parts or "dist-packages" in p.parts:
            continue
        if exclude is not None and p == Path(exclude).resolve():
            continue
        out[str(p.relative_to(ROOT))] = _sha(p)
    return dict(sorted(out.items()))


def provenance(script) -> dict:
    """git_sha, tree_dirty and script_sha256 as before, plus modules_sha256 and code_sha256.

    `script_sha256` keeps its 16-character truncation so existing manifests and the
    aggregators that compare it across seeds read the new ones unchanged.
    """
    script = Path(script).resolve()
    if not script.exists():
        raise ProvenanceError(f"script {script} does not exist")
    mods = local_modules(exclude=script)
    whole = hashlib.sha256()
    whole.update(_sha(script).encode())
    for rel, h in mods.items():
        whole.update(f"\0{rel}\0{h}".encode())
    return {"git_sha": git_sha(), "tree_dirty": tree_dirty(),
            "script_sha256": _sha(script)[:16],
            "modules_sha256": {k: v[:16] for k, v in mods.items()},
            "code_sha256": whole.hexdigest()[:16]}


def rel(path) -> str:
    """A path for a manifest: repo-relative where it is in the repo, absolute otherwise.

    WHY THIS IS SHARED. Four scripts recorded their inputs with a bare
    `path.relative_to(ROOT)`, which raises for any path outside the checkout. That held
    while every input was a file under `data/`, and stopped holding the first time a run
    directory was somewhere else: a stage that had already done its work died while
    writing down what it had used. The fix is not a try/except at each site, because the
    question - what do I call this path in a record - has one answer.

    The repo-relative form is kept where it applies, because it is what makes two
    machines' manifests comparable; an absolute path outside the repo is reported as it
    is, since there is nothing to make it relative to.
    """
    p = Path(path).resolve()
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
