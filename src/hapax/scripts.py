"""The one place the package starts a committed `experiments/` script.

WHY THESE STAY PROCESSES. Imo's call on 2026-11-24: the science stages remain the committed
scripts, byte-identical, run as child processes. Two reasons, both about provenance. A
script's `script_sha256` is what ties a committed result to the code that wrote it, so
lifting the code would leave every committed manifest pointing at a file the tool no longer
runs. And `experiments/provenance.py` records `modules_sha256` from `sys.modules` when the
manifest is written, so a stage run in-process would stamp every earlier stage's modules
into its own `code_sha256`. The process boundary is what scopes that hash to the stage.

WHAT THIS MODULE ADDS, which is PLAN.md 1.7. A subprocess hides what a script depends on.
Before this, `hapax.prepare`, `hapax.enrich` and `hapax.report` each built their own command
line with `sys.executable` and a path typed in place, and nothing said that
`classify_contigs.py` runs `subtract_known.py`'s code, or which file of it. Here every
script the package may start is declared, a name not declared is refused, and each call
returns the script's hash and the hash of every repo module it imports, transitively, read
from the files before the process starts. So a run's record names the code each stage
could have run, and `tests/test_scripts.py` fails when a package module starts an
`experiments/` script any other way.

The import closure is found by reading the source, not by importing it, because importing a
stage script runs its module-level code. A dynamic import (`importlib`, `__import__`) is
not seen; the closure is a floor on what the script can load, and says so.
"""
from __future__ import annotations

import ast
import hashlib
import subprocess
import sys
from pathlib import Path

from .stages import R

#: Scripts the package starts outside the stage registry: the assembly and the contig
#: split, which `hapax.prepare` runs before any stage. Every other script comes from the
#: registry's own `script` field, so a stage cannot be declared with a script this module
#: would refuse.
PREPARE_SCRIPTS = ("path_settings.py", "classify_contigs.py")


class ScriptError(RuntimeError):
    """A script the package has not declared, or one whose code is not where it should be."""


def declared() -> tuple[str, ...]:
    from_stages = [s.script for s in R.stages.values() if s.script]
    return tuple(dict.fromkeys([*from_stages, *PREPARE_SCRIPTS]))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _local(name: str, exp: Path) -> Path | None:
    head = name.split(".")[0]
    for cand in (exp / f"{head}.py", exp / head / "__init__.py"):
        if cand.exists():
            return cand
    return None


def import_closure(script: Path, exp: Path) -> dict[str, str]:
    """sha256 of every `experiments/` module the script imports, transitively, by name.

    Only top-level module names are resolved against `exp`, which is how every script here
    imports its neighbours (`sys.path.insert(0, <its own directory>)` then a bare import).
    """
    seen: dict[str, str] = {}
    todo = [script]
    while todo:
        path = todo.pop()
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            for name in names:
                found = _local(name, exp)
                if found is None or found == script:
                    continue
                rel = str(found.relative_to(exp.parent))
                if rel not in seen:
                    seen[rel] = _sha(found)
                    todo.append(found)
    return dict(sorted(seen.items()))


def command(repo: Path, name: str, args: list, isolated: bool = True) -> list[str]:
    """The argv for one declared script, by this interpreter, isolated (`-I`) by default."""
    if name not in declared():
        raise ScriptError(
            f"{name} is not a script the package declares. Add it to a stage's `script` "
            f"field or to PREPARE_SCRIPTS, so the run's record can name its code.")
    path = Path(repo) / "experiments" / name
    if not path.exists():
        raise ScriptError(f"{path} does not exist")
    return [sys.executable, *(["-I"] if isolated else []), str(path), *map(str, args)]


def describe(repo: Path, name: str) -> dict:
    """What a call to `name` could run: the script's hash and its import closure."""
    path = Path(repo) / "experiments" / name
    return {"script": f"experiments/{name}", "script_sha256": _sha(path),
            "imports_sha256": import_closure(path, path.parent),
            "imports_note": "read from the source before the process started; a dynamic "
                            "import is not seen, so this is a floor on what could load"}


def run(repo: Path, name: str, args: list, *, isolated: bool = True,
        cwd: Path | None = None) -> tuple[subprocess.CompletedProcess, dict]:
    """Start one declared script and return the finished process and its code record.

    The code record is taken BEFORE the process starts, because CLAUDE.md's rule about
    hashing applies here too: a hash taken afterwards would describe a file that may have
    been edited while the script ran.
    """
    cmd = command(repo, name, args, isolated)
    code = describe(repo, name)
    proc = subprocess.run(cmd, cwd=cwd or repo, capture_output=True, text=True)
    script_at = cmd.index(str(Path(repo) / "experiments" / name))
    return proc, {**code, "argv": cmd[script_at + 1:], "returncode": proc.returncode}
