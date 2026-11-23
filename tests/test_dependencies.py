"""Every third-party import in the repo is declared in requirements.txt.

The repo ran for weeks with no dependency file at all, so a fresh container could not run
the suite green unless someone already knew to `pip install mappy`, and the failure it got
instead was an assertion about an exit status, two layers away from the cause. This is the
tripwire that keeps the declaration honest: an import added without a line in
requirements.txt fails here, by name, rather than in the next fresh container.

Static (an AST walk), so it needs nothing installed and runs in the same container whose
missing packages it is about to name.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Import name -> distribution name, where they differ.
DIST = {"RNA": "ViennaRNA", "sklearn": "scikit-learn"}

#: Top-level names that are this repo's own code, not packages. `python` is the archived
#: codon assembler's package (experiments/codon_assembler/python/).
LOCAL = {"python"}


def _sources() -> list[Path]:
    return sorted([ROOT / "hapax.py", *(ROOT / "experiments").rglob("*.py"),
                   *(ROOT / "tests").glob("*.py")])


def _third_party() -> dict[str, list[str]]:
    local = LOCAL | {p.stem for p in ROOT.rglob("*.py") if ".git" not in p.parts}
    found: dict[str, list[str]] = {}
    for path in _sources():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for name in names:
                top = name.split(".")[0]
                if top not in sys.stdlib_module_names and top not in local:
                    found.setdefault(top, []).append(str(path.relative_to(ROOT)))
    return found


def _declared() -> set[str]:
    out = set()
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        line = line.split("#")[0].strip()
        if line:
            out.add(re.split(r"[<>=!~\[; ]", line, maxsplit=1)[0].lower())
    return out


def test_every_import_is_declared():
    declared = _declared()
    missing = {mod: paths for mod, paths in _third_party().items()
               if DIST.get(mod, mod).lower() not in declared}
    assert not missing, ("imported but not in requirements.txt:\n  " + "\n  ".join(
        f"{DIST.get(m, m)} (import {m}): {', '.join(sorted(set(p))[:3])}"
        for m, p in sorted(missing.items())))


def test_nothing_declared_is_unused():
    """A declared package nothing imports is a stale line, and a stale pin reads as a
    provenance claim nobody is making."""
    used = {DIST.get(m, m).lower() for m in _third_party()}
    assert not (_declared() - used), f"declared but never imported: {_declared() - used}"
