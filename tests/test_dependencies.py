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
    # src/ is the product since PLAN.md step 3, and it was the first place an import
    # could be added without this tripwire seeing it: the package is not under
    # experiments/ and not a test, so for one commit `jsonschema` was declared and, as
    # far as this file could tell, imported nowhere.
    return sorted([ROOT / "hapax.py", *(ROOT / "experiments").rglob("*.py"),
                   *(ROOT / "src").rglob("*.py"), *(ROOT / "tests").glob("*.py")])


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


#: The conda environment's pins, and TOOLS.md's. They are two statements of one fact.
#:
#: `environment.yml` arrived with PLAN.md step 2, because "the whole install is
#: `conda env create` plus `pip install .`" needs a file that asks for the tools rather
#: than a table that records what somebody happened to install. That makes two places a
#: version is written, which is the shape of every defaulted sample fact this repo has
#: had to unpick: the stage list was typed in three places and the launcher's copy is the
#: one that ran. So the two are compared here, and a pin moved in one place fails the
#: suite rather than letting a run resolve to a version no committed number was measured
#: with.
def _environment_pins() -> dict[str, str]:
    pins = {}
    for line in (ROOT / "environment.yml").read_text().splitlines():
        line = line.split("#")[0].strip()
        if line.startswith("- ") and "=" in line:
            name, _, version = line[2:].partition("=")
            pins[name.strip()] = version.strip()
    return pins


#: environment.yml name -> the name TOOLS.md's table uses for the same tool.
TOOLS_MD_NAME = {"megahit": "MEGAHIT", "diamond": "DIAMOND", "kraken2": "Kraken2",
                 "minimap2": "minimap2", "blast": "BLAST+"}


def test_the_conda_pins_match_tools_md():
    pins = _environment_pins()
    rows = (ROOT / "TOOLS.md").read_text()
    wrong = []
    for conda_name, md_name in TOOLS_MD_NAME.items():
        version = pins.get(conda_name)
        assert version, f"{conda_name} is not pinned in environment.yml"
        # The table's row for the tool, as `| Tool | **version** or version |`.
        for row in rows.splitlines():
            cells = [c.strip() for c in row.split("|")]
            if len(cells) > 3 and cells[1] == md_name:
                recorded = cells[2].replace("*", "")
                if version not in recorded:
                    wrong.append(f"{md_name}: environment.yml asks for {version}, "
                                 f"TOOLS.md pins {recorded}")
                break
        else:
            wrong.append(f"{md_name} is in environment.yml and not in TOOLS.md's table")
    assert not wrong, "\n  " + "\n  ".join(wrong)
