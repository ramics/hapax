"""No script the package runs may work out its own commit.

WHY THIS EXISTS, AND WHY IT IS WIDER THAN THE DRIVER'S VERSION. `tests/test_hapax_stages.py`
already fails if `hapax.py` grows a private `git()`: that copy wrote `git_sha: ""` and
`tree_dirty: false` into every manifest of the first production run, because `bool("")`
is False and the instance has a payload tarball and no git. The same copy existed in five
of the scripts the driver and the package run, where it caught every exception and
answered `"unknown"` and `True` - so the committed outputs of the second and third
production libraries carry `git_sha: "unknown-dirty"` and a positive claim that the tree
was dirty, about a tree nobody looked at. `architecture_call.py`'s copy did not even
catch: it raised FileNotFoundError and lost the stage after the per-gene workup it reads
had already run.

`experiments/provenance.py` is the one implementation. It reads HAPAX_GIT_SHA where there
is no checkout, which the launcher exports from the commit it built the payload from, and
it raises rather than answering where it can say nothing.

This is the shuffle tripwire applied to provenance, one layer out from the driver.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

#: `git` as the program of a subprocess call, in any of the spellings used here.
GIT_CALL = re.compile(r'\[\s*"git"|"git"\s*,\s*"|\(\s*"git"')


def test_no_declared_script_shells_out_to_git():
    from hapax import scripts

    offenders = []
    for name in sorted(scripts.declared()):
        path = ROOT / "experiments" / name
        if not path.exists():
            continue
        if GIT_CALL.search(path.read_text()):
            offenders.append(name)
    assert not offenders, (
        f"{offenders} work out their own commit. Import experiments/provenance.py, which "
        f"reads HAPAX_GIT_SHA on a machine with no checkout and raises where it cannot "
        f"say which code ran.")


def test_provenance_module_is_what_they_import():
    """The sweep above is satisfiable by writing no provenance at all, so check one."""
    src = (ROOT / "experiments/subtract_known.py").read_text()
    assert "import provenance as _provenance" in src
    assert "git_sha = _provenance.git_sha" in src
