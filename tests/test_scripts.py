"""Every committed script the package starts, it starts through `hapax.scripts`.

PLAN.md 1.7: a subprocess hides what a script depends on. The stage scripts stay processes
on purpose (Imo, 2026-11-24, so each stage's `code_sha256` is scoped to its own code), so
what these tests hold is the other half: that a call is declared, refused when it is not,
and recorded with the code it could run.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from hapax import scripts

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "src/hapax"


def test_no_package_module_builds_a_script_command_of_its_own():
    """The tripwire. Before this, three modules each typed `sys.executable, "-I", repo /
    "experiments/<name>.py"`, and none recorded what the script imported."""
    offenders = []
    for f in sorted(PKG.glob("*.py")):
        if f.name == "scripts.py":
            continue
        for i, line in enumerate(f.read_text().splitlines(), 1):
            # A command line, not a path: `_render_check` loads make_web.py as a module
            # to call one function, which is an import and is allowed.
            if re.search(r'(-I"|sys\.executable|\bpy\b,).*experiments/\w+\.py', line):
                offenders.append(f"{f.name}:{i}: {line.strip()}")
    assert not offenders, f"start experiments/ scripts through hapax.scripts: {offenders}"


def test_every_declared_script_exists_and_its_imports_resolve():
    for name in scripts.declared():
        rec = scripts.describe(ROOT, name)
        assert len(rec["script_sha256"]) == 64
        for rel in rec["imports_sha256"]:
            assert (ROOT / rel).exists(), f"{name} imports {rel}, which is not there"


def test_the_closure_is_transitive():
    """`known_strain_divergence.py` reaches `subtract_known.py` only through
    `public_health_profile.py`, which is the dependency a subprocess hid."""
    rec = scripts.describe(ROOT, "known_strain_divergence.py")
    assert "experiments/public_health_profile.py" in rec["imports_sha256"]
    assert "experiments/subtract_known.py" in rec["imports_sha256"]


def test_an_undeclared_script_is_refused():
    with pytest.raises(scripts.ScriptError, match="not a script the package declares"):
        scripts.command(ROOT, "shuffle_control.py", [])


def test_a_call_records_its_code_before_it_starts(tmp_path):
    proc, code = scripts.run(ROOT, "architecture_call.py", ["--help"])
    assert proc.returncode == 0, proc.stderr
    assert code["script"] == "experiments/architecture_call.py"
    assert code["argv"] == ["--help"]
    assert code["returncode"] == 0
