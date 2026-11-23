"""provenance.py: a manifest's code hash must move when imported shared code moves.

Run with `python3 -m pytest tests/test_provenance.py`.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import textwrap
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments"
sys.path.insert(0, str(EXP))

import provenance as pv  # noqa: E402

PROBE = textwrap.dedent("""
    import json, sys
    sys.path.insert(0, {exp!r})
    sys.path.insert(0, {tmp!r})
    import provenance, {mod}
    print(json.dumps(provenance.provenance(__file__)))
""")


def _run(tmp: Path, mod_body: str) -> dict:
    """A script inside the repo tree importing a module inside the repo tree."""
    import json
    (tmp / "probe_shared_mod.py").write_text(mod_body)
    script = tmp / "probe_script.py"
    script.write_text(PROBE.format(exp=str(EXP), tmp=str(tmp), mod="probe_shared_mod"))
    out = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                         check=True)
    return json.loads(out.stdout)


def test_existing_keys_keep_their_meaning():
    p = pv.provenance(EXP / "provenance.py")
    assert set(p) >= {"git_sha", "tree_dirty", "script_sha256"}
    assert p["script_sha256"] == hashlib.sha256(
        (EXP / "provenance.py").read_bytes()).hexdigest()[:16]
    assert p["git_sha"].endswith("-dirty") == p["tree_dirty"]


def test_imported_module_change_moves_code_hash_not_script_hash():
    # Under the repo root, because only repo modules are hashed; gitignored data/ is not
    # needed, a throwaway directory under tests/ is cleaned up below.
    tmp = ROOT / "tests" / "_provenance_probe"
    tmp.mkdir(exist_ok=True)
    try:
        a = _run(tmp, "X = 1\n")
        b = _run(tmp, "X = 2\n")
        assert a["script_sha256"] == b["script_sha256"]
        assert a["code_sha256"] != b["code_sha256"]
        key = "tests/_provenance_probe/probe_shared_mod.py"
        assert key in a["modules_sha256"]
        assert a["modules_sha256"][key] != b["modules_sha256"][key]
        assert "experiments/provenance.py" in a["modules_sha256"]
        assert "tests/_provenance_probe/probe_script.py" not in a["modules_sha256"]
    finally:
        shutil.rmtree(tmp)


def test_missing_script_raises():
    import pytest
    with pytest.raises(pv.ProvenanceError):
        pv.provenance(ROOT / "no_such_script.py")


def test_launcher_supplied_commit_is_marked_as_such():
    """An instance has no checkout, so the launcher's commit is the only true answer there.

    The test that matters is not that the field is carried but that its SOURCE is: "the
    launcher said so" and "the repo said so" are different strengths of claim, and a
    manifest where either could be meant is worse than one with no field. The content
    hashes must not move, because they pin the code whatever the commit says.
    """
    repo = pv.provenance(EXP / "provenance.py")
    given = pv.provenance(EXP / "provenance.py", "0" * 40)
    assert repo["git_sha_source"] == "repo"
    assert given["git_sha_source"] == "launcher"
    assert given["git_sha"] == "0" * 40
    assert given["tree_dirty"] is False
    assert given["code_sha256"] == repo["code_sha256"]
    assert given["script_sha256"] == repo["script_sha256"]


def test_a_dirty_or_empty_override_is_refused():
    """`-dirty` names HEAD, not the code that ran, and this is the one place it could slip
    in without a repo to contradict it."""
    for bad in ("", "   ", "unknown", "abc123-dirty"):
        with pytest.raises(pv.ProvenanceError):
            pv.provenance(EXP / "provenance.py", bad)
