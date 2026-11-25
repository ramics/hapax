"""Whether each stage's `--out` is a directory or a file, checked against the script.

WHY THIS TEST EXISTS. `--out` means a directory to some of these scripts and a file to
others, and the flag looks the same from the driver either way. Handed a file path, a
script that means a directory ran `out.mkdir(parents=True)` on it and wrote its own
per-seed filename inside, so the stage exited 0, the run directory had
`public_health_profile_s1.json/public_health_profile_s1.json`, and the only thing that
noticed was the driver failing to read its own expected output as JSON. Two of the stages
adopted in PLAN.md step 6 were wired that way.

HOW IT CHECKS. The script says which it means: `args.out.mkdir(...)` followed by
`args.out / "<name>"` is a directory, and `out.write_text` or `out.parent.mkdir` is a
file. The driver says what it passes: a `--out` whose expression ends in a filename
literal is a file. Both are read off the source, so a new stage, or a script that changes
its mind, fails here rather than in a run directory nobody inspects.

WHAT IT CANNOT CATCH. A script that takes a directory and writes a name the driver does
not expect. That is what the registry's declared outputs are for, and the end-to-end job
is what reads them.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "src/hapax/driver.py"
EXPERIMENTS = ROOT / "experiments"


def _driver_out_flags() -> dict[str, str]:
    """Stage method name -> the `--out` expression the driver passes it."""
    src = DRIVER.read_text()
    found: dict[str, str] = {}
    stage = None
    for line in src.splitlines():
        m = re.match(r"    def (\w+)\(self\)", line)
        if m:
            stage = m.group(1)
        hit = re.search(r'"--out", str\(([^)]*\)?[^)]*)\)\]?', line)
        if hit and stage:
            found.setdefault(stage, hit.group(1).strip())
    return found


def _script_of(stage: str) -> Path | None:
    """The experiment script a stage calls, read off its own `self.call`."""
    src = DRIVER.read_text()
    block = re.search(rf"    def {stage}\(self\).*?(?=\n    def |\nclass |\Z)", src,
                      re.S)
    if not block:
        return None
    name = re.search(r'self\.call\(\s*"[^"]+",\s*"([^"]+\.py)"', block.group(0))
    return EXPERIMENTS / name.group(1) if name else None


def _script_means(path: Path) -> str:
    """"dir" or "file", from what the script does with its own --out."""
    src = path.read_text()
    if re.search(r"\ba(?:rgs)?\.out\.mkdir\(", src):
        return "dir"
    if re.search(r"\ba(?:rgs)?\.out\.(?:write_text|open)\b", src) or \
            re.search(r"\ba(?:rgs)?\.out\.parent\.mkdir\(", src):
        return "file"
    return "unknown"


STAGES = sorted(_driver_out_flags())


@pytest.mark.parametrize("stage", STAGES)
def test_the_driver_passes_out_in_the_shape_the_script_means(stage):
    script = _script_of(stage)
    if script is None or not script.exists():
        pytest.skip(f"{stage} does not call an experiment script by a literal name")
    means = _script_means(script)
    if means == "unknown":
        pytest.skip(f"{script.name} does not say what its --out is in a form read here")
    passed = _driver_out_flags()[stage]
    looks_like_a_file = bool(re.search(r'\.json"?\)?\s*$', passed))
    assert (means == "file") == looks_like_a_file, (
        f"{stage} passes --out {passed!r} and {script.name} treats --out as a {means}. "
        f"A script that means a directory and is handed a file path makes a directory of "
        f"that name and writes inside it, which exits 0 and produces nothing the run "
        f"directory can find.")


def test_every_stage_that_passes_out_was_examined():
    """A stage list that silently shrinks would make this suite pass by testing nothing."""
    assert len(STAGES) >= 10, STAGES
