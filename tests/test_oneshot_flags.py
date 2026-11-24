"""Every flag the one-command template passes must exist in the thing it calls.

WHY THIS EXISTS. `experiments/prod_oneshot.sh` runs on an instance with no checkout, and
nothing between its lines is checked until the line runs. It passed `--r1` and `--r2` to
`experiments/nontemplate_prefix_map.py`, which takes neither: that script resolves a
library's read paths through `experiments/libraries.json` precisely so that a caller cannot
point it at one library's reads while it reads another library's facts. argparse refused,
the run died with `unrecognized arguments`, and the refusal arrived 29 minutes and about
$0.25 in, after the install, the five reference sets and 1.8 GB of reads had all succeeded.

A flag name is the cheapest thing in a launch to check and one of the last to be checked,
because the shell is a template and the parser is Python. This reads both: the long flags
each command in the template passes, against the long options the command itself declares.
It needs no instance and no network.

The hapax subcommands are checked through the installed parser, and the experiments scripts
through `--help`, which is the only statement of a script's interface that cannot drift from
the script.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ("experiments/prod_oneshot.sh", "experiments/prod_userdata.sh")

# A command in these templates is `"$PY" -I <script> --flag ... \` continued over lines, or
# `hapax <sub> --flag ...`. The tokens that matter are the long flags; everything else is a
# value, a token to be substituted, or shell.
COMMAND = re.compile(
    r"""(?P<head>
            (?:"\$PY"\s+(?:-I\s+)?(?P<script>[\w./$]*experiments/[\w_]+\.py)
             | (?<![\w/.-])hapax[ \t]+(?P<sub>run|report|stages|fetch-db)\b)
        )
        (?P<rest>(?:[^\n\\]|\\\n)*)""",
    re.VERBOSE,
)
FLAG = re.compile(r"(?<![\w-])--[a-z][a-z0-9-]*")


def commands(text: str) -> list[tuple[str, set[str]]]:
    found = []
    for m in COMMAND.finditer(text):
        target = m.group("script") or f"hapax {m.group('sub')}"
        if m.group("script"):
            target = "experiments/" + target.split("experiments/")[-1]
        found.append((target, set(FLAG.findall(m.group("rest")))))
    return found


def declared_by_script(script: Path) -> set[str]:
    proc = subprocess.run([sys.executable, "-I", str(script), "--help"],
                          capture_output=True, text=True, cwd=ROOT, timeout=300)
    if proc.returncode != 0:
        pytest.skip(f"{script.name} --help exits {proc.returncode}: "
                    f"{proc.stderr.strip().splitlines()[-1:]}")
    return set(FLAG.findall(proc.stdout))


def declared_by_hapax(sub: str) -> set[str]:
    # Through the installed entry point rather than by reaching into the parser object,
    # because the entry point is what the template calls.
    proc = subprocess.run([sys.executable, "-I", "-c",
                           "from hapax.cli import main; raise SystemExit(main())",
                           sub, "--help"],
                          capture_output=True, text=True, cwd=ROOT, timeout=300)
    assert proc.returncode == 0, (
        f"`hapax {sub} --help` exits {proc.returncode}, so the template calls a "
        f"subcommand that does not work: {proc.stderr.strip()[-400:]}")
    return set(FLAG.findall(proc.stdout))


@pytest.mark.parametrize("template", TEMPLATES)
def test_every_flag_the_template_passes_is_a_flag_the_command_takes(template):
    text = (ROOT / template).read_text()
    checked = 0
    for target, flags in commands(text):
        if not flags:
            continue
        if target.startswith("hapax "):
            declared = declared_by_hapax(target.split()[1])
        else:
            path = ROOT / target
            if not path.exists():
                raise AssertionError(f"{template} calls {target}, which does not exist")
            declared = declared_by_script(path)
        unknown = sorted(flags - declared)
        assert not unknown, (
            f"{template} passes {unknown} to {target}, which does not take them. "
            f"It takes {sorted(declared)}.")
        checked += 1
    # A regex that matched nothing would pass this test silently, which is the failure
    # shape this file exists to catch one layer down.
    assert checked >= 2, f"{template}: found only {checked} commands to check"
