"""A stage must be told where its external binaries are, never left to a script default.

WHY THIS EXISTS. Eight committed scripts declare `--diamond-bin` or `--classify-bin` with
a default that is one container's layout: `/home/user/tools/diamond`. The defaults are
grandfathered, because every committed command still reproduces with them. What must not
happen is the package calling such a script without the flag, which is a path written
where a parameter goes — `hapax.tools` exists to resolve it.

Measured: `divergent_contig_workup.py` was started without `--diamond-bin` and failed with
`FileNotFoundError: /home/user/tools/diamond` on an EC2 instance, 11 stages into a paid
run, after the whole read cascade had finished. `prepare.split` had the rule and `enrich`
did not, which is the shape of every lesson in this repo's conventions: the second writer
did not know which condition the first was fixing.

A static read of the declared scripts and of the package's own call sites. No binaries, no
network.
"""

from __future__ import annotations

import re
from pathlib import Path

from hapax import scripts

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src/hapax"

#: A `--something-bin` flag whose default is an ABSOLUTE path, i.e. one machine's layout.
#: A flag defaulting to a bare name is resolved on PATH and is not this fault.
BIN_FLAG = re.compile(
    r"""add_argument\(\s*["'](?P<flag>--[a-z0-9-]*-bin)["'][^)]*?default\s*=\s*"""
    r"""(?:Path\()?\s*["'](?P<default>/[^"']+)["']""",
    re.DOTALL)


def test_the_package_passes_every_absolute_default_binary_flag():
    package_source = "\n".join(p.read_text() for p in sorted(PACKAGE.glob("*.py")))
    offenders = []
    for name in scripts.declared():
        path = ROOT / "experiments" / name
        if not path.exists():
            continue
        source = path.read_text()
        for m in BIN_FLAG.finditer(source):
            flag, default = m.group("flag"), m.group("default")
            if f'"{flag}"' not in package_source:
                offenders.append(f"{name} defaults {flag} to {default} and no package "
                                 f"module passes it")
    assert not offenders, "\n".join(offenders)


def test_the_resolver_is_what_supplies_them():
    """Passing a path typed in the package would move the fault rather than remove it."""
    package_source = "\n".join(p.read_text() for p in sorted(PACKAGE.glob("*.py")))
    for flag in re.findall(r'"(--[a-z0-9-]*-bin)"', package_source):
        name = flag[2:-4].replace("-", "_")
        # Either form of the resolver: `resolve` for a stage that asks directly, or
        # `resolve_all`, which `cmd_run` uses so that every missing binary arrives in one
        # refusal rather than one per re-run.
        direct = f'tools.resolve(cfg, "{name}")' in package_source
        through_all = ("resolve_all" in package_source
                       and f'"{name}"' in package_source)
        assert direct or through_all, (
            f"{flag} is passed by the package but not from hapax.tools, which is the one "
            f"place that looks for a binary and refuses rather than guessing")
