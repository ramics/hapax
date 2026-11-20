"""Shared plumbing for experiments: run manifests, git SHA, hapax invocation.

Every experiment writes a manifest beside its results (seed, git SHA,
parameters, timings, tool versions) so that any reported number can be
reproduced from a committed script plus a seed.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
HAPAX = ROOT / "target" / "release" / "hapax"


def git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def ensure_release_build() -> None:
    """Benchmarks come from release builds, never debug."""
    subprocess.check_call(["cargo", "build", "--release", "--quiet"], cwd=ROOT)


def run_hapax(reads: Path, out: Path, args: list[str]) -> dict:
    t = time.perf_counter()
    subprocess.run([str(HAPAX), "assemble", "-i", str(reads), "-o", str(out), *args], check=True, stderr=subprocess.DEVNULL)
    wall = time.perf_counter() - t
    stats = json.loads((out / "stats.json").read_text())
    stats["wall"] = wall
    return stats


def tool_versions() -> dict:
    v = {"python": platform.python_version(), "platform": platform.platform()}
    try:
        v["rustc"] = subprocess.check_output(["rustc", "--version"], text=True).strip()
    except Exception:
        pass
    # the external assemblers live in a micromamba env, not on PATH: look there
    # too, or every manifest records "not installed" for tools that did run
    env = Path(os.environ.get("HAPAX_TOOLS", "/home/user/mamba/envs/asm/bin"))
    for tool in ("megahit", "metaspades.py", "spades.py", "haploflow", "art_illumina", "minimap2", "samtools", "plass", "penguin"):
        p = shutil.which(tool) or (str(env / tool) if (env / tool).exists() else None)
        if not p:
            v[tool] = "not installed"
            continue
        ver = ""
        for flag in ("--version", "-v"):
            try:
                ver = subprocess.run([p, flag], capture_output=True, text=True, timeout=30).stdout.strip().splitlines()[0]
            except Exception:
                ver = ""
            if ver:
                break
        v[tool] = f"{p} {ver}".strip()
    return v


def write_manifest(out: Path, seed, params: dict, timings: dict, extra: dict | None = None) -> None:
    m = {
        "seed": seed,
        "git_sha": git_sha(),
        "params": params,
        "timings_s": timings,
        "tools": tool_versions(),
        "command": " ".join(sys.argv),
        "date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if extra:
        m.update(extra)
    (out / "manifest.json").write_text(json.dumps(m, indent=2))
