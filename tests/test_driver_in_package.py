"""The driver is `hapax.driver`, run in-process by `hapax run`; `hapax.py` is a shim.

PLAN.md step 7. What these hold is that the move changed where the code lives and nothing
a caller sees: the committed `python hapax.py ...` command lines run the package's driver,
a module loaded from `hapax.py` by path still exposes the driver's tables (the EC2 launcher
reads its stage list that way), and the package entry point no longer starts a second
program to run it.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

from hapax import cli, driver

ROOT = Path(__file__).resolve().parents[1]


def test_the_shim_runs_the_package_driver_from_anywhere(tmp_path):
    """Run from another directory and isolated, so neither the cwd nor an installed
    `hapax` can decide which driver answers."""
    r = subprocess.run([sys.executable, "-I", str(ROOT / "hapax.py"), "run", "--help"],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "usage: hapax run" in r.stdout


def test_a_module_loaded_from_the_shim_carries_the_drivers_tables():
    spec = importlib.util.spec_from_file_location("hapax_shim_probe", ROOT / "hapax.py")
    shim = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(shim)
    assert shim.PRODUCTION_STAGES == driver.PRODUCTION_STAGES
    assert shim.STAGE_ORDER == driver.STAGE_ORDER


def test_the_manifest_hashes_the_driver_not_the_shim():
    src = (ROOT / "src/hapax/driver.py").read_text()
    assert 'sha256(ROOT / "hapax.py")' not in src, (
        "hashing the shim would pin nothing: the code that ran is src/hapax/driver.py")
    assert driver.DRIVER == (ROOT / "src/hapax/driver.py").resolve()


def test_hapax_run_does_not_start_a_second_program_for_the_driver():
    src = Path(cli.__file__).read_text()
    assert "subprocess" not in src
    assert cli.legacy_driver() is driver


def test_a_driver_refusal_comes_back_as_a_status(tmp_path, capsys):
    rc = cli._drive(["run", str(tmp_path / "no_1.fq"), str(tmp_path / "no_2.fq"),
                     "--seed", "1", "--out", str(tmp_path / "out")])
    assert rc != 0
