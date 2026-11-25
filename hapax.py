#!/usr/bin/env python3
"""`python hapax.py ...` - kept so every committed command line still runs.

The driver moved into the package as `hapax.driver` (PLAN.md step 7). This file is a
shim: it puts this checkout's `src/` first on the path, imports the driver from there and
runs it, and when it is loaded as a module by path (as `experiments/ec2_prod.py` does to
read the stage list) it exposes the driver's names.

WHY `src/` IS PUT FIRST AND THIS DIRECTORY TAKEN OFF. Run as a script, Python puts this
directory at the front of `sys.path`, so `import hapax` would find this file instead of
the package, which is the trap `experiments/prod_oneshot.sh` records. And a different
installed `hapax` must not be the one that runs: the driver that runs is the one beside
this file.

Committed run manifests before the move name this file, whose content was the driver;
manifests after it name `src/hapax/driver.py`.
"""
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
sys.path[:] = [p for p in sys.path if Path(p or ".").resolve() != _here]
sys.path.insert(0, str(_here / "src"))
for _name in [m for m in sys.modules if m == "hapax" or m.startswith("hapax.")]:
    _mod = sys.modules[_name]
    if not str(getattr(_mod, "__file__", "") or "").startswith(str(_here / "src")):
        del sys.modules[_name]

from hapax import driver as _driver  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(_driver.main())

globals().update({k: v for k, v in vars(_driver).items() if not k.startswith("__")})
