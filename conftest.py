"""Keep the repo root off `sys.path`, so that `hapax` means the package.

WHY THIS FILE EXISTS. There are two things called hapax in this checkout: the installed
package `src/hapax/`, which is the product, and `hapax.py` at the root, which is the
driver the package still shells out to (PLAN.md step 7 archives it). `python -m pytest`
puts the working directory at the front of `sys.path`, so from the repo root `import
hapax` reaches the SCRIPT and the package is invisible; plain `pytest` does not, so the
same suite passed or failed depending on how it was started.

A latent dependency on how a command was typed is the shape of defect this repo keeps
paying for, so it is removed rather than documented: the root comes off the path here,
once, for every entry point. Nothing needs it. The driver is loaded from its own file by
`tests/legacy.py`, and the scripts under `experiments/` are added to the path by the tests
that use them.

The driver's content is deliberately untouched by any of this. Its `script_sha256` pins
five committed seed manifests, so moving or editing it to resolve a name collision would
damage the one field that says which code produced those numbers.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = str(Path(__file__).resolve().parent)
sys.path[:] = [p for p in sys.path if p not in ("", ".", ROOT)]
