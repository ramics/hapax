"""Load the legacy driver, `hapax.py` at the repo root, from its own file.

NOT A TEST. A helper, because `import hapax` is now ambiguous: the product is the package
`src/hapax/` and the driver is still the script `hapax.py` beside it, so which one that
statement reaches depends on the order of `sys.path`, which depends on which test ran
first. That is the kind of defect this repo keeps paying for, so the name is not used:
the driver is loaded from its path and the package keeps the name.

The driver's content is deliberately untouched by any of this. Its `script_sha256` pins
five committed seed manifests, so editing it to make an import convenient would break the
one field that says which code produced those numbers.

PLAN.md step 7 archives the driver, and this helper goes with it.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def driver():
    """The root hapax.py, as a module, loaded once per process."""
    name = "hapax_legacy_driver"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "hapax.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
