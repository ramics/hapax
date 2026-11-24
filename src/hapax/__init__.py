"""hapax: subtract what is identifiable in a metagenome, then characterise what is left.

The product lives here and only here. `experiments/` keeps the research scripts and their
results and is never imported by this package; the package shells out to the ones it still
needs, and PLAN.md step 7 replaces those calls with imports.

The version is a provenance field. The per-script content hash rule (CLAUDE.md) applies to
`experiments/`, where a script is the unit that produced a number; product code is
versioned and edited normally, and a number in `results/` names the package version that
produced it.
"""
from __future__ import annotations

__all__ = ["__version__"]

try:                                    # installed
    from importlib.metadata import version as _version

    __version__ = _version("hapax")
except Exception:                       # a checkout with no install
    # Not a fallback to a plausible value: it says which of the two cases this is, so a
    # manifest written from a checkout cannot be read as one written from an install.
    __version__ = "0.3.0+checkout"
