"""Resolve a stage's sample-specific references from a committed file, refusing to guess.

WHY THIS EXISTS. The three stages promoted alongside it - the explained-fraction shortlist
over the removed pile, the cross-plant subtraction and the truth-panel check - each need to
know something about the library in front of them: which depositor genome panel is ground
truth for it, which cross-plant panel it can be subtracted against, and at what operating
point. Every one of those values existed as a literal inside a script or a shell template
before this file. A literal is the bug family this repo keeps being caught by: a default
where a refusal belongs. On a machine where one sample's panel exists, a stage with a
default scores the WRONG sample's genomes and the run directory looks complete.

THE THREE ANSWERS THIS LOADER CAN GIVE, AND WHY THERE ARE THREE. "Unknown library" and
"this library has no truth panel" are different facts and must not collapse into one:

  - the library has an entry and the reference is named -> resolve it;
  - the library has an entry and the reference is explicitly `null` -> a MEASURED absence
    with a committed reason. The caller records `skipped` with that reason, which is the
    honest output for Berlin, whose positive control is the spiked ladder and not a panel;
  - the library has no entry at all -> `UnknownLibrary`. The stage refuses to start.

THE FOLD-IN. `experiments/libraries.json` is meant to hold the per-library facts that
decide what gets read (trim_r2, block, read paths). It is not in the committed tree, so
this file holds the stage-reference half separately, and `resolve_library` is the single
place that changes when the two become one. Two files are the interim, not the design.

WHAT THIS DELIBERATELY DOES NOT DO. It does not sniff the library from a read path, and it
does not accept a library name that disagrees with the reads it was handed - that check
belongs to whatever resolves the reads, and guessing here would reintroduce the silence
this file exists to remove.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCES = Path(__file__).resolve().parent / "stage_references.json"


class UnknownLibrary(KeyError):
    """A library with no entry. Raised rather than defaulted, always."""


class MissingReference(FileNotFoundError):
    """A named reference whose file is absent. A refusal, not a zero."""


def load(path: Path | None = None) -> dict:
    p = Path(path) if path else REFERENCES
    if not p.exists():
        raise MissingReference(f"{p} is absent, so no stage reference can be resolved")
    return json.loads(p.read_text())


def libraries(refs: dict | None = None) -> list[str]:
    return sorted((refs or load())["libraries"])


def resolve_library(library: str, refs: dict | None = None) -> dict:
    """The entry for one library, or UnknownLibrary with the list of the ones that exist."""
    refs = refs or load()
    entry = refs["libraries"].get(library)
    if entry is None:
        raise UnknownLibrary(
            f"{library} has no entry in {REFERENCES.name}, so this stage does not know "
            f"which panel is ground truth for it and will not guess. Known libraries: "
            f"{', '.join(libraries(refs))}. Add an entry, or run without this stage.")
    return entry


def resolve(library: str, kind: str, refs: dict | None = None) -> dict:
    """Resolve one reference for one library.

    Returns `{"state": "resolved", ...}` with the panel's own record and its paths made
    absolute, or `{"state": "absent", "reason": ...}` when the entry says this library has
    none. Raises for an unknown library, an unknown panel name, or a named file that is not
    there - three conditions a caller must not be able to read as a zero.
    """
    refs = refs or load()
    entry = resolve_library(library, refs)
    if kind not in ("truth_panel", "crossplant_panel"):
        raise KeyError(f"no such stage reference: {kind}")
    name = entry.get(kind)
    if name is None:
        reason = entry.get(f"{kind}_absent_because")
        if not reason:
            raise MissingReference(
                f"{library}'s {kind} is null with no committed reason beside it. A null "
                f"here is a measured absence and has to say why; an empty one is "
                f"indistinguishable from a forgotten entry.")
        return {"state": "absent", "library": library, "kind": kind, "reason": reason}
    table = {"truth_panel": "truth_panels", "crossplant_panel": "crossplant_panels"}[kind]
    panel = refs.get(table, {}).get(name)
    if panel is None:
        raise UnknownLibrary(
            f"{library} names {kind} '{name}', which is not defined under {table} in "
            f"{REFERENCES.name}")
    out = {"state": "resolved", "library": library, "kind": kind, "panel": name,
           "record": panel, "files": {}}
    for key, value in panel.items():
        if key.startswith("_") or not isinstance(value, str):
            continue
        if value.startswith("results/") or value.startswith("data/"):
            p = ROOT / value
            out["files"][key] = p
    return out


def missing_files(resolved: dict) -> list[str]:
    """Which of a resolved reference's files are absent or empty.

    Empty counts as absent: a zero-byte panel is a failed download, not a panel with
    nothing in it, and this repo has had a stage report a number off one.
    """
    out = []
    for key, p in resolved.get("files", {}).items():
        if not p.exists():
            out.append(f"{key}: {p} is absent")
        elif p.is_file() and p.stat().st_size == 0:
            out.append(f"{key}: {p} is empty (0 bytes)")
    return out


def band(refs: dict | None = None) -> dict:
    """The shipped coverage x identity band for the removal shortlist, with its provenance.

    It is returned whole, including `chosen_after_seeing_the_positive`, so a caller cannot
    take the numbers without the sentence that qualifies them.
    """
    refs = refs or load()
    b = refs["removal_shortlist_band"]
    for key in ("coverage", "identity_percent", "significance_evalue", "measured_on",
                "chosen_after_seeing_the_positive"):
        if key not in b:
            raise MissingReference(f"removal_shortlist_band has no {key}")
    return b
