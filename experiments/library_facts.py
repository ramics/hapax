#!/usr/bin/env python3
"""Per-library constants, from one committed file, with a refusal instead of a default.

Why. A trim length and a block size are facts about a sequencing run. They were written
into roughly fifteen scripts as `default=14` and into five more as a bare `14` inside the
call. That is one library's biology applied to whatever FASTQ the next caller points at,
silently, in the direction of less recognition: on a library with no non-template prefix,
trimming 14 nt deletes real sequence from every mate and every downstream rate drops.

The shape of the fix is the one this repo keeps arriving at. There is no default. Ask for
a library by accession and get its values, or state them explicitly; a library this file
does not know raises, and the message says what has to be added. A value that is absent
is absent, never another library's.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FACTS = ROOT / "experiments/libraries.json"


class UnknownLibrary(KeyError):
    """Raised rather than falling back. A guess here is silent data loss."""


def library(accession: str) -> dict:
    d = json.loads(FACTS.read_text())["libraries"]
    if accession not in d:
        raise UnknownLibrary(
            f"{accession} is not in {FACTS.relative_to(ROOT)}. Add an entry giving at "
            f"least trim_r2 and block, each with how it was established, or pass the "
            f"values on the command line. Known: {', '.join(sorted(d))}")
    return d[accession]


def fact(accession: str, key: str):
    v = library(accession).get(key)
    if v is None:
        raise UnknownLibrary(
            f"{accession} has no {key} in {FACTS.relative_to(ROOT)}. It is a property of "
            f"that library and this code does not have one to fall back on.")
    return v


def resolve(accession: str | None, key: str, given):
    """The value the caller gave, else the library's, else a refusal naming both."""
    if given is not None:
        return given
    if accession is None:
        raise UnknownLibrary(
            f"--{key.replace('_', '-')} was not given and no --library was named, so "
            f"there is nothing to read it from. Name the library or pass the value.")
    return fact(accession, key)
