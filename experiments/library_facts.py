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
import re
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


# ------------------------------------------------------------------ the argparse shape
#
# Resolving the trim from a named library leaves one hole open, and it is the hole that
# matters. `--library` defaults to the discriminating dataset so that every committed
# command reproduces unchanged; a caller who points `--r1` at a different FASTQ and says
# nothing would therefore still get Berlin's trim. So the read paths are checked against
# the library: a path that does not carry the named accession means the two disagree, and
# a per-library value that was not stated explicitly is refused rather than guessed. The
# accession is the test because staged copies live under many directories - a container
# scratch dir, an EC2 mount, data/ww_berlin and data/berlin all appear in this repo's
# manifests - but every one of them keeps the accession in the file name.

DEFAULT_LIBRARY = "SRR23973193"

LIBRARY_HELP = ("accession in experiments/libraries.json, which supplies the values that "
                "are properties of that sequencing run. An unknown library raises rather "
                "than inheriting another library's")


def add_argument(ap, default: str | None = DEFAULT_LIBRARY) -> None:
    """Add --library. A script that cannot name one passes default=None and must be told."""
    ap.add_argument("--library", default=default, help=LIBRARY_HELP)


# Not \b: an underscore is a word character, so `ww_SRR23421544` would not match and
# the one case this check exists for would sail through.
RUN_IN_PATH = re.compile(r"(?<![A-Za-z0-9])[SED]RR\d{5,}")


def paths_name_library(accession: str, paths) -> bool:
    """True unless a given path names a run that is not this one.

    A path naming no run at all cannot disagree: `data/berlin` and a work directory are
    the library's reads under a name that does not carry the accession, and refusing those
    would refuse every committed command. A path naming a DIFFERENT run is the case this
    exists for, and it is the one that silently applied one library's prefix to another's
    reads.
    """
    for p in paths:
        if p is None:
            continue
        named = RUN_IN_PATH.findall(str(p))
        if named and accession not in named:
            return False
    return True


def resolve_args(args, keys, path_attrs=("r1", "r2"), library_attr="library") -> None:
    """Fill each of `keys` on `args` from the named library, in place.

    `keys` is a sequence of names that are both the attribute and the fact, or a mapping
    of attribute to fact where a script's flag is spelled differently or where it wants a
    different operating point of the same quantity (`{"trim_r2": "trim_r2_for_mapback"}`).

    A value the caller gave is kept and never checked against the library: stating it is
    how to override. A value that was not given is read from the library only if the read
    paths agree with the library named; otherwise this raises, because that is exactly the
    case where one library's biology is about to be applied to another library's reads.
    """
    acc = getattr(args, library_attr, None)
    paths = [getattr(args, a, None) for a in path_attrs]
    agree = acc is not None and paths_name_library(acc, paths)
    pairs = keys.items() if isinstance(keys, dict) else ((k, k) for k in keys)
    for attr, key in pairs:
        if getattr(args, attr, None) is not None:
            continue
        if not agree:
            named = ", ".join(str(p) for p in paths if p is not None)
            raise UnknownLibrary(
                f"--{attr.replace('_', '-')} was not given, --{library_attr} is {acc}, "
                f"and the read path does not carry that accession ({named}). Name the "
                f"library these reads come from, or pass the value: {key} is a property "
                f"of a sequencing run and this code has no right one to assume.")
        setattr(args, attr, fact(acc, key))
