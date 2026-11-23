#!/usr/bin/env python3
"""Reading a run directory's row dumps, gzipped or plain, from one place.

Why this is a module and not three lines in a reader. The first production run over an
untouched library gzipped its row dumps, because one of them is 79 MB a seed; every
earlier run directory in this repo writes them plain. A reader that knows only the plain
name does not fail on the other spelling - it finds no file and reports a zero, which is
how a page over 32 M pairs drew an empty residue table, an empty scatter and two empty
search arms while saying nothing was wrong.

Four scripts read these dumps by name (make_web.py, subject_names.py,
concern_ranking.py, holdout_gene_map.py). The same shuffle bug shipped twice in two days
because the second writer did not know which condition the first was fixing, so this is
imported rather than re-typed, and like the other shared modules it raises where a copy
would fall back: a directory carrying both spellings refuses, because nothing in it says
which one the run wrote.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path


def jsonl(rundir: Path, name: str):
    """Rows of a run directory dump, gzipped or plain; a refusal where it is both.

    The production run writes `reads.jsonl.gz`; the earlier runs in this repo write
    `reads.jsonl`. The reader looked for the plain name only, so the first real
    production directory built a page with **zero per-read rows and said nothing**: the
    residue table, the scatter and both search arms drew empty, which is the page
    reporting "nothing here" for a run of 32 M pairs. A default where a refusal belongs,
    the family this repo keeps paying for. So both spellings are read, and a directory
    that carries both refuses rather than picking one, because nothing says which is
    current.
    """
    plain, gz = rundir / name, rundir / (name + ".gz")
    if plain.exists() and gz.exists():
        raise SystemExit(f"{rundir} carries both {name} and {name}.gz and nothing says "
                         f"which the run wrote; remove one")
    if gz.exists():
        with gzip.open(gz, "rt") as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]
    if plain.exists():
        with open(plain) as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]
    return None
