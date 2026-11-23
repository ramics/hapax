#!/usr/bin/env python3
"""The run viewer's per-sample description, derived rather than written.

Why. Every page needs a title and one line saying what the water was. Those were typed
into `results/web/pages/<sample>.json` by whoever built the page, which meant a sample
nobody had touched could not have a page until somebody wrote a sentence for it - and the
sentence was unattributable English sitting beside numbers that all carry provenance.
That is the same asymmetry the vocabulary file was built to close: the convention that a
figure is generated from committed files was being applied to the numbers and not to the
claims beside them, and the claims are usually larger.

So the description is assembled here from facts that are properties of the library, in
`experiments/libraries.json` under `sample`, each with how it was established. This module
owns the grammar and holds no facts of its own: it joins the clauses whose facts exist and
leaves out the clauses whose facts do not. **A fact that is absent is left out; it is never
guessed**, and a library with no sample block gets a page that says that plainly rather
than a page with a confident-sounding sentence nobody can attribute.

A page descriptor may still state `title` or `lede` itself. That is a deliberate override,
not a default, so it is recorded as one: `derive()` returns what it would have said, and
`make_web.py` keeps the override and marks it. The rule the repo keeps arriving at applies
here too - the fallback is the bug - so there is no silent substitution in either
direction.

The test of the whole arrangement is operational, and it is SPEC 13.4's: a library nobody
has touched produces its page with nobody writing a sentence for it.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts  # noqa: E402

# The clause order of the sentence, and the only grammar in this file. Each entry is the
# fact's key in the library's `sample` block, how that fact is spoken, and whether it joins
# the clause before it with a space rather than a comma: "influent from Berlin" is one
# noun phrase, and the two clauses after it are a list.
CLAUSES = (
    ("matrix", "{}", " "),
    ("place", "from {}", " "),
    ("collected", "collected on {}", ", "),
    ("molecule", "sequenced as {}", ", "),
)


def sample_facts(accession: str) -> dict:
    """The library's `sample` block, or {} where it has none. Never another library's."""
    return library_facts.library(accession).get("sample") or {}


def lede(accession: str) -> str:
    s = sample_facts(accession)
    parts = [(t.format(s[k]), j) for k, t, j in CLAUSES if s.get(k)]
    if not parts:
        return (f"{accession} has no sample block in experiments/libraries.json, so this "
                f"page says nothing about what the water was. Add one, with how each "
                f"fact was established, and the sentence appears here.")
    out = parts[0][0][0].upper() + parts[0][0][1:]
    for text, join in parts[1:]:
        out += join + text
    return out + "."


def title(accession: str) -> str:
    """Named by its place where the library records one, by its accession otherwise.

    Imo's reading preference is names, not accessions. It is not applied where the name
    would have to be invented: a library whose deposit gives no locality is titled by the
    only identifier that is certainly its own.
    """
    s = sample_facts(accession)
    return f"{s.get('place') or accession} run viewer"


def derive(accession: str) -> dict:
    return {"title": title(accession), "lede": lede(accession),
            "facts": sample_facts(accession)}


if __name__ == "__main__":  # a one-line check for a library, before a page exists
    for a in sys.argv[1:]:
        d = derive(a)
        print(f"{a}\n  title: {d['title']}\n  lede:  {d['lede']}")
