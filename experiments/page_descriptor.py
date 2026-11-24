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
    ("country", "{}", ", "),
    ("collected", "collected on {}", ", "),
    ("molecule", "sequenced as {}", ", "),
)

#: Keys of a `sample` block that are not clauses. `basis` is how the facts were
#: established and is read by a person, not spoken by the sentence.
NOT_SPOKEN = ("basis",)


def sample_facts(accession: str, facts: dict | None = None) -> dict:
    """The library's `sample` block, or {} where it has none. Never another library's.

    `facts` lets the CALLER supply the block instead, for water that is not in
    `experiments/libraries.json` at all. That is the case this tool is for: a production
    run on an input nobody has deposited from has no registry entry, and until this
    argument existed such a run could not have a page, because the registry raises for a
    library it does not know - correctly, since the alternative is another library's
    facts. Supplied facts are an override with a provenance of their own (the run's
    config, recorded in the run directory), not a fallback: when they are given, nothing
    is read from the registry, and when they are not, the registry's refusal stands.
    """
    if facts is not None:
        return facts
    return library_facts.library(accession).get("sample") or {}


def check(accession: str, facts: dict | None = None) -> None:
    """Refuse a committed fact the grammar has no clause for.

    The sibling thread caught this on the library that the whole genericness test turns
    on: `country` was a committed fact with a basis in four sample blocks and there was no
    clause for it, so the sentence dropped it silently - losing "South Korea" on the one
    library that is not in Europe. A fact somebody wrote down, with its basis, that the
    renderer ignores is exactly the asymmetry this module exists to close, one layer in.
    So an unknown key raises rather than going quietly nowhere: either it gets a clause or
    it does not belong in the block.
    """
    spoken = {k for k, _, _ in CLAUSES} | set(NOT_SPOKEN)
    extra = sorted(set(sample_facts(accession, facts)) - spoken)
    if extra:
        raise library_facts.UnknownLibrary(
            f"{accession}: sample block carries {', '.join(extra)}, which "
            f"page_descriptor.CLAUSES has no clause for, so the sentence would drop it "
            f"without saying so. Add a clause or remove the fact.")


def lede(accession: str, facts: dict | None = None) -> str:
    check(accession, facts)
    s = sample_facts(accession, facts)
    parts = [(t.format(s[k]), j) for k, t, j in CLAUSES if s.get(k)]
    if not parts:
        where = ("the run's configuration" if facts is not None
                 else "experiments/libraries.json")
        return (f"{accession} has no sample facts in {where}, so this page says nothing "
                f"about what the water was. Add them, with how each fact was "
                f"established, and the sentence appears here.")
    out = parts[0][0][0].upper() + parts[0][0][1:]
    for text, join in parts[1:]:
        out += join + text
    return out + "."


def title(accession: str, facts: dict | None = None) -> str:
    """Named by its place where the library records one, by its accession otherwise.

    Imo's reading preference is names, not accessions. It is not applied where the name
    would have to be invented: a library whose deposit gives no locality is titled by the
    only identifier that is certainly its own.
    """
    s = sample_facts(accession, facts)
    return f"{s.get('place') or accession} run viewer"


def derive(accession: str, facts: dict | None = None) -> dict:
    return {"title": title(accession, facts), "lede": lede(accession, facts),
            "facts": sample_facts(accession, facts)}


if __name__ == "__main__":  # a one-line check for a library, before a page exists
    for a in sys.argv[1:]:
        d = derive(a)
        print(f"{a}\n  title: {d['title']}\n  lede:  {d['lede']}")
