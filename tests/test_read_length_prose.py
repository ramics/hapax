"""A read length written into prose must name the library it was measured on.

WHY THIS EXISTS. `experiments/libraries.json` exists because a trim length is a fact about
one sequencing run and had been written into twenty scripts as a default. The same is true
of a read length, and the same rule had not been applied to the WORDS: the driver wrote
`a 107 nt query against a large index finds about 40 bases at about 88% by chance` into
every run's `strains.json`, as a bare property of "a query". On SRR34217566, which is
2 x 150, that sentence appeared in the library's own output and told a reader its reads
were 107 nt. The measurement is sound and worth quoting; what it cannot do is describe
another run.

So a literal that names a read length in nucleotides has to name an accession, which is
what makes it a quotation rather than a claim about the library in hand. Docstrings are
NOT exempt: the reasoning is the research, and a reader of a docstring is as entitled to
know which library a number came from.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: A read length in prose. Bounds chosen to cover every library in the repo and the common
#: Illumina lengths, rather than any number followed by `nt`.
LENGTH = re.compile(r"\b(?:93|100|101|107|125|150|151|250|300) nt\b")
#: An accession, which is what turns the length into a quotation with a source.
ACCESSION = re.compile(r"\b[SED]RR\d+\b")

#: The package, plus the renderer. The renderer is here because it is the one file whose
#: strings a reader of a report actually sees: two sentences in its detail panels carried
#: SRR23973193's 107 nt to a 2 x 150 library's own page, which is the same defect this
#: test was written for one layer out. The rest of `experiments/` is not swept yet - 51
#: docstrings in the archived experiments quote a length without an accession - so the
#: line is drawn at what ships to a page.
FILES = tuple(sorted((ROOT / "src" / "hapax").glob("*.py"))) + \
    (ROOT / "experiments" / "make_web.py",)


def test_every_read_length_in_prose_names_its_library():
    unsourced = []
    for path in FILES:
        for node in ast.walk(ast.parse(path.read_text())):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            if LENGTH.search(node.value) and not ACCESSION.search(node.value):
                unsourced.append(f"{path.name}:{node.lineno}: {node.value[:120]!r}")
    assert not unsourced, (
        "a read length is a fact about one sequencing run, and these name one without "
        "saying which run it was measured on:\n" + "\n".join(unsourced))
