"""A file-valued flag must be passed a file whose NAME the script declares.

WHY THIS EXISTS. `divergent_contig_workup.py --strains` documents
`known_strain_divergence_s*.json` and reads it under a `taxa` key. The package passed it
`strains.json`, which is the `strains` stage's own output and a read-level table of
subjects with no such key. Nothing raised: selection simply found no contig in either
arm, the ORF file came out empty, and DIAMOND exited 1 on a blank query file after eleven
stages of a paid run. The line in the log read `0 divergent contigs, 0 near-reference
controls`, which is a statement about the water, and the cause was a filename two files
in one run directory both answer to.

A type checker cannot see this: both sides are `Path`. What is checkable is the name, and
every one of these scripts already states it in `--help` because a person has to type it.
So the test reads what the package passes, reads what the script declares, and compares.
Static plus `--help`, so it needs no instance, no network and no database.

It skips a flag whose help names no file, because a help that states nothing cannot be
contradicted - the test's own coverage is asserted below so that an unparseable call list
fails rather than silently checking nothing.
"""
from __future__ import annotations

import ast
import fnmatch
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: A filename in prose or in code, including the `{seed}` and `s*` spellings both sides use.
NAME = re.compile(r"[A-Za-z0-9_.*{}\-]+\.(?:json|jsonl|tsv|fa|faa|fastq|fq)(?:\.gz)?")

CALLERS = tuple(sorted((ROOT / "src" / "hapax").glob("*.py")))


def _basename(expr: str) -> str | None:
    """The last filename-shaped literal in an expression like `out / f"x_s{seed}.json"`."""
    for lit in reversed([a or b for a, b in re.findall(r"'([^']+)'|\"([^\"]+)\"", expr)]):
        found = NAME.findall(lit)
        if found:
            return found[-1]
    return None


def _help(script: str) -> str:
    done = subprocess.run([sys.executable, "-I", str(ROOT / "experiments" / script),
                           "--help"], capture_output=True, text=True, cwd="/tmp")
    return done.stdout or done.stderr


def _flag_block(text: str, flag: str) -> str:
    """The help text for one flag: from the flag to the next option at any indent."""
    found = re.search(rf"^(\s+){re.escape(flag)}\b(.*?)(?=^\s+-|\Z)", text, re.M | re.S)
    return found.group(0) if found else ""


def _pairs() -> list[tuple[str, str, str, str]]:
    """(caller, script, flag, basename) for every file-valued flag the package passes."""
    out = []
    for caller in CALLERS:
        tree = ast.parse(caller.read_text())
        env: dict[str, str] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and \
                    isinstance(node.targets[0], ast.Name):
                env.setdefault(node.targets[0].id, ast.unparse(node.value))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            scripts = [a.value for a in node.args
                       if isinstance(a, ast.Constant) and isinstance(a.value, str)
                       and a.value.endswith(".py")]
            if not scripts:
                continue
            for arglist in [a for a in node.args if isinstance(a, ast.List)]:
                els = arglist.elts
                for i, el in enumerate(els[:-1]):
                    if not (isinstance(el, ast.Constant) and isinstance(el.value, str)
                            and el.value.startswith("--")):
                        continue
                    value = els[i + 1]
                    expr = env.get(value.id, value.id) if isinstance(value, ast.Name) \
                        else ast.unparse(value)
                    name = _basename(expr)
                    if name:
                        out.append((caller.name, scripts[0], el.value, name))
    return out


def test_the_call_lists_are_still_readable():
    """If the glue is rewritten into a shape this cannot parse, say so rather than pass."""
    pairs = _pairs()
    assert len(pairs) >= 5, (
        f"found only {pairs}; the package passes file-valued flags to its stage scripts "
        f"and this test can no longer see them, so it is checking nothing")
    assert any(flag == "--strains" for _, _, flag, _ in pairs), \
        "the call this test was written for is gone; check it still covers the fault"


def test_every_file_valued_flag_gets_the_file_its_script_declares():
    wrong = []
    for caller, script, flag, passed in _pairs():
        declared = NAME.findall(_flag_block(_help(script), flag))
        if not declared:
            continue
        # `{seed}` on the package's side and `s*` on the script's side are the same
        # wildcard written in the two places' own spelling.
        got = passed.replace("{seed}", "*")
        if not any(fnmatch.fnmatch(got, d) or fnmatch.fnmatch(d, got) for d in declared):
            wrong.append(f"{caller} passes {script} {flag} {passed}, which is not one of "
                         f"{declared} as its own --help declares")
    assert not wrong, "\n".join(wrong)


#: The four facts that belong to a sequencing run and never to this code: the two read
#: paths, the mate-2 trim, and the accession the first three are checked against.
LIBRARY_FACTS = ("--r1", "--r2", "--trim-r2", "--library")


def _passed_flags() -> dict[str, set[str]]:
    """Every flag the package passes to each experiments script, across all its call sites.

    Read from the unparsed argument list rather than from its constant elements, because a
    flag can be passed conditionally - a trim that was not stated is left out so the script
    resolves it - and a starred conditional has no constant element to find.
    """
    passed: dict[str, set[str]] = {}
    for caller in CALLERS:
        for node in ast.walk(ast.parse(caller.read_text())):
            if not isinstance(node, ast.Call):
                continue
            scripts = [a.value for a in node.args
                       if isinstance(a, ast.Constant) and isinstance(a.value, str)
                       and a.value.endswith(".py")]
            lists = [a for a in node.args if isinstance(a, ast.List)]
            if not scripts or not lists:
                continue
            text = " ".join(ast.unparse(a) for a in lists)
            found = {f for f in LIBRARY_FACTS if f"'{f}'" in text or f'"{f}"' in text}
            passed.setdefault(scripts[0], set()).update(found)
    return passed


def test_a_script_the_package_runs_is_told_every_library_fact_it_accepts():
    """The driver guards its own stages against this; the package inherited no guard.

    `hapax.py` carries two tables for exactly this - LIBRARY_AWARE, which appends
    `--library` to every script that accepts one, and SAMPLE_SPECIFIC, which refuses to
    run a stage whose sample-specific flag was left to a built-in default. Both are
    asserted against the scripts' own source by `tests/test_hapax_stages.py`. Neither
    reaches the enrichments, which the PACKAGE runs: so the per-gene workup called and
    searched 76 genes of an unseen library and then opened
    `data/berlin/SRR23973193_1.fastq.gz`, its own default, after the expensive half.

    A flag the script does not accept is not required, and `--trim-r2` counts as passed
    where the package passes it conditionally, because the alternative is for the script
    to resolve it from the registry, which raises for a library it does not know.
    """
    missing = []
    for script, passed in sorted(_passed_flags().items()):
        source = (ROOT / "experiments" / script).read_text()
        accepts = [f for f in LIBRARY_FACTS
                   if f'"{f}"' in source
                   or (f == "--library" and "library_facts.add_argument" in source)]
        absent = [f for f in accepts if f not in passed]
        if absent:
            missing.append(f"{script} accepts {accepts} and the package passes none of "
                           f"{absent}")
    assert not missing, "\n".join(missing)
