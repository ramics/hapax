"""The shuffled null has one home, and these tests are why it can be trusted there.

Run with `python3 -m pytest tests/test_shuffle_control.py`. Pure Python, no external tools,
about two seconds.

WHY THIS FILE EXISTS. On 2026-11-21 and again on 2026-11-22 a local copy of the
dinucleotide shuffle dead-ended on real sequence and returned its input unchanged. Both
times the null quietly became a slice of the test, and published numbers had to be
withdrawn. Both times the copy had been read and passed by eye. So this file checks
three things by running code rather than reading it:

1. The canonical module (`experiments/shuffle_control.py`) holds its own invariant on
   every sequence of a corpus built to break walks: same length, same composition, not
   the input. Real sequence from the repo plus three low-complexity families.
2. `mono_shuffle` draws exactly what each grandfathered mononucleotide copy draws for the
   same seed. That is the proof that switching a script to the shared function changes
   none of its outputs, only its script hash.
3. No new shuffle appears in `experiments/` outside the module. Every existing copy is
   listed below with its measured behaviour. The list may shrink; it may not grow.
"""

from __future__ import annotations

import ast
import collections
import random
import sys
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments"
sys.path.insert(0, str(EXP))

import shuffle_control as sc  # noqa: E402

# ----------------------------------------------------------------------- the allowlist
# Every function in experiments/*.py whose name contains "shuffle", other than the module
# itself. A new one fails `test_no_new_shuffle_copies`: import from shuffle_control
# instead. Remove an entry when its script switches to the shared module.
#
# kind: "dinuc" or "mono" for a sequence shuffle that the behavioural tests run;
# "wrapper" for a function that calls a shuffle rather than being one;
# "audit" for the measuring harness that holds copies on purpose.
GRANDFATHERED = {
    # Dinucleotide walks. Measured on this file's corpus (see test_dinuc_copies_measured):
    # all but the canonical one hand back their input on dinucleotide-repeat sequence.
    ("holdout_run.py", "dinuc_shuffle"): ("dinuc", "tree condition, returns None on "
                                          "failure and its caller counts; publishes "
                                          "results/holdout_run/"),
    ("residue_structure.py", "dinuc_shuffle"): ("dinuc", "tree condition, 50 retries "
                                                "then returns input; 2 of 7,390 reads"),
    ("residue_contig_annotation.py", "dinuc_shuffle"): ("dinuc", "tree condition, 100 "
                                                        "retries then returns input"),
    # Mononucleotide permutations: no walk, so no dead end. Proven draw-identical to
    # sc.mono_shuffle by test_mono_copies_draw_identically.
    ("residue_artefacts.py", "shuffle"): ("mono", ""),
    ("readlevel_features.py", "shuffle_seq"): ("mono", ""),
    ("residue_contig_overlap.py", "shuffle"): ("mono", ""),
    ("residue_report.py", "shuffled"): ("mono", ""),
    ("profile_reach.py", "shuffle_seq"): ("mono", ""),
    ("recurrence_queries.py", "shuffled"): ("mono", ""),
    # Arriving with PR #32 (cross-plant), written the same day as this list: two more.
    ("crossplant_contig_queries.py", "shuffled"): ("mono", "PR #32"),
    ("crossplant_sample.py", "shuffled"): ("mono", "PR #32"),
    # Calls a shuffle per read of a FASTQ; its inline permutation is a mono shuffle.
    ("public_health_profile.py", "shuffle_fastq"): ("wrapper", ""),
    ("holdout_run.py", "shuffle_arm"): ("wrapper", "counts failed and returned_input"),
    # The harness that measures the copies above holds its own references to them.
    ("shuffle_fallback_audit.py", "*"): ("audit", ""),
}

# Functions whose name matches the tripwire but whose body is an import rather than a walk.
# They are not copies and are not measured by the behavioural tests, so they are kept out of
# GRANDFATHERED, whose entries record a copy's measured behaviour. The tripwire is
# name-based on purpose, and a loader is the one false positive that shape produces.
IMPORTERS = {
    ("reach_queries.py", "load_shuffle"): ("loads holdout_run.dinuc_shuffle by file path, "
                                           "because experiments/ is not a package; "
                                           "publishes results/reach/"),
}

CANONICAL = "shuffle_control.py"


def _top_level_functions(path: Path):
    with warnings.catch_warnings():               # other scripts' escape sequences
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(path.read_text())
    return [n for n in tree.body if isinstance(n, ast.FunctionDef)]


def _extract(fname: str, func: str):
    """Compile one top-level function out of a script without importing the script.

    Importing would pull in each script's heavy dependencies (mappy, pyrodigal, numpy)
    and run module-level code; the shuffles need only random and collections. Executing
    the function's own source is what makes this a test of the committed code rather
    than of a re-typed copy of it.
    """
    path = EXP / fname
    if not path.exists():
        return None
    for node in _top_level_functions(path):
        if node.name == func:
            ns = {"random": random, "collections": collections,
                  "defaultdict": collections.defaultdict, "Counter": collections.Counter}
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), ns)
            return ns[func]
    return None


# --------------------------------------------------------------------------- the corpus
def _fasta(path: Path) -> dict[str, str]:
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def _corpus() -> dict[str, list[str]]:
    """Real sequence in this library's two mate lengths and whole, plus walk-breakers.

    Real: the 78 committed holdout candidate genomes (`results/holdout_candidates/`),
    whole and cut into 107 and 93 nt windows. Synthetic, because "clean on real reads" was
    the claim that failed twice: poly-A flanks (homopolymer vertices), AT-only sequence
    (a two-letter graph), and a dinucleotide repeat with a short random tail, which is
    where every retry-limited walk in this repo gives up.
    """
    g = _fasta(ROOT / "results/holdout_candidates/candidates.fa")
    real = []
    for s in g.values():
        real.append(s)
        for w in (107, 93):
            real.extend(s[i:i + w] for i in range(0, len(s) - w, w * 7))
    r = random.Random(0)
    poly, at, rep = [], [], []
    for n in (20, 93, 107, 300):
        for _ in range(100):
            poly.append("A" * r.randint(5, 30) + "".join(r.choice("ACGT") for _ in range(n))
                        + "A" * r.randint(5, 30))
            at.append("".join(r.choice("AT") for _ in range(n)))
            rep.append("AC" * (n // 2) + "".join(r.choice("ACGTN") for _ in range(10)))
    return {"real": real, "polyA_flanked": poly, "AT_only": at, "dinuc_repeat": rep}


CORPUS = _corpus()


# ------------------------------------------------------------------- canonical module
@pytest.mark.parametrize("family", sorted(CORPUS))
def test_canonical_dinuc_holds_invariant_on_every_sequence(family):
    """Every output is a valid dinucleotide shuffle; nothing raises on this corpus; the
    only sequences returned unchanged are ones with no other arrangement."""
    rng = random.Random(1)
    identical = 0
    for s in CORPUS[family]:
        out = sc.dinuc_shuffle(s, rng)
        if out == s.upper():
            identical += 1
            continue
        sc.check_shuffle(s, out, "dinuc")
    # A handful of short sequences genuinely have one Eulerian path or a few; a rate like
    # the broken copies' (tens of percent on dinuc_repeat) would fail here.
    assert identical <= 0.01 * len(CORPUS[family]), (family, identical)


def test_canonical_arm_census_and_validation():
    seqs = {f"r{i}": s for i, s in enumerate(CORPUS["real"][:300])}
    seqs["homopolymer"] = "A" * 50
    res = sc.shuffle_arm(seqs, random.Random(3), kind="dinuc")
    c = res["census"]
    assert c["offered"] == len(seqs)
    assert c["shuffled"] + c["dropped_identical"] == c["offered"]
    assert "homopolymer" in res["dropped_identical"]
    for name, out in res["shuffled"].items():
        sc.check_shuffle(seqs[name], out, "dinuc")


def test_arm_is_reproducible_for_a_seed():
    seqs = {f"r{i}": s for i, s in enumerate(CORPUS["real"][:50])}
    a = sc.shuffle_arm(seqs, random.Random(7))
    b = sc.shuffle_arm(seqs, random.Random(7))
    assert a == b


def test_check_shuffle_rejects_the_bug():
    """The defect that cost two days: an 'output' that is the input."""
    s = CORPUS["real"][0]
    with pytest.raises(sc.ShuffleFailed):
        sc.check_shuffle(s, s, "dinuc")
    with pytest.raises(sc.ShuffleFailed):
        sc.check_shuffle(s, s[::-1], "dinuc")          # same bases, wrong doublets
    with pytest.raises(sc.ShuffleFailed):
        sc.check_shuffle(s, s[:-1], "mono")
    with pytest.raises(ValueError):
        sc.check_shuffle(s, s, "trinuc")


def test_short_sequence_raises_not_returns():
    with pytest.raises(sc.ShuffleFailed):
        sc.dinuc_shuffle("ACG", random.Random(1))


# ------------------------------------------------------- grandfathered mono copies
MONO = [k for k, v in GRANDFATHERED.items() if v[0] == "mono"]


@pytest.mark.parametrize("fname,func", MONO)
def test_mono_copies_draw_identically(fname, func):
    """Same seed, same draws, same outputs: switching to sc.mono_shuffle is a no-op on
    every number the script writes. Interleaved calls, as the scripts make them."""
    fn = _extract(fname, func)
    if fn is None:
        pytest.skip(f"{fname}:{func} is gone, which is the point; drop it from the list")
    seqs = CORPUS["real"][:200] + CORPUS["polyA_flanked"][:50]
    ra, rb = random.Random(11), random.Random(11)
    assert [fn(s, ra) for s in seqs] == [sc.mono_shuffle(s, rb) for s in seqs]
    assert ra.random() == rb.random()                  # generator left in the same state


# ------------------------------------------------------ grandfathered dinuc copies
def _measure(fn, seqs, seed=1):
    rng = random.Random(seed)
    returned_input = none = 0
    for s in seqs:
        out = fn(s, rng)
        if out is None:
            none += 1
        elif out.upper() == s.upper():
            returned_input += 1
        else:
            sc.check_shuffle(s, out, "dinuc")          # when they shuffle, they shuffle right
    return returned_input, none


def test_dinuc_copies_do_not_draw_identically():
    """Unlike the mono copies, a correct dinucleotide walk is NOT interchangeable with the
    canonical one, and this pins that so nobody switches one as a tidy-up.

    `holdout_run.dinuc_shuffle` carries the tree condition and is correct: over 400 random
    sequences it preserves length, dinucleotide composition and both endpoints on every one,
    and never fails. So does `sc.dinuc_shuffle`. They nevertheless agree on **0 of 400**,
    because they consume the generator differently - Altschul-Erikson admits many valid
    walks and the two pick different ones.

    The consequence is the one that matters for the record: switching `reach_queries.py` or
    `holdout_run.py` to the canonical module would redraw every shuffled sequence and so
    move the null in `results/reach/` and `results/holdout_run/`. Both nulls read 0.00%, so
    the conclusion is robust to the redraw, but the drawn sequences and any count are not.
    **That is a re-run, not a switch** - the opposite of the mono case above, where
    `test_mono_copies_draw_identically` proves the switch changes nothing but the hash.
    """
    fn = _extract("holdout_run.py", "dinuc_shuffle")
    if fn is None:
        pytest.skip("holdout_run.py:dinuc_shuffle is gone; drop this with its list entry")
    seqs = CORPUS["real"][:400]
    agree = valid = 0
    for seq in seqs:
        a = sc.dinuc_shuffle(seq, random.Random(7))
        b = fn(seq, random.Random(7))
        assert b is not None, "holdout_run's walk failed on real sequence; see the list"
        sc.check_shuffle(seq, a, "dinuc")
        sc.check_shuffle(seq, b, "dinuc")      # both correct, which is the whole point
        valid += 1
        agree += (a == b)
    assert valid == len(seqs)
    assert agree == 0, (f"{agree} of {valid} draws now agree; if a walk was changed to make "
                        "them agree, the committed nulls that used the old one need re-running")


@pytest.mark.parametrize("fname", ["holdout_run.py", "residue_contig_annotation.py",
                                   "residue_structure.py"])
def test_tree_condition_copies_are_clean_on_real_sequence(fname):
    """These three carry the Altschul-Erikson tree condition, and on real sequence of
    this repo none hands back its input. That is what keeps their committed results
    standing. Their behaviour on low-complexity sequence is reported, not asserted, by
    test_dinuc_copies_measured: there they all give up where the canonical walk does not.
    """
    fn = _extract(fname, "dinuc_shuffle")
    if fn is None:
        pytest.skip(f"{fname} no longer carries its own copy")
    returned_input, _none = _measure(fn, CORPUS["real"])
    assert returned_input == 0


def test_divergent_known_copy_is_the_known_defect():
    """Pins the defect so nobody mistakes main's copy for a working one. When the holdout
    thread's fix lands this copy is gone and the test skips."""
    fn = _extract("divergent_known.py", "dinuc_shuffle")
    if fn is None:
        pytest.skip("divergent_known.py imports the shared shuffle now")
    returned_input, _ = _measure(fn, CORPUS["real"])
    assert returned_input > 0


def test_dinuc_copies_measured(capsys):
    """Not an assertion: the table of what each copy does on each family, printed with
    `pytest -s`, so the next audit measures instead of reading."""
    rows = []
    fns = {"shuffle_control (canonical)": sc.dinuc_shuffle}
    for (fname, func), (kind, _) in GRANDFATHERED.items():
        if kind == "dinuc":
            fn = _extract(fname, func)
            if fn is not None:
                fns[fname] = fn
    for name, fn in fns.items():
        cells = []
        for fam in sorted(CORPUS):
            ri, none = _measure(fn, CORPUS[fam])
            cells.append(f"{fam} {ri}+{none}/{len(CORPUS[fam])}")
        rows.append(f"{name:32s} " + "  ".join(cells))
    with capsys.disabled():
        print("\nreturned-input + None, per family:\n" + "\n".join(rows))


# ---------------------------------------------------------------------- the tripwire
def test_no_new_shuffle_copies():
    """A sequence shuffle lives in shuffle_control.py and is imported, never re-typed.

    The second writer of a walk does not know which non-obvious condition the first was
    fixing; that is how the same bug shipped twice in two days. Name-based on purpose:
    cheap, and a function that shuffles sequence and is not called *shuffle* is a review
    problem this cannot solve.
    """
    allowed = set(GRANDFATHERED) | set(IMPORTERS)
    found = []
    for path in sorted(EXP.glob("*.py")):
        if path.name == CANONICAL or (path.name, "*") in allowed:
            continue
        for node in _top_level_functions(path):
            if "shuffle" in node.name.lower() and (path.name, node.name) not in allowed:
                found.append(f"{path.name}:{node.lineno} def {node.name}")
    assert not found, ("new shuffle implementation(s); import from experiments/"
                       "shuffle_control.py instead:\n  " + "\n  ".join(found))
