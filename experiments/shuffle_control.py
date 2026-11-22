"""The dinucleotide-preserving shuffle, in one place, because a silent fallback cost a day.

WHY THIS FILE EXISTS. Two scripts in this repo carried their own copy of a naive
dinucleotide shuffle: build the successor lists, shuffle each, then walk. That walk can
dead-end, because a random ordering of the successor lists is not generally an Eulerian
circuit, and both copies handled the dead end by retrying twenty times and then
**returning the input unchanged**. On SRR23973193 that happened for 324 of 11,361 contigs,
2.85% of them, and every single one of the 304 "shuffled" contigs that produced a
significant protein hit was one of those 324. The null was not a null. It was 2.85% of the
real library wearing the label, and it was used to publish a 33-fold existence ratio, a
claim that no property of a hit separates real from shuffled, and two kills of other
people's fields. None of those survive.

THE CORRECT ALGORITHM is Altschul and Erikson's: the walk succeeds exactly when the
last-used outgoing edge of every vertex except the final one forms a tree rooted at the
final vertex. So choose those last edges at random, test the tree property, and retry only
that choice; then shuffle the remaining edges freely. The walk after that cannot dead-end.

AND THERE IS NO FALLBACK. A sequence this cannot shuffle raises. A null that quietly
contains real sequence is worse than no null, because it reads as a passed control.
"""
from __future__ import annotations

import collections
import random


class ShuffleFailed(RuntimeError):
    """Raised rather than returning anything the caller might mistake for a shuffle."""


def dinuc_shuffle(seq: str, rng: random.Random, tries: int = 200) -> str:
    """Return a sequence with the same length and the same dinucleotide counts.

    Same length and same dinucleotide composition means the null differs from the real
    sequence in gene order and codon structure only, which is what makes it a control for
    homology rather than for base composition.
    """
    s = seq.upper()
    if len(s) < 4:
        raise ShuffleFailed(f"sequence of length {len(s)} is too short to shuffle")
    first, last = s[0], s[-1]
    edges = collections.defaultdict(list)
    for a, b in zip(s, s[1:]):
        edges[a].append(b)
    verts = sorted(edges)

    for _ in range(tries):
        # 1. One last-edge per vertex except the final one, chosen at random.
        last_edge = {}
        ok = True
        for v in verts:
            if v == last:
                continue
            if not edges[v]:
                ok = False
                break
            # Chosen among DISTINCT successors, not among edge instances. Weighting by
            # multiplicity makes a homopolymer vertex pick itself almost every time, the
            # tree test then fails on that self-cycle, and a poly-A contig exhausts every
            # retry. Measured: with multiplicity weighting a 1,343 nt contig from this
            # library failed 200 attempts; without it the whole library shuffles.
            last_edge[v] = rng.choice(sorted(set(edges[v])))
        if not ok:
            continue
        # 2. Those last edges must form a tree rooted at the final vertex: follow each
        #    vertex's last edge and it must reach `last` without cycling.
        good = True
        for v in last_edge:
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in last_edge:
                    good = False
                    break
                seen.add(cur)
                cur = last_edge[cur]
            if not good:
                break
        if not good:
            continue
        # 3. Shuffle the rest of each vertex's edges and put its last edge at the end.
        order = {}
        for v in verts:
            rest = list(edges[v])
            if v in last_edge:
                rest.remove(last_edge[v])
            rng.shuffle(rest)
            if v in last_edge:
                rest.append(last_edge[v])
            order[v] = rest
        # 4. Walk. With the tree property above this cannot dead-end.
        out, cur, used = [first], first, collections.defaultdict(int)
        for _ in range(len(s) - 1):
            nxt = order[cur][used[cur]]
            used[cur] += 1
            out.append(nxt)
            cur = nxt
        res = "".join(out)
        if len(res) != len(s):
            raise ShuffleFailed("walk produced the wrong length, which cannot happen")
        return res
    raise ShuffleFailed(f"no Eulerian walk after {tries} attempts on a {len(s)} nt "
                        f"sequence; this is not a sequence to skip quietly")


def shuffle_all(seqs: dict, rng: random.Random, allow_identical: bool = False) -> dict:
    """Shuffle every sequence, and refuse to hand back one that did not change.

    A shuffle that returns its input is the failure mode this module exists for, so it is
    checked here too rather than trusted. Homopolymers and very low-complexity sequence
    can be genuinely invariant under a dinucleotide shuffle; those are reported by count
    and excluded, never silently passed through as a null.
    """
    out, unchanged = {}, []
    for name, s in seqs.items():
        r = dinuc_shuffle(s, rng)
        if r == s.upper() and not allow_identical:
            unchanged.append(name)
            continue
        out[name] = r
    return {"shuffled": out, "unchanged_and_dropped": unchanged}


# ---------------------------------------------------------------------------------------
# Everything below this line was added by the consolidation pass (2026-11-22). The two
# functions above are byte-identical to the copy on the holdout thread's branch, so that
# whichever lands second merges as a superset rather than a rewrite.
# ---------------------------------------------------------------------------------------

def mono_shuffle(seq: str, rng: random.Random) -> str:
    """Return a permutation of the bases: same length, same base counts, nothing else kept.

    This is the null for "is this field just base composition". It destroys dinucleotide
    structure, so it inflates stop-codon density relative to real sequence and makes the
    null EASIER to beat; a kill of the form "indistinguishable from a mononucleotide
    shuffle" survives that, a claim of separation from one is weaker than it looks.

    It consumes the generator exactly as every local mononucleotide copy in experiments/ does
    (`list`, `rng.shuffle`, `join`), so a script that switches to it draws the same
    sequence for the same seed. `tests/test_shuffle_control.py` checks that equivalence
    against every grandfathered copy, which is what makes switching a no-op on outputs.

    It cannot dead-end, so it cannot fail the way the dinucleotide walk did. It can still
    hand back its input (a homopolymer, or by chance on a very short read), and that is
    the caller's to count: use `shuffle_arm`, which drops and counts those.
    """
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def _dinuc_counts(s: str) -> collections.Counter:
    return collections.Counter(zip(s, s[1:]))


def check_shuffle(original: str, shuffled: str, kind: str) -> None:
    """Raise unless `shuffled` is a valid `kind` shuffle of `original` and not a copy of it.

    The invariant a null is published on, checked on the sequence rather than trusted from
    the code: same length; same base counts (mono) or same dinucleotide counts and the same
    first and last base (dinuc, which an Eulerian walk from first to last preserves); and
    not identical to the input. A broken walk that returned its input passed the first two
    and failed only the third, which is why the third is not optional here.
    """
    o, s = original.upper(), shuffled.upper()
    if len(o) != len(s):
        raise ShuffleFailed(f"{kind} shuffle changed length {len(o)} -> {len(s)}")
    if kind == "mono":
        if collections.Counter(o) != collections.Counter(s):
            raise ShuffleFailed("mono shuffle changed base composition")
    elif kind == "dinuc":
        if _dinuc_counts(o) != _dinuc_counts(s):
            raise ShuffleFailed("dinuc shuffle changed dinucleotide composition")
        if o[:1] != s[:1] or o[-1:] != s[-1:]:
            raise ShuffleFailed("dinuc shuffle moved the first or last base")
    else:
        raise ValueError(f"unknown shuffle kind {kind!r}; expected 'mono' or 'dinuc'")
    if o == s:
        raise ShuffleFailed(f"{kind} shuffle returned its input unchanged")


def shuffle_arm(seqs, rng: random.Random, kind: str = "dinuc") -> dict:
    """Shuffle a whole control arm, validate every output, and return it with its census.

    `seqs` is a dict or an iterable of (name, sequence) pairs; order is kept, so the
    generator is consumed in the caller's order and a fixed seed gives a fixed arm.

    The census is what a manifest must carry: `offered`, `shuffled`, and
    `dropped_identical` (sequences genuinely invariant under the shuffle, such as a
    homopolymer, or a short read whose permutation happened to come back the same). On
    Berlin, 2.85% of contigs falling through a broken walk supplied 304 of 304 of the
    null's significant hits, so a drop rate is not a perturbation to ignore: report it,
    and read the null's hits against it.

    A sequence the dinucleotide walk cannot shuffle raises `ShuffleFailed` out of this
    function. That is deliberate. Skipping it would remove exactly the repetitive inputs
    that throw chance alignments, which biases the null in a direction nobody measured.
    A caller that wants to drop such sequences must catch the exception itself and count
    them, in its own code, where a reader can see it.
    """
    if kind not in ("mono", "dinuc"):
        raise ValueError(f"unknown shuffle kind {kind!r}; expected 'mono' or 'dinuc'")
    items = seqs.items() if isinstance(seqs, dict) else seqs
    fn = mono_shuffle if kind == "mono" else dinuc_shuffle
    out, dropped, offered = {}, [], 0
    for name, s in items:
        offered += 1
        r = fn(s, rng)
        if r.upper() == s.upper():
            dropped.append(name)
            continue
        check_shuffle(s, r, kind)
        out[name] = r
    return {"shuffled": out,
            "census": {"kind": kind, "offered": offered, "shuffled": len(out),
                       "dropped_identical": len(dropped)},
            "dropped_identical": dropped}
