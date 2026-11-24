"""The facts that are properties of the FILE, measured from the file.

LIFTED FROM `experiments/measure_library.py`, which stays byte-identical, for the reason
`hapax.assembly_input` gives: a script's hash is tied to what it produced, so the package
carries the logic and `tests/test_lifted_glue.py` asserts the two give the same answer on the
same reads. The reasoning behind each value (why the trim is not here, why the stride is
rounded down and kept off its ceiling, why one seed has no stride) is in that script's
docstring and comments, kept with the code they explain below.
"""
from __future__ import annotations

import gzip
import json
import time
from collections import Counter
from pathlib import Path

from .provenance import record

#: Fraction of the ceiling left unused, so a file that counts short still fits five seeds.
STRIDE_MARGIN = 0.10
#: Reads whose lengths are histogrammed. The length of a library's reads is a property of
#: the run, so a prefix of the file settles it; the COUNT is not, and is counted in full.
LENGTH_SAMPLE = 100_000


def _open(path: Path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def count_and_lengths(path: Path, sample: int = LENGTH_SAMPLE) -> tuple[int, Counter]:
    """Records in a FASTQ, and the length histogram of its first `sample` records.

    Counts lines and divides, which is what every FASTQ reader here does; a record whose
    sequence is wrapped over several lines would break it, and no instrument writes one.
    The file is read once, streaming, so the cost is one pass per mate.
    """
    lengths: Counter = Counter()
    lines = seen = 0
    with _open(path) as h:
        for i, line in enumerate(h):
            lines += 1
            if i % 4 == 1 and seen < sample:
                lengths[len(line.strip())] += 1
                seen += 1
    if lines % 4:
        raise ValueError(f"{path} has {lines} lines, not a multiple of 4: truncated?")
    return lines // 4, lengths


def stride_for(pairs: int, pairs_per_seed: int, seeds: int,
               margin: float = STRIDE_MARGIN) -> int:
    """The largest round stride that fits `seeds` disjoint blocks, with margin.

    Rounded down to 2 significant figures so the number in a manifest is one a reader can
    check against the pair count in their head, which is the only defence against a stride
    that is quietly off by an order of magnitude.
    """
    if seeds < 2:
        raise ValueError("a stride is only meaningful for two or more seeds")
    ceiling = (pairs - pairs_per_seed) // (seeds - 1)
    # DISJOINT, not merely in-range: a stride smaller than the block makes the seeds
    # OVERLAP, so five seeds would re-read the same reads and the spread between them
    # would be an artefact of the stride rather than of the flow cell. The arithmetic
    # passes happily in that case, which is why it is checked rather than assumed.
    if ceiling < pairs_per_seed:
        raise ValueError(
            f"{pairs} pairs cannot hold {seeds} DISJOINT blocks of {pairs_per_seed}: the "
            f"largest stride that fits is {ceiling}, below the block itself, so the seeds "
            f"would overlap. {seeds * pairs_per_seed} pairs are the minimum. Use fewer "
            f"seeds or a smaller block.")
    usable = int(ceiling * (1 - margin))
    digits = max(0, len(str(usable)) - 2)
    unit = 10 ** digits
    stride = max(unit, (usable // unit) * unit)
    if stride < pairs_per_seed:
        # The margin and the rounding together can cross the line the ceiling cleared.
        # Take the ceiling instead and say so, rather than returning overlapping blocks.
        return ceiling
    return stride


def measure(r1: Path, r2: Path | None, pairs_per_seed: int, seeds: int) -> dict:
    t0 = time.time()
    pairs, l1 = count_and_lengths(Path(r1))
    out = {"pairs": pairs, "pairs_basis": f"records counted in {Path(r1).name}",
           "read_length": [max(l1, key=l1.get)]}
    if r2 is not None:
        n2, l2 = count_and_lengths(Path(r2))
        if n2 != pairs:
            raise ValueError(
                f"mate 1 has {pairs} records and mate 2 has {n2}. They are not pairs, and "
                f"every stage here reads them as pairs.")
        out["read_length"].append(max(l2, key=l2.get))
        out["layout"] = "paired"
    else:
        out["layout"] = "single"
    out["read_length_basis"] = (
        f"modal sequence length over the first {LENGTH_SAMPLE} records of each mate. "
        f"Read length is the one recognition wall this project has measured moving "
        f"(SPEC 4.2), so it is recorded rather than assumed comparable to another "
        f"library's.")
    out["read_length_spread"] = {
        "r1": sorted(l1.items(), key=lambda kv: -kv[1])[:5],
        **({"r2": sorted(l2.items(), key=lambda kv: -kv[1])[:5]} if r2 is not None else {}),
    }
    # ONE SEED HAS NO STRIDE, AND THAT IS NOT A FAILURE TO MEASURE ONE. The stride exists
    # to keep several seeds reading disjoint blocks of the file, so with one seed there is
    # nothing for it to separate and `stride_for` raises rather than inventing a number.
    # Refusing the whole measurement on that basis would mean a one-seed instrument - a
    # fixture, a pilot - could not be measured at all, so the absence is recorded with its
    # reason instead. A headline result still needs five seeds (CLAUDE.md); this records
    # what this file can hold, not what a result may claim.
    out["seeds"] = seeds
    if seeds < 2:
        out["block"] = None
        out["block_absent_because"] = (
            f"one seed was asked for, and a stride only separates two or more. The single "
            f"block is the first {pairs_per_seed} pairs of {pairs}.")
    else:
        out["block"] = stride_for(pairs, pairs_per_seed, seeds)
        out["block_basis"] = (
            f"derived from the {pairs} pairs counted in this file: the largest round "
            f"stride that fits {seeds} disjoint contiguous blocks of {pairs_per_seed} "
            f"pairs with a {STRIDE_MARGIN:.0%} margin, so a file that counts shorter than "
            f"its metadata still reads five blocks rather than leaving the last seeds "
            f"empty.")
    out["measured_not_stated"] = (["pairs", "read_length", "layout"]
                                  + (["block"] if out["block"] is not None else []))
    out["not_measurable_from_the_file"] = {
        "trim_r2": "needs the reads mapped onto known sequence; see "
                   "experiments/nontemplate_prefix_map.py and trim_from_prefix.py",
        "sample": "the matrix, place, collection date and molecule are properties of the "
                  "water and of the deposit, not of the FASTQ",
        "crossplant_panel": "membership of a committed panel, not a property of the reads",
    }
    out["wall_seconds"] = round(time.time() - t0, 1)
    return out


def measure_to(r1: Path, r2: Path | None, pairs_per_seed: int, seeds: int,
               out: Path, accession: str | None = None) -> dict:
    """`measure`, plus what the script's main adds, written to `out`."""
    m = measure(r1, r2, pairs_per_seed, seeds)
    m["accession"] = accession
    m["r1"], m["r2"] = str(r1), (str(r2) if r2 else None)
    m["provenance"] = record(__file__, "hapax.measure.measure_to")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(m, indent=1) + "\n")
    return m
