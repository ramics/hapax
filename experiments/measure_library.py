#!/usr/bin/env python3
"""The facts that are properties of the FILE, measured from the file.

WHY THIS EXISTS, AND WHAT IT CORRECTS. `experiments/libraries.json` was built for
experiments, where pinning a value beside how it was established is how a result
reproduces. Production inherited the habit, and that is wrong in the one case production
is for: a sample nobody has touched. Running one therefore began by a person committing an
entry - pair count, read lengths, seed stride, read paths - before any stage would start.

Imo, 2026-11-23, on being shown that list: *"why are we not discovering the facts from the
sample itself?"* Three of those four are in the file. The pair count is a count. The read
lengths are lengths. The stride is arithmetic on the count. None of them needs a person,
and asking for them turns "new sample" into "new sample plus a commit".

So this measures them, writes them beside the run with provenance, and prints them for a
caller to state explicitly. The committed registry stays what it is good for: an override
with a basis, for the case where a measurement exists and you want to disagree with it,
and for the facts that genuinely are not in the file.

WHAT IS NOT IN THE FILE, and so is not here. Named rather than silently skipped, because
the point of the exercise is to know what is left:

- `trim_r2`, the non-template prefix. Needs the reads mapped onto something already known
  and the clip read off the front, which needs the sample to carry enough known sequence:
  `experiments/nontemplate_prefix_map.py`, and it REFUSED on SRR23702504 (0.13% of pairs
  anchored). Not a count, not derivable from read length, and the one fact on this library
  that had to be stated by a person.
- the matrix, place, date and molecule: properties of the water and of the deposit, not of
  the FASTQ. They go in the `sample` block and the page says so when they are absent.
- the cross-plant panel, which is a question about membership of a committed panel.

HOW THE STRIDE IS DERIVED, and the bug it exists to stop. Five seeds are five disjoint
contiguous blocks of one flow cell, so the stride has to satisfy
`(seeds - 1) * stride + pairs_per_seed <= pairs`. Berlin's 1.5 M written into eight scripts
left seeds 4 and 5 of a 4.29 M-read library reading NOTHING, and the directory still looked
like a five-seed result. Measured from the count there is no such case: the ceiling is
computed, a margin is taken off it, and the value is rounded down so it is a number a human
can check. It is deliberately NOT the ceiling itself - a file that counts shorter than its
metadata (a truncated transfer, a different filter) would then run off the end.

The spread of five such blocks is an observed range and never a dispersion estimate
(CLAUDE.md), which is a property of the design this file cannot change and does not try to.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import provenance  # noqa: E402

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
    out["block"] = stride_for(pairs, pairs_per_seed, seeds)
    out["block_basis"] = (
        f"derived from the {pairs} pairs counted in this file: the largest round stride "
        f"that fits {seeds} disjoint contiguous blocks of {pairs_per_seed} pairs with a "
        f"{STRIDE_MARGIN:.0%} margin, so a file that counts shorter than its metadata "
        f"still reads five blocks rather than leaving the last seeds empty.")
    out["measured_not_stated"] = ["pairs", "read_length", "layout", "block"]
    out["not_measurable_from_the_file"] = {
        "trim_r2": "needs the reads mapped onto known sequence; see "
                   "experiments/nontemplate_prefix_map.py and trim_from_prefix.py",
        "sample": "the matrix, place, collection date and molecule are properties of the "
                  "water and of the deposit, not of the FASTQ",
        "crossplant_panel": "membership of a committed panel, not a property of the reads",
    }
    out["wall_seconds"] = round(time.time() - t0, 1)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, default=None,
                    help="omit for a single-end library, which is recorded as single")
    ap.add_argument("--pairs-per-seed", type=int, required=True,
                    help="the block each seed reads; the stride has to fit all of them")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--accession", default=None,
                    help="recorded for the reader; nothing here is read from the registry")
    ap.add_argument("--out", type=Path, default=None, help="where to write the manifest")
    ap.add_argument("--print", dest="show", default=None,
                    choices=["pairs", "block", "read_length", "shell"],
                    help="print one value, or `shell` for eval-able assignments, and "
                         "nothing else, so a run script can read it")
    a = ap.parse_args()

    m = measure(a.r1, a.r2, a.pairs_per_seed, a.seeds)
    m["accession"] = a.accession
    m["r1"], m["r2"] = str(a.r1), (str(a.r2) if a.r2 else None)
    m["provenance"] = provenance.provenance(__file__)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(m, indent=1) + "\n")
    if a.show == "shell":
        print(f"PAIRS={m['pairs']}")
        print(f"BLOCK={m['block']}")
        print(f"READ_LENGTH={','.join(str(x) for x in m['read_length'])}")
    elif a.show == "read_length":
        print(",".join(str(x) for x in m["read_length"]))
    elif a.show:
        print(m[a.show])
    else:
        json.dump(m, sys.stdout, indent=1)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
