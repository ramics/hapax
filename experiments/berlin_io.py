"""Streaming reader for the Berlin wastewater run, shared by the subtraction experiments.

Mate 2 of SRR23973193 carries ~14 nt of non-template 5' sequence (measured in
results/ww_berlin_consensus/REPORT.md: 2,455 of 2,572 reads soft-clipped, mode 14).
Everything downstream is wrong if it is not trimmed first: it inflates insert estimates,
and here it would donate a constant 14 nt of unmatched sequence to every mate 2, which is
exactly the kind of artefact a subtraction experiment must not mistake for unknown
sequence.

`skip` selects a block of the file rather than a reshuffle: an Illumina run is laid out by
flow-cell tile, so distinct blocks are distinct regions of the flow cell. This is the same
seeding convention as experiments/insert_profile.py (PR #10), kept identical so results
from the two threads are comparable.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def read_pairs(f1: Path, f2: Path, limit: int, trim_r2: int = 14, skip: int = 0):
    """Yield (name, seq1, seq2, qual1, qual2) with mate 2's non-template prefix removed."""
    with gzip.open(f1, "rt") as h1, gzip.open(f2, "rt") as h2:
        name = None
        s1 = s2 = None
        for i, (l1, l2) in enumerate(zip(h1, h2)):
            rec = i // 4
            if rec < skip:
                continue
            phase = i % 4
            if phase == 0:
                name = l1[1:].split()[0]
            elif phase == 1:
                s1, s2 = l1.strip(), l2.strip()[trim_r2:]
            elif phase == 3:
                yield name, s1, s2, l1.strip(), l2.strip()[trim_r2:]
                if rec + 1 - skip >= limit:
                    return


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
