"""A synthetic contig-mode run with planted rates, so the contig report is checked by output.

`crossplant_contig_report.py` aggregates windows into per-contig fractions, unions libraries
in a shuffled order and scores an own-library ratio against a panel p99. Those are exactly
the shapes that go wrong silently - an axis averaged the wrong way reads plausible. So the
report is run first against a run whose answers are known: the identified set is witnessed
in 90% of libraries, tier 2 in 20%, the shuffles in none, the own library in all, and the
own-library depth is inflated 9x. The report should return those numbers, the nulls should
read zero and the growth swap control should stay flat while the test fires.

Written because the real run costs hours of EC2 and lands once; a shape bug found then is a
relaunch, and a relaunch has already cost this experiment three.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

#: planted witnessing rates per panel library
PLANTED = {"ident": 0.9, "unident": 0.2}
LENGTHS = (200, 300, 500, 1000, 250, 700)
OWN = "SRR23973193"
OWN_DEPTH_MULTIPLIER = 9


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--panel", type=Path, default=Path("results/crossplant/panel.tsv"))
    ap.add_argument("--subjects", type=int, default=40)
    ap.add_argument("--contigs", type=int, default=6)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    rows = list(csv.DictReader(args.panel.open(), delimiter="\t"))
    subs = sorted({r["run"] for r in rows[: args.subjects]} | {OWN})
    names, meta = [], []
    for cset in PLANTED:
        for pre in ("contig", "shuffled-contig"):
            for seed in (1, 2, 3, 4, 5):
                for ci in range(args.contigs):
                    ln = LENGTHS[ci % len(LENGTHS)]
                    for p in range(2):
                        names.append(f"{pre}|{cset}|{ci}|L{ln}|p{p}|s{seed}")
                        meta.append((cset, pre))
    nq, ns = len(names), len(subs)
    rng = np.random.default_rng(args.seed)
    f = np.zeros((ns, nq))
    for q, (cset, pre) in enumerate(meta):
        p = 0.0 if pre.startswith("shuffled") else PLANTED[cset]
        f[:, q] = (rng.random(ns) < p).astype(float)
    own = subs.index(OWN)
    f[own, [q for q, m in enumerate(meta) if not m[1].startswith("shuffled")]] = 1.0
    bits = {"any1": np.packbits(f > 0, axis=1), "any2": np.packbits(f > 0, axis=1),
            "most1": np.packbits(f >= 0.5, axis=1), "most2": np.packbits(f >= 0.5, axis=1)}
    med = {m: (f * 7).astype(np.uint16) for m in ("med1", "med2")}
    for m in med:
        med[m][own] = (med[m][own] * OWN_DEPTH_MULTIPLIER).astype(np.uint16)
    out = args.out / "packed"
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "packed.npz", subjects=np.array(subs), nq=nq,
                        library_pairs=np.full(ns, 2_000_000, dtype=np.int64), **bits, **med)
    (out / "query_names.txt").write_text("\n".join(names) + "\n")
    print(f"planted run at {out}: {ns} subjects, {nq} queries, rates {PLANTED}")


if __name__ == "__main__":
    main()
