"""Synthetic runs with planted rates, so both cross-plant reports are checked by output.

`crossplant_contig_report.py` aggregates windows into per-contig fractions, unions libraries
in a shuffled order and scores an own-library ratio against a panel p99. Those are exactly
the shapes that go wrong silently - an axis averaged the wrong way reads plausible. So the
report is run first against a run whose answers are known: the identified set is witnessed
in 90% of libraries, tier 2 in 20%, the shuffles in none, the own library in all, and the
own-library depth is inflated 9x. The report should return those numbers, the nulls should
read zero and the growth swap control should stay flat while the test fires.

The reads mode does the same for `crossplant_report.py`: each library's reads are witnessed
in their own plant at one planted rate and in other plants at another, so the plant matrix
has a known diagonal and a known off-diagonal, and the own-library and shuffled controls
have known values. A matrix that averages the wrong axis reads plausible; one that returns
the planted pair does not.

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
#: reads mode: a library's reads in its own plant, and in any other plant
PLANTED_READS = {"own_plant": 0.6, "other_plant": 0.15}
#: reads mode: the Berlin piles, witnessed at a rate each
PLANTED_BERLIN = {"identified": 0.8, "unident_contig": 0.3, "unplaced": 0.05}
LENGTHS = (200, 300, 500, 1000, 250, 700)
OWN = "SRR23973193"
OWN_DEPTH_MULTIPLIER = 9


def write_packed(out: Path, subs, names, f, boosts, mult):
    """`boosts` is (subject row, query columns) pairs whose depth is inflated `mult` times.

    Growth is scored against the row of the query's OWN library, which for library queries
    is that library's own subject row and not the held-out sample's. A fixture that inflates
    only the held-out row plants nothing the read report measures, and its growth test then
    reads exactly like its control - which is the failure this file exists to catch.
    """
    bits = {"any1": np.packbits(f > 0, axis=1), "any2": np.packbits(f > 0, axis=1),
            "most1": np.packbits(f >= 0.5, axis=1), "most2": np.packbits(f >= 0.5, axis=1)}
    med = {m: (f * 7).astype(np.uint16) for m in ("med1", "med2")}
    for m in med:
        for row, cols in boosts:
            med[m][row, cols] = (med[m][row, cols].astype(int) * mult).clip(0, 65535)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "packed.npz", subjects=np.array(subs), nq=len(names),
                        library_pairs=np.full(len(subs), 2_000_000, dtype=np.int64),
                        **bits, **med)
    (out / "query_names.txt").write_text("\n".join(names) + "\n")


def reads_run(rows, subs, args, rng):
    """Library reads at a planted own-plant and other-plant rate, plus the Berlin piles."""
    plant = {r["run"]: f'{r["project"]}:{r["plant"] or "untreated"}' for r in rows}
    # One source per distinct plant, so the matrix has a diagonal and an off-diagonal.
    seen, sources = set(), []
    for s in subs:
        if s in plant and plant[s] not in seen and len(sources) < args.query_libraries:
            seen.add(plant[s])
            sources.append(s)
    names, meta = [], []
    for acc in sources:
        for seed in (1, 2, 3, 4, 5):
            for i in range(args.per_seed):
                names.append(f"lib|{acc}|{i}|s{seed}")
                meta.append(("lib", acc))
                if i < args.per_seed // 10:
                    names.append(f"shuffled-lib|{acc}|{i}|s{seed}")
                    meta.append(("null", acc))
    for pile, rate in PLANTED_BERLIN.items():
        for seed in (1, 2, 3, 4, 5):
            for i in range(args.per_seed):
                names.append(f"berlin|{pile}|d0|{i}|s{seed}")
                meta.append(("berlin", pile))
    f = np.zeros((len(subs), len(names)))
    for q, (kind, who) in enumerate(meta):
        for j, sub in enumerate(subs):
            if kind == "null":
                p = 0.0
            elif kind == "berlin":
                p = PLANTED_BERLIN[who]
            elif sub == who:
                p = 1.0
            else:
                same = sub in plant and plant.get(who) == plant.get(sub)
                p = PLANTED_READS["own_plant" if same else "other_plant"]
            f[j, q] = float(rng.random() < p)
    # Each source library's own reads are planted deeper in its own library than anywhere
    # else, which is what the growth readout is for.
    boosts = []
    for acc in sources:
        cols = [q for q, (kind, who) in enumerate(meta) if kind == "lib" and who == acc]
        boosts.append((subs.index(acc), np.array(cols)))
    return names, f, boosts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", default="contigs", choices=("contigs", "reads"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--panel", type=Path, default=Path("results/crossplant/panel.tsv"))
    ap.add_argument("--subjects", type=int, default=40)
    ap.add_argument("--contigs", type=int, default=6)
    ap.add_argument("--query-libraries", type=int, default=6,
                    help="reads mode: how many libraries supply queries")
    ap.add_argument("--per-seed", type=int, default=40, help="reads mode: queries per seed")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    rows = list(csv.DictReader(args.panel.open(), delimiter="\t"))
    # A stride over the whole panel, not its first rows: the panel is sorted by project, so
    # a prefix is one plant and a matrix built on it has no off-diagonal to check.
    step = max(1, len(rows) // args.subjects)
    subs = sorted({r["run"] for r in rows[::step][: args.subjects]} | {OWN})
    if args.mode == "reads":
        rng = np.random.default_rng(args.seed)
        names, f, boosts = reads_run(rows, subs, args, rng)
        write_packed(args.out / "packed", subs, names, f, boosts, OWN_DEPTH_MULTIPLIER)
        (args.out / "sample").mkdir(parents=True, exist_ok=True)
        print(f"planted reads run at {args.out}: {len(subs)} subjects, {len(names)} "
              f"queries, rates {PLANTED_READS} and {PLANTED_BERLIN}")
        return
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
    write_packed(args.out / "packed", subs, names, f,
                 [(own, np.arange(nq))], OWN_DEPTH_MULTIPLIER)
    print(f"planted contig run at {args.out}: {ns} subjects, {nq} queries, rates {PLANTED}")


if __name__ == "__main__":
    main()
