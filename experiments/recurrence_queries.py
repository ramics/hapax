"""The query set for the cross-sample recurrence pass, stratified by depth.

WHAT THE RECURRENCE MEASUREMENT IS. For each read of this sample, is its sequence present in
each of the other libraries of PRJNA948850 (`experiments/panel_inventory.py` pins the panel),
by exact canonical 31-mer presence through `experiments/kmer_coverage`'s query-only path. One
row per read per subject library; occupancy is then the number of subjects that witness it,
and the deliverable is the distribution of that, never a threshold on it.

WHY THE QUERY SET IS STRATIFIED BY DEPTH RATHER THAN MATCHED ON IT. Recurrence across samples
and depth within this one are the same variable wearing a hat: an abundant molecule is
abundant everywhere, so a residue pile at median depth 1 compared against identified reads at
median depth 87 would show less recurrence for a reason that is arithmetic rather than
biological. The brainstorm called for an abundance-matched control. Matching throws away the
part of the data that answers the question, so this instead selects a FIXED NUMBER of pairs
in every (pile, depth bin) cell, which makes recurrence reportable as a curve against depth
with the identified piles' curve drawn beside it, per the repo's preference for continuous
readouts over thresholds. A matched comparison is then one weighted average away, and the
curve says whether the two piles differ at the same depth at all - which is the actual
question and is not visible in a single matched number.

THE DEPTH THE CELLS ARE BUILT FROM IS THIS LIBRARY'S OWN. It comes from the library-block
path of the same tool, which subtracts each pair's own occurrences, so it is "occurrences
elsewhere in this run" - the quantity `results/residue_census/` §3 reports. The cross-sample
scoring is the block-9 query-only path, which subtracts nothing because a foreign library has
no copy of the query to subtract. Those are different paths of one tool, and the rule learned
in `cross_library_witness.py` is that groups being COMPARED must go through the same path:
here every pile and the null go through block 9 against every subject, and depth is a
covariate computed identically for all of them. Nothing is compared across the two paths.

PILES, from the whole-library label table (`library_labels.py`, validated against the cascade
in `results/residue_census/` §2) rather than from the cascade's residue FASTAs, which live in
another thread's working store:
  - `unplaced`: neither mate matches SILVA or any of this sample's own contigs at >=40 bases.
    15.71% of the library, the pile the residue-only assembly was built from, and the one at
    median depth 1. STATED BIAS: it carries no read-level Kraken2 or DIAMOND, so it is a
    superset of the cascade's residue and holds some identifiable pairs; that can only make
    it look MORE like identified material, i.e. it biases against any difference found.
  - `unident_contig`: residue sitting on one of the sample's own unidentifiable contigs.
  - `identified`: a non-rRNA own contig that an index identified. The known-positive.
  - `rrna`: SILVA or an rRNA contig. Reported separately and never pooled, because shared
    ribosomal sequence is exactly what inflated the 71.2% figure this experiment replaces -
    it recurs in every library by conservation, so pooling it would manufacture recurrence.
  - `shuffled`: the null. Each selected read is shuffled with its mononucleotide composition
    preserved, so what it measures is how often composition alone finds a 31-mer witness in a
    foreign library of the same material. It is scored through the same block-9 path.

WHY A SHUFFLE AND NOT A RANDOM SEQUENCE. A uniform random 31-mer is absent from everything,
which would make the null trivially zero and tell us nothing; the composition-preserving
shuffle is the null this repo already uses for every field, and it reads 0.00 on the artefact
panel, so a non-zero reading here would be information about the library rather than noise.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

import library_facts

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from berlin_io import read_pairs  # noqa: E402
from subtract_known import provenance  # noqa: E402

#: Class codes as `library_label_classes.py` writes them.
CODES = ("none", "rrna", "identified", "unident", "other_ctg")

#: Depth bins on this library's own median depth for mate 1, left-closed. The edges are the
#: ones the census's table already distinguishes - 0 (unwitnessed), 1-2 (the unplaced pile's
#: median), 3-9, 10-99 (identified reads sit at 87), 100+ - so the cells line up with a
#: published table instead of introducing a second set of bins for the same axis.
DEPTH_EDGES = (0, 1, 3, 10, 100)
DEPTH_NAMES = ("d0", "d1-2", "d3-9", "d10-99", "d100+")


def bin_of(depth: int) -> str:
    i = 0
    for j, edge in enumerate(DEPTH_EDGES):
        if depth >= edge:
            i = j
    return DEPTH_NAMES[i]


def pile_of(c1: int, c2: int) -> str:
    """One pile per pair, from the two mates' classes.

    Order matters and is the census's: a pair with any ribosomal mate is ribosomal, because
    that is the conserved sequence whose recurrence would otherwise leak into every other
    pile. Then identified, then the unidentifiable-contig pile, then unplaced.
    """
    cls = {CODES[c1], CODES[c2]}
    if "rrna" in cls:
        return "rrna"
    if "identified" in cls:
        return "identified"
    if "other_ctg" in cls:
        return "other_ctg"
    if "unident" in cls:
        return "unident_contig"
    return "unplaced"


def shuffled(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def write_tsv(path: Path, rows: list[tuple[str, str, str]]) -> None:
    with open(path, "w") as h:
        for name, s1, s2 in rows:
            h.write(f"{name}\t{s1}\t{s2}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="which 100,000-pair block of the library, as everywhere else in "
                         "this repo: block 1 = seed 1, and blocks are flow-cell regions "
                         "rather than a reshuffle")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--classes", type=Path,
                    default=ROOT / "data/work/library_label_classes.npz")
    ap.add_argument("--block-pairs", type=int, default=100_000)
    ap.add_argument("--per-cell", type=int, default=2_000,
                    help="pairs per (pile, depth bin) cell. The cost of the cross-sample "
                         "pass is one scan of each subject library, which is independent of "
                         "the query count, so this is bounded by memory rather than by time")
    ap.add_argument("--null-per-pile", type=int, default=2_000)
    ap.add_argument("--tool", type=Path,
                    default=ROOT / "experiments/kmer_coverage/target/release/kmer_coverage")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/recurrence"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/recurrence")
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["trim_r2"])

    for path in (args.r1, args.r2, args.classes, args.tool):
        if not path.exists():
            raise SystemExit(f"{path} is missing; nothing can be selected without it")
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    rng = random.Random(args.seed)

    # ---------------------------------------------------------------- depth, own-subtracted
    # One pass of the tool over the whole library with this block as a library block, which
    # is the path that subtracts a pair's own occurrences. This is the census's instrument
    # and its numbers are reproduced below as a validation before anything is selected.
    skip = args.block_pairs * (args.seed - 1)
    cov = args.work / f"cov_block{args.seed}.tsv"
    t_cov = time.perf_counter()
    if not cov.exists():
        with open(cov, "w") as out, open(args.work / f"cov_s{args.seed}.err", "w") as err:
            subprocess.run([str(args.tool), str(args.r1), str(args.r2), str(args.trim_r2),
                            str(args.block_pairs), str(skip)],
                           stdout=out, stderr=err, check=True)
    seconds_depth = time.perf_counter() - t_cov

    depth: dict[int, int] = {}
    with open(cov) as h:
        hdr = h.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for need in ("name", "med1", "block"):
            if need not in ix:
                raise SystemExit(f"{cov} has no {need} column; the tool's format changed")
        for line in h:
            f = line.rstrip("\n").split("\t")
            depth[int(f[ix["name"]].split(".")[1])] = int(f[ix["med1"]])
    if not depth:
        raise SystemExit(f"{cov} holds no rows for block {args.seed}")

    # ---------------------------------------------------------------- piles, from the labels
    d = np.load(args.classes)
    spot = d["spot"].astype(np.int64)
    c1, c2 = d["class1"], d["class2"]
    order = np.argsort(spot)
    spot, c1, c2 = spot[order], c1[order], c2[order]

    cells: dict[tuple[str, str], list[int]] = {}
    for sp, depth_i in depth.items():
        j = int(np.searchsorted(spot, sp))
        if j >= spot.size or spot[j] != sp:
            continue                      # a pair the label table does not cover
        key = (pile_of(int(c1[j]), int(c2[j])), bin_of(depth_i))
        cells.setdefault(key, []).append(sp)

    # VALIDATION BEFORE SELECTION, per the repo's rule that the control is read first: the
    # pile sizes and the unplaced pile's depth must reproduce the published census. If they
    # do not, the label table and this block are not the same instrument and nothing below
    # means anything.
    pile_counts = Counter(k[0] for k, v in cells.items() for _ in v)
    n_block = sum(pile_counts.values())
    unplaced_depths = [depth[s] for k, v in cells.items() if k[0] == "unplaced" for s in v]
    check = {"block_pairs_labelled": n_block,
             "pile_fractions": {k: round(v / n_block, 4) for k, v in pile_counts.items()},
             "unplaced_median_depth": float(np.median(unplaced_depths)),
             "unplaced_frac_depth_0": round(float(np.mean(np.array(unplaced_depths) == 0)), 4),
             "unplaced_frac_depth_lt5": round(float(np.mean(np.array(unplaced_depths) < 5)), 4)}
    print(json.dumps(check, indent=2), flush=True)
    if not 0.12 <= pile_counts["unplaced"] / n_block <= 0.20:
        raise SystemExit(
            f"the unplaced pile is {pile_counts['unplaced'] / n_block:.1%} of this block "
            f"against 15.71% of the library in results/residue_census/ §2. That is not a "
            f"block-to-block difference; it means the label table is being indexed wrong.")
    if not 0.0 <= check["unplaced_median_depth"] <= 3.0:
        raise SystemExit(
            f"the unplaced pile's median depth is {check['unplaced_median_depth']} against "
            f"1.0 published. The depth column or the pile assignment is wrong.")

    # ---------------------------------------------------------------- the query set
    wanted: dict[int, tuple[str, str]] = {}
    for (pile, dbin), spots in sorted(cells.items()):
        rng.shuffle(spots)
        for sp in spots[:args.per_cell]:
            wanted[sp] = (pile, dbin)

    rows: list[tuple[str, str, str]] = []
    null_quota = Counter()
    for name, s1, s2, _, _ in read_pairs(args.r1, args.r2, args.block_pairs,
                                         trim_r2=args.trim_r2, skip=skip):
        sp = int(name.split(".")[1])
        if sp not in wanted:
            continue
        pile, dbin = wanted[sp]
        rows.append((f"{pile}|{dbin}|{depth[sp]}|{sp}|s{args.seed}", s1, s2))
        # The null is drawn from the same reads rather than from a separate sample, so it
        # carries the same composition, the same length and the same error content.
        if null_quota[pile] < args.null_per_pile:
            null_quota[pile] += 1
            rows.append((f"shuffled-{pile}|{dbin}|{depth[sp]}|{sp}|s{args.seed}",
                         shuffled(s1, rng), shuffled(s2, rng)))

    if len(rows) < 0.5 * len(wanted):
        raise SystemExit(
            f"only {len(rows)} of {len(wanted)} selected pairs were found in the reads. The "
            f"block offset and the label table disagree, and a silently short query set "
            f"would still produce a tidy distribution.")

    q = args.out / f"queries_s{args.seed}.tsv"
    write_tsv(q, rows)
    cell_sizes = Counter(name.split("|")[0] + "|" + name.split("|")[1] for name, _, _ in rows)
    manifest = {
        "script": "experiments/recurrence_queries.py", **provenance(__file__),
        "parameters": {k: str(v) for k, v in vars(args).items()},
        "query_pairs": len(rows),
        "real_pairs": sum(1 for n, _, _ in rows if not n.startswith("shuffled-")),
        "null_pairs": sum(1 for n, _, _ in rows if n.startswith("shuffled-")),
        "cells": dict(sorted(cell_sizes.items())),
        "validation_against_residue_census": check,
        "operating_points": {"k": 31, "class_split": "40 matching bases, SILVA or own contig",
                             "depth": "median over mate 1's scorable 31-mers, own "
                                      "occurrences subtracted",
                             "trim_r2": args.trim_r2},
        "timings": {"depth_pass_seconds": round(seconds_depth, 1),
                    "total_seconds": round(time.perf_counter() - t0, 1)},
    }
    (args.out / f"queries_manifest_s{args.seed}.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(f"{len(rows):,} query pairs -> {q}")
    for k, v in sorted(cell_sizes.items()):
        print(f"  {k:28s} {v:5d}")


if __name__ == "__main__":
    main()
