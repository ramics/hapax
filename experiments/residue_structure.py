"""Is the residue STRUCTURED or DIFFUSE in the joint space of its fields? (SPEC.md §4.2)

WHY THIS GATES EVERYTHING DOWNSTREAM. Every feature this project has killed, it killed as a
threshold on ONE axis: peptide absence, read-level homology, ORF-likeness, codon usage, GC,
every pair rule. That establishes that no single axis carries a call. It says nothing about
the joint distribution, which has never been measured. If the residue is one featureless
cloud in the joint space as well, the joint-space view in SPEC.md §8.3 has nothing to show
and should not be built; a clean negative here saves that work. If it separates, the
multidimensional measurement has something to measure.

THE AXES, and why these. A read drawn from a molecule that is present in the library at some
depth has 31-mers that recur elsewhere in the library, at a level that is consistent along
the read and shared with its mate. Each competing noise source breaks a different one of
those properties, which is the reason §4.2 gives for expecting joint separation where
thresholds failed:

  - a sequencing error destroys the 31 k-mers spanning it, so one mate carries a run of
    zero-witness k-mers while the other stays covered (high zeros, high mate disagreement,
    a coverage step);
  - a joined library molecule puts its two halves at different depths (step and mate
    disagreement without zeros);
  - a molecule present once in the library has no witnesses at all, evenly (zeros on both
    mates, flat, mates agreeing).

Composition (GC, CpG and UpA observed/expected, trinucleotide entropy, homopolymer) and mean
error probability ride along because they are what the per-read report already carries.

Coverage is EXACT, not sketched: experiments/kmer_coverage (Rust) counts every canonical
31-mer of every read in the five blocks across the whole 8,731,224-pair library and
subtracts each pair's own occurrences, so a count is "witnesses elsewhere in the library".
The sibling thread's `library_abundance` field is a count-min sketch over the first 400,000
pairs that includes the read itself; it is safe for "rare" but cannot resolve flatness or
mate agreement, so it is not reused here.

THE THREE POPULATIONS (SPEC.md §8.3 draws the latter two as reference clouds). Residue
membership is taken by read name from the committed re-cut cascade
(results/subtract_recut/residue_rows_s{1..5}.jsonl, 100,000-pair blocks at skip
(seed-1)x1.5M, the same blocks this script reads):

  - residue: the ~22.9% of pairs the cascade could not identify;
  - subtracted: every other pair in the same block, the known-positive cloud;
  - shuffled null: each residue mate dinucleotide-shuffled (Eulerian-walk, doublet
    frequencies preserved exactly), scored against the library but never added to it, so a
    null read cannot witness its siblings. It holds GC and CpG/UpA fixed and destroys
    sequence identity, which is what "no information" looks like on these axes.

THE PLANTED CONTROL. Reads simulated off a genome absent from this sample's assembly (human
rhinovirus A1B, D00239 - the one clean holdout of the spike ladder, PR for
results/spike_holdout), with this library's measured insert distribution (median 194, p95
363) and error rates (0.00089 R1, 0.00117 R2), at four abundances from 1e-3 to 1e-6 of
library pairs. Each abundance arm is an independently substituted copy of the genome so that
arms share no 31-mers and cannot witness each other. The arms are added to the library pass,
so a planted read is witnessed by exactly as many siblings as its abundance implies. This is
the only arm of the experiment with a known answer, and it is what says whether a cluster
this finds is the kind of thing a real low-abundance genome makes.

WHAT COUNTS AS STRUCTURE, decided before looking (three readouts, increasingly demanding):

  1. Separability. Cross-validated AUC for residue against subtracted, in the joint space
     versus the best single axis. The joint space earning more than the best threshold is
     the minimum claim of §4.2 and it is a number, not a picture.
  2. Clustering beyond a correlated Gaussian. k-means for k = 1..8 on z-scored axes, with
     the dispersion compared against a parametric bootstrap drawn with the SAME mean and
     covariance as the data (gap statistic). A dozen correlated fields produce apparent
     clusters when they are one measurement wearing a hat, and that null removes exactly
     that: only clustering beyond the covariance counts.
  3. Recovery of labels the clustering never saw. Does a cluster pick out the planted arms,
     or the residue piles the cascade defines, at better than chance?

A diffuse answer is: AUC at or near the best single axis, gap statistic flat in k, planted
arms indistinguishable from background. That is a publishable negative and is reported as
one.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import subprocess
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BLOCK = 1_500_000      # pairs skipped per seed, the cascade's convention
PAIRS = 100_000        # pairs per block, as the cascade used
TRIM_R2 = 14           # mate 2's non-template prefix (experiments/berlin_io.py)
INSERT_Q = [(0.50, 194), (0.75, 246), (0.95, 363), (0.99, 485), (1.0, 600)]
ERR1, ERR2 = 0.00089, 0.00117
READ_LEN = 107
LIBRARY_PAIRS = 8_731_224
ARMS = [("1e-3", 1e-3), ("1e-4", 1e-4), ("1e-5", 1e-5), ("1e-6", 1e-6)]
ALT = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}


# ---------------------------------------------------------------- sequence helpers

def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def dinuc_shuffle(s: str, rng: random.Random) -> str:
    """Doublet-preserving shuffle by an Eulerian walk on the nucleotide graph.

    A plain shuffle would change CpG and UpA observed/expected, which are two of the axes,
    so the null would be a different composition rather than the same composition without
    information. Altschul-Erickson: random Eulerian circuit on the doublet multigraph.
    """
    if len(s) < 4:
        return s
    edges: dict[str, list[str]] = {}
    for a, b in zip(s, s[1:]):
        edges.setdefault(a, []).append(b)
    last = s[-1]
    # Pick, for each vertex, a last-edge choice that keeps a spanning tree to `last`
    # (connectivity condition); retry on failure, which is cheap at this length.
    for _ in range(50):
        order = {v: rng.sample(e, len(e)) for v, e in edges.items()}
        reach = {last}
        ok = True
        for v in order:
            if v == last:
                continue
            # follow the chosen last edges toward `last`
            seen, cur = set(), v
            while cur not in reach:
                if cur in seen or cur not in order:
                    ok = False
                    break
                seen.add(cur)
                cur = order[cur][-1]
            if not ok:
                break
            reach |= seen | {v}
        if not ok:
            continue
        out, cur = [s[0]], s[0]
        avail = {v: list(e) for v, e in order.items()}
        for _ in range(len(s) - 1):
            if not avail.get(cur):
                break
            nxt = avail[cur].pop(0)
            out.append(nxt)
            cur = nxt
        if len(out) == len(s):
            return "".join(out)
    return s


def composition(s1: str, s2: str, q1: str, q2: str) -> dict:
    """GC, trinucleotide entropy, CpG/UpA o/e, homopolymer and mean error probability.

    Pair-level, over both mates concatenated, because the row the report emits is per pair.
    """
    s = s1 + s2
    n = len(s)
    counts = Counter(s)
    gc = (counts["G"] + counts["C"]) / max(n, 1)
    tri = Counter(s[i:i + 3] for i in range(n - 2) if set(s[i:i + 3]) <= set("ACGT"))
    tot = sum(tri.values()) or 1
    ent = -sum((c / tot) * math.log(c / tot, 2) for c in tri.values()) / 6.0
    di = Counter(s[i:i + 2] for i in range(n - 1))
    f = {b: counts[b] / max(n, 1) for b in "ACGT"}

    def oe(xy: str) -> float:
        exp = f[xy[0]] * f[xy[1]] * max(n - 1, 1)
        return di[xy] / exp if exp > 0 else 0.0

    hp, run, prev = 1, 1, ""
    for b in s:
        run = run + 1 if b == prev else 1
        prev = b
        hp = max(hp, run)
    qs = q1 + q2
    err = sum(10 ** (-(ord(c) - 33) / 10) for c in qs) / max(len(qs), 1)
    return {"gc": round(gc, 4), "entropy3": round(ent, 4), "CpG_oe": round(oe("CG"), 4),
            "UpA_oe": round(oe("TA"), 4), "homopolymer": hp,
            "mean_error_prob": round(err, 6)}


def read_pairs(f1: Path, f2: Path, limit: int, skip: int = 0):
    """Mate 2's non-template prefix trimmed, same convention as experiments/berlin_io.py."""
    with gzip.open(f1, "rt") as h1, gzip.open(f2, "rt") as h2:
        name = s1 = s2 = None
        for i, (l1, l2) in enumerate(zip(h1, h2)):
            rec = i // 4
            if rec < skip:
                continue
            ph = i % 4
            if ph == 0:
                name = l1[1:].split()[0]
            elif ph == 1:
                s1, s2 = l1.strip(), l2.strip()[TRIM_R2:]
            elif ph == 3:
                yield name, s1, s2, l1.strip(), l2.strip()[TRIM_R2:]
                if rec + 1 - skip >= limit:
                    return


def simulate_pairs(genome: str, n: int, rng: random.Random, tag: str):
    """Pairs off `genome` with this library's measured insert and error profile."""
    def errs(s: str, e: float) -> str:
        out = list(s)
        for i, b in enumerate(out):
            if b in ALT and rng.random() < e:
                out[i] = rng.choice(ALT[b])
        return "".join(out)

    out = []
    for i in range(n):
        u = rng.random()
        ins = INSERT_Q[-1][1]
        for q, v in INSERT_Q:
            if u <= q:
                ins = v
                break
        ins = min(max(ins, READ_LEN + 1), len(genome) - 1)
        st = rng.randrange(0, max(len(genome) - ins, 1))
        frag = genome[st:st + ins]
        if rng.random() < 0.5:
            frag = revcomp(frag)
        out.append((f"{tag}|{i}", errs(frag[:READ_LEN], ERR1),
                    errs(revcomp(frag)[:READ_LEN], ERR2)))
    return out


def substitute(seq: str, rate: float, rng: random.Random) -> str:
    return "".join(rng.choice(ALT[b]) if b in ALT and rng.random() < rate else b for b in seq)


# ---------------------------------------------------------------- prep

def cmd_prep(args: argparse.Namespace) -> None:
    """Write the planted-control and shuffled-null read files the Rust counter consumes."""
    rng = random.Random(args.seed)
    args.work.mkdir(parents=True, exist_ok=True)

    genome = "".join(l.strip() for l in open(args.genome) if not l.startswith(">")).upper()
    with open(args.work / "spike.tsv", "w") as h:
        for tag, rate in [(f"{t}_r{r}", v) for t, v in ARMS
                          for r in range(1, args.arm_replicates + 1)]:
            n = max(int(round(rate * LIBRARY_PAIRS)), 1)
            # Each arm is its own 25%-substituted copy: at 31 nt two arms share a k-mer with
            # probability 0.75**31 ~ 1e-4, so arms cannot witness each other and each arm's
            # coverage is its own abundance and nothing else.
            g = substitute(genome, 0.25, rng)
            for nm, r1, r2 in simulate_pairs(g, n, rng, f"spike_{tag}"):
                h.write(f"{nm}\t{r1}\t{r2}\n")

    if args.stage == "spike":
        print("spike arms written")
        return

    # Shuffled null: a sample of residue pairs per seed block, doublet-shuffled.
    with open(args.work / "null.tsv", "w") as h:
        for seed in range(1, 6):
            names = set(open(ROOT / f"data/residue_structure/residue_names_s{seed}.txt")
                        .read().split())
            kept = 0
            for nm, s1, s2, _q1, _q2 in read_pairs(args.r1, args.r2, PAIRS,
                                                   (seed - 1) * BLOCK):
                if nm in names and kept < args.null_per_seed:
                    h.write(f"null_s{seed}|{nm}\t{dinuc_shuffle(s1, rng)}"
                            f"\t{dinuc_shuffle(s2, rng)}\n")
                    kept += 1
    print("prep written to", args.work)


# ---------------------------------------------------------------- composition pass

def cmd_compose(args: argparse.Namespace) -> None:
    """Composition and quality for every pair of the five blocks, plus spike and null."""
    args.work.mkdir(parents=True, exist_ok=True)
    with open(args.work / "composition.tsv", "w") as h:
        keys = ["gc", "entropy3", "CpG_oe", "UpA_oe", "homopolymer", "mean_error_prob"]
        h.write("name\t" + "\t".join(keys) + "\n")
        for seed in range(1, 6):
            for nm, s1, s2, q1, q2 in read_pairs(args.r1, args.r2, PAIRS,
                                                 (seed - 1) * BLOCK):
                c = composition(s1, s2, q1, q2)
                h.write(nm + "\t" + "\t".join(str(c[k]) for k in keys) + "\n")
        for f in ("spike.tsv", "null.tsv"):
            for line in open(args.work / f):
                nm, r1, r2 = line.rstrip("\n").split("\t")
                c = composition(r1, r2, "I" * len(r1), "I" * len(r2))
                h.write(nm + "\t" + "\t".join(str(c[k]) for k in keys) + "\n")
    print("composition written")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--genome", type=Path, default=ROOT / "data/spike/D00239.fa")
    ap.add_argument("--work", type=Path, default=ROOT / "data/residue_structure")
    ap.add_argument("--null-per-seed", type=int, default=4000)
    ap.add_argument("--arm-replicates", type=int, default=5,
                    help="independent simulations per abundance, so the planted control "
                         "carries spread like every other headline in this repo")
    ap.add_argument("stage", choices=["prep", "spike", "compose"])
    args = ap.parse_args()
    t0 = time.perf_counter()
    {"prep": cmd_prep, "spike": cmd_prep, "compose": cmd_compose}[args.stage](args)
    print(f"{args.stage} {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
