"""Does a read find its own mate? Mate linkage as ground truth for "reads find their friends".

THE IDEA, AND WHY IT IS NOT THE EXPERIMENT WE ALREADY RAN. `results/residue_structure/` §4
measured whether one genome's reads find each other in the joint space, by planting a genome
at a known abundance and asking what fraction of a planted read's 10 nearest neighbours came
from the same source. It works to about 1 in 100,000 of the library and dies there. But it
needed a PLANTED arm to know the answer, it treated a PAIR as one point, and it therefore
never used the one piece of ground truth the real library already carries: **the two mates of
a pair came from the same molecule.**

So: make each MATE a point, and ask whether a read's own mate is among its k nearest
neighbours. That is a label-free "same molecule" test, on real sequence, at every abundance
present in the library, with no simulation anywhere in it. Imo's judgement is that mate
linkage is crucial; this is the cheapest measurement of whether the friendship geometry can
see it at all.

WHAT A POSITIVE WOULD MEAN AND WHAT IT WOULD NOT. Recovery above chance means the axes carry
molecule-level information, not merely abundance-level. It does NOT mean the geometry can
group reads from the same GENOME: two mates share an insert, so they are far more alike than
two reads a kilobase apart in the same genome. This measurement is therefore an UPPER BOUND
on friendship resolution, and the planted-genome purity in §4 is the corresponding lower
bound. Stating it as a bound is the point; quoting it as the grouping rate would be the
error.

THREE CONFOUNDS, EACH MEASURED RATHER THAN ASSUMED AWAY.

1. **26-27% of pairs in this library overlap their own mate** (`results/quality_profile/`).
   An overlapping pair's two mates contain literally the same bases, so their coverage
   profiles are correlated by construction and not by biology. Recovery is therefore reported
   SPLIT by whether the mates share a 31-mer, and the non-overlapping half is the honest
   number.
2. **The mates are not interchangeable.** Mate 2 is 93 nt against mate 1's 107 after the
   measured 14 nt trim (`results/nontemplate_prefix/`), so it carries 62 scorable 31-mers
   against 76 and its coverage axes are noisier. Recovery is reported per direction - mate 1
   finding mate 2, and mate 2 finding mate 1 - and beside a length-matched control that
   truncates mate 1 to 93 nt. CLAUDE.md requires this of any pair-symmetric field and the
   two directions here are expected to differ.
3. **Abundance does most of the work in the joint space already** (coverage-only purity 0.918
   against composition-only 0.292 at 1 in 10,000). So recovery is reported BY COVERAGE
   DECILE, which is what makes it comparable to the 1-in-100,000 wall, and with the axis
   groups split so the share carried by composition alone is visible.

TIES ARE BROKEN AT RANDOM, and that is not a detail. These axes are exactly degenerate for a
mate with no witnesses anywhere - coverage zero, spread zero, every k-mer unwitnessed - so
thousands of mates share a coordinate. sklearn returns neighbours in index order, and since
a pair's two mates are adjacent in any natural ordering, index-order ties would hand every
one of them its own partner and manufacture the entire result. This is the same failure that
inflated a nearest-neighbour control in this project once already.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts

from berlin_io import read_pairs
from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

#: Per-mate axes. The coverage four come from the Rust counter's per-mate columns; the
#: composition six are computed here on that mate's own sequence.
COVERAGE_AXES = ["log_cov", "cov_spread", "zero_frac", "cov_step"]
COMPOSITION_AXES = ["gc", "entropy3", "CpG_oe", "UpA_oe", "homopolymer", "mean_error_prob"]
AXES = COVERAGE_AXES + COMPOSITION_AXES

#: The tool's output header, asserted rather than assumed: reading this table by hard-coded
#: column index once produced a fraction of -0.71 in `cross_library_witness.py`.
HEADER = ("block name n1 n2 med1 med2 q10_1 q90_1 q10_2 q90_2 zero1 zero2 step1 step2 "
          "zrun1 zrun2 stepnz1 stepnz2 library_pairs").split()


def entropy3(seq: str) -> float:
    if len(seq) < 3:
        return 0.0
    counts: dict[str, int] = {}
    for i in range(len(seq) - 2):
        counts[seq[i:i + 3]] = counts.get(seq[i:i + 3], 0) + 1
    n = sum(counts.values())
    h = -sum((c / n) * math.log(c / n) for c in counts.values())
    return h / math.log(min(64, n)) if n > 1 else 0.0


def dinuc_oe(seq: str, a: str, b: str) -> float:
    """Observed over expected for one dinucleotide, the standard CpG-style ratio."""
    n = len(seq)
    if n < 2:
        return 0.0
    na = seq.count(a)
    nb = seq.count(b)
    obs = sum(1 for i in range(n - 1) if seq[i] == a and seq[i + 1] == b)
    exp = na * nb / n
    return obs / exp if exp > 0 else 0.0


def homopolymer(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def mean_error_prob(qual: str) -> float:
    if not qual:
        return 0.0
    return st.mean(10 ** (-(ord(c) - 33) / 10) for c in qual)


def composition(seq: str, qual: str) -> list[float]:
    acgt = sum(1 for c in seq if c in "ACGT") or 1
    return [sum(1 for c in seq if c in "GC") / acgt,
            entropy3(seq),
            dinuc_oe(seq, "C", "G"),
            dinuc_oe(seq, "T", "A"),            # UpA on the RNA convention
            homopolymer(seq) / max(len(seq), 1),
            mean_error_prob(qual)]


def coverage_axes(f: dict[str, str], mate: int) -> list[float] | None:
    """The four coverage axes for one mate, or None if it has no scorable k-mer."""
    try:
        n = int(f[f"n{mate}"])
        med = float(f[f"med{mate}"])
        q10 = float(f[f"q10_{mate}"])
        q90 = float(f[f"q90_{mate}"])
        zero = float(f[f"zero{mate}"])
        step = float(f[f"step{mate}"])
    except (KeyError, ValueError):
        return None
    if n <= 0 or zero != zero or med != med:
        return None
    if not 0.0 <= zero <= 1.0:
        raise SystemExit(f"zero fraction {zero} for mate {mate} is impossible, so a column "
                         f"is being read wrong")
    return [math.log1p(max(med, 0.0)),
            math.log2((q90 + 1.0) / (q10 + 1.0)),
            zero,
            math.log1p(max(step, 0.0))]


#: The k the mate-sharing split uses, named once so the bias statement and the
#: split cannot disagree about it.
KMER_K = 31


def kmers(seq: str, k: int = KMER_K) -> set[str]:
    """Canonical k-mers, for the mate-overlap test only."""
    comp = str.maketrans("ACGT", "TGCA")
    out = set()
    for i in range(len(seq) - k + 1):
        s = seq[i:i + k]
        if set(s) - set("ACGT"):
            continue
        out.add(min(s, s.translate(comp)[::-1]))
    return out


def recovery(X: np.ndarray, partner: np.ndarray, ks: list[int],
             rng: np.random.Generator, subset: np.ndarray | None = None) -> dict:
    """Is each point's partner among its k nearest neighbours?

    `partner[i]` is the index of i's own mate. Ties are permuted before the query for the
    reason in the module docstring: a pair's two mates are adjacent in index order, so
    index-order tie-breaking would return the partner first for every degenerate point and
    produce a perfect result from nothing.
    """
    n = len(X)
    # AN EMPTY SUBSET IS A CENSUS, NOT AN ERROR. `mates_overlap` and `mates_disjoint`
    # partition the pairs, and on a library whose insert is longer than twice the read
    # one side of that partition is empty - which is a fact about the library, measured.
    # Asked for zero query points, the neighbour search raised, so a legitimate reading
    # of a short-insert library stopped the stage. The count is reported with the reason,
    # because an empty table and a test that never ran must not look the same.
    rows_in = np.arange(n) if subset is None else np.flatnonzero(subset)
    if len(rows_in) == 0 or n < 2:
        out = {str(k): {"recovery": None, "chance": None, "enrichment": None}
               for k in ks}
        out["n_queries"] = int(len(rows_in))
        out["n_points"] = int(n)
        out["state"] = ("no query points in this subset, so the recovery rate is not "
                        "defined here" if n >= 2 else
                        "fewer than two points, so no point has a neighbour")
        return out
    perm = rng.permutation(n)
    inv = np.empty(n, dtype=np.int64)
    inv[perm] = np.arange(n)
    Xp = X[perm]
    kmax = min(max(ks), n - 1)
    nn = NearestNeighbors(n_neighbors=kmax + 1).fit(Xp)
    rows = np.arange(n) if subset is None else np.flatnonzero(subset)
    _d, idx = nn.kneighbors(Xp[inv[rows]])
    want = inv[partner[rows]][:, None]
    out = {}
    for k in ks:
        kk = min(k, kmax)
        hit = (idx[:, 1:kk + 1] == want).any(axis=1)
        rate = float(hit.mean())
        chance = kk / (n - 1)
        out[str(k)] = {"recovery": round(rate, 5), "chance": round(chance, 6),
                       "enrichment": round(rate / chance, 2) if chance else None}
    out["n_queries"] = int(len(rows))
    out["n_points"] = int(n)
    return out


def standardise(X: np.ndarray) -> np.ndarray:
    """Z-score the axes.

    Not cosmetic: the axes run from a fraction in [0,1] to a log coverage near 11, so a
    Euclidean neighbourhood on the raw matrix would be a statement about log coverage and
    nothing else, which is exactly the confound this experiment exists to separate.
    """
    mu, sd = X.mean(0), X.std(0)
    sd = sd.copy()
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def build_space(args, rows, tag: str):
    """Score one set of pairs against the whole library and return the per-mate matrix.

    Factored out because the length-matched control needs the identical pipeline on mate 1
    truncated to mate 2's length. Running that arm through a second code path would make the
    control and the test two instruments, which this project has a rule against.
    """
    q = args.work / f"query_{tag}_s{args.seed}.tsv"
    with open(q, "w") as h:
        for name, s1, s2, _q1, _q2 in rows:
            h.write(f"{name}\t{s1}\t{s2}\n")
    seqs = {r[0]: (r[1], r[2], r[3], r[4]) for r in rows}

    # Query-only (the tool's 7th argument): this path does NOT subtract a pair's own
    # occurrences from the library counts, which is what makes the coverage a statement
    # about the library rather than about the pair.
    tsv = args.work / f"coverage_{tag}_s{args.seed}.tsv"
    t = time.perf_counter()
    with open(tsv, "w") as out, open(args.work / f"tool_{tag}.err", "w") as err:
        subprocess.run([str(args.tool), str(args.r1), str(args.r2), str(args.trim_r2),
                        str(args.library_pairs), "0", "-", str(q)],
                       stdout=out, stderr=err, check=True)
    seconds = time.perf_counter() - t

    points: list[list[float]] = []
    owner: list[str] = []
    which: list[int] = []
    dropped = 0
    with open(tsv) as h:
        head = h.readline().rstrip("\n").split("\t")
        if head != HEADER:
            raise SystemExit(f"{tsv}: unexpected header {head[:6]}... - the tool's output "
                             f"format changed and every column below would be wrong")
        for line in h:
            f = dict(zip(HEADER, line.rstrip("\n").split("\t")))
            if f["block"] != "9":
                continue
            name = f["name"]
            if name not in seqs:
                continue
            s1, s2, q1, q2 = seqs[name]
            a = coverage_axes(f, 1)
            b = coverage_axes(f, 2)
            if a is None or b is None:
                dropped += 1
                continue
            points.append(a + composition(s1, q1))
            owner.append(name)
            which.append(1)
            points.append(b + composition(s2, q2))
            owner.append(name)
            which.append(2)

    if len(points) < 200:
        raise SystemExit(f"only {len(points)} mate points survived in {tag}; nothing can "
                         f"be measured")
    return np.asarray(points, dtype=np.float64), owner, which, dropped, seconds


def _bias(rows, k: int) -> str:
    """The bias statement, with THIS library's mate lengths in it.

    The sentence used to end "Mate 2 is 93 nt against mate 1's 107 and carries 62
    scorable 31-mers against 76", which is the discriminating dataset after its 14 nt
    non-template trim. Every library's linkage output carried it, and the web view prints
    it verbatim under the mate-recovery table, so a 2x151 library was told its mates were
    93 and 107 nt. The asymmetry is the point of the sentence and is real wherever it
    exists, so it is measured here and the clause is dropped where the two mates are the
    same length, which is itself the thing a reader needs to know.
    """
    r1 = st.median([len(r[1]) for r in rows]) if rows else 0
    r2 = st.median([len(r[2]) for r in rows]) if rows else 0
    head = ("a pair's two mates come from one insert and many of them overlap, so this "
            "is an UPPER BOUND on friendship resolution between reads of one genome, "
            "not an estimate of it. The planted-genome purity in "
            "results/residue_structure/ is the lower bound. ")
    if r1 == r2:
        return head + (
            f"Both mates are {int(r1)} nt here, carrying {max(0, int(r1) - k + 1)} "
            f"scorable {k}-mers each, so the two directions are the same measurement on "
            "this library; they are still reported apart, because a mate can map at a "
            "different rate without being shorter")
    return head + (
        f"Mate 2 is {int(r2)} nt against mate 1's {int(r1)} and carries "
        f"{max(0, int(r2) - k + 1)} scorable {k}-mers against "
        f"{max(0, int(r1) - k + 1)}, so the two directions are different measurements "
        "and are reported apart")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=20_000,
                    help="query pairs for this seed. 2x this many points go in the space")
    ap.add_argument("--block", type=int, default=None,
                    help="pairs skipped per seed, matching subtract_known.py, so the blocks "
                         "are the same ones every other experiment here uses")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="measured, not detected: results/nontemplate_prefix/")
    ap.add_argument("--ks", type=int, nargs="+", default=[1, 5, 10, 50])
    ap.add_argument("--tool", type=Path,
                    default=ROOT / "experiments/kmer_coverage/target/release/kmer_coverage")
    ap.add_argument("--library-pairs", type=int, default=2_000,
                    help="library pairs the tool also scores as its own block 1. Not a "
                         "control - the tool subtracts a library pair's own occurrences "
                         "and does not for the query-only path - and present only because "
                         "the tool requires a block specification")
    ap.add_argument("--length-match", action="store_true", default=True,
                    help="run the length-matched arm, which truncates mate 1 to mate 2's "
                         "length and pays a second whole-library pass. On by default "
                         "because CLAUDE.md requires a pair-symmetric field to ship either "
                         "per mate or with this control, and this ships both")
    ap.add_argument("--no-length-match", dest="length_match", action="store_false")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/matefriends"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/mate_friends")
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["trim_r2", "block"])

    if not args.tool.exists():
        raise SystemExit(f"{args.tool} is not built; cargo build --release in its directory")
    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    skip = args.block * (args.seed - 1)

    # ---- the query set, written in the tool's query-only format
    rows = list(read_pairs(args.r1, args.r2, args.pairs, trim_r2=args.trim_r2, skip=skip))
    if not rows:
        raise SystemExit(f"no pairs read from {args.r1}")

    X, owner, which, dropped, cover_seconds = build_space(args, rows, "full")
    which_arr = np.asarray(which)
    # partner[i] is the index of i's own mate. The two mates of a pair were appended
    # consecutively, so this is exact rather than a lookup.
    partner = np.arange(len(X)) ^ 1
    Z = standardise(X)
    seqs = {r[0]: (r[1], r[2], r[3], r[4]) for r in rows}

    rng = np.random.default_rng(args.seed)
    result: dict = {"all": recovery(Z, partner, args.ks, rng)}

    # ---- CONFOUND 1: mates that overlap share literal sequence
    pair_names = owner[0::2]
    overlap = {}
    for name in pair_names:
        s1, s2, _q1, _q2 = seqs[name]
        k1, k2 = kmers(s1), kmers(s2)
        overlap[name] = bool(k1 & k2)
    ov = np.asarray([overlap[n] for n in owner])
    result["overlap_rate_of_pairs"] = round(sum(overlap.values()) / len(overlap), 4)
    result["mates_overlap"] = recovery(Z, partner, args.ks, rng, subset=ov)
    result["mates_disjoint"] = recovery(Z, partner, args.ks, rng, subset=~ov)

    # ---- CONFOUND 2: the two directions are different measurements
    for m in (1, 2):
        result[f"from_mate{m}"] = recovery(Z, partner, args.ks, rng,
                                           subset=(which_arr == m))

    # ---- CONFOUND 3: abundance does most of the work already
    cov = X[:, AXES.index("log_cov")]
    dec = np.quantile(cov, np.linspace(0, 1, 11))
    by_decile = []
    for i in range(10):
        lo, hi = dec[i], dec[i + 1]
        sub = (cov >= lo) & (cov <= hi if i == 9 else cov < hi)
        if sub.sum() < 50:
            by_decile.append({"decile": i + 1, "n": int(sub.sum()), "skipped": True})
            continue
        r = recovery(Z, partner, [10], rng, subset=sub)
        by_decile.append({"decile": i + 1, "n": int(sub.sum()),
                          "log_cov_range": [round(float(lo), 3), round(float(hi), 3)],
                          "recovery_at_10": r["10"]["recovery"],
                          "enrichment_at_10": r["10"]["enrichment"]})
    result["by_coverage_decile"] = by_decile

    # Axis groups: how much of the recovery is composition alone, which is the half that is
    # not confounded by the cascade having removed what assembled.
    for gname, cols in (("coverage_only", COVERAGE_AXES),
                        ("composition_only", COMPOSITION_AXES)):
        sel = [AXES.index(a) for a in cols]
        result[gname] = recovery(Z[:, sel], partner, args.ks, rng)

    # ---- the null: a label permutation. This MUST land on chance, and if it does not the
    # geometry is leaking partner identity through something other than the axes.
    shuffled = rng.permutation(len(X))
    result["permuted_partner_null"] = recovery(Z, shuffled, args.ks, rng)

    # ---- DOES IT SURVIVE AT SCALE? Recovery at fixed k is a function of how crowded the
    # space is, so a figure from 40,000 points says nothing on its own about 35 million.
    # Subsampling whole PAIRS (never splitting one) gives the scaling for free, and it is the
    # scaling rather than the level that decides whether mate linkage is usable in a
    # production run: if recovery falls as 1/n the neighbourhood is a fixed-size lottery,
    # and if enrichment holds the axes genuinely locate the molecule.
    npairs = len(X) // 2
    ladder = []
    for frac in (0.125, 0.25, 0.5, 1.0):
        m = int(npairs * frac)
        if m < 100:
            continue
        keep = rng.permutation(npairs)[:m]
        idx = np.sort(np.concatenate([keep * 2, keep * 2 + 1]))
        # Standardised on the subsample, because that is what a run of that size would do.
        Zs = standardise(X[idx])
        ps = np.arange(len(idx)) ^ 1
        ovs = ov[idx]
        r = recovery(Zs, ps, [10], rng)
        rd = recovery(Zs, ps, [10], rng, subset=~ovs)
        ladder.append({"pairs": m, "points": len(idx),
                       "recovery_at_10": r["10"]["recovery"],
                       "enrichment_at_10": r["10"]["enrichment"],
                       "disjoint_recovery_at_10": rd["10"]["recovery"],
                       "disjoint_enrichment_at_10": rd["10"]["enrichment"]})
    result["scale_ladder"] = ladder

    # ---- CONFOUND 2, the control rather than the split: mate 1 truncated to mate 2's
    # length. The two mates differ by 14 nt and therefore by 14 scorable 31-mers, so any
    # asymmetry above could be length alone. This arm pays a second whole-library pass to
    # find out, which is the only honest way to ask it.
    if args.length_match:
        matched = [(n, s1[:len(s2)], s2, q1[:len(s2)], q2)
                   for n, s1, s2, q1, q2 in rows]
        Xm, owner_m, which_m, dropped_m, sec_m = build_space(args, matched, "lenmatch")
        Zm = standardise(Xm)
        pm = np.arange(len(Xm)) ^ 1
        wm = np.asarray(which_m)
        lm = {"all": recovery(Zm, pm, args.ks, rng),
              "pairs_dropped_no_kmers": dropped_m,
              "coverage_seconds": round(sec_m, 1),
              "mate1_truncated_to": "len(mate2)"}
        for m in (1, 2):
            lm[f"from_mate{m}"] = recovery(Zm, pm, args.ks, rng, subset=(wm == m))
        result["length_matched"] = lm

    out = {"script": "experiments/mate_friends.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "axes": AXES, "mate_points": int(len(X)), "pairs_used": int(len(X) // 2),
           "pairs_dropped_no_kmers": dropped, "file_offset_pairs": skip,
           "coverage_seconds": round(cover_seconds, 1),
           "result": result,
           "readout": "fraction of mates whose own partner is among their k nearest "
                      "neighbours in the standardised 10-axis per-mate space, against the "
                      "chance rate k/(n-1)",
           "bias": _bias(rows, KMER_K),
           "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    p = args.out / f"mate_friends_s{args.seed}.json"
    p.write_text(json.dumps(out, indent=2))
    r = result["all"]["10"] if "10" in result["all"] else {}
    print(f"seed {args.seed}: {len(X)} mate points; recovery@10 "
          f"{r.get('recovery')} against chance {r.get('chance')} "
          f"({r.get('enrichment')}x); disjoint-mate "
          f"{result['mates_disjoint']['10']['recovery']}; "
          f"null {result['permuted_partner_null']['10']['recovery']}", flush=True)
    print(f"written {p}", flush=True)


if __name__ == "__main__":
    main()
