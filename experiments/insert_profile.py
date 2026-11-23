"""What insert sizes a real wastewater RNA library actually delivers.

Why this exists. The wide-insert hypothesis under test in this thread is that a
library template with the two mates far apart maximises the chance that a pair
hits a coding region and the chance that it sees something novel. Before any of
that is modelled, one thing has to be measured rather than assumed: what insert
distribution the sample material can even produce. Wastewater RNA is already
fragmented, so the achievable insert may be capped by the molecules in the
sample rather than by the protocol, and in that case the hypothesis is not a
design choice at all. The project has twice quoted ~250 nt for SRR23973193's
inserts; that number was an ART simulation parameter in the abundance sweep, not
a measurement of this library (recorded as an open caveat in
results/readlevel_features/REPORT.md §6).

Two independent estimators, because each is biased in a known direction and the
pair of them bounds the truth:

1. **Mate overlap** (reference-free). If the insert is shorter than the sum of
   the two read lengths, the mates overlap and the overlap offset gives the
   insert exactly. This sees only inserts below ~200 nt, so on its own it
   measures the left tail and the mergeable fraction - which is also the number
   the dropped merged-mate mechanism turned on (NAO reports 58-76% of wastewater
   viral pairs merge at 2x150).
2. **Reference mapping** (mappy/minimap2 against the study authors' own SPAdes
   contigs for this sample). This sees the whole distribution but only for pairs
   whose template assembled, which biases towards abundant organisms.

Mate 2 of this library carries ~14 nt of non-template 5' sequence (measured in
results/ww_berlin_consensus/REPORT.md: 2,455 of 2,572 reads soft-clipped, mode
14). It is trimmed here before anything else; leaving it in inflates every
insert estimate by 14 nt and destroys the overlap estimator outright.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import statistics
import subprocess
import time
from collections import Counter
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

import library_facts

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_pairs(f1: Path, f2: Path, limit: int, trim_r2: int, skip: int = 0):
    """Stream (name, r1, r2) with mate 2's non-template prefix removed.

    `skip` discards whole pairs before yielding, which is how the seeds differ:
    a run is laid out by flow-cell tile, so distinct blocks of the file are
    distinct regions of the run rather than a reshuffle of the same reads.
    """
    with gzip.open(f1, "rt") as h1, gzip.open(f2, "rt") as h2:
        for i, (l1, l2) in enumerate(zip(h1, h2)):
            rec = i // 4
            if rec < skip:
                continue
            phase = i % 4
            if phase == 0:
                name = l1[1:].split()[0]
            elif phase == 1:
                yield name, l1.strip(), l2.strip()[trim_r2:]
                if rec + 1 - skip >= limit:
                    return


def overlap_insert(r1: str, r2rc: str, k: int = 16, min_ov: int = 20,
                   max_mm: float = 0.08) -> int | None:
    """Insert length from mate overlap, or None if the mates do not overlap.

    r2rc is the reverse complement of mate 2, so in template coordinates r1
    starts at 0 and r2rc ends at the template's 3' end. If r2rc starts at offset
    `o` inside r1 then the insert is o + len(r2rc). Candidate offsets come from
    shared k-mers (cheap) and are then verified over the whole overlap, keeping
    the longest verified overlap, which is the one consistent with a real
    template rather than a repeat.
    """
    pos: dict[str, int] = {}
    for i in range(len(r1) - k + 1):
        pos.setdefault(r1[i:i + k], i)
    cands = set()
    for j in range(0, min(len(r2rc) - k + 1, len(r2rc)), 4):
        p = pos.get(r2rc[j:j + k])
        if p is not None and p - j >= 0:
            cands.add(p - j)
    best = None
    for o in sorted(cands):
        ov = min(len(r1) - o, len(r2rc))
        if ov < min_ov:
            continue
        a, b = r1[o:o + ov], r2rc[:ov]
        mm = sum(1 for x, y in zip(a, b) if x != y)
        if mm <= max_mm * ov:
            ins = o + len(r2rc)
            if best is None or ov > best[1]:
                best = (ins, ov)
    return best[0] if best else None


def mapping_inserts(f1: Path, f2: Path, ref: Path, limit: int, trim_r2: int,
                    skip: int = 0) -> tuple[list[int], dict]:
    """Insert lengths for pairs whose template is in the authors' assembly."""
    import mappy

    aln = mappy.Aligner(str(ref), preset="sr", n_threads=4)
    if not aln:
        raise SystemExit(f"could not index {ref}")
    inserts: list[int] = []
    stats = Counter()
    for name, r1, r2 in read_pairs(f1, f2, limit, trim_r2, skip):
        stats["pairs"] += 1
        h1 = [h for h in aln.map(r1) if h.is_primary]
        if not h1:
            continue
        h2 = [h for h in aln.map(r2) if h.is_primary]
        if not h2:
            stats["one_mate_only"] += 1
            continue
        a, b = h1[0], h2[0]
        if a.ctg != b.ctg:
            stats["discordant_contig"] += 1
            continue
        if a.strand == b.strand:
            stats["same_strand"] += 1
            continue
        stats["concordant"] += 1
        span = max(a.r_en, b.r_en) - min(a.r_st, b.r_st)
        # clipped bases are template too: add back what the aligner left off
        span += a.q_st + (len(r1) - a.q_en) + b.q_st + (len(r2) - b.q_en)
        if 0 < span < 2000:
            inserts.append(span)
    return inserts, dict(stats)


def summarise(x: list[int]) -> dict:
    if not x:
        return {"n": 0}
    x = sorted(x)
    q = statistics.quantiles(x, n=100) if len(x) > 100 else [x[0]] * 99
    return {
        "n": len(x),
        "mean": round(statistics.fmean(x), 1),
        "sd": round(statistics.pstdev(x), 1) if len(x) > 1 else 0.0,
        "min": x[0],
        "p01": round(q[0]), "p05": round(q[4]), "p25": round(q[24]),
        "median": round(q[49]), "p75": round(q[74]), "p95": round(q[94]),
        "p99": round(q[98]), "max": x[-1],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    ap.add_argument("--ref", type=Path, default=None,
                    help="contigs for the mapping estimator (authors' SPAdes assembly)")
    ap.add_argument("--pairs", type=int, default=200_000)
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="non-template 5' prefix on mate 2 (measured: mode 14 nt)")
    ap.add_argument("--block", type=int, default=2_000_000,
                    help="pairs between seeds: seed s reads the block starting at (s-1)*block")
    ap.add_argument("--out", type=Path, default=Path("results/insert_profile"))
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["trim_r2"])

    random.seed(args.seed)
    skip = (args.seed - 1) * args.block
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    ov_inserts: list[int] = []
    lens1, lens2 = Counter(), Counter()
    n = 0
    for name, r1, r2 in read_pairs(args.r1, args.r2, args.pairs, args.trim_r2, skip):
        n += 1
        lens1[len(r1)] += 1
        lens2[len(r2)] += 1
        ins = overlap_insert(r1, revcomp(r2))
        if ins is not None:
            ov_inserts.append(ins)
    t_ov = time.perf_counter() - t0

    result = {
        "seed": args.seed,
        "pairs_examined": n,
        "read_len_mate1_mode": lens1.most_common(1)[0],
        "read_len_mate2_mode_after_trim": lens2.most_common(1)[0],
        "trim_r2": args.trim_r2,
        "overlap_estimator": {
            "pairs_with_detected_overlap": len(ov_inserts),
            "fraction_overlapping": round(len(ov_inserts) / n, 4) if n else None,
            "insert_summary": summarise(ov_inserts),
            "histogram_20nt_bins": dict(sorted(
                Counter(20 * (i // 20) for i in ov_inserts).items())),
            "wall_s": round(t_ov, 1),
        },
    }

    if args.ref:
        t1 = time.perf_counter()
        map_inserts, stats = mapping_inserts(args.r1, args.r2, args.ref,
                                             args.pairs, args.trim_r2, skip)
        result["mapping_estimator"] = {
            "reference": str(args.ref),
            "pair_disposition": stats,
            "insert_summary": summarise(map_inserts),
            "histogram_20nt_bins": dict(sorted(
                Counter(20 * (i // 20) for i in map_inserts).items())),
            "fraction_concordant_pairs_overlapping": round(
                sum(1 for i in map_inserts if i < 200) / len(map_inserts), 4)
            if map_inserts else None,
            # the cost of asking for width: what fraction of molecules survives a
            # minimum-insert requirement. This is the quantity the wide-insert
            # hypothesis has to pay, since at 1e-6 abundance molecules sampled is
            # the binding resource.
            "yield_at_min_insert": {
                str(t): round(sum(1 for i in map_inserts if i >= t) / len(map_inserts), 5)
                for t in (150, 200, 250, 300, 350, 400, 500, 600, 800)
            } if map_inserts else None,
            "wall_s": round(time.perf_counter() - t1, 1),
        }

    (args.out / f"insert_profile_s{args.seed}.json").write_text(json.dumps(result, indent=2))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/insert_profile.py",
        "git_sha": sha,
        "seed": args.seed,
        "params": {"pairs": args.pairs, "trim_r2": args.trim_r2, "skip_pairs": skip,
                   "r1": str(args.r1), "r2": str(args.r2), "ref": str(args.ref)},
        "wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2))
    print(json.dumps(result, indent=2)[:4000])


if __name__ == "__main__":
    main()
