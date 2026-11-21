#!/usr/bin/env python3
"""What deeper sequencing of this library would buy, priced from the library we have.

WHY THIS EXISTS
---------------
`results/residue_census/` §3 put the unplaced residue at median 31-mer depth 1, with 91%
of it occurring fewer than five times in 8.7 M pairs and about a third occurring nowhere
else. Everything that needs a second observation of a molecule — assembly, neighbourhood
grouping, mate-linkage kinship, within-run recurrence — therefore has nothing to work
with on that pile, and the census's own conclusion was that this is sampling rather than
tuning: *you cannot assemble sequence you hold one copy of*.

That makes depth the governing variable, and depth is the one variable money can change.
So the question "should we sequence this water deeper" is the most expensive open question
in the project, and it does not have to be answered by spending. It can be answered by
rarefaction: measure how the depth distribution of a FIXED set of molecules moves as the
library around them grows, and read the slope.

This is a decision instrument, not a novelty field. Its output is a curve with an
operating point, which is what SPEC.md §8.3 requires of anything shown.

DESIGN
------
A fixed query block of real pairs is scored against nested subsample libraries at
f = 1/8, 1/4, 1/2, 1 of the run. Nested, so the curve is a rarefaction of one library
rather than four unrelated samples: a pair included at f is included at every larger f.
The query pairs are EXCLUDED from every library at every level, so a depth is always
"occurrences elsewhere" and no level is contaminated by a query witnessing itself. That
differs from the census, which scored reads that were in the library and subtracted their
own occurrences; stated because it makes this instrument's depths slightly lower than the
census's at f = 1 by construction.

Five seeds: each draws its own query block AND its own nested subsample, so the spread
covers both.

WHAT CAN AND CANNOT BE READ OFF IT
----------------------------------
Readable: the slope. Whether the fully-unwitnessed fraction is falling fast or has
flattened between f = 1/8 and f = 1 is a measured fact about this library.

Extrapolation is where this kind of curve goes wrong, so the horizon is stated rather than
implied. Interpolation inside [1/8, 1] is measurement. A doubling is a short, fairly safe
extrapolation. Ten times is NOT: the accumulation estimators that are reliable to about
2x are not reliable to 10x, and this script therefore reports the measured points and the
doubling, and refuses to print a 10x number. The repo has already withdrawn one claim for
extrapolating a curve past where it was measured (`results/absence_curve/`), which is why.

BIAS, with direction
--------------------
Subsampling a finished library is not the same experiment as sequencing deeper. Deeper
sequencing of the SAME library re-samples the same molecules, which is what this models.
Sequencing new water samples a different population — and the same plant's rRNA fraction
went 38.0% to 67.1% in eight months (`results/dna_rna_pair/`), so a figure from new water
would be a figure about a different sample. Direction: this curve is optimistic about new
water and roughly right about re-sequencing the same library deeper.

Second bias: a subsample of a library with PCR duplicates in it will show duplicates as
depth, and this library's preparation is not duplicate-free. Direction: inflates depth at
every level, so the fully-unwitnessed fraction is a lower bound throughout.

USAGE
-----
  python3 experiments/depth_rarefaction.py --seeds 1,2,3,4,5 --out results/depth_rarefaction
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
LEVELS = (0.125, 0.25, 0.5, 1.0)
TRIM_R2 = 14
LIB_PAIRS = 8_731_224          # SRR23973193, from data/MANIFEST.md


def write_subsamples(r1: Path, r2: Path, outdir: Path, query_idx: set[int], seed: int,
                     log) -> tuple[dict[float, tuple[Path, Path]], Path, int]:
    """One streaming pass writing every nested level plus the query block.

    compresslevel=1 throughout: these files are an intermediate, and the time saved is
    the difference between this being a cheap experiment and an annoying one.
    """
    rng = random.Random(1_000_000 + seed)
    handles = {}
    for f in LEVELS:
        p1 = outdir / f"lib_f{f}_1.fastq.gz"
        p2 = outdir / f"lib_f{f}_2.fastq.gz"
        handles[f] = (gzip.open(p1, "wt", compresslevel=1),
                      gzip.open(p2, "wt", compresslevel=1), p1, p2)
    qpath = outdir / "query.tsv"
    qfh = open(qpath, "w")

    n = 0
    nq = 0
    t0 = time.time()
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline(); p1l = f1.readline(); q1 = f1.readline()
            h2 = f2.readline(); s2 = f2.readline(); p2l = f2.readline(); q2 = f2.readline()
            if n in query_idx:
                a = s1.strip()
                b = s2.strip()[TRIM_R2:]
                if "N" not in a and "N" not in b:
                    qfh.write(f"q{n}\t{a}\t{b}\n")
                    nq += 1
                n += 1
                continue
            u = rng.random()
            for f, (o1, o2, _, _) in handles.items():
                if u < f:
                    o1.write(h1); o1.write(s1); o1.write(p1l); o1.write(q1)
                    o2.write(h2); o2.write(s2); o2.write(p2l); o2.write(q2)
            n += 1
            if n % 2_000_000 == 0:
                log(f"  seed {seed}: {n:,} pairs in {time.time() - t0:.0f}s")
    paths = {}
    for f, (o1, o2, p1, p2) in handles.items():
        o1.close(); o2.close()
        paths[f] = (p1, p2)
    qfh.close()
    log(f"  seed {seed}: {n:,} pairs read, {nq:,} query pairs kept, "
        f"{time.time() - t0:.0f}s")
    return paths, qpath, n


def score_level(binary: Path, lib: tuple[Path, Path], query: Path, out: Path, log) -> dict:
    """Run the exact-31-mer instrument with the query block as QUERY_ONLY pairs.

    QUERY_ONLY (block 9 in that tool) means scored against the library but never added to
    it, which is the convention this experiment needs.
    """
    t0 = time.time()
    with open(out, "w") as fh:
        subprocess.run([str(binary), str(lib[0]), str(lib[1]), str(TRIM_R2),
                        "1", "0", "-", str(query)], check=True, stdout=fh)
    elapsed = time.time() - t0
    log(f"    scored in {elapsed:.0f}s -> {out.name}")
    return {"seconds": elapsed}


def summarise_level(path: Path) -> dict:
    """Per-mate depth summary for one (seed, level) from the instrument's TSV.

    Columns are taken from the file's own header rather than by position, because a column
    order assumed here and changed there is the silent-join failure this repo has already
    published once (`results/partial_eyeball/`).

    Only block 9 rows are kept: those are the QUERY_ONLY pairs, scored against the
    library but never added to it. The tool reports raw counts for block 9 and
    own-subtracted counts for every other block, so mixing blocks would mix two
    definitions of depth.

    Reported PER MATE, never pooled: mate 2 is 93 nt against mate 1's 107, which is 62
    scorable 31-mers against 76, and a pooled figure averages two different instruments
    (CLAUDE.md).
    """
    import csv
    med = {"1": [], "2": []}
    nk = {"1": [], "2": []}
    libpairs = None
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            if r["block"] != "9":
                continue
            libpairs = int(r["library_pairs"])
            for m in ("1", "2"):
                if int(r[f"n{m}"]) > 0:
                    med[m].append(int(r[f"med{m}"]))
                    nk[m].append(int(r[f"n{m}"]))
    out = {"library_pairs": libpairs}
    for m in ("1", "2"):
        v = np.array(med[m], dtype=float)
        if v.size == 0:
            out[f"mate{m}"] = {"n": 0}
            continue
        out[f"mate{m}"] = {
            "n": int(v.size),
            "median_depth": float(np.median(v)),
            "mean_depth": float(v.mean()),
            "frac_depth0": float((v == 0).mean()),
            "frac_depth_lt5": float((v < 5).mean()),
            "frac_depth_lt2": float((v < 2).mean()),
            "mean_scorable_kmers": float(np.mean(nk[m])),
        }
    return out


def curve_report(per_seed: dict, log) -> dict:
    """The curve, its slope, and a doubling — and deliberately no 10x number.

    The reported extrapolation is one doubling, by linear fit of the unwitnessed fraction
    against log2(effort) over the measured points. The fit's own residual is reported
    beside it: a curve that is not straight in log2 effort has no business being
    extrapolated at all, and saying so is cheaper than withdrawing it later.
    """
    rep = {"levels": list(LEVELS), "per_seed": per_seed, "mates": {}}
    for m in ("1", "2"):
        xs = np.log2(np.array(LEVELS, dtype=float))
        series = []
        for seed, lv in per_seed.items():
            ys = [lv[str(f)][f"mate{m}"]["frac_depth0"] for f in LEVELS]
            series.append(ys)
        arr = np.array(series)                      # seeds x levels
        mean = arr.mean(axis=0)
        slope, intercept = np.polyfit(xs, mean, 1)
        resid = float(np.max(np.abs(np.polyval([slope, intercept], xs) - mean)))
        rep["mates"][m] = {
            "frac_depth0_mean": [float(v) for v in mean],
            "frac_depth0_min": [float(v) for v in arr.min(axis=0)],
            "frac_depth0_max": [float(v) for v in arr.max(axis=0)],
            "slope_per_doubling": float(slope),
            "max_abs_residual_of_linear_fit": resid,
            "projected_one_doubling": float(mean[-1] + slope),
            "extrapolation_horizon": (
                "one doubling only; a 10x figure is not reported because the"
                " accumulation estimators behind it are not reliable that far"),
        }
        log(f"mate {m}: unwitnessed fraction {[round(float(v), 4) for v in mean]}"
            f" over f={list(LEVELS)}, slope {slope:+.4f}/doubling, fit residual {resid:.4f}")
    return rep


def git_sha() -> str:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, check=True,
                             capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO, check=True,
                               capture_output=True, text=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--out", default="results/depth_rarefaction")
    ap.add_argument("--query-pairs", type=int, default=10_000)
    ap.add_argument("--work", default=None)
    ap.add_argument("--keep-libs", action="store_true",
                    help="keep the intermediate subsample fastqs (they are large)")
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    outdir = REPO / args.out
    outdir.mkdir(parents=True, exist_ok=True)
    logfh = open(outdir / "run.log", "a")

    def log(msg: str) -> None:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        logfh.write(line + "\n")
        logfh.flush()

    r1 = REPO / "data" / "berlin" / "SRR23973193_1.fastq.gz"
    r2 = REPO / "data" / "berlin" / "SRR23973193_2.fastq.gz"
    binary = REPO / "experiments" / "kmer_coverage" / "target" / "release" / "kmer_coverage"
    for p in (r1, r2, binary):
        if not p.exists():
            log(f"MISSING INPUT: {p}")
            return 2

    work = Path(args.work) if args.work else outdir / "work"
    work.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    timings = {}
    per_seed: dict[int, dict] = {}
    for seed in seeds:
        sdir = work / f"s{seed}"
        sdir.mkdir(exist_ok=True)
        qidx = set(random.Random(seed).sample(range(LIB_PAIRS), args.query_pairs))
        log(f"seed {seed}: writing nested subsamples")
        libs, qpath, n_seen = write_subsamples(r1, r2, sdir, qidx, seed, log)
        if n_seen != LIB_PAIRS:
            log(f"  NOTE: library holds {n_seen:,} pairs, manifest says {LIB_PAIRS:,}")
        per_seed[seed] = {}
        for f in LEVELS:
            out = sdir / f"depth_f{f}.tsv"
            log(f"  seed {seed} level f={f}")
            timings[f"s{seed}_f{f}"] = score_level(binary, libs[f], qpath, out, log)
            per_seed[seed][str(f)] = summarise_level(out)
            shutil.copy(out, outdir / f"depth_s{seed}_f{f}.tsv")
        if not args.keep_libs:
            for f in LEVELS:
                for p in libs[f]:
                    p.unlink(missing_ok=True)
            log(f"  seed {seed}: intermediate libraries removed")

    report = curve_report(per_seed, log)
    with open(outdir / "curve.json", "w") as fh:
        json.dump(report, fh, indent=1)

    manifest = {
        "script": "experiments/depth_rarefaction.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "git_sha": git_sha(),
        "seeds": seeds,
        "params": {"levels": list(LEVELS), "query_pairs": args.query_pairs,
                   "trim_r2": TRIM_R2, "k": 31,
                   "query_excluded_from_library": True,
                   "library_pairs_expected": LIB_PAIRS},
        "timings_s": timings,
        "total_s": time.time() - t_start,
    }
    with open(outdir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    log(f"done in {time.time() - t_start:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
