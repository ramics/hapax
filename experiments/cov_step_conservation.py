"""Is within-read coverage step a conservation detector rather than a junction detector?

THE HYPOTHESIS, handed over by the thread measuring whether the residue is structured, and
the one explanation of its `cov_step` field still standing after two kills. The field is a
genuine change of exact 31-mer depth inside one 107 nt read (0.804 rank correlation with the
step recomputed after dropping unwitnessed k-mers, so it is not mostly sequencing error);
it is not a chimaeric junction (the contig anatomy fires on 0% of both unidentifiable and
identified contigs); and it is not an abundance boundary visible at contig level
(`cov_step_join.py`: +0.011 against the mates' depth ratio, interval across zero). What is
left is sequence SHARED with other molecules - a conserved or repeated segment having more
witnesses in the library than its neighbours - which would make the field a conservation
detector and much less interesting.

This repo is placed to test that, because the rRNA work already measured the mechanism it
would need: ribosomal sequence is reached 388 times per occurrence against 1.00 for a
shuffled control, so conserved material genuinely has more witnesses.

TWO MEASUREMENTS, both continuous and both unthresholded:

1. RIBOSOMAL SIMILARITY, as matching bases against SILVA with no cutoff applied. The
   cascade already removed every read matching SILVA over 40 bases, so what is left in the
   residue is weak similarity - exactly where a conserved segment shorter than the
   subtraction threshold would hide. If `cov_step` is conservation, it should rise with
   this.
2. SHARING WITHIN THE SAMPLE: how many DISTINCT contigs of the sample's own assembly a
   read aligns to, and the best alignment to a contig other than its own best. A segment
   repeated across molecules is the direct form of the hypothesis and does not depend on
   the segment being ribosomal or being in any database.

THE CONFOUND THAT HAS TO BE HELD. `cov_step` correlates with `log_cov`, and an abundant
read hits more contigs and more of SILVA for reasons that have nothing to do with a step.
So every correlation is also reported WITHIN quartiles of `log_cov`: if it survives inside
the strata it is about the step, and if it only exists pooled it is about abundance. The
shuffled column is the null for each cell.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cov_step_join import fisher_ci, read_fasta_mates, spearman
from subtract_known import git_sha

ROOT = Path(__file__).resolve().parents[1]
COV = re.compile(r"_cov_([0-9.]+)")


def contig_profile(idx, seq: str, min_mlen: int):
    """Best matching bases per distinct contig, so sharing is counted over MOLECULES rather
    than over alignments: one contig returning five overlapping alignments is one witness."""
    best: dict[str, int] = {}
    for h in idx.map(seq):
        if h.mlen < min_mlen:
            continue
        if h.mlen > best.get(h.ctg, 0):
            best[h.ctg] = h.mlen
    return best


def best_mlen(idx, seq: str) -> int:
    return max((h.mlen for h in idx.map(seq)), default=0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--axes", type=Path,
                    default=ROOT / "results/residue_structure/coverage_axes_s1.tsv")
    ap.add_argument("--reads", type=Path,
                    default=Path("/tmp/claude-0/wk/recut/residue_s1.fa"))
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--silva", type=Path, default=ROOT / "data/silva/silva_dna.fa")
    ap.add_argument("--share-min-mlen", type=int, default=30,
                    help="shorter than the cascade's 40, because the segment being tested "
                         "for is by hypothesis short enough to have survived subtraction")
    ap.add_argument("--step-high", type=float, default=2.0,
                    help="the sibling thread's own cut for a high-step read")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/covcons"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/cov_step_conservation")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    axes: dict[str, dict[str, float]] = {}
    with open(args.axes) as h:
        header = h.readline().rstrip("\n").split("\t")
        for line in h:
            f = line.rstrip("\n").split("\t")
            axes[f[0]] = {k: float(v) for k, v in zip(header[1:], f[1:])}
    mates = read_fasta_mates(args.reads)

    import mappy
    ctg_mmi, sil_mmi = args.work / "contigs.mmi", args.work / "silva.mmi"
    ctg = mappy.Aligner(str(args.contigs), preset="sr", n_threads=args.threads,
                        fn_idx_out=str(ctg_mmi)) if not ctg_mmi.exists() else \
        mappy.Aligner(str(ctg_mmi), preset="sr", n_threads=args.threads)
    sil = mappy.Aligner(str(args.silva), preset="sr", n_threads=args.threads,
                        fn_idx_out=str(sil_mmi)) if not sil_mmi.exists() else \
        mappy.Aligner(str(sil_mmi), preset="sr", n_threads=args.threads)
    if not ctg or not sil:
        raise SystemExit("minimap2 index failed")

    rows = []
    for name, d in mates.items():
        a = axes.get(name)
        if a is None:
            continue
        seqs = [s for s in (d.get("1"), d.get("2")) if s]
        if not seqs:
            continue
        prof: dict[str, int] = {}
        for s in seqs:
            for c, m in contig_profile(ctg, s, args.share_min_mlen).items():
                prof[c] = max(prof.get(c, 0), m)
        order = sorted(prof.values(), reverse=True)
        rows.append({
            "name": name,
            "silva_mlen": max(best_mlen(sil, s) for s in seqs),
            "distinct_contigs": len(prof),
            "second_contig_mlen": order[1] if len(order) > 1 else 0,
            **a,
        })

    rng = random.Random(args.seed)
    fields = ["cov_step", "cov_spread", "mate_disagree", "zero_frac", "log_cov"]
    targets = ["silva_mlen", "distinct_contigs", "second_contig_mlen"]

    def corr(sub, field, target):
        if len(sub) < 30:
            return None
        xs = [r[field] for r in sub]
        ys = [r[target] for r in sub]
        rho = spearman(xs, ys)
        perm = ys[:]
        rng.shuffle(perm)
        return {"n": len(sub), "spearman": round(rho, 4), "ci95": fisher_ci(rho, len(sub)),
                "shuffled": round(spearman(xs, perm), 4)}

    # strata on log_cov, so a correlation that is really abundance cannot hide in the pool
    lc = sorted(r["log_cov"] for r in rows)
    qs = [lc[int(len(lc) * f)] for f in (0.25, 0.5, 0.75)] if lc else []

    def stratum(i):
        if not qs:
            return []
        if i == 0:
            return [r for r in rows if r["log_cov"] <= qs[0]]
        if i == 3:
            return [r for r in rows if r["log_cov"] > qs[2]]
        return [r for r in rows if qs[i - 1] < r["log_cov"] <= qs[i]]

    hi = [r for r in rows if r["cov_step"] >= args.step_high]
    lo = [r for r in rows if r["cov_step"] < args.step_high]

    def group(sub, key):
        vals = [r[key] for r in sub]
        return {"n": len(sub), "mean": round(st.mean(vals), 4) if vals else None,
                "median": round(st.median(vals), 4) if vals else None,
                "frac_nonzero": round(sum(1 for v in vals if v > 0) / len(vals), 4)
                if vals else None}

    out = {
        "seed": args.seed,
        "residue_pairs": len(rows),
        "high_step_pairs": len(hi),
        "pooled": {t: {f: corr(rows, f, t) for f in fields} for t in targets},
        "within_log_cov_quartiles": {
            t: {f: [corr(stratum(i), f, t) for i in range(4)] for f in ["cov_step"]}
            for t in targets},
        "high_vs_low_step": {t: {"high": group(hi, t), "low": group(lo, t)}
                             for t in targets},
        "interpretation_rules": [
            "the hypothesis passes only if cov_step rises with sharing INSIDE the log_cov "
            "quartiles; a correlation that exists pooled and vanishes within strata is "
            "abundance, which cov_step is known to track",
            "silva_mlen is weak ribosomal similarity by construction: every read matching "
            "SILVA over 40 bases was already removed to form the residue, so this field "
            "can only see conserved segments shorter than that",
            "distinct_contigs counts molecules rather than alignments, and is still a "
            "property of the assembly: a shared segment whose other copies never assembled "
            "is invisible here, which biases the test towards finding nothing",
        ],
        "wall_seconds": time.perf_counter() - t0,
    }
    (args.out / f"conservation_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/cov_step_conservation.py", "git_sha": git_sha(),
        "axes_source": str(args.axes),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    print(f"{len(rows)} residue pairs, {len(hi)} with cov_step >= {args.step_high}")
    for t in targets:
        print(f"== {t}")
        for f in fields:
            c = out["pooled"][t][f]
            if c:
                print(f"   pooled {f:14} rho {c['spearman']:+.3f} ci {c['ci95']} "
                      f"shuffled {c['shuffled']:+.3f}")
        q = out["within_log_cov_quartiles"][t]["cov_step"]
        print("   cov_step within log_cov quartiles: " + ", ".join(
            f"{c['spearman']:+.3f}" if c else "n/a" for c in q))
        g = out["high_vs_low_step"][t]
        print(f"   high-step mean {g['high']['mean']} vs low-step {g['low']['mean']} "
              f"(nonzero {g['high']['frac_nonzero']} vs {g['low']['frac_nonzero']})")


if __name__ == "__main__":
    main()
