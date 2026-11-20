"""Does the coverage-shape family lose its separating power as the subtrahend grows?

THE MECHANISM TEST. The shape axes (`cov_spread`, `cov_step`, `mate_disagree`, measured by
PR #17/#18 over the whole library with each pair's own occurrences subtracted) separate
residue from subtracted reads at 0.76-0.78 AUC, and `cov_step_conservation.py` showed what
they are measuring: witnesses, i.e. other molecules in the library carrying the same
segment. If that is right, then removing identifiable material should remove the witnessed
material first and the shape axes should LOSE separating power as the subtrahend grows.
That is a prediction about mechanism rather than another AUC.

THE CONFOUND, as the sibling thread stated it: a falling AUC is also what you would see if
a bigger index simply removed the high-coverage reads, because `log_cov` is the dominant
axis and is partly definitional. The discriminating version is therefore to watch the shape
axes WITHIN A FIXED COVERAGE BAND, because they are ratios rather than levels:

- if the shape axes hold their AUC inside a fixed `log_cov` band while the overall number
  falls, the fall was abundance;
- if they fall inside the band too, it is witnesses.

The band is fixed ONCE, from the full population at rung zero, and reused at every rung, so
that the comparison is not made against a moving definition.

THE LADDER is nucleotide-only and cumulative: SILVA, then the sample's own ribosomal
contigs, then its own identified non-ribosomal contigs, then the 1 kb contigs of all 116
samples. Nucleotide rather than protein because this has to run over all ~100,000 pairs of
the block rather than a subsample, and because the ordering (ribosomal first, with a
nucleotide index) is a requirement of the cascade it mirrors. The subtrahend size in bytes
is reported beside each rung so "as the subtrahend grows" is a measured axis and not a
label.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import platform
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from cov_step_join import rank
from subtract_known import _mm_index, git_sha, map_filter

ROOT = Path(__file__).resolve().parents[1]


def auc(pos: list[float], neg: list[float]) -> float | None:
    """Mann-Whitney AUC from ranks, ties shared. pos = the population the axis should score
    high (here: NOT identified at this rung)."""
    if not pos or not neg:
        return None
    all_v = pos + neg
    r = rank(all_v)
    s = sum(r[:len(pos)])
    n1, n2 = len(pos), len(neg)
    return (s - n1 * (n1 + 1) / 2) / (n1 * n2)


def auc_ci(a: float | None, n1: int, n2: int):
    """Hanley-McNeil interval. Reported because the later rungs leave a small unidentified
    arm inside the band - 237 pairs at the last rung - and an AUC near 0.5 from a few
    hundred observations is not the same statement as one from twenty thousand. This repo
    has already been caught once by a fit resting on one unbalanced arm."""
    if a is None or n1 < 2 or n2 < 2:
        return [None, None]
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2)
    se = math.sqrt(max(var, 0.0))
    return [round(max(0.0, a - 1.96 * se), 4), round(min(1.0, a + 1.96 * se), 4)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--axes", type=Path,
                    default=ROOT / "results/residue_structure/coverage_axes_allpairs_s1.tsv.gz")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--band", type=float, nargs=2, default=[0.4, 0.6],
                    help="quantiles of log_cov, from rung zero, defining the fixed band")
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/shapeladder"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/shape_axis_ladder")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    axes: dict[str, dict[str, float]] = {}
    with gzip.open(args.axes, "rt") as h:
        header = h.readline().rstrip("\n").split("\t")
        for line in h:
            f = line.rstrip("\n").split("\t")
            axes[f[0]] = {k: float(v) for k, v in zip(header[1:], f[1:])}

    rows = list(read_pairs(args.r1, args.r2, args.pairs, 14, (args.seed - 1) * 1_500_000))
    seqs = {r[0]: (r[1], r[2]) for r in rows if r[0] in axes}
    print(f"{len(seqs)} of {len(rows)} pairs carry coverage axes", flush=True)

    rungs = [
        ("silva", ROOT / "data/silva/silva_dna.fa"),
        ("own_rrna_contigs", ROOT / "data/geo/sample_contigs_rrna.fa"),
        ("own_identified_contigs", ROOT / "data/geo/sample_contigs_identified_norrna.fa"),
        ("all_sample_contigs_1kb", ROOT / "data/geo/all_contigs_1kb.fa"),
    ]

    fields = ["cov_step", "cov_spread", "mate_disagree", "zero_frac", "log_cov"]
    lc = sorted(axes[k]["log_cov"] for k in seqs)
    band = (lc[int(len(lc) * args.band[0])], lc[int(len(lc) * args.band[1])])

    def measure(identified: set[str], label: str, cum_bytes: int):
        pos = [k for k in seqs if k not in identified]
        neg = [k for k in seqs if k in identified]
        in_band = [k for k in seqs if band[0] <= axes[k]["log_cov"] <= band[1]]
        bpos = [k for k in in_band if k not in identified]
        bneg = [k for k in in_band if k in identified]
        row = {"rung": label, "cumulative_subtrahend_bytes": cum_bytes,
               "identified": len(neg), "not_identified": len(pos),
               "band_identified": len(bneg), "band_not_identified": len(bpos),
               "auc_overall": {}, "auc_within_band": {}}
        row["auc_overall_ci95"], row["auc_within_band_ci95"] = {}, {}
        for f in fields:
            a = auc([axes[k][f] for k in pos], [axes[k][f] for k in neg])
            b = auc([axes[k][f] for k in bpos], [axes[k][f] for k in bneg])
            row["auc_overall"][f] = a
            row["auc_within_band"][f] = b
            row["auc_overall_ci95"][f] = auc_ci(a, len(pos), len(neg))
            row["auc_within_band_ci95"][f] = auc_ci(b, len(bpos), len(bneg))
        return row

    out_rungs = [measure(set(), "nothing_subtracted", 0)]
    identified: set[str] = set()
    cum = 0
    for name, fa in rungs:
        idx = _mm_index(fa, args.work / f"{name}.mmi", args.threads)
        hit = map_filter(idx, seqs, args.min_mlen)
        identified |= set(hit)
        cum += fa.stat().st_size
        out_rungs.append(measure(identified, name, cum))
        r = out_rungs[-1]
        print(f"{name:26} cum {cum / 1e6:7.1f} MB  removed {r['identified']:6}  "
              f"cov_step AUC {r['auc_overall']['cov_step']:.3f} overall / "
              f"{r['auc_within_band']['cov_step']:.3f} in band  "
              f"(log_cov {r['auc_overall']['log_cov']:.3f} / "
              f"{r['auc_within_band']['log_cov']:.3f})", flush=True)

    first, last = out_rungs[1], out_rungs[-1]
    verdict = {}
    for f in fields:
        d_all = (last["auc_overall"][f] or 0) - (first["auc_overall"][f] or 0)
        d_band = (last["auc_within_band"][f] or 0) - (first["auc_within_band"][f] or 0)
        verdict[f] = {"delta_auc_overall": round(d_all, 4),
                      "delta_auc_within_band": round(d_band, 4)}

    out = {
        "seed": args.seed, "pairs_with_axes": len(seqs),
        "log_cov_band": [round(band[0], 4), round(band[1], 4)],
        "band_quantiles": args.band,
        "axes_source": str(args.axes),
        "rungs": out_rungs,
        "change_from_first_to_last_rung": verdict,
        "reading_rule": [
            "witnesses: the shape axes' AUC falls INSIDE the fixed band as well as overall",
            "abundance: the shape axes hold their AUC inside the band while the overall "
            "number falls, because log_cov is the dominant axis and partly definitional",
            "the band is fixed once from the rung-zero population and reused, so the "
            "comparison is not against a moving definition",
            "an AUC below 0.5 means the axis scores the IDENTIFIED pairs higher, which is "
            "the familiarity polarity rather than a weak result",
            "the unidentified arm inside the band shrinks as the subtrahend grows, to a "
            "few hundred pairs at the last rung, so read the intervals: the rung to trust "
            "for the in-band claim is the last one whose smaller arm is still in the "
            "thousands",
        ],
        "wall_seconds": time.perf_counter() - t0,
    }
    (args.out / f"shape_ladder_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/shape_axis_ladder.py", "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    print("\nchange, first rung to last:")
    for f, v in verdict.items():
        print(f"  {f:15} overall {v['delta_auc_overall']:+.3f}   "
              f"within band {v['delta_auc_within_band']:+.3f}")
    print(f"band: log_cov {band[0]:.3f}-{band[1]:.3f} "
          f"(median {st.median(lc):.3f})")


if __name__ == "__main__":
    main()
