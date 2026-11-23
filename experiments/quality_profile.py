"""Do these reads get grittier towards one end, and does the mate asymmetry bias a field
that treats the two mates as interchangeable?

WHY THIS IS BEING MEASURED AGAIN. "Illumina reads get worse towards the 3' end" is true of
Illumina in general and it is the assumption a reader brings to any position-independent
field. This library was measured once as FLAT, which was put down to NovaSeq quality
binning, and every field that trusts base quality rests on that single measurement. One
measurement of a load-bearing fact is not enough, so this re-measures it across five
flow-cell blocks, and - more importantly - measures it a second way that does not trust the
quality strings at all.

TWO INSTRUMENTS, because binned quality scores are a report of the instrument's opinion
and not an observation:

1. THE QUALITY PROFILE: mean Q, and the fraction of bases below Q20 and Q30, per position,
   for each mate separately. This is what the usual claim is about, and if NovaSeq binning
   has flattened it this will look flat for a reason that has nothing to do with the real
   error rate.
2. THE OVERLAP PROFILE, which needs no reference and no trust in Q. 53% of this library's
   pairs overlap their own mate (median insert 194 nt against 2x107 reads), so where the
   mates overlap each base is observed twice, and a mismatch between them means one of the
   two observations is wrong. Counting mismatches by POSITION IN EACH MATE gives a direct
   per-position error signal: a ramp in a mate shows up as rising mismatch rate at that
   mate's late positions, whatever the quality strings say. The rate is the SUM of the two
   mates' error rates at the paired positions, which is why the two mates' profiles are
   reported separately rather than halved.

   Its bias is stated rather than corrected: only overlapping pairs enter, so the
   measurement is of short molecules, and a mismatch cannot be attributed to one mate
   without using the quality scores this instrument exists to avoid trusting.

3. THE CONSEQUENCE FOR MATE-INTERCHANGEABLE FIELDS. The error-floor field calibrated in
   report section 5.2 counts pairs where exactly one mate is unwitnessed, and treats the
   two mates as equivalent. If R2 really is ~1.3x noisier, that field should fire
   asymmetrically. Measured directly here as which mate is the unmapped one when exactly
   one mate maps to the sample's own assembly - a count the field itself never looks at.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts

from berlin_io import read_pairs, revcomp
from subtract_known import _mm_index, provenance

ROOT = Path(__file__).resolve().parents[1]


def q_of(ch: str) -> int:
    return ord(ch) - 33


def best_overlap(s1: str, s2rc: str, min_ov: int, max_mismatch_frac: float):
    """Best ungapped overlap of mate 1 with the reverse complement of mate 2.

    Returns (offset, overlap length, mismatch positions as (i in s1, j in s2rc)). The
    offset search is exhaustive over a 107 nt read, which is cheap, and a best-of-all-
    offsets choice is required rather than using the insert distribution: using an
    estimated insert to place the overlap would make the mismatch count depend on the
    estimate.
    """
    best = None
    for off in range(-len(s2rc) + min_ov, len(s1) - min_ov + 1):
        i0 = max(0, off)
        j0 = i0 - off
        n = min(len(s1) - i0, len(s2rc) - j0)
        if n < min_ov:
            continue
        mm = [(i0 + k, j0 + k) for k in range(n) if s1[i0 + k] != s2rc[j0 + k]]
        if len(mm) > max_mismatch_frac * n:
            continue
        if best is None or len(mm) / n < len(best[2]) / best[1]:
            best = (off, n, mm)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--pairs", type=int, default=20_000)
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="mate 2's non-template 5' prefix. From the library unless given")
    ap.add_argument("--block", type=int, default=None,
                    help="pairs skipped per seed, so the seeds read distinct regions of\n                         " "the flow cell. From the library unless given")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--min-overlap", type=int, default=40)
    ap.add_argument("--max-mismatch-frac", type=float, default=0.1,
                    help="above this an 'overlap' is two different molecules, not an error")
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--mate-asymmetry", action="store_true", default=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/qprof"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/quality_profile")
    args = ap.parse_args()
    # A trim length and a block stride are facts about a sequencing run, and both
    # were literals inside the read_pairs call below, where no help text and no
    # manifest ever showed them. They now come from experiments/libraries.json by
    # accession, which refuses a library it does not know rather than applying
    # Berlin's 14 nt and 1.5 M stride to another library's reads.
    library_facts.resolve_args(args, ["trim_r2", "block"])

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    import mappy
    idx = _mm_index(args.contigs, args.work / "contigs.mmi", args.threads) \
        if args.mate_asymmetry else None

    per_seed = []
    for seed in args.seeds:
        L = 120
        qsum = {1: [0] * L, 2: [0] * L}
        qn = {1: [0] * L, 2: [0] * L}
        below20 = {1: [0] * L, 2: [0] * L}
        below30 = {1: [0] * L, 2: [0] * L}
        ovn = {1: [0] * L, 2: [0] * L}
        ovmm = {1: [0] * L, 2: [0] * L}
        pairs = overlapping = 0
        only1 = only2 = both = neither = 0
        only1_t = only2_t = both_t = neither_t = 0

        for name, s1, s2, q1, q2 in read_pairs(args.r1, args.r2, args.pairs, args.trim_r2,
                                               (seed - 1) * args.block):
            pairs += 1
            for mate, q in ((1, q1), (2, q2)):
                for i, ch in enumerate(q[:L]):
                    v = q_of(ch)
                    qsum[mate][i] += v
                    qn[mate][i] += 1
                    if v < 20:
                        below20[mate][i] += 1
                    if v < 30:
                        below30[mate][i] += 1
            ov = best_overlap(s1, revcomp(s2), args.min_overlap, args.max_mismatch_frac)
            if ov:
                overlapping += 1
                off, n, mm = ov
                i0 = max(0, off)
                j0 = i0 - off
                for k in range(n):
                    ovn[1][i0 + k] += 1
                    # position in mate 2's own 5'->3' order, since the overlap was taken
                    # against its reverse complement
                    ovn[2][len(s2) - 1 - (j0 + k)] += 1
                for i, j in mm:
                    ovmm[1][i] += 1
                    ovmm[2][len(s2) - 1 - j] += 1
            if idx is not None:
                m1 = any(h.mlen >= args.min_mlen for h in idx.map(s1))
                m2 = any(h.mlen >= args.min_mlen for h in idx.map(s2))
                if m1 and m2:
                    both += 1
                elif m1:
                    only1 += 1
                elif m2:
                    only2 += 1
                else:
                    neither += 1
                # LENGTH-MATCHED CONTROL, and it is not optional: mate 2 loses 14 nt to
                # the non-template prefix trim, so it is 93 nt against 107, and a shorter
                # read maps less often for reasons that have nothing to do with quality.
                # Without this the asymmetry below is unreadable.
                m1t = any(h.mlen >= args.min_mlen for h in idx.map(s1[:len(s2)]))
                if m1t and m2:
                    both_t += 1
                elif m1t:
                    only1_t += 1
                elif m2:
                    only2_t += 1
                else:
                    neither_t += 1

        def prof(num, den):
            return [round(num[i] / den[i], 6) if den[i] >= 50 else None for i in range(L)]

        def ramp(values, frac=0.25):
            """Mean over the last quarter of the read minus mean over the first quarter.
            A ramp is a DIFFERENCE between ends, so report it as one rather than asking the
            reader to eyeball 107 numbers."""
            vals = [v for v in values if v is not None]
            if len(vals) < 20:
                return None
            k = max(1, int(len(vals) * frac))
            return round(st.mean(vals[-k:]) - st.mean(vals[:k]), 6)

        d = {"seed": seed, "pairs": pairs, "overlapping_pairs": overlapping,
             "overlap_frac": round(overlapping / pairs, 4) if pairs else None,
             "quality": {}, "overlap_mismatch": {}, "ramps": {}}
        for mate in (1, 2):
            mq = prof(qsum[mate], qn[mate])
            b20 = prof(below20[mate], qn[mate])
            mmr = prof(ovmm[mate], ovn[mate])
            d["quality"][f"mean_q_r{mate}"] = mq
            d["quality"][f"frac_below_q20_r{mate}"] = b20
            d["quality"][f"frac_below_q30_r{mate}"] = prof(below30[mate], qn[mate])
            d["overlap_mismatch"][f"rate_r{mate}"] = mmr
            d["overlap_mismatch"][f"observations_r{mate}"] = ovn[mate]
            d["ramps"][f"mean_q_r{mate}"] = ramp(mq)
            d["ramps"][f"frac_below_q20_r{mate}"] = ramp(b20)
            d["ramps"][f"overlap_mismatch_r{mate}"] = ramp(mmr)
            tot_n = sum(ovn[mate])
            d["overlap_mismatch"][f"overall_rate_r{mate}"] = \
                round(sum(ovmm[mate]) / tot_n, 6) if tot_n else None
        if idx is not None:
            d["mate_asymmetry"] = {
                "both_mates_map": both, "only_r1_maps": only1, "only_r2_maps": only2,
                "neither_maps": neither,
                "asymmetric_pairs": only1 + only2,
                "frac_of_asymmetric_where_r2_is_the_unmapped_one":
                    round(only1 / (only1 + only2), 4) if (only1 + only2) else None,
                "length_matched": {
                    "r1_truncated_to_mate_2_length": True,
                    "both_mates_map": both_t, "only_r1_maps": only1_t,
                    "only_r2_maps": only2_t, "neither_maps": neither_t,
                    "asymmetric_pairs": only1_t + only2_t,
                    "frac_of_asymmetric_where_r2_is_the_unmapped_one":
                        round(only1_t / (only1_t + only2_t), 4)
                        if (only1_t + only2_t) else None,
                },
            }
        per_seed.append(d)
        r = d["ramps"]
        print(f"seed {seed}: overlap {d['overlap_frac']}  "
              f"Q ramp r1 {r['mean_q_r1']} r2 {r['mean_q_r2']}  "
              f"mismatch ramp r1 {r['overlap_mismatch_r1']} r2 {r['overlap_mismatch_r2']}  "
              f"overall mismatch r1 {d['overlap_mismatch']['overall_rate_r1']} "
              f"r2 {d['overlap_mismatch']['overall_rate_r2']}", flush=True)

    def across(path):
        vals = []
        for d in per_seed:
            cur = d
            for p in path:
                cur = cur[p]
            if cur is not None:
                vals.append(cur)
        if not vals:
            return None
        return {"mean": round(st.mean(vals), 6), "min": min(vals), "max": max(vals),
                "seeds": len(vals)}

    summary = {
        "quality_ramp": {k: across(["ramps", k]) for k in per_seed[0]["ramps"]},
        "overall_overlap_mismatch_rate": {
            f"r{m}": across(["overlap_mismatch", f"overall_rate_r{m}"]) for m in (1, 2)},
        "mate_asymmetry": {
            k: across(["mate_asymmetry", k]) for k in
            ("only_r1_maps", "only_r2_maps",
             "frac_of_asymmetric_where_r2_is_the_unmapped_one")}
        if "mate_asymmetry" in per_seed[0] else None,
        "mate_asymmetry_length_matched": across(
            ["mate_asymmetry", "length_matched",
             "frac_of_asymmetric_where_r2_is_the_unmapped_one"])
        if "mate_asymmetry" in per_seed[0] else None,
    }

    out = {"seeds": args.seeds, "pairs_per_seed": args.pairs, "per_seed": per_seed,
           "summary": summary,
           "reading_rules": [
               "a ramp is the last quarter of the read minus the first quarter; positive "
               "means the end of the read is worse",
               "the overlap instrument does not trust the quality strings, which is the "
               "point: NovaSeq bins quality, so a flat Q profile is not by itself evidence "
               "of a flat error rate",
               "the overlap mismatch rate is the SUM of both mates' error rates at the "
               "paired positions, so the level is about twice a single mate's error rate; "
               "the shape and the comparison between mates are the readouts",
               "only overlapping pairs enter the overlap instrument, which selects short "
               "molecules, and a mismatch cannot be blamed on one mate without using the "
               "quality scores this instrument avoids",
               "the two mates' overlap profiles are ONE mismatch profile in two coordinate "
               "systems: mate 1's late positions pair with mate 2's early ones, so a ramp "
               "in one appears inverted in the other and the overall rates are equal by "
               "construction. Only the shape, and whether any positional structure exists "
               "at all, is a readout",
               "the mate-mapping asymmetry must be read from the length_matched block: "
               "mate 2 is 93 nt against mate 1's 107 after the non-template trim, and a "
               "shorter read maps less often for reasons unrelated to quality",
           ],
           "wall_seconds": time.perf_counter() - t0}
    (args.out / "quality_profile.json").write_text(json.dumps(out, indent=2))
    (args.out / "manifest.json").write_text(json.dumps({
        "script": "experiments/quality_profile.py", **provenance(__file__),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform()}, indent=2))

    print("\nacross seeds:")
    for k, v in summary["quality_ramp"].items():
        print(f"  ramp {k:28} {v}")
    for k, v in summary["overall_overlap_mismatch_rate"].items():
        print(f"  overall overlap mismatch {k}: {v}")
    if summary["mate_asymmetry"]:
        for k, v in summary["mate_asymmetry"].items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
