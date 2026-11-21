"""Does a mate's 5' prefix align with the reference the rest of the read aligns to?

THE CLAIM UNDER TEST. Every experiment in this repo trims 14 bases from the 5' end of mate 2
of SRR23973193, on the stated grounds that they are non-template sequence. That trim is
load-bearing: it is why R2 is 93 nt against R1's 107, which is the cause of every pair
asymmetry CLAUDE.md warns about. It has never been tested against a reference in this repo.
`dna_rna_pair.py` tried to test it by composition and could not: its detector fires on all
four mate-and-library combinations (DNA 14/14, the two RNA libraries 9/16 with near-identical
profiles), and the report records it as UNRESOLVED and the number as a flag.

WHY COMPOSITION CANNOT ANSWER IT AND MAPPING CAN. Illumina's 5' priming bias and a genuine
non-template prefix produce the same composition signature, because both make the first bases
unlike the body of the read. They differ in one respect that composition cannot see: priming
bias changes WHICH molecules get sampled, not whether their first bases came from the
template. So priming-biased bases still align to the genome, and non-template bases cannot.
Map the untrimmed read and the question becomes a clipping question.

THE CONTROL IS THE OTHER END OF THE SAME READ. A read that maps to a divergent reference
entry clips at both ends for divergence, and the level of clipping therefore says nothing on
its own. What a non-template prefix predicts is a DIFFERENCE: a sharp mode at a fixed length
on the 5' end of the affected mate, absent from the 3' end of that same read and from the
other mate. Three controls are reported beside every figure:

  1. the 3' clip of the same mate - the within-read baseline;
  2. the other mate of the same pair, which had the same library prep and the same reference;
  3. the 2022 RNA library beside the 2021 one, so an instrument property and a per-run
     property can be told apart.

The 3' control is conservative rather than neutral: adapter read-through clips 3' ends and
inflates it, which can only make a real 5' excess harder to see.

THE COORDINATE ASSUMPTION IS CHECKED, NOT ASSUMED. This reads mappy's q_st and q_en as
coordinates on the query AS SUPPLIED, so that the 5' clip is q_st for both strands. If that
is wrong - if they were reported on the aligned orientation - then minus-strand hits would
carry the mirror image and the two ends would average into looking alike, which is exactly
the failure mode CLAUDE.md's "read the control first" rule is about. So every distribution is
also reported split by strand, and a mirror between the two strands invalidates the run
rather than being averaged away.

WHAT THIS CANNOT ANSWER. The DNA library is 0.01% ribosomal, so SILVA gives it almost no
confidently mapped reads and its row here is a statement about the reference, not about the
library. No nucleotide reference this container holds identifies enough of a DNA virome to
replace it. That row is reported with its n so it cannot be read as a measurement.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics as st
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from subtract_known import _mm_index, provenance

ROOT = Path(__file__).resolve().parents[1]

LIBRARIES = {
    "rna_2021": {"run": "SRR23973193", "dir": "data/berlin",
                 "note": "the library whose R2 every experiment trims by 14"},
    "rna_2022": {"run": "SRR23973265", "dir": "data/berlin_pair",
                 "note": "same plant, eight months later, same naming series"},
    "dna_2022": {"run": "SRR23973375", "dir": "data/berlin_pair",
                 "note": "0.01% rRNA, so SILVA cannot anchor it; kept to show that"},
}


def clips(aligner, seq: str, min_mlen: int, min_mapq: int):
    """(five_prime_clip, three_prime_clip, strand) for the best anchored hit, or None.

    "Best" is the longest matching block, and the hit must be anchored before its clipping
    means anything: an alignment covering a third of the read has clipping that is a
    statement about the reference, not about the read's ends. Hence both floors.
    """
    best = None
    for h in aligner.map(seq):
        if h.mlen < min_mlen or h.mapq < min_mapq:
            continue
        if best is None or h.mlen > best.mlen:
            best = h
    if best is None:
        return None
    return best.q_st, len(seq) - best.q_en, best.strand


def summarise(v: list[int]) -> dict:
    if not v:
        return {"n": 0}
    c = Counter(v)
    return {"n": len(v), "mode": c.most_common(1)[0][0],
            "mode_share": round(c.most_common(1)[0][1] / len(v), 4),
            "median": st.median(v), "mean": round(st.mean(v), 2),
            "frac_ge_10": round(sum(1 for x in v if x >= 10) / len(v), 4),
            "frac_zero": round(sum(1 for x in v if x == 0) / len(v), 4),
            "top": [[k, round(n / len(v), 4)] for k, n in c.most_common(6)]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--min-mlen", type=int, default=60,
                    help="the body must be anchored before end clipping is readable")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--silva", type=Path, default=ROOT / "data/silva/silva_dna.fa")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--libraries", nargs="+", default=list(LIBRARIES))
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/ntprefix"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/nontemplate_prefix")
    args = ap.parse_args()

    if not args.silva.exists():
        raise SystemExit(f"{args.silva} is missing; data/ is gitignored and fetched locally")

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    idx = _mm_index(args.silva, args.work / "silva.mmi", args.threads)

    results, timings = {}, {}
    for name in args.libraries:
        spec = LIBRARIES[name]
        d = ROOT / spec["dir"]
        r1, r2 = d / f"{spec['run']}_1.fastq.gz", d / f"{spec['run']}_2.fastq.gz"
        for p in (r1, r2):
            if not p.exists():
                raise SystemExit(f"{p} is missing")
        t = time.perf_counter()
        # SEEDS ARE DISJOINT BLOCKS OF THE FILE, NOT A SHUFFLE. read_pairs streams from the
        # start, so a seed that changed only an RNG would re-measure the same reads and five
        # of them would be one result reported five times. Offsetting instead makes the
        # spread a statement about position in the run as well, which is the axis a
        # flowcell-dependent artefact would move along.
        skip = (args.seed - 1) * args.pairs
        # trim_r2=0 is the whole point: the prefix under test must still be on the read.
        five = {1: [], 2: []}
        three = {1: [], 2: []}
        by_strand = {(m, s): [] for m in (1, 2) for s in (1, -1)}
        anchored = {1: 0, 2: 0}
        total = 0
        for _, s1, s2, _, _ in read_pairs(r1, r2, args.pairs, trim_r2=0, skip=skip):
            total += 1
            for mate, seq in ((1, s1), (2, s2)):
                c = clips(idx, seq, args.min_mlen, args.min_mapq)
                if c is None:
                    continue
                f, th, strand = c
                anchored[mate] += 1
                five[mate].append(f)
                three[mate].append(th)
                by_strand[(mate, strand)].append(f)
        timings[name] = time.perf_counter() - t

        row = {"run": spec["run"], "note": spec["note"], "pairs_read": total,
               "file_offset_pairs": skip,
               "anchored": {f"r{m}": anchored[m] for m in (1, 2)},
               "anchored_frac": {f"r{m}": round(anchored[m] / total, 4) if total else None
                                 for m in (1, 2)}}
        for m in (1, 2):
            row[f"r{m}_five_prime_clip"] = summarise(five[m])
            row[f"r{m}_three_prime_clip"] = summarise(three[m])
            row[f"r{m}_five_minus_three_mean"] = (
                round(st.mean(five[m]) - st.mean(three[m]), 2) if five[m] else None)
            # The coordinate check. Equal modes across strands means q_st is on the query as
            # supplied; a mirror means it is not and nothing above can be read.
            row[f"r{m}_five_prime_by_strand"] = {
                ("plus" if s == 1 else "minus"): summarise(by_strand[(m, s)])
                for s in (1, -1)}
        results[name] = row
        print(f"{name}: anchored r1 {anchored[1]} r2 {anchored[2]} of {total}; "
              f"5' mode r1 {row['r1_five_prime_clip'].get('mode')} "
              f"r2 {row['r2_five_prime_clip'].get('mode')}; "
              f"5'-3' mean r1 {row['r1_five_minus_three_mean']} "
              f"r2 {row['r2_five_minus_three_mean']}", flush=True)

    out = {"script": "experiments/nontemplate_prefix_map.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "libraries": results,
           "readout": "soft-clip length at each end of the best anchored SILVA alignment, "
                      "on UNTRIMMED reads. A non-template prefix predicts a sharp mode at a "
                      "fixed length on one mate's 5' end, absent from its own 3' end and "
                      "from the other mate",
           "bias": "adapter read-through clips 3' ends and inflates the within-read "
                   "control, which makes a real 5' excess harder to see, not easier. "
                   "Divergence from the SILVA entry clips both ends and cancels in the "
                   "difference. The DNA library is 0.01% ribosomal so its row is a "
                   "statement about the reference, not the library",
           "coordinate_check": "every 5' distribution is also split by strand; equal modes "
                               "confirm q_st is on the query as supplied, a mirror "
                               "invalidates the run",
           "timings": timings, "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    p = args.out / f"nontemplate_prefix_s{args.seed}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"written {p}", flush=True)


if __name__ == "__main__":
    main()
