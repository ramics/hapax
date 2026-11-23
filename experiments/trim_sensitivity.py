"""How much does the mate-2 trim move the numbers this project reports?

WHY THIS EXISTS, AND IT IS NOT THE QUESTION I SET OUT TO ASK. `nontemplate_prefix_map.py`
established that mate 2's non-template prefix has a 5' soft-clip mode of exactly 14 in 5/5
seeds on both RNA libraries. `dna_rna_pair.py` had been trimming 16, taken from a composition
detector, so I corrected it to the measured 14 and re-ran - and the Kraken2 removal fraction
went DOWN, by 1.66 points on the 2022 RNA library and 1.06 on the 2021 one, from being given
two MORE bases of what the mapping says is template.

That is the wrong direction for a correction, so "16 was an over-trim" was too strong a
reading and is withdrawn. Two bases of real template cannot make a read harder to classify;
what can is that **the prefix length is a distribution, not a constant**. Its mode is 14 at
51% of anchored reads, with 13 at 16%, 15 at 6.5%, 16 at 2.2% and a tail. Kraken2's
confidence is the FRACTION of a read's minimizers falling in the called clade, so two extra
minimizers at the end of the read where the prefix variability lives lower the confidence of
every read a little and push marginal ones below the threshold. Trimming at the mode is the
best estimate of the prefix; trimming above it is the better operating point for a
fraction-thresholded classifier. Those are different objectives and the mode does not settle
the second one.

SO THE TRIM IS AN OPERATING POINT AND IT GETS MEASURED LIKE ONE. This repo's rule is that
every number ships with its operating point, and 1.66 points is the same size as effects it
reports as findings - the Kraken2 confidence sweep moves 42.8/21.7/3.5% across three
settings, and the eight-month drift between two samples of one plant is 18 points. A 2-base
parameter choice worth 1.7 points belongs on the record as a curve rather than as a decision
somebody made once.

WHAT IS CONTROLLED. Mate 1 is untouched at every rung, so the curve isolates mate 2. The same
block of the same library is used at every rung, so nothing varies but the trim. Both the
identified fraction and the mean read length are reported, because a trim that removes real
sequence and a trim that removes junk both shorten the read and only the classification
distinguishes them: removing junk should raise the identified fraction while shortening the
read, and removing template should lower it.

WHAT IS NOT CONTROLLED. One library per run and one confidence setting, so this says how
sensitive the number is here and not what the trim should be in general. The ribosomal
fraction is not measured at each rung, which would need the SILVA index and roughly ten times
the wall time; the classification fraction is the cheaper proxy and the one the cascade
reports.
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

from berlin_io import read_pairs
from subtract_known import kraken2, provenance, write_fastq

ROOT = Path(__file__).resolve().parents[1]

LIBRARIES = {
    "rna_2021": ("data/berlin", "SRR23973193"),
    "rna_2022": ("data/berlin_pair", "SRR23973265"),
    "dna_2022": ("data/berlin_pair", "SRR23973375"),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--library", default="rna_2022", choices=sorted(LIBRARIES))
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--block", type=int, default=None,
                    help="pairs skipped per seed, matching dna_rna_pair.py. From "
                         "experiments/libraries.json for the chosen library unless given")
    ap.add_argument("--trims", type=int, nargs="+", default=[0, 10, 12, 14, 16, 18, 20, 24])
    ap.add_argument("--k2-db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--k2-confidence", type=float, default=0.1,
                    help="the cascade's operating confidence, so the curve is comparable "
                         "to every removal figure in results/")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/"
                                 "libexec/classify"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/trimsens"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/nontemplate_prefix")
    args = ap.parse_args()
    # The stride is a property of the library: sized to 8.73 M pairs, it leaves the
    # last seeds of a smaller library reading nothing, and the directory still looks
    # like a five-seed result. Read from experiments/libraries.json, which raises
    # for a library it does not know.
    args.block = library_facts.resolve(LIBRARIES[args.library][1], "block",
                                       args.block)

    d, run = LIBRARIES[args.library]
    r1, r2 = ROOT / d / f"{run}_1.fastq.gz", ROOT / d / f"{run}_2.fastq.gz"
    for p in (r1, r2, args.k2_db, args.classify_bin):
        if not p.exists():
            raise SystemExit(f"{p} is missing")

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    skip = args.block * (args.seed - 1)

    rungs = []
    for trim in sorted(args.trims):
        rows = list(read_pairs(r1, r2, args.pairs, trim_r2=trim, skip=skip))
        if not rows:
            raise SystemExit(f"no pairs read from {r1}")
        fq1 = args.work / f"t{trim}_1.fq"
        fq2 = args.work / f"t{trim}_2.fq"
        write_fastq(fq1, rows, 1)
        write_fastq(fq2, rows, 2)
        calls, wall, rss = kraken2(args.k2_db, fq1, fq2, args.classify_bin, args.threads,
                                   args.work / f"t{trim}.k2", args.k2_confidence)
        n = len(rows)
        # taxid 0 is Kraken2's "unclassified", and counting it as a call once produced
        # 1.0000 for two libraries at once. Count tx != 0.
        hit = sum(1 for name, tx in calls.items() if tx != 0)
        if n and hit / n >= 0.999:
            raise SystemExit(
                f"trim {trim}: {hit}/{n} classified, which is not a measurement - taxid 0 "
                f"is being counted as a call or the classifier got the wrong file.")
        rung = {"trim": trim, "pairs": n,
                "identified": round(hit / n, 5),
                "mean_len_r1": round(st.mean(len(r[1]) for r in rows), 2),
                "mean_len_r2": round(st.mean(len(r[2]) for r in rows), 2),
                "kraken2_seconds": round(wall, 1), "peak_rss_kb": rss}
        rungs.append(rung)
        print(f"trim {trim:3d}: identified {rung['identified']:.5f}  "
              f"r1 {rung['mean_len_r1']} r2 {rung['mean_len_r2']} nt  "
              f"({wall:.0f}s)", flush=True)
        fq1.unlink(missing_ok=True)
        fq2.unlink(missing_ok=True)

    best = max(rungs, key=lambda r: r["identified"])
    at14 = next((r for r in rungs if r["trim"] == 14), None)
    at16 = next((r for r in rungs if r["trim"] == 16), None)
    out = {"script": "experiments/trim_sensitivity.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "library": {"name": args.library, "run": run},
           "rungs": rungs,
           "best_trim": best["trim"], "best_identified": best["identified"],
           "span_points": round(100 * (best["identified"]
                                       - min(r["identified"] for r in rungs)), 3),
           "measured_mode_14_vs_16": (round(100 * (at16["identified"] - at14["identified"]), 3)
                                      if at14 and at16 else None),
           "readout": "fraction of pairs Kraken2 PlusPF-16 classifies at confidence 0.1, "
                      "against bases removed from mate 2's 5' end. Mate 1 is untouched at "
                      "every rung and the same block is used throughout, so nothing varies "
                      "but the trim",
           "bias": "removing junk and removing template both shorten the read, so length "
                   "cannot tell them apart; the classified fraction can, because junk "
                   "removal should raise it while template removal lowers it. Kraken2's "
                   "confidence is a FRACTION of a read's minimizers, so this curve is a "
                   "property of a fraction-thresholded classifier and not of the prefix",
           "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    p = args.out / f"trim_sensitivity_{args.library}_s{args.seed}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"best trim {best['trim']} at {best['identified']:.5f}; "
          f"span {out['span_points']} points; 16 minus 14 = "
          f"{out['measured_mode_14_vs_16']} points", flush=True)
    print(f"written {p}", flush=True)


if __name__ == "__main__":
    main()
