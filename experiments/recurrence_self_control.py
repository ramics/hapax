"""The instrument's own controls, scored against the library the queries came from.

WHY THIS EXISTS, AND WHY IT IS NOT A RESULT. Every number in the recurrence pass is a
fraction of a read's 31-mers found in some other library, and that quantity has two ends
whose values are known in advance. A read scored against the library it was sequenced from
must be fully witnessed: its own k-mers are in there by construction. A
composition-preserving shuffle of that same read must be witnessed hardly at all, because a
shuffled 107-mer's 31-mers are new sequence. If either end does not read what it must, the
plumbing is wrong - the query file and the library have drifted apart, a column is being read
off by one, the trim differs between the two sides - and every cross-sample figure downstream
would be a measurement of that fault.

This is the control the earlier cross-library witness figure did not have. It is cheap: one
scan of our own library, the same cost as one panel subject.

WHAT WOULD BE A BUG RATHER THAN A RESULT. Both ends reading the SAME extreme value. This
repo's own rule: an identical extreme number in a pile and its control is a bug. Here the
two ends must come out at opposite extremes, and the gap between them is the dynamic range
every later figure lives inside.

NOTE ON THE 1.0: it is not evidence that the pass works on FOREIGN libraries, only that the
join, the parse and the k-mer arithmetic are right. The foreign-library behaviour has no
known answer, which is the whole point of measuring it.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

import library_facts

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from subtract_known import provenance  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded for the convention; the query file names its own seeds")
    ap.add_argument("--queries", type=Path,
                    default=ROOT / "results/recurrence/queries_s1-5.tsv")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
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

    args.work.mkdir(parents=True, exist_ok=True)
    tsv = args.work / "self_witness.tsv"
    t0 = time.perf_counter()
    if not tsv.exists() or tsv.stat().st_size == 0:
        with open(tsv, "w") as o, open(args.work / "self_witness.err", "w") as e:
            subprocess.run([str(args.tool), str(args.r1), str(args.r2), str(args.trim_r2),
                            "2000", "0", "-", str(args.queries)],
                           stdout=o, stderr=e, check=True)
    seconds = round(time.perf_counter() - t0, 1)

    per_pile: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    with open(tsv) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["block"] != "9":
                continue
            pile = r["name"].split("|")[0]
            for mate in ("1", "2"):
                z = r[f"zero{mate}"]
                if z in ("NaN", "nan", ""):
                    continue
                per_pile[pile][f"mate{mate}"].append(1.0 - float(z))

    out: dict = {"seed": args.seed, "piles": {}}
    for pile, mates in sorted(per_pile.items()):
        out["piles"][pile] = {
            m: {"n": len(v), "mean": round(st.mean(v), 6),
                "median": round(st.median(v), 6),
                "frac_fully_witnessed": round(sum(1 for x in v if x > 0.999) / len(v), 6),
                "min": round(min(v), 6), "max": round(max(v), 6)}
            for m, v in sorted(mates.items()) if v}

    real = [p for p in out["piles"] if not p.startswith("shuffled-")]
    null = [p for p in out["piles"] if p.startswith("shuffled-")]
    worst_real = min(out["piles"][p]["mate1"]["mean"] for p in real)
    worst_null = max(out["piles"][p]["mate1"]["mean"] for p in null)
    out["verdict"] = {
        "lowest_real_pile_mean": worst_real,
        "highest_null_pile_mean": worst_null,
        "dynamic_range": round(worst_real - worst_null, 6),
        "reading": "a real read must be fully witnessed in its own library and a shuffle "
                   "hardly at all; the two ends must be at OPPOSITE extremes, because the "
                   "same extreme on both is this repo's signature for a bug",
    }
    if worst_real < 0.99 or worst_null > 0.01:
        raise SystemExit(
            f"the controls do not read what they must: real piles down to {worst_real} "
            f"(must be ~1) and nulls up to {worst_null} (must be ~0). Fix the pipeline "
            f"before any cross-sample number is read.")

    out["manifest"] = {"script": "experiments/recurrence_self_control.py",
                       **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()},
                       "seconds": seconds,
                       "operating_points": {"k": 31, "trim_r2": args.trim_r2}}
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "self_control.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    for pile in sorted(out["piles"]):
        m = out["piles"][pile]["mate1"]
        print(f"{pile:26s} n={m['n']:6d} mean={m['mean']:.4f} min={m['min']:.3f} "
              f"full={m['frac_fully_witnessed']:.4f}")
    print(f"\ndynamic range {out['verdict']['dynamic_range']}  ->  {path}")


if __name__ == "__main__":
    main()
