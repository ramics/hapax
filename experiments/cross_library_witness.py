"""Does one library's sequence explain the other's unidentified reads?

THE QUESTION, AND WHY IT IS THE ONE WORTH ASKING OF A MATCHED PAIR. `dna_rna_pair.py`
establishes that the DNA member of the pair is a virome and that 99% of it is unidentifiable
by any index we have. That is interesting but it is not a reason to sequence both. The reason
would be this: if the DNA library contains the molecules our RNA residue cannot account for,
then a second extraction from the same water is a cheaper subtrahend than any public
database, which is exactly what this project's thesis predicts for sample-local data
(1.4 MB of the sample's own rRNA contigs beat all 272 MB of SILVA).

HOW IT IS MEASURED. `experiments/kmer_coverage/` indexes a whole library's exact canonical
31-mers and scores query pairs against it without adding them to the library, which is what
its block-9 path exists for. So the query set is reads from the 2022 RNA library, split by
whether Kraken2 identified them, and the subject is a whole library. The readout is the
fraction of a read's 31-mers that occur at least once in the subject library - a deliberately
weak criterion, because the question here is whether the material is PRESENT at all, not how
deeply.

THE CONTROL IS THE WHOLE EXPERIMENT. "A quarter of unidentified reads are witnessed in the
DNA library" means nothing on its own: unidentified reads might be witnessed in any library,
or in none, for reasons that have nothing to do with nucleic acid. So the same query set is
scored against a second subject - our 2021 RNA library, a different sample of the same plant
- and what is reported is the DIFFERENCE between what the two subjects explain. Without it
this measurement cannot distinguish "the DNA library holds the residue" from "the residue is
ordinary sequence that recurs everywhere".

WHAT IS NOT CONTROLLED. The query reads are RNA, so a same-nucleic-acid subject shares
ribosomal sequence with them that a virome cannot. That inflates the 2021 subject, and the
direction is stated with every number rather than corrected, because the correction would
need a third library this project does not have.

TWO COLUMN BUGS, ONE CAUGHT BY ARITHMETIC. The first version of this analysis read the
coverage TSV by hard-coded column index and was off by one, which produced a "fraction
present" of -0.71. An impossible value is the cheapest possible bug report; a plausible one
would have shipped. So the header is parsed into a name map here and asserted against the
tool's own header line, and every fraction is checked to be in [0, 1] before it is used.
The same run also showed why the tool's own library blocks cannot serve as the control: for
those it subtracts each pair's own occurrences and for block 9 it does not, so a library
read is scored more harshly than a foreign one and the control moves against the test.
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

#: The tool's output header, asserted rather than assumed. See the docstring.
HEADER = ("block name n1 n2 med1 med2 q10_1 q90_1 q10_2 q90_2 zero1 zero2 step1 step2 "
          "zrun1 zrun2 stepnz1 stepnz2 library_pairs").split()


def parse(path: Path) -> dict[str, list[float]]:
    """{group: [fraction of the pair's 31-mers present in the subject library]}."""
    cols = {k: i for i, k in enumerate(HEADER)}
    out: dict[str, list[float]] = {}
    with open(path) as h:
        head = h.readline().rstrip("\n").split("\t")
        if head != HEADER:
            raise SystemExit(f"{path}: unexpected header {head[:6]}... - the tool's output "
                             f"format changed and every column index below would be wrong")
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[cols["block"]] != "9":  # only the query-only pairs are comparable
                continue
            group = f[cols["name"]].split("|", 1)[0]
            try:
                z1 = float(f[cols["zero1"]])
                z2 = float(f[cols["zero2"]])
            except ValueError:
                continue
            if z1 != z1 or z2 != z2:  # a mate with no valid k-mer
                continue
            present = 1.0 - (z1 + z2) / 2
            if not 0.0 <= present <= 1.0:
                raise SystemExit(f"{path}: fraction present {present} is impossible, so a "
                                 f"column is being read wrong")
            out.setdefault(group, []).append(present)
    return out


def summarise(v: list[float]) -> dict:
    n = len(v)
    return {"pairs": n, "mean_frac_present": round(st.mean(v), 4),
            "median_frac_present": round(st.median(v), 4),
            "frac_with_half_present": round(sum(1 for w in v if w >= 0.5) / n, 4),
            "frac_with_none_present": round(sum(1 for w in v if w <= 0.0) / n, 4)}


def build_queries(reads1: Path, reads2: Path, calls: Path, out: Path, pairs: int,
                  per_group: int, trim_r2: int, seed: int) -> dict:
    """Query set: pairs of one library, split by whether Kraken2 identified them.

    The split comes from the classifier's own per-read output rather than being recomputed,
    so the groups are exactly the ones `dna_rna_pair.py` reports a rate for.
    """
    classified = {}
    with open(calls) as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            classified[f[1]] = f[2] != "0"
    rows = list(read_pairs(reads1, reads2, pairs, trim_r2=trim_r2, skip=0))
    groups = {"identified": [r for r in rows if classified.get(r[0])],
              "unidentified": [r for r in rows if not classified.get(r[0])]}
    rng = random.Random(seed)
    counts = {}
    with open(out, "w") as h:
        for tag, group in groups.items():
            rng.shuffle(group)
            take = group[:per_group]
            counts[tag] = len(take)
            for r in take:
                h.write(f"{tag}|{r[0]}\t{r[1]}\t{r[2]}\n")
    return {"available": {k: len(v) for k, v in groups.items()}, "written": counts}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--query-r1", type=Path,
                    default=ROOT / "data/berlin_pair/SRR23973265_1.fastq.gz")
    ap.add_argument("--query-r2", type=Path,
                    default=ROOT / "data/berlin_pair/SRR23973265_2.fastq.gz")
    ap.add_argument("--query-trim-r2", type=int, default=16)
    ap.add_argument("--calls", type=Path, required=True,
                    help="Kraken2 per-read output for the query library, from the same "
                         "block: dna_rna_pair.py leaves it in its work dir")
    ap.add_argument("--pairs", type=int, default=200_000)
    ap.add_argument("--per-group", type=int, default=10_000)
    ap.add_argument("--subject", action="append", nargs=4,
                    metavar=("LABEL", "R1", "R2", "TRIM_R2"), required=True,
                    help="a subject library to score against; give it at least twice, "
                         "because a single subject's number cannot be read")
    ap.add_argument("--tool", type=Path,
                    default=ROOT / "experiments/kmer_coverage/target/release/kmer_coverage")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/crosswitness"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/dna_rna_pair")
    args = ap.parse_args()

    if len(args.subject) < 2:
        raise SystemExit("one subject library is not a measurement; pass --subject twice")
    if not args.tool.exists():
        raise SystemExit(f"{args.tool} is not built; cargo build --release in its directory")

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    q = args.work / f"query_s{args.seed}.tsv"
    qmeta = build_queries(args.query_r1, args.query_r2, args.calls, q, args.pairs,
                          args.per_group, args.query_trim_r2, args.seed)
    print(f"query set: {qmeta}", flush=True)

    results, timings = {}, {}
    for label, r1, r2, trim in args.subject:
        tsv = args.work / f"witness_{label}_s{args.seed}.tsv"
        t = time.perf_counter()
        with open(tsv, "w") as out, open(args.work / f"{label}.err", "w") as err:
            # 2000 library pairs are scored as block 1 as well. They are NOT the control -
            # the tool subtracts their own occurrences and does not for block 9 - and are
            # kept only because the tool requires a block specification.
            subprocess.run([str(args.tool), r1, r2, str(trim), "2000", "0", "-", str(q)],
                           stdout=out, stderr=err, check=True)
        timings[label] = time.perf_counter() - t
        groups = parse(tsv)
        results[label] = {g: summarise(v) for g, v in sorted(groups.items())}
        results[label]["gap_unidentified_minus_identified"] = round(
            results[label]["unidentified"]["mean_frac_present"]
            - results[label]["identified"]["mean_frac_present"], 4)
        print(f"{label}: " + ", ".join(
            f"{g} {d['mean_frac_present']:.4f}" for g, d in results[label].items()
            if isinstance(d, dict)) +
            f", gap {results[label]['gap_unidentified_minus_identified']:+.4f}", flush=True)

    out = {"script": "experiments/cross_library_witness.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "query": qmeta, "subjects": results,
           "readout": "fraction of a pair's exact canonical 31-mers occurring at least "
                      "once in the subject library, averaged over the two mates",
           "bias": "the query reads are RNA, so a same-nucleic-acid subject shares "
                   "ribosomal sequence with them that a virome cannot. That inflates the "
                   "RNA subject; no third library exists here to correct it",
           "timings": timings, "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    (args.out / f"cross_witness_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    print(f"written {args.out}/cross_witness_s{args.seed}.json", flush=True)


if __name__ == "__main__":
    main()
