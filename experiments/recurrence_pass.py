"""Score one query set against many panel libraries: the cross-sample recurrence pass.

WHAT RUNS WHERE. This is the on-instance half. The query set is built in the container by
`recurrence_queries.py` and uploaded; the panel is 130 GB of FASTQ in `sra-pub-run-odp`,
which is in-region on EC2 and about 1 MB/s from the container, so the scan has to happen
beside the data. Per subject library: stage the SRA object, make FASTQ, scan it once with
`experiments/kmer_coverage` through the query-only path, write the per-query result, delete
the reads. Disk never holds more than one subject.

WHAT IS MEASURED PER QUERY AND SUBJECT. Not a yes or no. For each mate separately - because
mate 2 is 93 nt against mate 1's 107 after the trim and a pooled figure averages two
instruments (CLAUDE.md) - the fraction of the mate's 31-mers that occur at all in the
subject, and the median number of times they occur. Presence answers "is this molecule
there"; the median count answers "how abundant is it there", which is what makes occupancy
a shape rather than a bitmask, and SPEC.md §4.2 asks for the shape.

WHY THE COUNT MATTERS AS MUCH AS THE PRESENCE. A molecule witnessed once in forty libraries
and a molecule abundant in two are both "present in n samples" under a threshold, and they
are different organisms. The pass keeps both, and nothing here thresholds either: the
analysis reports distributions with the identified and ribosomal piles drawn beside the
residue, per this repo's rule that a reference cloud travels with every figure.

THE SUBJECT'S OWN TRIM. Every panel library is the same protocol as ours (SMARTer Stranded
Total RNA-Seq v3 Pico), whose mate 2 carries a non-template prefix measured at mode 14 on two
of these libraries (`results/nontemplate_prefix/`). It is trimmed on the subject side too.
This matters less than on the query side - an untrimmed subject donates a few extra k-mers to
the index and can only make a query look MORE witnessed - but the direction is stated rather
than left to the reader, and it is the same operating point on both sides.

IDEMPOTENCE IS A COST PROPERTY HERE. A subject whose output already exists is skipped, so a
run that dies at library 60 of 116 resumes rather than paying for the first 59 again. The
cost of this pass is one library scan each and nothing else, which is why the query set can
be generous.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HEADER = ("block name n1 n2 med1 med2 q10_1 q90_1 q10_2 q90_2 zero1 zero2 step1 step2 "
          "zrun1 zrun2 stepnz1 stepnz2 library_pairs").split()


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    print("+ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def stage_sra(accession: str, work: Path) -> Path:
    """Copy one run's SRA object out of the public bucket.

    Unsigned first and signed as the fallback, because which one works depends on where this
    is running: from the container the regional and signed paths are refused and the global
    unsigned one works, and on EC2 there is no egress proxy signing our principal onto a
    public bucket. Both are tried rather than chosen, since the failure is an `AccessDenied`
    that reads exactly like a missing grant (`ec2_common.py`'s module docstring).
    """
    dest = work / f"{accession}.sra"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    uri = f"s3://sra-pub-run-odp/sra/{accession}/{accession}"
    for extra in (["--no-sign-request"], []):
        try:
            run(["aws", "s3", "cp", uri, str(dest), *extra, "--only-show-errors"])
            return dest
        except subprocess.CalledProcessError as exc:
            last = exc
    raise SystemExit(f"{accession} did not stage from {uri}: {last}")


def to_fastq(sra: Path, work: Path, threads: int) -> tuple[Path, Path]:
    """`fasterq-dump` to split FASTQ, with the split asserted.

    A run that comes back single-ended would leave mate 2 missing, and the scan would then
    score every query's mate 2 against nothing and report it as absent - a zero column
    beside a varying one, which this project has already learned to read as absence rather
    than as a measurement of zero.
    """
    acc = sra.stem
    run(["fasterq-dump", "--split-files", "--threads", str(threads), "--outdir", str(work),
         "--temp", str(work), str(sra)])
    r1, r2 = work / f"{acc}_1.fastq", work / f"{acc}_2.fastq"
    if not (r1.exists() and r2.exists()):
        raise SystemExit(f"{acc} did not produce two mate files; it is not paired as the "
                         f"panel inventory recorded, so it does not belong in the panel")
    return r1, r2


def scan(tool: Path, r1: Path, r2: Path, trim: int, queries: Path, out_tsv: Path) -> None:
    # 2,000 subject pairs are also scored, as block 1, because the tool's interface requires
    # a library block. They are NOT a control: that path subtracts a pair's own occurrences
    # and the query-only path does not, so a subject read is scored more harshly than a
    # query (`cross_library_witness.py`). They are ignored on the way out.
    with open(out_tsv, "w") as o, open(str(out_tsv) + ".err", "w") as e:
        subprocess.run([str(tool), str(r1), str(r2), str(trim), "2000", "0", "-",
                        str(queries)], stdout=o, stderr=e, check=True)


def read_rows(tsv: Path, names: list[str]) -> dict[str, np.ndarray]:
    """Parse the tool's output into arrays in the query file's order.

    The header is asserted, not assumed, and the join is by name with a completeness check:
    a name-keyed join that silently drops non-matches is what produced `partial_eyeball/`'s
    small, plausible, wrong file.
    """
    ix = {k: i for i, k in enumerate(HEADER)}
    pos = {n: i for i, n in enumerate(names)}
    n = len(names)
    frac1 = np.full(n, np.nan, dtype=np.float32)
    frac2 = np.full(n, np.nan, dtype=np.float32)
    med1 = np.zeros(n, dtype=np.uint32)
    med2 = np.zeros(n, dtype=np.uint32)
    kmers1 = np.zeros(n, dtype=np.uint16)
    kmers2 = np.zeros(n, dtype=np.uint16)
    seen = 0
    library_pairs = 0
    with open(tsv) as h:
        head = h.readline().rstrip("\n").split("\t")
        if head != HEADER:
            raise SystemExit(f"{tsv}: unexpected header {head[:6]}...; the tool's output "
                             f"format changed and every column below would be wrong")
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[ix["block"]] != "9":
                continue
            i = pos.get(f[ix["name"]])
            if i is None:
                continue
            seen += 1
            library_pairs = int(f[ix["library_pairs"]])
            for tag, (frac, med, kmers) in {"1": (frac1, med1, kmers1),
                                            "2": (frac2, med2, kmers2)}.items():
                z = f[ix[f"zero{tag}"]]
                kmers[i] = int(f[ix[f"n{tag}"]])
                if z in ("NaN", "nan", ""):
                    continue
                value = 1.0 - float(z)
                if not -1e-6 <= value <= 1 + 1e-6:
                    raise SystemExit(f"{tsv}: present fraction {value} is impossible, so a "
                                     f"column is being read wrong")
                frac[i] = min(max(value, 0.0), 1.0)
                med[i] = int(f[ix[f"med{tag}"]])
    if seen < 0.99 * n:
        raise SystemExit(f"{tsv} scored {seen} of {n} queries. The query file and this "
                         f"output are not the same set, and a short join would still "
                         f"produce a tidy distribution.")
    return {"frac1": frac1, "frac2": frac2, "med1": med1, "med2": med2,
            "kmers1": kmers1, "kmers2": kmers2,
            "library_pairs": np.array([library_pairs], dtype=np.int64)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--subjects", nargs="+", required=True,
                    help="run accessions, in the order they should be scanned")
    ap.add_argument("--trim-r2", type=int, required=True,
                    help="mate 2 bases trimmed from every subject library scanned. "
                         "Required, with no default: this pass runs over panel "
                         "libraries that experiments/libraries.json has no entries "
                         "for, so the caller states the value. "
                         "experiments/ec2_recurrence.py resolves it from the query "
                         "run and passes it in the UserData")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("THREADS", "8")))
    ap.add_argument("--tool", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("/mnt/work"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep-tsv", action="store_true",
                    help="keep the tool's raw output per subject; off by default because it "
                         "is ~10 MB per subject of what the npz already holds")
    ap.add_argument("--timings-name", default="pass_timings.json",
                    help="the progress file's name. Shards write their own: one shared file "
                         "would have each shard's write erase the others' progress, and the "
                         "run would look like it had done a fraction of what it had")
    ap.add_argument("--budget-hours", type=float, default=0.0,
                    help="stop starting new subjects after this long. 0 means no limit; the "
                         "pilot sets it so the measurement of the rate cannot itself run "
                         "away, which is the failure the big-index run's ordering exists to "
                         "prevent")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    names = [line.split("\t", 1)[0] for line in open(args.queries) if line.strip()]
    if not names:
        raise SystemExit(f"{args.queries} is empty")
    (args.out / "query_names.txt").write_text("\n".join(names) + "\n")
    print(f"{len(names):,} queries, {len(args.subjects)} subjects", flush=True)

    t_start = time.perf_counter()
    done, timings = [], {}
    for acc in args.subjects:
        npz = args.out / f"witness_{acc}.npz"
        if npz.exists():
            print(f"= {acc} already done", flush=True)
            done.append(acc)
            continue
        if args.budget_hours and (time.perf_counter() - t_start) / 3600 > args.budget_hours:
            print(f"! budget of {args.budget_hours}h reached before {acc}; stopping with "
                  f"{len(done)} subjects done", flush=True)
            break
        t = {}
        t0 = time.perf_counter()
        try:
            sra = stage_sra(acc, args.work)
            t["stage"] = round(time.perf_counter() - t0, 1)
            t1 = time.perf_counter()
            r1, r2 = to_fastq(sra, args.work, args.threads)
            t["fastq"] = round(time.perf_counter() - t1, 1)
            t2 = time.perf_counter()
            tsv = args.work / f"{acc}.cov.tsv"
            scan(args.tool, r1, r2, args.trim_r2, args.queries, tsv)
            t["scan"] = round(time.perf_counter() - t2, 1)
            arrays = read_rows(tsv, names)
            np.savez_compressed(npz, **arrays)
            if args.keep_tsv:
                shutil.copy(tsv, args.out / f"{acc}.cov.tsv")
            t["library_pairs"] = int(arrays["library_pairs"][0])
            t["total"] = round(time.perf_counter() - t0, 1)
            done.append(acc)
            print(f"+ {acc} {t}", flush=True)
        finally:
            # Disk is the silent killer on a fixed volume, and a failure that leaves 20 GB of
            # FASTQ behind turns one broken subject into a broken run.
            for leftover in args.work.glob(f"{acc}*"):
                if leftover.is_file():
                    leftover.unlink()
                else:
                    shutil.rmtree(leftover, ignore_errors=True)
        timings[acc] = t
        (args.out / args.timings_name).write_text(json.dumps(
            {"queries": len(names), "subjects_done": done, "timings": timings,
             "seconds_elapsed": round(time.perf_counter() - t_start, 1)}, indent=2) + "\n")

    print(json.dumps({"subjects_done": len(done),
                      "seconds": round(time.perf_counter() - t_start, 1)}, indent=2))
    if not done:
        raise SystemExit("no subject completed; the pass produced nothing")


if __name__ == "__main__":
    main()
