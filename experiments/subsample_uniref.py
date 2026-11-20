"""Stream UniRef50 over HTTP and keep nested deterministic subsamples, without
ever storing the 9.2 GB file.

Why this exists. The open item left by `mate_rho.py` and by the index-coverage
arm of `results/mate_width/REPORT.md` is the same measurement against a
*comprehensive public* protein index: against RefSeq viral + Swiss-Prot, 86% of
Berlin background reads have every peptide 10-mer absent, and indexing the
sample's own predicted proteins showed that absence is database incompleteness
rather than novelty. The honest version of that curve needs real public indexes
of several sizes.

Two things blocked it and neither needed to. Storing the file failed: this
environment's egress proxy would not carry 9.2 GB to disk (corrupt once, refused
once, a ranged download dead on chunk 3 of 18). And a UniRef50 10-mer index is
~160 GB of sorted hashes, which does not fit. Both dissolve if the file is never
stored and only a fraction is kept: the stream itself is reliable, and a
subsample is a representative fraction of protein space whose k-mer count scales
with it.

The subsamples are nested (every 100th is a subset of every 10th) so the points
on the curve differ only in coverage. Every Nth record rather than a random
sample, so the choice reproduces from the file alone with no seed state, and
because UniRef50's record order is by cluster identifier, which is unrelated to
the taxonomy or length of its members.

Completeness is checkable, which matters here: a first attempt was discarded as
prefix-biased on the strength of a `BadGzipFile` raised after the last record.
UniProt's REST API reports the release's exact cluster count
(`x-total-results` for `query=identity:0.5`), 41,243,524 for this release, so
the record count printed below either matches it or the run is short. Trailing
bytes after the final member are therefore reported, not treated as failure.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

URL = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=URL)
    ap.add_argument("--every", type=int, nargs="+", default=[10, 100],
                    help="keep one record in every N, one output per value")
    ap.add_argument("--out-dir", type=Path, default=Path("data/prot"))
    ap.add_argument("--expect-records", type=int, default=41243524,
                    help="release cluster count from UniProt's x-total-results")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    outs = {n: gzip.open(args.out_dir / f"uniref50_every{n}.fasta.gz", "wt")
            for n in args.every}
    residues = {n: 0 for n in args.every}
    kept = {n: 0 for n in args.every}

    proc = subprocess.Popen(["curl", "-sS", "--fail", "--retry", "8",
                             "--retry-delay", "5", "-L", args.url],
                            stdout=subprocess.PIPE)
    seen = 0
    writing: list = []
    trailing = None
    try:
        with gzip.open(proc.stdout, "rt") as fh:  # type: ignore[arg-type]
            for line in fh:
                if line[0] == ">":
                    writing = [outs[n] for n in args.every if seen % n == 0]
                    for n in args.every:
                        if seen % n == 0:
                            kept[n] += 1
                    seen += 1
                    if seen % 2_000_000 == 0:
                        print(f"  {seen} records", flush=True)
                for oh in writing:
                    oh.write(line)
                if line[0] != ">":
                    for n in args.every:
                        if seen and (seen - 1) % n == 0:
                            residues[n] += len(line) - 1
    except (gzip.BadGzipFile, EOFError) as e:
        trailing = repr(e)
    finally:
        for oh in outs.values():
            oh.close()
        proc.stdout.close()  # type: ignore[union-attr]
        proc.wait()

    complete = seen == args.expect_records
    print(f"records {seen} expected {args.expect_records} "
          f"complete {complete} trailing {trailing}")
    for n in args.every:
        print(f"  every{n}: {kept[n]} records, {residues[n]} residues "
              f"-> {args.out_dir / f'uniref50_every{n}.fasta.gz'}")
    sys.exit(0 if complete else 1)


if __name__ == "__main__":
    main()
