"""Stream UniRef50 and keep nested random subsets of it, because the whole thing does not fit.

A comprehensive protein index is the far end of the subtraction curve and it is out of
reach here by a wide margin: UniRef50 is ~13 GB compressed, a DIAMOND database of it is
~30 GB, and this container has ~31 GB of writable disk left after the nucleotide indexes.
The sibling thread hit the same wall from the other side (a UniRef50 peptide 10-mer index
would want ~160 GB).

So index size becomes the independent variable rather than a blocker. This script streams
the file over HTTP without ever storing it and writes nested subsets: the 10% subset is a
subset of the 25% subset, which is a subset of UniRef50. Nested rather than independent so
that the removal-per-GB curve is monotone by construction and the three points differ only
in how much of the same catalogue was subtracted - an independent draw at each size would
add sampling noise to the very quantity being measured.

Sampling is by hash of the accession, not by a random number generator, so a re-run picks
the same records without storing a record list, and --seed changes which records are
chosen rather than only the order.

What a subset licenses and what it does not: a 25% subset of UniRef50 bounds the removal
of full UniRef50 from below, and the shape of removal against size over several subsets is
what the extrapolation rests on. It is not a substitute for the full index, and the
extrapolation is reported as one.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

URL = ("https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--url", default=URL)
    ap.add_argument("--fractions", type=float, nargs="+", default=[0.02, 0.10, 0.25])
    ap.add_argument("--out", type=Path, default=Path("data/prot"))
    ap.add_argument("--prefix", default="uniref50")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    fracs = sorted(args.fractions)
    handles = [(f, gzip.open(args.out / f"{args.prefix}_{int(f * 100):02d}pct.fasta.gz",
                             "wt", compresslevel=1)) for f in fracs]
    counts = {f: 0 for f in fracs}
    seqs_seen = 0
    aa_kept = {f: 0 for f in fracs}
    t0 = time.perf_counter()

    proc = subprocess.Popen(["curl", "-sL", "--retry", "20", "--retry-all-errors", args.url],
                            stdout=subprocess.PIPE)
    keep: list = []
    name = None
    buf: list[str] = []

    def flush():
        nonlocal seqs_seen
        if name is None:
            return
        seqs_seen += 1
        seq = "".join(buf)
        u = int(hashlib.blake2b(f"{args.seed}:{name}".encode(), digest_size=8).hexdigest(),
                16) / 2 ** 64
        for f, h in handles:
            if u < f:
                h.write(f">{name}\n{seq}\n")
                counts[f] += 1
                aa_kept[f] += len(seq)

    with gzip.open(proc.stdout, "rt") as g:
        for line in g:
            if line.startswith(">"):
                flush()
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        flush()
    for _f, h in handles:
        h.close()
    proc.stdout.close()
    proc.wait()

    meta = {
        "url": args.url,
        "seed": args.seed,
        "records_streamed": seqs_seen,
        "subsets": [{"fraction": f, "records": counts[f], "residues": aa_kept[f],
                     "file": f"{args.prefix}_{int(f * 100):02d}pct.fasta.gz"}
                    for f in fracs],
        "wall_seconds": time.perf_counter() - t0,
    }
    (args.out / f"{args.prefix}_subsets.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
