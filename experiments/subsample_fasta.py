"""Deterministic every-Nth-record subsample of a (gzipped) protein FASTA.

Why this exists. The open question left by the rho_fp measurement is whether the
high false-flag correlation is a property of the flags or an artefact of an
incomplete protein index: against RefSeq viral + Swiss-Prot, 86% of Berlin
background reads have *every* peptide 10-mer absent, so there is little room for
the score to vary. The clean answer would be the same measurement against a
comprehensive index, and that does not fit here - UniRef50 is ~20 Gaa, i.e.
~2e10 peptide 10-mers, about 160 GB as sorted hashes against 50 GB of free disk.

What fits is a scaling curve: measure the absence rate and the mate correlation
at several index sizes and see how fast they move. Subsamples give the extra
points. Every Nth record rather than a random sample, so the choice is
reproducible from the file alone with no seed state, and because UniRef50's
record order is by cluster identifier, which is unrelated to the taxonomy or
length of the member sequences.

The subsample is a *representative fraction of protein space*, not a smaller
complete index: it keeps the diversity of the full file while reducing the
k-mer count roughly in proportion. That is what the curve needs - the question
is how the absence rate falls as coverage of protein space grows.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--every", type=int, required=True,
                    help="keep one record in every N")
    args = ap.parse_args()

    opener_in = gzip.open if args.inp.suffix == ".gz" else open
    opener_out = gzip.open if args.out.suffix == ".gz" else open
    kept = seen = residues = 0
    keep = False
    with opener_in(args.inp, "rt") as fh, opener_out(args.out, "wt") as oh:
        for line in fh:
            if line[0] == ">":
                keep = (seen % args.every) == 0
                seen += 1
                kept += keep
            if keep:
                oh.write(line)
                if line[0] != ">":
                    residues += len(line) - 1
    print(f"{args.inp.name}: {seen} records -> kept {kept} ({residues} residues) "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
