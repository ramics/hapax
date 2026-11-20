"""What the assembler actually sees: shared exact k-mers, not mean divergence.

Section 5 plots pooling against whole-genome divergence, but a de Bruijn
assembler joins two lineages only where they share an exact k-mer, and that
depends on how the differences are spaced rather than on how many there are.
Two pairs at the same mean divergence can differ several-fold in 21-mer survival
if one pair's differences sit in a hypervariable region and the other's are
spread evenly; clustered differences leave long identical stretches that pool
freely. Pair identity is confounded with divergence across the rungs, so this
measurement is also the control for that confound.

Reported per pair: the fraction of the smaller genome's canonical 21-mers also
present in the other (canonical because MEGAHIT canonicalises), and the longest
exactly identical substring, found by binary search over length with a hash set
per step rather than by alignment.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import edlib


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.open():
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


COMP = str.maketrans("ACGT", "TGCA")


def canon_kmers(s: str, k: int) -> set[str]:
    out = set()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        rc = w.translate(COMP)[::-1]
        out.add(min(w, rc))
    return out


def longest_common_substring(a: str, b: str) -> int:
    """Binary search on length; a length is feasible when any substring of that
    length in a also occurs in b. Linear per step, so ~log2(len) steps total."""
    lo, hi = 0, min(len(a), len(b))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        bs = {b[i:i + mid] for i in range(len(b) - mid + 1)}
        if any(a[i:i + mid] in bs for i in range(len(a) - mid + 1)):
            lo = mid
        else:
            hi = mid - 1
    return lo


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, nargs="+", required=True,
                    help="one or more fasta files holding every accession named in --pairs")
    ap.add_argument("--pairs", default="",
                    help="tag:accA:accB, comma separated; default is the rung ladder")
    ap.add_argument("--k", type=int, default=21, help="MEGAHIT's smallest k, where joining starts")
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep/kmer_survival.json"))
    args = ap.parse_args()

    seqs: dict[str, str] = {}
    for f in args.fasta:
        seqs |= read_fasta(f)

    if args.pairs:
        pairs = [tuple(p.split(":")) for p in args.pairs.split(",")]
    else:
        from two_lineage_rung import PAIRS  # noqa
        pairs = [(t, a, b) for t, a, b, _ in PAIRS]

    rows = []
    for tag, acc_a, acc_b in pairs:
        ka = next(k for k in seqs if k.startswith(acc_a.split(".")[0]))
        kb = next(k for k in seqs if k.startswith(acc_b.split(".")[0]))
        a, b = seqs[ka], seqs[kb]
        small, big = (a, b) if len(a) <= len(b) else (b, a)
        sa, sb = canon_kmers(small, args.k), canon_kmers(big, args.k)
        shared = len(sa & sb)
        div = edlib.align(a, b, task="distance")["editDistance"] / max(len(a), len(b))
        row = {
            "pair": tag, "a": ka, "b": kb,
            "divergence": round(div, 4),
            "k": args.k,
            "kmers_smaller": len(sa),
            "shared_kmers": shared,
            "shared_fraction": round(shared / len(sa), 4),
            "longest_identical": longest_common_substring(a, b),
            # uniform-substitution expectation, for contrast with the measurement
            "expected_if_uniform": round((1 - div) ** args.k, 4),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"k": args.k, "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
