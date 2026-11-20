"""What a gene caller does when handed pure ribosomal RNA.

`results/rrna_content/` found that ribosomal Berlin reads are less than half as likely to
have every peptide 10-mer absent as non-ribosomal ones (0.125 against 0.286), which chance
cannot produce against an index of real proteins, and that the non-absent ribosomal reads
sit overwhelmingly inside Prodigal-predicted genes. The inference was that the sample's own
predicted proteins include spurious genes called on ribosomal contigs, so the protein index
"explains" ribosomal reads that no real protein could explain.

The other branch is testing that on this sample's assembly, which is the right place to fix
it. This tests the general claim instead, and it is the more portable one: hand pyrodigal
nothing but curated ribosomal RNA from SILVA and measure what it calls. No assembly, no
coverage, no mapping. If a gene caller invents substantial coding sequence on pure rRNA,
then every sample-self protein index built this way carries the same pollution, and the
background-removal curve has to mask rRNA or count those reads separately. If it calls
almost nothing, the Berlin result needs another explanation and the inference was wrong.

Two readouts, because they answer different questions. The fraction of sequences with any
call says how often the trap is sprung; the fraction of nucleotide length inside a call is
the coding density, directly comparable to the 0.92 the other branch measured over the
whole Berlin assembly.

A composition control runs beside it: each sequence is shuffled, preserving length and base
composition but destroying every real motif and the actual rRNA secondary structure. If
shuffled sequence yields the same coding density then the calls are what any sequence of
this GC content would produce and the finding is about gene calling in general; if real rRNA
yields more, ribosomal sequence is specifically prone to it. Either outcome leaves the
consequence for the index intact, so this is a mechanism control rather than a check on the
headline.

pyrodigal is used in meta mode, which is what a metagenomic assembly is annotated with and
what the other branch used, so the comparison is like for like.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import statistics as st
import subprocess
import time
from pathlib import Path

import numpy as np
import pyrodigal

ROOT = Path(__file__).resolve().parents[1]


def sample_sequences(fastas: list[Path], n: int, min_len: int, rng) -> list[tuple[str, str]]:
    """Reservoir-sample sequences, so one pass suffices and no taxon is favoured."""
    res: list[tuple[str, str]] = []
    seen = 0
    for fa in fastas:
        name, cur = None, []
        with gzip.open(fa, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name and cur:
                        seen = _offer(res, n, seen, name, "".join(cur), min_len, rng)
                    name, cur = line[1:].strip(), []
                else:
                    cur.append(line.strip())
            if name and cur:
                seen = _offer(res, n, seen, name, "".join(cur), min_len, rng)
    return res


def _offer(res, n, seen, name, seq, min_len, rng) -> int:
    seq = seq.upper().replace("U", "T")
    if len(seq) < min_len:
        return seen
    if len(res) < n:
        res.append((name, seq))
    else:
        j = rng.integers(0, seen + 1)
        if j < n:
            res[int(j)] = (name, seq)
    return seen + 1


def call_genes(seqs: list[tuple[str, str]], orf_finder) -> dict:
    with_call = 0
    coding_bp = 0
    total_bp = 0
    aa = 0
    genes = 0
    per_seq_density = []
    for _, seq in seqs:
        preds = orf_finder.find_genes(seq)
        total_bp += len(seq)
        cbp = sum(abs(g.end - g.begin) + 1 for g in preds)
        coding_bp += cbp
        aa += sum((abs(g.end - g.begin) + 1) // 3 for g in preds)
        genes += len(preds)
        if len(preds):
            with_call += 1
        per_seq_density.append(cbp / len(seq))
    return {
        "sequences": len(seqs),
        "total_bp": total_bp,
        "fraction_with_any_call": with_call / len(seqs) if seqs else 0.0,
        "coding_density": coding_bp / total_bp if total_bp else 0.0,
        "median_per_sequence_density": st.median(per_seq_density) if per_seq_density else 0.0,
        "genes_called": genes,
        "amino_acids_called": aa,
        "mean_gene_aa": aa / genes if genes else 0.0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sequences", type=int, default=3000)
    ap.add_argument("--min-length", type=int, default=500,
                    help="SILVA holds partial entries; 500 keeps sequences long enough for a "
                         "gene call to be possible at all")
    ap.add_argument("--silva", type=Path, default=ROOT / "data" / "silva")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    fastas = sorted(args.silva.glob("*.fasta.gz"))
    seqs = sample_sequences(fastas, args.sequences, args.min_length, rng)
    print(f"{len(seqs):,} SILVA sequences, "
          f"{sum(len(s) for _, s in seqs):,} bp, median length "
          f"{st.median([len(s) for _, s in seqs]):.0f}")

    finder = pyrodigal.GeneFinder(meta=True)
    real = call_genes(seqs, finder)

    shuffled = []
    for name, seq in seqs:
        arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
        rng.shuffle(arr)
        shuffled.append((name, arr.tobytes().decode()))
    ctrl = call_genes(shuffled, finder)

    out = {
        "seed": args.seed,
        "pyrodigal": pyrodigal.__version__,
        "mode": "meta",
        "reference": "SILVA 138.1 SSURef_NR99 + LSURef_NR99",
        "real_rrna": real,
        "composition_shuffled_control": ctrl,
        "comparison": {
            "berlin_assembly_coding_density": 0.92,
            "note": "0.92 is the other branch's figure over the whole Berlin assembly, "
                    "quoted for scale only; it was measured elsewhere, not here",
        },
    }
    (args.out / f"prodigal_on_rrna_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"prodigal_on_rrna_manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/prodigal_on_rrna.py",
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "pyrodigal": pyrodigal.__version__,
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))

    for label, d in (("real rRNA", real), ("shuffled", ctrl)):
        print(f"{label:>10}: any call {d['fraction_with_any_call']:.4f}  "
              f"coding density {d['coding_density']:.4f}  "
              f"median per-seq {d['median_per_sequence_density']:.4f}  "
              f"{d['genes_called']:,} genes, {d['amino_acids_called']:,} aa, "
              f"mean {d['mean_gene_aa']:.0f} aa")


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
