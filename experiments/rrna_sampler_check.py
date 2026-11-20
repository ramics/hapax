"""Does the k-mer subsampler actually retain the fraction it claims?

The other branch's contig-level labelling hit a real bug of this kind: an FNV-style hash
over a 4-letter alphabet has non-uniform low bits, so masking the low bits retained 1 k-mer
in 65 rather than the intended 1 in 8. The same mistake would be invisible here — the
headline would still look reasonable, because the sampling applies to reference and query
alike, and only the per-read hit count would quietly collapse.

`experiments/rrna_content.py` masks the low bits of a splitmix64 finalizer applied to the
k-mer's own 62-bit code rather than the code itself, which should avalanche. "Should" is
not a measurement, so this measures it, on both the reference and the reads, and also
reports how many reads carry too few retained k-mers to be callable at all. That last
number is a ceiling on sensitivity that no threshold choice can recover.

Cheap enough to re-run whenever the sampler or k changes, which is the point.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rrna_content as rc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def retention(seqs: list[str], rate_log2: int) -> dict:
    allc, keptc = [], []
    for s in seqs:
        b = rc._MAP[np.frombuffer(s.encode(), dtype=np.uint8)]
        c = rc.codes_from_bases(b)
        allc.append(c)
        keptc.append(rc.keep(c, rate_log2))
    a = np.unique(np.concatenate(allc))
    k = np.unique(np.concatenate(keptc))
    return {"distinct_kmers": int(a.size), "retained": int(k.size),
            "realised_rate": float(k.size / a.size) if a.size else 0.0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="unused: the check is deterministic over file prefixes, but the "
                         "flag is kept so every experiment takes --seed")
    ap.add_argument("--subsample-log2", type=int, default=3)
    ap.add_argument("--min-hits", type=int, default=3)
    ap.add_argument("--silva-sequences", type=int, default=2000)
    ap.add_argument("--reads", type=int, default=20000)
    ap.add_argument("--silva", type=Path,
                    default=ROOT / "data" / "silva" / "SSURef_NR99.fasta.gz")
    ap.add_argument("--fastq", type=Path,
                    default=ROOT / "data" / "berlin" / "prefix_R1.fastq")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    target = 1.0 / (1 << args.subsample_log2)

    seqs, cur = [], []
    with gzip.open(args.silva, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper().replace("U", "T"))
                    cur = []
                if len(seqs) >= args.silva_sequences:
                    break
            else:
                cur.append(line.strip())

    reads = []
    with args.fastq.open() as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                reads.append(line.strip())
            if len(reads) >= args.reads:
                break

    ref = retention(seqs, args.subsample_log2)
    qry = retention(reads, args.subsample_log2)

    per_read = np.array([rc.hits(r, np.empty(0, dtype=np.uint64),
                                 args.subsample_log2)[0] for r in reads])
    out = {
        "target_rate": target,
        "reference": ref,
        "reads": qry,
        "retained_kmers_per_read": {
            "mean": float(per_read.mean()),
            "median": float(np.median(per_read)),
            "min": int(per_read.min()),
            "expected": (args.reads and 0) or 0,
        },
        "uncallable_read_fraction": float(np.mean(per_read < args.min_hits)),
        "note": "uncallable reads carry fewer than min_hits retained k-mers, so they can "
                "never be called ribosomal whatever the threshold; it is a ceiling on "
                "sensitivity, not a tuning parameter",
    }
    out["retained_kmers_per_read"]["expected"] = round(
        (len(reads[0]) - rc.K + 1) * target, 2) if reads else 0.0

    (args.out / "rrna_sampler_check.json").write_text(json.dumps(out, indent=2))
    (args.out / "rrna_sampler_check_manifest.json").write_text(json.dumps({
        "script": "experiments/rrna_sampler_check.py",
        "git_sha": rc._git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))

    print(f"target rate {target:.5f}")
    print(f"  reference: {ref['distinct_kmers']:,} distinct, {ref['retained']:,} retained, "
          f"realised {ref['realised_rate']:.5f}")
    print(f"  reads    : {qry['distinct_kmers']:,} distinct, {qry['retained']:,} retained, "
          f"realised {qry['realised_rate']:.5f}")
    print(f"  per read : mean {out['retained_kmers_per_read']['mean']:.2f}, expected "
          f"{out['retained_kmers_per_read']['expected']:.2f}; uncallable "
          f"{out['uncallable_read_fraction']:.4f}")


if __name__ == "__main__":
    main()
