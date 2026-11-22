"""Pack the per-subject witness files into query x subject bit matrices for download.

Per mate: `any` (a 31-mer present) and `most` (half or more present) as packed bits, plus the
median 31-mer count per query and subject capped at 65535. A run's 540 npz files are a few
GB; this is what the report reads.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--witness", type=Path, required=True)
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    op = gzip.open if args.queries.suffix == ".gz" else open
    with op(args.queries, "rt") as h:
        names = [line.split("\t", 1)[0] for line in h if line.strip()]
    subs = sorted(p.stem.removeprefix("witness_") for p in args.witness.glob("witness_*.npz"))
    nq, ns = len(names), len(subs)
    bits = {k: np.zeros((ns, (nq + 7) // 8), dtype=np.uint8)
            for k in ("any1", "any2", "most1", "most2")}
    med = {k: np.zeros((ns, nq), dtype=np.uint16) for k in ("med1", "med2")}
    pairs = np.zeros(ns, dtype=np.int64)
    for j, acc in enumerate(subs):
        z = np.load(args.witness / f"witness_{acc}.npz")
        for m in ("1", "2"):
            f = np.nan_to_num(z[f"frac{m}"], nan=0.0)
            bits[f"any{m}"][j] = np.packbits(f > 0)
            bits[f"most{m}"][j] = np.packbits(f >= 0.5)
            med[f"med{m}"][j] = np.minimum(z[f"med{m}"], 65535)
        pairs[j] = int(z["library_pairs"][0])
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out / "packed.npz", subjects=np.array(subs), nq=nq,
                        library_pairs=pairs, **bits, **med)
    (args.out / "query_names.txt").write_text("\n".join(names) + "\n")
    print(json.dumps({"queries": nq, "subjects": ns}))


if __name__ == "__main__":
    main()
