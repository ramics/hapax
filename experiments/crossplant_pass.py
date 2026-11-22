"""On-instance phase 2 of the cross-plant pass: scan every fixed-depth subject once.

Reuses `recurrence_pass.py`'s scan and parser unchanged, against the FASTQ subsets that
`crossplant_sample.py` kept, so the only differences from the Berlin recurrence pass are the
panel, the fixed subject depth and the per-library R2 trim.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from recurrence_pass import read_rows, scan


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--subjects", nargs="+", required=True, help="ACC:TRIM_R2")
    ap.add_argument("--subject-dir", type=Path, required=True)
    ap.add_argument("--tool", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--timings-name", default="pass_timings.json")
    ap.add_argument("--budget-hours", type=float, default=0.0)
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    names = [line.split("\t", 1)[0] for line in open(args.queries) if line.strip()]
    t_start = time.perf_counter()
    done, timings = [], {}
    for spec in args.subjects:
        acc, trim = spec.split(":")
        npz = args.out / f"witness_{acc}.npz"
        if npz.exists():
            done.append(acc)
            continue
        if args.budget_hours and (time.perf_counter() - t_start) / 3600 > args.budget_hours:
            print(f"! budget reached before {acc}", flush=True)
            break
        r1 = args.subject_dir / f"{acc}_1.fastq.gz"
        r2 = args.subject_dir / f"{acc}_2.fastq.gz"
        if not (r1.exists() and r2.exists()):
            print(f"! {acc} has no subset; skipped", flush=True)
            continue
        t0 = time.perf_counter()
        tsv = args.work / f"{acc}.cov.tsv"
        scan(args.tool, r1, r2, int(trim), args.queries, tsv)
        arrays = read_rows(tsv, names)
        np.savez_compressed(npz, **arrays)
        tsv.unlink()
        timings[acc] = {"scan": round(time.perf_counter() - t0, 1),
                        "library_pairs": int(arrays["library_pairs"][0])}
        done.append(acc)
        print(f"+ {acc} {timings[acc]}", flush=True)
        (args.out / args.timings_name).write_text(json.dumps(
            {"queries": len(names), "subjects_done": done, "timings": timings,
             "seconds_elapsed": round(time.perf_counter() - t_start, 1)}, indent=2) + "\n")
    if not done:
        raise SystemExit("no subject completed")


if __name__ == "__main__":
    main()
