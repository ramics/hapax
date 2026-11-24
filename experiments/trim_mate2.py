#!/usr/bin/env python3
"""Write a copy of a paired library with mate 2's non-template prefix removed.

WHY A FILE ON DISK EXISTS AT ALL. Every read-level stage trims at read time:
`berlin_io.read_pairs` cuts `trim_r2` bases off mate 2 as it streams, so no trimmed copy
is needed for the cascade, the per-read report or the linkage stages. The assembler is the
exception. MEGAHIT reads FASTQ from disk, so the only way to keep the non-template bases
out of the graph is to hand it reads that no longer carry them.

WHY IT MATTERS FOR AN ASSEMBLY AND NOT ONLY FOR A SEARCH. The prefix is kit sequence, not
template (`results/nontemplate_prefix/`), so in an assembly graph it is either a tip on
every mate-2 read, which pruning removes along with whatever thin real sequence shares the
tip, or - if the kit's sequence is constant - a motif that appears at one end of every
mate-2 read in the library and joins unrelated molecules through it. The first costs
recall, the second manufactures joins. Neither is visible in a contig afterwards.

WHAT THIS DOES NOT DO. It is not adapter or quality trimming and it does not touch mate 1,
which has ~2 bases of 5' excess and no mode - the difference that separates priming bias
from a non-template prefix (SPEC 4.1). Mate 1 is copied through byte for byte so the two
output files stay in lockstep for the assembler, and the record count is asserted equal.

THE TRIM IS A PROPERTY OF THE LIBRARY. It is resolved through `library_facts`, which
raises for a library this repo has no entry for rather than lending it another library's
value, and `--trim-r2 0` is a legitimate instruction that copies both mates unchanged -
written to the manifest as an explicit zero, which is a different fact from "not trimmed".

OUTPUT IS UNCOMPRESSED on purpose: it is an intermediate the assembler reads once and the
caller deletes, and gzipping 10 Gbase costs more wall clock than the assembly saves.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

import library_facts
from provenance import provenance


def _open(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    library_facts.add_argument(ap, default=None)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="bases to cut from the 5' end of mate 2. Resolved from "
                         "experiments/libraries.json when not given; a library with no "
                         "entry is refused rather than given another library's prefix")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--prefix", required=True,
                    help="output basename; writes <prefix>_1.fq and <prefix>_2.fq")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="where to write the manifest; defaults beside the output")
    args = ap.parse_args()
    library_facts.resolve_args(args, ["trim_r2"])
    if args.trim_r2 < 0:
        raise SystemExit("--trim-r2 cannot be negative")
    for p in (args.r1, args.r2):
        if not p.exists():
            raise SystemExit(f"{p} is missing")

    t0 = time.perf_counter()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    f1 = args.out_dir / f"{args.prefix}_1.fq"
    f2 = args.out_dir / f"{args.prefix}_2.fq"
    n = short = 0
    with _open(args.r1) as h1, _open(args.r2) as h2, \
            open(f1, "w") as o1, open(f2, "w") as o2:
        while True:
            a = [h1.readline() for _ in range(4)]
            b = [h2.readline() for _ in range(4)]
            if not a[0]:
                if b[0]:
                    raise SystemExit(f"{args.r2} has more records than {args.r1}")
                break
            if not b[0]:
                raise SystemExit(f"{args.r1} has more records than {args.r2}")
            seq, qual = b[1].rstrip("\n"), b[3].rstrip("\n")
            if len(seq) <= args.trim_r2:
                # A mate shorter than the prefix would become an empty record, which some
                # readers accept and silently drop: counted and reported rather than
                # written, because a changing record count between the two files is the
                # one thing the assembler cannot be told about.
                short += 1
                seq, qual = "N", "!"
            else:
                seq, qual = seq[args.trim_r2:], qual[args.trim_r2:]
            o1.write(a[0] + a[1] + a[2] + a[3])
            o2.write(f"{b[0]}{seq}\n{b[2]}{qual}\n")
            n += 1
    man = {"script": "experiments/trim_mate2.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "pairs": n, "trim_r2": args.trim_r2,
           "mate2_shorter_than_the_prefix": short,
           "mate2_shorter_note": "written as a single N so the record count matches mate 1",
           "outputs": {str(f1): f1.stat().st_size, str(f2): f2.stat().st_size},
           "wall_seconds": round(time.perf_counter() - t0, 1)}
    mp = args.manifest or args.out_dir / f"{args.prefix}_trim_manifest.json"
    mp.write_text(json.dumps(man, indent=2) + "\n")
    print(json.dumps({k: v for k, v in man.items() if k != "parameters"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
