#!/usr/bin/env python3
"""Relabel MEGAHIT contig names into the SPAdes-style form the committed scripts parse.

`classify_contigs.py`, `subtract_known.py` and `residue_report.py` read a contig's length
and coverage out of its NAME, because that is what the authors' SPAdes assembly carried
(`NODE_1_length_12710_cov_727.4`). MEGAHIT writes the same two numbers in its own words
(`>k141_12 flag=1 multi=3.7261 len=486`). Rewriting the name is a pure relabelling of the
same contig and keeps every downstream script unmodified; the alternative is a second
header parser, which is a second instrument reading the same quantity.

WHY A SCRIPT AND NOT A HEREDOC. This is the third place in the repo that needs it:
`holdout_cascade_run.sh` and `path_settings_run.sh` each carry their own copy inline. Those
two are grandfathered - editing them would change the provenance of results nobody is
re-running (experiments/SHARED.md's rule) - and nothing new should re-type it. The same
argument as `shuffle_control.py`: the second writer of a copy does not know which condition
the first was fixing.

IT REFUSES rather than guessing. A header with no `len=` or no `multi=` is not a MEGAHIT
contig, and a file of them would otherwise come out named `length_None`, which parses to
nothing and reads downstream as a contig of unknown length. The original name is kept as a
suffix so a contig can be traced back to the assembly it came from.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

from provenance import provenance

LEN = re.compile(r"len=(\d+)")
COV = re.compile(r"multi=([0-9.]+)")


def rename(src: Path, dst: Path, tag: str) -> dict:
    n = 0
    with open(src) as fh, open(dst, "w") as w:
        for line in fh:
            if not line.startswith(">"):
                w.write(line)
                continue
            n += 1
            L, C = LEN.search(line), COV.search(line)
            if not L or not C:
                raise SystemExit(
                    f"{src}: contig {n} has no len= or multi= in its header "
                    f"({line.strip()!r}). That is not a MEGAHIT contig name, and renaming "
                    f"it would invent a length and a coverage.")
            orig = line[1:].split()[0]
            w.write(f">{tag}_NODE_{n}_length_{L.group(1)}_cov_{C.group(1)}"
                    f"_megahit_{orig}\n")
    if n == 0:
        raise SystemExit(f"{src} holds no contigs")
    return {"contigs": n}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--contigs", type=Path, required=True,
                   help="MEGAHIT final.contigs.fa")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tag", required=True,
                   help="prefix for every contig name, normally the run accession, so a "
                        "contig carries which library it was assembled from")
    ap.add_argument("--manifest", type=Path, default=None)
    args = ap.parse_args()
    if not args.contigs.exists():
        raise SystemExit(f"{args.contigs} is missing")
    t0 = time.perf_counter()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stats = rename(args.contigs, args.out, args.tag)
    man = {"script": "experiments/megahit_rename.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           **stats, "wall_seconds": round(time.perf_counter() - t0, 1),
           "what_this_is": "a pure relabelling: sequence is unchanged and the MEGAHIT name "
                           "is kept as a suffix"}
    mp = args.manifest or args.out.with_suffix(".rename_manifest.json")
    mp.write_text(json.dumps(man, indent=2) + "\n")
    print(f"{stats['contigs']} contigs relabelled -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
