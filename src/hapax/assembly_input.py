"""The two steps either side of the assembler: trim mate 2 for it, relabel what it wrote.

LIFTED, NOT IMPORTED. These were `experiments/trim_mate2.py` and
`experiments/megahit_rename.py`, run as subprocesses by `hapax.prepare`. Both scripts are
the committed producers of the files under `results/prod_SRR23702504/`, and their
`script_sha256` is what ties those files to the code that wrote them, so they stay in
`experiments/` byte-identical and this module carries the same logic for the package. That
is a deliberate duplicate, and the risk a duplicate carries here is the one CLAUDE.md
records for the shuffle: a second copy that differs from the first without anyone knowing.
So `tests/test_lifted_glue.py` runs each script and its lifted function on the same input
and asserts the outputs are byte-identical. A change to either copy fails that test, and
the change then has to be made to both or the script retired with its results.

WHY THE TRIM IS A REQUIRED ARGUMENT HERE. The script resolves it from `libraries.json` when
not given. The package resolves it once, in `hapax.prepare.resolve_trim`, which records
where it came from, so a second resolution here could only disagree with the first.

What each step does and why is documented in the script it was lifted from, and that
reasoning is not repeated here: see the two scripts' docstrings.
"""
from __future__ import annotations

import gzip
import json
import re
import time
from pathlib import Path

from .provenance import record

LEN = re.compile(r"len=(\d+)")
COV = re.compile(r"multi=([0-9.]+)")


class GlueError(ValueError):
    """An input the step refuses, with the reason the script gave."""


def _open(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def trim_mate2(r1: Path, r2: Path, trim_r2: int, out_dir: Path, prefix: str,
               manifest: Path | None = None) -> dict:
    """Write <prefix>_1.fq unchanged and <prefix>_2.fq with `trim_r2` bases cut from 5'.

    A mate 2 no longer than the prefix is written as a single N so the record counts stay
    equal, and counted in the manifest; unequal record counts are refused.
    """
    r1, r2, out_dir = Path(r1), Path(r2), Path(out_dir)
    if trim_r2 < 0:
        raise GlueError("--trim-r2 cannot be negative")
    for p in (r1, r2):
        if not p.exists():
            raise GlueError(f"{p} is missing")
    t0 = time.perf_counter()
    out_dir.mkdir(parents=True, exist_ok=True)
    f1, f2 = out_dir / f"{prefix}_1.fq", out_dir / f"{prefix}_2.fq"
    n = short = 0
    with _open(r1) as h1, _open(r2) as h2, open(f1, "w") as o1, open(f2, "w") as o2:
        while True:
            a = [h1.readline() for _ in range(4)]
            b = [h2.readline() for _ in range(4)]
            if not a[0]:
                if b[0]:
                    raise GlueError(f"{r2} has more records than {r1}")
                break
            if not b[0]:
                raise GlueError(f"{r1} has more records than {r2}")
            seq, qual = b[1].rstrip("\n"), b[3].rstrip("\n")
            if len(seq) <= trim_r2:
                short += 1
                seq, qual = "N", "!"
            else:
                seq, qual = seq[trim_r2:], qual[trim_r2:]
            o1.write(a[0] + a[1] + a[2] + a[3])
            o2.write(f"{b[0]}{seq}\n{b[2]}{qual}\n")
            n += 1
    man = {**record(__file__, "hapax.assembly_input.trim_mate2"),
           "parameters": {"r1": str(r1), "r2": str(r2), "trim_r2": trim_r2,
                          "out_dir": str(out_dir), "prefix": prefix},
           "pairs": n, "trim_r2": trim_r2,
           "mate2_shorter_than_the_prefix": short,
           "mate2_shorter_note": "written as a single N so the record count matches mate 1",
           "outputs": {str(f1): f1.stat().st_size, str(f2): f2.stat().st_size},
           "wall_seconds": round(time.perf_counter() - t0, 1)}
    mp = Path(manifest) if manifest else out_dir / f"{prefix}_trim_manifest.json"
    mp.write_text(json.dumps(man, indent=2) + "\n")
    return man


def rename_megahit(contigs: Path, out: Path, tag: str,
                   manifest: Path | None = None) -> dict:
    """Relabel MEGAHIT headers as `<tag>_NODE_<n>_length_<len>_cov_<multi>_megahit_<orig>`.

    Refuses a header with no `len=` or `multi=`, and a file with no contigs, rather than
    inventing a length or a coverage.
    """
    contigs, out = Path(contigs), Path(out)
    if not contigs.exists():
        raise GlueError(f"{contigs} is missing")
    t0 = time.perf_counter()
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(contigs) as fh, open(out, "w") as w:
        for line in fh:
            if not line.startswith(">"):
                w.write(line)
                continue
            n += 1
            L, C = LEN.search(line), COV.search(line)
            if not L or not C:
                raise GlueError(
                    f"{contigs}: contig {n} has no len= or multi= in its header "
                    f"({line.strip()!r}). That is not a MEGAHIT contig name, and renaming "
                    f"it would invent a length and a coverage.")
            orig = line[1:].split()[0]
            w.write(f">{tag}_NODE_{n}_length_{L.group(1)}_cov_{C.group(1)}"
                    f"_megahit_{orig}\n")
    if n == 0:
        raise GlueError(f"{contigs} holds no contigs")
    man = {**record(__file__, "hapax.assembly_input.rename_megahit"),
           "parameters": {"contigs": str(contigs), "out": str(out), "tag": tag},
           "contigs": n, "wall_seconds": round(time.perf_counter() - t0, 1),
           "what_this_is": "a pure relabelling: sequence is unchanged and the MEGAHIT name "
                           "is kept as a suffix"}
    mp = Path(manifest) if manifest else out.with_suffix(".rename_manifest.json")
    mp.write_text(json.dumps(man, indent=2) + "\n")
    return man
