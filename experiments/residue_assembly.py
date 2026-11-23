"""Assemble the library's unplaced pairs on their own.

WHY THIS EXISTS. About a third of the residue has no other read covering it, and most of
the rest sits among reads nobody can place. A 107 nt read is too short for a profile search
(Pfam, Rfam) to have any power, so the composition question - is this coding, is it a known
RNA family, is it a virus-like protein - can only be asked of longer sequence. The authors'
SPAdes assembly was built from the whole library, where 38% rRNA and a few very abundant
transcripts dominate the graph. Assembling only what that assembly and SILVA could not place
removes that crowding, so pieces too thin or too tangled to come out before may come out now.

INPUT. Every library pair whose two mates both read class 0 ("none") in
`library_label_classes.py`: no SILVA match and no own-contig match at >= 40 matching bases.
STATED BIAS: those labels have no read-level Kraken2 or DIAMOND, so this set is a superset of
the cascade's "maps nothing" residue - it also holds pairs Kraken2 or DIAMOND would identify.
That inflates "known" calls on these contigs rather than hiding them, so it errs against a
finding of novelty.

OUTPUT. `data/work/resasm/final.contigs.fa` plus a manifest. Annotation is a separate step
(`residue_contig_annotate.py`), with the sample's own identified and unidentified contigs run
through the same code as controls.
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

import library_facts

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance, _peak_rss  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="convention; MEGAHIT is run with 1 thread-independent setting")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--classes", type=Path, default=ROOT / "data/work/library_label_classes.npz")
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--megahit", default="/home/user/tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=ROOT / "data/work")
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["trim_r2"])
    t0 = time.perf_counter()
    d = np.load(args.classes)
    keep = set(d["spot"][(d["class1"] == 0) & (d["class2"] == 0)].tolist())
    f1, f2 = args.work / "unplaced_1.fq", args.work / "unplaced_2.fq"
    n = 0
    with gzip.open(args.r1, "rt") as h1, gzip.open(args.r2, "rt") as h2, open(f1, "w") as o1, open(f2, "w") as o2:
        while True:
            a = [h1.readline() for _ in range(4)]
            b = [h2.readline() for _ in range(4)]
            if not a[0]:
                break
            spot = int(a[0][1:].split()[0].split(".")[1])
            if spot in keep:
                o1.write(f"@{spot}/1\n{a[1]}+\n{a[3]}")
                o2.write(f"@{spot}/2\n{b[1][args.trim_r2:]}+\n{b[3][args.trim_r2:]}")
                n += 1
    t_extract = time.perf_counter() - t0
    out = args.work / "resasm"
    t1 = time.perf_counter()
    subprocess.run([args.megahit, "-1", str(f1), "-2", str(f2), "-o", str(out), "-t", str(args.threads),
                    "--presets", "meta-sensitive", "--min-contig-len", "200"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t_asm = time.perf_counter() - t1
    lens = []
    for line in open(out / "final.contigs.fa"):
        if line[0] == ">":
            lens.append(int(line.split("len=")[1].split()[0]))
    lens = np.array(lens)
    man = {"script": "experiments/residue_assembly.py", **provenance(__file__),
           "parameters": {k: str(v) for k, v in vars(args).items()},
           "pairs_in": n, "pairs_in_frac_of_library": n / len(d["spot"]),
           "contigs": int(lens.size), "bases": int(lens.sum()),
           "len_p50": float(np.median(lens)), "ge500": int((lens >= 500).sum()), "ge1000": int((lens >= 1000).sum()),
           "max": int(lens.max()), "seconds_extract": t_extract, "seconds_assemble": t_asm,
           "peak_rss_bytes": _peak_rss(),
           "megahit": subprocess.run([args.megahit, "--version"], capture_output=True, text=True).stdout.strip()}
    (out / "manifest.json").write_text(json.dumps(man, indent=2))
    print(json.dumps({k: v for k, v in man.items() if k not in ("parameters",)}, indent=1))


if __name__ == "__main__":
    main()
