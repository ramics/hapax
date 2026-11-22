"""How often does each dinucleotide shuffle in this repo hand back its own input?

WHY. A dinucleotide shuffle whose walk dead-ends and then `return s` turns the null into a
copy of the test, and a null that reads like the test is how a real field gets killed.
That happened twice (holdout_run.py, fixed in PR #29; divergent_known.py, 324 of 11,361
Berlin contigs, fixed on the eswrkj branch). This measures every implementation that is
still in the tree on the same real sequence, so "affected" is a measured rate rather than
a reading of the code. Mononucleotide shuffles (`rng.shuffle` of a list) and
partner permutations have no walk and no fallback and are not measured.

Input is any FASTA of real assembled sequence, cut into read-sized windows (107 and 93 nt,
the two mate lengths of SRR23973193 after trimming) and also used whole (300-5,000 nt).
Calibrated against: `holdout_run.dinuc_shuffle`, the fixed walk, which must read 0.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from subtract_known import provenance  # noqa: E402

IMPLEMENTATIONS = ["holdout_run", "divergent_known", "residue_structure",
                   "residue_contig_annotation"]


def read_fasta(p: Path) -> list[str]:
    out, buf = [], []
    for line in open(p):
        if line.startswith(">"):
            if buf:
                out.append("".join(buf))
            buf = []
        else:
            buf.append(line.strip().upper())
    if buf:
        out.append("".join(buf))
    return [s for s in out if set(s) <= set("ACGT")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--reads", type=int, default=2000)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parents[1] / "results/fail_loudly_audit")
    a = ap.parse_args()
    seqs = read_fasta(a.fasta)
    rng = random.Random(a.seed)
    reads = []
    for s in seqs:
        for L in (107, 93):
            if len(s) > L:
                i = rng.randrange(len(s) - L)
                reads.append(s[i:i + L])
    reads = reads[:a.reads]
    whole = [s for s in seqs if 300 <= len(s) <= 5000]
    dn = lambda s: collections.Counter(zip(s, s[1:]))
    res = {}
    for mod in IMPLEMENTATIONS:
        f = importlib.import_module(mod).dinuc_shuffle
        res[mod] = {}
        for label, pool in (("read_windows", reads), ("whole_contigs", whole)):
            r = random.Random(a.seed)
            same = none = broken = 0
            for s in pool:
                o = f(s, r)
                if o is None:
                    none += 1
                elif o == s:
                    same += 1
                elif len(o) != len(s) or dn(o) != dn(s):
                    broken += 1
            res[mod][label] = {"n": len(pool), "returned_input": same,
                               "returned_input_frac": round(same / len(pool), 5),
                               "returned_none": none, "broke_dinucleotide_counts": broken}
    out = {"script": "experiments/shuffle_fallback_audit.py", **provenance(__file__),
           "seed": a.seed, "fasta": str(a.fasta),
           "fasta_sha256_16": hashlib.sha256(a.fasta.read_bytes()).hexdigest()[:16],
           "sequences": len(seqs), "results": res}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / f"shuffle_fallback_s{a.seed}.json").write_text(json.dumps(out, indent=2))
    for m, d in res.items():
        print(m, {k: v["returned_input_frac"] for k, v in d.items()})


if __name__ == "__main__":
    main()
