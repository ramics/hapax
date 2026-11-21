"""How deep in the library the residue's molecules actually are.

WHY THIS EXISTS. The residue-only assembly (`residue_assembly.py`) produced 1.31 Mb of
short contigs from 1.37 M pairs and recruited only a few per cent of them back, which says
the unplaced fraction is not a few abundant unknown genomes. That is an inference from an
assembler's behaviour. This measures the quantity directly: for each pair in one block, how
many times its 31-mers occur ELSEWHERE in the whole library, exactly, with the pair's own
occurrences subtracted (`experiments/kmer_coverage`, the repo's existing instrument).

Depth is the quantity the whole project's walls are stated in - detection dies below about
1 in a million, a genome's reads find each other only to about 1 in 100,000 - so saying
where the residue sits on that axis is saying how much room any method has.

PILES, all from the same block and the same instrument run, so no two numbers come from
different passes:
  - `identified_rrna`, `identified_other`: pairs the cascade subtracted, split by the
    library label table. The known-positive controls.
  - `residue_on_unident_contig`: residue pairs that map one of the sample's own contigs
    that nothing could identify (tier 2's read-level population).
  - `residue_unplaced`: residue pairs where neither mate matches SILVA or any own contig at
    >= 40 bases. The pile the residue-only assembly was built from.
  - `residue_on_unsplit_contig`: residue pairs whose only own-contig match is one of the
    7,944 contigs the cascade's three-way split omits. Reported separately rather than
    folded in, because it turns out to be a small, very high-depth pile - repetitive
    sequence - and averaging it into either neighbour would move that neighbour's median.

POLARITY, stated before the numbers: depth here is a description, not a novelty score. The
coverage-shape family was measured to have FAMILIARITY polarity
(`results/cov_step_conservation/`), so a high-depth read is not a more interesting read;
this field is reported as abundance and nothing else.

STATED BIAS. (a) The label table has no read-level Kraken2 or DIAMOND, so `residue_unplaced`
is a superset of the cascade's "maps nothing" pile and holds some identifiable pairs; that
pulls its depth UP, against the finding. (b) Mates are reported separately: R2 is 93 nt
after the 14 nt trim against R1's 107, so it carries 63 scorable 31-mers against 76, and a
pooled figure would average two different instruments (CLAUDE.md).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance  # noqa: E402


def summarise(med, zero) -> dict:
    med = np.asarray(med, dtype=float)
    return {"n": int(med.size),
            "median_depth_p50": float(np.median(med)),
            "median_depth_p10": float(np.quantile(med, 0.1)),
            "median_depth_p90": float(np.quantile(med, 0.9)),
            "frac_depth_0": float((med == 0).mean()),
            "frac_depth_lt_5": float((med < 5).mean()),
            "frac_depth_ge_50": float((med >= 50).mean()),
            "unwitnessed_kmer_frac_mean": float(np.mean(zero))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="which 100,000-pair block; block 1 = seed 1")
    ap.add_argument("--cov", type=Path, default=ROOT / "data/work/cov_blocks5.tsv")
    ap.add_argument("--block", type=int, default=0,
                    help="block column to read; 0 means use --seed as the block number")
    ap.add_argument("--classes", type=Path, default=ROOT / "data/work/library_label_classes.npz")
    ap.add_argument("--recut", type=Path, default=ROOT / "data/recut")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    d = np.load(args.classes)
    cls = {int(s): (int(a), int(b)) for s, a, b in zip(d["spot"], d["class1"], d["class2"])}
    residue = {int(l[1:].split("/")[0].split(".")[1])
               for l in open(args.recut / f"residue_s{args.seed}.fa") if l[0] == ">"}

    piles: dict[str, dict[str, list]] = {}
    with open(args.cov) as h:
        hdr = h.readline().split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        want = str(args.block or args.seed)
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[0] != want:
                continue
            spot = int(f[ix["name"]].split(".")[1])
            c1, c2 = cls.get(spot, (0, 0))
            if spot in residue:
                pile = "residue_unplaced" if (c1 == 0 and c2 == 0) else (
                    "residue_on_unident_contig" if 3 in (c1, c2) else "residue_on_unsplit_contig")
            elif 1 in (c1, c2):
                pile = "identified_rrna"
            else:
                pile = "identified_other"
            p = piles.setdefault(pile, {"med1": [], "med2": [], "z1": [], "z2": []})
            p["med1"].append(int(f[ix["med1"]]))
            p["med2"].append(int(f[ix["med2"]]))
            p["z1"].append(float(f[ix["zero1"]]))
            p["z2"].append(float(f[ix["zero2"]]) if f[ix["zero2"]] != "NaN" else np.nan)

    out = {"seed": args.seed, "block_pairs": sum(len(p["med1"]) for p in piles.values()),
           "piles": {k: {"mate1": summarise(v["med1"], v["z1"]),
                         "mate2": summarise(v["med2"], np.nan_to_num(v["z2"]))}
                     for k, v in sorted(piles.items())},
           "manifest": {"script": "experiments/residue_depth.py", **provenance(__file__),
                        "parameters": {k: str(v) for k, v in vars(args).items()},
                        "instrument": "experiments/kmer_coverage, exact canonical 31-mers over all "
                                      "8,731,224 pairs, each pair's own occurrences subtracted",
                        "operating_points": {"class_split": "40 matching bases, SILVA or own contig",
                                             "depth": "median over a mate's scorable 31-mers"}}}
    (args.out / f"depth_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    for k, v in out["piles"].items():
        print(f"{k:28s} n={v['mate1']['n']:6d}  m1 depth p50={v['mate1']['median_depth_p50']:7.1f} "
              f"0={v['mate1']['frac_depth_0']*100:5.2f}% <5={v['mate1']['frac_depth_lt_5']*100:5.2f}% "
              f"| m2 p50={v['mate2']['median_depth_p50']:7.1f} <5={v['mate2']['frac_depth_lt_5']*100:5.2f}%")


if __name__ == "__main__":
    main()
