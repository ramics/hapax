"""Are the residue-only contigs new sequence, or the sample's own contigs re-assembled?

WHY THIS EXISTS. `residue_assembly.py` built 3,695 contigs from the 1,372,021 pairs that
match neither SILVA nor any of this sample's own contigs at >= 40 bases. Those pairs were
*selected* for not matching a contig, so the contigs built from them ought to be sequence the
sample's own assembly does not hold. Ought to is not measured. A read can fail to map a
contig at 40 bases and still come from the same molecule - the assembly can stop short of it,
or the read can sit over a variant - and in that case this assembly is re-deriving the ends of
contigs that already exist and the annotation work has one pile, not two.

So this aligns the residue-only contigs against the sample's own, and reports the covered
fraction as a continuous number rather than a boolean: a contig sharing 60 of 300 bases with
an existing contig is a different statement from one sharing 290.

CONTROLS, same code, same pass:
  - `own_unidentified`: the 46,157 contigs the cascade could not identify, aligned against
    the rest of the sample's own assembly with themselves excluded by name. This is the
    known-positive for "how much does one contig of this assembly share with the others
    anyway" - without it, any overlap figure is unreadable, because contigs of one assembly
    share repeats and conserved regions with each other.
  - `shuffled`: each residue-only contig mononucleotide-shuffled. The null; a shuffle keeps
    length and GC and destroys homology, so anything it aligns is a chance match at this
    operating point.

Also reported, because it decides what the annotation can reach: the length and own-assembly
depth of the unidentifiable contigs that residue reads actually map, from the whole-library
label table. SPEC.md §3 says only 3.87% of unidentifiable contigs clear 500 nt and 2x, which
is where assembly and classification begin to work; this says how much of the residue's
read mass sits on the ones that do.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance  # noqa: E402


def read_fa(p: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name:
        out[name] = "".join(buf)
    return out


def covered_fraction(aligner, seq: str, exclude: str | None = None) -> tuple[float, str | None]:
    """Fraction of the query's bases inside any alignment, and the best target.

    Union of alignment intervals rather than the single best hit: a contig re-derived in two
    pieces, or one matching a repeat twice, is still covered, and the single-best-hit reading
    would call it new.
    """
    mask = np.zeros(len(seq), dtype=bool)
    best, which = 0, None
    for h in aligner.map(seq):
        if exclude is not None and h.ctg == exclude:
            continue
        mask[h.q_st:h.q_en] = True
        if h.mlen > best:
            best, which = h.mlen, h.ctg
    return float(mask.mean()), which


def shuffle(s: str, rng: random.Random) -> str:
    x = list(s)
    rng.shuffle(x)
    return "".join(x)


def summarise(v) -> dict:
    a = np.asarray(v, dtype=float)
    return {"n": int(a.size), "mean": float(a.mean()), "p50": float(np.median(a)),
            "frac_ge_0.5": float((a >= 0.5).mean()), "frac_ge_0.9": float((a >= 0.9).mean()),
            "frac_lt_0.1": float((a < 0.1).mean())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--resasm", type=Path, default=ROOT / "data/work/resasm/final.contigs.fa")
    ap.add_argument("--own-all", type=Path, default=ROOT / "data/work/sample_contigs_all.fa")
    ap.add_argument("--own-unident", type=Path, default=ROOT / "data/geo/sets/sample_contigs_unidentified.fa")
    ap.add_argument("--labels", type=Path, default=ROOT / "data/work/library_labels.npz")
    ap.add_argument("--classes", type=Path, default=ROOT / "data/work/library_label_classes.npz")
    ap.add_argument("--control-n", type=int, default=2000, help="own-unidentified contigs to score")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    ap.add_argument("--write-split", action="store_true",
                    help="write resasm_{new,partial,rederived}.fa and covered_fraction.tsv beside the assembly")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    import mappy

    res = read_fa(args.resasm)
    own_unident = read_fa(args.own_unident)
    al = mappy.Aligner(str(args.own_all), preset="sr")
    if not al:
        raise SystemExit("could not index the sample's own assembly")

    out: dict = {"seed": args.seed, "operating_point": "minimap2 sr, union of alignment intervals"}

    cov, hits, per = [], Counter(), {}
    for k, s in res.items():
        f, w = covered_fraction(al, s)
        cov.append(f)
        hits[w] += 1
        per[k] = f
    # Split the contigs by that fraction, because the annotation downstream wants the pile
    # the sample's own assembly does NOT already hold. The thresholds are an operating point
    # and the continuous values travel beside them in covered_fraction.tsv.
    if args.write_split:
        (args.resasm.parent / "covered_fraction.tsv").write_text(
            "contig\tcovered_fraction_of_own_assembly\tlen\n"
            + "".join(f"{k}\t{v:.4f}\t{len(res[k])}\n" for k, v in per.items()))
        for tag, pick in (("new", lambda f: f < 0.1), ("partial", lambda f: 0.1 <= f < 0.5),
                          ("rederived", lambda f: f >= 0.5)):
            sel = [k for k, f in per.items() if pick(f)]
            (args.resasm.parent / f"resasm_{tag}.fa").write_text(
                "".join(f">{k} covered={per[k]:.4f}\n{res[k]}\n" for k in sel))
            out[f"split_{tag}"] = {"contigs": len(sel), "bases": sum(len(res[k]) for k in sel)}
    out["residue_contigs_vs_own_assembly"] = summarise(cov)
    out["residue_contigs"] = {"n": len(res), "bases": sum(len(s) for s in res.values())}

    sh = [covered_fraction(al, shuffle(s, rng))[0] for s in res.values()]
    out["shuffled_null"] = summarise(sh)

    # Control: an own-unidentified contig against the rest of the same assembly, itself out.
    keys = sorted(own_unident)
    rng.shuffle(keys)
    ctrl = [covered_fraction(al, own_unident[k], exclude=k)[0] for k in keys[:args.control_n]]
    out["own_unidentified_vs_rest_of_own_assembly"] = summarise(ctrl)

    # Where the residue's read mass sits among unidentifiable contigs, by contig length and
    # the coverage SPAdes recorded in the contig name.
    d = np.load(args.labels)
    cl = np.load(args.classes)
    names = d["contig_names"]
    unid = {n.split()[0] for n in own_unident}
    resid_mass, lens, covs, distinct = 0, [], [], set()
    sel = (cl["class1"] == 3) | (cl["class2"] == 3)
    for m in (1, 2):
        ctg, cm = d[f"ctg{m}"], d[f"cm{m}"]
        ok = sel & (cm >= 40) & (ctg >= 0)
        for i in np.flatnonzero(ok):
            n = str(names[ctg[i]])
            if n not in unid:
                continue
            resid_mass += 1
            distinct.add(n)
            lens.append(int(n.split("_length_")[1].split("_")[0]))
            covs.append(float(n.split("_cov_")[1]))
    lens, covs = np.array(lens), np.array(covs)
    workable = (lens >= 500) & (covs >= 2)
    out["residue_mass_on_unidentifiable_contigs"] = {
        "mate_observations": int(resid_mass),
        "contig_len_p50": float(np.median(lens)), "contig_cov_p50": float(np.median(covs)),
        "frac_on_contigs_ge500nt_and_2x": float(workable.mean()),
        "frac_on_contigs_ge500nt": float((lens >= 500).mean()),
        "frac_on_contigs_ge1kb": float((lens >= 1000).mean()),
        "distinct_contigs": len(distinct),
    }
    out["manifest"] = {"script": "experiments/residue_contig_overlap.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()}}
    (args.out / f"contig_overlap_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    for k in ("residue_contigs_vs_own_assembly", "own_unidentified_vs_rest_of_own_assembly", "shuffled_null"):
        v = out[k]
        print(f"{k:42s} n={v['n']:5d} mean={v['mean']:.3f} p50={v['p50']:.3f} "
              f">=0.5 {v['frac_ge_0.5']*100:5.1f}% >=0.9 {v['frac_ge_0.9']*100:5.1f}% <0.1 {v['frac_lt_0.1']*100:5.1f}%")
    print(json.dumps(out["residue_mass_on_unidentifiable_contigs"], indent=1))


if __name__ == "__main__":
    main()
