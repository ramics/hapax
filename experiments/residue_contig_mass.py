"""How the residue's read mass is spread over the contigs it sits on, and how few carry it.

WHY THIS EXISTS. The residue splits into two halves an order of magnitude apart in depth
(`results/residue_census/` §3): about 8.4 points of the library sit on contigs the cascade
could not identify, at median depth 6, and about 11 points are unplaced at depth 1. The
depth-6 half is the one contig-level work can reach, and the question that decides what to do
with it is not how many such contigs exist but **how concentrated the read mass is**. 46,157
contigs is a count nobody can act on; if half the mass sits on a few hundred, that is a list a
person can read, which is what Imo's recut of the target asks for - the eyes of a skilled
bioinformatician over a shortlist, not a per-read verdict.

So this reports the Lorenz curve of residue mass over contigs: how many contigs carry 50% and
90% of it, and what those contigs look like (length, the coverage SPAdes recorded, residue
mass per kb).

CONTROL, same code, same pass: the identified non-rRNA contigs with the mass of reads the
cascade identified. Concentration is a property of any metagenome - abundance is heavy-tailed
whatever the sequence is - so a concentration figure without that reference says nothing. If
the residue's mass is no more concentrated than the identified mass, then "a readable
shortlist" is just what a wastewater library looks like and is not a property of the residue.

STATED BIAS. Read mass here is mate observations at >= 40 matching bases from the whole-library
label table, which holds no read-level Kraken2 or DIAMOND, so a contig's residue mass includes
mates those indexes would have identified; that inflates mass on contigs resembling known
material, i.e. against a finding that the residue concentrates somewhere interesting. Mates
are counted separately and reported per mate as well as pooled, since R2 is 93 nt against
R1's 107 and maps at a different rate (CLAUDE.md).
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


def lorenz(mass: np.ndarray) -> dict:
    """Contigs needed for 50% and 90% of the mass, and the Gini coefficient."""
    m = np.sort(mass)[::-1]
    tot = m.sum()
    c = np.cumsum(m)
    out = {"contigs_with_mass": int((m > 0).sum()), "total_mass": int(tot)}
    for f in (0.5, 0.9):
        out[f"contigs_for_{int(f * 100)}pct"] = int(np.searchsorted(c, f * tot) + 1)
    asc = np.sort(m[m > 0])
    n = asc.size
    out["gini"] = float((2 * np.arange(1, n + 1) - n - 1).dot(asc) / (n * asc.sum()))
    out["top10_share"] = float(m[:10].sum() / tot)
    out["top100_share"] = float(m[:100].sum() / tot)
    return out


def contig_len(name: str) -> int:
    return int(name.split("_length_")[1].split("_")[0])


def contig_cov(name: str) -> float:
    return float(name.split("_cov_")[1])


def describe(names, mass, k: int = 200) -> dict:
    """Length, SPAdes coverage and mass per kb of the k highest-mass contigs, against all."""
    order = np.argsort(mass)[::-1]
    top = [i for i in order[:k] if mass[i] > 0]
    allm = [i for i in range(len(mass)) if mass[i] > 0]

    def block(ix):
        L = np.array([contig_len(str(names[i])) for i in ix], dtype=float)
        C = np.array([contig_cov(str(names[i])) for i in ix], dtype=float)
        M = np.array([mass[i] for i in ix], dtype=float)
        return {"n": len(ix), "len_p50": float(np.median(L)), "cov_p50": float(np.median(C)),
                "mass_p50": float(np.median(M)), "mass_per_kb_p50": float(np.median(M / (L / 1000))),
                "frac_ge500_and_2x": float(np.mean((L >= 500) & (C >= 2)))}

    return {f"top{k}": block(top), "all_with_mass": block(allm)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded for the convention. THIS PASS IS DETERMINISTIC: it counts "
                         "every mate of all 8,731,224 pairs, so there is no sampling and no "
                         "spread to report. Five identical seeds would be five copies of one "
                         "number, which is not a five-seed result - the census is the whole "
                         "population, and that is stated rather than dressed up.")
    ap.add_argument("--labels", type=Path, default=ROOT / "data/work/library_labels.npz")
    ap.add_argument("--classes", type=Path, default=ROOT / "data/work/library_label_classes.npz")
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    d = np.load(args.labels)
    cl = np.load(args.classes)
    names = d["contig_names"]
    nc = len(names)

    # Mass per contig, per mate, for the two piles. A mate contributes to the contig it best
    # matches at >= 40 bases; class 3 is an unidentifiable contig, class 2 an identified one.
    piles = {"residue_on_unidentifiable": 3, "identified_contigs_control": 2}
    out: dict = {"seed": args.seed, "contigs_in_assembly": nc,
                 "operating_point": "mate's best own-contig match >= 40 matching bases"}
    for pile, want in piles.items():
        per_mate = {}
        pooled = np.zeros(nc, dtype=np.int64)
        for m in (1, 2):
            ctg, cm, c = d[f"ctg{m}"], d[f"cm{m}"], cl[f"class{m}"]
            ok = (c == want) & (cm >= 40) & (ctg >= 0)
            mass = np.bincount(ctg[ok], minlength=nc).astype(np.int64)
            per_mate[f"mate{m}"] = lorenz(mass)
            pooled += mass
        out[pile] = {**per_mate, "pooled": lorenz(pooled), "shape": describe(names, pooled, args.top)}
        # The shortlist itself, so the number is actionable rather than quoted.
        if want == 3:
            order = np.argsort(pooled)[::-1][:args.top]
            (args.out / "top_unidentifiable_contigs.tsv").write_text(
                "contig\tresidue_mate_observations\tlen\tspades_cov\n"
                + "".join(f"{names[i]}\t{pooled[i]}\t{contig_len(str(names[i]))}\t"
                          f"{contig_cov(str(names[i]))}\n" for i in order if pooled[i] > 0))

    out["manifest"] = {"script": "experiments/residue_contig_mass.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()}}
    (args.out / f"contig_mass_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    for pile in piles:
        p = out[pile]["pooled"]
        s = out[pile]["shape"]
        print(f"{pile:32s} contigs={p['contigs_with_mass']:6d} mass={p['total_mass']:9d} "
              f"50%={p['contigs_for_50pct']:5d} 90%={p['contigs_for_90pct']:6d} gini={p['gini']:.3f} "
              f"top100={p['top100_share']*100:5.1f}%  top{args.top} len={s[f'top{args.top}']['len_p50']:.0f} "
              f"cov={s[f'top{args.top}']['cov_p50']:.1f}")


if __name__ == "__main__":
    main()
