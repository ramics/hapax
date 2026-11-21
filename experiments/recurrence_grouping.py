"""Can reads be grouped by WHERE they occur instead of by what they are?

THE IDEA, from `/mnt/project-files/single_molecule/BRAINSTORM.md` §2. Two reads at depth 1
can never be grouped by sequence: they share no overlap, no neighbourhood and no second
observation, which is the floor `results/residue_structure/` §4 and `results/mate_friends/`
measured. But if each read is present in the same subset of the plant's other samples, that
pattern is evidence they came from the same organism - the phylogenetic-profile trick from
comparative genomics, applied to reads. It needs no second copy in this run, which is why it
is the one grouping idea the sampling wall does not already answer.

THE POSITIVE CONTROL IS THE WHOLE EXPERIMENT, AND THIS LIBRARY SUPPLIES IT FREE. A read's
own mate came off the same molecule. So the question "can co-occurrence group reads of one
organism" has a label-free ground truth: rank every pair's mate 2 by how similar its
occupancy vector is to a given pair's mate 1, and ask whether that pair's OWN mate is near
the top. That is exactly the instrument `results/mate_friends/` used for sequence
neighbourhoods, reused on a different space, so the two walls are measured the same way and
the numbers are comparable.

Chance is the first thing reported, not the last: with n candidates and a top-k list, chance
is k/n, and any multiplier is read against that.

THE NULL FOR STRUCTURE IS A CURVEBALL SWAP, NOT A SHUFFLE. Occupancy vectors have two
marginals that make reads look alike for uninteresting reasons: an abundant read is present
nearly everywhere, and a deeply sequenced subject witnesses nearly everything. A null that
destroys both would make any observed similarity look like structure. So the null preserves
BOTH - each read's number of samples and each sample's number of reads - by repeated
trial swaps of 2x2 submatrices (Strona et al.'s curveball algorithm), which is the standard
null for presence-absence matrices in ecology for exactly this reason.

WHAT WOULD KILL IT, stated before the run, per the brainstorm: occupancy vectors with too
little entropy to separate anything, i.e. reads that are in all samples or in one. The
entropy is reported first for that reason, and if the informative middle band is empty, the
grouping question is answered without running the grouping.

POLARITY. Co-occurrence similarity is a FAMILIARITY-FLAVOURED field like everything else
here: reads that recur together recur, and recurrence is the opposite of novel. What this
measures is not novelty but whether the unplaceable pile has internal structure - whether it
is many organisms or one, and which of its reads belong together. That is the grouping
question SPEC.md defers, not a score.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from subtract_known import provenance  # noqa: E402


def load(run_dir: Path, panel: Path, group_by: str) -> tuple[list[str], dict, np.ndarray]:
    """Occupancy as (reads x groups) booleans per mate, at the 'any shared 31-mer' point."""
    names = (run_dir / "query_names.txt").read_text().split()
    rows = {r["run"]: r for r in csv.DictReader(open(panel), delimiter="\t")}
    files = sorted(run_dir.glob("witness_*.npz"))
    if not files:
        raise SystemExit(f"no witness_*.npz in {run_dir}")
    runs = [f.name[len("witness_"):-len(".npz")] for f in files]
    buckets: dict[str, list[int]] = defaultdict(list)
    for i, run in enumerate(runs):
        buckets[rows[run][group_by] if group_by != "run" else run].append(i)
    keys = sorted(buckets)

    present = {1: np.zeros((len(keys), len(names)), dtype=bool),
               2: np.zeros((len(keys), len(names)), dtype=bool)}
    raw = {1: np.zeros((len(files), len(names)), dtype=bool),
           2: np.zeros((len(files), len(names)), dtype=bool)}
    for i, f in enumerate(files):
        d = np.load(f)
        for mate in (1, 2):
            raw[mate][i] = np.nan_to_num(d[f"frac{mate}"], nan=0.0) > 0
    for j, k in enumerate(keys):
        for mate in (1, 2):
            present[mate][j] = raw[mate][buckets[k]].any(axis=0)
    return names, present, np.array(keys)


def curveball(mat: np.ndarray, rng: np.random.Generator, rounds: int = 8) -> np.ndarray:
    """A presence-absence matrix with both marginals preserved (Strona et al. 2014).

    Rows are reads, columns samples. Each trial takes two rows, keeps what they share, and
    deals the rest between them at random, which cannot change either row's total and cannot
    change any column's total either.
    """
    m = mat.copy()
    n = m.shape[0]
    sets = [set(np.flatnonzero(r).tolist()) for r in m]
    for _ in range(rounds * n):
        i, j = rng.integers(0, n, size=2)
        if i == j:
            continue
        a, b = sets[i], sets[j]
        shared = a & b
        swappable = list((a | b) - shared)
        if len(swappable) < 2:
            continue
        rng.shuffle(swappable)
        take = len(a) - len(shared)
        sets[i] = shared | set(swappable[:take])
        sets[j] = shared | set(swappable[take:])
    out = np.zeros_like(m)
    for i, s in enumerate(sets):
        if s:
            out[i, list(s)] = True
    return out


def jaccard_topk(query: np.ndarray, candidates: np.ndarray, k: int) -> np.ndarray:
    """For each query row, the indices of the k most similar candidate rows (Jaccard).

    Ties are a real hazard here - occupancy vectors are low-dimensional and many reads share
    one - and `results/residue_structure/` lost a result to ties broken by array order. So
    ties are broken at random by adding a tiny random key before the sort, which makes the
    chance rate the honest one rather than an artefact of how the rows were stacked.
    """
    q = query.astype(np.float32)
    c = candidates.astype(np.float32)
    inter = q @ c.T
    qs = q.sum(axis=1, keepdims=True)
    cs = c.sum(axis=1, keepdims=True).T
    union = qs + cs - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        sim = np.where(union > 0, inter / np.maximum(union, 1e-9), 0.0)
    jitter = np.random.default_rng(0).random(sim.shape).astype(np.float32) * 1e-6
    return np.argsort(-(sim + jitter), axis=1)[:, :k], sim


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--panel", type=Path, default=ROOT / "results/recurrence/panel.tsv")
    ap.add_argument("--group-by", default="water", choices=("run", "water", "date"),
                    help="the occupancy vector's dimensions. `water` by default, because "
                         "preparation variants of one water are not independent samples")
    ap.add_argument("--piles", nargs="+",
                    default=["unplaced", "unident_contig", "identified"])
    ap.add_argument("--k", type=int, default=10,
                    help="top-k for the mate-finding control, as results/mate_friends/ used")
    ap.add_argument("--max-reads", type=int, default=4000,
                    help="per pile per seed. The similarity matrix is quadratic, and the "
                         "control's chance rate depends on this, so it is recorded with "
                         "every number rather than tuned")
    ap.add_argument("--out", type=Path, default=ROOT / "results/recurrence")
    ap.add_argument("--label", default="grouping")
    args = ap.parse_args()

    t0 = time.perf_counter()
    names, present, keys = load(args.run_dir, args.panel, args.group_by)
    pile = np.array([x.split("|")[0] for x in names])
    dbin = np.array([x.split("|")[1] for x in names])
    seed_of = np.array([x.split("|")[4] if len(x.split("|")) > 4 else "s1" for x in names])
    rng = np.random.default_rng(args.seed)

    out: dict = {"seed": args.seed, "group_by": args.group_by,
                 "dimensions": int(keys.size), "piles": {}}

    for p in args.piles:
        for seed_tag in sorted(set(seed_of)):
            sel = np.flatnonzero((pile == p) & (seed_of == seed_tag))
            if sel.size < 200:
                continue
            if sel.size > args.max_reads:
                sel = rng.choice(sel, args.max_reads, replace=False)
            m1 = present[1][:, sel].T            # reads x groups
            m2 = present[2][:, sel].T

            # ENTROPY FIRST. If the vectors are all-in or all-out there is nothing to group
            # by, and the brainstorm's own kill criterion applies before any grouping runs.
            occ = m1.sum(axis=1)
            d = keys.size
            informative = (occ > 0) & (occ < d)
            frac_mid = float(informative.mean())
            bits = 0.0
            if informative.any():
                q = m1[informative].mean(axis=0)
                q = np.clip(q, 1e-9, 1 - 1e-9)
                bits = float(-(q * np.log2(q) + (1 - q) * np.log2(1 - q)).sum())

            # THE POSITIVE CONTROL: can a read's own mate be found by occupancy alone?
            ok = (m1.sum(axis=1) > 0) & (m2.sum(axis=1) > 0)
            res_control: dict = {"usable_pairs": int(ok.sum())}
            if ok.sum() >= 200:
                q, c = m1[ok], m2[ok]
                k = min(args.k, c.shape[0] - 1)
                top, _ = jaccard_topk(q, c, k)
                own = np.arange(ok.sum())[:, None]
                hit = float((top == own).any(axis=1).mean())
                chance = k / ok.sum()
                # The null keeps both marginals, so a multiplier over it is structure and
                # not the two things that make vectors look alike for free.
                cnull = curveball(c, rng)
                topn, _ = jaccard_topk(q, cnull, k)
                hit_null = float((topn == own).any(axis=1).mean())
                res_control.update({
                    "k": k, "own_mate_in_topk": hit, "chance": chance,
                    "own_mate_in_topk_curveball_null": hit_null,
                    "multiplier_over_chance": hit / chance if chance else float("nan"),
                    "multiplier_over_null": hit / hit_null if hit_null else float("inf"),
                })

            # STRUCTURE: how similar is a read to its ten nearest neighbours, against the
            # marginal-preserving null on the same reads.
            res_struct: dict = {}
            if informative.sum() >= 200:
                mm = m1[informative]
                k = min(args.k, mm.shape[0] - 1)
                _, sim = jaccard_topk(mm, mm, k)
                np.fill_diagonal(sim, -1.0)
                knn = -np.sort(-sim, axis=1)[:, :k]
                null = curveball(mm, rng)
                _, nsim = jaccard_topk(null, null, k)
                np.fill_diagonal(nsim, -1.0)
                nknn = -np.sort(-nsim, axis=1)[:, :k]
                res_struct = {"n": int(informative.sum()), "k": k,
                              "mean_knn_jaccard": float(knn.mean()),
                              "mean_knn_jaccard_null": float(nknn.mean()),
                              "excess": float(knn.mean() - nknn.mean())}

            out["piles"][f"{p}|{seed_tag}"] = {
                "n": int(sel.size),
                "occupancy_mean": float(occ.mean()),
                "frac_vectors_informative": frac_mid,
                "entropy_bits_of_informative_vectors": bits,
                "max_possible_bits": float(d),
                "mate_control": res_control,
                "structure": res_struct,
            }

    out["manifest"] = {"script": "experiments/recurrence_grouping.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()},
                       "operating_points": {"presence": "at least one shared 31-mer",
                                            "k": 31, "topk": args.k},
                       "seconds": round(time.perf_counter() - t0, 1)}
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / f"{args.label}_s{args.seed}.json"
    path.write_text(json.dumps(out, indent=1) + "\n")

    print(f"{keys.size} {args.group_by} dimensions")
    for key, v in out["piles"].items():
        c, s = v["mate_control"], v["structure"]
        print(f"{key:26s} n={v['n']:5d} informative={v['frac_vectors_informative']*100:5.1f}% "
              f"mate@{c.get('k', 0)}={c.get('own_mate_in_topk', float('nan')):.4f} "
              f"(chance {c.get('chance', float('nan')):.4f}, "
              f"null {c.get('own_mate_in_topk_curveball_null', float('nan')):.4f})  "
              f"knn J={s.get('mean_knn_jaccard', float('nan')):.3f} vs null "
              f"{s.get('mean_knn_jaccard_null', float('nan')):.3f}")
    print(f"\n{path}")


if __name__ == "__main__":
    main()
