"""Structured or diffuse: the three readouts, on the joint table built by
experiments/residue_structure.py and experiments/kmer_coverage.

The axes are listed once, in AXES, and every readout uses the same list, so no readout can
quietly use a different space from the figure. Each is z-scored on the residue population
of the seed being analysed, because the question is about the shape of the residue cloud and
the reference clouds are drawn into the same frame.

1. SEPARABILITY. Logistic regression, residue against subtracted, five-fold cross-validated
   AUC, compared against the best single axis's AUC on the same folds. This is the minimum
   claim of SPEC.md §4.2: if the joint space does not beat the best threshold, the
   multidimensional framing buys nothing and should be reported as buying nothing. The
   single-axis arm is the honest comparator precisely because every killed feature in this
   project was a single-axis threshold.

2. CLUSTERING BEYOND THE COVARIANCE. k-means for k = 1..8, with log within-cluster
   dispersion compared against a parametric bootstrap drawn with the same mean and
   covariance as the residue cloud (the gap statistic of Tibshirani et al. 2001, with the
   multivariate-normal reference rather than a uniform box). The normal reference is the
   point: a dozen correlated fields make apparent clusters when they are one measurement
   wearing a hat, and the covariance-matched null removes exactly that illusion. Only
   clustering the covariance cannot produce counts as structure.

3. LABEL RECOVERY. Three labels the clustering never sees: the planted arms (known answer),
   the residue piles the cascade defines (maps the sample's own contigs or nothing at all),
   and the mate-disagreement extremes. Reported as the enrichment of each label in the
   cluster that holds most of it, against its background share.

Every readout runs per seed and is reported with spread across the five seeds, per the
repo's five-seed rule, and the per-seed values are written out so the spread is auditable
rather than asserted.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
MIN_KMERS = 40   # a pair below this has no within-read coverage shape to measure

# (axis name, how it is built from the coverage and composition columns). Coverage axes are
# logged because library depth spans five orders of magnitude; the mate and step axes are
# already ratios in log space.
AXES = [
    "log_cov",          # pair's median witnesses elsewhere in the library
    "cov_spread",       # within-read log2 spread of coverage, q90 over q10
    "zero_frac",        # fraction of k-mers with no witness anywhere else
    "mate_disagree",    # |log2| ratio of the two mates' medians
    "cov_step",         # largest within-read coverage step, either mate
    "gc",
    "entropy3",
    "CpG_oe",
    "UpA_oe",
    "homopolymer",
    "mean_error_prob",
]


def load(work: Path) -> dict:
    comp: dict[str, list[float]] = {}
    with open(work / "composition.tsv") as h:
        head = h.readline().rstrip("\n").split("\t")[1:]
        for line in h:
            f = line.rstrip("\n").split("\t")
            comp[f[0]] = [float(x) for x in f[1:]]
    ci = {k: i for i, k in enumerate(head)}

    # A pair with almost no scorable k-mers (a short mate, or N-rich sequence) has no
    # coverage shape to measure, so it is dropped rather than imputed and the count is
    # reported: 408 pairs in 529,700 carry no 31-mer at all, 1,422 carry fewer than 40.
    dropped = 0
    rows: list[dict] = []
    with open(work / "coverage.tsv") as h:
        head2 = h.readline().rstrip("\n").split("\t")
        for line in h:
            f = line.rstrip("\n").split("\t")
            r = dict(zip(head2, f))
            name = r["name"]
            c = comp.get(name)
            if c is None:
                continue
            m1, m2 = float(r["med1"]), float(r["med2"])
            n1, n2 = float(r["n1"]), float(r["n2"])
            if n1 + n2 < MIN_KMERS:
                dropped += 1
                continue
            z1, z2 = float(r["zero1"]), float(r["zero2"])
            q10 = float(r["q10_1"]) + float(r["q10_2"])
            q90 = float(r["q90_1"]) + float(r["q90_2"])
            row = {
                "block": int(r["block"]),
                "name": name,
                "log_cov": math.log1p((m1 + m2) / 2),
                "cov_spread": math.log2((q90 + 1) / (q10 + 1)),
                "zero_frac": (z1 * n1 + z2 * n2) / max(n1 + n2, 1),
                "mate_disagree": abs(math.log2((m1 + 1) / (m2 + 1))),
                "cov_step": max(float(r["step1"]), float(r["step2"])),
                **{k: c[ci[k]] for k in ("gc", "entropy3", "CpG_oe", "UpA_oe",
                                         "homopolymer")},
                "mean_error_prob": math.log10(c[ci["mean_error_prob"]] + 1e-6),
            }
            if all(math.isfinite(row[a]) for a in AXES):
                rows.append(row)
            else:
                dropped += 1
    return rows, dropped


def matrix(rows: list[dict]) -> np.ndarray:
    return np.array([[r[a] for a in AXES] for r in rows], dtype=float)


def gap_statistic(X: np.ndarray, ks: list[int], rng: np.random.Generator,
                  b: int = 10) -> list[dict]:
    """Tibshirani gap with a covariance-matched multivariate-normal reference.

    The reference is drawn with the data's own mean and covariance, so a gap above zero
    means clustering that the correlation structure alone does not produce. That is the
    distinction the coordinator's brief asked for: independence and structure are different
    properties, and a dozen fields tracking one thing must not read as a dozen clusters.
    """
    mu, cov = X.mean(0), np.cov(X, rowvar=False)
    refs = [rng.multivariate_normal(mu, cov, size=X.shape[0]) for _ in range(b)]

    def wk(A: np.ndarray, k: int) -> float:
        if k == 1:
            return math.log(((A - A.mean(0)) ** 2).sum() + 1e-12)
        km = KMeans(n_clusters=k, n_init=4, random_state=0).fit(A)
        return math.log(km.inertia_ + 1e-12)

    out = []
    for k in ks:
        w = wk(X, k)
        rs = [wk(R, k) for R in refs]
        sk = float(np.std(rs)) * math.sqrt(1 + 1 / b)
        out.append({"k": k, "log_W": w, "ref_log_W": float(np.mean(rs)),
                    "gap": float(np.mean(rs) - w), "s_k": sk})
    return out


COVERAGE_AXES = ["log_cov", "cov_spread", "zero_frac", "mate_disagree", "cov_step"]
COMPOSITION_AXES = ["gc", "entropy3", "CpG_oe", "UpA_oe", "homopolymer"]
# The simulator writes no quality string, so every planted and shuffled read carries the
# same mean error probability. That constant would hand the neighbourhood test a free label
# identifying planted reads, so the planted readouts run in spaces that exclude it. Stated
# here rather than fixed by inventing a quality model this library's profile does not
# justify (flat, NovaSeq-binned: results/berlin_quality).
ARMS_TAGS = [("spike_1e-3", 1e-3), ("spike_1e-4", 1e-4), ("spike_1e-5", 1e-5),
             ("spike_1e-6", 1e-6)]
PLANTED_SPACES = {
    "all_but_quality": [a for a in AXES if a != "mean_error_prob"],
    "coverage_only": COVERAGE_AXES,
    "composition_only": COMPOSITION_AXES,
}


def knn_purity(Xq: np.ndarray, Xbg: np.ndarray, k: int, rng: np.random.Generator) -> dict:
    """Do reads from ONE source find each other in the joint space?

    This is the operational form of "a cluster that behaves like a genome": pool a set of
    reads that genuinely share an origin with a background of residue reads, and ask what
    fraction of each query read's k nearest neighbours share its origin. The chance rate is
    the query's share of the pool, so purity at the chance rate means a clustering of any
    kind could not recover that genome, however many axes it is given. Reported for the
    planted arms, where the answer is known, and for same-sized random residue subsets,
    which is the no-structure control measured rather than assumed.
    """
    X = np.vstack([Xq, Xbg])
    nq = len(Xq)
    is_q = np.r_[np.ones(nq, bool), np.zeros(len(Xbg), bool)]
    # Ties must be broken at random. These axes are degenerate for a pair with no witnesses
    # anywhere (coverage zero, spread zero, every k-mer unwitnessed), so thousands of reads
    # share a coordinate exactly; sklearn returns neighbours in index order, which would
    # hand every tie to whichever population was stacked first and manufacture the result.
    perm = rng.permutation(len(X))
    Xp, qp = X[perm], is_q[perm]
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(Xp))).fit(Xp)
    _d, idx = nn.kneighbors(Xp[qp])
    same = [qp[row[1:]].mean() for row in idx]
    return {"n_query": int(nq), "k": k, "purity": round(float(np.mean(same)), 4),
            "chance": round(nq / len(X), 4),
            "enrichment": round(float(np.mean(same)) / (nq / len(X)), 3)}


def separability(Xr: np.ndarray, Xs: np.ndarray, seed: int) -> dict:
    """Joint AUC against the best single axis, on identical folds.

    Also reported with the axis groups split, because the comparison is confounded by
    construction and the split is what makes the confound visible: the cascade subtracted
    the sample's OWN assembled contigs, and assembly needs depth, so residue reads are
    low-coverage partly by definition. The composition-only arm is free of that, and the
    coverage-only arm is what the confound can reach at most.
    """
    X = np.vstack([Xr, Xs])
    y = np.r_[np.ones(len(Xr)), np.zeros(len(Xs))]
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    joint, single = [], {a: [] for a in AXES}
    for tr, te in skf.split(X, y):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=2000, C=1.0))
        m.fit(X[tr], y[tr])
        joint.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
        for i, a in enumerate(AXES):
            s = roc_auc_score(y[te], X[te][:, i])
            single[a].append(max(s, 1 - s))
    groups = {}
    for gname, cols in (("coverage_only", COVERAGE_AXES),
                        ("composition_only", COMPOSITION_AXES + ["mean_error_prob"])):
        sel = [AXES.index(a) for a in cols]
        gs = []
        for tr, te in skf.split(X, y):
            m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
            m.fit(X[tr][:, sel], y[tr])
            gs.append(roc_auc_score(y[te], m.predict_proba(X[te][:, sel])[:, 1]))
        groups[gname] = round(float(np.mean(gs)), 4)
    best = max(single, key=lambda a: float(np.mean(single[a])))
    return {
        "group_auc": groups,"joint_auc": float(np.mean(joint)), "joint_auc_sd": float(np.std(joint)),
            "best_single_axis": best,
            "best_single_auc": float(np.mean(single[best])),
            "single_auc": {a: round(float(np.mean(v)), 4) for a, v in single.items()},
            "n_residue": len(Xr), "n_subtracted": len(Xs)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1, help="analysis RNG; blocks are 1..5")
    ap.add_argument("--work", type=Path, default=ROOT / "data/residue_structure")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_structure")
    ap.add_argument("--k-max", type=int, default=20)
    ap.add_argument("--knn", type=int, default=10)
    ap.add_argument("--planted-only", action="store_true",
                    help="recompute the planted and null readouts, reusing the committed "
                         "separability and gap results for the five seeds")
    ap.add_argument("--subtracted-sample", type=int, default=20000)
    args = ap.parse_args()
    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows, dropped = load(args.work)
    by_block: dict[int, list[dict]] = {}
    for r in rows:
        by_block.setdefault(r["block"], []).append(r)

    resid_names = {s: set(open(ROOT / f"data/residue_structure/residue_names_s{s}.txt")
                          .read().split()) for s in range(1, 6)}

    result = {"axes": AXES, "pairs_dropped_no_kmers": dropped,
              "min_kmers": MIN_KMERS, "seeds": {}, "planted": {}, "null": {}}
    prev_path = args.out / f"structure_s{args.seed}.json"
    if args.planted_only and prev_path.exists():
        result["seeds"] = json.loads(prev_path.read_text())["seeds"]

    # ---- the two reference clouds, pooled (they are reference, not per-seed claims)
    nulls = matrix(by_block.get(9, []))
    spike = by_block.get(0, [])

    for seed in ([] if args.planted_only and result["seeds"] else range(1, 6)):
        blk = by_block.get(seed, [])
        res = [r for r in blk if r["name"] in resid_names[seed]]
        sub = [r for r in blk if r["name"] not in resid_names[seed]]
        sub = sub[:args.subtracted_sample]
        Xr, Xs = matrix(res), matrix(sub)
        sep = separability(Xr, Xs, seed)
        # z-score the residue cloud on itself; clustering is about its own shape
        mu, sd = Xr.mean(0), Xr.std(0) + 1e-9
        gap = gap_statistic((Xr - mu) / sd, list(range(1, args.k_max + 1)), rng)
        best_k = max(
            [g for g in gap],
            key=lambda g: g["gap"] - (0 if g["k"] == gap[-1]["k"] else 0),
        )["k"]
        # Tibshirani's own rule: smallest k whose gap >= gap(k+1) - s(k+1)
        rule_k = gap[-1]["k"]
        for i in range(len(gap) - 1):
            if gap[i]["gap"] >= gap[i + 1]["gap"] - gap[i + 1]["s_k"]:
                rule_k = gap[i]["k"]
                break
        result["seeds"][seed] = {
            "separability": sep,
            "gap": gap,
            "gap_max_k": best_k,
            "gap_rule_k": rule_k,
            "axis_means": {
                "residue": {a: round(float(Xr[:, i].mean()), 4) for i, a in enumerate(AXES)},
                "subtracted": {a: round(float(Xs[:, i].mean()), 4) for i, a in enumerate(AXES)},
            },
        }

    # ---- the planted control, five independent simulations per abundance
    #
    # Each replicate is its own substituted copy of the genome and its own reads, so the
    # spread below is simulation spread and not a single lucky draw - the repo's five-seed
    # rule applied to the one arm of this experiment that has a known answer. The
    # neighbourhood test is run against the residue cloud of block 1, which is the cloud a
    # real planted genome would have to be found inside.
    res1 = [r for r in by_block[1] if r["name"] in resid_names[1]]
    Xr1 = matrix(res1)
    mu, sd = Xr1.mean(0), Xr1.std(0) + 1e-9
    bg = (Xr1 - mu) / sd
    by_tag: dict[str, list[dict]] = {}
    for r in spike:
        by_tag.setdefault(r["name"].split("|")[0], []).append(r)
    for rate_tag, _rate in ARMS_TAGS:
        reps = sorted(t for t in by_tag if t.startswith(rate_tag + "_r"))
        if not reps:
            continue
        per_rep, per_ctrl, means = [], [], []
        for t in reps:
            arm = by_tag[t]
            Xa = (matrix(arm) - mu) / sd
            ctrl_idx = rng.choice(len(bg), size=min(len(arm), len(bg)), replace=False)
            sub_bg = np.delete(bg, ctrl_idx, axis=0)
            per_rep.append({sp: knn_purity(Xa[:, [AXES.index(a) for a in cols]],
                                           bg[:, [AXES.index(a) for a in cols]],
                                           args.knn, rng)
                            for sp, cols in PLANTED_SPACES.items()})
            per_ctrl.append({sp: knn_purity(bg[ctrl_idx][:, [AXES.index(a) for a in cols]],
                                            sub_bg[:, [AXES.index(a) for a in cols]],
                                            args.knn, rng)
                             for sp, cols in PLANTED_SPACES.items()})
            means.append(matrix(arm).mean(0))
        def agg(rows, key):
            v = [r[key]["purity"] for r in rows]
            return {"mean": round(float(np.mean(v)), 4), "min": round(min(v), 4),
                    "max": round(max(v), 4),
                    "chance": round(float(np.mean([r[key]["chance"] for r in rows])), 5)}
        result["planted"][rate_tag] = {
            "replicates": len(reps),
            "n_pairs_each": len(by_tag[reps[0]]),
            "knn_purity": {sp: agg(per_rep, sp) for sp in PLANTED_SPACES},
            "knn_purity_residue_control": {sp: agg(per_ctrl, sp) for sp in PLANTED_SPACES},
            "axis_means": {a: round(float(np.mean([m[i] for m in means])), 4)
                           for i, a in enumerate(AXES)},
        }
    if len(nulls):
        # The shuffled null is a reference cloud for composition and a DEGENERATE one for
        # coverage: any shuffle destroys every 31-mer, so a shuffled read has no witnesses
        # by construction and collapses to the corner of the space. Reported, because a
        # calibration that cannot discriminate is still the field's honest calibration.
        result["null"] = {
            "n_pairs": int(len(nulls)),
            "axis_means": {a: round(float(nulls[:, i].mean()), 4)
                           for i, a in enumerate(AXES)},
        }

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                         text=True).stdout.strip()
    (args.out / f"structure_s{args.seed}.json").write_text(json.dumps(result, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/residue_structure_analyse.py",
        "git_sha": sha,
        "parameters": vars(args) | {"work": str(args.work), "out": str(args.out)},
        "inputs": {
            "coverage": str(args.work / "coverage.tsv"),
            "composition": str(args.work / "composition.tsv"),
            "residue_membership": "results/subtract_recut/residue_rows_s{1..5}.jsonl "
                                  "(claude/project-thread-eswrkj)",
        },
        "wall_seconds": round(time.perf_counter() - t0, 1),
    }, indent=2))
    print(json.dumps({s: {"joint": v["separability"]["joint_auc"],
                          "best_single": v["separability"]["best_single_auc"],
                          "axis": v["separability"]["best_single_axis"],
                          "gap_rule_k": v["gap_rule_k"]}
                      for s, v in result["seeds"].items()}, indent=2))


if __name__ == "__main__":
    main()
