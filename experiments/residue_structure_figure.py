"""The joint space as SPEC.md §8.3 describes it, from committed result files only.

Four panels, each one of the readouts rather than decoration:

  A  the joint space itself, the two axes that carry the separation, with the residue, the
     subtracted cloud and the shuffled null drawn underneath as the two reference clouds,
     and the planted arms over the top at four abundances;
  B  the gap statistic against k, with the covariance-matched reference, which is where
     "structured or diffuse" is actually decided;
  C  separability: the joint space against the best single axis, and against each axis
     group, so the confound in the comparison is visible in the figure and not only in the
     text;
  D  nearest-neighbour purity for the planted arms against chance, the test of whether a
     clustering could recover a genome's reads at all.

Every number drawn here is read from results/residue_structure/*.json.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
C = {"residue": "#2a4d8f", "subtracted": "#9aa4b2", "null": "#c86a4a"}
ARM_C = ["#1b7f5c", "#2e9e74", "#58b894", "#9ed6bf"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path, default=ROOT / "results/residue_structure")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    out = args.out or args.results / "joint_space.png"

    st = json.loads((args.results / "structure_s1.json").read_text())
    pts = json.loads((args.results / "joint_points_s1.json").read_text())
    axes = pts["axes"]
    ix, iy = axes.index("log_cov"), axes.index("zero_frac")

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    fig.suptitle("The residue in the joint space of its fields: Berlin wastewater "
                 "SRR23973193, block 1 of 5", fontsize=13)

    # ---- A: the joint space
    a = ax[0][0]
    for pop in ("subtracted", "null", "residue"):
        P = np.array(pts["populations"][pop]["points"])
        a.scatter(P[:, ix], P[:, iy], s=3, alpha=0.25, c=C[pop], linewidths=0,
                  label=f"{pop} (n={pts['populations'][pop]['n_total']:,})")
    for i, rate in enumerate(["1e-3", "1e-4", "1e-5", "1e-6"]):
        # Replicate 1 of each abundance; the other four are in the committed points file
        # and plotting all five would only overdraw the same cloud.
        p = pts["populations"].get(f"spike_{rate}_r1")
        if not p or not p["points"]:
            continue
        P = np.array(p["points"])
        a.scatter(P[:, ix], P[:, iy], s=14, c=ARM_C[i], edgecolors="black",
                  linewidths=0.3, label=f"planted {rate} (n={p['n_total']})")
    a.set_xlabel("log(1 + median 31-mer witnesses elsewhere in the library)")
    a.set_ylabel("fraction of 31-mers with no witness")
    a.set_title("A. The joint space, with both reference clouds")
    a.legend(fontsize=7, markerscale=2, loc="upper right")

    # ---- B: the gap statistic
    b = ax[0][1]
    for seed, v in st["seeds"].items():
        ks = [g["k"] for g in v["gap"]]
        b.plot(ks, [g["gap"] for g in v["gap"]], marker="o", ms=3, lw=1,
               label=f"seed {seed}")
        b.fill_between(ks, [g["gap"] - g["s_k"] for g in v["gap"]],
                       [g["gap"] + g["s_k"] for g in v["gap"]], alpha=0.1)
    b.axhline(0, color="black", lw=0.8)
    b.set_xlabel("k (k-means clusters)")
    b.set_ylabel("gap: log W of covariance-matched reference − log W of data")
    b.set_title("B. Clustering beyond the covariance\n"
                f"Tibshirani rule picks k = "
                f"{', '.join(str(v['gap_rule_k']) for v in st['seeds'].values())}")
    b.legend(fontsize=7)

    # ---- C: separability
    c = ax[1][0]
    seeds = sorted(st["seeds"])
    bars = {
        "joint (all axes)": [st["seeds"][s]["separability"]["joint_auc"] for s in seeds],
        "best single axis": [st["seeds"][s]["separability"]["best_single_auc"]
                             for s in seeds],
        "coverage axes only": [st["seeds"][s]["separability"]["group_auc"]["coverage_only"]
                               for s in seeds],
        "composition axes only":
            [st["seeds"][s]["separability"]["group_auc"]["composition_only"]
             for s in seeds],
    }
    xs = np.arange(len(bars))
    c.bar(xs, [np.mean(v) for v in bars.values()],
          yerr=[np.std(v) for v in bars.values()], capsize=4,
          color=["#2a4d8f", "#6f87b5", "#9aa4b2", "#c86a4a"])
    c.axhline(0.5, color="black", lw=0.8, ls=":")
    c.set_xticks(xs)
    c.set_xticklabels(bars.keys(), fontsize=8, rotation=12)
    c.set_ylim(0.45, 1.0)
    c.set_ylabel("cross-validated AUC, residue against subtracted")
    best = st["seeds"][seeds[0]]["separability"]["best_single_axis"]
    c.set_title(f"C. Separability, 5 seeds (best single axis: {best})")

    # ---- D: nearest-neighbour purity
    d = ax[1][1]
    tags = [t for t in ("spike_1e-3", "spike_1e-4", "spike_1e-5", "spike_1e-6")
            if t in st["planted"]]
    w = 0.38
    xs = np.arange(len(tags))
    sp = "all_but_quality"
    pur = [st["planted"][t]["knn_purity"][sp]["mean"] for t in tags]
    pe = [[st["planted"][t]["knn_purity"][sp]["mean"]
           - st["planted"][t]["knn_purity"][sp]["min"] for t in tags],
          [st["planted"][t]["knn_purity"][sp]["max"]
           - st["planted"][t]["knn_purity"][sp]["mean"] for t in tags]]
    ctr = [st["planted"][t]["knn_purity_residue_control"][sp]["mean"] for t in tags]
    ce = [[st["planted"][t]["knn_purity_residue_control"][sp]["mean"]
           - st["planted"][t]["knn_purity_residue_control"][sp]["min"] for t in tags],
          [st["planted"][t]["knn_purity_residue_control"][sp]["max"]
           - st["planted"][t]["knn_purity_residue_control"][sp]["mean"] for t in tags]]
    d.bar(xs - w / 2, pur, w, yerr=pe, capsize=3, color="#1b7f5c",
          label="planted genome, 5 simulations")
    d.bar(xs + w / 2, ctr, w, yerr=ce, capsize=3, color="#9aa4b2",
          label="same-sized random residue subset")
    for i, t in enumerate(tags):
        d.plot([xs[i] - w, xs[i] + w],
               [st["planted"][t]["knn_purity"][sp]["chance"]] * 2,
               color="black", lw=1, ls=":")
    d.set_xticks(xs)
    d.set_xticklabels([t.split("_")[1] + f"\nn={st['planted'][t]['n_pairs_each']}"
                       for t in tags], fontsize=8)
    d.set_ylabel(f"fraction of the 10 "
                 "nearest neighbours from the same source")
    d.set_title("D. Could a neighbourhood recover one genome's reads?\n"
                "dotted line: chance, the query's share of the pool")
    d.legend(fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=170)
    print("wrote", out)


if __name__ == "__main__":
    main()
