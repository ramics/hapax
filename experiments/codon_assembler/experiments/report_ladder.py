"""Figure and paired-difference table for the diversity ladder.

The four hapax configurations run on the *same* simulated worlds, so the
codon-vs-nucleotide comparison is paired by seed. Paired differences are far
more sensitive than comparing two means with overlapping spreads, and they
are what the write-up quotes.

Usage: python experiments/report_ladder.py results/ladder
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = {"hapax-codon": "#1f5fa8", "hapax-nucleotide": "#d1495b", "hapax-codon-spaced": "#6fa8dc", "hapax-nucleotide-spaced": "#f4a3a8"}


def main(d: Path):
    rows = list(csv.DictReader(open(d / "runs.tsv"), delimiter="\t"))
    divs = sorted({float(r["diversity"]) for r in rows})
    taus = sorted({r["tau"] for r in rows if r["tau"]})
    tau = "0.9" if "0.9" in taus else taus[len(taus) // 2]
    metrics = [
        ("mean_contigs_per_found_entity", "contigs per found entity (target 1)"),
        ("chimaera_rate", "chimaera rate (contigs)"),
        ("mean_recovered", "recovered fraction (largest contig)"),
        ("hmm_fraction", "fraction of candidates sent to codon HMM"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4))
    for ax, (m, label) in zip(axes, metrics):
        for asm, col in COLORS.items():
            ys, es = [], []
            for dv in divs:
                v = [float(r[m]) for r in rows if r["assembler"] == asm and r["tau"] == tau and float(r["diversity"]) == dv]
                ys.append(st.mean(v) if v else float("nan"))
                es.append(st.stdev(v) if len(v) > 1 else 0)
            ax.errorbar(divs, ys, yerr=es, label=asm, color=col, marker="o", capsize=3, ls="--" if "spaced" in asm else "-")
        ax.set_xlabel("within-population diversity (subs/site)")
        ax.set_title(label, fontsize=10)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"hapax diversity ladder, tau={tau}, mean ± sd over seeds")
    fig.tight_layout()
    fig.savefig(d / "ladder.png", dpi=130)

    idx = {(r["seed"], r["diversity"], r["assembler"], r["tau"]): r for r in rows}
    seeds = sorted({r["seed"] for r in rows})
    lines = ["| seeds | diversity | tau | Δ contigs/entity (codon − nt), mean [min, max] | seeds codon better / worse | Δ chimeric contigs (sum) |", "|---|---|---|---|---|---|"]
    for shape in ("", "-spaced"):
        for dv in sorted({r["diversity"] for r in rows}, key=float):
            for t in taus:
                diffs, dch = [], 0
                for s in seeds:
                    c = idx.get((s, dv, "hapax-codon" + shape, t))
                    n = idx.get((s, dv, "hapax-nucleotide" + shape, t))
                    if not c or not n:
                        continue
                    diffs.append(float(c["mean_contigs_per_found_entity"]) - float(n["mean_contigs_per_found_entity"]))
                    dch += int(c["chimeric_contigs"]) - int(n["chimeric_contigs"])
                if diffs:
                    better = sum(x < 0 for x in diffs)
                    worse = sum(x > 0 for x in diffs)
                    lines.append(f"| {'spaced' if shape else 'contiguous'} | {dv} | {t} | {st.mean(diffs):+.3f} [{min(diffs):+.2f}, {max(diffs):+.2f}] | {better} / {worse} | {dch:+d} |")
    (d / "paired.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
