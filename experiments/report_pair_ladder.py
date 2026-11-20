"""Summarise results/pair_ladder/rows.jsonl into the tables the write-up quotes.

Why the readouts are these. The question is whether a tool returns each lineage as its
own contigs, so the headline is *separated breadth*: the fraction of a lineage's genome
covered by contigs whose every called window is that lineage. The weaker lineage of the
two is what matters (min over a and b), since returning one lineage and burying the other
is the failure mode. Mosaic mass fraction is reported beside it from the same runs, as
CLAUDE.md requires for any contiguity number. Spread is given two ways: SD over all
pair x seed runs, and the range of per-pair means inside a rung, which is the number the
three-pairs-per-rung design exists to expose.

Usage: python experiments/report_pair_ladder.py results/pair_ladder
"""
from __future__ import annotations

import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def mean_sd(v):
    return (st.mean(v), st.stdev(v) if len(v) > 1 else 0.0) if v else (float("nan"), 0.0)


def main() -> None:
    d = Path(sys.argv[1] if len(sys.argv) > 1 else "results/pair_ladder")
    rows = [json.loads(x) for x in (d / "rows.jsonl").open()]
    mix = [r for r in rows if r["mode"] == "mix"]
    ctl = [r for r in rows if r["mode"] != "mix"]
    tools = sorted({r["tool"] for r in mix})
    rungs = sorted({(r["divergence_rung"] if "divergence_rung" in r else r["rung"]) for r in mix},
                   key=lambda t: min(r["divergence"] for r in mix if r["rung"] == t))
    covs = sorted({r["cov_per_lineage"] for r in mix})
    out = []
    out.append(f"{len(mix)} mixture runs, {len(ctl)} single-lineage control runs, "
               f"{sum(1 for r in rows if r['error'])} runs with a tool error.\n")

    def sep_min(r):
        return min(r["sep_breadth_a"], r["sep_breadth_b"])

    def mos_frac(r):
        return r["mosaic_mass"] / r["mass"] if r["mass"] else 0.0

    for metric, fn, label in [("sep", sep_min, "separated breadth of the weaker lineage"),
                              ("mos", mos_frac, "mosaic fraction of assembled mass"),
                              ("any", lambda r: min(r["breadth_a"], r["breadth_b"]),
                               "breadth of the weaker lineage by any contig")]:
        for cov in covs:
            out.append(f"\n### {label}, {cov:g}x per lineage\n")
            out.append("mean ± SD over 3 pairs x 5 seeds; [range of the 3 per-pair means]\n")
            out.append("| rung | divergence | " + " | ".join(tools) + " |")
            out.append("|---|---|" + "---|" * len(tools))
            for rung in rungs:
                divs = sorted({r["divergence"] for r in mix if r["rung"] == rung})
                cells = []
                for t in tools:
                    sel = [r for r in mix if r["rung"] == rung and r["tool"] == t and r["cov_per_lineage"] == cov]
                    v = [fn(r) for r in sel]
                    m, s = mean_sd(v)
                    per = defaultdict(list)
                    for r in sel:
                        per[(r["a"], r["b"])].append(fn(r))
                    pm = [st.mean(x) for x in per.values()]
                    cells.append(f"{m:.2f} ± {s:.2f} [{min(pm):.2f}-{max(pm):.2f}]" if pm else "-")
                out.append(f"| {rung} | {divs[0]:.1%}-{divs[-1]:.1%} | " + " | ".join(cells) + " |")

    out.append("\n### Single-lineage controls (labeller null on real assembler output)\n")
    out.append("| tool | runs | runs with any mosaic contig | mean breadth of the present lineage |")
    out.append("|---|---|---|---|")
    for t in sorted({r["tool"] for r in ctl}):
        sel = [r for r in ctl if r["tool"] == t]
        own = [r["breadth_a"] if r["mode"] == "a_alone" else r["breadth_b"] for r in sel]
        out.append(f"| {t} | {len(sel)} | {sum(r['mosaic_contigs'] > 0 for r in sel)} | {st.mean(own):.3f} |")

    out.append("\n### Output size (redundancy)\n")
    out.append("| tool | median contigs per mixture run | median assembled mass / (len a + len b) |")
    out.append("|---|---|---|")
    for t in tools:
        sel = [r for r in mix if r["tool"] == t]
        out.append(f"| {t} | {st.median(r['contigs'] for r in sel):g} | "
                   f"{st.median(r['mass'] / 23800 for r in sel):.2f} |")
    out += switch_tables(mix, tools, rungs, covs, d)
    text = "\n".join(out) + "\n"
    (d / "summary.md").write_text(text)
    print(text)


def switch_tables(mix, tools, rungs, covs, d: Path) -> list[str]:
    """Lineage switches per 10 kb of labelled contig, from the kept per-contig labels.

    The strict 'pure contig' readout scores a full-length contig with one swapped 300 nt
    window the same as a 50:50 fusion. Switch density separates the two: a near-perfect
    separation with rare local swaps sits near 1 per 10 kb, a collapse much higher.
    """
    import gzip
    keep = Path("/mnt/project-files/strain_ladder/contigs")
    out = []
    for cov in covs:
        out.append(f"\n### lineage switches per 10 kb of labelled contig, {cov:g}x per lineage\n")
        out.append("| rung | " + " | ".join(tools) + " |")
        out.append("|---|" + "---|" * len(tools))
        for rung in rungs:
            cells = []
            for t in tools:
                v = []
                for r in mix:
                    if r["rung"] != rung or r["tool"] != t or r["cov_per_lineage"] != cov:
                        continue
                    tag = f"{r['rung']}_{r['a']}_{r['b']}_mix_c{cov}_s{r['seed']}"
                    f = keep / tag / f"{t}.labels.json.gz"
                    if not f.exists():
                        continue
                    per = json.load(gzip.open(f, "rt"))
                    nt = sum(len(c["calls"]) for c in per) * 300
                    if nt:
                        v.append(sum(c["switches"] for c in per) / nt * 1e4)
                m, sd = mean_sd(v)
                cells.append(f"{m:.1f} ± {sd:.1f}" if v else "-")
            out.append(f"| {rung} | " + " | ".join(cells) + " |")
    return out


if __name__ == "__main__":
    main()
