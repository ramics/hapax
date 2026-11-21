#!/usr/bin/env python3
"""Tables for results/weak_hits/ from its committed cell files.

Generated, never hand-written: a number typed into a table is a number without
provenance.

Usage: python3 experiments/report_weak_hits.py [--dir results/weak_hits]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
LABEL = {("background", "dNA"): "real library reads (known material)",
         ("shuffle", "dNA"): "**composition-matched shuffle (null)**"}


def order(key):
    pile, rung = key[0], key[1]
    rank = {"spike": 0, "background": 1, "shuffle": 2}.get(pile, 9)
    return (rank, int(rung[1:]) if rung[1:].isdigit() else -1)


def label(key):
    return LABEL.get(key, f"holdout, {key[1][1:]}% divergence")


def agg(rows, key, cond, field):
    v = [r[field] for r in rows
         if (r["pile"], r["rung"]) == key and r["condition"] == cond and field in r]
    return (float(np.mean(v)), float(np.min(v)), float(np.max(v))) if v else None


def fmt(t, pct=True, dp=1):
    if t is None:
        return "n/a"
    m, lo, hi = t
    if pct:
        return f"{100 * m:.{dp}f} ({100 * lo:.{dp}f}–{100 * hi:.{dp}f})"
    return f"{m:.{dp}f} ({lo:.{dp}f}–{hi:.{dp}f})"


def sci(t):
    if t is None:
        return "n/a"
    return f"{t[0]:.2g} ({t[1]:.2g}–{t[2]:.2g})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/weak_hits")
    args = ap.parse_args()
    d = REPO / args.dir
    bands = json.load(open(d / "bands.json"))
    conv = json.load(open(d / "conversion.json"))
    man = json.load(open(d / "manifest.json"))
    sig = man["params"]["significance_cutoff_evalue"]
    perm = man["params"]["reporting_cutoff_evalue"]
    keys = sorted({(r["pile"], r["rung"]) for r in bands}, key=order)

    out = []
    w = out.append
    w("# Is there headroom in weak protein hits?")
    w("")
    w(f"Generated from `bands.json` and `conversion.json` by "
      f"`experiments/report_weak_hits.py`. {len(man['seeds'])} seeds, "
      f"{man['params']['query_count']:,} queries, DIAMOND "
      f"`{man['params']['diamond_mode']}` reported down to E <= {perm:g} so that a "
      f"weak-hit population exists at all, with E <= {sig:g} as the call line. Mean over "
      f"seeds with the min–max range.")
    w("")

    w("## 1. The weak band reads the same on real sequence as on shuffles")
    w("")
    w("Mate 1. 'Any hit' is at the reporting cutoff, which is why it is near 100% "
      "everywhere **including the null** — that column is the chance band, not evidence.")
    w("")
    w("| pile | any hit (E<=100) | significant (E<=1e-3) | median E of hits | median bits | median %id |")
    w("|---|---|---|---|---|---|")
    for k in keys:
        w(f"| {label(k)} | " + " | ".join([
            fmt(agg(bands, k, "m1", "any_hit")),
            fmt(agg(bands, k, "m1", "sig")),
            sci(agg(bands, k, "m1", "median_evalue_of_hits")),
            fmt(agg(bands, k, "m1", "median_bitscore_of_hits"), pct=False),
            fmt(agg(bands, k, "m1", "median_pident_of_hits"), pct=False),
        ]) + " |")
    w("")
    w("The weak hits themselves, with the null beside them:")
    w("")
    w("| pile | weak hits (1e-3 < E <= 100) | median bits of those | median %id of those |")
    w("|---|---|---|---|")
    for k in keys:
        w(f"| {label(k)} | " + " | ".join([
            fmt(agg(bands, k, "m1", "weak_frac")),
            fmt(agg(bands, k, "m1", "median_bitscore_of_weak"), pct=False),
            fmt(agg(bands, k, "m1", "median_pident_of_weak"), pct=False),
        ]) + " |")
    w("")

    w("## 2. Merging converts weak hits, and converts none of the noise")
    w("")
    w("Among pairs whose mate 1 produced only a weak hit, the fraction whose merged read "
      "produces a significant one. The control column is pairs already significant on "
      "mate 1, which merging should keep. The last column is the length-arithmetic null: "
      "a merged read's bitscore divided by its mate's scaled by the length ratio, where "
      "1.0 means a longer query simply accumulated more of the same score.")
    w("")
    w("| pile | n weak on mate 1 | converted | kept (already significant) | observed/length-predicted bits |")
    w("|---|---|---|---|---|")
    ckeys = sorted({(r["pile"], r["rung"]) for r in conv}, key=order)
    for k in ckeys:
        rs = [r for r in conv if (r["pile"], r["rung"]) == k]

        def mm(f, pct=True, dp=1):
            v = [r[f] for r in rs if f in r]
            if not v:
                return "n/a"
            a = np.array(v, dtype=float)
            if pct:
                return f"{100 * a.mean():.{dp}f} ({100 * a.min():.{dp}f}–{100 * a.max():.{dp}f})"
            return f"{a.mean():.3f} ({a.min():.3f}–{a.max():.3f})"
        nweak = int(np.mean([r["n_weak_m1"] for r in rs]))
        w(f"| {label(k)} | {nweak:,} | {mm('converted')} | {mm('kept')} | "
          f"{mm('bitscore_observed_over_length_prediction', pct=False)} |")
    w("")

    w("## 3. The merged gain at a fixed E-value, not a bitscore floor")
    w("")
    w("`../profile_reach/` §2 measured the merged advantage at a bitscore floor, where a "
      "longer query has an arithmetic edge. An E-value penalises query length, so if the "
      "gain were arithmetic it would shrink here. Same overlapping pairs, merged against "
      "mate 1 of the same molecule.")
    w("")
    w("| pile | E<=1e-3 merged / mate 1 | E<=1e-5 merged / mate 1 | E<=1e-10 merged / mate 1 |")
    w("|---|---|---|---|")
    for k in ckeys:
        rs = [r for r in conv if (r["pile"], r["rung"]) == k]
        cells = []
        for tag in ("1e-3", "1e-5", "1e-10"):
            def one(cond):
                v = [r[f"sig_E{tag}_{cond}"] for r in rs if f"sig_E{tag}_{cond}" in r]
                if not v:
                    return "n/a"
                a = np.array(v, dtype=float)
                return f"{100 * a.mean():.1f} ({100 * a.min():.1f}–{100 * a.max():.1f})"
            cells.append(f"{one('merged')} / {one('m1')}")
        w(f"| {label(k)} | " + " | ".join(cells) + " |")
    w("")
    w("The advantage survives at every threshold and **grows** at the strictest, which is "
      "the opposite of what a length artefact does.")
    w("")

    w("## Provenance")
    w("")
    w(f"- generator `{man['script']}`, content hash `{man['script_sha256'][:16]}`")
    w(f"- this report `experiments/report_weak_hits.py`, content hash "
      f"`{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}`")
    w(f"- git SHA `{man['git_sha']}`")
    w(f"- {man['params']['note']}")
    w(f"- DIAMOND time {man['timings_s']['diamond']:.0f}s "
      f"(0 means the search was reused from an earlier identical run)")
    w("")

    (d / "TABLES.md").write_text("\n".join(out) + "\n")
    print(f"wrote {d / 'TABLES.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
