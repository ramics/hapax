#!/usr/bin/env python3
"""Tables for results/profile_reach/ from its committed cell files.

Figures and tables are generated from committed result files, never hand-written: a
number typed into a table is a number without provenance, and this repo has had to
withdraw claims whose prose drifted from its outputs (CLAUDE.md).

Usage: python3 experiments/report_profile_reach.py [--dir results/profile_reach]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent


def agg(rows, pile, cond, field, rung=None):
    """Mean over seeds with the min-max range, for one cell of the design."""
    vals = [r[field] for r in rows
            if r["pile"] == pile and r["condition"] == cond
            and (rung is None or r["rung"] == rung) and field in r]
    if not vals:
        return None
    return float(np.mean(vals)), float(np.min(vals)), float(np.max(vals))


def fmt(t, pct=True):
    if t is None:
        return "n/a"
    m, lo, hi = t
    if pct:
        return f"{100 * m:.1f} ({100 * lo:.1f}–{100 * hi:.1f})"
    return f"{m:.2f} ({lo:.2f}–{hi:.2f})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/profile_reach")
    args = ap.parse_args()
    d = REPO / args.dir
    rows = json.load(open(d / "cells.json"))
    paired = json.load(open(d / "paired_merged.json"))
    man = json.load(open(d / "manifest.json"))
    rungs = [f"d{r}" for r in man["params"]["rungs_pct"]]
    seeds = man["seeds"]

    out = []
    w = out.append
    w("# Profile reach, and the pair as one molecule")
    w("")
    w(f"Generated from `cells.json` and `paired_merged.json` by "
      f"`experiments/report_profile_reach.py`. {len(seeds)} seeds "
      f"({man['params']['n_sim_pairs_per_seed_per_rung']} simulated pairs per seed per "
      f"rung, {man['params']['n_background_pairs_per_seed']} real pairs per seed). "
      f"Every cell is mean over seeds with the min–max range.")
    w("")
    w("Operating points, stated because thresholds dominate everything else here: "
      f"DIAMOND `{man['params']['diamond_mode']}` against "
      f"{' + '.join(man['params']['diamond_db'])}, reported at a bitscore floor; Pfam-A "
      "searched with no score cutoff and reported at per-domain E-value floors, with the "
      "families' own curated gathering thresholds (`GA`) as one point on that curve.")
    w("")

    # ---- table 1: reach of each instrument along the ladder, mate 1
    w("## 1. Reach along the divergence ladder (mate 1, 107 nt)")
    w("")
    w("| pile | DIAMOND bs>=50 | DIAMOND bs>=40 | Pfam E<=1e-2 | Pfam E<=1e-5 | Pfam at GA |")
    w("|---|---|---|---|---|---|")
    for rung in rungs:
        label = f"holdout, {rung[1:]}% divergence"
        w(f"| {label} | " + " | ".join([
            fmt(agg(rows, "spike", "m1", "diamond_bs50", rung)),
            fmt(agg(rows, "spike", "m1", "diamond_bs40", rung)),
            fmt(agg(rows, "spike", "m1", "pfam_e0.01", rung)),
            fmt(agg(rows, "spike", "m1", "pfam_e1e-05", rung)),
            fmt(agg(rows, "spike", "m1", "pfam_ga", rung)),
        ]) + " |")
    for pile, label in (("background", "real library reads (known material)"),
                        ("shuffle", "composition-matched shuffle (null)")):
        w(f"| {label} | " + " | ".join([
            fmt(agg(rows, pile, "m1", "diamond_bs50")),
            fmt(agg(rows, pile, "m1", "diamond_bs40")),
            fmt(agg(rows, pile, "m1", "pfam_e0.01")),
            fmt(agg(rows, pile, "m1", "pfam_e1e-05")),
            fmt(agg(rows, pile, "m1", "pfam_ga")),
        ]) + " |")
    w("")

    # ---- table 2: the gap, which is the actual deliverable
    w("## 2. The gap: no translated hit, yes profile hit")
    w("")
    w("The reportable quantity, because a profile hit rate on a coding background is a "
      "familiarity signal. Read the background and shuffle rows FIRST: if the gap reads "
      "the same on known material, it describes the database rather than the query.")
    w("")
    w("| pile | gap (bs<50, E<=1e-2) | gap (bs<50, E<=1e-5) | gap (bs<40, E<=1e-2) | gap (bs<50, GA) |")
    w("|---|---|---|---|---|")
    for rung in rungs:
        w(f"| holdout, {rung[1:]}% divergence | " + " | ".join([
            fmt(agg(rows, "spike", "m1", "gap_bs50_e0.01", rung)),
            fmt(agg(rows, "spike", "m1", "gap_bs50_e1e-05", rung)),
            fmt(agg(rows, "spike", "m1", "gap_bs40_e0.01", rung)),
            fmt(agg(rows, "spike", "m1", "gap_bs50_ga", rung)),
        ]) + " |")
    for pile, label in (("background", "real library reads (known material)"),
                        ("shuffle", "composition-matched shuffle (null)")):
        w(f"| {label} | " + " | ".join([
            fmt(agg(rows, pile, "m1", "gap_bs50_e0.01")),
            fmt(agg(rows, pile, "m1", "gap_bs50_e1e-05")),
            fmt(agg(rows, pile, "m1", "gap_bs40_e0.01")),
            fmt(agg(rows, pile, "m1", "gap_bs50_ga")),
        ]) + " |")
    w("")

    # ---- table 3: per mate, never pooled
    w("## 3. Per mate, never pooled")
    w("")
    w("Mate 2 is 93 nt against mate 1's 107 after the non-template prefix trim, so a "
      "pooled figure averages two instruments.")
    w("")
    w("| pile | mate | mean query nt | peptide segments | DIAMOND bs>=50 | Pfam E<=1e-2 |")
    w("|---|---|---|---|---|---|")
    for rung in rungs:
        for cond in ("m1", "m2"):
            w(f"| holdout {rung[1:]}% | {cond} | " + " | ".join([
                fmt(agg(rows, "spike", cond, "mean_query_len", rung), pct=False),
                fmt(agg(rows, "spike", cond, "mean_peptide_segments", rung), pct=False),
                fmt(agg(rows, "spike", cond, "diamond_bs50", rung)),
                fmt(agg(rows, "spike", cond, "pfam_e0.01", rung)),
            ]) + " |")
    for cond in ("m1", "m2"):
        w(f"| real library | {cond} | " + " | ".join([
            fmt(agg(rows, "background", cond, "mean_query_len"), pct=False),
            fmt(agg(rows, "background", cond, "mean_peptide_segments"), pct=False),
            fmt(agg(rows, "background", cond, "diamond_bs50")),
            fmt(agg(rows, "background", cond, "pfam_e0.01")),
        ]) + " |")
    w("")

    # ---- table 4: the paired merged comparison
    w("## 4. The pair as one molecule: merged against mate 1 of the SAME pair")
    w("")
    w("Restricted to pairs whose mates overlap, so the comparison is within a molecule "
      "rather than between insert lengths. The merged read carries no error-correction "
      "credit for its overlap (stated bias 2), so these are lower bounds.")
    w("")
    w("| pile | mean merged nt | mean mate-1 nt | DIAMOND bs>=50 merged | same, mate 1 | "
      "Pfam E<=1e-2 merged | same, mate 1 |")
    w("|---|---|---|---|---|---|---|")
    by = defaultdict(list)
    for r in paired:
        by[(r["pile"], r["rung"])].append(r)
    # ordered by pile then by rung as a NUMBER: sorting the keys as strings puts the 5%
    # rung after the 30% one, and a ladder out of order invites being read as a trend
    order = {"spike": 0, "background": 1, "shuffle": 2}
    labels = {"spike": None, "background": "real library reads (known material)",
              "shuffle": "composition-matched shuffle (null)"}

    def rung_num(r):
        return int(r[1:]) if r[1:].isdigit() else -1

    for (pile, rung), rs in sorted(by.items(),
                                   key=lambda kv: (order.get(kv[0][0], 9),
                                                   rung_num(kv[0][1]))):
        def mm(f):
            v = [r[f] for r in rs]
            return float(np.mean(v)), float(np.min(v)), float(np.max(v))
        label = labels[pile] or f"holdout, {rung[1:]}% divergence"
        w(f"| {label} | " + " | ".join([
            fmt(mm("mean_merged_len"), pct=False), fmt(mm("mean_m1_len"), pct=False),
            fmt(mm("diamond_bs50_merged")), fmt(mm("diamond_bs50_m1")),
            fmt(mm("pfam_e0.01_merged")), fmt(mm("pfam_e0.01_m1")),
        ]) + " |")
    w("")
    w("## Provenance")
    w("")
    w(f"- generator `{man['script']}`, content hash `{man['script_sha256'][:16]}`")
    w(f"- this report `experiments/report_profile_reach.py`, content hash "
      f"`{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}`")
    w(f"- git SHA `{man['git_sha']}`")
    w(f"- holdout: {man['inputs']['holdout']}, {man['inputs']['holdout_len_nt']} nt")
    w(f"- library: {man['inputs']['library']}")
    w(f"- timings: DIAMOND {man['timings_s']['diamond']:.0f}s, "
      f"Pfam {man['timings_s']['pfam']:.0f}s, total {man['timings_s']['total']:.0f}s")
    w("")

    (d / "TABLES.md").write_text("\n".join(out) + "\n")
    print(f"wrote {d / 'TABLES.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
