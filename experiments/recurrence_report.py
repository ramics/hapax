"""Turn the per-subject witness arrays into the recurrence distribution.

WHAT COMES IN. One `witness_<run>.npz` per panel library from `recurrence_pass.py`, each
holding, per query and per mate, the fraction of that mate's 31-mers present in that library
and the median number of times they occur there. Plus the panel table, which says which runs
share a water and a date.

WHAT GOES OUT, AND WHAT DELIBERATELY DOES NOT. Out: for every (pile, depth bin, mate) cell,
the distribution of how many subjects witness a read - as a histogram and as quantiles, with
the shuffled null and the identified and ribosomal piles beside it. Not out: a threshold, a
score, or a per-read verdict. The occupancy *shape* is the deliverable (SPEC.md §4.2), and
the project's standing finding is that single axes here are continuous and do not separate,
so a cut point would be inventing a decision the data has not supported yet.

THREE DENOMINATORS, because they answer different questions and differ by up to 2.8x.
`runs` counts libraries; `waters` counts date-plus-window composites, so preparation
variants of one water count once; `dates` pools an evening's two windows. A molecule present
in two preparations of the same water has recurred in the laboratory, not in the sewer, and
only the water and date denominators say so.

PRESENCE AT TWO OPERATING POINTS, since there is no single right one and the repo's rule is
that no number travels without its operating point. `any` is at least one 31-mer of the mate
present, the weakest criterion there is and the one the earlier cross-library witness figure
used. `most` is at least half the mate's 31-mers, which is close to "the molecule is there"
rather than "something that shares 31 bases is there". Both are reported for every cell; they
behave differently and that difference is information.

THE BIAS DIRECTIONS, stated rather than corrected.
  - The subject libraries differ in depth by a factor of ~3 (read counts 3-17 M), and a
    deeper library witnesses more of everything. So occupancy mixes recurrence with subject
    depth. The control for it is in the output: `subject_depth_spearman` correlates each
    subject's witnessed fraction with its read count, and the per-subject table lets the
    reader see whether a pile's pattern is anything other than that.
  - Mates are reported separately throughout: mate 2 is 93 nt against mate 1's 107 after the
    trim, so it carries 62 scorable 31-mers against 76 and is witnessed less often for a
    reason that has nothing to do with the molecule (CLAUDE.md).
  - The ribosomal pile is the familiarity reference: it should recur nearly everywhere by
    conservation, and if it does not, the instrument is wrong rather than the sewage
    interesting. Read that row first.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from subtract_known import provenance  # noqa: E402

#: Fraction of a mate's 31-mers that must be present for the mate to count as witnessed.
POINTS = {"any": 1e-9, "most": 0.5}


def quantiles(v: np.ndarray) -> dict:
    if v.size == 0:
        return {}
    return {"p10": float(np.quantile(v, 0.1)), "p50": float(np.median(v)),
            "p90": float(np.quantile(v, 0.9)), "mean": float(v.mean())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="directory of witness_<run>.npz plus query_names.txt")
    ap.add_argument("--panel", type=Path, default=ROOT / "results/recurrence/panel.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "results/recurrence")
    ap.add_argument("--label", default="pilot",
                    help="names the output, so a pilot's table never overwrites a full "
                         "pass's and the two can be compared instead")
    args = ap.parse_args()

    names = (args.run_dir / "query_names.txt").read_text().split()
    if not names:
        raise SystemExit(f"{args.run_dir}/query_names.txt is empty")
    panel = {r["run"]: r for r in csv.DictReader(open(args.panel), delimiter="\t")}

    files = sorted(args.run_dir.glob("witness_*.npz"))
    if not files:
        raise SystemExit(f"no witness_*.npz in {args.run_dir}")
    runs = [f.name[len("witness_"):-len(".npz")] for f in files]
    missing = [r for r in runs if r not in panel]
    if missing:
        raise SystemExit(f"{missing} are not in the panel table; the subject list and the "
                         f"panel are not the same panel")

    # ---------------------------------------------------------------- load, per mate
    n = len(names)
    frac = {1: np.zeros((len(files), n), dtype=np.float32),
            2: np.zeros((len(files), n), dtype=np.float32)}
    med = {1: np.zeros((len(files), n), dtype=np.float32),
           2: np.zeros((len(files), n), dtype=np.float32)}
    kmers = {1: None, 2: None}
    for i, f in enumerate(files):
        d = np.load(f)
        for mate in (1, 2):
            frac[mate][i] = np.nan_to_num(d[f"frac{mate}"], nan=0.0)
            med[mate][i] = d[f"med{mate}"]
            if kmers[mate] is None:
                kmers[mate] = d[f"kmers{mate}"]
    # A mate with no scorable 31-mer cannot be witnessed or not witnessed; it is excluded
    # rather than counted as absent, which is the difference between a measurement of zero
    # and absence (the repo's own distinction).
    scorable = {mate: np.asarray(kmers[mate]) > 0 for mate in (1, 2)}

    # ---------------------------------------------------------------- piles and depth bins
    pile = np.array([x.split("|")[0] for x in names])
    dbin = np.array([x.split("|")[1] for x in names])
    own_depth = np.array([int(x.split("|")[2]) for x in names])
    # The five seeds are five flow-cell blocks of this library scanned in ONE query file,
    # because the cost of this pass is one scan per subject and is independent of how many
    # queries that scan carries. So the seed travels in the read's name and the spread across
    # seeds is computed here rather than by five runs of the instrument (CLAUDE.md: a
    # headline number reports its spread, and five runs would also be five different passes).
    seed_of = np.array([x.split("|")[4] if len(x.split("|")) > 4 else f"s{args.seed}"
                        for x in names])

    # ---------------------------------------------------------------- the three denominators
    groups = {"runs": [[i] for i in range(len(runs))], "waters": None, "dates": None}
    for key, field in (("waters", "water"), ("dates", "date")):
        buckets = defaultdict(list)
        for i, run in enumerate(runs):
            buckets[panel[run][field]].append(i)
        groups[key] = list(buckets.values())

    out: dict = {"seed": args.seed, "label": args.label,
                 "subjects": {r: {"water": panel[r]["water"], "date": panel[r]["date"],
                                  "read_count": int(panel[r]["read_count"])} for r in runs},
                 "denominators": {k: len(v) for k, v in groups.items()},
                 "cells": {}}

    # The subject-depth control, read before the piles: if occupancy is mostly a statement
    # about how deeply each subject was sequenced, every row below inherits that and the
    # pattern is not about the sewage.
    witnessed_per_subject = {}
    for i, run in enumerate(runs):
        ok = scorable[1]
        witnessed_per_subject[run] = float((frac[1][i][ok] >= POINTS["most"]).mean())
    reads = np.array([int(panel[r]["read_count"]) for r in runs], dtype=float)
    wit = np.array([witnessed_per_subject[r] for r in runs], dtype=float)
    if len(runs) >= 3:
        rank = lambda v: np.argsort(np.argsort(v)).astype(float)  # noqa: E731
        rho = float(np.corrcoef(rank(reads), rank(wit))[0, 1])
    else:
        rho = float("nan")
    out["subject_depth_control"] = {
        "witnessed_fraction_per_subject_mate1_most": witnessed_per_subject,
        "spearman_with_subject_read_count": rho,
        "reading": "occupancy mixes recurrence with how deeply each subject was sequenced; "
                   "a high rho means the per-read occupancy is partly a property of the "
                   "panel rather than of the molecule",
    }

    for mate in (1, 2):
        ok = scorable[mate]
        for point, cut in POINTS.items():
            present = frac[mate] >= cut                     # subjects x queries, boolean
            for gname, gidx in groups.items():
                # A query is present in a group if any of that group's runs witnesses it.
                occ = np.zeros(n, dtype=np.int32)
                for members in gidx:
                    occ += present[members].any(axis=0).astype(np.int32)
                total = len(gidx)
                for p in sorted(set(pile)):
                    for b in sorted(set(dbin[pile == p])):
                        sel = (pile == p) & (dbin == b) & ok
                        if sel.sum() < 50:
                            continue
                        v = occ[sel]
                        key = f"{p}|{b}|mate{mate}|{point}|{gname}"
                        depth_here = med[mate][:, sel]
                        seen = depth_here[present[:, sel]]
                        # Per-seed, so a cell's headline numbers carry their spread across
                        # the five blocks rather than one pooled figure.
                        per_seed = {}
                        for s in sorted(set(seed_of[sel])):
                            m = sel & (seed_of == s)
                            if m.sum() < 20:
                                continue
                            w = occ[m].astype(float)
                            per_seed[s] = {"n": int(m.sum()),
                                           "occupancy_p50": float(np.median(w)),
                                           "occupancy_mean": float(w.mean()),
                                           "frac_in_none": float((w == 0).mean()),
                                           "frac_in_all": float((w == total).mean())}
                        out["cells"][key] = {
                            "per_seed": per_seed,
                            "seed_spread": {
                                field: [round(min(x[field] for x in per_seed.values()), 4),
                                        round(max(x[field] for x in per_seed.values()), 4)]
                                for field in ("occupancy_mean", "frac_in_none",
                                              "frac_in_all")} if per_seed else {},
                            "n": int(sel.sum()),
                            "own_depth_median": float(np.median(own_depth[sel])),
                            "denominator": total,
                            "occupancy": quantiles(v.astype(float)),
                            "frac_in_none": float((v == 0).mean()),
                            "frac_in_all": float((v == total).mean()),
                            "frac_in_middle": float(((v > 0) & (v < total)).mean()),
                            "subject_depth_when_present": quantiles(seen.astype(float))
                                                          if seen.size else {},
                            "occupancy_histogram": np.bincount(
                                v, minlength=total + 1).tolist(),
                        }

    out["manifest"] = {"script": "experiments/recurrence_report.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()},
                       "operating_points": {
                           "presence": {k: f"fraction of a mate's 31-mers >= {v}"
                                        for k, v in POINTS.items()},
                           "k": 31,
                           "piles": "library label table, 40 matching bases",
                       },
                       "subjects_scanned": len(runs)}
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / f"recurrence_{args.label}_s{args.seed}.json"
    path.write_text(json.dumps(out, indent=1) + "\n")

    # ---------------------------------------------------------------- the readable summary
    print(f"{len(runs)} subjects, {out['denominators']} denominators, "
          f"subject-depth rho {rho:+.2f}")
    print(f"{'cell':34s} {'n':>6s} {'own d':>6s} {'occ p50':>8s} {'none':>7s} {'all':>7s}")
    for key, c in out["cells"].items():
        p, b, mate, point, gname = key.split("|")
        if mate != "mate1" or point != "most" or gname != "waters":
            continue
        print(f"{p + ' ' + b:34s} {c['n']:6d} {c['own_depth_median']:6.0f} "
              f"{c['occupancy']['p50']:8.1f} {c['frac_in_none'] * 100:6.1f}% "
              f"{c['frac_in_all'] * 100:6.1f}%")
    print(f"\n{path}")


if __name__ == "__main__":
    main()
