"""Cross-plant contig report: what a panel of other plants' water removes from tier 2.

Reads `crossplant_reduce.py`'s packed matrices for a contig-mode run. The unit is a contig,
not a query row: a contig's windows are aggregated into one witnessed-window fraction per
subject, which is what makes this arm stronger than the read arm - a read gives one bit per
library, a contig gives as many as it has windows.

Every number is per set (`unident` is tier 2, `ident` is the known-background comparator
from the same assembly), per seed, median and min-max over the five seeds, at two operating
points (`any`: one 31-mer of a window present; `most`: half or more). Both windows of a
query pair are treated as two independent windows of the same molecule rather than as two
mates, because here they are the same length and the same kit - the read arm's mate
asymmetry does not exist in this material.

1. Plant x set: the witnessed-window fraction of a contig in a library of plant Q, averaged
   over that plant's libraries. The identified set is the ceiling and the shuffles the floor,
   both drawn in the same pass.
2. Saturation: the fraction of a contig's windows still unwitnessed after a union of k panel
   libraries, as a curve in k, and the share of contigs left wholly untouched. This is the
   subtraction question: does adding plants keep removing tier 2, or does it stop.
3. Length strata: the same readouts split by contig length, because length drives both how
   many windows a contig offers and how findable it is, and matching on it would hide the
   question rather than answer it.
4. Growth: a contig whose own-library window depth stands far above the panel's p99 for the
   same windows is reported rather than subtracted, with a swap control in which a panel
   library stands in for the own library.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from subtract_known import provenance  # noqa: E402

QUERY_RUN = "SRR23973193"
KS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
SEEDS = [1, 2, 3, 4, 5]
#: Length strata in nt. A contig's window count, and so its resolution, follows its length.
STRATA = [(180, 270), (270, 450), (450, 900), (900, 10 ** 9)]
#: A contig counts as removed when this share of its windows is witnessed in the union.
REMOVED_AT = 0.5


def load(packed: Path):
    z = np.load(packed / "packed.npz")
    nq = int(z["nq"])
    names = (packed / "query_names.txt").read_text().split("\n")[:nq]
    subs = [str(s) for s in z["subjects"]]
    mats = {k: np.unpackbits(z[k], axis=1, count=nq).astype(bool)
            for k in ("any1", "any2", "most1", "most2")}
    meds = {k: z[k] for k in ("med1", "med2")}
    return names, subs, mats, meds


def parse(names):
    """(set, kind, seed) -> {contig index -> [query rows]}, plus each contig's length."""
    out = defaultdict(lambda: defaultdict(list))
    length = {}
    for i, n in enumerate(names):
        f = n.split("|")
        if f[0] not in ("contig", "shuffled-contig"):
            continue
        kind = "null" if f[0].startswith("shuffled") else "real"
        cset, ci, ln, seed = f[1], int(f[2]), int(f[3][1:]), int(f[5][1:])
        out[(cset, kind, seed)][ci].append(i)
        length[(cset, ci)] = ln
    return out, length


def spread(vals):
    v = [float(x) for x in vals if x == x]
    if not v:
        return "", "", ""
    return f"{np.median(v):.4f}", f"{min(v):.4f}", f"{max(v):.4f}"


def stratum(ln: int) -> str:
    for lo, hi in STRATA:
        if lo <= ln < hi:
            return f"{lo}-{hi if hi < 10 ** 9 else 'max'}"
    return "short"


def group_blocks(groups, key):
    """Query columns of one (set, kind, seed) group, grouped by contig, plus block starts.

    The per-contig readouts are aggregations over a contig's window columns, and doing them
    with a Python loop over (contig, subject) pairs is 1.7 billion iterations for this panel -
    the first version of this report ran past forty minutes and was killed. Laying the columns
    out contiguously per contig lets `np.add.reduceat` do all of them in one call per matrix.
    """
    items = sorted(groups.get(key, {}).items())
    if not items:
        return None
    cols, starts, ids, widths = [], [], [], []
    for ci, rows in items:
        starts.append(len(cols))
        cols.extend(rows)
        ids.append(ci)
        widths.append(len(rows))
    return (np.array(cols), np.array(starts), ids, np.array(widths, dtype=float))


def contig_fractions(mats, op, cols, starts, widths):
    """(subjects x contigs) witnessed-window fraction: both windows of every pair."""
    a = np.add.reduceat(mats[f"{op}1"][:, cols], starts, axis=1)
    b = np.add.reduceat(mats[f"{op}2"][:, cols], starts, axis=1)
    return (a + b) / (2.0 * widths)


def plant_table(subs, info, mats, groups, length):
    rows = []
    plants = sorted({info[x]["key"] for x in subs if x in info})
    cols_by_plant = {p: [j for j, x in enumerate(subs)
                         if x in info and info[x]["key"] == p and x != QUERY_RUN]
                     for p in plants}
    for op in ("most", "any"):
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                per_seed = {}
                n_contigs = 0
                for seed in SEEDS:
                    blocks = group_blocks(groups, (cset, kind, seed))
                    if blocks is None:
                        continue
                    cols, starts, ids, widths = blocks
                    n_contigs = len(ids)
                    F = contig_fractions(mats, op, cols, starts, widths)
                    per_seed[seed] = F
                for p in plants:
                    cj = cols_by_plant[p]
                    if not cj:
                        continue
                    vals = [float(F[cj].mean()) for F in per_seed.values()]
                    med, lo, hi = spread(vals)
                    rows.append({"op": op, "set": cset, "kind": kind, "subject_plant": p,
                                 "libraries": len(cj), "contigs_per_seed": n_contigs,
                                 "witnessed_windows": med, "seed_min": lo, "seed_max": hi})
    return rows


def saturation(subs, info, mats, groups, length):
    """Leftover windows, untouched contigs and removed contigs against union size."""
    pool = [j for j, x in enumerate(subs) if x in info and x != QUERY_RUN]
    rows = []
    for op in ("most", "any"):
        acc = defaultdict(lambda: defaultdict(list))
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                for seed in SEEDS:
                    blocks = group_blocks(groups, (cset, kind, seed))
                    if blocks is None:
                        continue
                    cols, starts, ids, widths = blocks
                    strata = np.array([stratum(length[(cset, ci)]) for ci in ids])
                    order = list(pool)
                    random.Random(f"{cset}|{kind}|{seed}").shuffle(order)
                    c1 = np.logical_or.accumulate(mats[f"{op}1"][order][:, cols], axis=0)
                    c2 = np.logical_or.accumulate(mats[f"{op}2"][order][:, cols], axis=0)
                    ks = [x for x in KS if x <= len(order)] + [len(order)]
                    for k in ks:
                        seen = (np.add.reduceat(c1[k - 1], starts)
                                + np.add.reduceat(c2[k - 1], starts))
                        left = 1.0 - seen / (2.0 * widths)
                        for s_ in np.unique(np.append(strata, "all")):
                            m = slice(None) if s_ == "all" else (strata == s_)
                            v = left[m]
                            if not len(v):
                                continue
                            acc[(cset, kind, str(s_), k, "left")][seed].append(float(v.mean()))
                            acc[(cset, kind, str(s_), k, "untouched")][seed].append(
                                float((v == 1.0).mean()))
                            acc[(cset, kind, str(s_), k, "removed")][seed].append(
                                float((v <= 1 - REMOVED_AT).mean()))
                            acc[(cset, kind, str(s_), k, "contigs")][seed].append(float(len(v)))
        # Sorted on k as a number, not as text: a committed table that reads
        # 1, 16, 2, 32 is a table nobody checks against the prose.
        for key, by_seed in sorted(acc.items(),
                                   key=lambda kv: (kv[0][0], kv[0][1], kv[0][2],
                                                   kv[0][4], kv[0][3])):
            cset, kind, s_, k, field = key
            if field == "contigs":
                continue
            meds = [float(np.mean(v)) for v in by_seed.values()]
            med, lo, hi = spread(meds)
            n = acc[(cset, kind, s_, k, "contigs")]
            rows.append({"op": op, "set": cset, "kind": kind, "stratum": s_, "k": k,
                         "field": field,
                         "contigs": int(np.mean([v[0] for v in n.values()])) if n else 0,
                         "value": med, "seed_min": lo, "seed_max": hi})
    return rows


def _per_contig_median(v, starts, widths, w):
    """Median over a contig's own window columns, in both mates."""
    out = np.empty(len(starts))
    for t, (st, wd) in enumerate(zip(starts, widths)):
        c = np.arange(int(st), int(st + wd))
        out[t] = np.median(v[np.r_[c, c + w]])
    return out


def growth(subs, info, mats, meds, groups, length, top: int):
    """Own-assembly depth against the panel's p99 over the same windows, per contig.

    Imo's third question: a background that stops being a background gets reported rather
    than subtracted. The swap control answers the only question that makes the ratio
    legible - what growth reads when the query library is just another panel draw - and the
    two mates are kept apart because they are different lengths and map at different rates,
    so a pooled ratio averages two leans (CLAUDE.md).
    """
    pool = [j for j, x in enumerate(subs) if x in info and x != QUERY_RUN]
    own = next((j for j, x in enumerate(subs) if x == QUERY_RUN), None)
    if own is None:
        return [], {"note": f"{QUERY_RUN} was not among the subjects; growth needs it"}
    out, summary = [], {}
    for cset in ("unident", "ident"):
        for kind in ("real", "null"):
            ratios, swapped, carried, by_mate = [], [], [], {1: [], 2: []}
            for seed in SEEDS:
                blocks = group_blocks(groups, (cset, kind, seed))
                if blocks is None:
                    continue
                cols, starts, ids, widths = blocks
                w = len(cols)
                # One panel library stands in for the own library and is scored against the
                # rest: the ratio a library with no local bloom returns.
                swap = random.Random(f"swap|{cset}|{seed}").randrange(len(pool))
                rest = [j for t, j in enumerate(pool) if t != swap]
                own_d = np.concatenate([meds["med1"][own, cols],
                                        meds["med2"][own, cols]]).astype(float)
                p99 = np.percentile(np.concatenate(
                    [meds["med1"][pool][:, cols], meds["med2"][pool][:, cols]],
                    axis=1).astype(float), 99, axis=0)
                g = (own_d + 1) / (p99 + 1)
                swap_d = np.concatenate([meds["med1"][pool[swap], cols],
                                         meds["med2"][pool[swap], cols]]).astype(float)
                p99_s = np.percentile(np.concatenate(
                    [meds["med1"][rest][:, cols], meds["med2"][rest][:, cols]],
                    axis=1).astype(float), 99, axis=0)
                gs = (swap_d + 1) / (p99_s + 1)
                per = _per_contig_median(g, starts, widths, w)
                per_s = _per_contig_median(gs, starts, widths, w)
                ratios.append(per)
                swapped.append(per_s)
                carried.append(_per_contig_median(p99, starts, widths, w) >= 1.0)
                for m, sl in ((1, slice(0, w)), (2, slice(w, 2 * w))):
                    gm, cm = g[sl], np.arange(w)
                    by_mate[m].append(np.array([
                        np.median(gm[int(st):int(st + wd)])
                        for st, wd in zip(starts, widths)]))
                if kind == "real":
                    F = contig_fractions(mats, "most", cols, starts, widths)
                    occ = F[pool].mean(axis=0)
                    own_med = np.add.reduceat(own_d[:w], starts) / widths
                    p99_med = np.add.reduceat(p99[:w], starts) / widths
                    for t, ci in enumerate(ids):
                        c = slice(int(starts[t]), int(starts[t] + widths[t]))
                        out.append({"set": cset, "contig": ci,
                                    "length": length[(cset, ci)], "seed": seed,
                                    "windows": int(2 * widths[t]),
                                    "growth": round(float(per[t]), 2),
                                    "growth_mate1": round(float(by_mate[1][-1][t]), 2),
                                    "growth_mate2": round(float(by_mate[2][-1][t]), 2),
                                    "own_mean_depth_mate1": round(float(own_med[t]), 1),
                                    "panel_p99_mean_mate1": round(float(p99_med[t]), 1),
                                    "panel_witnessed_windows": round(float(occ[t]), 4),
                                    "panel_carries_contig": int(carried[-1][t]),
                                    "own_count_capped": int(
                                        (meds["med1"][own, cols][c] >= 65535).any())})

            def agg(parts, mask=None):
                a = np.concatenate(parts) if parts else np.array([])
                if mask is not None:
                    a = a[np.concatenate(mask)]
                return {"contigs": int(a.size),
                        **{f"frac_ge_{t}": round(float((a >= t).mean()), 6) if a.size else None
                           for t in (2, 5, 10, 100)},
                        "p99": round(float(np.percentile(a, 99)), 2) if a.size else None}

            summary[f"{cset}|{kind}"] = agg(ratios)
            summary[f"{cset}|{kind}_swap_control"] = agg(swapped)
            summary[f"{cset}|{kind}_mate1"] = agg(by_mate[1])
            summary[f"{cset}|{kind}_mate2"] = agg(by_mate[2])
            # Where the panel carries no depth at all the ratio's denominator is the +1
            # alone, so "growth" there is own depth wearing a ratio's clothes - and that is
            # most of the top of the list. Split it: carried contigs give a real ratio,
            # uncarried ones are an absence result and belong with the occupancy numbers,
            # not among the blooms. The swap control is split the same way or the comparison
            # moves with the stratum.
            summary[f"{cset}|{kind}_panel_carried"] = agg(ratios, carried)
            summary[f"{cset}|{kind}_panel_carried_swap_control"] = agg(swapped, carried)
            summary[f"{cset}|{kind}_panel_absent"] = agg(ratios, [~m for m in carried])
            summary[f"{cset}|{kind}_panel_absent_frac"] = round(
                float(np.concatenate([~m for m in carried]).mean()), 6) if carried else None
    out.sort(key=lambda d: -d["growth"])
    return out[:top], summary


def controls(subs, info, mats, groups):
    """The floor and the ceiling of this arm, both measured in the same pass."""
    out = {}
    own = next((j for j, x in enumerate(subs) if x == QUERY_RUN), None)
    panel = [j for j, x in enumerate(subs) if x in info and x != QUERY_RUN]
    for op in ("most", "any"):
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                o_vals, p_vals = [], []
                for seed in SEEDS:
                    blocks = group_blocks(groups, (cset, kind, seed))
                    if blocks is None:
                        continue
                    cols, starts, ids, widths = blocks
                    F = contig_fractions(mats, op, cols, starts, widths)
                    if own is not None:
                        o_vals.append(float(F[own].mean()))
                    if panel:
                        p_vals.append(float(F[panel].mean()))
                if own is not None:
                    out[f"own_assembly_library|{op}|{cset}|{kind}"] = spread(o_vals)
                out[f"panel_mean|{op}|{cset}|{kind}"] = spread(p_vals)
    return out


def write(path: Path, rows):
    if not rows:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="unused; seeds 1-5 are all read")
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--panel", type=Path, default=Path("results/crossplant/panel.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/crossplant"))
    ap.add_argument("--top", type=int, default=500)
    args = ap.parse_args()
    t0 = time.time()
    info = {}
    for r in csv.DictReader(args.panel.open(), delimiter="\t"):
        plant = r["plant"] or "untreated"
        info[r["run"]] = {"project": r["project"], "key": f"{r['project']}:{plant}",
                          "water": r["water"], "date": r["date"]}
    names, subs, mats, meds = load(args.run_dir / "packed")
    groups, length = parse(names)
    if not groups:
        raise SystemExit("no contig queries in this run; is it a reads-mode run?")
    args.out.mkdir(parents=True, exist_ok=True)
    write(args.out / "contig_plants.tsv", plant_table(subs, info, mats, groups, length))
    write(args.out / "contig_saturation.tsv", saturation(subs, info, mats, groups, length))
    top, gsum = growth(subs, info, mats, meds, groups, length, args.top)
    write(args.out / "contig_growth_top.tsv", top)
    ctl = {"rates": controls(subs, info, mats, groups), "growth": gsum,
           "contigs_per_seed": {f"{k[0]}|{k[1]}|s{k[2]}": len(v)
                                for k, v in sorted(groups.items(), key=str)},
           "subjects": len(subs), "panel_libraries": len(
               [s for s in subs if s in info and s != QUERY_RUN])}
    (args.out / "contig_controls.json").write_text(json.dumps(ctl, indent=2) + "\n")
    manifest = {"script": "experiments/crossplant_contig_report.py", **provenance(__file__),
                "parameters": {"run_dir": str(args.run_dir), "seeds": SEEDS, "ks": KS,
                               "strata": STRATA, "removed_at": REMOVED_AT, "top": args.top},
                "queries": len(names), "subjects": len(subs),
                "seconds": round(time.time() - t0, 1)}
    (args.out / "contig_report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(ctl["rates"], indent=2))


if __name__ == "__main__":
    main()
