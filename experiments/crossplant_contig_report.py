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


def per_contig_fraction(M, rows, j):
    """Witnessed-window fraction of one contig in subject j: both windows of every pair."""
    return (M["w1"][j, rows].sum() + M["w2"][j, rows].sum()) / (2 * len(rows))


def windows_view(mats, op):
    return {"w1": mats[f"{op}1"], "w2": mats[f"{op}2"]}


def plant_table(subs, info, mats, groups, length):
    rows = []
    for op in ("most", "any"):
        M = windows_view(mats, op)
        plants = sorted({info[s]["key"] for s in subs if s in info})
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                for p in plants:
                    cols = [j for j, s in enumerate(subs)
                            if s in info and info[s]["key"] == p and s != QUERY_RUN]
                    if not cols:
                        continue
                    vals = []
                    for seed in SEEDS:
                        g = groups.get((cset, kind, seed), {})
                        if not g:
                            continue
                        f = [per_contig_fraction(M, r, j) for ci, r in g.items()
                             for j in cols]
                        vals.append(np.mean(f) if f else float("nan"))
                    med, lo, hi = spread(vals)
                    n = len(groups.get((cset, kind, SEEDS[0]), {}))
                    rows.append({"op": op, "set": cset, "kind": kind, "subject_plant": p,
                                 "libraries": len(cols), "contigs_per_seed": n,
                                 "witnessed_windows": med, "seed_min": lo, "seed_max": hi})
    return rows


def saturation(subs, info, mats, groups, length):
    """Leftover windows and untouched contigs against union size, by set and stratum."""
    pool = [j for j, s in enumerate(subs) if s in info and s != QUERY_RUN]
    rows = []
    for op in ("most", "any"):
        M = windows_view(mats, op)
        acc = defaultdict(lambda: defaultdict(list))
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                for seed in SEEDS:
                    g = groups.get((cset, kind, seed), {})
                    if not g:
                        continue
                    order = list(pool)
                    random.Random(f"{cset}|{kind}|{seed}").shuffle(order)
                    idx = np.array(sorted({i for r in g.values() for i in r}))
                    cum1 = np.logical_or.accumulate(M["w1"][order][:, idx], axis=0)
                    cum2 = np.logical_or.accumulate(M["w2"][order][:, idx], axis=0)
                    pos = {q: t for t, q in enumerate(idx)}
                    for k in [x for x in KS if x <= len(order)] + [len(order)]:
                        left, untouched = [], []
                        for ci, r in g.items():
                            t = [pos[q] for q in r]
                            seen = int(cum1[k - 1, t].sum() + cum2[k - 1, t].sum())
                            frac = 1.0 - seen / (2 * len(r))
                            left.append(frac)
                            untouched.append(frac == 1.0)
                            s = stratum(length[(cset, ci)])
                            acc[(cset, kind, s, k, "left")][seed].append(frac)
                            acc[(cset, kind, s, k, "untouched")][seed].append(float(frac == 1.0))
                            acc[(cset, kind, s, k, "removed")][seed].append(
                                float(frac <= 1 - REMOVED_AT))
                        for field, v in (("left", left), ("untouched", untouched)):
                            acc[(cset, kind, "all", k, field)][seed] += [float(x) for x in v]
                        acc[(cset, kind, "all", k, "removed")][seed] += [
                            float(x <= 1 - REMOVED_AT) for x in left]
        # Sorted on k as a number, not as text: a committed table that reads
        # 1, 16, 2, 32 is a table nobody checks against the prose.
        for key, by_seed in sorted(acc.items(),
                                   key=lambda kv: (kv[0][0], kv[0][1], kv[0][2],
                                                   kv[0][4], kv[0][3])):
            cset, kind, s, k, field = key
            meds = [np.mean(v) for v in by_seed.values()]
            med, lo, hi = spread(meds)
            rows.append({"op": op, "set": cset, "kind": kind, "stratum": s, "k": k,
                         "field": field, "contigs": len(next(iter(by_seed.values()))),
                         "value": med, "seed_min": lo, "seed_max": hi})
    return rows


def growth(subs, info, mats, meds, groups, length, top: int):
    pool = [j for j, s in enumerate(subs) if s in info and s != QUERY_RUN]
    own = next((j for j, s in enumerate(subs) if s == QUERY_RUN), None)
    if own is None:
        return [], {"note": f"{QUERY_RUN} was not among the subjects; growth needs it"}
    out, summary = [], {}
    M = windows_view(mats, "most")
    for cset in ("unident", "ident"):
        for kind in ("real", "null"):
            ratios, swapped = [], []
            for seed in SEEDS:
                g = groups.get((cset, kind, seed), {})
                if not g:
                    continue
                swap = random.Random(f"swap|{cset}|{seed}").choice(range(len(pool)))
                rest = [j for t, j in enumerate(pool) if t != swap]
                for ci, r in g.items():
                    o = np.concatenate([meds["med1"][own, r], meds["med2"][own, r]]).astype(float)
                    panel = np.concatenate([meds["med1"][pool][:, r],
                                            meds["med2"][pool][:, r]], axis=1).astype(float)
                    p99 = np.percentile(panel, 99, axis=0)
                    ratio = float(np.median((o + 1) / (p99 + 1)))
                    o_s = np.concatenate([meds["med1"][pool[swap], r],
                                          meds["med2"][pool[swap], r]]).astype(float)
                    rest_p99 = np.percentile(np.concatenate(
                        [meds["med1"][rest][:, r], meds["med2"][rest][:, r]], axis=1
                    ).astype(float), 99, axis=0)
                    swapped.append(float(np.median((o_s + 1) / (rest_p99 + 1))))
                    ratios.append(ratio)
                    occ = float(np.mean([per_contig_fraction(M, r, j) for j in pool]))
                    if kind == "real":
                        out.append({"set": cset, "contig": ci,
                                    "length": length[(cset, ci)], "seed": seed,
                                    "windows": 2 * len(r), "growth": round(ratio, 2),
                                    "own_median_depth": float(np.median(o)),
                                    "panel_p99_median": float(np.median(p99)),
                                    "panel_witnessed_windows": round(occ, 4)})
            for label, v in ((kind, ratios), (f"{kind}_swap_control", swapped)):
                a = np.array(v)
                summary[f"{cset}|{label}"] = {
                    "contigs": int(a.size),
                    **{f"frac_ge_{t}": round(float((a >= t).mean()), 6) if a.size else None
                       for t in (2, 5, 10, 100)},
                    "p99": round(float(np.percentile(a, 99)), 2) if a.size else None}
    out.sort(key=lambda d: -d["growth"])
    return out[:top], summary


def controls(subs, info, mats, groups):
    """The floor and the ceiling of this arm, both measured in the same pass."""
    out = {}
    own = next((j for j, s in enumerate(subs) if s == QUERY_RUN), None)
    for op in ("most", "any"):
        M = windows_view(mats, op)
        for cset in ("unident", "ident"):
            for kind in ("real", "null"):
                if own is not None:
                    vals = []
                    for seed in SEEDS:
                        g = groups.get((cset, kind, seed), {})
                        if g:
                            vals.append(np.mean([per_contig_fraction(M, r, own)
                                                 for r in g.values()]))
                    out[f"own_assembly_library|{op}|{cset}|{kind}"] = spread(vals)
                vals = []
                for seed in SEEDS:
                    g = groups.get((cset, kind, seed), {})
                    cols = [j for j, s in enumerate(subs) if s in info and s != QUERY_RUN]
                    if g and cols:
                        vals.append(np.mean([per_contig_fraction(M, r, j)
                                             for r in g.values() for j in cols]))
                out[f"panel_mean|{op}|{cset}|{kind}"] = spread(vals)
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
