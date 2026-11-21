"""Cross-plant report: how different the backgrounds are, what a union leaves, what grew.

Reads `crossplant_reduce.py`'s packed matrices. Every rate is per mate, at two operating
points (`any`: one 31-mer present; `most`: half or more), per seed, reported as the median
and the min-max over seeds. Subjects are all cut to the same depth, so a difference between
plants is not a difference in library size.

1. Plant x plant: the fraction of plant P's query reads witnessed in a library of plant Q,
   averaged over library pairs at different waters. The diagonal (same plant, different
   water) is the within-plant comparator.
2. Saturation: a held-out library's reads left unwitnessed by a union of k other libraries,
   drawn from its own plant, other plants, or everything except its own water. The Berlin
   piles (identified, rRNA, unplaced and their shuffles) go through the same curves, with
   SRR23973193 as the held-out library.
3. Growth: reads whose 31-mer count in their own library is far above that read's count in
   every panel library, which is reported instead of subtracted. Count ratios are comparable
   because depth is fixed. Library queries come from spots outside the subject subset, so the
   own-library count excludes the read itself and is not size-biased by the draw.
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
#: A read counts as background if both mates are witnessed in this share of the panel.
BACKGROUND_MIN = 0.10


def load(packed: Path):
    z = np.load(packed / "packed.npz")
    names = (packed / "query_names.txt").read_text().split("\n")[: int(z["nq"])]
    subs = [str(s) for s in z["subjects"]]
    nq = int(z["nq"])
    mats = {k: np.unpackbits(z[k], axis=1, count=nq).astype(bool)
            for k in ("any1", "any2", "most1", "most2")}
    meds = {k: z[k] for k in ("med1", "med2")}
    return names, subs, mats, meds, z["library_pairs"]


def parse(names):
    """Index the query rows: library queries by (acc, seed, real|null), Berlin by pile."""
    lib, berlin = defaultdict(list), defaultdict(list)
    for i, n in enumerate(names):
        f = n.split("|")
        seed = int(f[-1][1:])
        if f[0] in ("lib", "shuffled-lib"):
            lib[(f[1], seed, "null" if f[0].startswith("shuffled") else "real")].append(i)
        elif f[0] == "berlin":
            berlin[(f[1], f[2], seed)].append(i)
    return ({k: np.array(v) for k, v in lib.items()},
            {k: np.array(v) for k, v in berlin.items()})


def spread(vals):
    v = [x for x in vals if x == x]
    if not v:
        return "", "", ""
    return f"{np.median(v):.4f}", f"{min(v):.4f}", f"{max(v):.4f}"


def plant_matrix(subs, info, mats, lib):
    plants = sorted({info[s]["key"] for s in subs if s in info})
    rows = []
    for op in ("most1", "most2", "any1", "any2"):
        M = mats[op]
        per_seed = {seed: defaultdict(list) for seed in SEEDS}
        for a_j, a in enumerate(subs):
            if a not in info:
                continue
            for seed in SEEDS:
                idx = lib.get((a, seed, "real"))
                if idx is None:
                    continue
                rate = M[:, idx].mean(axis=1)
                for b_j, b in enumerate(subs):
                    if b == a or b not in info or info[b]["water"] == info[a]["water"]:
                        continue
                    per_seed[seed][(info[a]["key"], info[b]["key"])].append(rate[b_j])
        for p in plants:
            for q in plants:
                vals = [np.mean(per_seed[s][(p, q)]) if per_seed[s][(p, q)] else float("nan")
                        for s in SEEDS]
                n_pairs = len(per_seed[SEEDS[0]][(p, q)])
                med, lo, hi = spread(vals)
                rows.append({"op": op, "query_plant": p, "subject_plant": q,
                             "within": int(p == q), "library_pairs": n_pairs,
                             "witnessed": med, "seed_min": lo, "seed_max": hi})
    return rows


def null_rates(subs, mats, lib):
    out = {}
    for op in mats:
        vals = []
        for seed in SEEDS:
            hits = tot = 0
            for a in subs:
                idx = lib.get((a, seed, "null"))
                if idx is None:
                    continue
                others = [j for j, b in enumerate(subs) if b != a]
                hits += int(mats[op][others][:, idx].sum())
                tot += len(others) * len(idx)
            vals.append(hits / tot if tot else float("nan"))
        out[op] = spread(vals)
    return out


def self_rates(subs, mats, lib):
    out = {}
    for op in mats:
        vals = []
        for seed in SEEDS:
            got = [mats[op][j, lib[(a, seed, "real")]].mean() for j, a in enumerate(subs)
                   if (a, seed, "real") in lib]
            vals.append(float(np.min(got)) if got else float("nan"))
        out[op] = spread(vals)
    return out


def pools(a, subs, info):
    own_water = info[a]["water"] if a in info else None
    key = info[a]["key"] if a in info else None
    ok = [j for j, b in enumerate(subs) if b in info and b != a and info[b]["water"] != own_water]
    return {"own_plant": [j for j in ok if info[subs[j]]["key"] == key],
            "other_plants": [j for j in ok if info[subs[j]]["key"] != key],
            "all": ok}


def curve(M, pool, idx, rng):
    """Leftover fraction after a union of the first k libraries of a random order."""
    order = list(pool)
    rng.shuffle(order)
    cum = np.logical_or.accumulate(M[order][:, idx], axis=0)
    left = 1.0 - cum.mean(axis=1)
    return {k: float(left[k - 1]) for k in KS if k <= len(order)} | {"max": float(left[-1])}


def saturation(subs, info, mats, lib, berlin):
    rows = []
    for op in ("most1", "most2", "any1"):
        M = mats[op]
        acc_curves = defaultdict(lambda: defaultdict(list))
        for a in subs:
            if a not in info:
                continue
            P = pools(a, subs, info)
            for seed in SEEDS:
                idx = lib.get((a, seed, "real"))
                if idx is None:
                    continue
                for mode, pool in P.items():
                    if not pool:
                        continue
                    c = curve(M, pool, idx, random.Random(f"{a}|{seed}|{mode}"))
                    for k, v in c.items():
                        acc_curves[(info[a]["project"], mode, k)][seed].append(v)
                        acc_curves[("ALL", mode, k)][seed].append(v)
        for (proj, mode, k), by_seed in sorted(acc_curves.items(), key=str):
            meds = [float(np.median(v)) for v in by_seed.values()]
            p10 = float(np.median([np.percentile(v, 10) for v in by_seed.values()]))
            p90 = float(np.median([np.percentile(v, 90) for v in by_seed.values()]))
            n = len(next(iter(by_seed.values())))
            med, lo, hi = spread(meds)
            rows.append({"op": op, "heldout": proj, "pool": mode, "k": k,
                         "heldout_libraries": n, "left_median": med, "seed_min": lo,
                         "seed_max": hi, "left_p10": f"{p10:.4f}", "left_p90": f"{p90:.4f}"})
        # Berlin piles: SRR23973193 held out; pools are the 6 Berlin comparator libraries,
        # every other plant, and both.
        bpools = {"berlin_panel": [j for j, b in enumerate(subs) if b in info
                                   and info[b]["project"] == "PRJNA948850"],
                  "other_plants": [j for j, b in enumerate(subs) if b in info
                                   and info[b]["project"] != "PRJNA948850"]}
        bpools["all"] = bpools["berlin_panel"] + bpools["other_plants"]
        piles = sorted({(p, d) for p, d, _ in berlin})
        for pile, dbin in piles:
            for mode, pool in bpools.items():
                per = defaultdict(list)
                for seed in SEEDS:
                    idx = berlin.get((pile, dbin, seed))
                    if idx is None or not pool:
                        continue
                    for k, v in curve(M, pool, idx, random.Random(f"berlin|{seed}|{mode}")).items():
                        per[k].append(v)
                for k, vals in per.items():
                    med, lo, hi = spread(vals)
                    rows.append({"op": op, "heldout": f"berlin:{pile}|{dbin}", "pool": mode,
                                 "k": k, "heldout_libraries": 1, "left_median": med,
                                 "seed_min": lo, "seed_max": hi, "left_p10": "",
                                 "left_p90": ""})
    return rows


def growth(names, subs, info, meds, mats, lib, top: int):
    """Own-library count against the panel's spread for the same read, both mates."""
    out, summary, n_all = [], {}, {"real": 0, "null": 0}
    swaps = {"real": [], "null": []}
    for kind in ("real", "null"):
        ratios = []
        for a_j, a in enumerate(subs):
            if a not in info:
                continue
            P = pools(a, subs, info)["all"]
            for seed in SEEDS:
                idx = lib.get((a, seed, kind))
                if idx is None:
                    continue
                g, own, pan, gs = [], [], [], []
                # Exchangeability control: a random panel library standing in for the own
                # library, scored against the rest of the panel the same way. If the own
                # library were just another panel draw, this is what growth would read.
                swap = random.Random(f"swap|{a}|{seed}").choice(range(len(P)))
                rest = [j for t, j in enumerate(P) if t != swap]
                for m in ("1", "2"):
                    o = meds[f"med{m}"][a_j, idx].astype(float)
                    panel = meds[f"med{m}"][P][:, idx].astype(float)
                    p99 = np.percentile(panel, 99, axis=0)
                    g.append((o + 1) / (p99 + 1))
                    o_s = meds[f"med{m}"][P[swap], idx].astype(float)
                    p99_s = np.percentile(meds[f"med{m}"][rest][:, idx].astype(float), 99, axis=0)
                    gs.append((o_s + 1) / (p99_s + 1))
                    own.append(o)
                    pan.append((np.median(panel, axis=0), p99, panel.max(axis=0)))
                gmin = np.minimum(g[0], g[1])
                # Only reads that ARE background can outgrow it: witnessed (most, both
                # mates) in at least BACKGROUND_MIN of the panel.
                occ = mats["any1"][P][:, idx].sum(axis=0)
                bg = ((mats["most1"][P][:, idx] & mats["most2"][P][:, idx]).mean(axis=0)
                      >= BACKGROUND_MIN)
                ratios.append(gmin[bg])
                swaps[kind].append(np.minimum(gs[0], gs[1])[bg])
                n_all[kind] += len(idx)
                if kind == "real":
                    for t, qi in enumerate(idx):
                        if not bg[t]:
                            continue
                        out.append({"query": names[qi], "library": a,
                                    "plant": info[a]["key"], "date": info[a]["date"],
                                    "seed": seed, "growth": round(float(gmin[t]), 2),
                                    "own1": int(own[0][t]), "own2": int(own[1][t]),
                                    "panel_med1": float(pan[0][0][t]),
                                    "panel_p99_1": float(pan[0][1][t]),
                                    "panel_max1": int(pan[0][2][t]),
                                    "panel_libraries_any1": int(occ[t]),
                                    "panel_size": len(P)})
        for label, parts in ((kind, ratios), (f"{kind}_swap_control", swaps[kind])):
            r = np.concatenate(parts) if parts else np.array([])
            summary[label] = {"queries": n_all[kind], "background_reads": int(r.size),
                              **{f"frac_ge_{t}": round(float((r >= t).mean()), 6)
                                 if r.size else None for t in (2, 5, 10, 100)},
                              "p99": round(float(np.percentile(r, 99)), 2) if r.size else None}
    out.sort(key=lambda d: -d["growth"])
    return out[:top], summary


def composition(sample_dir: Path):
    rows = []
    for p in sorted(sample_dir.glob("sample_*.json")):
        d = json.loads(p.read_text())
        row = {"run": d["acc"], "trim_r2": d["trim_r2"], "spots": d["spots"],
               "mean_len1": d["mean_len"]["r1"], "mean_len2": d["mean_len"]["r2"]}
        for m in ("r1", "r2"):
            c = np.array(d["composition_5p"][m], dtype=float)[:, :4]
            top = c.max(axis=1) / np.maximum(c.sum(axis=1), 1)
            row[f"{m}_maxbase_0_14"] = round(float(top[:15].mean()), 3)
            row[f"{m}_maxbase_20_29"] = round(float(top[20:30].mean()), 3)
        rows.append(row)
    return rows


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
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="downloaded run: packed/, sample/")
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
    names, subs, mats, meds, pairs = load(args.run_dir / "packed")
    lib, berlin = parse(names)
    args.out.mkdir(parents=True, exist_ok=True)
    write(args.out / "plant_matrix.tsv", plant_matrix(subs, info, mats, lib))
    write(args.out / "saturation.tsv", saturation(subs, info, mats, lib, berlin))
    top, gsum = growth(names, subs, info, meds, mats, lib, args.top)
    write(args.out / "growth_top.tsv", top)
    write(args.out / "composition_5p.tsv", composition(args.run_dir / "sample"))
    controls = {"own_library_other_spots_min_over_libraries": self_rates([s for s in subs if s in info], mats, lib),
                "shuffled_in_other_libraries": null_rates([s for s in subs if s in info],
                                                          mats, lib),
                "growth": gsum,
                "subject_pairs_min_max": [int(pairs.min()), int(pairs.max())]}
    (args.out / "controls.json").write_text(json.dumps(controls, indent=2) + "\n")
    manifest = {"script": "experiments/crossplant_report.py", **provenance(__file__),
                "parameters": {"run_dir": str(args.run_dir), "seeds": SEEDS, "ks": KS,
                               "background_min": BACKGROUND_MIN,
                               "top": args.top},
                "subjects": len(subs), "queries": len(names),
                "seconds": round(time.time() - t0, 1)}
    (args.out / "report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(controls, indent=2))


if __name__ == "__main__":
    main()
