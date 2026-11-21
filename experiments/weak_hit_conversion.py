#!/usr/bin/env python3
"""Is there headroom in turning weak protein hits into identifications?

WHY THIS EXISTS
---------------
The protein-index ladder established that bigger databases buy hits rather than
identifications: across four cumulative indexes the residue's "any hit at all" climbed to
near 100% while identifications at bitscore >= 50 stayed near 1% and were still tripling.
The standing reading of that was **headroom**: almost every residue read already has a
protein hit, so the prize is converting those weak hits into calls.

This tests that reading, and it needs two things the ladder did not have in the same pass:
a **shuffled null**, and a query set that includes the same molecule at two lengths.

It is also a correction. `results/profile_reach/` searched with DIAMOND's *default*
`--evalue 0.001`, so no weak hit was ever reported and its "any hit" columns were all
already significant. Nothing published from that run depended on a weak-hit population —
the reach and merged-versus-mate results are all at or above significance — but the weak
hit question could not be asked of it at all, and that is why this search is re-run
permissively at `--evalue 100` over the identical query set.

THE TWO QUESTIONS
-----------------
1. **Is the weak-hit band distinguishable from noise?** The same fields on real library
   reads, on the holdout at each divergence rung, and on composition-matched shuffles. If
   the band reads the same on shuffles as on real sequence, it is the chance band and
   there is nothing in it to convert. This project has already learned exactly that lesson
   for nucleotide BLAST, where the chance band read identically on residue and
   already-identified reads and a claim was withdrawn for it. Here it is measured for
   translated search against a protein index, with the null in the same pass.

2. **Does merging the mates convert weak hits into identifications?** Among pairs whose
   mate 1 produced only a weak hit, the fraction whose merged read produces a significant
   one. **The shuffle row decides whether this is a gain or arithmetic**: a longer query
   accumulates more bitscore, so if noise converts at the same rate the conversion means
   nothing. An E-value is the right readout here precisely because it penalises query
   length, so a chance hit cannot be promoted by making the query longer.

READ THE CONTROL FIRST.

Usage: python3 experiments/weak_hit_conversion.py --out results/weak_hits
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
SIG = 1e-3          # DIAMOND's own default significance cutoff, used as the call line
PERMISSIVE = 100.0  # the reporting cutoff, so a weak-hit population exists at all


def parse(out: Path) -> dict:
    best: dict[str, dict] = {}
    with open(out) as fh:
        for line in fh:
            q, subj, pid, ln, ev, bs = line.rstrip("\n").split("\t")
            score = float(bs)
            if q not in best or score > best[q]["bitscore"]:
                best[q] = {"subject": subj, "pident": float(pid), "length": int(ln),
                           "evalue": float(ev), "bitscore": score}
    return best


def run_diamond(db: Path, qfa: Path, out: Path, threads: int, log) -> float:
    if out.exists():
        log(f"reusing {out.name}")
        return 0.0
    dmnd = REPO / "bin_local" / "diamond"
    t0 = time.time()
    subprocess.run([str(dmnd), "blastx", "-d", str(db), "-q", str(qfa), "-o", str(out),
                    "--ultra-sensitive", "--max-target-seqs", "1",
                    "--evalue", str(PERMISSIVE), "--threads", str(threads), "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue",
                    "bitscore"], check=True)
    el = time.time() - t0
    log(f"diamond blastx --evalue {PERMISSIVE:g}: {el:.0f}s")
    return el


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/weak_hits")
    ap.add_argument("--work", default="results/profile_reach/work")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    outdir = REPO / args.out
    outdir.mkdir(parents=True, exist_ok=True)
    work = REPO / args.work
    logfh = open(outdir / "run.log", "a")

    def log(m):
        line = f"[{time.strftime('%H:%M:%S')}] {m}"
        print(line, flush=True)
        logfh.write(line + "\n")
        logfh.flush()

    qfa = work / "queries.fa"
    db = work / "prot.dmnd"
    for p in (qfa, db):
        if not p.exists():
            log(f"MISSING INPUT: {p} — run experiments/profile_reach.py first")
            return 2

    t_start = time.time()
    tsv = work / "diamond_permissive.tsv"
    elapsed = run_diamond(db, qfa, tsv, args.threads, log)
    best = parse(tsv)

    qlen, name = {}, None
    for line in open(qfa):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            qlen[name] = len(line.strip())
    log(f"{len(qlen)} queries, {len(best)} with a hit at E<={PERMISSIVE:g}")

    # ---- question 1: is the weak band distinguishable from noise? per seed, per mate
    bands = defaultdict(lambda: defaultdict(list))
    for q in qlen:
        seed, pile, rung, cond, _ = q.split("|")
        bands[(pile, rung, cond)][seed].append(q)

    rows = []
    for (pile, rung, cond), per_seed in bands.items():
        for seed, qs in per_seed.items():
            hits = [best[q] for q in qs if q in best]
            n = len(qs)
            row = {"seed": seed, "pile": pile, "rung": rung, "condition": cond, "n": n,
                   "mean_query_len": float(np.mean([qlen[q] for q in qs])),
                   "any_hit": len(hits) / n,
                   "sig": sum(1 for h in hits if h["evalue"] <= SIG) / n}
            if hits:
                row["median_evalue_of_hits"] = float(np.median([h["evalue"] for h in hits]))
                row["median_bitscore_of_hits"] = float(np.median([h["bitscore"] for h in hits]))
                row["median_pident_of_hits"] = float(np.median([h["pident"] for h in hits]))
            weak = [h for h in hits if h["evalue"] > SIG]
            row["weak_frac"] = len(weak) / n
            if weak:
                row["median_bitscore_of_weak"] = float(np.median([h["bitscore"] for h in weak]))
                row["median_pident_of_weak"] = float(np.median([h["pident"] for h in weak]))
            rows.append(row)

    # ---- question 2: conversion, paired within a molecule
    allq = defaultdict(dict)
    for q in qlen:
        seed, pile, rung, cond, pid = q.split("|")
        allq[(seed, pile, rung, pid)][cond] = q
    conv_rows = []
    groups = defaultdict(list)
    for key, c in allq.items():
        if "merged" in c and "m1" in c:
            groups[(key[0], key[1], key[2])].append(c)
    for (seed, pile, rung), items in groups.items():
        weak = [c for c in items
                if (h := best.get(c["m1"])) is not None and h["evalue"] > SIG]
        already = [c for c in items
                   if (h := best.get(c["m1"])) is not None and h["evalue"] <= SIG]
        r = {"seed": seed, "pile": pile, "rung": rung,
             "n_pairs": len(items), "n_weak_m1": len(weak), "n_sig_m1": len(already)}
        if weak:
            r["converted"] = sum(
                1 for c in weak
                if (h := best.get(c["merged"])) is not None and h["evalue"] <= SIG) / len(weak)
        if already:
            r["kept"] = sum(
                1 for c in already
                if (h := best.get(c["merged"])) is not None and h["evalue"] <= SIG) / len(already)
        # the length-arithmetic null: does a merged read score more than its mate's
        # bitscore scaled by the length ratio? both-hit pairs only
        ratios = []
        for c in items:
            a, b = best.get(c["m1"]), best.get(c["merged"])
            if a and b:
                pred = a["bitscore"] * qlen[c["merged"]] / qlen[c["m1"]]
                if pred > 0:
                    ratios.append(b["bitscore"] / pred)
        if ratios:
            r["bitscore_observed_over_length_prediction"] = float(np.median(ratios))
            r["n_both_hit"] = len(ratios)
        # The decisive control for results/profile_reach/ §2, which measured the merged
        # gain at a BITSCORE floor. An E-value penalises a longer query, so if the gain
        # were a longer query clearing a fixed floor it would shrink here. Reported as a
        # curve over thresholds, within the same overlapping pairs.
        for t, tag in ((1e-3, "1e-3"), (1e-5, "1e-5"), (1e-10, "1e-10")):
            for cond in ("merged", "m1", "m2"):
                have = [c for c in items if cond in c]
                if not have:
                    continue
                r[f"sig_E{tag}_{cond}"] = sum(
                    1 for c in have
                    if (h := best.get(c[cond])) is not None and h["evalue"] <= t
                ) / len(have)
        conv_rows.append(r)

    json.dump(rows, open(outdir / "bands.json", "w"), indent=1)
    json.dump(conv_rows, open(outdir / "conversion.json", "w"), indent=1)

    man = {
        "script": "experiments/weak_hit_conversion.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "git_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                  capture_output=True, text=True).stdout.strip()
        + ("-dirty" if subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                                      capture_output=True, text=True).stdout.strip()
           else ""),
        "seeds": sorted({r["seed"] for r in rows}),
        "params": {
            "significance_cutoff_evalue": SIG,
            "reporting_cutoff_evalue": PERMISSIVE,
            "diamond_mode": "--ultra-sensitive --max-target-seqs 1",
            "queries": str(qfa.relative_to(REPO)),
            "query_count": len(qlen),
            "note": "the query set is byte-identical to results/profile_reach/'s; only "
                    "the reporting E-value cutoff differs, so the two runs are the same "
                    "instrument at two operating points",
        },
        "timings_s": {"diamond": elapsed, "total": time.time() - t_start},
    }
    json.dump(man, open(outdir / "manifest.json", "w"), indent=1)
    log(f"done in {time.time() - t_start:.0f}s -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
