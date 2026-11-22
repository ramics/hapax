#!/usr/bin/env python3
"""Tables for results/holdout_cascade/REPORT.md, computed from the committed run files.

Separate from the scripts that produce those files, for the reason the repo already has a
rule about: a figure typed into prose is a figure without provenance, and this project has
had to withdraw claims whose text drifted from its outputs. Everything the report quotes is
printed here from `holdout_cascade_s1.json`, the five `run_s*/` directories and
`contig_identification.json`, so re-running this is how the report is checked.

THE TABLE THIS EXISTS FOR. The question is whether the withheld genome is anywhere a reader
would land. It is not in the residue - the cascade classified its contig as identified and
subtracted it - so the interesting distribution is over the pile nobody looks at: for every
contig the cascade REMOVED, how far was the relative that removed it? If removals are
overwhelmingly near-identical and the withheld genome sits alone in a sparse tail, then a
distance-with-coverage readout on the removed pile would have surfaced it, and the report
can say what would justify a flag. If that tail is crowded, it would not, and the report
says that instead. Which of those is true is measured here and not chosen.
"""

from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BANDS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95), (95, 99), (99, 100.01)]


def band_of(x: float) -> str:
    for lo, hi in BANDS:
        if lo <= x < hi:
            return f"{lo}-{hi if hi <= 100 else 100}"
    return "?"


def pct(x) -> str:
    return "n/a" if x is None else f"{100 * x:.2f}%"


def minmax(vals) -> str:
    """Five seeds give an honest range and a useless variance, so only the range is shown."""
    v = [x for x in vals if x is not None]
    if not v:
        return "not measured"
    return f"{100 * min(v):.2f}-{100 * max(v):.2f}%" if max(v) <= 1.0 else \
        f"{min(v):.4g}-{max(v):.4g}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=ROOT / "results/holdout_cascade")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    a = ap.parse_args()
    D = a.dir

    sc = json.loads((D / "holdout_cascade_s1.json").read_text())
    ci = json.loads((D / "contig_identification.json").read_text())

    print("=" * 78)
    print("1. THE CASCADE, FIVE SEEDS (each 200,000 reads; one read pool, so the spread")
    print("   below is subsampling noise and nothing else)")
    print("=" * 78)
    runs, tiers_by_seed = {}, {}
    for s in a.seeds:
        p = D / f"run_s{s}" / "tiers.json"
        if p.exists():
            tiers_by_seed[s] = json.loads(p.read_text())
        m = D / f"run_s{s}" / "manifest.json"
        if m.exists():
            runs[s] = json.loads(m.read_text())
    if tiers_by_seed:
        order = [c["tier"] for c in next(iter(tiers_by_seed.values()))["curve"]]
        print(f"{'tier':<42}{'alone':>10}{'marginal':>12}{'cumulative':>14}")
        for i, t in enumerate(order):
            al, mg, cu = [], [], []
            for d in tiers_by_seed.values():
                c = next((x for x in d["curve"] if x["tier"] == t), None)
                if c:
                    al.append(c.get("tier_alone_removed"))
                    mg.append(c.get("marginal_removed"))
                    cu.append(c.get("cumulative_removed"))
            print(f"{t:<42}{minmax(al):>10}{minmax(mg):>12}{minmax(cu):>14}")
        print()
        for label, key in (("identified (tier 1)", "union_removed"),
                           ("residue", None), ("tier 3 (unplaced by anything)", None)):
            if key:
                print(f"  {label:<34}{minmax([d[key] for d in tiers_by_seed.values()])}")
        res = [d["residue_pairs"] / d["pairs"] for d in tiers_by_seed.values()]
        t3 = [d["tier3_pairs"] / d["pairs"] for d in tiers_by_seed.values()]
        t2 = [(d["residue_pairs"] - d["tier3_pairs"]) / d["pairs"]
              for d in tiers_by_seed.values()]
        print(f"  {'residue':<34}{minmax(res)}")
        print(f"  {'tier 2 (assembled, unidentified)':<34}{minmax(t2)}")
        print(f"  {'tier 3 (unplaced by anything)':<34}{minmax(t3)}")
        print("\n  stages, per seed:")
        for s, m in runs.items():
            st_ = {x["stage"]: x["status"] for x in m["stages"]}
            print(f"    seed {s}: complete={m['complete']}  "
                  + " ".join(f"{k}={v}" for k, v in st_.items()))

    print()
    print("=" * 78)
    print("2. THE ASSEMBLY SPLIT (whole library, one assembly, deterministic)")
    print("=" * 78)
    print(f"  contigs {ci['contigs']:,}  bases {ci['bases']:,}")
    print(f"  identified {pct(ci['identified_frac'])} of contigs, "
          f"{pct(ci['identified_bases_frac'])} of bases")
    print(f"  operating point: {json.dumps(ci['operating_point'])}")
    print(f"  identified by protein only: {pct(ci['route_protein_only'])}; "
          f"by nucleotide only: {pct(ci['route_nucleotide_only'])}")
    print("  Kraken2 confidence sweep (fraction of contigs identified):")
    for k, v in ci["kraken2_confidence_sweep"].items():
        print(f"    {k}: {pct(v)}")
    print("  protein bitscore sweep:")
    for db, sw in ci["protein_bitscore_sweep"].items():
        print(f"    {db}: " + "  ".join(f"{k}={pct(v)}" for k, v in sw.items()))

    print()
    print("=" * 78)
    print("3. THE TWO ARMS, CONTIG BY CONTIG")
    print("=" * 78)
    idx = sc["indexes"]
    for arm, rows in (("candidate (withheld, deposited after the cutoff)",
                       sc["candidate_rows"]),
                      ("leak (already in the reference)", sc["leak_rows"])):
        print(f"\n--- {arm}: {len(rows)} contigs ---")
        for r in sorted(rows, key=lambda r: -r["length"])[:6]:
            print(f"  {r['contig'].split('_megahit_')[0]}")
            print(f"    {r['length']} nt, {pct(r['attributed_covered_fraction'])} of the "
                  f"contig attributed, tier = {r['tier']}")
            print(f"    rank by length: {r['rank_by_length_all_contigs']} of "
                  f"{sum(v['contigs'] for v in sc['summary'].values())}")
            for label in idx:
                h = r.get(label, {})
                if "subject" in h:
                    print(f"    {label:<14} {h['identity_percent']}% over "
                          f"{h['aligned']} (union cov {h['query_coverage_union']}) "
                          f"E={h['evalue']:.3g} bits={h['bitscore']:.0f} -> "
                          f"{h['subject']}")
                else:
                    print(f"    {label:<14} {h.get('state')}")

    print()
    print("=" * 78)
    print("4. WHERE THE REMOVED PILE'S IDENTITIES SIT - the control that decides whether")
    print("   a distance readout on tier 1 would have surfaced this")
    print("=" * 78)
    rows = []
    tsv = D / "contig_rows_s1.tsv"
    import csv
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["length"] = int(r["length"])
            rows.append(r)
    for label in idx:
        print(f"\n--- {label} ---")
        for tier in ("identified", "assembled_unidentifiable"):
            pile = [r for r in rows
                    if r["tier"].startswith(tier) and r["arm"] == "background"]
            sig = [r for r in pile if r[f"{label}_significant"] == "True"
                   and r[f"{label}_identity"]]
            print(f"  {tier}: {len(pile):,} background contigs, "
                  f"{len(sig):,} with a significant hit "
                  f"({100 * len(sig) / max(1, len(pile)):.1f}%)")
            if not sig:
                continue
            ids = [float(r[f"{label}_identity"]) for r in sig]
            counts = {}
            for x in ids:
                counts[band_of(x)] = counts.get(band_of(x), 0) + 1
            print(f"    median identity {st.median(ids):.1f}%  "
                  f"range {min(ids):.1f}-{max(ids):.1f}")
            for lo, hi in BANDS:
                b = f"{lo}-{hi if hi <= 100 else 100}"
                if b in counts:
                    print(f"      {b:>8}% : {counts[b]:>6,}  "
                          f"({100 * counts[b] / len(sig):.2f}% of significant)")
        # The candidate's own band, named against that distribution.
        for r in rows:
            if r["arm"] == "candidate" and r[f"{label}_identity"]:
                x = float(r[f"{label}_identity"])
                pile = [float(q[f"{label}_identity"]) for q in rows
                        if q["arm"] == "background" and q["tier"] == "identified"
                        and q[f"{label}_significant"] == "True"
                        and q[f"{label}_identity"]]
                if pile:
                    below = sum(1 for y in pile if y <= x)
                    print(f"    CANDIDATE sits at {x:.1f}% ({band_of(x)}% band): "
                          f"{below:,} of {len(pile):,} identified background contigs "
                          f"({100 * below / len(pile):.2f}%) are at or below it")
                # Distance says how far the relative was; bitscore says how marginal the
                # removal was. A contig removed far ABOVE the tier's threshold cannot be
                # spared by raising it, which is the obvious first fix and is why this
                # column is here.
                bs = r.get(f"{label}_bitscore")
                bpile = [float(q[f"{label}_bitscore"]) for q in rows
                         if q["arm"] == "background" and q["tier"] == "identified"
                         and q.get(f"{label}_bitscore")]
                if bs and bpile:
                    bs = float(bs)
                    under = sum(1 for y in bpile if y <= bs)
                    print(f"    CANDIDATE bitscore {bs:.0f} against the tier threshold "
                          f"of 50: {100 * under / len(bpile):.2f}% of identified "
                          f"background contigs score at or below it, so raising the "
                          f"threshold to spare it would also spare that share")
    print()
    print("=" * 78)
    print("4b. COVERAGE x IDENTITY JOINTLY, over the pile the cascade removed, and the")
    print("    SAME cell counts for the leak arm. Identity alone is uninterpretable, so")
    print("    the question is whether adding the axis that has calibration behind it")
    print("    isolates the candidate. The leak arm is the direction check: it is known")
    print("    material from this same library and must land top-right.")
    print("=" * 78)
    COV = [(0.70, 1.01), (0.50, 0.70), (0.25, 0.50), (0.00, 0.25)]
    IDB = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100.01)]
    for label in ("holdout_prot", "holdout_nt"):
        for arm, tier in (("background", "identified"), ("leak", None)):
            sel = [r for r in rows
                   if r["arm"] == arm
                   and (tier is None or r["tier"] == tier)
                   and r[f"{label}_significant"] == "True"
                   and r[f"{label}_identity"] and r[f"{label}_coverage"]]
            if not sel:
                continue
            print(f"\n  --- {label}, arm={arm}"
                  + (f", tier={tier}" if tier else "") + f": n={len(sel):,} ---")
            hdr = "  cov \\ id      " + "".join(
                f"{lo}-{min(hi, 100):.0f}".rjust(9) for lo, hi in IDB)
            print(hdr)
            for clo, chi in COV:
                cells = []
                for ilo, ihi in IDB:
                    n = sum(1 for r in sel
                            if clo <= float(r[f"{label}_coverage"]) < chi
                            and ilo <= float(r[f"{label}_identity"]) < ihi)
                    cells.append(f"{n:,}".rjust(9))
                print(f"  {clo:.2f}-{chi:.2f}  " + "".join(cells))
            # Where the candidate falls, named against the cell it falls in.
            for r in rows:
                if r["arm"] != "candidate" or not r[f"{label}_identity"]:
                    continue
                x, c = float(r[f"{label}_identity"]), float(r[f"{label}_coverage"])
                cb = next((p for p in COV if p[0] <= c < p[1]), None)
                ib = next((p for p in IDB if p[0] <= x < p[1]), None)
                if cb and ib and arm == "background":
                    n = sum(1 for q in sel
                            if cb[0] <= float(q[f"{label}_coverage"]) < cb[1]
                            and ib[0] <= float(q[f"{label}_identity"]) < ib[1])
                    print(f"    CANDIDATE: identity {x:.1f}%, union coverage {c:.4f} -> "
                          f"cell cov {cb[0]:.2f}-{cb[1]:.2f} x id {ib[0]}-"
                          f"{min(ib[1], 100):.0f}, which holds {n:,} contigs "
                          f"({100 * n / len(sel):.2f}% of this pile)")
                    # The cell's own top rows, by bitscore: if the candidate is not the
                    # most striking thing in its own cell, a box drawn round the cell is
                    # not a flag for this genome.
                    peers = sorted(
                        (q for q in sel
                         if cb[0] <= float(q[f"{label}_coverage"]) < cb[1]
                         and ib[0] <= float(q[f"{label}_identity"]) < ib[1]),
                        key=lambda q: -float(q[f"{label}_bitscore"]))[:4]
                    for q in peers:
                        print(f"      peer: {q['length']:>6} nt  "
                              f"{float(q[f'{label}_identity']):.1f}%  "
                              f"cov {float(q[f'{label}_coverage']):.3f}  "
                              f"bits {float(q[f'{label}_bitscore']):.0f}")

    print()
    print("=" * 78)
    print("4c. THE COVERAGE GAP: how much more of a contig is explained by known")
    print("    proteins TOGETHER than by the single best-covering one. The one field in")
    print("    this experiment whose direction is pinned on known material from the same")
    print("    water, so it is printed with that arm beside it and with a length-matched")
    print("    comparison, because the gap can only open on a contig with several genes.")
    print("=" * 78)
    for label in ("holdout_prot", "holdout_nt"):
        if not rows or f"{label}_coverage_gap" not in rows[0]:
            print(f"\n  {label}: not measured - this contig_rows file predates the "
                  f"any-subject union column; re-run holdout_cascade.py")
            continue
        cand = next((r for r in rows if r["arm"] == "candidate"
                     and r.get(f"{label}_coverage_gap")), None)
        if not cand:
            print(f"\n  {label}: the candidate has no significant hit, so no gap")
            continue
        g = float(cand[f"{label}_coverage_gap"])
        print(f"\n  --- {label} ---")
        print(f"  CANDIDATE {cand['length']} nt: per-subject "
              f"{float(cand[f'{label}_coverage']):.4f}, any-subject "
              f"{float(cand[f'{label}_coverage_any']):.4f}, GAP +{100 * g:.1f} points "
              f"over {cand[f'{label}_significant_subjects']} significant subjects")
        lo, hi = cand["length"] * 0.7, cand["length"] * 1.3
        def gaps(pred):
            return [100 * float(r[f"{label}_coverage_gap"]) for r in rows
                    if pred(r) and r.get(f"{label}_coverage_gap") not in (None, "")]
        band = gaps(lambda r: r["arm"] == "background" and r["tier"] == "identified"
                    and lo <= r["length"] <= hi)
        if band:
            print(f"  length-matched identified background ({lo:.0f}-{hi:.0f} nt), "
                  f"n={len(band):,}: median {st.median(band):.1f}, "
                  f"mean {st.mean(band):.1f}, max {max(band):.1f}")
            print(f"    candidate percentile: "
                  f"{100 * sum(1 for x in band if x <= 100 * g) / len(band):.2f}%")
            for th in (10, 20, 30):
                n = sum(1 for x in band if x > th)
                print(f"    gap > {th:>2} points: {n:,} ({100 * n / len(band):.2f}%)")
        leak = gaps(lambda r: r["arm"] == "leak")
        if leak:
            print(f"  LEAK ARM, known material from this library, n={len(leak)}: "
                  f"median {st.median(leak):.1f}, mean {st.mean(leak):.1f}, "
                  f"max {max(leak):.1f}")
        # The confound: a gap needs several genes, so it could be a gene-count readout.
        byn: dict[str, list] = {}
        for r in rows:
            if (r["arm"] != "background" or r["tier"] != "identified"
                    or r.get(f"{label}_coverage_gap") in (None, "")):
                continue
            k = r.get(f"{label}_significant_subjects") or "?"
            byn.setdefault(str(k), []).append(
                100 * float(r[f"{label}_coverage_gap"]))
        print("  confound check, gap by number of significant subjects:")
        for k in sorted(byn, key=lambda x: (not x.isdigit(), x.zfill(3))):
            v = byn[k]
            print(f"    {k:>3} subjects: n={len(v):>6,}  median {st.median(v):>5.1f}  "
                  f"mean {st.mean(v):>5.1f}")

    print()
    print("=" * 78)
    print("5. WHAT THE PIPELINE ACTUALLY REPORTED - the `strains` and `layers` output")
    print("   a reader would have in front of them, with the withheld genome removed")
    print("   in tier 1 and therefore absent from both by construction")
    print("=" * 78)
    for s in a.seeds:
        sp = D / f"run_s{s}" / "strains.json"
        lp = D / f"run_s{s}" / "layers.json"
        if not sp.exists():
            continue
        sj = json.loads(sp.read_text())
        subs = sj.get("subjects", [])
        reads = [r["residue_reads"] for r in subs]
        ids = [r["median_identity_percent"] for r in subs]
        aln = [r["median_aligned_aa"] for r in subs]
        # The subject the withheld genome's own contig matched. If `strains` had any
        # purchase on this genome it would name this subject; it is searched for rather
        # than assumed absent, because "we did not see it" and "we did not look" are the
        # failure this experiment already found a bug for.
        wanted = ("074562", "010771404")
        hit = [r["subject"] for r in subs if any(w in r["subject"] for w in wanted)]
        print(f"\n  seed {s}: {sj.get('subjects_total')} subjects over "
              f"{sj.get('known_pile_searched_reads')} searched residue reads")
        print(f"    residue reads per subject: max {max(reads) if reads else 0}, "
              f"subjects with >=2 {sum(1 for x in reads if x >= 2)}, "
              f"with >=3 {sum(1 for x in reads if x >= 3)}")
        if ids:
            print(f"    median aligned {st.median(aln):.0f} aa, identity "
                  f"{min(ids):.1f}-{max(ids):.1f}%")
        print(f"    the withheld genome's own subject present: "
              f"{hit if hit else 'NO'}")
        if lp.exists():
            lj = json.loads(lp.read_text())
            for arm in ("residue", "known"):
                c = lj["counts"].get(arm, {})
                tot = lj["totals"].get(arm)
                print(f"    layers[{arm}] of {tot}: "
                      + "  ".join(f"{k}={v}" for k, v in sorted(c.items())))

    print()
    print("=" * 78)
    print("6. DRIFT CONTROL: held-out index against today's, same pass")
    print("=" * 78)
    print(f"  {sc['drift_control']}")
    for arm in ("candidate", "leak"):
        for r in (sc[f"{arm}_rows"] or [])[:3]:
            pairs = [(l, r.get(l, {}).get("identity_percent"),
                      r.get(l, {}).get("query_coverage_union"))
                     for l in idx]
            print(f"  [{arm}] {r['length']} nt: "
                  + "  ".join(f"{l}={v}%/cov {c}" for l, v, c in pairs if v))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
