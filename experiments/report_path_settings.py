#!/usr/bin/env python3
"""Print the tables of results/path_settings/ from its committed outputs.

Nothing here computes a measurement; every number is read from a file written by
`experiments/path_settings.py` or `experiments/subtract_known.py`, because a number typed
into a report is a number without provenance and this repo has had to withdraw claims whose
prose drifted from its outputs.

THE READING ORDER IS THE POINT. The ladder is read against two rows and not on its own:
`baseline_repeat` is the same setting run again, so any arm that moves a column by less
than that row does bounds nothing; and the absent arm of the panel (53 genomes with zero
reads mapped in this library) is the false-join control. A setting that recovers more
genome while its purity falls has absorbed sequence into longer contigs, which is the
failure this experiment exists to catch, so purity and recovery are printed side by side
and never one without the other.
"""
from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def q(xs, p):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    k = max(0, min(len(xs) - 1, int(round(p * (len(xs) - 1)))))
    return xs[k]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", type=Path, default=ROOT / "results/path_settings")
    ap.add_argument("--seeds", type=int, nargs="*", default=[1, 2, 3, 4, 5])
    args = ap.parse_args()
    d = args.dir
    arms = [p.name[len("assembly_"):-len(".json")] for p in sorted(d.glob("assembly_*.json"))]
    man = json.loads(next(d.glob("manifest_s*.json")).read_text())
    ladder = man.get("ladder", {})

    print("## 1. What each setting is, and which way it was expected to move\n")
    print("| arm | flags added to the baseline | stated direction |")
    print("|---|---|---|")
    for a in arms:
        spec = ladder.get(a, {})
        flags = " ".join(spec.get("flags", [])) or "(none)"
        print(f"| `{a}` | `{flags}` | {spec.get('direction', '')} |")

    print("\n## 2. The assemblies\n")
    print("| arm | contigs | bases | N50 | max | >=1 kb | >=3 kb | minutes | threads |")
    print("|---|---|---|---|---|---|---|---|---|")
    for a in arms:
        s = json.loads((d / f"assembly_{a}.json").read_text())
        mins = f"{s['seconds'] / 60:.1f}" if s.get("seconds") else "-"
        print(f"| `{a}` | {s['contigs']:,} | {s['bases']:,} | {s['n50']:,} | {s['max']:,} "
              f"| {s['ge1000']:,} | {s['ge3000']:,} | {mins} | {s['threads']} |")

    recs = {}
    for a in arms:
        p = d / f"recovery_{a}.json"
        if p.exists():
            recs[a] = json.loads(p.read_text())
    if recs:
        print("\n## 3. Ground truth: 25 genomes the depositors assembled from this same water\n")
        print("| arm | >=90% recovered | >=50% | mean breadth | median purity | "
              "p10 purity | contig joins | false joins (control) |")
        print("|---|---|---|---|---|---|---|---|")
        for a, s in recs.items():
            pres = [r for r in s["genomes"].values() if r["present"]]
            pur = [r["best_contig_purity"] for r in pres]
            print(f"| `{a}` | {s['recovered_90pct']}/{s['panel_present']} "
                  f"| {s['recovered_50pct']}/{s['panel_present']} "
                  f"| {s['mean_covered_frac']:.3f} "
                  f"| {q(pur, 0.5):.3f} | {q(pur, 0.1):.3f} "
                  f"| {s['inter_group_mosaics']} | {len(s['false_joins_absent_arm'])}/"
                  f"{s['panel_absent']} |")

        base = recs.get("baseline")
        if base:
            print("\n### 3a. The trade, which the aggregate hides\n")
            print("Genomes whose breadth moves by at least 10 points against the "
                  "`baseline` arm. The `baseline_repeat` row is the null: it is the same "
                  "setting, so whatever it shows is what a difference has to beat.\n")
            print("| arm | gained >=10 pts | lost >=10 pts | net | "
                  "recovered but embedded (breadth >=0.9, purity <0.9) |")
            print("|---|---|---|---|---|")
            for a, s in recs.items():
                gain = loss = 0
                for acc, r in s["genomes"].items():
                    if not r["present"]:
                        continue
                    delta = r["covered_frac"] - base["genomes"][acc]["covered_frac"]
                    if delta >= 0.10:
                        gain += 1
                    elif delta <= -0.10:
                        loss += 1
                emb = sum(1 for r in s["genomes"].values()
                          if r["present"] and r["covered_frac"] >= 0.9
                          and r["best_contig_purity"] is not None
                          and r["best_contig_purity"] < 0.9)
                print(f"| `{a}` | {gain} | {loss} | {gain - loss:+d} | {emb} |")

        print("\n### 3b. Per genome, baseline against every other arm\n")
        base = recs.get("baseline")
        if base:
            order = sorted(((acc, r) for acc, r in base["genomes"].items() if r["present"]),
                           key=lambda kv: -kv[1]["mean_depth"])
            print("| genome | depth | divergence % | " +
                  " | ".join(f"{a}: breadth/purity" for a in recs) + " |")
            print("|---" * (3 + len(recs)) + "|")
            for acc, r in order:
                cells = []
                for a, s in recs.items():
                    g = s["genomes"][acc]
                    pur = g["best_contig_purity"]
                    cells.append(f"{g['covered_frac']:.2f}/" +
                                 (f"{pur:.2f}" if pur is not None else "-"))
                dv = r["divergence_pct"]
                print(f"| {acc} | {r['mean_depth']:.1f}x | "
                      f"{dv if dv is not None else '-'} | " + " | ".join(cells) + " |")

    ids = {}
    for a in arms:
        p = d / f"identification_{a}.json"
        if p.exists():
            ids[a] = json.loads(p.read_text())
    if ids:
        print("\n## 3c. What each arm leaves unidentifiable, which is tier 2\n")
        print("The absorption question the truth panel cannot answer. A setting that fuses "
              "novel sequence into abundant known contigs shows up here as identified bases "
              "rising while tier-2 mass falls, with no change in the panel. Identification "
              "is against indexes rolled back to 2022-05-10, the day this water was "
              "collected, at the committed operating point (Kraken2 confidence 0.1, "
              "protein bitscore 50).\n")
        print("| arm | identified contigs % | identified bases % | tier-2 contigs | "
              "tier-2 bases | tier-2 bases >=500 nt and >=2x |")
        print("|---|---|---|---|---|---|")
        for a, s in ids.items():
            t2 = s["tier2_population"]
            print(f"| `{a}` | {s['identified_frac'] * 100:.2f} "
                  f"| {s['identified_bases_frac'] * 100:.2f} | {t2['contigs']:,} "
                  f"| {t2['bases']:,} | {t2.get('frac_over_500nt', 0) * 100:.1f}% / "
                  f"{t2.get('frac_over_2x', 0) * 100:.1f}% |")

    cas = {}
    for a in arms:
        rows = []
        for sd in args.seeds:
            p = d / f"cascade_{a}_s{sd}.json"
            if p.exists():
                rows.append(json.loads(p.read_text()))
        if rows:
            cas[a] = rows
    if cas:
        print("\n## 4. Residue, as an output of each arm\n")
        print("NOTE: these fractions are not comparable to the 77.10% of "
              "`results/subtract_known/`. The two 8 GB Kraken2 tiers are not staged here "
              "(fixed across arms, marginal 0.27 and 0.00 points there), and this is a "
              "different library. They are comparable across the arms below and nowhere else.\n")
        print("| arm | removed % | residue % (mean) | residue min-max | tier 2 % "
              "(assembled, unidentifiable) | tier 3 % (maps nothing) | seeds |")
        print("|---|---|---|---|---|---|---|")
        for a, rows in cas.items():
            res = [(1 - r["union_removed"]) * 100 for r in rows]
            rem = [r["union_removed"] * 100 for r in rows]
            t2 = [sum(r.get("residue_composition", {}).values()) * 100 for r in rows]
            t3 = [r["tier3_pairs"] / r["pairs"] * 100 for r in rows
                  if r.get("tier3_pairs") is not None]
            print(f"| `{a}` | {st.mean(rem):.2f} | {st.mean(res):.2f} "
                  f"| {min(res):.2f}-{max(res):.2f} | {st.mean(t2):.2f} "
                  f"| {st.mean(t3):.2f} | {len(res)} |" if t3 else
                  f"| `{a}` | {st.mean(rem):.2f} | {st.mean(res):.2f} "
                  f"| {min(res):.2f}-{max(res):.2f} | {st.mean(t2):.2f} | - | {len(res)} |")

    if cas and "baseline" in cas:
        print("\n### 4b. Residue against baseline, PAIRED BY SEED\n")
        print("Five seeds are five disjoint blocks of one read pool, so the same seed is "
              "the same reads and the paired difference is the comparison with the "
              "subsampling noise taken out. The `baseline_repeat` row is the null.\n")
        bl = {r["seed"]: r for r in cas["baseline"]}
        print("| arm | residue change (pts, mean) | per-seed range | tier 2 change | "
              "tier 3 change | seeds paired |")
        print("|---|---|---|---|---|---|")
        for a, rows in cas.items():
            d_res, d_t2, d_t3 = [], [], []
            for r in rows:
                b = bl.get(r["seed"])
                if not b:
                    continue
                d_res.append(((1 - r["union_removed"]) - (1 - b["union_removed"])) * 100)
                d_t2.append((sum(r.get("residue_composition", {}).values())
                             - sum(b.get("residue_composition", {}).values())) * 100)
                if r.get("tier3_pairs") is not None and b.get("tier3_pairs") is not None:
                    d_t3.append((r["tier3_pairs"] / r["pairs"]
                                 - b["tier3_pairs"] / b["pairs"]) * 100)
            if d_res:
                print(f"| `{a}` | {st.mean(d_res):+.2f} | {min(d_res):+.2f} to "
                      f"{max(d_res):+.2f} | {st.mean(d_t2):+.2f} "
                      f"| {st.mean(d_t3):+.2f} | {len(d_res)} |"
                      if d_t3 else
                      f"| `{a}` | {st.mean(d_res):+.2f} | {min(d_res):+.2f} to "
                      f"{max(d_res):+.2f} | {st.mean(d_t2):+.2f} | - | {len(d_res)} |")

    absf = d / "absorption_s1.json"
    if absf.exists():
        ab = json.loads(absf.read_text())
        print("\n## 4c. Absorption at the contig level\n")
        print("How much of the BASELINE's unidentifiable pile now sits inside a contig the "
              "cascade calls identified, per arm. `baseline` is 0 by construction (the two "
              "piles are disjoint) and `baseline_repeat` is the null: the same setting at a "
              "different thread count, so it is what this reads when nothing changed.\n")
        print("| arm | contigs >=90% inside an identified contig | of bases | "
              "coverage-weighted | of which the single largest contig |")
        print("|---|---|---|---|---|")
        for a, r in ab.items():
            t1 = r.get("top1_share_of_absorbed_mass")
            print(f"| `{a}` | {r['contigs_ninety_covered']} of "
                  f"{r['baseline_unidentified_contigs']:,} "
                  f"| {r['bases_absorbed_frac'] * 100:.2f}% "
                  f"| {r['coverage_weighted_mass_frac'] * 100:.2f}% "
                  f"| {t1 * 100:.0f}% |" if t1 is not None else
                  f"| `{a}` | {r['contigs_ninety_covered']} of "
                  f"{r['baseline_unidentified_contigs']:,} "
                  f"| {r['bases_absorbed_frac'] * 100:.2f}% "
                  f"| {r['coverage_weighted_mass_frac'] * 100:.2f}% | - |")

    print("\n## 5. Provenance\n")
    print(f"- script `{man['script']}` at `{man['script_sha256'][:12]}`, git `{man['git_sha']}`")
    print(f"- {man.get('megahit_version', '?')}, minimap2 {man.get('minimap2_version', '?')}")
    print(f"- reads sha256 `{str(man.get('reads_sha256'))[:12]}`")


if __name__ == "__main__":
    main()
