"""Summarise multi-swarm metagenome runs across seeds into committed tables.

    python experiments/report_swarm_metagenome.py results/swarm_metagenome

Writes summary.md, per_swarm.tsv and runs.tsv beside the work directories.
Every number is mean +- sd over the seeds present, and contiguity numbers are
printed in the same row as the chimaera number from the same run, as
CLAUDE.md requires.
"""

from __future__ import annotations

import json
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def ms(v):
    v = [x for x in v if x is not None]
    if not v:
        return "-"
    if len(v) == 1:
        return f"{v[0]:.2f}"
    return f"{st.mean(v):.2f} ± {st.stdev(v):.2f}"


def main() -> None:
    root = Path(sys.argv[1])
    seeds = sorted(root.glob("work/seed*/scores.json"))
    if not seeds:
        raise SystemExit(f"no scores.json under {root}/work")
    per_tool = defaultdict(lambda: defaultdict(list))
    per_swarm = defaultdict(lambda: defaultdict(list))
    for p in seeds:
        sc = json.loads(p.read_text())
        for tool, r in sc.items():
            if not r.get("n_contigs"):
                per_tool[tool]["status_fail"].append(1)
                continue
            swarms = r["swarms"]
            found = [e for e in swarms.values() if e["contigs"] > 0]
            per_tool[tool]["swarms_found"].append(len(found))
            per_tool[tool]["contigs_per_found"].append(st.mean([e["contigs"] for e in found]) if found else None)
            per_tool[tool]["completeness"].append(st.mean([e["completeness"] for e in swarms.values()]))
            per_tool[tool]["cross_swarm_chim"].append(sum(1 for c in r["chimaera_detail"].values()
                                                          if sum(1 for k in c if not k.startswith("NC_")) >= 2))
            per_tool[tool]["chimaeras_all"].append(r["chimaeras"])
            per_tool[tool]["resolved"].append(sum(e["resolved"] for e in swarms.values()))
            per_tool[tool]["n_contigs"].append(r["n_contigs"])
            per_tool[tool]["wall"].append(r["wall"])
            per_tool[tool]["rss"].append(r.get("max_rss_gb"))
            vr = [e["variants_remap"] for e in swarms.values() if e["variants_remap"]["recall"] is not None]
            per_tool[tool]["var_recall"].append(st.mean([x["recall"] for x in vr]) if vr else None)
            per_tool[tool]["var_prec"].append(st.mean([x["precision"] for x in vr if x["precision"] is not None]) if vr else None)
            bg = r["background"]
            per_tool[tool]["bg_frac"].append(st.mean([v["len_frac"] for v in bg.values()]))
            for sw, e in swarms.items():
                for k in ("contigs", "completeness", "resolved", "chimaeric_contigs"):
                    per_swarm[(tool, sw)][k].append(e[k])
                per_swarm[(tool, sw)]["var_recall"].append(e["variants_remap"]["recall"])
                per_swarm[(tool, sw)]["var_prec"].append(e["variants_remap"]["precision"])
                if "variants_own" in e:
                    per_swarm[(tool, sw)]["var_recall_own"].append(e["variants_own"]["recall"])
    n = len(seeds)
    lines = [f"# Multi-swarm metagenome: {n} seed(s)", "",
             "Eight real within-host swarms (see experiments/fetch_swarm_panels.py) on a",
             "bacterial and phage background. Contiguity and chimaeras from the same runs.", "",
             "| tool | swarms found /8 | contigs per found swarm | completeness | cross-swarm chimaeras | all chimaeras | haplotypes resolved | total contigs | variant recall | variant precision | wall (s) | peak RSS (GB) |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for tool in sorted(per_tool):
        d = per_tool[tool]
        if "status_fail" in d and not d.get("n_contigs"):
            lines.append(f"| {tool} | failed / no contigs | | | | | | | | | | |")
            continue
        lines.append("| " + " | ".join([tool, ms(d["swarms_found"]), ms(d["contigs_per_found"]), ms(d["completeness"]),
                                        ms(d["cross_swarm_chim"]), ms(d["chimaeras_all"]), ms(d["resolved"]),
                                        ms(d["n_contigs"]), ms(d["var_recall"]), ms(d["var_prec"]),
                                        ms(d["wall"]), ms(d["rss"])]) + " |")
    (root / "summary.md").write_text("\n".join(lines) + "\n")
    # work/ is gitignored, so copy each seed's manifest and per-seed numbers out
    # beside the tables: a committed result has to carry its own provenance.
    man = root / "manifests"
    man.mkdir(exist_ok=True)
    with open(root / "runs.tsv", "w") as f:
        f.write("seed\ttool\tstatus\tswarms_found\tcontigs_per_found\tcross_swarm_chimaeras\tcompleteness\tvariant_recall\tn_contigs\twall_s\tmax_rss_gb\n")
        for pth in seeds:
            seed_dir = pth.parent
            seed = seed_dir.name.replace("seed", "")
            src = seed_dir / "manifest.json"
            if src.exists():
                (man / f"seed{seed}.json").write_text(src.read_text())
            for tool, r in json.loads(pth.read_text()).items():
                if not r.get("n_contigs"):
                    f.write(f"{seed}\t{tool}\t{r.get('status', 'no contigs')}\t-\t-\t-\t-\t-\t0\t{r.get('wall', 0):.0f}\t-\n")
                    continue
                sw = r["swarms"]
                found = [e for e in sw.values() if e["contigs"] > 0]
                cpf = st.mean([e["contigs"] for e in found]) if found else 0
                cross = sum(1 for c in r["chimaera_detail"].values() if sum(1 for k in c if not k.startswith("NC_")) >= 2)
                vr = [e["variants_remap"]["recall"] for e in sw.values() if e["variants_remap"]["recall"] is not None]
                f.write(f"{seed}\t{tool}\t{r['status']}\t{len(found)}\t{cpf:.2f}\t{cross}\t"
                        f"{st.mean([e['completeness'] for e in sw.values()]):.3f}\t{st.mean(vr) if vr else 0:.3f}\t"
                        f"{r['n_contigs']}\t{r['wall']:.0f}\t{r.get('max_rss_gb') or 0:.2f}\n")
    with open(root / "per_swarm.tsv", "w") as f:
        f.write("tool\tswarm\tcontigs\tcompleteness\tchimaeric_contigs\tresolved\tvar_recall\tvar_precision\tvar_recall_own\n")
        for (tool, sw), d in sorted(per_swarm.items()):
            f.write("\t".join([tool, sw] + [ms(d[k]) for k in ("contigs", "completeness", "chimaeric_contigs", "resolved", "var_recall", "var_prec", "var_recall_own")]) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
