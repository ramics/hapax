"""The diversity ladder (METHOD.md 4): does synonymous collapse buy anything?

For each seed and each within-population diversity (0.05-0.20 subs/site, the
thesis ladder), simulate a world of four entities -- two sister populations
plus two unrelated families -- and assemble it with:

* hapax codon mode (the method),
* hapax nucleotide mode (identical pipeline, codon states switched off),

each with contiguous and codon-spaced seeds, across a sweep of the tier-2
acceptance threshold tau. The codon-vs-nucleotide comparison at matched
pipeline is the cleanest test of the thesis: anything else that differs
between hapax and MEGAHIT confounds it. The tau sweep is there because a
single threshold could favour either mode; what matters is the trade-off
curve between contigs-per-entity and chimaera rate, so every contiguity
number is reported with the chimaera number from the same run.

External assemblers (MEGAHIT, metaSPAdes, PLASS, PenguiN) run when they are on
PATH and are scored with the contig-only k-mer evaluator; they are not
installed in the development sandbox, so `--external` is opt-in.

Usage:
    python experiments/diversity_ladder.py --seed 1 --n-seeds 5 --out results/ladder
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import time
from pathlib import Path

from common import ROOT, ensure_release_build, run_hapax, write_manifest

from python import evaluate, swarm

CONFIGS = {
    "codon": ["--mode", "codon", "--seed", "contiguous"],
    "nucleotide": ["--mode", "nucleotide", "--seed", "contiguous"],
    "codon-spaced": ["--mode", "codon", "--seed", "spaced", "-k", "16"],
    "nucleotide-spaced": ["--mode", "nucleotide", "--seed", "spaced", "-k", "16"],
    # dispatch ablations (experiments/README: "dispatch")
    "codon-spaced-v0": ["--mode", "codon", "--seed", "spaced", "-k", "16", "--bucket-cap", "0", "--lazy-t2", "0"],
    "codon-spaced-cap": ["--mode", "codon", "--seed", "spaced", "-k", "16", "--lazy-t2", "0"],
    "codon-spaced-t1groups": ["--mode", "codon", "--seed", "spaced", "-k", "16", "--stage1", "tier1"],
    "nucleotide-spaced-t1groups": ["--mode", "nucleotide", "--seed", "spaced", "-k", "16", "--stage1", "tier1"],
    "nucleotide-spaced-v0": ["--mode", "nucleotide", "--seed", "spaced", "-k", "16", "--bucket-cap", "0", "--lazy-t2", "0"],
}
DEFAULT_CONFIGS = ["codon", "nucleotide", "codon-spaced", "nucleotide-spaced"]


def run_external(tool: str, reads: Path, out: Path, threads: int) -> Path | None:
    """Run one external assembler; returns its contig FASTA or None."""
    if out.exists():
        shutil.rmtree(out)
    if tool == "megahit" and shutil.which("megahit"):
        subprocess.run(["megahit", "-r", str(reads), "-o", str(out), "-t", str(threads)], check=True, capture_output=True)
        return out / "final.contigs.fa"
    if tool == "metaspades" and shutil.which("spades.py"):
        # metaSPAdes needs paired-end input; single-end reads run as --rnaviral
        subprocess.run(["spades.py", "--rnaviral", "-s", str(reads), "-o", str(out), "-t", str(threads)], check=True, capture_output=True)
        return out / "contigs.fasta"
    if tool == "penguin" and shutil.which("penguin"):
        out.mkdir(parents=True)
        subprocess.run(["penguin", "nuclassemble", str(reads), str(out / "contigs.fa"), str(out / "tmp"), "--threads", str(threads)], check=True, capture_output=True)
        return out / "contigs.fa"
    # PLASS emits protein contigs; it is scored separately (protein recovery),
    # not with the nucleotide evaluator, so it is not run here.
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True, help="first seed")
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--diversities", type=float, nargs="+", default=[0.05, 0.10, 0.15, 0.20])
    ap.add_argument("--taus", type=float, nargs="+", default=[0.7, 0.9, 1.0])
    ap.add_argument("--configs", nargs="+", default=DEFAULT_CONFIGS, choices=list(CONFIGS))
    ap.add_argument("--sister-div", type=float, default=0.30)
    ap.add_argument("--coverage", type=float, default=50.0)
    ap.add_argument("--external", nargs="*", default=[], help="megahit metaspades penguin")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    ensure_release_build()
    out = a.out
    work = out / "work"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    t_start = time.perf_counter()
    t_sim = 0.0
    for s in range(a.seed, a.seed + a.n_seeds):
        for div in a.diversities:
            world_seed = s * 1000 + int(round(div * 100))
            wdir = work / f"s{s}_d{div:.2f}"
            t = time.perf_counter()
            w = swarm.build_world(world_seed, div, sister_div=a.sister_div, coverage=a.coverage)
            swarm.write_world(w, wdir, {"seed": world_seed, "diversity": div})
            t_sim += time.perf_counter() - t
            for cfg in a.configs:
                for tau in a.taus:
                    adir = wdir / f"{cfg}_tau{tau}"
                    st = run_hapax(wdir / "reads.fa", adir, CONFIGS[cfg] + ["--tau", str(tau), "-t", str(a.threads)])
                    ev = evaluate.membership(wdir, adir)
                    rows.append(
                        {
                            "seed": s,
                            "world_seed": world_seed,
                            "diversity": div,
                            "assembler": f"hapax-{cfg}",
                            "tau": tau,
                            "contigs": ev["contigs"],
                            "mean_contigs_per_found_entity": ev["mean_contigs_per_found_entity"],
                            "entities_missed": ev["entities_missed"],
                            "chimeric_contigs": ev["chimeric_contigs"],
                            "chimaera_rate": ev["chimaera_rate"],
                            "chimeric_base_fraction": ev["chimeric_base_fraction"],
                            "mean_recovered": ev["mean_recovered"],
                            "contigs_per_entity": ",".join(map(str, ev["contigs_per_entity"])),
                            "hmm_fraction": st["hmm_fraction"],
                            "wall_s": st["wall"],
                        }
                    )
            for tool in a.external:
                fa = run_external(tool, wdir / "reads.fa", wdir / tool, a.threads)
                if fa is None or not fa.exists():
                    continue
                ev = evaluate.kmer(wdir, fa)
                rows.append({"seed": s, "world_seed": world_seed, "diversity": div, "assembler": tool, "tau": "", **{k: ev[k] for k in ("contigs", "mean_contigs_per_found_entity", "entities_missed", "chimeric_contigs", "chimaera_rate", "chimeric_base_fraction", "mean_recovered")}, "contigs_per_entity": ",".join(map(str, ev["contigs_per_entity"])), "hmm_fraction": "", "wall_s": ""})
            print(f"seed {s} div {div} done ({time.perf_counter() - t_start:.0f}s)", flush=True)

    with open(out / "runs.tsv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)

    # summary: mean and sd across seeds
    keys = sorted({(r["diversity"], r["assembler"], r["tau"]) for r in rows}, key=lambda x: (x[0], x[1], str(x[2])))
    metrics = ["mean_contigs_per_found_entity", "entities_missed", "chimaera_rate", "chimeric_base_fraction", "mean_recovered", "hmm_fraction", "wall_s"]
    summ = []
    for k in keys:
        rs = [r for r in rows if (r["diversity"], r["assembler"], r["tau"]) == k]
        row = {"diversity": k[0], "assembler": k[1], "tau": k[2], "n_seeds": len(rs)}
        for m in metrics:
            vals = [float(r[m]) for r in rs if r[m] != ""]
            if vals:
                row[m + "_mean"] = round(statistics.mean(vals), 4)
                row[m + "_sd"] = round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0
        summ.append(row)
    with open(out / "summary.tsv", "w", newline="") as f:
        fields = sorted({k for r in summ for k in r}, key=lambda x: (x not in ("diversity", "assembler", "tau", "n_seeds"), x))
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        wr.writeheader()
        wr.writerows(summ)

    write_manifest(
        out,
        {"first": a.seed, "n": a.n_seeds},
        {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()},
        {"total": time.perf_counter() - t_start, "simulation": t_sim},
        {"configs": CONFIGS},
    )
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
