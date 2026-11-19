"""Why hapax merges two swarms of the same virus, and what would stop it.

On the multi-swarm metagenome hapax returned one contig for several swarms but
merged each sister pair: the two HIV *env* hosts (86% identical consensuses),
HBV p03 into p10 (91%), HXB2 into the *gag-pol* swarm (89%). This isolates
that failure: one pair at a time, no background, so the only question is the
grouping rule.

The suspect is `layout::group`'s stage-2 merge rule. Two stage-1 clusters are
merged when the number of overlap edges bridging them is at least `min_bridge`
AND at least `bridge_frac` of the smaller cluster's read count. Both terms
grow with coverage: at 300x two sister swarms share enough conserved stretches
to clear an absolute threshold of 5 edges many times over, and `bridge_frac`
compares bridges to *reads*, not to the within-cluster edge density that the
same coverage also inflates. The prediction is that the merge happens at high
coverage and not at low, for the same pair, which is what this sweeps.

For each pair and coverage it reports whether the two swarms ended up in one
read group (merged), and the bridge ratio actually observed: bridging edges
between the two swarms' reads divided by edges within the sparser swarm. That
ratio is the quantity a coverage-independent rule would threshold.

    python experiments/sister_pairs.py --seed 1 --n-seeds 3 --out results/sister_pairs
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import random
import shutil
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HAPAX, ROOT, ensure_release_build, write_manifest  # noqa: E402
from score_swarm_metagenome import assign_windows  # noqa: E402
from swarm_metagenome import DATA, ENV, read_fasta, safe  # noqa: E402

PAIRS = {
    "hiv_env": ("hiv_env_ZM247F", "hiv_env_ZM246F"),      # 86% identical consensuses
    "hbv": ("hbv_p10", "hbv_p03"),                        # 91%
    "hcv_superinfection": ("hcv_Pt3P", "hcv_Pt3S"),       # one host, two strains
}
# The three ways of trying to keep two hosts' swarms apart, plus the knobs that
# were suspected first. Every one of them fails, each for its own reason; the
# reasons are the result (see results/swarm_metagenome/README.md).
GRID = [
    ("default", []),
    ("bridge-frac-0.6", ["--bridge-frac", "0.6"]),
    ("min-bridge-50", ["--min-bridge", "50"]),
    ("t1groups", ["--stage1", "tier1"]),
    ("split-phase", ["--split-phase", "1"]),
    ("split-cut-0.05", ["--split-cut", "0.05"]),
    ("split-lineage", ["--split-lineage", "1"]),
    ("split-lineage-0.075", ["--split-lineage", "1", "--lineage-split-div", "0.075"]),
]


def simulate_pair(pair: tuple[str, str], cov: float, seed: int, work: Path) -> Path:
    """ART reads for both swarms, log-normal haplotype abundance, single-end pool."""
    rng = random.Random(seed)
    src = work / "src"
    src.mkdir(parents=True, exist_ok=True)
    reads = work / "reads.fq"
    jobs = []
    for sw in pair:
        haps = read_fasta(DATA / "panels" / f"{sw}.fa")
        w = [math.exp(rng.gauss(0, 1.0)) for _ in haps]
        tot = sum(w)
        for (n, s), wi in zip(haps, w):
            fn = src / f"{safe(n)}.fa"
            fn.write_text(f">{safe(n)}\n{s}\n")
            jobs.append((safe(n), cov * wi / tot))
    tmp = work / "art"
    tmp.mkdir(exist_ok=True)
    with open(reads, "w") as out:
        for i, (n, c) in enumerate(jobs):
            subprocess.run([str(ENV / "art_illumina"), "-ss", "HS25", "-na", "-q", "-l", "150",
                            "-f", f"{max(c, 0.01):.4f}", "-rs", str(seed * 1000 + i),
                            "-i", str(src / f"{n}.fa"), "-o", str(tmp / f"{n}_")],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            fq = tmp / f"{n}_.fq"
            if not fq.exists():
                continue
            lines = open(fq).read().splitlines()
            for k in range(0, len(lines), 4):
                out.write(f"@{n}__r{k}\n{lines[k + 1]}\n+\n{lines[k + 3]}\n")
    shutil.rmtree(tmp)
    return reads


def swarm_of(read_name: str) -> str:
    return read_name.split("__")[0]


def analyse(out: Path, pair: tuple[str, str], work: Path) -> dict:
    """Score the contig *sequences* against the two swarms, and read edges.tsv
    for the bridge ratio.

    Sequence scoring, not read membership: a split contig shares the reads from
    conserved stretches with its sibling by design (src/phase.rs), so counting
    reads per contig would call every split contig a chimaera. What matters is
    whether the sequence emitted is one host's or a mixture, which is the same
    window test the metagenome scorer uses.
    """
    allfa = work / "pair_truth.fa"
    entity = {}
    if not allfa.exists():
        with open(allfa, "w") as f:
            for sw in pair:
                for n, s in read_fasta(DATA / "panels" / f"{sw}.fa"):
                    f.write(f">{safe(n)}\n{s}\n")
    for sw in pair:
        for n, _ in read_fasta(DATA / "panels" / f"{sw}.fa"):
            entity[safe(n)] = sw
    contigs = [(n, s) for n, s in read_fasta(out / "contigs.fa") if len(s) >= 500]
    asg = assign_windows(contigs, allfa, entity, out) if contigs else {}
    merged = sum(1 for a in asg.values() if a["chimaera"])
    own = {sw: sum(1 for a in asg.values() if a["entity"] == sw) for sw in pair}
    within = collections.Counter()
    between = 0
    for line in open(out / "edges.tsv").readlines()[1:]:
        f = line.split("\t")
        a, b = swarm_of(f[0]), swarm_of(f[1])
        if a == b:
            within[a] += 1
        else:
            between += 1
    sparser = min(within.values()) if within else 0
    return {"chimaeric_contigs": merged, "contigs_per_swarm": own,
            "edges_within": dict(within), "edges_between": between,
            "bridge_ratio": between / sparser if sparser else None}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=1)
    ap.add_argument("--covs", default="30,100,300")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", default=str(ROOT / "results" / "sister_pairs"))
    ap.add_argument("--configs", nargs="+", default=[c for c, _ in GRID], choices=[c for c, _ in GRID])
    a = ap.parse_args()
    ensure_release_build()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.perf_counter()
    for seed in range(a.seed, a.seed + a.n_seeds):
        for pname, pair in PAIRS.items():
            for cov in [float(x) for x in a.covs.split(",")]:
                work = out / "work" / f"{pname}_cov{cov:g}_seed{seed}"
                work.mkdir(parents=True, exist_ok=True)
                reads = simulate_pair(pair, cov, seed, work)
                for cfg, extra in [g for g in GRID if g[0] in a.configs]:
                    o = work / cfg
                    o.mkdir(exist_ok=True)
                    t = time.perf_counter()
                    subprocess.run([str(HAPAX), "assemble", "-i", str(reads), "-o", str(o), "-t", str(a.threads),
                                    "--seed", "spaced", "--dump-edges", "1", *extra],
                                   check=True, stderr=subprocess.DEVNULL)
                    r = analyse(o, pair, work)
                    r.update({"pair": pname, "cov": cov, "seed": seed, "config": cfg, "wall": time.perf_counter() - t})
                    rows.append(r)
                    print(json.dumps(r), flush=True)
                    (o / "edges.tsv").unlink()
                reads.unlink()
                shutil.rmtree(work / "src")
    with open(out / "runs.tsv", "w") as f:
        f.write("pair\tcov\tseed\tconfig\tchimaeric_contigs\tcontigs_a\tcontigs_b\tedges_between\tbridge_ratio\twall\n")
        for r in rows:
            ca, cb = r["contigs_per_swarm"].values()
            f.write(f"{r['pair']}\t{r['cov']:g}\t{r['seed']}\t{r['config']}\t{r['chimaeric_contigs']}\t{ca}\t{cb}\t{r['edges_between']}\t{r['bridge_ratio'] or float('nan'):.4f}\t{r['wall']:.1f}\n")
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r["pair"], r["cov"], r["config"])].append(r)
    lines = ["| pair | coverage | config | chimaeric contigs | contigs per swarm | bridge ratio |", "|---|---|---|---|---|---|"]
    for (p, c, cfg), v in sorted(agg.items()):
        br = [x["bridge_ratio"] for x in v if x["bridge_ratio"] is not None]
        lines.append(f"| {p} | {c:g} | {cfg} | {st.mean(x['chimaeric_contigs'] for x in v):.1f} | "
                     f"{st.mean(sum(x['contigs_per_swarm'].values()) / 2 for x in v):.2f} | "
                     f"{st.mean(br):.3f} |" if br else "")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    write_manifest(out, a.seed, vars(a), {"total": time.perf_counter() - t0})
    print("\n".join(lines))


if __name__ == "__main__":
    main()
