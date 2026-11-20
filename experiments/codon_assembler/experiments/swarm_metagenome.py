"""Multi-swarm metagenome: seven real within-host viral swarms on a bacterial
and phage background, assembled by hapax and by the incumbents.

Why this experiment: the stress-test community (another thread) showed that a
plain metagenome assembler (MEGAHIT) already finds novel and divergent viruses
in a mixed sample, but returns 0 of 10 co-existing HIV env haplotypes, with
chimaeras between them. The open question is what each tool does with
*several* swarms in one sample, some of them sister populations of each other
(two subtype-C env hosts, two HBV hosts), with real bacterial reads around
them. hapax's design target is one variant-aware contig per swarm, so the
questions are: one contig per swarm or many; any joins across swarms or into
the background; and does the variant structure it reports match the swarm's
real haplotype variation.

Truth: experiments/fetch_swarm_panels.py (real SGA/clonal panels). Haplotype
abundances inside each swarm are log-normal (sigma 1), so a few haplotypes
dominate and many sit near 1x, as in real within-host samples. That skew is
drawn per seed, as is the ART read sampling.

Reads: ART HiSeq 2500, 2x150, 350 nt inserts. hapax is single-end, so it gets
both mates as independent reads; the other tools get pairs.

    python experiments/swarm_metagenome.py --seed 1 --n-seeds 5 --out results/swarm_metagenome
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HAPAX, ROOT, ensure_release_build, git_sha, write_manifest  # noqa: E402

DATA = ROOT / "data" / "swarm_metagenome"
ENV = Path(os.environ.get("HAPAX_TOOLS", "/home/user/mamba/envs/asm/bin"))

# fold coverage per background genome
BG_COV = {"NC_000913.3": 5, "NC_003228.3": 5, "NC_004461.1": 5,
          "NC_001416.1": 30, "NC_000866.4": 30, "NC_024711.1": 30, "NC_001422.1": 30, "K03455.1": 30}

HAPAX_CONFIGS = {
    "hapax-codon": ["--seed", "spaced"],
    "hapax-nucleotide": ["--seed", "spaced", "--mode", "nucleotide"],
    "hapax-codon-t1groups": ["--seed", "spaced", "--stage1", "tier1"],
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        out.append((name, "".join(seq)))
    return out


def safe(name: str) -> str:
    return name.replace("|", "__")


def build_community(seed: int, swarm_cov: float, sigma: float, work: Path) -> dict:
    """Write one FASTA per source with its fold coverage. Returns the truth table."""
    rng = random.Random(seed)
    srcdir = work / "sources"
    srcdir.mkdir(parents=True, exist_ok=True)
    truth = {"swarms": {}, "background": {}}
    for p in sorted((DATA / "panels").glob("*.fa")):
        swarm = p.stem
        haps = read_fasta(p)
        w = [math.exp(rng.gauss(0, sigma)) for _ in haps]
        tot = sum(w)
        entry = {}
        for (n, s), wi in zip(haps, w):
            frac = wi / tot
            cov = swarm_cov * frac
            fn = srcdir / f"{safe(n)}.fa"
            fn.write_text(f">{safe(n)}\n{s}\n")
            entry[safe(n)] = {"frac": frac, "cov": cov, "len": len(s)}
        truth["swarms"][swarm] = entry
    for acc, cov in BG_COV.items():
        (n, s), = read_fasta(DATA / "background" / f"{acc}.fa")
        fn = srcdir / f"{acc}.fa"
        fn.write_text(f">{acc}\n{s}\n")
        truth["background"][acc] = {"cov": cov, "len": len(s)}
    return truth


def simulate_reads(seed: int, truth: dict, work: Path) -> tuple[Path, Path, Path]:
    """ART per source at its own fold coverage, then pooled and shuffled per pair."""
    src = work / "sources"
    tmp = work / "art"
    tmp.mkdir(exist_ok=True)
    jobs = []
    for swarm, haps in truth["swarms"].items():
        for n, h in haps.items():
            jobs.append((n, h["cov"]))
    for acc, b in truth["background"].items():
        jobs.append((acc, b["cov"]))
    for i, (n, cov) in enumerate(jobs):
        subprocess.run([str(ENV / "art_illumina"), "-ss", "HS25", "-p", "-na", "-q", "-l", "150", "-m", "350", "-s", "50",
                        "-f", f"{cov:.4f}", "-rs", str(seed * 1000 + i), "-i", str(src / f"{n}.fa"), "-o", str(tmp / f"{n}_")],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pairs = []
    for n, _ in jobs:
        r1, r2 = tmp / f"{n}_1.fq", tmp / f"{n}_2.fq"
        if not r1.exists():
            continue
        a, b = open(r1).read().splitlines(), open(r2).read().splitlines()
        for k in range(0, len(a), 4):
            pairs.append((a[k:k + 4], b[k:k + 4]))
    random.Random(seed).shuffle(pairs)
    p1, p2, pse = work / "reads_1.fq", work / "reads_2.fq", work / "reads_se.fq"
    with open(p1, "w") as f1, open(p2, "w") as f2, open(pse, "w") as fs:
        for k, (x, y) in enumerate(pairs):
            # ART names reads <source>-<n>/1; keep the source for scoring, make names unique and tool-safe
            src_name = x[0][1:].rsplit("-", 1)[0]
            nm = f"r{k}_{src_name}"
            f1.write(f"@{nm}/1\n{x[1]}\n+\n{x[3]}\n")
            f2.write(f"@{nm}/2\n{y[1]}\n+\n{y[3]}\n")
            fs.write(f"@{nm}_1\n{x[1]}\n+\n{x[3]}\n@{nm}_2\n{y[1]}\n+\n{y[3]}\n")
    shutil.rmtree(tmp)
    return p1, p2, pse


def timed(cmd: list[str], log: Path, timeout: int) -> tuple[float, str]:
    t = time.perf_counter()
    with open(log, "w") as f:
        try:
            # run through a tiny wrapper so peak RSS is per tool, not cumulative over the experiment
            wrap = "import resource,subprocess,sys;r=subprocess.call(sys.argv[1:]);print('Maximum resident set size (kbytes):',resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);sys.exit(r)"
            r = subprocess.run([sys.executable, "-c", wrap, *cmd], stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
            status = "ok" if r.returncode == 0 else f"exit {r.returncode}"
        except subprocess.TimeoutExpired:
            status = "timeout"
    return time.perf_counter() - t, status


def run_tools(work: Path, p1: Path, p2: Path, pse: Path, tools: list[str], threads: int, timeout: int) -> dict:
    res = {}
    asm = work / "asm"
    asm.mkdir(exist_ok=True)
    for name, extra in HAPAX_CONFIGS.items():
        if name not in tools:
            continue
        out = asm / name
        out.mkdir(exist_ok=True)
        wall, st = timed([str(HAPAX), "assemble", "-i", str(pse), "-o", str(out), "-t", str(threads), *extra], out / "log.txt", timeout)
        res[name] = {"wall": wall, "status": st, "contigs": str(out / "contigs.fa")}
        if (out / "stats.json").exists():
            res[name]["stats"] = json.loads((out / "stats.json").read_text())
    if "megahit" in tools:
        out = asm / "megahit"
        shutil.rmtree(out, ignore_errors=True)
        wall, st = timed([str(ENV / "megahit"), "-1", str(p1), "-2", str(p2), "-o", str(out), "-t", str(threads)], asm / "megahit.log", timeout)
        res["megahit"] = {"wall": wall, "status": st, "contigs": str(out / "final.contigs.fa")}
    if "metaspades" in tools:
        out = asm / "metaspades"
        # --only-assembler: read correction on 540k reads took longer than the whole
        # rest of the experiment, and the reads are already simulated with known error
        wall, st = timed([str(ENV / "metaspades.py"), "--only-assembler", "-1", str(p1), "-2", str(p2), "-o", str(out), "-t", str(threads), "-m", "12"], asm / "metaspades.log", timeout)
        res["metaspades"] = {"wall": wall, "status": st, "contigs": str(out / "contigs.fasta")}
    if "haploflow" in tools:
        out = asm / "haploflow"
        both = work / "reads_cat.fq"
        with open(both, "w") as f:
            for p in (p1, p2):
                f.write(open(p).read())
        wall, st = timed([str(ENV / "haploflow"), "--read-file", str(both), "--out", str(out), "--log", str(out) + ".log"], asm / "haploflow.time", timeout)
        both.unlink()
        res["haploflow"] = {"wall": wall, "status": st, "contigs": str(out / "contigs.fa")}
    for name, r in res.items():
        log = asm / f"{name}.log" if name in ("megahit", "metaspades") else (asm / name / "log.txt" if name.startswith("hapax") else asm / "haploflow.time")
        try:
            for line in open(log):
                if "Maximum resident set size" in line:
                    r["max_rss_gb"] = int(line.split()[-1]) / 1e6
        except FileNotFoundError:
            pass
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=1)
    ap.add_argument("--swarm-cov", type=float, default=300.0, help="total fold coverage per swarm, split across its haplotypes")
    ap.add_argument("--sigma", type=float, default=1.0, help="log-normal sigma of within-swarm haplotype abundance")
    ap.add_argument("--tools", default="hapax-codon,hapax-nucleotide,hapax-codon-t1groups,megahit,metaspades,haploflow")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--out", default=str(ROOT / "results" / "swarm_metagenome"))
    ap.add_argument("--keep-reads", action="store_true")
    a = ap.parse_args()
    ensure_release_build()
    tools = a.tools.split(",")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for seed in range(a.seed, a.seed + a.n_seeds):
        work = out / "work" / f"seed{seed}"
        work.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()
        truth = build_community(seed, a.swarm_cov, a.sigma, work)
        (work / "truth.json").write_text(json.dumps(truth, indent=1))
        p1, p2, pse = simulate_reads(seed, truth, work)
        t_sim = time.perf_counter() - t0
        res = run_tools(work, p1, p2, pse, tools, a.threads, a.timeout)
        (work / "runs.json").write_text(json.dumps(res, indent=1))
        subprocess.run([sys.executable, str(ROOT / "experiments" / "score_swarm_metagenome.py"), str(work)], check=True)
        write_manifest(work, seed, vars(a), {"simulate": t_sim, **{k: v["wall"] for k, v in res.items()}},
                       {"status": {k: v["status"] for k, v in res.items()}})
        if not a.keep_reads:
            for p in (p1, p2, pse):
                p.unlink(missing_ok=True)
    print("done; summarise with experiments/report_swarm_metagenome.py", out)


if __name__ == "__main__":
    main()
