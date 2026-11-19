"""Where does assembly stop and read-level evidence have to start?

Both open directions for this project rest on the same unmeasured number. The
contig-classifier argument says the assembler sets the sensitivity ceiling,
because viral classifiers need long contigs. The read-level argument says that
at surveillance abundance there is no contig at all. Lander-Waterman says they
are an order of magnitude apart in abundance and nobody has measured the
crossover on real reads.

So: spike a held-out real viral genome into a real wastewater background at a
series of absolute read counts, assemble with several assemblers and several
minimum-overlap settings, and report the contig length distribution of the
spike-derived contigs at each rung. The deliverable is a curve, not a
threshold, because the other thread's check showed the 1500 bp figure is a
benchmark harmonisation filter rather than a tool requirement (geNomad states
no minimum; VIBRANT's four-ORF rule is the real floor).

Why absolute read counts rather than relative abundance: coverage over the
target genome is what decides whether an overlap exists, and the Berlin run is
only 8.7M pairs, so a purely relative design would put one-in-a-million at nine
reads. Relative abundance is reported alongside, against the background size
actually used.

Why a real background: SecureBio's 4.9e-7 spike-in was about 10x by Poisson and
still stopped at a zero-coverage gap, so over-dispersion of real coverage is
what decides outcomes, and a simulated background cannot show it.

Why k is the minimum-overlap arm: in a de Bruijn assembler the k-mer size *is*
the minimum overlap two reads must share to be joined. Relaxing it is the lever
the sister thread measured as moving contig mass over 1500 bp from 20% to 52%
at 3x coverage, so it belongs in the sweep rather than being left at defaults.
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import git_sha  # noqa: E402

# Spike sizes in read PAIRS. The low end is the SecureBio regime (hundreds of
# pairs); the high end is where contigs are expected to cross every classifier
# floor, so the curve should bracket the crossover rather than sit on one side.
SPIKES = [10, 30, 100, 300, 1000, 3000]

# Assembler arms. The k-list arms are the minimum-overlap lever on MEGAHIT:
# default sweeps k=21..141, 'lowk' refuses nothing above 21, 'highk' demands a
# 71 nt overlap. metaSPAdes needs --only-assembler (read correction alone ran
# over an hour in this sandbox) and asks for 12 GB, so it never runs beside
# another job.
ARMS = {
    "megahit": ["--min-contig-len", "200"],
    "megahit-lowk": ["--min-contig-len", "200", "--k-list", "21,29,39"],
    "megahit-highk": ["--min-contig-len", "200", "--k-list", "71,99,127"],
}


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def fastq_pair_count(path: Path) -> int:
    out = sh(["seqkit", "stats", "-T", str(path)]).stdout.strip().splitlines()[1]
    return int(out.split("\t")[3])


def make_background(r1: Path, r2: Path, pairs: int, out: Path, seed: int) -> tuple[Path, Path]:
    """A fixed background subsample, identical across every rung and arm.

    Held fixed on purpose: if the background moved between rungs, a difference
    in contig length could come from the background rather than from the spike.
    """
    out.mkdir(parents=True, exist_ok=True)
    b1, b2 = out / "bg_1.fq.gz", out / "bg_2.fq.gz"
    if b1.exists() and b2.exists():
        return b1, b2
    for src, dst in ((r1, b1), (r2, b2)):
        # head, not random sampling: mates must stay in register, and seqkit
        # sample with the same seed on each mate is not guaranteed to agree
        # across files of unequal record order.
        sh(["seqkit", "head", "-n", str(pairs), "-o", str(dst), str(src)])
    return b1, b2


def simulate_spike_pool(genome: Path, out: Path, read_len: int, coverage: int, seed: int) -> tuple[Path, Path]:
    """One deep pool of spike reads per seed; rungs are nested subsamples of it.

    Nesting matters: a rung of 300 pairs then contains the 100-pair rung's
    reads, so differences along the curve are added coverage and not a
    different draw.
    """
    out.mkdir(parents=True, exist_ok=True)
    pref = out / "spike"
    p1, p2 = Path(str(pref) + "1.fq"), Path(str(pref) + "2.fq")
    if not (p1.exists() and p2.exists()):
        sh([
            "art_illumina", "-ss", "MSv3", "-i", str(genome), "-p",
            "-l", str(read_len), "-f", str(coverage), "-m", "250", "-s", "40",
            "-rs", str(seed), "-na", "-o", str(pref),
        ])
    return p1, p2


def take_pairs(p1: Path, p2: Path, n: int, out: Path, tag: str) -> tuple[Path, Path]:
    a, b = out / f"{tag}_1.fq", out / f"{tag}_2.fq"
    sh(["seqkit", "head", "-n", str(n), "-o", str(a), str(p1)])
    sh(["seqkit", "head", "-n", str(n), "-o", str(b), str(p2)])
    return a, b


def assemble(arm: str, r1: Path, r2: Path, out: Path, threads: int) -> tuple[Path | None, float]:
    t = time.perf_counter()
    if out.exists():
        shutil.rmtree(out)
    if arm.startswith("megahit"):
        sh(["megahit", "-1", str(r1), "-2", str(r2), "-o", str(out),
            "-t", str(threads), *ARMS[arm]])
        contigs = out / "final.contigs.fa"
    elif arm == "metaspades":
        sh(["spades.py", "--meta", "--only-assembler", "-1", str(r1), "-2", str(r2),
            "-o", str(out), "-t", str(threads), "-m", "12"])
        contigs = out / "contigs.fasta"
    else:
        raise ValueError(arm)
    wall = time.perf_counter() - t
    return (contigs if contigs.exists() else None), wall


def score(contigs: Path, genome: Path, min_ident: float = 0.90) -> dict:
    """Spike-derived contigs, their length distribution, breadth and chimaeras.

    A contig counts as spike-derived when an alignment to the spike genome at
    >=min_ident covers at least 80% of the contig. A contig with >=500 nt
    aligned to the spike and >=500 nt not aligned is a candidate chimaera: real
    background sequence welded to the target. Contiguity never ships without
    that number (CLAUDE.md).
    """
    paf = sh(["minimap2", "-x", "asm10", "-c", "--secondary=no", str(genome), str(contigs)]).stdout
    glen = int(sh(["seqkit", "fx2tab", "-nl", str(genome)]).stdout.split("\t")[1])
    by_contig: dict[str, dict] = {}
    covered: set[int] = set()
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        name, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        ts, te = int(f[7]), int(f[8])
        matches, blocklen = int(f[9]), int(f[10])
        if blocklen == 0 or matches / blocklen < min_ident:
            continue
        rec = by_contig.setdefault(name, {"qlen": qlen, "aligned": 0})
        rec["aligned"] += qe - qs
        covered.update(range(ts, te))
    spike, chimaera = [], 0
    for name, rec in by_contig.items():
        if rec["aligned"] >= 0.8 * rec["qlen"]:
            spike.append(rec["qlen"])
        elif rec["aligned"] >= 500 and rec["qlen"] - rec["aligned"] >= 500:
            chimaera += 1
    spike.sort(reverse=True)
    total = sum(spike)
    return {
        "spike_contigs": len(spike),
        "spike_mass": total,
        "longest": spike[0] if spike else 0,
        "mass_over_500": sum(x for x in spike if x >= 500),
        "mass_over_1000": sum(x for x in spike if x >= 1000),
        "mass_over_1500": sum(x for x in spike if x >= 1500),
        "n_over_500": sum(1 for x in spike if x >= 500),
        "n_over_1500": sum(1 for x in spike if x >= 1500),
        "breadth": len(covered) / glen,
        "genome_len": glen,
        "chimaera_candidates": chimaera,
        "lengths": spike[:50],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=1)
    ap.add_argument("--bg1", type=Path, required=True)
    ap.add_argument("--bg2", type=Path, required=True)
    ap.add_argument("--bg-pairs", type=int, default=1_000_000)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--read-len", type=int, default=107, help="matches SRR23973193")
    ap.add_argument("--arms", default="megahit,megahit-lowk,megahit-highk")
    ap.add_argument("--spikes", default=",".join(str(s) for s in SPIKES))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/home/user/sweepdata/work"))
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    spikes = [int(s) for s in args.spikes.split(",")]
    arms = args.arms.split(",")

    bg1, bg2 = make_background(args.bg1, args.bg2, args.bg_pairs, args.work / "bg", args.seed)
    rows, timings = [], {}
    t0 = time.perf_counter()

    for seed in range(args.seed, args.seed + args.n_seeds):
        random.seed(seed)
        pool1, pool2 = simulate_spike_pool(
            args.genome, args.work / f"pool_s{seed}", args.read_len, coverage=200, seed=seed
        )
        for n in spikes:
            s1, s2 = take_pairs(pool1, pool2, n, args.work, f"spk_s{seed}_n{n}")
            m1 = args.work / f"mix_s{seed}_n{n}_1.fq.gz"
            m2 = args.work / f"mix_s{seed}_n{n}_2.fq.gz"
            for mix, bg, spk in ((m1, bg1, s1), (m2, bg2, s2)):
                with open(mix, "wb") as fh:
                    subprocess.run(["bash", "-c", f"cat {bg} <(gzip -c {spk})"],
                                   stdout=fh, check=True)
            for arm in arms:
                tag = f"s{seed}_n{n}_{arm}"
                try:
                    contigs, wall = assemble(arm, m1, m2, args.work / f"asm_{tag}", args.threads)
                except subprocess.CalledProcessError as e:
                    print(f"[{tag}] assembler failed: {e.stderr[-400:]}", flush=True)
                    continue
                if contigs is None:
                    print(f"[{tag}] no contigs", flush=True)
                    continue
                sc = score(contigs, args.genome)
                cov = n * 2 * args.read_len / sc["genome_len"]
                row = {"seed": seed, "spike_pairs": n, "arm": arm,
                       "coverage": round(cov, 2),
                       "rel_abundance": n / args.bg_pairs,
                       "wall": round(wall, 1), **sc}
                rows.append(row)
                timings[tag] = round(wall, 1)
                print(json.dumps({k: v for k, v in row.items() if k != "lengths"}), flush=True)
                (args.out / "rows.jsonl").open("a").write(json.dumps(row) + "\n")
                shutil.rmtree(args.work / f"asm_{tag}", ignore_errors=True)
            for p in (m1, m2, s1, s2):
                p.unlink(missing_ok=True)

    (args.out / "manifest.json").write_text(json.dumps({
        "git_sha": git_sha(),
        "params": vars(args) | {"bg1": str(args.bg1), "bg2": str(args.bg2),
                                "genome": str(args.genome), "work": str(args.work),
                                "out": str(args.out)},
        "spikes": spikes, "arms": arms,
        "wall_total": round(time.perf_counter() - t0, 1),
        "timings": timings,
    }, indent=2, default=str))
    print(f"{len(rows)} rows -> {args.out}/rows.jsonl", flush=True)


if __name__ == "__main__":
    main()
