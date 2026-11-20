"""Two related lineages at fixed total coverage: does pooling buy length, and
what does it cost in mosaicism?

This closes the loop between this project's two measurements. The abundance
sweep located where contigs cross a classifier's useful length (1,500 nt between
5x and 7.5x for one lineage). The stress test measured the opposite failure: a
minor member erased below 3% divergence, its reads absorbed into the majority's
contig. The sister thread's point is that these are the same event seen twice.
An assembler that joins reads across a 2% difference sees two lineages at 3.75x
each as one sequence at 7.5x, and 7.5x is exactly where length arrives.

So pooling is not a bug to be fixed or a feature to be enabled. It converts
coverage into length at the price of a mosaic, and whether that is a gain or a
loss depends entirely on the question being asked. That argues for labelling it,
not preventing it, which is what the per-contig novelty/collapse flag would do.

Prediction being tested, stated before running: at fixed 7.5x total, an
assembler unwilling to pool returns about 25-30% breadth per lineage (each
lineage behaving like the sweep's 3.35x rung), while one willing to pool returns
about 91% breadth of a mosaic (the sweep's 7.54x rung).

Truth is real published sequence throughout, as the repo requires: four pairs of
complete rabies genomes whose divergence is what it is, measured with edlib and
reported, never simulated evolution. Lineage assignment is per 300 nt window, so
a contig that switches lineage is counted as the mosaic it is rather than
assigned wholesale to its majority.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import edlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abundance_sweep import make_background, sh  # noqa: E402
from common import git_sha  # noqa: E402

# Real pairs, divergence measured with edlib over the whole genome. No pair sits
# at 5%: the deposited rabies genomes cluster, so the ladder is 2.0, 3.8, 7.3,
# 9.7% rather than the 2/5/10 asked for, and the gap is reported rather than
# filled by mutating a genome to order.
PAIRS = [
    ("div02", "PZ227608.1", "PZ227621.1", 0.0198),
    ("div04", "PV491557.1", "PZ324370.1", 0.0381),
    ("div07", "PZ227608.1", "AB044824.1", 0.0729),
    ("div10", "NC_001542.1", "AB009663.2", 0.0965),
]
WINDOW = 300


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.open():
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def write_one(seqs: dict[str, str], acc: str, out: Path) -> Path:
    key = next(k for k in seqs if k.startswith(acc.split(".")[0]))
    out.write_text(f">{key}\n{seqs[key]}\n")
    return out


def simulate(genome: Path, out: Path, pairs_wanted: int, read_len: int, seed: int) -> tuple[Path, Path]:
    """ART at the coverage that yields the wanted number of pairs, then truncate.

    Generating a surplus and taking a prefix keeps rungs nested and makes the
    read count exact rather than approximate, which matters when the whole point
    is a fixed coverage budget.
    """
    out.mkdir(parents=True, exist_ok=True)
    pref = out / "r"
    p1, p2 = Path(str(pref) + "1.fq"), Path(str(pref) + "2.fq")
    glen = len(next(iter(read_fasta(genome).values())))
    cov = max(2.0, (pairs_wanted * 2 * read_len / glen) * 1.4)
    if not (p1.exists() and p2.exists()):
        sh(["art_illumina", "-ss", "MSv3", "-i", str(genome), "-p", "-l", str(read_len),
            "-f", str(round(cov, 2)), "-m", "250", "-s", "40", "-rs", str(seed),
            "-na", "-o", str(pref)])
    return p1, p2


def lineage_windows(contig: str, a: str, b: str) -> tuple[int, int, int]:
    """Per-window lineage call: (windows_a, windows_b, switches).

    Each window is aligned to both lineages in infix mode and assigned to the
    closer one; ties and windows that match neither are not counted. Switches
    count transitions along the contig, so a clean contig of either lineage has
    zero and a mosaic has at least one.
    """
    calls = []
    for i in range(0, len(contig) - WINDOW + 1, WINDOW):
        w = contig[i:i + WINDOW]
        da = edlib.align(w, a, mode="HW", task="distance")["editDistance"]
        db = edlib.align(w, b, mode="HW", task="distance")["editDistance"]
        if da < db:
            calls.append("a")
        elif db < da:
            calls.append("b")
    switches = sum(1 for x, y in zip(calls, calls[1:]) if x != y)
    return calls.count("a"), calls.count("b"), switches


def score(contigs: Path, ref2: Path, a: str, b: str) -> dict:
    """Spike-derived contigs, per-lineage breadth, and mosaicism.

    Breadth is computed per lineage from the alignments to that lineage, so two
    lineages recovered separately and one mosaic covering both read differently
    here, which is the distinction the whole rung exists to make.
    """
    paf = sh(["minimap2", "-x", "asm20", "-c", "--secondary=no", str(ref2), str(contigs)]).stdout
    aligned: dict[str, int] = {}
    qlen: dict[str, int] = {}
    covered = {"a": set(), "b": set()}
    names = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        target, ts, te = f[5], int(f[7]), int(f[8])
        matches, blocklen = int(f[9]), int(f[10])
        if blocklen == 0 or matches / blocklen < 0.85:
            continue
        qlen[q] = ql
        aligned[q] = aligned.get(q, 0) + (qe - qs)
        side = "a" if target.startswith(a.split(".")[0]) else "b"
        covered[side].update(range(ts, te))
        names[q] = side
    spike = [q for q in aligned if aligned[q] >= 0.8 * qlen[q]]
    seqs = read_fasta(contigs)
    refs = read_fasta(ref2)
    sa = next(v for k, v in refs.items() if k.startswith(a.split(".")[0]))
    sb = next(v for k, v in refs.items() if k.startswith(b.split(".")[0]))
    lens = sorted((qlen[q] for q in spike), reverse=True)
    mosaic, switches = 0, 0
    for q in spike:
        if qlen[q] < 2 * WINDOW:
            continue
        wa, wb, sw = lineage_windows(seqs[q], sa, sb)
        if wa and wb and min(wa, wb) >= 1 and sw:
            mosaic += 1
        switches += sw
    return {
        "spike_contigs": len(spike),
        "spike_mass": sum(lens),
        "longest": lens[0] if lens else 0,
        "mass_over_1500": sum(x for x in lens if x >= 1500),
        "breadth_a": round(len(covered["a"]) / len(sa), 4),
        "breadth_b": round(len(covered["b"]) / len(sb), 4),
        "mosaic_contigs": mosaic,
        "lineage_switches": switches,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=3)
    ap.add_argument("--bg1", type=Path, required=True)
    ap.add_argument("--bg2", type=Path, required=True)
    ap.add_argument("--bg-pairs", type=int, default=1_000_000)
    ap.add_argument("--pool-fa", type=Path, required=True, help="multi-fasta holding the lineages")
    ap.add_argument("--extra-fa", type=Path, help="second fasta searched for accessions too")
    ap.add_argument("--total-pairs", type=int, default=450,
                    help="total spike pairs per run; 450 is about 7.5x of a 11.9 kb genome")
    ap.add_argument("--read-len", type=int, default=107)
    ap.add_argument("--r2-trim-start", type=int, default=15)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--controls", action="store_true",
                    help="also run each lineage alone at half the budget")
    ap.add_argument("--only", default="",
                    help="comma-separated pair tags to run; the rest are skipped, so a run "
                         "interrupted partway can resume without duplicating rows")
    ap.add_argument("--work", type=Path, default=Path("/home/user/sweepdata/twolin"))
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep"))
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    r2_len = args.read_len - args.r2_trim_start + 1
    bg1, bg2 = make_background(args.bg1, args.bg2, args.bg_pairs,
                               Path("/home/user/sweepdata/work/bg"), args.seed,
                               args.r2_trim_start, args.read_len)

    seqs = read_fasta(args.pool_fa)
    if args.extra_fa:
        seqs |= read_fasta(args.extra_fa)
    rows, t0 = [], time.perf_counter()
    half = args.total_pairs // 2

    only = {t for t in args.only.split(",") if t}
    for tag, acc_a, acc_b, div in PAIRS:
        if only and tag not in only:
            continue
        fa_a = write_one(seqs, acc_a, args.work / f"{tag}_a.fa")
        fa_b = write_one(seqs, acc_b, args.work / f"{tag}_b.fa")
        ref2 = args.work / f"{tag}_both.fa"
        ref2.write_text(fa_a.read_text() + fa_b.read_text())
        # measured here rather than trusted from the table above
        sa = next(iter(read_fasta(fa_a).values()))
        sb = next(iter(read_fasta(fa_b).values()))
        div_measured = edlib.align(sa, sb, task="distance")["editDistance"] / max(len(sa), len(sb))

        for seed in range(args.seed, args.seed + args.n_seeds):
            a1, a2 = simulate(fa_a, args.work / f"{tag}_a_s{seed}", half, args.read_len, seed)
            b1, b2 = simulate(fa_b, args.work / f"{tag}_b_s{seed}", half, args.read_len, seed + 500)
            modes = [("both", [(a1, a2, half), (b1, b2, half)])]
            if args.controls:
                modes += [("a_alone", [(a1, a2, half)]), ("b_alone", [(b1, b2, half)])]
            for mode, parts in modes:
                m1 = args.work / f"mix_{tag}_{mode}_s{seed}_1.fq.gz"
                m2 = args.work / f"mix_{tag}_{mode}_s{seed}_2.fq.gz"
                for mix, bg, idx in ((m1, bg1, 0), (m2, bg2, 1)):
                    cmd = f"cat {bg}"
                    for p, n in [((p1, p2)[idx], n) for p1, p2, n in parts]:
                        cut = f"seqkit head -n {n} {p}"
                        if idx == 1:
                            cut += f" | seqkit subseq -r 1:{r2_len}"
                        cmd += f" <({cut} | gzip -c)"
                    with open(mix, "wb") as fh:
                        subprocess.run(["bash", "-c", cmd], stdout=fh, check=True)
                asm = args.work / f"asm_{tag}_{mode}_s{seed}"
                if asm.exists():
                    subprocess.run(["rm", "-rf", str(asm)], check=True)
                try:
                    sh(["megahit", "-1", str(m1), "-2", str(m2), "-o", str(asm),
                        "-t", str(args.threads), "--min-contig-len", "200"])
                except subprocess.CalledProcessError as e:
                    print(f"[{tag} {mode} s{seed}] megahit failed: {e.stderr[-300:]}", flush=True)
                    continue
                sc = score(asm / "final.contigs.fa", ref2, acc_a, acc_b)
                total = sum(n for _, _, n in parts)
                row = {"pair": tag, "lineage_a": acc_a, "lineage_b": acc_b,
                       "divergence": round(div_measured, 4), "mode": mode, "seed": seed,
                       "spike_pairs_total": total,
                       "coverage_total": round(total * (args.read_len + r2_len) / len(sa), 2),
                       **sc}
                rows.append(row)
                print(json.dumps(row), flush=True)
                (args.out / "two_lineage_rows.jsonl").open("a").write(json.dumps(row) + "\n")
                subprocess.run(["rm", "-rf", str(asm)], check=True)
                m1.unlink(missing_ok=True)
                m2.unlink(missing_ok=True)

    (args.out / "two_lineage_manifest.json").write_text(json.dumps({
        "git_sha": git_sha(),
        "params": {k: str(v) for k, v in vars(args).items()},
        "pairs": PAIRS, "window": WINDOW,
        "wall_total": round(time.perf_counter() - t0, 1),
    }, indent=2))
    print(f"{len(rows)} rows", flush=True)


if __name__ == "__main__":
    main()
