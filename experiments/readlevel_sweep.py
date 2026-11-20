"""Read-level detection across the same coverage rungs as the assembly sweep.

This is the arm that matters below about 2x spike coverage, where the assembly
sweep returns no contig at all and geNomad therefore has nothing to call.

Two efficiencies make it cheap, and both are exact rather than approximate.
The background's translated search is run once, because the background is held
fixed across every rung; and the spike pool's search is run once, because the
rungs are nested subsamples of that pool, so a rung's spike hits are the hits of
its first n read pairs. Composing rungs from those two searches gives the same
answer as searching each mixture, for a fraction of the compute.

The metric is sensitivity against flagged background reads per million, not a
false-positive rate: this background is a real wastewater virome and genuinely
contains divergent viruses, so a background read with a viral hit is usually a
true viral read. The spike-free rate is the number any claim has to beat, and
reporting it as a rate keeps that honest. The companion finding from the
classifier arm makes the point sharply: geNomad calls 98% of background contigs
over 3 kb from this sample viral, so "viral" is nearly uninformative here and
the question that matters is "new".
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

SPIKES = [10, 30, 100, 300, 1000, 3000]


def interleave_to_fasta(r1: Path, r2: Path, out: Path, limit_pairs: int | None = None) -> int:
    """Both mates as one fasta, mate 2 kept at its trimmed length.

    Read-level detection treats each read separately, so the pair structure is
    not used here; it is kept only so that the read counts match the assembly
    sweep's rungs exactly.
    """
    if out.exists():
        return sum(1 for line in out.open() if line.startswith(">"))
    n = 0
    cmd = f"seqkit fq2fa {r1} ; seqkit fq2fa {r2}"
    if limit_pairs is not None:
        cmd = (f"seqkit head -n {limit_pairs} {r1} | seqkit fq2fa ; "
               f"seqkit head -n {limit_pairs} {r2} | seqkit fq2fa")
    with out.open("w") as fh:
        subprocess.run(["bash", "-c", cmd], stdout=fh, check=True)
    n = sum(1 for line in out.open() if line.startswith(">"))
    return n


def diamond(query: Path, db: Path, out: Path, threads: int, sensitivity: str) -> Path:
    if out.exists():
        return out
    cmd = ["diamond", "blastx", "-q", str(query), "-d", str(db), "-o", str(out),
           "--threads", str(threads), "--outfmt", "6", "qseqid", "sseqid",
           "pident", "length", "evalue", "bitscore", "--max-target-seqs", "1",
           "--quiet"]
    if sensitivity != "default":
        cmd.append(f"--{sensitivity}")
    subprocess.run(cmd, check=True)
    return out


def best_scores(tsv: Path) -> dict[str, float]:
    best: dict[str, float] = {}
    with tsv.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, bits = f[0], float(f[5])
            if bits > best.get(q, 0.0):
                best[q] = bits
    return best


def read_order(fasta: Path) -> list[str]:
    """Spike read names in pool order, so a rung is a prefix of this list."""
    return [line[1:].split()[0] for line in fasta.open() if line.startswith(">")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bg1", type=Path, required=True)
    ap.add_argument("--bg2", type=Path, required=True)
    ap.add_argument("--pool1", type=Path, required=True)
    ap.add_argument("--pool2", type=Path, required=True)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--sensitivity", default="sensitive")
    ap.add_argument("--bits", type=float, default=50.0,
                    help="bitscore threshold for a flag; the curve reports several")
    ap.add_argument("--spikes", default=",".join(str(s) for s in SPIKES))
    ap.add_argument("--work", type=Path, default=Path("/home/user/sweepdata/readlevel"))
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep"))
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    bg_fa = args.work / "bg.fa"
    n_bg = interleave_to_fasta(args.bg1, args.bg2, bg_fa)
    pool_fa = args.work / f"pool_s{args.seed}.fa"
    interleave_to_fasta(args.pool1, args.pool2, pool_fa)
    print(f"background reads {n_bg}", flush=True)

    bg_hits = best_scores(diamond(bg_fa, args.db, args.work / "bg_hits.tsv",
                                  args.threads, args.sensitivity))
    pool_hits = best_scores(diamond(pool_fa, args.db, args.work / f"pool_s{args.seed}_hits.tsv",
                                    args.threads, args.sensitivity))
    print(f"background reads with a hit {len(bg_hits)}; pool reads with a hit {len(pool_hits)}",
          flush=True)

    names = read_order(pool_fa)
    half = len(names) // 2  # mate 1 block then mate 2 block
    thresholds = [40, 50, 60, 80, 100, 150, 200]
    rows = []
    for n in [int(s) for s in args.spikes.split(",")]:
        # a rung is the first n pairs of the pool: the first n of each mate block
        rung = names[:n] + names[half:half + n]
        for t in thresholds:
            tp = sum(1 for q in rung if pool_hits.get(q, 0) >= t)
            fp = sum(1 for b in bg_hits.values() if b >= t)
            rows.append({
                "seed": args.seed, "spike_pairs": n, "spike_reads": len(rung),
                "bitscore": t,
                "spike_flagged": tp,
                "sensitivity": round(tp / len(rung), 4),
                "background_flagged": fp,
                "background_per_million": round(fp * 1e6 / n_bg, 1),
            })
    result = {
        "seed": args.seed,
        "background_reads": n_bg,
        "background_with_hit": len(bg_hits),
        "pool_reads_with_hit": len(pool_hits),
        "pool_reads": len(names),
        "sensitivity_setting": args.sensitivity,
        "wall": round(time.perf_counter() - t0, 1),
        "rows": rows,
    }
    (args.out / f"readlevel_s{args.seed}.json").write_text(json.dumps(result, indent=2))
    for r in rows:
        if r["bitscore"] == 50:
            print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main()
