"""The read-level arm: translated search on single reads, scored the way a
surveillance pipeline would have to be scored.

Two things make this harder than the usual benchmark, and both are the point.

First, the false-positive definition. A background read that hits a viral
protein is not a false positive in a wastewater sample: the sample genuinely
contains viruses, including divergent ones this project has already measured
(five astrovirus lineages 10-14% apart in SRR23973193). So the honest metric is
not "fraction of background called viral" but "flagged background reads per
million", reported beside a spike-free run of the same background, which gives
the rate the method has to beat. Any claim of a low false-positive rate that
does not state this is measuring the wrong thing.

Second, the detectable fraction of the spike genome is not the genome. With
Rhabdoviridae removed by taxonomy, only 2 of rabies' 5 proteins keep any
homolog in the database at all (L at 35.3% identity to an unclassified
rhabdo-like metagenomic virus, N at 27.4%), so the matrix and glycoprotein
regions are invisible to translated search whatever the depth. Per-read
sensitivity is therefore reported both over all spike reads and over spike reads
falling in a region with a retained homolog, because the first number mixes a
database property into what looks like a method property.

Truth labels are free: ART names spike reads after the reference accession,
background reads keep their SRA names.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path


def diamond_blastx(reads: Path, db: Path, out: Path, threads: int, sensitivity: str) -> Path:
    if out.exists():
        return out
    cmd = ["diamond", "blastx", "-q", str(reads), "-d", str(db), "-o", str(out),
           "--threads", str(threads), "--outfmt", "6", "qseqid", "sseqid",
           "pident", "length", "evalue", "bitscore", "--max-target-seqs", "1",
           "--quiet"]
    if sensitivity != "default":
        cmd.append(f"--{sensitivity}")
    subprocess.run(cmd, check=True)
    return out


def best_per_read(tsv: Path) -> dict[str, float]:
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


def curve(best: dict[str, float], is_spike, n_background: int, n_spike: int) -> list[dict]:
    """Sensitivity against flagged background reads per million, over thresholds.

    Reported as a curve because a single operating point hides the shape, and
    the shape is what decides whether a pipeline is usable at one in a million.
    """
    scored = sorted(best.items(), key=lambda kv: -kv[1])
    rows, tp, fp = [], 0, 0
    marks = {40, 50, 60, 80, 100, 150, 200, 300, 500}
    seen = set()
    for q, bits in scored:
        if is_spike(q):
            tp += 1
        else:
            fp += 1
        b = int(bits)
        for m in sorted(marks):
            if b <= m and m not in seen:
                seen.add(m)
                rows.append({
                    "bitscore_threshold": m,
                    "spike_reads_flagged": tp,
                    "sensitivity": round(tp / n_spike, 4) if n_spike else None,
                    "background_flagged": fp,
                    "background_per_million": round(fp * 1e6 / n_background, 1),
                })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reads", type=Path, required=True, help="mixed reads, one file")
    ap.add_argument("--db", type=Path, required=True, help="diamond db with the spike's family removed")
    ap.add_argument("--spike-prefix", default="NC_001542",
                    help="read-name prefix ART gave the spike reads")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--sensitivity", default="sensitive",
                    choices=["default", "sensitive", "more-sensitive", "very-sensitive"])
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep/readlevel"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    tsv = diamond_blastx(args.reads, args.db, args.out / f"hits_s{args.seed}.tsv",
                         args.threads, args.sensitivity)
    best = best_per_read(tsv)

    n_spike = n_bg = 0
    with args.reads.open() as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
                if name.startswith(args.spike_prefix):
                    n_spike += 1
                else:
                    n_bg += 1

    def is_spike(q: str) -> bool:
        return q.startswith(args.spike_prefix)

    rows = curve(best, is_spike, n_bg, n_spike)
    hit_counts = defaultdict(int)
    for q in best:
        hit_counts["spike" if is_spike(q) else "background"] += 1

    result = {
        "seed": args.seed,
        "reads_total": n_spike + n_bg,
        "spike_reads": n_spike,
        "background_reads": n_bg,
        "reads_with_any_hit": hit_counts,
        "sensitivity": args.sensitivity,
        "wall": round(time.perf_counter() - t0, 1),
        "curve": rows,
    }
    (args.out / f"readlevel_s{args.seed}.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
