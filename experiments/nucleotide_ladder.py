"""Does the nucleotide side of the removal curve flatten? A content ladder over reference
taxa, from a classification run already on disk.

WHY THIS EXISTS. The cascade's nucleotide tiers gave three index sizes - Kraken2 viral at
0.67 GB, Standard/PlusPF-8 at 8 GB, PlusPF-16 at 16 GB - and the last doubling bought 0.42
points of marginal removal. Read as a curve that says public nucleotide catalogues are
nearly exhausted. But the protein thread made exactly that mistake in the other direction
(two small virus-weighted points, read as "protein databases barely help", overturned when
a tenth of UniRef50 moved the same quantity 23 points), so three points on a line are not a
curve here either. This script puts more points on the nucleotide axis.

WHAT IT MEASURES, AND WHAT IT CANNOT. Building real Kraken2 databases from nested genome
subsets needs the source library (hundreds of GB) and a build step this environment cannot
run, so the ladder is constructed the other way round: classify once against the largest
database available, then ask what removal would have been if only a nested random subset of
REFERENCE TAXA had been present. The axis is therefore reference content, not bytes.

Two biases, both in the same direction, both stated rather than corrected:

1. A real database built from a taxon subset would still classify many reads whose own
   taxon was dropped, by assigning them to a retained relative at a higher rank. Deleting
   calls, as done here, loses those. So each point UNDERSTATES what a database of that
   content would remove, which makes the curve here steeper at the top than the truth and
   biases this script TOWARDS the conclusion it was written to challenge.
2. Kraken2's hash is capped, and minimizers are shared between related genomes, so bytes
   grow sublinearly with taxa. A bytes axis would compress the right-hand side.

Read the output for shape only: whether marginal removal per added fraction of reference
content is still rising at the largest content available, which is the question the protein
curve answered with "not flattened yet".
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from subtract_known import (_mm_index, dir_bytes, git_sha, kraken2, map_filter, spread,
                            write_fastq)

ROOT = Path(__file__).resolve().parents[1]


def taxon_kept(taxid: int, frac: float, salt: int) -> bool:
    """Deterministic nested membership: a taxon kept at frac is kept at every larger frac.

    Nested rather than independently drawn, because the question is what ADDING reference
    content buys; independent draws per rung would mix in the variance of which taxa were
    picked and make the marginal column meaningless.
    """
    h = (taxid * 2654435761 + salt * 40503) & 0xFFFFFFFF
    h ^= h >> 15
    h = (h * 2246822519) & 0xFFFFFFFF
    h ^= h >> 13
    return (h % 1_000_000) < frac * 1_000_000


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--fractions", type=float, nargs="+",
                    default=[0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1.0])
    ap.add_argument("--confidence", type=float, nargs="+", default=[0.0, 0.1])
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--silva", type=Path, default=ROOT / "data/silva/silva_dna.fa")
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/nucladder"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/nucleotide_ladder")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    rows = list(read_pairs(args.r1, args.r2, args.pairs, 14,
                           (args.seed - 1) * 1_500_000))
    seqs = {r[0]: (r[1], r[2]) for r in rows}
    n = len(rows)

    # rRNA first, as a design requirement rather than an ordering accident: against any
    # large index ribosomal reads stop looking novel for annotation reasons (automatically
    # annotated entries that are really translated ribosomal sequence), so a third of this
    # library would otherwise be scoring database annotation quality.
    rrna = set(map_filter(_mm_index(args.silva, args.work / "silva.mmi", args.threads),
                          seqs, args.min_mlen))

    fq1, fq2 = args.work / "r1.fq", args.work / "r2.fq"
    write_fastq(fq1, rows, 1)
    write_fastq(fq2, rows, 2)

    out_conf = {}
    for conf in args.confidence:
        calls, wall, _rss = kraken2(args.db, fq1, fq2, args.classify_bin, args.threads,
                                    args.work / "k2.out", conf)
        called = {k: t for k, t in calls.items() if t != 0}
        taxa = sorted(set(called.values()))
        rungs, prev_all, prev_nr = [], 0.0, 0.0
        for frac in args.fractions:
            keep = {t for t in taxa if taxon_kept(t, frac, args.seed)}
            removed = {k for k, t in called.items() if t in keep}
            all_frac = len(removed) / n
            nr_frac = len(removed - rrna) / max(n - len(rrna), 1)
            rungs.append({
                "content_fraction": frac,
                "reference_taxa_retained": len(keep),
                "removed_frac_of_library": all_frac,
                "marginal_frac_of_library": all_frac - prev_all,
                "removed_frac_of_non_rrna": nr_frac,
                "marginal_frac_of_non_rrna": nr_frac - prev_nr,
            })
            prev_all, prev_nr = all_frac, nr_frac
        out_conf[str(conf)] = {
            "distinct_reference_taxa_called": len(taxa),
            "classified_frac": len(called) / n,
            "rungs": rungs,
            "kraken2_wall_seconds": wall,
        }

    out = {
        "seed": args.seed,
        "pairs": n,
        "rrna_frac": len(rrna) / n,
        "db": str(args.db),
        "db_bytes": dir_bytes(args.db),
        "by_confidence": out_conf,
        "caveats": [
            "content axis, not bytes: a real DB of this content would be smaller than "
            "db_bytes times the fraction, because minimizers are shared between relatives",
            "deleting calls understates a smaller DB: reads whose taxon was dropped would "
            "often still be classified to a retained relative at a higher rank",
        ],
        "wall_seconds": time.perf_counter() - t0,
    }
    (args.out / f"nucladder_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/nucleotide_ladder.py",
        "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }, indent=2))

    for conf, d in out_conf.items():
        print(f"confidence {conf}: {d['distinct_reference_taxa_called']} taxa, "
              f"classified {d['classified_frac'] * 100:.2f}%")
        for r in d["rungs"]:
            print(f"  content {r['content_fraction']:<6} taxa {r['reference_taxa_retained']:>6}"
                  f"  removed {r['removed_frac_of_library'] * 100:6.2f}%"
                  f"  marginal {r['marginal_frac_of_library'] * 100:6.2f}%"
                  f"  | non-rRNA removed {r['removed_frac_of_non_rrna'] * 100:6.2f}%"
                  f"  marginal {r['marginal_frac_of_non_rrna'] * 100:6.2f}%")


if __name__ == "__main__":
    main()
