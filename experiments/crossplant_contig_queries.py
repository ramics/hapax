"""Build the contig arm's query file: tier 2 contigs as windows, with an identified control.

A single 90 nt read at depth one is a weak probe of whether a molecule is background
elsewhere. A contig is a strong one, because it offers many independent windows of the same
molecule and the fraction of them witnessed in another plant is a readout rather than a
coin flip. So the same subtraction question is asked of the assembled-but-unidentifiable
pile (tier 2) as of the reads, over the same subject panel at the same fixed depth.

Each contig is cut into non-overlapping windows of `--window` nt, which are paired
consecutively so the scan sees two real windows of one molecule per query row. Both
"mates" are the same length and come from the same molecule, so unlike a read pair they
*are* interchangeable here - the mate asymmetry that every read-level number in this repo
carries does not apply, and the two mates are a replicate rather than a second instrument.
Long contigs are capped at `--max-pairs` pairs taken at an even stride along the whole
contig, not the first few, so a 6 kb contig is not represented only by its 5' end.

Identified contigs ride along as the comparator: they are the same assembly, the same
lengths and the same instrument, and they are known shared background, so they say what
the witnessing rate reads on material that certainly recurs. A mononucleotide shuffle of
part of each draw is the null at the other end.

Contigs are sampled per seed rather than taken whole, because the five-seed range is the
only honest dispersion this repo reports and a whole population has none.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from subtract_known import provenance  # noqa: E402

SETS = {"unident": "sample_contigs_unidentified.fa", "ident": "sample_contigs_identified.fa"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def windows(seq: str, w: int, max_pairs: int) -> list[tuple[str, str]]:
    """Consecutive window pairs, at an even stride along the contig, ACGT only."""
    n = len(seq) // w
    cells = [seq[i * w:(i + 1) * w] for i in range(n)]
    cells = [c for c in cells if set(c) <= set("ACGT")]
    pairs = [(cells[i], cells[i + 1]) for i in range(0, len(cells) - 1, 2)]
    if len(pairs) > max_pairs:
        step = (len(pairs) - 1) / (max_pairs - 1) if max_pairs > 1 else 0
        pairs = [pairs[round(i * step)] for i in range(max_pairs)]
    return pairs


def shuffled(seq: str, rng: random.Random) -> str:
    c = list(seq)
    rng.shuffle(c)
    return "".join(c)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--contig-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--window", type=int, default=90,
                    help="same width as the read queries, so a contig window and a read "
                         "offer the same number of 31-mers and neither gets more chances")
    ap.add_argument("--max-pairs", type=int, default=8)
    ap.add_argument("--unident", type=int, default=6000, help="contigs per seed")
    ap.add_argument("--ident", type=int, default=3000, help="control contigs per seed")
    ap.add_argument("--null-frac", type=float, default=0.1)
    ap.add_argument("--seeds", default="1,2,3,4,5")
    args = ap.parse_args()
    t0 = time.time()

    pools, lens = {}, {}
    for key, fn in SETS.items():
        got = read_fasta(args.contig_dir / fn)
        pools[key] = got
        lens[key] = [len(s) for _, s in got]

    want = {"unident": args.unident, "ident": args.ident}
    rows, stats = [], {}
    for seed in (int(s) for s in args.seeds.split(",")):
        for key, pool in pools.items():
            rng = random.Random(f"{key}|{seed}")
            eligible = [i for i, (_, s) in enumerate(pool)
                        if len(s) >= 2 * args.window and set(s) <= set("ACGTN")]
            take = min(want[key], len(eligible))
            picked = rng.sample(eligible, take)
            n_pairs = 0
            for ci in picked:
                name, seq = pool[ci]
                ws = windows(seq, args.window, args.max_pairs)
                for p, (a, b) in enumerate(ws):
                    tag = f"contig|{key}|{ci}|L{len(seq)}|p{p}|s{seed}"
                    rows.append((tag, a, b))
                    if rng.random() < args.null_frac:
                        rows.append((f"shuffled-{tag}", shuffled(a, rng), shuffled(b, rng)))
                    n_pairs += 1
            st = stats.setdefault(key, {"eligible": len(eligible), "contigs_total": len(pool),
                                        "per_seed": {}})
            st["per_seed"][seed] = {"contigs": take, "pairs": n_pairs}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    manifest = {"script": "experiments/crossplant_contig_queries.py", **provenance(__file__),
                "parameters": {"window": args.window, "max_pairs": args.max_pairs,
                               "unident_per_seed": args.unident,
                               "ident_per_seed": args.ident,
                               "null_frac": args.null_frac, "seeds": args.seeds},
                "source": "s3 shared/contig_sets/, the cascade's own contig sets",
                "length_summary": {k: {"contigs": len(v), "total_nt": sum(v),
                                       "median_nt": sorted(v)[len(v) // 2]}
                                   for k, v in lens.items()},
                "sets": stats, "rows": len(rows),
                "seconds": round(time.time() - t0, 1)}
    (args.out.parent / "contig_queries_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"rows": len(rows), "sets": stats}, indent=2))


if __name__ == "__main__":
    main()
