"""Does the SILVA k-mer rule actually identify rRNA, and only rRNA?

`experiments/rrna_content.py` calls a Berlin read ribosomal when at least `--min-hits` of
its retained exact 31-mers are present in a SILVA index. That rule has two failure modes
and the headline number means nothing until both are measured on data where the answer is
known in advance.

**Sensitivity, and the reason it cannot be assumed.** SILVA NR99 is clustered at 99%
identity, so an environmental rRNA lineage is represented only by its nearest curated
relative. Exact 31-mer matching needs a 31 nt stretch that survived whatever divergence
separates the read from that relative, and a few percent divergence breaks most 31-mers.
The honest test therefore holds sequences OUT of the index and samples reads from the
held-out ones: every such read is genuinely rRNA, and whether it is detected depends on a
different sequence being close enough, which is exactly the situation a real read is in.
Holding out whole sequences rather than individual k-mers is the point; a k-mer-level
holdout would leave the read's own neighbours in the index and measure nothing.

**Specificity.** Reads from sequence that is certainly not rRNA must not be called. The
negatives are the six genomes already on disk from `experiments/peptide_absence_noise.py`
(two of them bacteriophage, so they carry no rRNA of their own and no host sequence
either). A rule that flags these is picking up low-complexity or repetitive matches rather
than ribosomal identity.

Both read sets carry substitutions at Berlin's own measured rate (0.00089 on R1, 0.00117
on R2, `results/berlin_quality/`), because error degrades exact k-mer matching and a
control run on perfect sequence would overstate sensitivity.

The number this produces is a correction factor with a known direction: observed rRNA
fraction divided by sensitivity is the floor on the true fraction. If sensitivity is poor
the headline stays a lower bound and says so, rather than being quietly rescaled.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rrna_content as rc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def read_silva_with_holdout(fastas: list[Path], every: int):
    """Sequences for the index, and the held-out ones, split by position in the file."""
    kept: list[bytes] = []
    held: list[bytes] = []
    n = 0
    for fa in fastas:
        cur: list[bytes] = []
        with gzip.open(fa, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    if cur:
                        (held if n % every == 0 else kept).append(b"".join(cur))
                        n += 1
                        cur = []
                else:
                    cur.append(line.strip())
            if cur:
                (held if n % every == 0 else kept).append(b"".join(cur))
                n += 1
    return kept, held


def index_from_sequences(seqs: list[bytes], rate_log2: int, chunk_bp: int = 5_000_000):
    parts, buf, buflen = [], [], 0
    for s in seqs:
        buf.append(s)
        buflen += len(s)
        if buflen >= chunk_bp:
            bases = rc._MAP[np.frombuffer(b"N".join(buf), dtype=np.uint8)]
            parts.append(rc.keep(rc.codes_from_bases(bases), rate_log2))
            buf, buflen = [], 0
    if buf:
        bases = rc._MAP[np.frombuffer(b"N".join(buf), dtype=np.uint8)]
        parts.append(rc.keep(rc.codes_from_bases(bases), rate_log2))
    return np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.uint64)


def genbank_sequence(path: Path) -> str:
    """The ORIGIN block of a GenBank flat file, uppercased."""
    out, on = [], False
    for line in path.read_text().splitlines():
        if line.startswith("ORIGIN"):
            on = True
            continue
        if on:
            if line.startswith("//"):
                break
            out.append("".join(c for c in line if c.isalpha()))
    return "".join(out).upper()


def sample_reads(seqs: list[str], n: int, length: int, err: float, rng) -> list[str]:
    """n reads of `length`, uniform over the sequences long enough to hold one."""
    pool = [s for s in seqs if len(s) >= length]
    if not pool:
        return []
    weights = np.array([len(s) - length + 1 for s in pool], dtype=float)
    weights /= weights.sum()
    bases = np.array(list("ACGT"))
    reads = []
    for _ in range(n):
        s = pool[rng.choice(len(pool), p=weights)]
        i = rng.integers(0, len(s) - length + 1)
        read = np.array(list(s[i : i + length].replace("U", "T")))
        hit = rng.random(length) < err
        if hit.any():
            read[hit] = bases[rng.integers(0, 4, hit.sum())]
        reads.append("".join(read))
    return reads


def flag_rate(reads: list[str], idx: np.ndarray, rate_log2: int, min_hits: int) -> dict:
    hits = [rc.hits(r, idx, rate_log2)[1] for r in reads]
    h = np.array(hits) if hits else np.zeros(0)
    return {
        "reads": len(reads),
        "flag_rate": float(np.mean(h >= min_hits)) if h.size else 0.0,
        "mean_hits": float(h.mean()) if h.size else 0.0,
        "median_hits": float(np.median(h)) if h.size else 0.0,
        "zero_hit_fraction": float(np.mean(h == 0)) if h.size else 0.0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--holdout-every", type=int, default=50,
                    help="every Nth SILVA sequence is held out of the index and used as a "
                         "positive; 50 leaves the index at 98% of full size")
    ap.add_argument("--reads", type=int, default=20000)
    ap.add_argument("--length", type=int, default=107, help="Berlin's modal read length")
    ap.add_argument("--error-rate", type=float, default=0.001,
                    help="Berlin's measured substitution rate, results/berlin_quality/")
    ap.add_argument("--subsample-log2", type=int, default=3)
    ap.add_argument("--min-hits", type=int, nargs="+", default=[1, 2, 3, 5, 8],
                    help="swept, so the rule is chosen from the curve rather than asserted")
    ap.add_argument("--silva", type=Path, default=ROOT / "data" / "silva")
    ap.add_argument("--negatives", type=Path, default=ROOT / "data" / "peptide_noise")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    fastas = sorted(args.silva.glob("*.fasta.gz"))
    print("splitting SILVA ...")
    kept, held = read_silva_with_holdout(fastas, args.holdout_every)
    print(f"{len(kept):,} sequences in the index, {len(held):,} held out")
    idx = index_from_sequences(kept, args.subsample_log2)
    print(f"holdout index: {idx.size:,} k-mers")

    pos = sample_reads([s.decode() for s in held], args.reads, args.length,
                       args.error_rate, rng)
    neg_seqs = [genbank_sequence(p) for p in sorted(args.negatives.glob("*.gb"))]
    neg = sample_reads(neg_seqs, args.reads, args.length, args.error_rate, rng)
    print(f"{len(pos):,} held-out rRNA reads, {len(neg):,} non-rRNA reads")

    rows = []
    for T in args.min_hits:
        p = flag_rate(pos, idx, args.subsample_log2, T)
        n = flag_rate(neg, idx, args.subsample_log2, T)
        rows.append({
            "min_hits": T,
            "sensitivity_heldout_rrna": p["flag_rate"],
            "false_flag_rate_non_rrna": n["flag_rate"],
            "positive_detail": p,
            "negative_detail": n,
        })
        print(f"T>={T}: sensitivity {p['flag_rate']:.4f}  "
              f"false flags {n['flag_rate']:.5f}  "
              f"(positive mean hits {p['mean_hits']:.2f}, zero-hit {p['zero_hit_fraction']:.3f})")

    out = {
        "seed": args.seed,
        "holdout_every": args.holdout_every,
        "sequences_in_index": len(kept),
        "sequences_held_out": len(held),
        "kmers_in_holdout_index": int(idx.size),
        "read_length": args.length,
        "error_rate": args.error_rate,
        "negatives": [p.name for p in sorted(args.negatives.glob("*.gb"))],
        "thresholds": rows,
        "note": "sensitivity here is against a 99%-clustered reference with the read's own "
                "sequence removed, which is the situation a real environmental read is in; "
                "it is the correction factor for the headline, not a tuning knob",
    }
    (args.out / f"rrna_control_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"rrna_control_manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/rrna_control.py",
        "reference": "SILVA 138.1 SSURef_NR99 + LSURef_NR99",
        "git_sha": rc._git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))


if __name__ == "__main__":
    main()
