"""How much of the Berlin library is ribosomal RNA, measured positively.

`results/peptide_absence_noise/REPORT.md` §11 withdrew the claim that peptide-k-mer
absence on real background means "not a novel virus", and left open what the still-absent
residue actually is. Two mechanisms were on the table and they imply opposite things:

1. Rare and unsequenced. A known-but-undersequenced organism, or a genuinely novel one.
   No protein index reachable from here would help, because the peptide is not in any
   database yet.
2. Abundant and non-coding. A read off ribosomal RNA, a UTR or a structural RNA can never
   match a *protein* index however complete, and in a total-RNA library rRNA is the most
   abundant thing present, not the rarest. If this dominates, the instrument is wrong for
   most of the residue rather than the residue being interesting.

The mate-separation thread is separating these by mapping the residue to the sample's own
assembly and asking whether a read falls inside a Prodigal-predicted CDS. That readout is
an *absence*: "Prodigal called no gene here". Prodigal is a prokaryotic gene caller run
over a mixed metatranscriptome, so it will miss spliced eukaryotic genes, short or
atypical viral ORFs, and anything with unusual coding statistics, and every miss inflates
the non-coding bin. This project has already lost one claim to a gene caller's miss rate.

This script measures the same branch the other way round, by *positive* identification: a
read is called rRNA when it shares exact nucleotide k-mers with SILVA, the curated rRNA
reference. "This read IS rRNA" does not depend on any gene caller, so the two readouts
fail independently.

Direction of the error, which is what makes the number usable. Exact k=31 matching
demands near-identity, and SILVA NR99 is clustered at 99% so divergent environmental
lineages are represented only by their nearest curated relative. Reads off an rRNA
molecule whose closest SILVA entry is a few percent away will therefore be missed. Every
such miss lowers the count, so the figure here is a LOWER BOUND on the library's rRNA
content. That is the useful direction: if the lower bound is already large, the abundant
non-coding branch is established without needing the bound to be tight.

The index is subsampled by a hash of the k-mer rather than truncated or clustered, so the
retained k-mers are a uniform random sample of SILVA's and the per-read hit count scales
by the sampling rate instead of being biased toward any taxon. Membership within the
retained set is exact — the 31-mer's own 62-bit code is stored, not a hash of it — so
there are no false positives from the data structure, only the vanishing chance that a
background read shares a true 31-mer with rRNA.

Reads are selected to match `experiments/absence_residual.py` on the other branch exactly
(seed s takes pairs [400000*(s-1), 400000*(s-1)+20000) in file order, mate 2 trimmed of
its 14 nt non-template prefix), and the per-read TSV is keyed on read_id and mate, so the
two tables join row for row and rRNA status can be crossed with mapping, coverage and CDS
status.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
K = 31

# 2-bit codes; everything else is 255 and invalidates every window covering it. SILVA is
# RNA, so U must map alongside T, and the files carry lower case in masked regions.
_MAP = np.full(256, 255, dtype=np.uint8)
for ch, code in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
    _MAP[ord(ch)] = code
    _MAP[ord(ch.lower())] = code


def codes_from_bases(bases: np.ndarray) -> np.ndarray:
    """Canonical 62-bit codes of every valid K-mer in a 2-bit base array (255 = invalid).

    Rolling accumulation rather than a sliding-window view: a window view of a 5 Mbp chunk
    is 31x the memory of the chunk, and the accumulation is K passes over one array.
    """
    n = bases.size - K + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)

    valid = bases != 255
    ok = np.ones(n, dtype=bool)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    b64 = np.where(valid, bases, 0).astype(np.uint64)

    for j in range(K):
        sl = b64[j : j + n]
        # forward: base j is the most significant pair
        fwd = (fwd << np.uint64(2)) | sl
        # reverse complement: base j is the least significant pair, complemented
        rev = rev | ((np.uint64(3) - sl) << np.uint64(2 * j))
        ok &= valid[j : j + n]

    canon = np.minimum(fwd, rev)
    return canon[ok]


def _mix(x: np.ndarray) -> np.ndarray:
    """splitmix64 finalizer, so the subsampling decision is unbiased across k-mer codes."""
    x = x.copy()
    x ^= x >> np.uint64(30)
    x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27)
    x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return x


def keep(codes: np.ndarray, rate_log2: int) -> np.ndarray:
    if rate_log2 == 0:
        return codes
    mask = np.uint64((1 << rate_log2) - 1)
    return codes[(_mix(codes) & mask) == np.uint64(0)]


def build_index(fastas: list[Path], rate_log2: int, chunk_bp: int = 5_000_000) -> tuple[np.ndarray, dict]:
    """Sorted unique canonical K-mer codes retained from the SILVA references."""
    parts: list[np.ndarray] = []
    stats = {"files": {}, "total_bp": 0, "sequences": 0}

    for fa in fastas:
        seqs = 0
        bp = 0
        buf: list[bytes] = []
        buflen = 0

        def flush() -> None:
            nonlocal buf, buflen
            if not buf:
                return
            # 'N' between sequences so no K-mer spans two records
            joined = b"N".join(buf)
            bases = _MAP[np.frombuffer(joined, dtype=np.uint8)]
            parts.append(keep(codes_from_bases(bases), rate_log2))
            buf, buflen = [], 0

        with gzip.open(fa, "rb") as fh:
            cur: list[bytes] = []
            for line in fh:
                if line.startswith(b">"):
                    if cur:
                        s = b"".join(cur)
                        buf.append(s)
                        buflen += len(s)
                        bp += len(s)
                        seqs += 1
                        cur = []
                        if buflen >= chunk_bp:
                            flush()
                else:
                    cur.append(line.strip())
            if cur:
                s = b"".join(cur)
                buf.append(s)
                bp += len(s)
                seqs += 1
            flush()

        stats["files"][fa.name] = {"sequences": seqs, "bases": bp}
        stats["sequences"] += seqs
        stats["total_bp"] += bp

    idx = np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.uint64)
    stats["kmers_retained"] = int(idx.size)
    stats["subsample_rate"] = 1.0 / (1 << rate_log2)
    return idx, stats


def read_pairs(r1: Path, r2: Path, start: int, count: int, trim2: int):
    """Yield (name, seq1, seq2) for pairs [start, start+count) in file order."""
    with r1.open() as f1, r2.open() as f2:
        for i in range(start):
            for _ in range(4):
                if not f1.readline():
                    return
                f2.readline()
        for _ in range(count):
            h1 = f1.readline()
            if not h1:
                return
            s1 = f1.readline().strip()
            f1.readline()
            f1.readline()
            f2.readline()
            s2 = f2.readline().strip()
            f2.readline()
            f2.readline()
            name = h1[1:].split()[0]
            yield name, s1, s2[trim2:]


def hits(seq: str, idx: np.ndarray, rate_log2: int) -> tuple[int, int]:
    """(retained K-mers queried, how many were in the index) for one read."""
    bases = _MAP[np.frombuffer(seq.encode(), dtype=np.uint8)]
    q = keep(codes_from_bases(bases), rate_log2)
    if q.size == 0 or idx.size == 0:
        return int(q.size), 0
    pos = np.searchsorted(idx, q)
    pos[pos >= idx.size] = 0
    return int(q.size), int(np.count_nonzero(idx[pos] == q))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="selects the read block, matching experiments/absence_residual.py: "
                         "pairs [400000*(seed-1), +20000)")
    ap.add_argument("--pairs", type=int, default=20000)
    ap.add_argument("--block", type=int, default=400000)
    ap.add_argument("--trim-mate2", type=int, default=14,
                    help="the non-template prefix found in PR #4; matches the other branch")
    ap.add_argument("--subsample-log2", type=int, default=3,
                    help="keep 1 in 2^this of SILVA's K-mers. 3 leaves ~9 expected hits in a "
                         "107 nt rRNA read, enough for a >=3 rule with room to spare")
    ap.add_argument("--min-hits", type=int, default=3,
                    help="retained K-mer hits needed to call a read rRNA")
    ap.add_argument("--silva", type=Path, default=ROOT / "data" / "silva")
    ap.add_argument("--reads", type=Path, default=ROOT / "data" / "berlin")
    ap.add_argument("--index-cache", type=Path,
                    default=ROOT / "data" / "silva" / "index_k31.npy")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)

    cache_meta = args.index_cache.with_suffix(".meta.json")
    if args.index_cache.exists() and cache_meta.exists():
        idx = np.load(args.index_cache)
        istats = json.loads(cache_meta.read_text())
        print(f"index from cache: {idx.size:,} k-mers")
    else:
        fastas = sorted(args.silva.glob("*.fasta.gz"))
        if not fastas:
            raise SystemExit(f"no SILVA fasta.gz under {args.silva}")
        print(f"building index from {[f.name for f in fastas]} ...")
        idx, istats = build_index(fastas, args.subsample_log2)
        np.save(args.index_cache, idx)
        cache_meta.write_text(json.dumps(istats, indent=2))
        print(f"index built: {idx.size:,} k-mers from {istats['total_bp']:,} bp")

    r1 = args.reads / "prefix_R1.fastq"
    r2 = args.reads / "prefix_R2.fastq"
    start = args.block * (args.seed - 1)

    tsv = args.out / f"rrna_s{args.seed}.tsv"
    counts = {1: [0, 0], 2: [0, 0]}   # mate -> [reads, rrna]
    hit_hist: dict[int, int] = {}
    n_pairs = 0
    with tsv.open("w") as out:
        out.write("read_id\tmate\tkmers_queried\tkmer_hits\trrna\n")
        for name, s1, s2 in read_pairs(r1, r2, start, args.pairs, args.trim_mate2):
            n_pairs += 1
            rid = name.decode() if isinstance(name, bytes) else name
            for mate, seq in ((1, s1), (2, s2)):
                q, h = hits(seq, idx, args.subsample_log2)
                flag = int(h >= args.min_hits)
                counts[mate][0] += 1
                counts[mate][1] += flag
                hit_hist[h] = hit_hist.get(h, 0) + 1
                out.write(f"{rid}\t{mate}\t{q}\t{h}\t{flag}\n")

    summary = {
        "seed": args.seed,
        "pairs": n_pairs,
        "read_block": [start, start + n_pairs],
        "min_hits": args.min_hits,
        "subsample_rate": 1.0 / (1 << args.subsample_log2),
        "k": K,
        "index": istats,
        "rrna_fraction": {
            f"mate{m}": (counts[m][1] / counts[m][0] if counts[m][0] else 0.0)
            for m in (1, 2)
        },
        "rrna_fraction_all_reads": (
            (counts[1][1] + counts[2][1]) / (counts[1][0] + counts[2][0])
            if counts[1][0] else 0.0),
        "hit_count_histogram": dict(sorted(hit_hist.items())[:30]),
        "note": "a LOWER BOUND: exact k=31 against a 99%-clustered reference misses "
                "divergent environmental rRNA lineages",
    }
    (args.out / f"rrna_s{args.seed}.json").write_text(json.dumps(summary, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/rrna_content.py",
        "accession": "SRR23973193",
        "reference": "SILVA 138.1 SSURef_NR99 + LSURef_NR99",
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))

    print(f"seed {args.seed}: {n_pairs:,} pairs, rRNA "
          f"mate1 {summary['rrna_fraction']['mate1']:.4f} "
          f"mate2 {summary['rrna_fraction']['mate2']:.4f} "
          f"all {summary['rrna_fraction_all_reads']:.4f} "
          f"({time.perf_counter() - t0:.0f} s)")


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
