"""The absence rate, and the mate correlation, against protein indexes of
increasing coverage - including real public ones.

Why this exists. `results/mate_width/REPORT.md` §4.3 showed that the 86% of
Berlin background reads with no known peptide 10-mer collapses to 22% once the
sample's own predicted proteins are indexed, i.e. that the absence was mostly
public-database incompleteness. The obvious objection is that the comparison ran
against a small public index (RefSeq viral + Swiss-Prot, 0.22 Gaa), and the
honest version needs real public indexes several times larger. UniRef50 was
first reported here as unobtainable; that was wrong. The file cannot be *stored*
(the proxy will not carry 9.2 GB to disk) but it *streams* reliably, and
`subsample_uniref.py` keeps nested every-10th and every-100th subsamples whose
record count matches UniProt's own release total exactly.

That leaves the index too big to hold: ~1.2e9 distinct 10-mers for the every-10th
subsample, ~10 GB as sorted hashes, against 15 GB of RAM. So this script inverts
the lookup. Instead of building an index and probing it per read, it hashes the
reads' peptide k-mers once, sorts them, and then *streams* each proteome,
marking which of the (few million) query k-mers it contains. Memory is set by
the reads, not by the database, and nothing of the index is ever written. The
cost is one pass per proteome per run.

The sources are applied cumulatively, so each level is a superset of the last
and the levels differ only in coverage. The sample's own assembly goes last, so
the final step isolates what 8 Maa taken from the sample itself adds on top of
every public index reachable here.

Readouts per level: the fraction of reads with every k-mer absent, the mean
best-frame absence (continuous, because the fully-absent rate saturates), and
the mates' correlation, which is the quantity the pair rules' merit depends on.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts
from insert_profile import read_pairs                        # noqa: E402
from peptide_absence import AA_CODE, hash_kmers, proteome_chunks, six_frames  # noqa: E402


POP = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def distinct_residues(codes: np.ndarray, k: int) -> np.ndarray:
    """Distinct residues in every length-k window, as a complexity readout.

    Why. A peptide k-mer translated from ribosomal RNA has no true protein
    behind it, yet reads called rRNA are measurably *less* likely to be fully
    absent than other reads. One live explanation is low complexity: a
    repetitive translated window can match the low-complexity stretches of real
    proteins, which no amount of cleaning the index of spurious gene calls would
    remove. So complexity is carried per k-mer and the absence statistic is
    reported a second time over complex k-mers only, which separates the two
    explanations instead of arguing them.

    Implemented as an OR of residue bitmasks across the window, then a popcount,
    so it costs k passes over the array like the hash itself.
    """
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint8)
    bits = np.where(codes > 0, np.uint32(1) << (codes.astype(np.uint32) - 1),
                    np.uint32(0))
    acc = np.zeros(n, dtype=np.uint32)
    ok = np.ones(n, dtype=bool)
    for off in range(k):
        acc |= bits[off:off + n]
        ok &= codes[off:off + n] != 0
    b = acc.view(np.uint8).reshape(-1, 4)
    return POP[b].sum(axis=1).astype(np.uint8)[ok]


def read_queries(r1: Path, r2: Path, pairs: int, trim_r2: int, skip: int, k: int):
    """Hash every read's six-frame peptide k-mers once.

    Returns the concatenated hashes plus, for each read, the per-frame slices,
    so the best-frame statistic can be recomputed at every index level without
    re-hashing. A read is "known" if *any* one frame is known, which is why the
    frames cannot be pooled.
    """
    hs: list[np.ndarray] = []
    ds: list[np.ndarray] = []
    names: list[tuple[str, int]] = []
    frames: list[list[tuple[int, int]]] = []   # per read: (start, stop) per frame
    pos = 0
    for _name, a, b in read_pairs(r1, r2, pairs, trim_r2, skip):
        for mate, read in ((1, a), (2, b)):
            names.append((_name, mate))
            per_frame: list[tuple[int, int]] = []
            for f in six_frames(read):
                start = pos
                for seg in f.split("*"):
                    if len(seg) < k:
                        continue
                    codes = AA_CODE[np.frombuffer(seg.encode(), dtype=np.uint8)]
                    h = hash_kmers(codes, k)
                    if h.size:
                        hs.append(h)
                        ds.append(distinct_residues(codes, k))
                        pos += h.size
                per_frame.append((start, pos))
            frames.append(per_frame)
    q = np.concatenate(hs) if hs else np.zeros(0, dtype=np.uint64)
    d = np.concatenate(ds) if ds else np.zeros(0, dtype=np.uint8)
    return q, d, names, frames


def mark_present(uq: np.ndarray, present: np.ndarray, faas: list[Path],
                 k: int) -> int:
    """Stream proteomes, marking which query k-mers they contain. Never stores
    the index; returns the number of (non-unique) database k-mers seen.

    `uq` must be the *unique* sorted query k-mers. Marking duplicates directly
    would silently under-count presence, because searchsorted returns only the
    leftmost of a run of equal values, and duplicate k-mers across reads from the
    same organism are the common case rather than the exception.""" 
    seen = 0
    for codes in proteome_chunks(faas):
        h = hash_kmers(codes, k)
        seen += int(h.size)
        if not h.size:
            continue
        idx = np.searchsorted(uq, h)
        np.clip(idx, 0, uq.size - 1, out=idx)
        hit = uq[idx] == h
        if hit.any():
            present[idx[hit]] = True
    return seen


def best_frame_absence(frames, present_unsorted: np.ndarray,
                       mask: np.ndarray | None = None) -> np.ndarray:
    """Per read: the lowest absent fraction across its six frames.

    `mask` restricts the k-mers counted (used for the complex-k-mers-only
    variant). A frame with no eligible k-mer left contributes nothing, and a
    read with none at all scores 1.0, which is the same convention as a read
    whose every k-mer is absent.
    """
    out = np.empty(len(frames), dtype=np.float64)
    for i, per_frame in enumerate(frames):
        best = np.nan
        for s, e in per_frame:
            if e <= s:
                continue
            pres = present_unsorted[s:e]
            if mask is not None:
                m = mask[s:e]
                if not m.any():
                    continue
                pres = pres[m]
            frac = 1.0 - float(pres.mean())
            best = frac if np.isnan(best) else min(best, frac)
        out[i] = 1.0 if np.isnan(best) else best
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=20000)
    ap.add_argument("--block", type=int, default=400000)
    # No default: a read path is a property of the sequencing run, so it resolves from
    # --library like the trim does. The defaults these two flags used to carry named
    # `data/ww_berlin/`, which no manifest documents (results/PROVENANCE.md §5).
    ap.add_argument("--r1", type=Path, default=None)
    ap.add_argument("--r2", type=Path, default=None)
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--min-distinct", type=int, default=5,
                    help="a k-mer is 'complex' at this many distinct residues")
    ap.add_argument("--out", type=Path, default=Path("results/absence_curve"))
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["r1", "r2", "trim_r2"])

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    np.random.seed(args.seed)

    levels = [
        ("viral", [Path("data/prot/viral.1.protein.faa.gz")]),
        ("+sprot", [Path("data/prot/uniprot_sprot.fasta.gz")]),
        ("+uniref50_1pct", [Path("data/prot/uniref50_every100.fasta.gz")]),
        ("+uniref50_10pct", [Path("data/prot/uniref50_every10.fasta.gz")]),
        # the sample's own proteins, split by whether their contig is ribosomal
        # (experiments/rrna_contigs.py). The split is cumulative, so the first
        # of the two is the cleaned index and the pair together is the whole
        # sample-self index measured before.
        ("+sample_nonrrna", [Path("data/prot/berlin_sample_proteins_nonrrna.faa")]),
        ("+sample_rrna_contigs",
         [Path("data/prot/berlin_sample_proteins_rrnacontigs.faa")]),
    ]

    skip = (args.seed - 1) * args.block
    q, dist, names, frames = read_queries(args.r1, args.r2, args.pairs,
                                          args.trim_r2, skip, args.k)
    complex_kmer = dist >= args.min_distinct
    uq, back = np.unique(q, return_inverse=True)
    present_u = np.zeros(uq.size, dtype=bool)

    rows = []
    per_level_absence: dict[str, np.ndarray] = {}
    per_level_absence_cx: dict[str, np.ndarray] = {}
    cum_kmers = 0
    for name, faas in levels:
        cum_kmers += mark_present(uq, present_u, faas, args.k)
        pres = present_u[back]
        absence = best_frame_absence(frames, pres)
        absence_cx = best_frame_absence(frames, pres, complex_kmer)
        per_level_absence[name] = absence
        per_level_absence_cx[name] = absence_cx
        s1, s2 = absence[0::2], absence[1::2]
        a1, a2 = s1 >= 1.0, s2 >= 1.0
        single = float(np.mean(np.concatenate([a1, a2])))
        rows.append({
            "level": name,
            "db_kmers_cumulative": cum_kmers,
            "fully_absent_single_read": round(single, 4),
            "mean_absence": round(float(absence.mean()), 4),
            "both_mates_fully_absent": round(float(np.mean(a1 & a2)), 4),
            "either_mate_fully_absent": round(float(np.mean(a1 | a2)), 4),
            "excess_over_independence": round(float(np.mean(a1 & a2)) / single ** 2, 3)
            if single > 0 else None,
            "r_mates_continuous": round(float(np.corrcoef(s1, s2)[0, 1]), 4),
            "fully_absent_complex_kmers_only": round(
                float(np.mean(absence_cx >= 1.0)), 4),
            "mean_absence_complex_kmers_only": round(float(absence_cx.mean()), 4),
            "seconds": round(time.perf_counter() - t0, 1),
        })
        print(json.dumps(rows[-1]), flush=True)

    # per-read dump, so P(absent | rRNA) can be computed at every index level by
    # joining the sibling thread's SILVA calls on read_id + mate
    with (args.out / f"reads_s{args.seed}.tsv").open("w") as fh:
        cols = [lv for lv, _ in levels]
        fh.write("read_id\tmate\t" + "\t".join(
            [f"abs_{c}" for c in cols] + [f"abscx_{c}" for c in cols]) + "\n")
        for i, (rid, mate) in enumerate(names):
            fh.write(f"{rid}\t{mate}\t" + "\t".join(
                [f"{per_level_absence[c][i]:.4f}" for c in cols]
                + [f"{per_level_absence_cx[c][i]:.4f}" for c in cols]) + "\n")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    res = {"seed": args.seed, "git_sha": sha, "k": args.k,
           "pairs": len(frames) // 2, "query_kmers": int(q.size),
           "query_kmers_unique": int(uq.size),
           "levels": rows, "seconds": round(time.perf_counter() - t0, 1)}
    (args.out / f"absence_curve_s{args.seed}.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
