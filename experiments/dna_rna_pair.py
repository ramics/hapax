"""The matched DNA and RNA libraries of one wastewater sample, compared rung by rung.

WHY THIS PAIR AND NOT A SECOND CITY. Every number this project has is from one library:
SRR23973193, Berlin, 2021-10-13, RNA. The standing worry has been that none of it
generalises. A second city would answer that but changes everything at once - catchment,
season, protocol, operator. PRJNA948850 turns out to hold 24 `GENOMIC` libraries beside its
215 RNA ones, and 22 of them have an RNA library from the same date and the same collection
window, so the pair used here (`DNA_220616_20-22hmix` / `RNA_WW_220616_20-22hmix`) differs in
ONE variable: which nucleic acid was extracted from the same water, on the same morning, by
the same group. That is a tighter instrument than a second city for the question Imo asked -
whether comparing DNA to RNA costs us the power we need - and it is also the first time any
number here can be checked against a different month of the same plant.

WHAT IS CONFOUNDED AND MUST BE SAID EVERY TIME. The RNA member is 2022-06-16 and our
calibration library is 2021-10-13, eight months apart. So a DNA-RNA difference measured here
is clean, and a comparison of EITHER against our 2021 baseline mixes nucleic acid with eight
months of sewage. The script therefore reports the pair against each other as the result,
and against 2021 only as a separately labelled drift check.

WHAT IT MEASURES, IN THE ORDER THAT THE ANSWERS DEPEND ON EACH OTHER.

1. **Are the libraries even shaped the same?** Read length per mate, and per-position base
   composition. This is not boilerplate: the 2021 library carries ~14 nt of non-template
   prefix on mate 2, which every experiment here trims, and that number is a property of one
   prep rather than of sequencing. Assuming it holds for a library that was never checked is
   how a 93 nt mate gets silently cut to 79. So the prefix is MEASURED per library, by the
   composition anomaly that reveals it, and the measured value is what the later rungs use.
2. **rRNA fraction.** Our single largest subtrahend is ribosomal: 34.5% of the 2021 RNA
   library. A DNA library should be nearly free of it, since rDNA is a few copies of a
   genome rather than most of a transcriptome. If that holds, the DNA library skips the
   cascade's most expensive rung for free, and the comparison says what that is worth.
3. **Identifiable fraction**, Kraken2 at the cascade's own operating confidence (0.1) with
   the same database. This is the thesis quantity - how much of the library is already known
   - measured on two nucleic acids from one sample.
4. **Do they see the same organisms?** The overlap of called taxa, weighted by reads, and
   the taxa each library finds alone. A DNA library that merely re-finds the RNA library's
   organisms adds coverage; one that finds different ones adds content, and only the second
   would change what the production tool should sequence.

Every count is a read PAIR, the occurrence-weighted unit used everywhere in this repo, and
every proportion ships with a Wilson interval because these are single-block samples.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import platform
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from subtract_known import (_mm_index, _silva_to_dna, dir_bytes, kraken2, map_filter,
                            provenance, write_fastq)

ROOT = Path(__file__).resolve().parents[1]

#: The pair. Both are `cer` extraction and `BIOMcombi` library prep, 2022-06-16, the 20-22h
#: composite, from the same plant - the only pair-level difference is DNA against RNA.
LIBRARIES = {
    "dna": {"run": "SRR23973375", "title": "DNA_220616_20-22hmix_cer_BIOMcombi",
            "source": "GENOMIC"},
    "rna": {"run": "SRR23973265", "title": "RNA_WW_220616_20-22hmix", "source":
            "TRANSCRIPTOMIC"},
}


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p, z = k / n, 1.959963985
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def base_composition(seqs: list[str], upto: int) -> list[dict[str, float]]:
    cols = []
    for i in range(upto):
        c = collections.Counter(s[i] for s in seqs if len(s) > i)
        n = sum(c.values()) or 1
        cols.append({b: c.get(b, 0) / n for b in "ACGTN"})
    return cols


def nontemplate_prefix(cols: list[dict[str, float]], window: int = 25,
                       tolerance: float = 0.04) -> dict:
    """How many leading positions have a composition unlike the body of the read.

    The 2021 library's mate 2 was found to carry ~14 nt that do not come from the template,
    and the tell is compositional: the first positions sit far from the per-base frequencies
    the rest of the read settles at. So the body is taken as positions `window`..end, and the
    prefix is the run of leading columns whose maximum per-base deviation from the body
    exceeds `tolerance`. Reported with the deviation profile rather than as a bare number,
    because this decides how many bases every later rung throws away and a wrong value is
    invisible downstream.
    """
    if not cols:
        return {"prefix_len": 0, "deviation": []}
    body = {b: st.mean([c[b] for c in cols[window:]]) for b in "ACGT"} \
        if len(cols) > window else {b: 0.25 for b in "ACGT"}
    dev = [max(abs(c[b] - body[b]) for b in "ACGT") for c in cols]
    n = 0
    for d in dev:
        if d <= tolerance:
            break
        n += 1
    return {"prefix_len": n, "tolerance": tolerance,
            "body_composition": {b: round(v, 4) for b, v in body.items()},
            "deviation": [round(d, 4) for d in dev[:40]]}


def gc(seq: str) -> float:
    n = sum(1 for c in seq if c in "ACGT")
    return (sum(1 for c in seq if c in "GC") / n) if n else 0.0


def profile_library(name: str, r1: Path, r2: Path, args, work: Path) -> dict:
    """Everything that can be measured from the reads alone, before any index is consulted."""
    raw = list(read_pairs(r1, r2, args.pairs, trim_r2=0, skip=args.block * (args.seed - 1)))
    if not raw:
        raise SystemExit(f"no pairs read from {r1}")
    s1 = [p[1] for p in raw]
    s2 = [p[2] for p in raw]
    len1 = collections.Counter(len(s) for s in s1)
    len2 = collections.Counter(len(s) for s in s2)
    mode1, mode2 = len1.most_common(1)[0][0], len2.most_common(1)[0][0]
    comp1 = base_composition(s1, min(mode1, 60))
    comp2 = base_composition(s2, min(mode2, 60))
    pre1, pre2 = nontemplate_prefix(comp1), nontemplate_prefix(comp2)
    # Duplicates on the first 50 bases of both mates: a crude but prep-independent rate, and
    # the only one available without mapping. Compared BETWEEN libraries, not against a
    # published figure, so the crudeness cancels.
    seen = collections.Counter((a[:50], b[:50]) for a, b in zip(s1, s2))
    dup = sum(v - 1 for v in seen.values() if v > 1)
    return {"run": LIBRARIES[name]["run"], "title": LIBRARIES[name]["title"],
            "source": LIBRARIES[name]["source"], "pairs": len(raw),
            "modal_length": {"r1": mode1, "r2": mode2},
            "length_spread": {"r1": dict(len1.most_common(4)),
                              "r2": dict(len2.most_common(4))},
            "nontemplate_prefix": {"r1": pre1, "r2": pre2},
            "gc": {"r1": round(st.mean([gc(s) for s in s1]), 4),
                   "r2": round(st.mean([gc(s) for s in s2]), 4)},
            "gc_sd": {"r1": round(st.pstdev([gc(s) for s in s1]), 4),
                      "r2": round(st.pstdev([gc(s) for s in s2]), 4)},
            "duplicate_frac": round(dup / len(raw), 5),
            "_pairs": raw}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=200_000)
    ap.add_argument("--block", type=int, default=1_500_000,
                    help="pairs skipped per seed, as in subtract_known.py")
    ap.add_argument("--data", type=Path, default=ROOT / "data/berlin_pair")
    ap.add_argument("--silva", type=Path, nargs="*",
                    default=[ROOT / "data/silva/SSU.fa.gz", ROOT / "data/silva/LSU.fa.gz"])
    ap.add_argument("--k2-db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--k2-confidence", type=float, default=0.1,
                    help="the cascade's operating confidence; comparability demands the "
                         "same value, not a re-tuned one")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/"
                                 "libexec/classify"))
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/dnarna"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/dna_rna_pair")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    libs: dict[str, dict] = {}
    for name, meta in LIBRARIES.items():
        r1 = args.data / f"{meta['run']}_1.fastq.gz"
        r2 = args.data / f"{meta['run']}_2.fastq.gz"
        for p in (r1, r2):
            if not p.exists():
                raise SystemExit(f"{p} is missing")
        t = time.perf_counter()
        libs[name] = profile_library(name, r1, r2, args, args.work)
        timings[f"profile_{name}"] = time.perf_counter() - t
        print(f"{name}: {libs[name]['pairs']} pairs, modal "
              f"{libs[name]['modal_length']['r1']}/{libs[name]['modal_length']['r2']} nt, "
              f"non-template prefix r1 {libs[name]['nontemplate_prefix']['r1']['prefix_len']}"
              f" r2 {libs[name]['nontemplate_prefix']['r2']['prefix_len']}, GC "
              f"{libs[name]['gc']['r1']:.3f}, duplicates "
              f"{libs[name]['duplicate_frac']:.4f}", flush=True)

    # ---- rRNA, the cascade's largest subtrahend, with the same index and threshold
    silva_fa = args.work / "silva.fa"
    if not silva_fa.exists():
        t = time.perf_counter()
        _silva_to_dna(args.silva, silva_fa)
        timings["silva_prepare"] = time.perf_counter() - t
    t = time.perf_counter()
    idx = _mm_index(silva_fa, args.work / "silva.mmi", args.threads)
    timings["silva_index"] = time.perf_counter() - t

    for name, lib in libs.items():
        pre = lib["nontemplate_prefix"]["r2"]["prefix_len"]
        seqs = {r[0]: (r[1], r[2][pre:]) for r in lib["_pairs"]}
        t = time.perf_counter()
        hits = map_filter(idx, seqs, args.min_mlen)
        timings[f"rrna_{name}"] = time.perf_counter() - t
        k, n = len(hits), len(seqs)
        lib["rrna"] = {"pairs_hit": k, "frac": k / n,
                       "ci95": wilson(k, n), "min_mlen": args.min_mlen,
                       "r2_trimmed_by": pre,
                       "index_bytes": sum(dir_bytes(p) for p in args.silva)}
        print(f"{name}: rRNA {k / n:.4f} of pairs", flush=True)

    # ---- identifiable fraction, same database and confidence as the cascade's rung
    for name, lib in libs.items():
        pre = lib["nontemplate_prefix"]["r2"]["prefix_len"]
        rows = [(r[0], r[1], r[2][pre:], r[3], r[4][pre:]) for r in lib["_pairs"]]
        fq1 = args.work / f"{name}_1.fq"
        fq2 = args.work / f"{name}_2.fq"
        write_fastq(fq1, rows, 1)
        write_fastq(fq2, rows, 2)
        t = time.perf_counter()
        calls, wall, rss = kraken2(args.k2_db, fq1, fq2, args.classify_bin, args.threads,
                                   args.work / f"{name}.k2", args.k2_confidence)
        timings[f"kraken2_{name}"] = time.perf_counter() - t
        n = len(rows)
        lib["kraken2"] = {"db": str(args.k2_db), "db_bytes": dir_bytes(args.k2_db),
                          "confidence": args.k2_confidence,
                          "classified": len(calls), "frac": len(calls) / n,
                          "ci95": wilson(len(calls), n),
                          "wall_seconds": wall, "peak_rss_kb": rss}
        lib["_taxa"] = collections.Counter(calls.values())
        print(f"{name}: identified {len(calls) / n:.4f} at confidence "
              f"{args.k2_confidence} against {args.k2_db.name}", flush=True)
        fq1.unlink(missing_ok=True)
        fq2.unlink(missing_ok=True)

    # ---- do the two nucleic acids see the same organisms?
    tdna, trna = libs["dna"]["_taxa"], libs["rna"]["_taxa"]
    shared = set(tdna) & set(trna)
    def weighted(c, keys):
        tot = sum(c.values()) or 1
        return sum(c[k] for k in keys) / tot
    overlap = {
        "taxa_dna": len(tdna), "taxa_rna": len(trna), "taxa_shared": len(shared),
        "jaccard": len(shared) / len(set(tdna) | set(trna)) if (tdna or trna) else None,
        "dna_reads_on_shared_taxa": weighted(tdna, shared),
        "rna_reads_on_shared_taxa": weighted(trna, shared),
        "top_dna_only": [t for t, _ in tdna.most_common(400) if t not in trna][:15],
        "top_rna_only": [t for t, _ in trna.most_common(400) if t not in tdna][:15],
        "note": "taxids, not names; a taxid common to both says the organism is seen by "
                "both nucleic acids, and the read-weighted columns say how much of each "
                "library those organisms account for",
    }
    print(f"shared taxa {len(shared)} of {len(set(tdna) | set(trna))}; "
          f"{overlap['dna_reads_on_shared_taxa']:.3f} of DNA reads and "
          f"{overlap['rna_reads_on_shared_taxa']:.3f} of RNA reads sit on them", flush=True)

    out = {"script": "experiments/dna_rna_pair.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else
                              [str(x) for x in v] if isinstance(v, list) else v)
                          for k, v in vars(args).items()},
           "pair": {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
                    for k, v in libs.items()},
           "taxon_overlap": overlap,
           "confound": "the RNA member is 2022-06-16 and the project's calibration library "
                       "is 2021-10-13, so DNA-vs-RNA here is clean but either against 2021 "
                       "mixes nucleic acid with eight months",
           "timings": timings, "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    (args.out / f"dna_rna_pair_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    print(f"written {args.out}/dna_rna_pair_s{args.seed}.json", flush=True)


if __name__ == "__main__":
    main()
