"""What are the reads that still have no known peptide 10-mer, after the sample's
own proteins are in the index?

Why this exists. The index-coverage arm (`results/mate_width/REPORT.md` §4.3)
showed that 86.3% of Berlin background reads being "fully absent" against public
protein databases collapses to 22.2% once the proteins predicted from the
sample's own assembly are indexed. That settles where the 86% came from and
leaves the residual unexplained, and the residual is what any absence-based
detector would actually flag. Two mechanisms predict it, and they have opposite
consequences:

1. **Rare and unsequenced.** Reads from organisms neither in public databases nor
   abundant enough to assemble here. If this dominates, the feature separates
   "abundant or known" from "rare and unknown" - and the detection target sits on
   the same side as the false positives, so no index fixes it.
2. **Abundant but non-coding.** A read off ribosomal RNA, a UTR, or a structural
   RNA can never match a protein index however complete, and in a total-RNA
   library rRNA is the most abundant thing present, not the rarest.

These are distinguishable with data already on disk, which is why this is cheap.
Every read is mapped to the sample's own assembly, and SPAdes writes per-contig
coverage into the contig name, so abundance comes free with the mapping:

- mechanism 1 predicts residual reads are mostly **unmapped**, or mapped to
  **low-coverage** contigs;
- mechanism 2 predicts they are mapped, to **ordinary or high-coverage** contigs,
  **outside** any predicted CDS.

Two artefact classes are measured alongside, because a residue defined by
matching nothing is enriched for them and they are neither mechanism: adapter
content (an exact Illumina adapter prefix) and low sequence complexity (the
dinucleotide-entropy and longest-homopolymer readouts, continuous rather than a
flag, so the comparison is a distribution rather than a threshold).

Readouts are reported for the fully-absent reads and for the rest side by side;
the contrast between the two groups is the measurement, since any absolute rate
here is a property of this library.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts
from insert_profile import read_pairs            # noqa: E402
from peptide_absence import Scorer, build_index, read_fasta  # noqa: E402

ADAPTER = "AGATCGGAAGAGC"   # Illumina TruSeq/Nextera read-through prefix


def contig_cov(name: str) -> float | None:
    """SPAdes writes `..._cov_48.913574` into the contig name: free abundance."""
    for part in name.split("_"):
        pass
    toks = name.split("_")
    for i, t in enumerate(toks):
        if t == "cov" and i + 1 < len(toks):
            try:
                return float(toks[i + 1])
            except ValueError:
                return None
    return None


def dinuc_entropy(s: str) -> float:
    """Shannon entropy over dinucleotides, in bits, max 4.

    Chosen over a mononucleotide entropy because the low-complexity sequence that
    collects in an unmatched residue is usually repeat-structured rather than
    base-biased: `ATATATAT` has maximal mononucleotide entropy.
    """
    if len(s) < 3:
        return 0.0
    c = Counter(s[i:i + 2] for i in range(len(s) - 1))
    n = sum(c.values())
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def longest_homopolymer(s: str) -> int:
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def cds_intervals(ref: Path, cache: Path) -> dict[str, list[tuple[int, int]]]:
    """Predicted CDS intervals per contig, cached.

    The same pyrodigal meta call that produced the indexed proteins, kept as
    coordinates so a mapped read can be placed inside or outside a gene.
    """
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        return {k: [tuple(x) for x in v] for k, v in json.loads(cache.read_text()).items()}
    import pyrodigal
    gf = pyrodigal.GeneFinder(meta=True)
    out: dict[str, list[tuple[int, int]]] = {}
    for name, seq in read_fasta(ref):
        genes = gf.find_genes(seq.encode())
        if len(genes):
            out[name.split()[0]] = [(g.begin, g.end) for g in genes]
    cache.write_text(json.dumps(out))
    return out


def in_cds(iv: list[tuple[int, int]], start: int, end: int) -> bool:
    return any(not (end < b or start > e) for b, e in iv)


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
    ap.add_argument("--ref", type=Path, default=Path("data/ww_berlin/spades_sample.fa"))
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--faa", type=Path, nargs="+",
                    default=[Path("data/prot/viral.1.protein.faa.gz"),
                             Path("data/prot/uniprot_sprot.fasta.gz"),
                             Path("data/prot/berlin_sample_proteins.faa")])
    ap.add_argument("--gpff", type=Path, default=Path("data/prot/viral.1.protein.gpff.gz"))
    ap.add_argument("--cache", type=Path, default=Path("data/pepidx"))
    ap.add_argument("--out", type=Path, default=Path("results/absence_residual"))
    args = ap.parse_args()
    # A trim length is a fact about a sequencing run, not a default. It is resolved
    # from experiments/libraries.json so a library this repo has no entry for is
    # refused by name rather than silently given Berlin's 14 nt, which on a library
    # with no non-template prefix deletes real bases from every mate, in the
    # direction of less recognition.
    library_facts.resolve_args(args, ["r1", "r2", "trim_r2"])

    random.seed(args.seed)
    np.random.seed(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    idx = build_index(args.faa, args.k, args.cache, None, args.gpff)
    scorer = Scorer(idx, args.k)
    iv = cds_intervals(args.ref, args.cache / "cds_intervals_spades_sample.json")

    import mappy
    aln = mappy.Aligner(str(args.ref), preset="sr", n_threads=4)
    if not aln:
        raise SystemExit(f"could not index {args.ref}")

    # Why the mate is carried here. A residual read that does not map to the
    # sample's own assembly is ambiguous on its own: it can come from an
    # organism too rare to assemble (the hypothesis under test) or be an
    # unalignable fragment of an abundant molecule - a junction, an
    # adapter-bearing end, a chimaera. The mate disambiguates, because both
    # mates come off the same molecule: if the mate maps to a high-coverage
    # contig, the molecule was abundant whatever the read itself does.
    groups = {"absent": [], "rest": []}
    skip = (args.seed - 1) * args.block
    n = 0
    for name, r1, r2 in read_pairs(args.r1, args.r2, args.pairs, args.trim_r2, skip):
        hits = {}
        for lbl, read in (("r1", r1), ("r2", r2)):
            hits[lbl] = next(iter(aln.map(read)), None) if len(read) >= 60 else None
        for lbl, read in (("r1", r1), ("r2", r2)):
            if len(read) < 60:
                continue
            mate = hits["r2" if lbl == "r1" else "r1"]
            n += 1
            frac, _, _ = scorer.score_nt_best_frame(read)
            g = "absent" if frac >= 1.0 else "rest"
            hit = hits[lbl]
            rec = {
                "mapped": hit is not None,
                "cov": None, "in_cds": None, "mapq": None,
                "entropy": dinuc_entropy(read),
                "homopolymer": longest_homopolymer(read),
                "adapter": ADAPTER in read,
                "mate_mapped": mate is not None,
                "mate_cov": contig_cov(mate.ctg) if mate is not None else None,
                "score": frac,
            }
            if hit is not None:
                rec["cov"] = contig_cov(hit.ctg)
                rec["mapq"] = hit.mapq
                rec["in_cds"] = in_cds(iv.get(hit.ctg, []), hit.r_st, hit.r_en)
            rec["read_id"] = name
            rec["mate"] = 1 if lbl == "r1" else 2
            rec["contig"] = hit.ctg if hit is not None else ""
            groups[g].append(rec)

    def summarise(rows: list[dict]) -> dict:
        if not rows:
            return {"n": 0}
        mapped = [r for r in rows if r["mapped"]]
        covs = [r["cov"] for r in mapped if r["cov"] is not None]
        cds = [r["in_cds"] for r in mapped if r["in_cds"] is not None]
        return {
            "n": len(rows),
            "mapped_frac": round(len(mapped) / len(rows), 4),
            "cov_median": round(float(np.median(covs)), 2) if covs else None,
            "cov_p10": round(float(np.percentile(covs, 10)), 2) if covs else None,
            "cov_p90": round(float(np.percentile(covs, 90)), 2) if covs else None,
            "cov_lt_2_frac": round(float(np.mean([c < 2 for c in covs])), 4) if covs else None,
            "in_cds_frac_of_mapped": round(float(np.mean(cds)), 4) if cds else None,
            "entropy_mean": round(float(np.mean([r["entropy"] for r in rows])), 4),
            "entropy_p10": round(float(np.percentile([r["entropy"] for r in rows], 10)), 4),
            "homopolymer_mean": round(float(np.mean([r["homopolymer"] for r in rows])), 2),
            "homopolymer_ge10_frac": round(
                float(np.mean([r["homopolymer"] >= 10 for r in rows])), 4),
            "adapter_frac": round(float(np.mean([r["adapter"] for r in rows])), 4),
            "mate_mapped_frac": round(float(np.mean([r["mate_mapped"] for r in rows])), 4),
            "unmapped_but_mate_mapped_frac": round(
                float(np.mean([(not r["mapped"]) and r["mate_mapped"] for r in rows])), 4),
            "mate_cov_median_where_self_unmapped": round(float(np.median(
                [r["mate_cov"] for r in rows
                 if not r["mapped"] and r["mate_cov"] is not None])), 2)
            if any(not r["mapped"] and r["mate_cov"] is not None for r in rows) else None,
        }

    # Per-read dump, so the rRNA arm in the sibling thread can join positive
    # "this read IS rRNA" calls onto these columns: "outside any Prodigal call"
    # is the absence of a gene call, not the presence of non-coding sequence,
    # and a prokaryotic caller over a mixed metatranscriptome will miss spliced
    # eukaryotic genes and atypical viral ORFs straight into the non-coding bin.
    cols = ["read_id", "mate", "absent", "mapped", "contig", "cov", "in_cds",
            "entropy", "homopolymer", "adapter"]
    dump = args.out / f"reads_s{args.seed}.tsv"
    with dump.open("w") as fh:
        fh.write("read_id\tmate\tabsent\tmapped\tcontig\tcontig_cov\tin_cds"
                 "\tentropy\thomopolymer\tadapter\n")
        for g, rows in groups.items():
            for r in rows:
                r = dict(r, absent=(g == "absent"))
                fh.write("\t".join(
                    "" if r[c] is None else
                    (str(int(r[c])) if isinstance(r[c], bool) else str(r[c]))
                    for c in cols) + "\n")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    res = {
        "seed": args.seed, "git_sha": sha, "reads_scored": n,
        "index_faa": [str(f) for f in args.faa], "k": args.k,
        "absent_frac_of_reads": round(len(groups["absent"]) / max(n, 1), 4),
        "absent": summarise(groups["absent"]),
        "rest": summarise(groups["rest"]),
        "seconds": round(time.perf_counter() - t0, 1),
    }
    (args.out / f"absence_residual_s{args.seed}.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "index_faa"}, indent=2))


if __name__ == "__main__":
    main()
