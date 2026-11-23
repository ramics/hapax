"""How much of a real wastewater library can be subtracted as already known, per gigabyte
of index, and what survives.

THE QUESTION. The project's thesis, as Imo put it on 2026-11-20: find the most efficient
way to extract the known, then score the unknown. "Efficient" has two readings and this
script measures both at once - cheapest (bytes of index, seconds of wall time) and most
complete (fraction of the library actually removed). The headline is the curve between
them, and the residue at the far end.

WHY THE POLARITY MATTERS. Every earlier read-level feature in this project failed as a
novelty score for the same reason: a wastewater background is coding too, so "looks like a
gene" cannot separate novel-coding from known-coding
(results/readlevel_features/REPORT.md, PR #8). Subtraction does not have that problem,
because the question asked of the index is "is this read already known", which is what an
index can answer. The failure mode moves somewhere else, and this script is built to
expose it: an incomplete index makes ordinary sewage look novel. PR #11 measured 86% of
Berlin reads with no known peptide 10-mer against a viral + Swiss-Prot index, which is a
statement about the index, not about the sewage. So nothing here reports an absence rate
without the index that produced it, in bytes.

WHAT IS MEASURED. Each tier is run INDEPENDENTLY over the same reads, not cascaded, so
the output supports three different readings:

  - per-tier removal: what this index alone would remove, order-free;
  - union removal: what all tiers together remove, which is the honest "known" fraction;
  - marginal removal in a stated order, which is what a pipeline builder actually buys by
    adding the next tier.

Only the protein tier is run on the survivors of the nucleotide tiers, because DIAMOND on
the full set costs more than the rest of the cascade together; that conditioning is
recorded in the output rather than smoothed over.

THE TIERS, cheapest first.

  artefact   ~0 bytes. Adapter and primer-dimer sequence, low-complexity and
             homopolymer runs, N-rich reads. Load-bearing rather than a footnote: a
             residue defined by matching nothing is enriched for exactly these, and a
             dengue/Nextera study found most or all host-virus chimaeras were artificial
             (Haas 2011 reports amplicon chimaera rates above 70% for the rarest taxa,
             which is our target by construction).
  rRNA       SILVA SSU+LSU NR99. SRR23973193 is library_selection=cDNA with no rRNA
             depletion recorded in ENA, so ribosomal RNA is expected to be a large
             fraction - and it is invisible to a protein index by construction, which is
             the first candidate explanation for the 86% figure.
  k2_*       Kraken2 at four index sizes (viral 0.6 GB, standard-8, PlusPF-8,
             PlusPF-16). Four points on the size axis with the same algorithm, which is
             what makes the removal-per-GB curve a curve rather than a pair of anecdotes.
             Confidence is swept rather than fixed, because the classified fraction of a
             wastewater library is known to depend on it strongly, and a single number
             would be a choice disguised as a measurement.
  selfasm    The study authors' own SPAdes contigs for this sample, and for all 116
             samples. Not a public catalogue at all: it is the same environment's own
             assembled sequence, and it tests whether the residue is unknown to
             databases or merely unassembled here.
  protein    DIAMOND blastx against RefSeq viral proteins + Swiss-Prot, with the bitscore
             threshold swept. Nucleotide identity decays with divergence faster than
             peptide identity, so this tier is the one that can reach a relative the
             nucleotide tiers cannot see.

CALIBRATION, stated up front because every number here is one-sided. Removal measured
against an index is a lower bound on what is known: it cannot distinguish sequence that is
genuinely unrepresented from sequence whose relative exists but is below the threshold.
The script therefore reports, for every tier, the removal rate as a continuous function of
its threshold, never a single yes/no count, and the residue is reported as a function of
bytes subtracted so that the extrapolation to a comprehensive index is visible as an
extrapolation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import resource
import shutil
import statistics as st
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

from berlin_io import read_pairs, revcomp
from library_facts import fact, resolve
from loud import (Consulted, Preflight, cache_is_current, mark_cache, require_output,
                  run_tool)

ROOT = Path(__file__).resolve().parents[1]

# Illumina adapter and primer-dimer seeds. 20 nt prefixes are used as literal seeds: a
# read that begins its adapter anywhere carries the first 20 nt of it unless the insert
# ends mid-adapter, and a 20-mer is long enough that a chance hit in 107 nt of sewage is
# ~1e-7 per read per adapter.
ADAPTERS = {
    "truseq_r1": "AGATCGGAAGAGCACACGTC",
    "truseq_r2": "AGATCGGAAGAGCGTCGTGT",
    "nextera": "CTGTCTCTTATACACATCT",
    "smallrna_3p": "TGGAATTCTCGGGTGCCAAGG",
    "polyA": "A" * 20,
    "polyG": "G" * 20,
}


# ---------------------------------------------------------------------------- artefacts

def shannon_entropy(seq: str, k: int = 3) -> float:
    """Entropy of k-mer usage, in bits per k-mer, normalised to [0, 1].

    DUST-style triplet counting rather than per-base entropy: per-base entropy misses
    tandem repeats, which are the dominant low-complexity class in environmental reads and
    the class most likely to produce a spurious "no known relative" verdict.
    """
    if len(seq) < k:
        return 1.0
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    n = sum(counts.values())
    h = -sum((c / n) * math.log2(c / n) for c in counts.values())
    hmax = math.log2(min(n, 4 ** k))
    return h / hmax if hmax > 0 else 1.0


def longest_homopolymer(seq: str) -> int:
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def artefact_fields(s1: str, s2: str) -> dict:
    """Per-pair artefact evidence. Continuous where possible, so thresholds stay external."""
    out = {}
    hit = None
    for name, seed in ADAPTERS.items():
        if seed in s1 or seed in s2 or revcomp(seed) in s1 or revcomp(seed) in s2:
            hit = name
            break
    out["adapter"] = hit
    # A primer dimer is an adapter-only molecule: the adapter starts at or very near the
    # read start, leaving no insert to speak of.
    dimer = False
    for seed in ADAPTERS.values():
        for s in (s1, s2):
            p = s.find(seed[:15])
            if 0 <= p <= 5:
                dimer = True
    out["primer_dimer"] = dimer
    out["entropy3"] = round(min(shannon_entropy(s1), shannon_entropy(s2)), 4)
    out["homopolymer"] = max(longest_homopolymer(s1), longest_homopolymer(s2))
    out["n_frac"] = round(max(s1.count("N") / max(len(s1), 1),
                              s2.count("N") / max(len(s2), 1)), 4)
    out["gc"] = round((s1 + s2).count("G") / max(len(s1) + len(s2), 1)
                      + (s1 + s2).count("C") / max(len(s1) + len(s2), 1), 4)
    return out


def artefact_flag(f: dict, entropy_min: float, homopolymer_max: int) -> bool:
    return bool(f["primer_dimer"] or f["adapter"] or f["n_frac"] > 0.1
                or f["entropy3"] < entropy_min or f["homopolymer"] > homopolymer_max)


# ------------------------------------------------------------------------------ kraken2

def kraken2(db: Path, fq1: Path, fq2: Path, classify_bin: Path, threads: int,
            out: Path, confidence: float = 0.0,
            memory_map: bool = True) -> tuple[dict[str, int], float, int]:
    """Classify a pair file at one confidence threshold; return {read name: taxid}.

    The threshold is passed to Kraken2 rather than reconstructed afterwards. Kraken2's
    confidence is the fraction of a read's minimizers that fall in the called taxon's
    CLADE, which needs the taxonomy tree to evaluate; a hits-string approximation
    computed here would be a different quantity wearing the same name. The sweep is
    therefore several runs of the classifier, which on a memory-mapped database is cheap
    relative to getting the definition wrong.

    Per-read output rather than Kraken2's own report, because this experiment has to
    intersect tiers read by read.
    """
    t0 = time.perf_counter()
    # -M memory-maps the hash instead of loading it. It was hardcoded, and on an 8 GB
    # database that turns every minimizer lookup into a random disk read: this project
    # measured 100k pairs unfinished in 9 minutes mapped against 27 seconds loaded. The
    # flag changes wall time and peak RSS only - the calls and confidences are identical
    # either way - so it is a parameter rather than a constant, and `timings` and
    # `peak_rss_mb` in the output are only comparable between runs that chose the same
    # value. Default stays -M so previously committed runs reproduce unchanged.
    cmd = [str(classify_bin), "-H", str(db / "hash.k2d"), "-t", str(db / "taxo.k2d"),
           "-o", str(db / "opts.k2d"), *(["-M"] if memory_map else []),
           "-p", str(threads), "-P",
           "-T", str(confidence), "-O", str(out), str(fq1), str(fq2)]
    run_tool(cmd, f"Kraken2 {db.name} at confidence {confidence}",
             stdout=subprocess.DEVNULL)
    require_output(out, f"Kraken2 {db.name}")
    wall = time.perf_counter() - t0
    calls: dict[str, int] = {}
    with open(out) as h:
        for line in h:
            _st, name, taxid, _len, _hits = line.rstrip("\n").split("\t")
            calls[name] = int(taxid)
    return calls, wall, _peak_rss()


# ------------------------------------------------------------------------------- mappy

def map_filter(ref_index, seqs: dict[str, tuple[str, str]], min_mlen: int):
    """Best nucleotide alignment per pair against a mappy index.

    Returns {name: (max matching bases over the two mates, reference hit)}. The matching
    base count is kept rather than a boolean, because it is the continuous readout: a read
    sharing 30 of 107 bases with a reference is a different statement from one sharing
    100, and a threshold applied inside this function would throw that away.
    """
    out = {}
    for name, (s1, s2) in seqs.items():
        best, which = 0, None
        for s in (s1, s2):
            for h in ref_index.map(s):
                if h.mlen > best:
                    best, which = h.mlen, h.ctg
        if best >= min_mlen:
            out[name] = (best, which)
    return out


# ------------------------------------------------------------------------------ helpers

def _peak_rss() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024


def dir_bytes(p: Path) -> int:
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def tree_dirty() -> bool:
    try:
        return subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
    except Exception:
        return True


def provenance(script) -> dict:
    """The three facts a manifest needs to identify the code that produced it.

    `git_sha()` alone is not enough and a bare HEAD SHA is worse than not enough: it looks
    resolved while the tree was dirty. So the flag is explicit rather than hidden in a
    suffix, and the script's own content hash is carried beside it, which is the only one
    of the three that pins the code when the tree was not clean.
    """
    return {"git_sha": git_sha(), "tree_dirty": tree_dirty(),
            "script_sha256": script_sha256(script)}


def script_sha256(path) -> str:
    """Content hash of the script that wrote a result.

    `git_sha()` returns HEAD plus "-dirty", and "-dirty" is where provenance ends: a run
    made halfway through an editing session records a SHA that does not describe the code
    that produced it. That is not hypothetical here - a codon-table fix landed between two
    seeds of the same multi-seed run, leaving one results directory with some seeds right
    and some wrong, both tagged with a dirty SHA. Hashing the file itself pins the code
    whatever the working tree looked like, and makes two runs comparable by inspection.
    """
    import hashlib
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]
    except Exception:
        return "unknown"


def spread(xs) -> dict:
    xs = list(xs)
    return {"mean": st.mean(xs), "min": min(xs), "max": max(xs),
            "sd": st.stdev(xs) if len(xs) > 1 else 0.0, "n": len(xs)}


def write_fastq(path: Path, rows, mate: int) -> None:
    with open(path, "w") as h:
        for name, s1, s2, q1, q2 in rows:
            s, q = (s1, q1) if mate == 1 else (s2, q2)
            h.write(f"@{name}\n{s}\n+\n{q}\n")


def write_fasta(path: Path, rows) -> None:
    with open(path, "w") as h:
        for name, s1, s2, _q1, _q2 in rows:
            h.write(f">{name}/1\n{s1}\n>{name}/2\n{s2}\n")


# --------------------------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=200_000)
    # Three values that are properties of a sequencing run and were written here as
    # Berlin's. --block sized to 8.73 M pairs runs off the end of a smaller library, and
    # --trim-r2 14 deletes 14 real bases from every mate of a library with no
    # non-template prefix, silently, in the direction of less recognition. They now come
    # from experiments/libraries.json by accession, and there is no default: an unknown
    # library raises rather than inheriting Berlin's biology.
    ap.add_argument("--library", default="SRR23973193",
                    help="accession in experiments/libraries.json, supplying --block, "
                         "--trim-r2 and the read paths unless those are given")
    ap.add_argument("--block", type=int, default=None,
                    help="pairs skipped per seed: distinct regions of the flow cell. "
                         "From the library unless given")
    ap.add_argument("--r1", type=Path, default=None)
    ap.add_argument("--r2", type=Path, default=None)
    ap.add_argument("--trim-r2", type=int, default=None)
    ap.add_argument("--k2-dbs", type=Path, nargs="+",
                    default=[ROOT / "data/k2_viral", ROOT / "data/k2_standard08",
                             ROOT / "data/k2_pluspf08", ROOT / "data/k2_pluspf16"])
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--silva", type=Path, nargs="*",
                    default=[ROOT / "data/silva/SSU.fa.gz", ROOT / "data/silva/LSU.fa.gz"])
    ap.add_argument("--selfasm", type=Path, nargs="*",
                    default=[ROOT / "data/geo/sample_contigs.fa"])
    ap.add_argument("--observe", type=Path, nargs="*", default=[],
                    help="references mapped and reported but NOT subtracted: the "
                         "sample's own contigs that nothing identified, which are tier "
                         "2's population rather than background")
    ap.add_argument("--dmnd", type=Path, nargs="*",
                    default=[ROOT / "data/prot/viral.dmnd", ROOT / "data/prot/sprot.dmnd"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--k2-confidence", type=float, nargs="+",
                    default=[0.0, 0.05, 0.1, 0.2],
                    help="Kraken2 confidence thresholds to sweep; each is a separate run")
    ap.add_argument("--no-k2-memory-map", dest="k2_memory_map", action="store_false",
                    help="load each Kraken2 hash into RAM instead of memory-mapping it. "
                         "The calls and confidences are identical; only wall time and "
                         "peak RSS change, so `timings` and `peak_rss_mb` are comparable "
                         "only across runs that made the same choice. On an 8 GB database "
                         "memory-mapping makes every lookup a random disk read - measured "
                         "here at 100k pairs unfinished in 9 minutes mapped against 27 "
                         "seconds loaded - so this is the setting to use when the machine "
                         "has the memory and the database is large")
    ap.add_argument("--k2-operating-confidence", type=float, default=0.1,
                    help="the threshold whose removals enter the cascade. 0.0 accepts a "
                         "single matching minimizer, which on a wastewater library is "
                         "mostly noise; 0.1 is the lowest value in the sweep at which the "
                         "classified fraction stops being dominated by one-minimizer calls")
    ap.add_argument("--diamond-sensitivity", default="--sensitive",
                    help="DIAMOND mode for the subtraction tier. --very-sensitive is kept "
                         "for the residue report, where the query set is 10x smaller")
    ap.add_argument("--diamond-bitscore", type=float, default=50.0,
                    help="bitscore at which a protein hit counts as subtracted; swept in output")
    ap.add_argument("--entropy-min", type=float, default=0.5)
    ap.add_argument("--homopolymer-max", type=int, default=25)
    ap.add_argument("--min-mlen", type=int, default=40,
                    help="matching bases needed to call a nucleotide hit; swept in output")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/subtract"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/subtract_known")
    ap.add_argument("--keep-known", type=int, default=20_000,
                    help="subtracted pairs kept as the known-positive comparison set")
    ap.add_argument("--rows-out", action="store_true",
                    help="write the per-pair evidence table for the residue")
    ap.add_argument("--preflight", action="store_true",
                    help="check every input, index and tool, then exit: 0 if the run "
                         "could start, 2 if not. hapax.py calls this before any stage runs")
    args = ap.parse_args()

    # Resolve the per-library facts before the preflight, so an unknown library is
    # refused by name rather than reported as four missing paths.
    args.block = resolve(args.library, "block", args.block)
    args.trim_r2 = resolve(args.library, "trim_r2", args.trim_r2)
    if args.r1 is None or args.r2 is None:
        if args.r1 is not None or args.r2 is not None:
            raise SystemExit("--r1 and --r2 are given together or not at all: one of "
                             "each from two libraries is a pairing that does not exist")
        args.r1 = ROOT / fact(args.library, "r1")
        args.r2 = ROOT / fact(args.library, "r2")


    # Every index named on the command line is one the result depends on, so an absent,
    # empty or half-built one is a refusal, not a tier quietly dropped from the curve. The
    # old behaviour was `continue`: the curve lost a row, the manifest still listed the
    # index, and nothing said so. To run without a tier, do not name it.
    pf = Preflight("subtract_known.py")
    pf.file(args.r1, "--r1")
    pf.file(args.r2, "--r2")
    for db in args.k2_dbs:
        pf.kraken2_db(db)
    if args.k2_dbs:
        pf.exe(args.classify_bin, "--classify-bin (Kraken2 classify)")
        pf.require(args.k2_operating_confidence in args.k2_confidence,
                   f"--k2-operating-confidence {args.k2_operating_confidence} is not in "
                   f"the --k2-confidence sweep {args.k2_confidence}, so every Kraken2 tier "
                   f"would remove nothing and report it as a measurement")
    for p in args.silva:
        pf.file(p, "--silva")
    for p in args.selfasm:
        pf.file(p, "--selfasm")
    for p in args.observe:
        pf.file(p, "--observe")
    for p in args.dmnd:
        pf.file(p, "--dmnd")
    if args.dmnd:
        pf.exe(args.diamond_bin, "--diamond-bin")
    stems = [p.stem for p in [*args.selfasm, *args.observe]]
    pf.require(len(stems) == len(set(stems)),
               f"--selfasm/--observe references share a file stem ({stems}); tiers and "
               f"cached indexes are keyed by stem, so one would overwrite the other")
    names = [d.name for d in args.k2_dbs]
    pf.require(len(names) == len(set(names)),
               f"--k2-dbs share a directory name ({names}); tiers are keyed by it")
    pf.module("mappy", "nucleotide tiers")
    pf.done(args.preflight)
    consulted = Consulted()

    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()
    timings: dict[str, float] = {}
    index_bytes: dict[str, int] = {}

    # ---- read a block
    t0 = time.perf_counter()
    rows = list(read_pairs(args.r1, args.r2, args.pairs, args.trim_r2,
                           (args.seed - 1) * args.block))
    timings["read_block"] = time.perf_counter() - t0
    names = [r[0] for r in rows]
    seqs = {r[0]: (r[1], r[2]) for r in rows}
    n = len(rows)
    if n == 0:
        raise SystemExit("no pairs read")
    if n < args.pairs:
        # Not fatal - a short library is a real library - but it changes every
        # denominator, so it is said out loud and recorded in the result.
        print(f"WARNING: asked for {args.pairs} pairs and the block held {n}", flush=True)

    fq1, fq2 = args.work / f"s{args.seed}_1.fq", args.work / f"s{args.seed}_2.fq"
    write_fastq(fq1, rows, 1)
    write_fastq(fq2, rows, 2)

    removed: dict[str, set[str]] = {}
    tier_meta: dict[str, dict] = {}

    # ---- tier: artefact
    t0 = time.perf_counter()
    art = {r[0]: artefact_fields(r[1], r[2]) for r in rows}
    removed["artefact"] = {k for k, f in art.items()
                           if artefact_flag(f, args.entropy_min, args.homopolymer_max)}
    timings["artefact"] = time.perf_counter() - t0
    index_bytes["artefact"] = sum(len(s) for s in ADAPTERS.values())
    tier_meta["artefact"] = {
        "adapter_any": sum(1 for f in art.values() if f["adapter"]) / n,
        "primer_dimer": sum(1 for f in art.values() if f["primer_dimer"]) / n,
        "low_entropy": sum(1 for f in art.values() if f["entropy3"] < args.entropy_min) / n,
        "homopolymer": sum(1 for f in art.values() if f["homopolymer"] > args.homopolymer_max) / n,
        "n_rich": sum(1 for f in art.values() if f["n_frac"] > 0.1) / n,
        "entropy_sweep": {f"{thr:.2f}": sum(1 for f in art.values() if f["entropy3"] < thr) / n
                          for thr in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)},
    }

    # ---- tier: rRNA (SILVA), nucleotide
    import mappy

    rrna_hits = {}
    if args.silva:
        silva_fa = args.work / "silva.fa"
        # Keyed to its sources: a work directory reused with a different --silva would
        # otherwise serve the old reference under the new one's name.
        if not cache_is_current(silva_fa, args.silva):
            t0 = time.perf_counter()
            _silva_to_dna(args.silva, silva_fa)
            mark_cache(silva_fa, args.silva)
            timings["silva_prepare"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        idx = _mm_index(silva_fa, args.work / "silva.mmi", args.threads)
        timings["silva_index"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        rrna_hits = map_filter(idx, seqs, args.min_mlen)
        timings["rrna"] = time.perf_counter() - t0
        index_bytes["rrna"] = sum(dir_bytes(p) for p in args.silva)
        removed["rrna"] = set(rrna_hits)
        for p in args.silva:
            consulted.searched("rrna", p, queries=n, hits=len(rrna_hits))
        tier_meta["rrna"] = {"mlen_sweep": _mlen_sweep(rrna_hits, n)}
        del idx

    # ---- tiers: kraken2 at several index sizes
    for db in args.k2_dbs:
        tier = db.name
        out = args.work / f"{tier}_s{args.seed}.k2"
        sweep, taxa, wall_total = {}, None, 0.0
        for conf in args.k2_confidence:
            calls, wall, rss = kraken2(db, fq1, fq2, args.classify_bin, args.threads,
                                       out, conf,
                                       memory_map=args.k2_memory_map)
            if len(calls) != n:
                raise SystemExit(f"Kraken2 {tier} at {conf} returned {len(calls)} calls "
                                 f"for {n} pairs; Kraken2 writes one line per pair")
            hit = {k for k, tx in calls.items() if tx != 0}
            sweep[f"{conf:.2f}"] = len(hit) / n
            wall_total += wall
            if conf == args.k2_operating_confidence:
                removed[tier] = hit
                taxa = Counter(tx for tx in calls.values() if tx != 0).most_common(15)
            out.unlink(missing_ok=True)
        timings[tier] = wall_total
        index_bytes[tier] = dir_bytes(db / "hash.k2d") + dir_bytes(db / "taxo.k2d")
        if tier not in removed:
            raise SystemExit(f"internal: Kraken2 {tier} never ran at the operating "
                             f"confidence {args.k2_operating_confidence}")
        consulted.searched(tier, db, queries=n, hits=len(removed[tier]),
                           confidence=args.k2_operating_confidence)
        tier_meta[tier] = {
            "peak_rss_bytes": rss,
            "wall_seconds_per_run": wall_total / max(len(args.k2_confidence), 1),
            "operating_confidence": args.k2_operating_confidence,
            "confidence_sweep": sweep,
            "top_taxa": taxa,
        }

    # ---- tier: the environment's own assembled sequence
    for asm in args.selfasm:
        tier = f"selfasm_{asm.stem}"
        t0 = time.perf_counter()
        idx = _mm_index(asm, args.work / f"{asm.stem}.mmi", args.threads)
        timings[f"{tier}_index"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        hits = map_filter(idx, seqs, args.min_mlen)
        timings[tier] = time.perf_counter() - t0
        index_bytes[tier] = dir_bytes(asm)
        removed[tier] = set(hits)
        consulted.searched(tier, asm, queries=n, hits=len(hits))
        tier_meta[tier] = {"mlen_sweep": _mlen_sweep(hits, n)}
        del idx

    # ---- observed but NOT subtracted: the sample's own contigs that nothing identified
    #
    # Imo's three-tier architecture (2026-11-20): tier 1 removes what we can identify,
    # tier 2 is what assembled but could not be identified, tier 3 is the residual reads.
    # Subtracting the assembly wholesale - what the first cut of this script did - deletes
    # tier 2 before anyone looks at it, because a novel virus abundant enough to assemble
    # would be removed as background. So these references are mapped and reported, and
    # deliberately left out of the subtracted union: reads that are explained ONLY by an
    # unidentifiable contig are tier 2's read-level population, not residue.
    observed: dict[str, set[str]] = {}
    for ref in args.observe:
        tier = f"observe_{ref.stem}"
        t0 = time.perf_counter()
        idx = _mm_index(ref, args.work / f"{ref.stem}.mmi", args.threads)
        hits = map_filter(idx, seqs, args.min_mlen)
        timings[tier] = time.perf_counter() - t0
        observed[tier] = set(hits)
        consulted.searched(tier, ref, queries=n, hits=len(hits), subtracted=False)
        index_bytes[tier] = dir_bytes(ref)
        tier_meta[tier] = {"subtracted": False, "mlen_sweep": _mlen_sweep(hits, n),
                           "tier_alone": len(hits) / n}
        del idx

    # ---- tier: protein, on the survivors of everything above
    nucleotide_union = set().union(*removed.values()) if removed else set()
    survivors = [r for r in rows if r[0] not in nucleotide_union]
    prot_rows = {}
    if args.dmnd and survivors:
        fa = args.work / f"residue_s{args.seed}.fa"
        write_fasta(fa, survivors)
        for db in args.dmnd:
            tier = f"prot_{db.stem}"
            t0 = time.perf_counter()
            tsv = args.work / f"{tier}_s{args.seed}.tsv"
            tsv.unlink(missing_ok=True)
            run_tool([str(args.diamond_bin), "blastx", "--db", str(db), "--query", str(fa),
                      "--out", str(tsv), "--threads", str(args.threads),
                      "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                      "evalue", "bitscore", "--max-target-seqs", "1",
                      "--evalue", "10", args.diamond_sensitivity, "--quiet"],
                     f"DIAMOND blastx against {db.name}")
            # An empty table is a real answer (searched, zero hits); a missing one is not.
            require_output(tsv, f"DIAMOND {db.name}", allow_empty=True)
            timings[tier] = time.perf_counter() - t0
            index_bytes[tier] = dir_bytes(db)
            best = {}
            with open(tsv) as h:
                for line in h:
                    q, s, pid, length, ev, bits = line.split("\t")
                    name = q.rsplit("/", 1)[0]
                    bits = float(bits)
                    if bits > best.get(name, (0.0,))[0]:
                        best[name] = (bits, s, float(pid), int(length), float(ev))
            prot_rows[tier] = best
            removed[tier] = {k for k, v in best.items() if v[0] >= args.diamond_bitscore}
            consulted.searched(tier, db, queries=len(survivors), hits=len(best),
                               removed=len(removed[tier]),
                               conditioned_on="survivors of all nucleotide tiers")
            tier_meta[tier] = {
                "conditioned_on": "survivors of all nucleotide tiers",
                "survivors_in": len(survivors),
                "bitscore_sweep": {str(b): sum(1 for v in best.values() if v[0] >= b) / len(survivors)
                                   for b in (30, 40, 50, 60, 80, 100)},
                "bitscore_sweep_of_library": {str(b): sum(1 for v in best.values() if v[0] >= b) / n
                                              for b in (30, 40, 50, 60, 80, 100)},
                "removed_at_what_distance": _removal_distance(best, args.diamond_bitscore),
            }
            tsv.unlink(missing_ok=True)

    elif args.dmnd:
        for db in args.dmnd:
            consulted.not_searched(f"prot_{db.stem}", db,
                                   "no pair survived the nucleotide tiers, so there was "
                                   "nothing to search; this is not a zero-hit search")

    # ---- cumulative curve in cheapest-first order
    #
    # rRNA FIRST, with a NUCLEOTIDE index, is a design requirement and not an artefact of
    # sorting by price. Against a large protein index ribosomal reads stop looking novel
    # faster than anything else, because public databases carry automatically annotated
    # entries that are really translated ribosomal sequence (the conservation multiplier
    # measured database-independently by the protein thread; the wastewater entry counts
    # that confirm it, 63 entries each reaching >=1% of ribosomal reads, from PR #10). In
    # a total-RNA library that is a third of the reads, so a cascade that consults a
    # protein index first measures public annotation quality rather than the sample.
    order = [t for t in ["artefact", "rrna"] if t in removed]
    assert not any(t.startswith("prot_") for t in order), "rRNA must precede the proteins"
    order += [d.name for d in args.k2_dbs if d.name in removed]
    order += [t for t in removed if t.startswith("selfasm_")]
    order += [t for t in removed if t.startswith("prot_")]
    seen: set[str] = set()
    curve = []
    cum_bytes = 0
    for tier in order:
        before = len(seen)
        seen |= removed[tier]
        cum_bytes += index_bytes.get(tier, 0)
        curve.append({
            "tier": tier,
            "index_bytes": index_bytes.get(tier, 0),
            "cumulative_index_bytes": cum_bytes,
            "tier_alone_removed": len(removed[tier]) / n,
            "marginal_removed": (len(seen) - before) / n,
            "cumulative_removed": len(seen) / n,
            "residue": 1 - len(seen) / n,
            "wall_seconds": timings.get(tier),
        })

    result = {
        "seed": args.seed,
        "pairs": n,
        "block_skipped": (args.seed - 1) * args.block,
        "curve": curve,
        "tiers": tier_meta,
        "union_removed": len(seen) / n,
        "residue_pairs": n - len(seen),
        # The residue, split the way Imo's architecture asks for rather than reported as
        # one number. "tier2_reads" are residual pairs that DO map to one of this
        # sample's own contigs that nothing could identify: present and often abundant
        # here, but unplaceable - the read-level population of tier 2. The rest is the
        # tier-3 population, which is by construction neither in the catalogues nor
        # assembled here, i.e. rare AND unsequenced.
        "residue_composition": {
            tier: len({x for x in hits if x not in seen}) / n
            for tier, hits in observed.items()
        },
        "tier3_pairs": (n - len(seen) - len({x for s in observed.values() for x in s
                                             if x not in seen})) if observed else None,
        "timings_seconds": timings,
        "index_bytes": index_bytes,
    }
    (args.out / f"subtract_s{args.seed}.json").write_text(json.dumps(result, indent=2))

    # The residue and a matched sample of SUBTRACTED pairs are both kept, because every
    # field in the residue report is calibrated against the subtracted reads: a field that
    # reads the same on both carries no information about the residue.
    residue_rows = [r for r in rows if r[0] not in seen]
    if len(residue_rows) != result["residue_pairs"]:
        raise SystemExit(f"internal: {len(residue_rows)} residue rows against "
                         f"residue_pairs {result['residue_pairs']}")
    write_fasta(args.work / f"residue_s{args.seed}.fa", residue_rows)
    import random as _random
    known_pairs = [r for r in rows if r[0] in seen]
    _random.Random(args.seed).shuffle(known_pairs)
    write_fasta(args.work / f"subtracted_s{args.seed}.fa", known_pairs[:args.keep_known])

    if args.rows_out:
        resid = args.out / f"residue_rows_s{args.seed}.jsonl"
        written = 0
        with open(resid, "w") as h:
            for name, s1, s2, q1, q2 in rows:
                if name in seen:
                    continue
                # Base quality travels with the row because a read can be unique to this
                # library for two quite different reasons - a rare organism, or errors -
                # and the abundance field alone cannot tell them apart. Mean error
                # probability rather than mean Phred: averaging Phred scores averages
                # logarithms and understates the error rate.
                errs = [10 ** (-(ord(c) - 33) / 10) for c in (q1 + q2)]
                row = {"name": name, "len1": len(s1), "len2": len(s2), **art[name],
                       "mean_error_prob": round(sum(errs) / len(errs), 6) if errs else None,
                       "frac_bases_under_q20": round(
                           sum(1 for c in (q1 + q2) if ord(c) - 33 < 20) / len(q1 + q2), 4)
                       if (q1 + q2) else None}
                for tier, best in prot_rows.items():
                    if name in best:
                        bits, sid, pid, length, ev = best[name]
                        row[tier] = {"subject": sid, "bitscore": bits, "pident": pid,
                                     "aln_len": length, "evalue": ev}
                if name in rrna_hits:
                    row["rrna_mlen"] = rrna_hits[name][0]
                h.write(json.dumps(row) + "\n")
                written += 1
        if written != result["residue_pairs"]:
            raise SystemExit(f"internal: wrote {written} residue rows against "
                             f"residue_pairs {result['residue_pairs']}")

    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/subtract_known.py",
        **provenance(__file__),
        "pairs_requested": args.pairs,
        "pairs_read": n,
        # What was opened, with its size at the time, beside `parameters` (what was asked
        # for). hapax.py refuses a stage whose two lists disagree.
        "consulted": consulted.as_list(),
        "parameters": {k: (str(v) if isinstance(v, (Path,)) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "peak_rss_bytes": _peak_rss(),
        "python": platform.python_version(),
        "mappy": mappy.__version__,
        "diamond": _version([str(args.diamond_bin), "version"]),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
    }, indent=2))

    for c in curve:
        print(f"{c['tier']:>18}  {c['index_bytes'] / 1e9:7.3f} GB  alone "
              f"{c['tier_alone_removed'] * 100:6.2f}%  marginal {c['marginal_removed'] * 100:6.2f}%  "
              f"cumulative {c['cumulative_removed'] * 100:6.2f}%  residue {c['residue'] * 100:6.2f}%")
    print(f"union removed {len(seen) / n * 100:.2f}%, residue {n - len(seen)} of {n} pairs")

    for f in (fq1, fq2):
        f.unlink(missing_ok=True)


def _silva_to_dna(paths, out: Path) -> None:
    """SILVA ships RNA with U, and mappy wants DNA. Also drops the long taxonomy string
    from the header, which otherwise appears in every hit record."""
    import gzip
    with open(out, "w") as h:
        for p in paths:
            op = gzip.open if str(p).endswith(".gz") else open
            with op(p, "rt") as g:
                for line in g:
                    if line.startswith(">"):
                        h.write(">" + line[1:].split()[0] + "\n")
                    else:
                        h.write(line.strip().replace("U", "T").upper() + "\n")


def _mm_index(fa: Path, mmi: Path, threads: int):
    """Build a minimap2 index once and reuse it across seeds.

    The index is a function of the reference alone, so rebuilding it per seed would put
    tens of seconds of constant cost into every tier's wall time and make the
    cost-per-tier comparison - the point of the experiment - read wrong.

    The cache is keyed to the reference's path, size and mtime. It used to be keyed to the
    file STEM alone, so a work directory reused with a different reference of the same name
    served the old index under the new name - the index named in the provenance not being
    the index used.
    """
    import mappy
    if cache_is_current(mmi, [fa]):
        a = mappy.Aligner(str(mmi), preset="sr", n_threads=threads)
        if not a:
            raise SystemExit(f"could not load cached index {mmi}")
        return a
    mmi.unlink(missing_ok=True)
    a = mappy.Aligner(str(fa), preset="sr", n_threads=threads, fn_idx_out=str(mmi))
    if not a:
        raise SystemExit(f"could not index {fa}")
    mark_cache(mmi, [fa])
    return a


def _removal_distance(best: dict, bitscore: float) -> dict:
    """The distance at which the protein tier removed each read, not just that it did.

    The tier's output is one bit per read, and the holdout thread showed what that costs:
    a contig 30.3% divergent from its withheld reference is removed at bitscore 840
    against a threshold of 50, so the removal is not marginal and no threshold recovers
    it. The identity was computed and thrown away. Keeping its distribution is free -
    `best` already carries pident, alignment length and E-value per read.

    Bias, stated rather than corrected: per-read identity leans towards "familiar". On
    the holdout's known distance, 145 nt reads report about 83% where the genome-level
    truth is 69.72%, up to 96.9% on the shortest alignments, and the error falls from
    +17.2 points below 50% coverage to +7.1 above 90% and is zero at zero distance. So
    every number here is an UPPER bound on identity and a lower bound on divergence, and
    the direction never reverses. Identity is also quantised at short alignment lengths,
    so the histogram runs only over alignments of at least 30 aa and the shorter ones are
    counted separately rather than binned.
    """
    kept = [v for v in best.values() if v[0] >= bitscore]
    scorable = [v for v in kept if v[3] >= 30]
    edges = [0, 30, 40, 50, 60, 70, 80, 90, 95, 100]
    hist = []
    for lo, hi in zip(edges, edges[1:]):
        hist.append({"from_percent": lo, "to_percent": hi,
                     "reads": sum(1 for v in scorable if lo <= v[2] < hi)})
    hist.append({"from_percent": 100, "to_percent": 100.1,
                 "reads": sum(1 for v in scorable if v[2] >= 100)})
    ids = sorted(v[2] for v in scorable)
    aln = sorted(v[3] for v in scorable)
    med = lambda z: z[len(z) // 2] if z else None
    return {
        "unit": "read pairs removed by this tier at the cascade's bitscore",
        "bitscore_threshold": bitscore,
        "removed": len(kept),
        "alignment_at_least_30_aa": len(scorable),
        "alignment_under_30_aa_unscorable": len(kept) - len(scorable),
        "identity_histogram": hist,
        "median_identity_percent": med(ids),
        "median_alignment_aa": med(aln),
        "barely_recognisable_share": (
            round(sum(1 for v in scorable if v[2] < 70) / len(scorable), 4)
            if scorable else None),
        "barely_recognisable_means": ("under 70% identity over at least 30 aa: removed as "
                                      "known, but at a distance where the subject lends "
                                      "its name rather than establishing an identity"),
        "bias": ("per-read identity reads high towards familiar, by +17.2 points below "
                 "50% alignment coverage and +7.1 above 90% on the holdout's known "
                 "distance, so these are upper bounds on identity"),
    }


def _mlen_sweep(hits: dict, n: int) -> dict:
    return {str(m): sum(1 for v in hits.values() if v[0] >= m) / n
            for m in (30, 40, 50, 60, 80, 100)}


def _version(cmd) -> str:
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT).strip().splitlines()[0]
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
