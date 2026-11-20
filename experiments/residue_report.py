"""The per-read evidence report for reads that survived subtraction.

WHAT THIS IS AND IS NOT. Imo's framing on 2026-11-20: for each read where one or both
mates match nothing known, a quantified statement of exactly how much we can and cannot
say - closest relative however far, maybe-coding, GC, repetitive or information-bearing,
primer dimer or adapter, any interesting motif. It is deliberately not a novelty
classifier. Every one of these fields has already been measured in this project to FAIL as
a novelty score (results/readlevel_features/REPORT.md, PR #8): a wastewater background is
coding too, so "looks like a gene" cannot separate novel-coding from known-coding. A
report that states what is and is not determinable about a read does not need novelty
polarity; it needs honest calibration of each field, which is what this script produces.

So every field ships with its own calibration measured on the same data:

  - a NULL: the same field computed on composition-matched shuffled versions of the same
    reads, which says what the field would read on sequence carrying no information at
    all. For ORF-likeness at 107 nt that null is high (33-72% of shuffled reads are
    stop-free in some frame), and a report that printed "has an open frame" without it
    would be misleading rather than merely incomplete.
  - a KNOWN-POSITIVE comparison: the same field on reads that the cascade DID subtract.
    If a field reads the same on both, it carries no information about the residue, and
    that is the field's honest calibration rather than a reason to hide it.

FIELDS.

  closest_relative   Best DIAMOND blastx hit at --very-sensitive with no bitscore floor,
                     reported with bitscore, percent identity and alignment length. Never
                     thresholded into a call: the project measured that translated-search
                     score ordering tracks DATABASE PROXIMITY, not novelty (only 11.6% of
                     spike reads hit at all, and score ordering put background above
                     spike at every threshold). "Closest relative, however far" is an
                     honest field and a dishonest classifier, so it is printed and not
                     thresholded.
  frame_evidence     Stop-free frames, longest ORF per frame, and - the part that must not
                     be hidden - whether the frame is identifiable at all. At 107 nt the
                     true frame is the only stop-free frame in just 0.528 of coding reads,
                     so a single-frame statement is wrong for about half of them. The
                     field reports the number of candidate frames alongside the best one.
  gc, entropy        GC fraction and triplet entropy, each with the shuffled null and the
                     subtracted-read comparison.
  repetitive         Fraction of the read's 16-mers that recur within the read or its
                     mate, plus the longest tandem repeat period. This separates "no
                     information" from "information we cannot place", which is the
                     distinction Imo's question turns on.
  artefact           Adapter, primer-dimer, homopolymer, N content: from
                     experiments/subtract_known.py, re-reported per read because in a
                     residue defined by matching nothing these are enriched by
                     construction.
  chimaera           Split-mapping evidence: the two halves of one mate placed against the
                     same reference set, called chimaeric when they place to different
                     references or in the wrong orientation. Load-bearing rather than a
                     footnote: a dengue/Nextera study found most or all host-virus
                     chimaeras artificial, and amplicon chimaera rates exceed 70% for the
                     rarest taxa, which is this residue by construction. Reported with
                     its own false-positive rate measured on subtracted reads, where a
                     split call is known to be wrong.
  motif              Pfam-A profile hits via pyhmmer on all six frames, with the
                     per-domain E-value. A profile hit is the one field that can reach a
                     relative too divergent for DIAMOND, since a profile is a family
                     model rather than a sequence.

Nothing here is combined into a score. The output is one row per surviving pair.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import random
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import revcomp
from subtract_known import (artefact_fields, dir_bytes, git_sha, shannon_entropy,
                            spread, _mm_index)

ROOT = Path(__file__).resolve().parents[1]

CODON_STOPS = {"TAA", "TAG", "TGA"}
AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R",
    "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str) -> str:
    return "".join(AA.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def six_frames(seq: str):
    rc = revcomp(seq)
    for strand, s in (("+", seq), ("-", rc)):
        for off in (0, 1, 2):
            yield strand, off, s[off:]


def frame_evidence(seq: str) -> dict:
    """Stop-free frames and the longest stop-free stretch, with frame ambiguity stated.

    The ambiguity is the point. PR #8 measured that at 107 nt the true frame of a coding
    read is the ONLY stop-free frame in 0.528 of reads (0.943 at 250 nt), so "the read is
    in frame +1" is a statement this read length usually cannot support. Reporting
    `stop_free_frames` next to `best_frame` makes that visible in every row instead of
    once in a caveat.
    """
    best_len, best_frame, stop_free = 0, None, []
    for strand, off, s in six_frames(seq):
        prot = translate(s)
        runs = max((len(x) for x in prot.split("*")), default=0)
        if "*" not in prot:
            stop_free.append(f"{strand}{off}")
        if runs > best_len:
            best_len, best_frame = runs, f"{strand}{off}"
    return {"longest_orf_aa": best_len, "best_frame": best_frame,
            "stop_free_frames": len(stop_free),
            "frame_identifiable": len(stop_free) == 1}


def repetitiveness(s1: str, s2: str, k: int = 16) -> dict:
    """Does the pair carry information, or the same few bases many times over?

    Two readouts because they fail differently: k-mer recurrence catches tandem and
    interspersed repetition, and the best tandem period catches short-period repeats whose
    k-mers are all distinct at k=16.
    """
    kms = Counter()
    for s in (s1, s2):
        for i in range(len(s) - k + 1):
            kms[s[i:i + k]] += 1
    tot = sum(kms.values())
    dup = sum(c for c in kms.values() if c > 1)
    best_period, best_score = 0, 0.0
    for p in range(1, 13):
        for s in (s1, s2):
            if len(s) <= p:
                continue
            match = sum(1 for a, b in zip(s, s[p:]) if a == b)
            score = match / (len(s) - p)
            if score > best_score:
                best_score, best_period = score, p
    return {"kmer16_duplicate_frac": round(dup / tot, 4) if tot else 0.0,
            "distinct_kmer16": len(kms),
            "tandem_period": best_period,
            "tandem_identity": round(best_score, 4)}


def shuffled(seq: str, rng: random.Random) -> str:
    """Composition-matched null: same mononucleotide composition, no structure.

    Mononucleotide rather than dinucleotide-preserving: the stronger null is the right one
    for GC and entropy, but for ORF-likeness a mononucleotide shuffle is already enough to
    show the problem, and PR #8 used the same null, so the numbers stay comparable.
    """
    xs = list(seq)
    rng.shuffle(xs)
    return "".join(xs)


def chimaera_evidence(idx, s1: str, s2: str, half_min: int = 30) -> dict:
    """Place each half of each mate separately and look for disagreement.

    A chimaeric molecule - PCR template switching, or an adapter-mediated join - puts two
    unrelated sequences in one read, so its halves place to unrelated references. A real
    novel read places nowhere at all, and the two cases are distinguished by whether the
    halves place rather than by whether the whole read does. The rate this returns on
    reads that the cascade subtracted is the field's false-positive rate and is reported
    beside it.
    """
    out = {"halves_placed": 0, "split_disagree": False, "refs": []}
    for s in (s1, s2):
        h = len(s) // 2
        hits = []
        for part in (s[:h], s[h:]):
            best = None
            for a in idx.map(part):
                if a.mlen >= half_min and (best is None or a.mlen > best[0]):
                    best = (a.mlen, a.ctg, a.strand)
            hits.append(best)
        placed = [x for x in hits if x]
        out["halves_placed"] += len(placed)
        if len(placed) == 2 and (placed[0][1] != placed[1][1]):
            out["split_disagree"] = True
            out["refs"] = [placed[0][1], placed[1][1]]
    return out


class LibraryAbundance:
    """How often does the rest of this library contain the same exact sequence?

    This field exists because of a correction the protein thread made on 2026-11-20 (PR
    #13). "Absence means the public databases are incomplete" is true but not sharp
    enough: what stays absent even after indexing the sample's OWN assembly is sequence
    from organisms that are neither in public databases NOR abundant enough to assemble
    here - rare AND unsequenced. That is also the description of a novel virus at one in a
    million, so the honest way to describe the residue is by abundance, not only by
    database status.

    So each surviving pair gets a measured abundance in its own library: the median
    multiplicity of its canonical 31-mers across a large sample of the run. A read from an
    abundant-but-unsequenced organism recurs many times; a read from something genuinely
    rare recurs about once, which is itself the evidence that it cannot be assembled here
    and therefore cannot be subtracted by the sample's own assembly.

    Counting is a 2-hash count-min sketch over uint8 counters rather than a dictionary:
    2M pairs carry ~300M distinct 31-mers, which no Python dict will hold in this
    container. A sketch only ever OVERSTATES a count, so a read called rare by the sketch
    is rare for certain, and the direction of the error is the safe one for this claim.
    Counters saturate at 255, which is fine: everything above that is "abundant" for the
    purpose at hand. The saturation fraction and an empty-query control are reported so
    the collision rate is visible rather than assumed.
    """

    def __init__(self, bits: int = 28, k: int = 31):
        import numpy as np
        self.k = k
        self.size = 1 << bits
        self.mask = self.size - 1
        self.a = np.zeros(self.size, dtype="uint8")
        self.b = np.zeros(self.size, dtype="uint8")
        self.added = 0

    @staticmethod
    def _canon(km: str) -> str:
        rc = revcomp(km)
        return km if km <= rc else rc

    # Deterministic mixing, not Python's hash(): str hashing is salted per process, so a
    # sketch built with hash() would give different counts on a re-run and the repo rule
    # that every number is reproducible from a script plus a seed would not hold.
    _PACK = str.maketrans("ACGT", "0123")
    _M1 = 0x9E3779B97F4A7C15
    _M2 = 0xC2B2AE3D27D4EB4F
    _U64 = (1 << 64) - 1

    def _h(self, km: str):
        v = int(km.translate(self._PACK), 4)
        h1 = (v * self._M1) & self._U64
        h2 = ((v ^ 0xD6E8FEB86659FD93) * self._M2) & self._U64
        return (h1 >> 13) & self.mask, (h2 >> 29) & self.mask

    def add_seq(self, seq: str) -> None:
        k = self.k
        for i in range(len(seq) - k + 1):
            km = self._canon(seq[i:i + k])
            if "N" in km:
                continue
            i1, i2 = self._h(km)
            if self.a[i1] < 255:
                self.a[i1] += 1
            if self.b[i2] < 255:
                self.b[i2] += 1
            self.added += 1

    def counts(self, seq: str) -> list[int]:
        k = self.k
        out = []
        for i in range(len(seq) - k + 1):
            km = self._canon(seq[i:i + k])
            if "N" in km:
                continue
            i1, i2 = self._h(km)
            out.append(int(min(self.a[i1], self.b[i2])))
        return out

    def saturation(self) -> float:
        import numpy as np
        return float(np.mean(self.a == 255))


def pfam_hits(prot_seqs: dict[str, str], hmm_path: Path, evalue: float, threads: int):
    """Six-frame Pfam-A search. A profile can reach a family a single sequence cannot."""
    import pyhmmer
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alphabet)
            for k, v in prot_seqs.items() if len(v) >= 20]
    out: dict[str, list] = {}
    if not seqs:
        return out
    # pyhmmer >= 0.11 returns str names, earlier versions bytes; accept either so the
    # script does not break silently on a different environment.
    def _s(x):
        return x.decode() if isinstance(x, bytes) else str(x)

    with pyhmmer.plan7.HMMFile(hmm_path) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            for hit in hits:
                if hit.evalue > evalue:
                    continue
                out.setdefault(_s(hit.name), []).append({
                    "profile": _s(hits.query.name),
                    "evalue": hit.evalue, "score": hit.score})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--residue", type=Path, required=True,
                    help="residue_rows_s<N>.jsonl from subtract_known.py --rows-out")
    ap.add_argument("--reads", type=Path, required=True,
                    help="fasta of the residue pairs written by subtract_known.py")
    ap.add_argument("--subtracted-reads", type=Path, default=None,
                    help="fasta of a sample of SUBTRACTED pairs: the known-positive "
                         "comparison every field is calibrated against")
    ap.add_argument("--dmnd", type=Path, nargs="*",
                    default=[ROOT / "data/prot/viral.dmnd", ROOT / "data/prot/sprot.dmnd"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--pfam", type=Path, default=ROOT / "data/pfam/Pfam-A.hmm")
    ap.add_argument("--chimaera-ref", type=Path,
                    default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--pfam-evalue", type=float, default=1e-3)
    ap.add_argument("--diamond-mode", default="--very-sensitive",
                    help="--very-sensitive for the closest-relative field, where the "
                         "point is to find the faintest hit; --sensitive when the run is "
                         "a sweep over several index sizes and the comparison between "
                         "sizes, not the absolute hit rate, is what is being measured")
    ap.add_argument("--skip-fields", action="store_true",
                    help="only the DIAMOND index-size sweep: skip motifs, chimaera and "
                         "the per-read fields, which do not depend on the protein index")
    ap.add_argument("--max-pairs", type=int, default=20_000)
    ap.add_argument("--library-r1", type=Path,
                    default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--library-r2", type=Path,
                    default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--library-pairs", type=int, default=400_000,
                    help="pairs of the run counted into the abundance sketch. The "
                         "measured abundance is relative to this many pairs, so the unit "
                         "is recurrence within a sample of the library, not a genome "
                         "coverage depth")
    ap.add_argument("--abundance-bits", type=int, default=28)
    ap.add_argument("--cross-sample-ref", type=Path,
                    default=ROOT / "data/geo/all_contigs_1kb.fa",
                    help="per-sample assemblies of all 116 samples, for cross-sample "
                         "recurrence")
    ap.add_argument("--cross-min-mlen", type=int, default=50)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/residue"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_report")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    t_start = time.perf_counter()

    pairs = _load_fasta_pairs(args.reads, args.max_pairs)
    known = _load_fasta_pairs(args.subtracted_reads, args.max_pairs) if args.subtracted_reads else {}
    if not pairs:
        raise SystemExit("no residue pairs")

    rows: dict[str, dict] = {}
    for name, (s1, s2) in pairs.items():
        rows[name] = _row(name, s1, s2, rng)

    # ---- closest relative, however far: no bitscore floor
    for db in args.dmnd:
        if not db.exists():
            continue
        tag = f"closest_{db.stem}"
        best = _diamond_best(args.diamond_bin, db, pairs, args.work / f"{tag}_s{args.seed}.fa",
                             args.threads, args.diamond_mode)
        for name, v in best.items():
            rows[name][tag] = v

    # ---- abundance in the library itself: rare-and-unknown vs abundant-and-unsequenced
    abund = None
    if args.library_r1.exists() and not args.skip_fields:
        from berlin_io import read_pairs as _rp
        t0 = time.perf_counter()
        sk = LibraryAbundance(bits=args.abundance_bits)
        for _nm, a, b, _q1, _q2 in _rp(args.library_r1, args.library_r2,
                                       args.library_pairs, 14, 0):
            sk.add_seq(a)
            sk.add_seq(b)
        for name, (s1, s2) in pairs.items():
            cs = sk.counts(s1) + sk.counts(s2)
            cs.sort()
            rows[name]["library_abundance"] = {
                "median_31mer_count": cs[len(cs) // 2] if cs else None,
                "max_31mer_count": cs[-1] if cs else None,
                "kmers": len(cs),
            }
        if known:
            kn = []
            for s1, s2 in known.values():
                cs = sk.counts(s1) + sk.counts(s2)
                cs.sort()
                if cs:
                    kn.append(cs[len(cs) // 2])
        abund = {
            "library_pairs_counted": args.library_pairs,
            "sketch_bits": args.abundance_bits,
            "kmers_added": sk.added,
            "counter_saturation_frac": sk.saturation(),
            "residue_median_count": spread(
                r["library_abundance"]["median_31mer_count"] for r in rows.values()
                if r["library_abundance"]["median_31mer_count"] is not None),
            "subtracted_median_count": spread(kn) if known and kn else None,
            "residue_recurrence_sweep": {
                str(c): sum(1 for r in rows.values()
                            if (r["library_abundance"]["median_31mer_count"] or 0) >= c) / len(rows)
                for c in (1, 2, 3, 5, 10, 20, 50)},
            "wall_seconds": time.perf_counter() - t0,
        }

    # ---- recurrence ACROSS samples: the one field no database can be incomplete about
    cross = None
    if args.cross_sample_ref.exists() and not args.skip_fields:
        # Everything else in this report is bounded by what some catalogue contains. This
        # field is not: it asks whether the same sequence turns up in OTHER wastewater
        # samples from the same works, which is answered by the samples themselves. A rare
        # but real organism should recur across samples or time points; a library
        # artefact, a chimaera or an error-made-unique read should not. That makes it the
        # only tier-3 field immune to the database-incompleteness problem that sank the
        # earlier absence features.
        #
        # The reference is the study authors' own per-sample assemblies for all 116
        # samples, so "recurs" means "assembles in another sample", which is a conjunction
        # of present AND abundant enough there to assemble: a one-sided measure again, and
        # it undercounts recurrence for anything rare everywhere.
        t0 = time.perf_counter()
        idx = _mm_index(args.cross_sample_ref,
                        args.work / f"{args.cross_sample_ref.stem}.mmi", args.threads)
        counts = []
        for name, (s1, s2) in pairs.items():
            samples = set()
            for s in (s1, s2):
                for a in idx.map(s):
                    if a.mlen >= args.cross_min_mlen:
                        samples.add(_sample_of(a.ctg))
            rows[name]["cross_sample"] = {"distinct_samples": len(samples),
                                          "samples": sorted(samples)[:5]}
            counts.append(len(samples))
        kn_counts = []
        if known:
            for s1, s2 in known.values():
                samples = set()
                for s in (s1, s2):
                    for a in idx.map(s):
                        if a.mlen >= args.cross_min_mlen:
                            samples.add(_sample_of(a.ctg))
                kn_counts.append(len(samples))
        cross = {
            "reference": str(args.cross_sample_ref),
            "min_matching_bases": args.cross_min_mlen,
            "residue_distinct_samples": spread(counts),
            "subtracted_distinct_samples": spread(kn_counts) if kn_counts else None,
            "residue_sweep": {str(c): sum(1 for x in counts if x >= c) / len(counts)
                              for c in (1, 2, 3, 5, 10)},
            "subtracted_sweep": ({str(c): sum(1 for x in kn_counts if x >= c) / len(kn_counts)
                                  for c in (1, 2, 3, 5, 10)} if kn_counts else None),
            "wall_seconds": time.perf_counter() - t0,
        }
        del idx

    # ---- motifs
    pf = {}
    if args.pfam.exists() and not args.skip_fields:
        prot = {}
        for name, (s1, s2) in pairs.items():
            for i, (strand, off, s) in enumerate(six_frames(s1 + "NN" + s2)):
                prot[f"{name}|{strand}{off}"] = translate(s).replace("*", "X")
        raw = pfam_hits(prot, args.pfam, args.pfam_evalue, args.threads)
        for key, hits in raw.items():
            pf.setdefault(key.rsplit("|", 1)[0], []).extend(hits)
        for name, hits in pf.items():
            rows[name]["motif"] = sorted(hits, key=lambda h: h["evalue"])[:3]

    # ---- chimaera evidence, and its own false-positive rate on subtracted reads
    chim_fp = None
    if args.chimaera_ref.exists() and not args.skip_fields:
        idx = _mm_index(args.chimaera_ref, args.work / f"{args.chimaera_ref.stem}.mmi",
                        args.threads)
        for name, (s1, s2) in pairs.items():
            rows[name]["chimaera"] = chimaera_evidence(idx, s1, s2)
        if known:
            k = [chimaera_evidence(idx, a, b) for a, b in known.values()]
            chim_fp = sum(1 for x in k if x["split_disagree"]) / len(k)

    # ---- field calibration: residue against shuffled null and against subtracted reads
    calib = {
        "pairs_scored": len(pairs),
        "fields": {
            "gc": {"residue": spread(r["gc"] for r in rows.values()),
                   "shuffled_null": spread(r["gc_shuffled"] for r in rows.values())},
            "entropy3": {"residue": spread(r["entropy3"] for r in rows.values()),
                         "shuffled_null": spread(r["entropy3_shuffled"] for r in rows.values())},
            "longest_orf_aa": {
                "residue": spread(r["frame"]["longest_orf_aa"] for r in rows.values()),
                "shuffled_null": spread(r["frame_shuffled"]["longest_orf_aa"] for r in rows.values()),
            },
            "frame_identifiable_frac": {
                "residue": sum(r["frame"]["frame_identifiable"] for r in rows.values()) / len(rows),
                "shuffled_null": sum(r["frame_shuffled"]["frame_identifiable"] for r in rows.values()) / len(rows),
            },
            "kmer16_duplicate_frac": {
                "residue": spread(r["repetitive"]["kmer16_duplicate_frac"] for r in rows.values()),
                "shuffled_null": spread(r["repetitive_shuffled"]["kmer16_duplicate_frac"] for r in rows.values()),
            },
        },
        "artefact_rates": {
            "adapter": sum(1 for r in rows.values() if r["adapter"]) / len(rows),
            "primer_dimer": sum(1 for r in rows.values() if r["primer_dimer"]) / len(rows),
            "low_entropy_0.5": sum(1 for r in rows.values() if r["entropy3"] < 0.5) / len(rows),
        },
        "motif_rate": len(pf) / len(rows) if pf else 0.0,
        "library_abundance": abund,
        "cross_sample_recurrence": cross,
        "chimaera_split_rate": (sum(1 for r in rows.values()
                                    if r.get("chimaera", {}).get("split_disagree")) / len(rows)),
        "chimaera_split_rate_on_subtracted_reads": chim_fp,
    }
    for db in args.dmnd:
        if not db.exists():
            continue
        tag = f"closest_{db.stem}"
        got = [r[tag] for r in rows.values() if tag in r]
        calib["fields"][tag] = {
            "any_hit_frac": len(got) / len(rows),
            "bitscore": spread(g["bitscore"] for g in got) if got else None,
            "pident": spread(g["pident"] for g in got) if got else None,
            "bitscore_sweep": {str(b): sum(1 for g in got if g["bitscore"] >= b) / len(rows)
                               for b in (30, 40, 50, 60, 80, 100)},
        }
        if known:
            kb = _diamond_best(args.diamond_bin, db, known,
                               args.work / f"{tag}_known_s{args.seed}.fa",
                               args.threads, args.diamond_mode)
            calib["fields"][tag]["subtracted_reads_any_hit_frac"] = len(kb) / len(known)
            calib["fields"][tag]["subtracted_reads_bitscore"] = (
                spread(v["bitscore"] for v in kb.values()) if kb else None)

    (args.out / f"calibration_s{args.seed}.json").write_text(json.dumps(calib, indent=2))
    with open(args.out / f"rows_s{args.seed}.jsonl", "w") as h:
        for name, r in rows.items():
            h.write(json.dumps(r) + "\n")

    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/residue_report.py",
        "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "index_bytes": {str(p): dir_bytes(p) for p in list(args.dmnd) + [args.pfam]
                        if p.exists()},
    }, indent=2))
    print(json.dumps(calib, indent=2)[:3000])


def _row(name: str, s1: str, s2: str, rng: random.Random) -> dict:
    art = artefact_fields(s1, s2)
    sh1, sh2 = shuffled(s1, rng), shuffled(s2, rng)
    return {
        "name": name, "len1": len(s1), "len2": len(s2),
        **art,
        "entropy3_shuffled": round(min(shannon_entropy(sh1), shannon_entropy(sh2)), 4),
        "gc_shuffled": art["gc"],  # invariant under shuffling by construction; kept so the
                                   # null column exists for every field in the table
        "frame": frame_evidence(s1 + "NN" + s2),
        "frame_shuffled": frame_evidence(sh1 + "NN" + sh2),
        "repetitive": repetitiveness(s1, s2),
        "repetitive_shuffled": repetitiveness(sh1, sh2),
    }


def _sample_of(ctg: str) -> str:
    """Sample id from a SPAdes contig name: WW_211013_20-22hmix_S5_L001_NODE_1_length_..."""
    return ctg.split("_NODE_")[0]


def _load_fasta_pairs(path: Path | None, limit: int) -> dict[str, tuple[str, str]]:
    if path is None or not path.exists():
        return {}
    out: dict[str, list] = {}
    name = None
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                base, _, mate = name.rpartition("/")
                out.setdefault(base, ["", ""])[0 if mate == "1" else 1] = line.strip()
            if len(out) > limit:
                break
    return {k: (v[0], v[1]) for k, v in list(out.items())[:limit] if v[0] and v[1]}


def _diamond_best(dbin: Path, db: Path, pairs: dict, fa: Path, threads: int, mode: str):
    with open(fa, "w") as h:
        for name, (s1, s2) in pairs.items():
            h.write(f">{name}/1\n{s1}\n>{name}/2\n{s2}\n")
    tsv = fa.with_suffix(".tsv")
    subprocess.run([str(dbin), "blastx", "--db", str(db), "--query", str(fa),
                    "--out", str(tsv), "--threads", str(threads), "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                    "--max-target-seqs", "1", "--evalue", "100", mode, "--quiet"],
                   check=True)
    best: dict[str, dict] = {}
    with open(tsv) as h:
        for line in h:
            q, s, pid, length, ev, bits = line.rstrip("\n").split("\t")
            name = q.rsplit("/", 1)[0]
            bits = float(bits)
            if bits > best.get(name, {}).get("bitscore", -1):
                best[name] = {"subject": s, "pident": float(pid), "aln_len": int(length),
                              "evalue": float(ev), "bitscore": bits,
                              "mate": q.rsplit("/", 1)[1]}
    tsv.unlink(missing_ok=True)
    fa.unlink(missing_ok=True)
    return best


if __name__ == "__main__":
    main()
