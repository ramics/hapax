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
from loud import Consulted, Preflight, require_output, run_tool
from subtract_known import (artefact_fields, dir_bytes, provenance, shannon_entropy,
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
    # The three stops. Their absence here was a real bug, caught on 2026-11-20 by a
    # calibration reading that could not be true: longest_orf_aa equalled the read length
    # for every read and stop_free_frames was 6 for random sequence. Without them
    # translate() mapped TAA/TAG/TGA to "X", no frame ever contained a stop, and every
    # coding field silently reported "this read is one long open reading frame".
    "TAA": "*", "TAG": "*", "TGA": "*",
}

# One source of truth for the stops, taken from the script whose frame work is already
# under test (PR #8, experiments/readlevel_features.py). A second private copy of a codon
# table is how the bug above happened, so this import plus the assertion below is the
# structural fix rather than the patch: if the two ever disagree, the import fails loudly
# at start-up instead of producing plausible numbers.
def _reference_stops() -> set[str]:
    """The stops as PR #8's tested frame code defines them.

    Imported through importlib with that script's own dependency directory on the path,
    because readlevel_features.py lives beside this file but imports a `common` module
    from the archived assembler's tree. A hard import at module scope would therefore
    couple this script's start-up to that layout; this reads the one constant and fails
    loudly if it cannot.
    """
    import importlib.util
    here = Path(__file__).resolve().parent
    sys.path.insert(0, str(here / "codon_assembler" / "experiments"))
    try:
        spec = importlib.util.spec_from_file_location(
            "_readlevel_features", here / "readlevel_features.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return set(mod.STOPS)
    finally:
        sys.path.pop(0)


assert {c for c, a in AA.items() if a == "*"} == _reference_stops(), (
    "codon table disagrees with readlevel_features.STOPS")


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


DINUCS = [a + b for a in "ACGT" for b in "ACGT"]


def dinucleotide_bias(seq: str) -> dict:
    """Observed-over-expected dinucleotide ratios, with CpG and UpA called out.

    Valid at 107 nt, which is why it is here: a read carries ~106 dinucleotides, so the
    16 ratios are estimated from enough observations to mean something, unlike codon usage
    (61 categories from ~35 observations). CpG and UpA depression is the classic
    host-adaptation signal in RNA viruses - vertebrate-infecting viruses track their
    host's suppression - so it carries information about WHERE a sequence has been
    replicating rather than about whether it is novel. It is reported for that reason, not
    as a novelty score: a bacteriophage and a vertebrate virus differ here, and both are
    equally unknown to a database.
    """
    n = len(seq)
    comp = {b: seq.count(b) / n for b in "ACGT"}
    obs = {d: 0 for d in DINUCS}
    tot = 0
    for i in range(n - 1):
        d = seq[i:i + 2]
        if d in obs:
            obs[d] += 1
            tot += 1
    out = {}
    for d in DINUCS:
        exp = comp[d[0]] * comp[d[1]]
        out[d] = round((obs[d] / tot) / exp, 4) if tot and exp > 0 else None
    return {"oe": out, "CpG_oe": out["CG"], "UpA_oe": out["TA"]}


def folding(seq: str, rng: random.Random, shuffles: int = 20) -> dict:
    """Minimum free energy of the RNA fold, against a dinucleotide-preserving control.

    The raw energy is almost useless on its own: it tracks GC content and length, so a
    GC-rich read folds "well" with no structure worth the name. What carries information
    is the z-score against sequences with the same dinucleotide composition, which is the
    standard control for structured RNA and needs NO database - the reason this field is
    worth having at all, given that database incompleteness has defeated nearly every
    other feature in this project.

    The control is a first-order Markov resample, which preserves dinucleotide frequencies
    in expectation rather than exactly. That is weaker than an Altschul-Erikson shuffle
    and is stated rather than glossed: it makes the z-score slightly optimistic, so a
    small positive result here should not be believed without the stronger control.
    """
    import RNA
    _, mfe = RNA.fold(seq.replace("T", "U"))
    nulls = []
    trans: dict[str, list[str]] = {}
    for i in range(len(seq) - 1):
        trans.setdefault(seq[i], []).append(seq[i + 1])
    for _ in range(shuffles):
        if not trans:
            break
        cur = seq[0]
        out = [cur]
        for _i in range(len(seq) - 1):
            nxt = rng.choice(trans.get(cur) or list(trans))
            out.append(nxt)
            cur = nxt
        _, m = RNA.fold("".join(out).replace("T", "U"))
        nulls.append(m)
    if len(nulls) < 2:
        return {"mfe": round(mfe, 2), "z": None}
    mu = sum(nulls) / len(nulls)
    sd = (sum((x - mu) ** 2 for x in nulls) / (len(nulls) - 1)) ** 0.5
    return {"mfe": round(mfe, 2), "null_mfe_mean": round(mu, 2),
            "z": round((mfe - mu) / sd, 3) if sd > 0 else None}


class SpacerIndex:
    """CRISPR spacers from this environment's own assembled bacteria.

    The best field in the non-coding ladder, because it carries POSITIVE identity with no
    database in it: a bacterium records the mobile elements that infected it, so a
    residual read matching a spacer found in this sewer's own assembled bacterial contigs
    is host-recorded evidence that the read is a phage or other mobile element, even when
    no catalogue anywhere holds a relative.

    Matching is by exact shared 20-mer in either orientation. Spacers are 26-50 nt, so a
    20-mer is most of one, and an exact 20-mer match between a 107 nt read and a spacer
    set this small is not a chance event: with S spacers of length L the expected number
    of chance matches per read is about 2 * S * (L - 19) * 88 / 4^20, which the report
    states numerically for the actual S rather than calling it negligible.
    """

    def __init__(self, fasta: Path, k: int = 20):
        self.k = k
        self.map: dict[str, str] = {}
        self.spacers = 0
        self.bases = 0
        name = None
        buf: list[str] = []
        with open(fasta) as h:
            for line in h:
                if line.startswith(">"):
                    if name:
                        self._add(name, "".join(buf))
                    name, buf = line[1:].strip(), []
                else:
                    buf.append(line.strip())
        if name:
            self._add(name, "".join(buf))

    def _add(self, name: str, seq: str) -> None:
        self.spacers += 1
        self.bases += len(seq)
        for s in (seq, revcomp(seq)):
            for i in range(len(s) - self.k + 1):
                self.map.setdefault(s[i:i + self.k], name)

    def chance_per_read(self, read_len: int = 107) -> float:
        kmers = max(read_len - self.k + 1, 0)
        return len(self.map) * kmers / (4 ** self.k)

    def match(self, s1: str, s2: str):
        for s in (s1, s2):
            for i in range(len(s) - self.k + 1):
                hit = self.map.get(s[i:i + self.k])
                if hit:
                    return hit
        return None


def rfam_scan(seqs: dict[str, str], cm: Path, cmscan: Path, threads: int, work: Path,
              evalue: float) -> dict[str, dict]:
    """Infernal structure search against Rfam, per sequence.

    Worth doing at read length, which was not obvious: Rfam's 4,227 covariance models have
    a median length of 88 and 66% are 107 nt or shorter, so most families CAN match a read
    from this library. Structure is conserved long after sequence is not, so a covariance
    model reaches relatives that neither a nucleotide nor a protein search can - and for
    the non-coding half of the residue it is the only homology route there is.

    Run at a stated E-value rather than Rfam's curated GA cutoffs, because GA thresholds
    are tuned for full-length sequence and a 107 nt fragment rarely clears them; the price
    is false positives, which is why the same scan is run on composition-matched shuffled
    reads to measure them directly instead of trusting the E-value.
    """
    fa = work / "rfam_query.fa"
    with open(fa, "w") as h:
        for k, v in seqs.items():
            h.write(f">{k}\n{v}\n")
    tbl = work / "rfam.tbl"
    subprocess.run([str(cmscan), "--rfam", "--noali", "--cpu", str(threads),
                    "-E", str(evalue), "--tblout", str(tbl), str(cm), str(fa)],
                   check=True, stdout=subprocess.DEVNULL)
    best: dict[str, dict] = {}
    with open(tbl) as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 17:
                continue
            model, query, ev = f[0], f[2], float(f[15])
            if query not in best or ev < best[query]["evalue"]:
                best[query] = {"model": model, "evalue": ev, "score": float(f[14])}
    fa.unlink(missing_ok=True)
    tbl.unlink(missing_ok=True)
    return best


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
    ap.add_argument("--known-motifs", action="store_true", default=True,
                    help="also run the Pfam pass on the subtracted pile, which is that "
                         "field's known-positive comparison")
    ap.add_argument("--no-known-motifs", dest="known_motifs",
                    action="store_false")
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
    ap.add_argument("--rfam", type=Path, default=ROOT / "data/rfam/Rfam.cm")
    ap.add_argument("--cmscan-bin", type=Path,
                    default=Path("/home/user/tools/inf/bin/cmscan"))
    ap.add_argument("--rfam-evalue", type=float, default=1e-3)
    ap.add_argument("--spacers", type=Path,
                    default=ROOT / "data/crispr/all_arrays_spacers.fa",
                    help="CRISPR spacers extracted from the assemblies by minced")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/residue"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_report")
    ap.add_argument("--without", nargs="+", default=[], choices=sorted(FIELD_INPUTS),
                    help="fields to leave out ON PURPOSE. Each is recorded as 'not "
                         "searched: --without' in the calibration file. A field whose "
                         "input is absent and is not named here is a refusal to start, "
                         "not a field quietly dropped")
    ap.add_argument("--preflight", action="store_true",
                    help="check every input, index and tool, then exit: 0 if the run "
                         "could start, 2 if not")
    args = ap.parse_args()

    # A field whose reference is absent used to be skipped by `if path.exists()`, and its
    # absence then read downstream as a measured zero (`motif_rate` was 0.0 with no Pfam
    # on disk). Now every field is either computed, or left out by name and recorded so.
    pf = Preflight("residue_report.py")
    if not args.preflight:
        # The read piles are an earlier stage's products, so a preflight run before that
        # stage cannot check them; everything else it can.
        pf.file(args.reads, "--reads")
        if args.subtracted_reads is not None:
            pf.file(args.subtracted_reads, "--subtracted-reads")
    for db in args.dmnd:
        pf.file(db, "--dmnd")
    if args.dmnd:
        pf.exe(args.diamond_bin, "--diamond-bin")
    for field, inputs in FIELD_INPUTS.items():
        if not field_on(args, field):
            continue
        for attr, kind in inputs:
            v = attr if kind == "module" else getattr(args, attr)
            if kind == "exe":
                pf.exe(v, f"{field} field (--{attr.replace('_', '-')}); or pass "
                          f"--without {field}")
            elif kind == "module":
                pf.module(v, f"{field} field; or pass --without {field}")
            else:
                pf.file(v, f"{field} field (--{attr.replace('_', '-')}); or pass "
                           f"--without {field}")
    pf.module("RNA", "fold field (ViennaRNA), computed on every row")
    pf.done(args.preflight)
    consulted = Consulted()
    field_state = {f: ("computed" if field_on(args, f) else
                       "not searched: --skip-fields" if args.skip_fields else
                       "not searched: --without")
                   for f in FIELD_INPUTS}

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
        rows[name]["set"] = "residue"

    # The same per-read fields on reads the cascade DID identify. Written out beside the
    # residue rows, tagged by `set`, because every statement of the form "the residue has
    # property X" needs the contrast population in the same file at the same seed - and
    # because a sibling thread joins on these rows by read name. The expensive per-read
    # work (DIAMOND, Rfam, Pfam) stays residue-only; what is written here are the fields
    # that cost nothing extra once the read is in hand.
    known_rows: dict[str, dict] = {}
    for name, (s1, s2) in known.items():
        known_rows[name] = _row(name, s1, s2, rng)
        known_rows[name]["set"] = "known"

    # ---- closest relative, however far: no bitscore floor
    for db in args.dmnd:
        tag = f"closest_{db.stem}"
        best = _diamond_best(args.diamond_bin, db, pairs, args.work / f"{tag}_s{args.seed}.fa",
                             args.threads, args.diamond_mode)
        for name, v in best.items():
            rows[name][tag] = v
        consulted.searched(f"{tag} residue", db, queries=len(pairs), hits=len(best))
    if not args.dmnd:
        consulted.not_searched("closest_*", None, "no --dmnd given: no read in this run "
                               "was searched for a relative")

    # ---- abundance in the library itself: rare-and-unknown vs abundant-and-unsequenced
    abund = None
    if field_on(args, "abundance"):
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
            for kname, (s1, s2) in known.items():
                cs = sk.counts(s1) + sk.counts(s2)
                cs.sort()
                if cs:
                    kn.append(cs[len(cs) // 2])
                known_rows[kname]["library_abundance"] = {
                    "median_31mer_count": cs[len(cs) // 2] if cs else None,
                    "max_31mer_count": cs[-1] if cs else None,
                    "kmers": len(cs),
                }
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
    if field_on(args, "cross_sample"):
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
        consulted.searched("cross_sample", args.cross_sample_ref, queries=len(pairs),
                           hits=sum(1 for x in counts if x))
        del idx

    # ---- Rfam structure search, with a shuffled control for its false-positive rate
    rfam = None
    if field_on(args, "rfam"):
        t0 = time.perf_counter()
        q = {}
        for name, (s1, s2) in pairs.items():
            q[f"{name}/1"] = s1
            q[f"{name}/2"] = s2
        hits = rfam_scan(q, args.rfam, args.cmscan_bin, args.threads, args.work,
                         args.rfam_evalue)
        per_pair = {}
        for key, v in hits.items():
            nm = key.rsplit("/", 1)[0]
            if nm not in per_pair or v["evalue"] < per_pair[nm]["evalue"]:
                per_pair[nm] = v
        for nm, v in per_pair.items():
            rows[nm]["rfam"] = v
        null_q = {f"{nm}/1": shuffled(s1, rng) for nm, (s1, _s2) in pairs.items()}
        null_hits = rfam_scan(null_q, args.rfam, args.cmscan_bin, args.threads, args.work,
                              args.rfam_evalue)
        kn_hits = {}
        if known:
            kq = {f"{nm}/1": a for nm, (a, _b) in known.items()}
            kn_hits = rfam_scan(kq, args.rfam, args.cmscan_bin, args.threads, args.work,
                                args.rfam_evalue)
        rfam = {
            "evalue": args.rfam_evalue,
            "models": _cm_count(args.rfam),
            "residue_hit_frac": len(per_pair) / len(pairs),
            "shuffled_control_hit_frac": len(null_hits) / len(pairs),
            "subtracted_hit_frac": (len(kn_hits) / len(known)) if known else None,
            "top_models": Counter(v["model"] for v in per_pair.values()).most_common(10),
            "wall_seconds": time.perf_counter() - t0,
        }

    # ---- CRISPR spacers: positive identity with no database in it
    crispr = None
    if field_on(args, "crispr"):
        si = SpacerIndex(args.spacers)
        hits = 0
        for name, (s1, s2) in pairs.items():
            m = si.match(s1, s2)
            if m:
                rows[name]["crispr_spacer"] = m
                hits += 1
        kn_hits = sum(1 for a, b in known.values() if si.match(a, b)) if known else None
        crispr = {
            "spacer_file": str(args.spacers),
            "spacers": si.spacers,
            "spacer_bases": si.bases,
            "distinct_20mers": len(si.map),
            "expected_chance_matches_per_read": si.chance_per_read(),
            "residue_hit_frac": hits / len(pairs),
            "subtracted_hit_frac": (kn_hits / len(known)) if known else None,
        }

    # ---- motifs
    pf = {}
    pf_known_rate = None
    if field_on(args, "pfam"):
        prot = {}
        for name, (s1, s2) in pairs.items():
            for i, (strand, off, s) in enumerate(six_frames(s1 + "NN" + s2)):
                prot[f"{name}|{strand}{off}"] = translate(s).replace("*", "X")
        raw = pfam_hits(prot, args.pfam, args.pfam_evalue, args.threads)
        for key, hits in raw.items():
            pf.setdefault(key.rsplit("|", 1)[0], []).extend(hits)
        for name, hits in pf.items():
            rows[name]["motif"] = sorted(hits, key=lambda h: h["evalue"])[:3]
        # The same profiles on the known pile. A profile hit is the one field that can reach
        # a relative too divergent for DIAMOND, which makes it the field whose
        # known-positive rate matters most: if a Pfam model hits identified wastewater reads
        # at the same rate, "a profile matches" says nothing about the residue. It costs a
        # second pyhmmer pass over the control pile, which is why it is optional, but an
        # uncalibrated field is worse than a slower run.
        if known and args.known_motifs:
            kprot = {}
            for name, (s1, s2) in known.items():
                for strand, off, s in six_frames(s1 + "NN" + s2):
                    kprot[f"{name}|{strand}{off}"] = translate(s).replace("*", "X")
            kraw = pfam_hits(kprot, args.pfam, args.pfam_evalue, args.threads)
            kpf: dict[str, list] = {}
            for key, hits in kraw.items():
                kpf.setdefault(key.rsplit("|", 1)[0], []).extend(hits)
            for name, hits in kpf.items():
                if name in known_rows:
                    known_rows[name]["motif"] = sorted(hits,
                                                       key=lambda h: h["evalue"])[:3]
            pf_known_rate = len(kpf) / len(known)

    # ---- chimaera evidence, and its own false-positive rate on subtracted reads
    chim_fp = None
    if field_on(args, "chimaera"):
        idx = _mm_index(args.chimaera_ref, args.work / f"{args.chimaera_ref.stem}.mmi",
                        args.threads)
        for name, (s1, s2) in pairs.items():
            rows[name]["chimaera"] = chimaera_evidence(idx, s1, s2)
        if known:
            # Per row as well as pooled: the pooled rate is the test's false-positive rate,
            # which is what makes a read-level chimaera call readable, and the per-row copy
            # is what lets anything downstream compare the two piles read by read.
            for name, (a, b) in known.items():
                known_rows[name]["chimaera"] = chimaera_evidence(idx, a, b)
            k = [known_rows[name]["chimaera"] for name in known]
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
        # 0.0 only when Pfam ran and matched nothing; otherwise the reason it did not run.
        "motif_rate": (len(pf) / len(rows) if field_on(args, "pfam")
                       else field_state["pfam"]),
        "motif_rate_on_subtracted_reads": pf_known_rate,
        "library_abundance": abund,
        "crispr": crispr,
        "rfam": rfam,
        "dinucleotide": {
            "CpG_oe_residue": spread(r["dinuc"]["CpG_oe"] for r in rows.values()
                                     if r["dinuc"]["CpG_oe"] is not None),
            "UpA_oe_residue": spread(r["dinuc"]["UpA_oe"] for r in rows.values()
                                     if r["dinuc"]["UpA_oe"] is not None),
        },
        "folding": {
            "mfe_residue": spread(r["fold"]["mfe"] for r in rows.values()),
            "z_residue": spread(r["fold"]["z"] for r in rows.values()
                                if r["fold"]["z"] is not None),
            "frac_z_under_minus2": (sum(1 for r in rows.values()
                                        if (r["fold"]["z"] or 0) <= -2) / len(rows)),
        },
        "cross_sample_recurrence": cross,
        "chimaera_split_rate": (sum(1 for r in rows.values()
                                    if r.get("chimaera", {}).get("split_disagree")) / len(rows)),
        "chimaera_split_rate_on_subtracted_reads": chim_fp,
    }
    for db in args.dmnd:
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
            consulted.searched(f"{tag} known", db, queries=len(known), hits=len(kb))
            calib["fields"][tag]["subtracted_reads_any_hit_frac"] = len(kb) / len(known)
            calib["fields"][tag]["subtracted_reads_bitscore"] = (
                spread(v["bitscore"] for v in kb.values()) if kb else None)
            # The hit is attached to the KNOWN ROW as well as summarised, because a summary
            # cannot be compared per read and a missing field in the control pile reads as
            # a zero rather than as "not measured". Building `hapax run` on these rows
            # surfaced exactly that: five of six residue layers were reporting a perfect
            # separation from the known pile that was only the search never having been
            # written down for it. The search already ran; this costs nothing but the key.
            for kname, v in kb.items():
                if kname in known_rows:
                    known_rows[kname][tag] = v

    # Which fields ran, and which piles each search covered. A reader of rows_s*.jsonl
    # cannot tell a missing `closest_*` key (searched, no hit) from an index never searched
    # without this; `hapax.py strains` reads it rather than inferring from hits.
    calib["field_state"] = field_state
    calib["consulted"] = consulted.as_list()
    calib["known_pile_pairs"] = len(known)
    (args.out / f"calibration_s{args.seed}.json").write_text(json.dumps(calib, indent=2))
    with open(args.out / f"rows_s{args.seed}.jsonl", "w") as h:
        for name, r in rows.items():
            h.write(json.dumps(r) + "\n")
        for name, r in known_rows.items():
            h.write(json.dumps(r) + "\n")

    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/residue_report.py",
        **provenance(__file__),
        "consulted": consulted.as_list(),
        "field_state": field_state,
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


#: The optional per-read fields and the inputs each one reads. Every input listed is
#: checked before the run starts; a field is computed unless named in --without (or
#: --skip-fields is given), and never dropped because its input happened to be absent.
FIELD_INPUTS = {
    "abundance": [("library_r1", "file"), ("library_r2", "file")],
    "cross_sample": [("cross_sample_ref", "file"), ("mappy", "module")],
    "rfam": [("rfam", "file"), ("cmscan_bin", "exe")],
    "crispr": [("spacers", "file")],
    "pfam": [("pfam", "file"), ("pyhmmer", "module")],
    "chimaera": [("chimaera_ref", "file"), ("mappy", "module")],
}


def field_on(args, field: str) -> bool:
    return not args.skip_fields and field not in args.without


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
        "dinuc": dinucleotide_bias(s1 + s2),
        "fold": folding(s1, rng),
    }


def _cm_count(cm: Path) -> int:
    try:
        return sum(1 for line in open(cm) if line.startswith("CLEN"))
    except Exception:
        return -1


def _sample_of(ctg: str) -> str:
    """Sample id from a SPAdes contig name: WW_211013_20-22hmix_S5_L001_NODE_1_length_..."""
    return ctg.split("_NODE_")[0]


def _load_fasta_pairs(path: Path | None, limit: int) -> dict[str, tuple[str, str]]:
    """Pairs from a FASTA of /1 /2 records. An absent file is fatal, not an empty pile:
    returning {} made a missing control pile read as "no known-positive in this run"."""
    if path is None or not path.exists():
        raise SystemExit(f"{path}: pair FASTA is absent")
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
    kept = list(out.items())[:limit]
    broken = [k for k, v in kept if not (v[0] and v[1])]
    if broken:
        # The cascade writes both mates of every pair, so a pair with one mate is a
        # truncated or foreign file, and dropping it would shrink the pile unannounced.
        raise SystemExit(f"{path}: {len(broken)} of {len(kept)} pairs lack a mate "
                         f"(first: {broken[0]})")
    return {k: (v[0], v[1]) for k, v in kept}


def _diamond_best(dbin: Path, db: Path, pairs: dict, fa: Path, threads: int, mode: str):
    with open(fa, "w") as h:
        for name, (s1, s2) in pairs.items():
            h.write(f">{name}/1\n{s1}\n>{name}/2\n{s2}\n")
    tsv = fa.with_suffix(".tsv")
    tsv.unlink(missing_ok=True)
    run_tool([str(dbin), "blastx", "--db", str(db), "--query", str(fa),
              "--out", str(tsv), "--threads", str(threads), "--outfmt", "6",
              "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
              "--max-target-seqs", "1", "--evalue", "100", mode, "--quiet"],
             f"DIAMOND blastx against {db.name}")
    require_output(tsv, f"DIAMOND {db.name}", allow_empty=True)
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
