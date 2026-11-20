"""The noise floor under peptide k-mer absence, measured before any index is built.

Why this exists
---------------
`results/readlevel_features/PRIOR_ART.md` left one quantity unmeasured and called it
"structurally worse than the nucleotide case": how far a single sequencing error
propagates in peptide k-mer space. It is the noise floor under the whole peptide-absence
proposal, because absence is the signal we intend to read as novelty and sequencing error
manufactures absence for free. One substitution changes at most one codon, and a changed
amino acid corrupts every peptide k-mer overlapping it — up to k_aa of them. At the 107 nt
read length of SRR23973193 a read carries only ~26 peptide 10-mers in its true frame, so a
handful of errors could in principle corrupt most of them.

If error alone drives a read from a *known* organism to look mostly absent, there is no
headroom left for real novelty and no proteome index is worth downloading. That is the kill
criterion, and it costs nothing but CPU to check, which is why it runs before the index.

What this can and cannot settle
-------------------------------
It bounds the FALSE-POSITIVE side only, and it does so optimistically on purpose.

The idealisation is that a read from a known organism has every one of its error-free
peptide k-mers present in the database, so that the k-mers which come back absent are
exactly the ones error corrupted. Real life adds absence that this does not model:
undersequenced known strains, low-complexity peptides, reads from non-coding sequence, and
the five wrong frames. So the false-positive rates reported here are a FLOOR — the real
feature can only be noisier. That makes a failure here decisive and a pass here merely
permissive.

It says nothing about sensitivity. A genuinely novel read's absent fraction depends on how
far its proteins sit from the nearest retained relative, which needs the index. So the
readout is deliberately framed as HEADROOM: the absence threshold at which error-induced
false flags drop below a given per-million budget, and how many of the read's k-mers are
left above that threshold for the real signal to use.

Sections
--------
1. Synonymous shielding. The genetic code is redundant, so a substantial fraction of
   substitutions leave the peptide untouched. This is the effect that could make peptide
   space *more* error-robust than nucleotide space per base, and it runs counter to the
   worry recorded in PRIOR_ART.md. Measured on real CDS rather than taken from the codon
   table, because codon usage is not uniform.
2. Corruption per read, peptide vs nucleotide at matched k (k_nt = 3 * k_aa). The direct
   comparison. Both are expressed as the fraction of a read's k-mers that differ from the
   error-free truth.
3. Indels. A substitution is local; an indel frameshifts everything downstream and destroys
   every peptide k-mer after it while leaving nucleotide k-mers only locally disturbed.
   Swept from Illumina rates to nanopore rates, because this is the term that decides
   whether the feature is Illumina-only.
4. Headroom. Section 2 turned into the number that matters: error-induced false-flag rate
   per million as a function of the absence threshold, read against the 1,489/M floor that
   translated search already set in results/abundance_sweep/REPORT.md.

Error model
-----------
Substitutions at a flat per-base rate, and separately on a 3'-ramped profile with the same
mean, because real Illumina error concentrates toward the end of the read. Clustering
matters in the direction people do not expect: errors bunched together corrupt fewer
DISTINCT k-mers than the same number of errors spread out, so a flat model is the
pessimistic one and the ramp is reported beside it rather than instead of it.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STOPS = {"TAA", "TAG", "TGA"}
BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

CODON_TABLE = {}
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _aa in enumerate(_AA):
    _c = "TCAG"[_i >> 4] + "TCAG"[(_i >> 2) & 3] + "TCAG"[_i & 3]
    CODON_TABLE[_c] = _aa

# Diverse small genomes with CDS annotation, spanning GC content and genome type.
# Kept small so the whole section runs in minutes; the quantity being measured is a
# property of the genetic code and the error model, not of any one organism, so breadth
# across codon usage matters more than depth in any genome.
GENOMES = {
    "hiv1_hxb2": "K03455.1",
    "rabies": "NC_001542.1",
    "astrovirus": "NC_001943.1",
    "hbv": "NC_003977.2",
    "lambda": "NC_001416.1",
    "phix174": "NC_001422.1",
}


# ---------------------------------------------------------------- fetching / parsing

def fetch_gb(accession: str, dest: Path) -> Path:
    """Fetch a GenBank record, caching in data/. data/ is gitignored, so an experiment
    that needs sequence fetches it rather than assuming a previous run left it behind."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype=gbwithparts&retmode=text")
    # E-utilities rate-limits shared egress IPs hard, so back off rather than hammer.
    # No contact address is sent: the sandbox's outbound identity is not the user's.
    req = urllib.request.Request(url, headers={"User-Agent": "hapax-experiment/1.0"})
    last = None
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=120) as fh:
                body = fh.read()
            if body.strip():
                dest.write_bytes(body)
                time.sleep(0.5)
                return dest
            last = "empty body"
        except Exception as exc:  # noqa: BLE001 - retried below, raised if persistent
            last = exc
        time.sleep(2 ** attempt)
    raise RuntimeError(f"could not fetch {accession}: {last}")


def parse_gb(path: Path) -> tuple[str, list[tuple[int, int, bool]]]:
    """Return (sequence, [(start, end, is_complement)]) for simple CDS features.

    Joined and fuzzy locations are skipped rather than guessed at: a spliced CDS would
    give a wrong frame, and a wrong frame silently poisons every number downstream. The
    CDS invariant check in main() is what catches a parsing error of this kind, following
    the lesson recorded in results/readlevel_features/REPORT.md.
    """
    text = path.read_text()
    seq_part = text.split("ORIGIN", 1)
    seq = ""
    if len(seq_part) == 2:
        seq = "".join(re.findall(r"[acgtnACGTN]", seq_part[1].split("//")[0])).upper()
    feats = []
    for m in re.finditer(r"^     CDS             (\S.*)$", text, re.M):
        loc = m.group(1).strip()
        comp = loc.startswith("complement(")
        if comp:
            loc = loc[len("complement("):-1]
        if "join" in loc or "<" in loc or ">" in loc or ".." not in loc:
            continue
        try:
            a, b = loc.split("..")
            feats.append((int(a) - 1, int(b), comp))
        except ValueError:
            continue
    return seq, feats


def revcomp(s: str) -> str:
    return "".join(COMP.get(c, "N") for c in reversed(s))


def translate(nt: str) -> str:
    return "".join(CODON_TABLE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def clean_cds(seq: str, feats: list[tuple[int, int, bool]]) -> list[str]:
    """CDS sequences that are usable as ground truth: in frame, no ambiguity, no internal
    stop. A CDS failing these is an annotation we cannot trust to define a true frame."""
    out = []
    for a, b, comp in feats:
        s = seq[a:b]
        if comp:
            s = revcomp(s)
        if len(s) < 120 or len(s) % 3 or set(s) - set(BASES):
            continue
        prot = translate(s)
        if "*" in prot[:-1] or "X" in prot:
            continue
        out.append(s)
    return out


# ---------------------------------------------------------------- error model

def ramp_profile(read_len: int, mean_rate: float) -> list[float]:
    """Per-position substitution probability rising linearly toward the 3' end with the
    requested mean. Shape taken as linear rather than fitted: the point is to show the
    direction clustering moves the answer, not to model a specific instrument."""
    lo, hi = 0.25, 1.75  # mean 1.0, so scaling by mean_rate preserves the mean
    return [mean_rate * (lo + (hi - lo) * i / max(1, read_len - 1)) for i in range(read_len)]


def mutate(read: str, rng: random.Random, profile: list[float],
           indel_rate: float = 0.0) -> tuple[str, int, int]:
    """Apply substitutions per the profile and indels at a flat rate.

    Returns (mutated read, n_substitutions, n_indels). The read is NOT re-trimmed to the
    original length after an indel: downstream code works in k-mer counts, and padding or
    truncating here would hide the frameshift that is the whole point of section 3.
    """
    out = []
    n_sub = n_indel = 0
    for i, base in enumerate(read):
        if indel_rate and rng.random() < indel_rate:
            n_indel += 1
            if rng.random() < 0.5:
                continue  # deletion
            out.append(rng.choice(BASES))  # insertion, then the original base follows
        if rng.random() < profile[i]:
            out.append(rng.choice([b for b in BASES if b != base]))
            n_sub += 1
        else:
            out.append(base)
    return "".join(out), n_sub, n_indel


# ---------------------------------------------------------------- k-mer corruption

def kmers(s: str, k: int) -> list[str]:
    return [s[i:i + k] for i in range(len(s) - k + 1)] if len(s) >= k else []


def corrupted_fraction(clean: str, dirty: str, k: int) -> tuple[int, int]:
    """(corrupted k-mer count, clean k-mer count), positionally compared.

    Position i of the dirty sequence is compared with position i of the clean one. Under
    an indel the sequences fall out of register and essentially every downstream k-mer
    differs, which is the intended behaviour: in a database lookup a frameshifted peptide
    really is absent. A set-membership comparison would charitably re-align them and
    understate the damage.
    """
    ck = kmers(clean, k)
    dk = kmers(dirty, k)
    if not ck:
        return 0, 0
    bad = 0
    for i, km in enumerate(ck):
        if i >= len(dk) or dk[i] != km:
            bad += 1
    return bad, len(ck)


def read_metrics(clean_read: str, dirty_read: str, frame_off: int,
                 k_aa: int) -> dict:
    """Peptide and nucleotide corruption for one read, in its known true frame.

    frame_off is the read's offset into the source CDS reading frame, so the true frame is
    known by construction rather than inferred — which matters because 47% of coding reads
    at 107 nt cannot identify their own frame (results/readlevel_features/REPORT.md). Every
    number here is therefore an upper bound on what a frame-blind scorer could achieve.
    """
    k_nt = 3 * k_aa
    c_pep = translate(clean_read[frame_off:])
    d_pep = translate(dirty_read[frame_off:])
    pep_bad, pep_n = corrupted_fraction(c_pep, d_pep, k_aa)
    nt_bad, nt_n = corrupted_fraction(clean_read, dirty_read, k_nt)
    return {"pep_bad": pep_bad, "pep_n": pep_n, "nt_bad": nt_bad, "nt_n": nt_n}


# ---------------------------------------------------------------- sections

def sample_cds_reads(cdss: list[str], read_len: int, n: int,
                     rng: random.Random) -> list[tuple[str, int]]:
    """Reads drawn inside a CDS, each with its offset into the reading frame.

    Reads start at arbitrary positions, not codon boundaries, because real reads do; the
    offset is carried so the true frame stays known.
    """
    usable = [c for c in cdss if len(c) > read_len + 3]
    out = []
    if not usable:
        return out
    while len(out) < n:
        cds = usable[rng.randrange(len(usable))]
        i = rng.randrange(len(cds) - read_len)
        out.append((cds[i:i + read_len], (-i) % 3))
    return out


def section_synonymous(cdss: list[str], rng: random.Random, n: int) -> dict:
    """What fraction of random substitutions in real CDS leave the peptide unchanged.

    This is the term that could make peptide absence MORE error-robust than nucleotide
    absence, contrary to the worry in PRIOR_ART.md. Measured on real codon usage because
    the uniform-codon expectation is not what a genome does.
    """
    silent = total = 0
    by_pos = {0: [0, 0], 1: [0, 0], 2: [0, 0]}
    for _ in range(n):
        cds = cdss[rng.randrange(len(cdss))]
        i = rng.randrange(len(cds))
        orig = cds[i]
        new = rng.choice([b for b in BASES if b != orig])
        c_start = (i // 3) * 3
        codon = cds[c_start:c_start + 3]
        if len(codon) < 3:
            continue
        mut = list(codon)
        mut[i - c_start] = new
        mut = "".join(mut)
        same = CODON_TABLE.get(codon, "X") == CODON_TABLE.get(mut, "X")
        total += 1
        by_pos[i - c_start][1] += 1
        if same:
            silent += 1
            by_pos[i - c_start][0] += 1
    return {
        "n_substitutions": total,
        "synonymous_fraction": silent / total if total else 0.0,
        "synonymous_fraction_by_codon_position": {
            str(p): (v[0] / v[1] if v[1] else 0.0) for p, v in by_pos.items()
        },
    }


def section_corruption(cdss: list[str], read_len: int, k_aa: int, rates: list[float],
                       n_reads: int, rng: random.Random) -> list[dict]:
    rows = []
    for rate in rates:
        for shape in ("flat", "ramped"):
            profile = ([rate] * read_len if shape == "flat"
                       else ramp_profile(read_len, rate))
            pep_fracs, nt_fracs, pep_counts, subs, pep_totals = [], [], [], [], []
            for clean, off in sample_cds_reads(cdss, read_len, n_reads, rng):
                dirty, n_sub, _ = mutate(clean, rng, profile)
                m = read_metrics(clean, dirty, off, k_aa)
                if m["pep_n"] == 0 or m["nt_n"] == 0:
                    continue
                pep_fracs.append(m["pep_bad"] / m["pep_n"])
                nt_fracs.append(m["nt_bad"] / m["nt_n"])
                pep_counts.append(m["pep_bad"])
                pep_totals.append(m["pep_n"])
                subs.append(n_sub)
            rows.append({
                "error_rate": rate, "profile": shape, "read_len": read_len, "k_aa": k_aa,
                "n_reads": len(pep_fracs),
                "mean_substitutions_per_read": _mean(subs),
                "peptide_corrupted_fraction": _mean(pep_fracs),
                "nucleotide_corrupted_fraction": _mean(nt_fracs),
                "peptide_minus_nucleotide": _mean(pep_fracs) - _mean(nt_fracs),
                "peptide_corrupted_count_mean": _mean(pep_counts),
                "reads_with_zero_corrupted": sum(1 for c in pep_counts if c == 0) / len(pep_counts) if pep_counts else 0.0,
                "kmers_per_read_mean": _mean(pep_totals),
                "kmers_per_read_min": min(pep_totals) if pep_totals else 0,
                "_pep_counts": pep_counts,
            })
    return rows


def section_indels(cdss: list[str], read_len: int, k_aa: int, sub_rate: float,
                   indel_rates: list[float], n_reads: int, rng: random.Random) -> list[dict]:
    rows = []
    profile = [sub_rate] * read_len
    for ir in indel_rates:
        pep_fracs, nt_fracs, n_ind = [], [], []
        for clean, off in sample_cds_reads(cdss, read_len, n_reads, rng):
            dirty, _, ni = mutate(clean, rng, profile, indel_rate=ir)
            m = read_metrics(clean, dirty, off, k_aa)
            if m["pep_n"] == 0 or m["nt_n"] == 0:
                continue
            pep_fracs.append(m["pep_bad"] / m["pep_n"])
            nt_fracs.append(m["nt_bad"] / m["nt_n"])
            n_ind.append(ni)
        rows.append({
            "indel_rate": ir, "sub_rate": sub_rate, "n_reads": len(pep_fracs),
            "mean_indels_per_read": _mean(n_ind),
            "reads_with_an_indel": sum(1 for x in n_ind if x) / len(n_ind) if n_ind else 0.0,
            "peptide_corrupted_fraction": _mean(pep_fracs),
            "nucleotide_corrupted_fraction": _mean(nt_fracs),
        })
    return rows


def tail_distribution(n_aa: int, k_aa: int, read_len: int, sub_rate: float,
                      p_nonsyn: float, n_samples: int, rng: random.Random) -> list[int]:
    """Corrupted-k-mer counts sampled fast enough to resolve a per-million tail.

    The direct simulation in section 2 cannot see the tail that matters. At 20,000 reads
    the smallest non-zero rate observable is 50 per million, so a threshold reported as
    meeting a 1-per-million budget would really be saying "no read in the sample did
    this" — exactly the kind of zero-count claim the project's process rules warn about.

    So the tail is sampled through the structure of the problem instead of through string
    manipulation. Given the set P of amino-acid positions that actually changed, the
    corrupted k-mer count is determined: a k-mer is corrupted exactly when its window
    intersects P. That makes a sample cheap, and the only inputs are the substitution
    count (Poisson in the read length) and the measured probability that a substitution is
    non-synonymous. Validated against the direct simulation by the caller.
    """
    n_pep = max(0, n_aa - k_aa + 1)
    lam = read_len * sub_rate
    counts = []
    for _ in range(n_samples):
        n_sub = _poisson(lam, rng)
        changed = set()
        for _ in range(n_sub):
            if rng.random() < p_nonsyn:
                changed.add(rng.randrange(n_aa))
        if not changed:
            counts.append(0)
            continue
        bad = sum(1 for i in range(n_pep)
                  if any(i <= c < i + k_aa for c in changed))
        counts.append(bad)
    return counts


def _poisson(lam: float, rng: random.Random) -> int:
    # Knuth's method; lam is well under 2 here, so the loop is short.
    el, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= el:
            return k
        k += 1


def section_tail(read_len: int, k_aa: int, rates: list[float], p_nonsyn: float,
                 budgets: list[float], n_samples: int, rng: random.Random,
                 direct: list[dict]) -> dict:
    """Headroom from the fast sampler, with a validation against the direct simulation."""
    n_aa = read_len // 3
    rows, checks = [], []
    for rate in rates:
        counts = tail_distribution(n_aa, k_aa, read_len, rate, p_nonsyn, n_samples, rng)
        n = len(counts)
        n_pep = max(0, n_aa - k_aa + 1)
        entry = {"error_rate": rate, "kmers_per_read": n_pep,
                 "n_samples": n, "resolution_per_million": 1e6 / n,
                 "mean_corrupted_fraction": _mean([c / n_pep for c in counts]) if n_pep else 0.0,
                 "thresholds": {}}
        for budget in budgets:
            t_needed, fp_at_t = None, None
            for T in range(0, n_pep + 2):
                fp = sum(1 for c in counts if c >= T) / n
                if fp * 1e6 <= budget:
                    t_needed, fp_at_t = T, fp * 1e6
                    break
            entry["thresholds"][str(budget)] = {
                "min_threshold": t_needed,
                "headroom_kmers": (n_pep - t_needed) if t_needed is not None else None,
                "achieved_fp_per_million": fp_at_t,
                "resolved": (fp_at_t is not None and fp_at_t > 0),
            }
        # What a mate pair buys. Sequencing errors in two mates are separate events on
        # separate molecules of the library, so for the ERROR component of a false flag
        # the mate correlation rho_fp is zero by construction — the one component of the
        # absence score whose correlation we do not have to measure. Requiring both mates
        # to clear the threshold therefore squares the error-induced rate. This is the
        # only aggregation guaranteed to be available at one-in-a-million abundance, where
        # a rare organism contributes about one fragment and read-to-read grouping has
        # nothing to group (results/readlevel_features/REPORT.md).
        #
        # It bounds the error term alone. Absence that is a property of the ORGANISM — an
        # undersequenced known strain, a low-complexity protein, a non-coding fragment —
        # appears in both mates together and does not square. That correlated term is
        # unmeasured and is what step two has to find on real background.
        entry["mate_conjunction"] = {}
        for budget in budgets:
            t_needed, achieved = None, None
            for T in range(0, n_pep + 2):
                fp = sum(1 for c in counts if c >= T) / n
                if (fp ** 2) * 1e6 <= budget:
                    t_needed, achieved = T, (fp ** 2) * 1e6
                    break
            entry["mate_conjunction"][str(budget)] = {
                "min_threshold": t_needed,
                "headroom_kmers": (n_pep - t_needed) if t_needed is not None else None,
                "achieved_fp_per_million_pairs": achieved,
                "per_read_rate_per_million": (
                    (sum(1 for c in counts if c >= t_needed) / n) * 1e6
                    if t_needed is not None else None),
            }

        rows.append(entry)
        # validation: the fast sampler must reproduce the mean the direct simulation saw
        for d in direct:
            if d["profile"] == "flat" and d["error_rate"] == rate:
                checks.append({"error_rate": rate,
                               "direct_mean_fraction": d["peptide_corrupted_fraction"],
                               "sampler_mean_fraction": entry["mean_corrupted_fraction"]})
    return {"rows": rows, "validation_against_direct_simulation": checks}


def section_headroom(corruption_rows: list[dict], budgets: list[float]) -> list[dict]:
    """Turn the corruption distribution into the number that decides the direction.

    A read from a known organism is flagged as novel if at least T of its peptide k-mers
    come back absent. Under the idealisation stated in the module docstring, the absent
    ones are exactly the corrupted ones, so the error-induced false-flag rate at threshold
    T is P(corrupted >= T), reported per million. Headroom is how many k-mers remain above
    the threshold for a genuinely novel read to spend — if the budget forces T up to the
    full k-mer count, the feature has nothing left to say.
    """
    out = []
    for row in corruption_rows:
        counts = row["_pep_counts"]
        if not counts:
            continue
        n = len(counts)
        # The conservative count: a read at frame offset 1 or 2 carries one fewer
        # k-mer than one at offset 0, and headroom should be quoted for the worst case.
        total_kmers = row["kmers_per_read_min"]
        entry = {"error_rate": row["error_rate"], "profile": row["profile"],
                 "kmers_per_read": total_kmers, "thresholds": {}}
        for budget in budgets:
            t_needed = None
            for T in range(0, (total_kmers or 40) + 2):
                fp = sum(1 for c in counts if c >= T) / n
                if fp * 1e6 <= budget:
                    t_needed = T
                    break
            entry["thresholds"][str(budget)] = {
                "min_threshold": t_needed,
                "headroom_kmers": (total_kmers - t_needed) if (t_needed is not None and total_kmers) else None,
            }
        out.append(entry)
    return out


def _mean(xs) -> float:
    return sum(xs) / len(xs) if xs else 0.0


def _spread(xs) -> dict:
    if not xs:
        return {"mean": 0.0, "min": 0.0, "max": 0.0}
    return {"mean": _mean(xs), "min": min(xs), "max": max(xs)}


# ---------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=5,
                    help="repo rule: a headline number runs across at least 5 seeds")
    ap.add_argument("--read-len", type=int, default=107,
                    help="107 nt is the SRR23973193 read length")
    ap.add_argument("--k-aa", type=int, default=10,
                    help="k_aa=10 is the published specificity floor; below ~7 most "
                         "peptide k-mers are shared across phyla")
    ap.add_argument("--n-reads", type=int, default=20000)
    ap.add_argument("--error-rates", default="0.0005,0.001,0.002,0.005,0.01")
    ap.add_argument("--indel-rates", default="0,1e-5,1e-4,1e-3,0.01,0.05")
    ap.add_argument("--budgets", default="93478,1489,100,1")
    ap.add_argument("--tail-samples", type=int, default=20000000,
                    help="must exceed the tightest budget's reciprocal by a margin or "
                         "the threshold it reports is a zero count, not a measurement")
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "peptide_noise")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "peptide_absence_noise")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)

    # ---- sequence
    cdss, per_genome = [], {}
    for name, acc in GENOMES.items():
        gb = fetch_gb(acc, args.data / f"{acc}.gb")
        seq, feats = parse_gb(gb)
        got = clean_cds(seq, feats)
        per_genome[name] = {"accession": acc, "n_cds_usable": len(got),
                            "total_cds_nt": sum(len(c) for c in got)}
        cdss.extend(got)
    if not cdss:
        sys.exit("no usable CDS parsed; refusing to report numbers from nothing")

    # ---- invariant: a read drawn inside a CDS must be stop-free in its true frame.
    # The same check caught two frame bugs in the previous read-level experiment, so it
    # runs before any measurement rather than after.
    rng = random.Random(args.seed)
    bad = 0
    for clean, off in sample_cds_reads(cdss, args.read_len, 2000, rng):
        if "*" in translate(clean[off:])[:-1]:
            bad += 1
    invariant = {"reads_checked": 2000, "stop_in_true_frame": bad,
                 "passed": bad == 0}
    if not invariant["passed"]:
        sys.exit(f"CDS invariant failed: {bad}/2000 reads have a stop in the true frame")

    rates = [float(x) for x in args.error_rates.split(",")]
    indel_rates = [float(x) for x in args.indel_rates.split(",")]
    budgets = [float(x) for x in args.budgets.split(",")]

    per_seed = []
    for s in range(args.n_seeds):
        rng = random.Random(args.seed + s)
        syn = section_synonymous(cdss, rng, 200000)
        corr = section_corruption(cdss, args.read_len, args.k_aa, rates,
                                  args.n_reads, rng)
        ind = section_indels(cdss, args.read_len, args.k_aa, 0.001, indel_rates,
                             args.n_reads, rng)
        head = section_headroom(corr, budgets)
        tail = section_tail(args.read_len, args.k_aa, rates,
                            1.0 - syn["synonymous_fraction"], budgets,
                            args.tail_samples, rng, corr)
        for row in corr:
            row.pop("_pep_counts", None)
        per_seed.append({"seed": args.seed + s, "synonymous": syn,
                         "corruption": corr, "indels": ind, "headroom": head,
                         "tail": tail})

    # ---- aggregate across seeds, reporting spread as the repo requires
    agg = {"synonymous_fraction": _spread([p["synonymous"]["synonymous_fraction"]
                                           for p in per_seed])}
    agg["corruption"] = []
    for i in range(len(per_seed[0]["corruption"])):
        ref = per_seed[0]["corruption"][i]
        agg["corruption"].append({
            "error_rate": ref["error_rate"], "profile": ref["profile"],
            "peptide_corrupted_fraction": _spread([p["corruption"][i]["peptide_corrupted_fraction"] for p in per_seed]),
            "nucleotide_corrupted_fraction": _spread([p["corruption"][i]["nucleotide_corrupted_fraction"] for p in per_seed]),
            "peptide_minus_nucleotide": _spread([p["corruption"][i]["peptide_minus_nucleotide"] for p in per_seed]),
            "reads_with_zero_corrupted": _spread([p["corruption"][i]["reads_with_zero_corrupted"] for p in per_seed]),
        })
    agg["indels"] = []
    for i in range(len(per_seed[0]["indels"])):
        ref = per_seed[0]["indels"][i]
        agg["indels"].append({
            "indel_rate": ref["indel_rate"],
            "peptide_corrupted_fraction": _spread([p["indels"][i]["peptide_corrupted_fraction"] for p in per_seed]),
            "nucleotide_corrupted_fraction": _spread([p["indels"][i]["nucleotide_corrupted_fraction"] for p in per_seed]),
        })

    payload = {
        "parameters": vars(args) | {"data": str(args.data), "out": str(args.out)},
        "genomes": per_genome,
        "n_cds": len(cdss),
        "cds_invariant": invariant,
        "aggregate": agg,
        "per_seed": per_seed,
    }
    (args.out / "peptide_absence_noise.json").write_text(json.dumps(payload, indent=2))

    manifest = {
        "script": "experiments/peptide_absence_noise.py",
        "seed": args.seed, "n_seeds": args.n_seeds,
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(agg, indent=2)[:4000])


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
