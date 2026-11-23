"""Read simulators for the incumbent stress test: paired-end Illumina and nanopore.

`swarm.simulate_reads` makes single-end, fixed-length short reads, which is
what hapax itself takes. The stress test also has to feed tools that need
paired-end input (VStrains, and SAVAGE/Haploflow do better with it) and the
long-read haplotype assemblers (Strainline, RVHaplo, HaploDMF, viralFlye), so
this module adds both read types on plain nucleotide strings.

Why a hand-rolled nanopore model rather than Badread or NanoSim: the point of
the long-read arm is to see how the tools cope with indel error that lands in
homopolymers and so can shift a reading frame. That needs the error mix to be
explicit and settable, and the per-base truth (where each error went) to be
recorded so a "wrong frameshift" in a contig can be traced to read error
rather than a real codon indel in the swarm. Calibration, chosen to match
published per-read accuracy rather than tuned:

* **R9.4.1** (Guppy/Dorado HAC era): modal read identity ~94%, i.e. ~6% error.
  Errors split roughly 35% substitution, 40% deletion, 25% insertion
  (Delahaye & Nicolas 2021 PLoS ONE; Wick et al. 2019 Genome Biol), deletions
  concentrated in homopolymers.
* **R10.4.1** (Dorado SUP): modal identity ~98-99%, ~2% error, with a smaller
  but still homopolymer-dominated indel share (Sereika et al. 2022 Nat
  Methods; Ni et al. 2023). R10's dual reader head helps homopolymers but runs
  of 5+ still lose or gain bases.
* **Homopolymer length error.** For a run of length h >= 3, the called length
  is h + delta with P(delta != 0) = hp_base * (h - 2), capped at 0.6, and
  delta = -1 three times as often as +1 (under-calling dominates). This term
  is applied on top of the background indel rate.
* **Per-read identity** varies: each read draws its error rate from a gamma
  around the mean (shape 8), because real runs have a long tail of poor reads
  that the tools must survive.

HIV env long-read data is almost always amplicon sequencing, so reads default
to full-length amplicons with a fraction truncated (broken molecules), rather
than a random fragment-length distribution.
"""

from __future__ import annotations

import numpy as np

BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

NANOPORE = {
    # mean error, fractions (sub, del, ins) of background error, homopolymer base rate
    "r9": dict(err=0.06, sub=0.35, dele=0.40, ins=0.25, hp_base=0.08),
    "r10": dict(err=0.02, sub=0.40, dele=0.35, ins=0.25, hp_base=0.03),
}


def revcomp(s: str) -> str:
    return "".join(COMP[b] for b in reversed(s))


def _illumina_errors(rng, frag: str, sub_err: float, indel_err: float, qramp: float = 3.0) -> str:
    """Substitution error rising along the read (Illumina 3' quality decay)."""
    n = len(frag)
    out = []
    for i, b in enumerate(frag):
        r = rng.random()
        if r < indel_err / 2:
            continue
        if r < indel_err:
            out.append(BASES[int(rng.integers(4))])
        p = sub_err * (1.0 + qramp * (i / max(n - 1, 1)) ** 2)
        if rng.random() < p:
            b = BASES[(BASES.index(b) + 1 + int(rng.integers(3))) % 4]
        out.append(b)
    return "".join(out)


def paired_end(rng, seqs, abundances, coverage, read_len=150, insert_mean=350, insert_sd=50, sub_err=0.002, indel_err=5e-5):
    """Fragments sampled across haplotypes by abundance; returns (r1, r2, truth).

    truth rows: (pair_id, haplotype index, fragment start, fragment end, strand).
    """
    L = float(np.mean([len(s) for s in seqs]))
    n_pairs = int(coverage * L / (2 * read_len))
    p = np.asarray(abundances, float)
    p /= p.sum()
    r1, r2, truth = [], [], []
    for k in range(n_pairs):
        h = int(rng.choice(len(seqs), p=p))
        s = seqs[h]
        ins = int(np.clip(rng.normal(insert_mean, insert_sd), read_len, len(s)))
        a = int(rng.integers(0, len(s) - ins + 1))
        frag = s[a : a + ins]
        strand = "+"
        if rng.random() < 0.5:
            frag, strand = revcomp(frag), "-"
        rid = f"p{k}"
        r1.append((rid + "/1", _illumina_errors(rng, frag[:read_len], sub_err, indel_err)))
        r2.append((rid + "/2", _illumina_errors(rng, revcomp(frag)[:read_len], sub_err, indel_err)))
        truth.append((rid, h, a, a + ins, strand))
    return r1, r2, truth


def _runs(s: str):
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        yield s[i], j - i
        i = j


def nanopore_read(rng, frag: str, chem: str):
    """Apply a nanopore error profile; returns (read, n_sub, n_del, n_ins, n_hp)."""
    m = NANOPORE[chem]
    err = rng.gamma(8.0, m["err"] / 8.0)
    # background error excludes the homopolymer term, which is added separately
    ps, pd, pi = err * m["sub"], err * m["dele"], err * m["ins"]
    out = []
    ns = nd = ni = nh = 0
    for base, h in _runs(frag):
        if h >= 3 and rng.random() < min(0.6, m["hp_base"] * (h - 2) * (err / m["err"])):
            h += -1 if rng.random() < 0.75 else 1
            nh += 1
        for _ in range(h):
            r = rng.random()
            if r < pd:
                nd += 1
                continue
            if r < pd + pi:
                out.append(BASES[int(rng.integers(4))])
                ni += 1
            b = base
            if rng.random() < ps:
                b = BASES[(BASES.index(b) + 1 + int(rng.integers(3))) % 4]
                ns += 1
            out.append(b)
    return "".join(out), ns, nd, ni, nh


def nanopore(rng, seqs, abundances, n_reads, chem="r9", truncated_frac=0.2):
    """Amplicon-style nanopore reads; returns (reads, truth).

    truth rows: (read id, haplotype, start, end, strand, subs, dels, ins, hp_errors).
    """
    p = np.asarray(abundances, float)
    p /= p.sum()
    reads, truth = [], []
    for k in range(n_reads):
        h = int(rng.choice(len(seqs), p=p))
        s = seqs[h]
        a, b = 0, len(s)
        if rng.random() < truncated_frac:
            span = int(rng.uniform(0.3, 0.9) * len(s))
            a = int(rng.integers(0, len(s) - span + 1))
            b = a + span
        frag = s[a:b]
        strand = "+"
        if rng.random() < 0.5:
            frag, strand = revcomp(frag), "-"
        read, ns, nd, ni, nh = nanopore_read(rng, frag, chem)
        rid = f"n{k}"
        reads.append((rid, read))
        truth.append((rid, h, a, b, strand, ns, nd, ni, nh))
    return reads, truth


def write_fastq(path, reads, q="5"):
    """FASTQ with flat quality; several tools refuse FASTA input."""
    with open(path, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{q * len(s)}\n")
