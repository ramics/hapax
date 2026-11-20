"""Synthetic clade generator for leave-one-clade-out novelty evaluation.

Produces coding-like genomes grouped into families. Each family has its own
codon-usage signature, so family membership leaves a real compositional trace
(which is what alignment-free novelty scorers rely on). Members within a family
are related by substitution and indel events applied under a weak coding
constraint, so ORF structure survives -- as it does under selection in real
genomes.

Also produces "junk" controls matched to a genome at increasing orders of
composition. The strongest control is `markov3`, which reproduces the source
genome's 4-mer statistics while destroying its reading-frame structure. Any
novelty score built purely on composition is blind to the difference between a
real divergent genome and its markov3 control; that blindness is the thing the
evaluation is designed to expose.
"""

from __future__ import annotations

import numpy as np

BASES = np.array(list("ACGT"))
B2I = {b: i for i, b in enumerate("ACGT")}
STOPS = {"TAA", "TAG", "TGA"}

SENSE_CODONS = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in STOPS
]
assert len(SENSE_CODONS) == 61


def _codon_usage(rng: np.random.Generator, concentration: float) -> np.ndarray:
    """Draw a family-specific codon usage vector.

    Low concentration gives a strongly skewed, highly distinctive usage;
    high concentration gives usage close to uniform and so families that are
    compositionally hard to tell apart.
    """
    return rng.dirichlet(np.full(61, concentration))


def _make_gene(rng: np.random.Generator, usage: np.ndarray, n_codons: int) -> str:
    body = rng.choice(SENSE_CODONS, size=n_codons, p=usage)
    stop = rng.choice(sorted(STOPS))
    return "ATG" + "".join(body) + stop


def family_root(
    rng: np.random.Generator,
    n_genes: int = 4,
    gene_codons: tuple[int, int] = (250, 600),
    concentration: float = 0.35,
) -> tuple[str, np.ndarray]:
    """Build the ancestral genome of one family, plus its codon usage."""
    usage = _codon_usage(rng, concentration)
    genes = []
    for _ in range(n_genes):
        n_cod = int(rng.integers(*gene_codons))
        genes.append(_make_gene(rng, usage, n_cod))
        # short intergenic spacer, composition drawn from the same usage so the
        # family signature is genome-wide rather than confined to genes
        spacer_len = int(rng.integers(20, 80))
        spacer = rng.choice(BASES, size=spacer_len)
        genes.append("".join(spacer))
    return "".join(genes), usage


def _would_make_stop(seq_list: list[str], pos: int, new_base: str) -> bool:
    """Check whether substituting `new_base` at `pos` creates an in-frame stop.

    Frame is taken relative to the start of the sequence. This is a deliberately
    crude stand-in for purifying selection: it is enough to keep long ORFs
    intact across evolutionary distance, which is all the coherence scorers
    need in order to be testable.
    """
    frame_start = (pos // 3) * 3
    if frame_start + 3 > len(seq_list):
        return False
    codon = list(seq_list[frame_start : frame_start + 3])
    codon[pos - frame_start] = new_base
    return "".join(codon) in STOPS


def evolve(
    rng: np.random.Generator,
    root: str,
    sub_rate: float,
    indel_rate: float = 2e-4,
) -> str:
    """Apply substitutions and in-frame indels to a root genome."""
    seq = list(root)
    n = len(seq)

    n_subs = rng.binomial(n, sub_rate)
    positions = rng.choice(n, size=n_subs, replace=False) if n_subs else []
    for pos in positions:
        current = seq[pos]
        options = [b for b in "ACGT" if b != current]
        rng.shuffle(options)
        for new_base in options:
            if not _would_make_stop(seq, pos, new_base):
                seq[pos] = new_base
                break

    # indels in multiples of three, so reading frames downstream are preserved
    n_indels = rng.binomial(n // 3, indel_rate)
    for _ in range(n_indels):
        if len(seq) < 400:
            break
        codon_pos = int(rng.integers(0, len(seq) // 3)) * 3
        if rng.random() < 0.5:
            del seq[codon_pos : codon_pos + 3]
        else:
            ins = list(rng.choice(BASES, size=3))
            seq[codon_pos:codon_pos] = ins

    return "".join(seq)


def build_families(
    n_families: int = 6,
    members_per_family: int = 8,
    seed: int = 0,
    concentration: float = 0.35,
    sub_rate_range: tuple[float, float] = (0.02, 0.22),
) -> dict:
    """Build the full synthetic world.

    Returns a dict with `genomes` (list of (family_id, member_id, sequence))
    and `roots` (family_id -> ancestral sequence).
    """
    rng = np.random.default_rng(seed)
    genomes, roots = [], {}

    for fam in range(n_families):
        root, _usage = family_root(rng, concentration=concentration)
        roots[fam] = root
        # graded divergence: member 0 sits close to the root, later members far
        rates = np.linspace(*sub_rate_range, members_per_family)
        for member, rate in enumerate(rates):
            genomes.append((fam, member, evolve(rng, root, float(rate))))

    return {"genomes": genomes, "roots": roots}


# ---------------------------------------------------------------- junk controls


def shuffle_mono(rng: np.random.Generator, seq: str) -> str:
    arr = np.array(list(seq))
    rng.shuffle(arr)
    return "".join(arr)


def random_matched_gc(rng: np.random.Generator, seq: str) -> str:
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
    return "".join(rng.choice(BASES, size=len(seq), p=p))


def markov_surrogate(rng: np.random.Generator, seq: str, order: int = 3) -> str:
    """Generate a sequence reproducing `seq`'s (order+1)-mer statistics.

    This is the control that matters. It matches the source genome on exactly
    the features that alignment-free composition scorers read, so a scorer that
    only looks at composition cannot separate the two. What it does not
    reproduce is long-range structure: reading frames, ORF lengths, consistent
    codon usage across a whole gene.
    """
    idx = np.array([B2I[c] for c in seq if c in B2I], dtype=np.int64)
    n = len(idx)
    if n <= order + 1:
        return seq

    # transition counts over contexts of length `order`
    n_ctx = 4**order
    counts = np.ones((n_ctx, 4), dtype=np.float64)  # Laplace smoothing
    powers = 4 ** np.arange(order - 1, -1, -1)
    ctx_ids = (idx[np.arange(n - order)[:, None] + np.arange(order)] * powers).sum(1)
    np.add.at(counts, (ctx_ids, idx[order:]), 1.0)
    probs = counts / counts.sum(1, keepdims=True)

    # sample via precomputed cumulatives and one batch of uniforms: the chain is
    # sequential, but per-step rng calls dominate the cost otherwise
    cum = probs.cumsum(1)
    cum[:, -1] = 1.0
    uniforms = rng.random(n - order)

    out = list(idx[:order])
    ctx = int((idx[:order] * powers).sum())
    for u in uniforms:
        row = cum[ctx]
        nxt = 0 if u < row[0] else 1 if u < row[1] else 2 if u < row[2] else 3
        out.append(nxt)
        ctx = (ctx * 4 + nxt) % n_ctx
    return "".join(BASES[out])


def codon_permutation_surrogate(rng: np.random.Generator, seq: str) -> str:
    """Permute codons within frame 0, preserving codon usage exactly.

    Complementary to `markov_surrogate`, and hard in a different place. Codon
    usage in frame 0 is untouched, so any statistic reading codon-usage
    concentration sees the real genome's signature. What moves is arrangement:
    stop codons that sat in intergenic regions are scattered uniformly, so long
    open spans disappear.

    A statistic that beats this null is reading gene architecture rather than
    codon composition.
    """
    n_codons = len(seq) // 3
    codons = [seq[i * 3 : i * 3 + 3] for i in range(n_codons)]
    rng.shuffle(codons)
    return "".join(codons) + seq[n_codons * 3 :]


def fragment(rng: np.random.Generator, seq: str, length: int) -> str:
    """Draw one fragment of the requested length from a random position.

    Used by the Phase 0 length sweep. Fragments shorter than the genome are
    drawn uniformly; a request at or above genome length returns the genome.
    """
    if length >= len(seq):
        return seq
    start = int(rng.integers(0, len(seq) - length))
    return seq[start : start + length]


JUNK_MAKERS = {
    "random_gc": random_matched_gc,
    "shuffle_mono": shuffle_mono,
    "markov3": lambda rng, seq: markov_surrogate(rng, seq, order=3),
    "markov5": lambda rng, seq: markov_surrogate(rng, seq, order=5),
    "codon_perm": codon_permutation_surrogate,
}
