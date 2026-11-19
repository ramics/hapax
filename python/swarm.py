"""Within-population swarm simulator for the diversity ladder.

`simulate.py` builds separate families for a novelty-scoring experiment. This
module reuses its codon-usage and gene primitives to build what the assembler
actually faces: a few viral *entities*, each sampled as a swarm of related
variants from one population, sequenced as short reads.

Design choices and what they are calibrated against:

* **Population structure.** Each entity is a random 32-tip tree (the author's
  thesis used a random 32-taxon guide tree with EvolveAGene), scaled so that
  the mean pairwise distance between tips equals the target diversity
  (0.05-0.20 substitutions per site, the thesis ladder).
* **Selection.** Substitutions are proposed at random and accepted under a
  codon-level filter: synonymous always, non-synonymous with probability
  omega (default 0.23, intra-patient HIV-1 env backbone dN/dS, Lemey et al.
  2007), stops never. 10% of codons are under positive selection (omega 1.5)
  and 10% under strong purifying selection (omega 0.05), as in the thesis
  EvolveAGene runs. EvolveAGene itself is not used: it cannot be fetched in
  this sandbox and its site model is not settable, so the selection filter
  here makes the parameters explicit instead.
* **Dual coding.** Each genome carries a region where a second gene overlaps
  the first in the +1 frame (as tat/rev overlap env in HIV-1). A substitution
  there must pass the filter in *both* frames, which is what makes such
  regions conserved and is the regime METHOD.md section 5 worries about.
* **Untranslated sequence.** 5' and 3' UTRs and intergenic spacers evolve
  neutrally and take single-base indels.
* **Sister entities.** Optionally one entity has a sister population whose
  root diverged from it by `sister_div`, so that "same entity" and "related
  but different entity" overlaps can be confused. Joins between them are the
  chimaeras that matter. Unrelated families make chimaeras easy to avoid and
  would flatter every assembler.
* **Reads.** Single-end, fixed length, both strands, abundance across tips
  log-normal (swarms are skewed), substitution and indel error at modern
  Illumina rates (0.2% and 0.005% per base) rather than the 2014 platform
  table -- chemistry has moved on since then.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:  # package-relative when imported from tests, plain when run as a script
    from . import simulate as sim
except ImportError:  # pragma: no cover
    import simulate as sim

BASES = "ACGT"
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
CODON_AA = {}
_AA = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
for _i, _c in enumerate(a + b + c for a in BASES for b in BASES for c in BASES):
    CODON_AA[_c] = _AA[_i]


@dataclass
class Genome:
    """A genome plus per-position annotation carried through evolution.

    `phase[t][p]` is the codon phase (0, 1, 2) of position p in coding track t,
    or -1 if p is not coding in that track. Track 0 is the main gene set;
    track 1 holds the overlapping +1-frame gene. `omega[p]` is the site's
    dN/dS multiplier class. `orig[p]` is the root coordinate of p (-1 for
    inserted bases), used to score recovered length against the root.
    """

    seq: list
    phase: list  # two lists of ints
    omega: list
    orig: list
    # relative substitution rate per site (hypervariable loops > 1); None = uniform
    rate: list | None = None
    # per-site codon-indel hotspot flag (1 = indels of 1-4 codons allowed here)
    hot: list | None = None

    def copy(self) -> "Genome":
        return Genome(
            list(self.seq),
            [list(self.phase[0]), list(self.phase[1])],
            list(self.omega),
            list(self.orig),
            None if self.rate is None else list(self.rate),
            None if self.hot is None else list(self.hot),
        )

    def arrays(self):
        """All per-position arrays, for splicing on indels."""
        out = [self.seq, self.phase[0], self.phase[1], self.omega, self.orig]
        if self.rate is not None:
            out.append(self.rate)
        if self.hot is not None:
            out.append(self.hot)
        return out

    def region_of(self, p: int) -> str:
        n = (self.phase[0][p] >= 0) + (self.phase[1][p] >= 0)
        return ["noncoding", "coding", "dual"][n]


def _sample_codons(rng, usage, n, prev_tail=None, avoid_plus1_stop=False):
    """Sample n sense codons; optionally keep the +1 frame free of stops."""
    out = []
    tail = prev_tail
    for _ in range(n):
        while True:
            c = sim.SENSE_CODONS[rng.choice(61, p=usage)]
            if avoid_plus1_stop and tail is not None and (tail + c[0]) in sim.STOPS:
                continue
            break
        out.append(c)
        tail = c[1:]
    return out


def build_root(rng, concentration=0.35, utr=(200, 200), genes=(400, 600, 300), dual_codons=100, spacer=50) -> Genome:
    """Ancestral genome: UTR, three genes, a dual-coding overlap inside gene 2, UTR."""
    usage = sim._codon_usage(rng, concentration)
    seq, ph0, ph1 = [], [], []

    def noncoding(n):
        s = list(rng.choice(list(BASES), size=n))
        seq.extend(s)
        ph0.extend([-1] * n)
        ph1.extend([-1] * n)

    def gene(n_codons, dual=False):
        start = len(seq)
        codons = ["ATG"]
        if dual:
            # dual region sits in the middle third of the gene
            a = n_codons // 3
            codons += _sample_codons(rng, usage, a)
            codons += _sample_codons(rng, usage, dual_codons, prev_tail=codons[-1][1:], avoid_plus1_stop=True)
            codons += _sample_codons(rng, usage, n_codons - a - dual_codons)
            dual_start = start + 3 * (1 + a) + 1  # +1 frame begins one base in
            dual_end = dual_start + 3 * (dual_codons - 1)
        else:
            codons += _sample_codons(rng, usage, n_codons)
        codons.append(str(rng.choice(sorted(sim.STOPS))))
        s = "".join(codons)
        seq.extend(s)
        ph0.extend([i % 3 for i in range(len(s))])
        ph1.extend([-1] * len(s))
        if dual:
            for p in range(dual_start, dual_end):
                ph1[p] = (p - dual_start) % 3

    noncoding(utr[0])
    gene(genes[0])
    noncoding(spacer)
    gene(genes[1], dual=True)
    noncoding(spacer)
    gene(genes[2])
    noncoding(utr[1])

    n = len(seq)
    # site classes per codon of track 0: 10% positive, 10% purifying (thesis setup)
    omega = [1.0] * n
    cls = rng.random(n // 3 + 1)
    for p in range(n):
        if ph0[p] >= 0 or ph1[p] >= 0:
            c = cls[p // 3]
            omega[p] = 1.5 if c < 0.1 else (0.05 / 0.23 if c < 0.2 else 1.0)
    return Genome(seq, [ph0, ph1], omega, list(range(n)))


def _accept(g: Genome, p: int, new: str, rng, omega: float) -> bool:
    for t in (0, 1):
        ph = g.phase[t][p]
        if ph < 0:
            continue
        cs = p - ph
        if cs + 3 > len(g.seq):
            return False
        old = "".join(g.seq[cs : cs + 3])
        cod = list(old)
        cod[ph] = new
        cod = "".join(cod)
        a_old, a_new = CODON_AA.get(old, "X"), CODON_AA.get(cod, "X")
        if a_new == "*":
            return False
        if a_new != a_old and rng.random() >= omega * g.omega[p]:
            return False
    return True


def mutate(rng, g: Genome, n_subs: int, omega: float, kappa: float = 2.0, indel_per_sub: float = 0.01) -> Genome:
    """Apply `n_subs` *accepted* substitutions plus a proportional number of indels."""
    g = g.copy()
    p_ts = kappa / (kappa + 2.0)
    done = 0
    tries = 0
    cum = None
    if g.rate is not None:
        cum = np.cumsum(g.rate)
    while done < n_subs and tries < 200 * (n_subs + 1):
        tries += 1
        p = int(rng.integers(len(g.seq))) if cum is None else int(np.searchsorted(cum, rng.random() * cum[-1], side="right"))
        cur = g.seq[p]
        if rng.random() < p_ts:
            new = TRANSITION[cur]
        else:
            new = str(rng.choice([b for b in BASES if b != cur and b != TRANSITION[cur]]))
        if _accept(g, p, new, rng, omega):
            g.seq[p] = new
            done += 1
    n_indel = rng.poisson(n_subs * indel_per_sub)
    for _ in range(n_indel):
        p = int(rng.integers(len(g.seq) - 3)) if cum is None else int(np.searchsorted(cum, rng.random() * cum[-1], side="right"))
        p = min(p, len(g.seq) - 4)
        _indel_at(rng, g, p)
    return g


def _splice_insert(g: Genome, p: int, bases: list, ph0: list) -> None:
    n = len(bases)
    g.seq[p:p] = bases
    g.phase[0][p:p] = ph0
    g.phase[1][p:p] = [-1] * n
    g.omega[p:p] = [1.0] * n
    g.orig[p:p] = [-1] * n
    if g.rate is not None:
        g.rate[p:p] = [g.rate[p - 1] if p > 0 else 1.0] * n
    if g.hot is not None:
        g.hot[p:p] = [g.hot[p - 1] if p > 0 else 0] * n


def _indel_at(rng, g: Genome, p: int) -> None:
    """One indel at p: frame-preserving in coding sequence, free in UTRs.

    Dual-coding sites take none. Hotspot sites (hypervariable loops) take
    indels of 1-4 codons, as env V1/V2/V4/V5 do; other coding sites take one
    codon at most.
    """
    if g.phase[1][p] >= 0:
        return
    if g.phase[0][p] >= 0:
        p -= g.phase[0][p]  # codon boundary
        if p < 0 or g.phase[1][p] >= 0 or g.phase[0][p] != 0:
            return
        n_cod = int(rng.integers(1, 5)) if (g.hot is not None and g.hot[p]) else 1
        if rng.random() < 0.5:
            span = 3 * n_cod
            if any(g.phase[0][q] < 0 or g.phase[1][q] >= 0 for q in range(p, min(p + span, len(g.seq)))):
                return
            for arr in g.arrays():
                del arr[p : p + span]
        else:
            ins = []
            for _ in range(n_cod):
                ins += list(sim.SENSE_CODONS[int(rng.integers(61))])
            _splice_insert(g, p, ins, [0, 1, 2] * n_cod)
    else:
        if rng.random() < 0.5:
            for arr in g.arrays():
                del arr[p]
        else:
            _splice_insert(g, p, [str(rng.choice(list(BASES)))], [-1])


def random_tree(rng, n_tips: int):
    """Kingman-style random coalescent tree. Returns (parent, branch_len, tips)."""
    nodes = list(range(n_tips))
    parent = {}
    blen = {}
    height = {i: 0.0 for i in nodes}
    t = 0.0
    nxt = n_tips
    while len(nodes) > 1:
        k = len(nodes)
        t += rng.exponential(1.0 / (k * (k - 1) / 2))
        i, j = rng.choice(len(nodes), size=2, replace=False)
        a, b = nodes[i], nodes[j]
        for c in (a, b):
            parent[c] = nxt
            blen[c] = t - height[c]
        height[nxt] = t
        nodes = [x for x in nodes if x not in (a, b)] + [nxt]
        nxt += 1
    return parent, blen, list(range(n_tips)), nodes[0]


def mean_pairwise(parent, blen, tips):
    def path(x):
        out = {}
        d = 0.0
        while x in parent:
            out[x] = d
            d += blen[x]
            x = parent[x]
        out[x] = d
        return out

    paths = {t: path(t) for t in tips}
    tot, n = 0.0, 0
    for i, a in enumerate(tips):
        for b in tips[i + 1 :]:
            pa, pb = paths[a], paths[b]
            # lowest common ancestor = first ancestor of b present in a's path
            x = b
            while x not in pa:
                x = parent[x]
            tot += pa[x] + pb[x]
            n += 1
    return tot / n


def evolve_swarm(rng, root: Genome, diversity: float, n_tips: int, omega: float) -> list:
    parent, blen, tips, top = random_tree(rng, n_tips)
    scale = diversity / mean_pairwise(parent, blen, tips)
    children = {}
    for c, p in parent.items():
        children.setdefault(p, []).append(c)
    genomes = {top: root}
    stack = [top]
    while stack:
        u = stack.pop()
        for c in children.get(u, []):
            n_subs = rng.poisson(len(genomes[u].seq) * blen[c] * scale)
            genomes[c] = mutate(rng, genomes[u], n_subs, omega)
            stack.append(c)
    return [genomes[t] for t in tips]


def simulate_reads(rng, variants, abundances, coverage, read_len, sub_err, indel_err, entity, start_id):
    L = np.mean([len(v.seq) for v in variants])
    n_reads = int(coverage * L / read_len)
    probs = np.asarray(abundances) / np.sum(abundances)
    reads, truth = [], []
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for k in range(n_reads):
        vi = int(rng.choice(len(variants), p=probs))
        v = variants[vi]
        s = int(rng.integers(0, len(v.seq) - read_len + 1))
        frag = v.seq[s : s + read_len]
        out = []
        for b in frag:
            r = rng.random()
            if r < indel_err / 2:
                continue  # deletion
            if r < indel_err:
                out.append(str(rng.choice(list(BASES))))  # insertion before b
            if rng.random() < sub_err:
                b = str(rng.choice([x for x in BASES if x != b]))
            out.append(b)
        strand = "+"
        if rng.random() < 0.5:
            out = [comp[b] for b in reversed(out)]
            strand = "-"
        origs = [o for o in v.orig[s : s + read_len] if o >= 0]
        rid = f"r{start_id + k}"
        reads.append((rid, "".join(out)))
        truth.append((rid, entity, vi, strand, s, s + read_len, min(origs) if origs else -1, max(origs) + 1 if origs else -1))
    return reads, truth


@dataclass
class World:
    roots: dict = field(default_factory=dict)  # entity -> Genome
    variants: dict = field(default_factory=dict)  # entity -> [Genome]
    family: dict = field(default_factory=dict)  # entity -> family id
    reads: list = field(default_factory=list)
    truth: list = field(default_factory=list)


def build_world(
    seed: int,
    diversity: float,
    n_families: int = 3,
    sister_div: float | None = 0.30,
    n_tips: int = 32,
    omega: float = 0.23,
    coverage: float = 50.0,
    read_len: int = 150,
    sub_err: float = 0.002,
    indel_err: float = 5e-5,
    abundance_sigma: float = 1.0,
) -> World:
    rng = np.random.default_rng(seed)
    w = World()
    roots = []
    for f in range(n_families):
        roots.append((f, build_root(rng)))
    if sister_div is not None:
        r0 = roots[0][1]
        sister = mutate(rng, r0, rng.poisson(len(r0.seq) * sister_div), omega)
        sister.orig = list(range(len(sister.seq)))  # own coordinate system
        roots.insert(1, (0, sister))
    rid = 0
    for e, (fam, root) in enumerate(roots):
        w.roots[e] = root
        w.family[e] = fam
        vs = evolve_swarm(rng, root, diversity, n_tips, omega)
        w.variants[e] = vs
        ab = rng.lognormal(0.0, abundance_sigma, size=len(vs))
        reads, truth = simulate_reads(rng, vs, ab, coverage, read_len, sub_err, indel_err, e, rid)
        rid += len(reads)
        w.reads += reads
        w.truth += truth
    # shuffle read order so input order carries no information
    order = rng.permutation(len(w.reads))
    w.reads = [w.reads[i] for i in order]
    w.truth = [w.truth[i] for i in order]
    return w


def write_world(w: World, out: Path, meta: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "reads.fa", "w") as f:
        for rid, s in w.reads:
            f.write(f">{rid}\n{s}\n")
    with open(out / "truth.tsv", "w") as f:
        f.write("read\tentity\tvariant\tstrand\tstart\tend\troot_start\troot_end\n")
        for t in w.truth:
            f.write("\t".join(map(str, t)) + "\n")
    with open(out / "genomes.fa", "w") as f:
        for e, vs in w.variants.items():
            for i, v in enumerate(vs):
                f.write(f">e{e}_v{i}\n{''.join(v.seq)}\n")
    with open(out / "roots.fa", "w") as f:
        for e, r in w.roots.items():
            f.write(f">e{e} family={w.family[e]}\n{''.join(r.seq)}\n")
    with open(out / "annotation.tsv", "w") as f:
        f.write("entity\tstart\tend\tregion\n")
        for e, r in w.roots.items():
            s = 0
            n = len(r.seq)
            while s < n:
                t = s
                while t < n and r.region_of(t) == r.region_of(s):
                    t += 1
                f.write(f"{e}\t{s}\t{t}\t{r.region_of(s)}\n")
                s = t
    (out / "world.json").write_text(json.dumps({**meta, "entities": {str(e): {"family": w.family[e], "root_len": len(w.roots[e].seq)} for e in w.roots}}, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--diversity", type=float, default=0.10)
    ap.add_argument("--families", type=int, default=3)
    ap.add_argument("--sister-div", type=float, default=0.30)
    ap.add_argument("--no-sister", action="store_true")
    ap.add_argument("--coverage", type=float, default=50.0)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    w = build_world(a.seed, a.diversity, a.families, None if a.no_sister else a.sister_div, coverage=a.coverage)
    write_world(w, a.out, vars(a) | {"out": str(a.out)})


if __name__ == "__main__":
    main()
