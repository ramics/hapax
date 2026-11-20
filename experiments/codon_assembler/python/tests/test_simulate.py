"""Keep the Markov surrogate sampler under test (CLAUDE.md: it is the fiddly part)."""

from collections import Counter

import numpy as np

from python import simulate as sim


def _kmer_freq(s, k):
    c = Counter(s[i : i + k] for i in range(len(s) - k + 1))
    n = sum(c.values())
    return {x: v / n for x, v in c.items()}


def test_markov3_preserves_4mer_statistics():
    rng = np.random.default_rng(0)
    root, _ = sim.family_root(rng)
    sur = sim.markov_surrogate(np.random.default_rng(1), root, order=3)
    assert len(sur) == len(root)
    a, b = _kmer_freq(root, 4), _kmer_freq(sur, 4)
    keys = set(a) | set(b)
    l1 = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    assert l1 < 0.25, l1


def test_markov3_destroys_reading_frames():
    rng = np.random.default_rng(2)
    root, _ = sim.family_root(rng)
    sur = sim.markov_surrogate(np.random.default_rng(3), root, order=3)

    def longest_open(s):
        best = 0
        for f in range(3):
            run = 0
            for i in range(f, len(s) - 2, 3):
                run = 0 if s[i : i + 3] in sim.STOPS else run + 1
                best = max(best, run)
        return best

    assert longest_open(root) > 3 * longest_open(sur)


def test_codon_permutation_keeps_usage():
    rng = np.random.default_rng(4)
    root, _ = sim.family_root(rng)
    p = sim.codon_permutation_surrogate(np.random.default_rng(5), root)
    n = len(root) // 3
    assert Counter(root[i * 3 : i * 3 + 3] for i in range(n)) == Counter(p[i * 3 : i * 3 + 3] for i in range(n))


def test_evolve_keeps_frames_free_of_new_stops():
    rng = np.random.default_rng(6)
    g = sim._make_gene(rng, sim._codon_usage(rng, 0.35), 300)
    e = sim.evolve(rng, g, 0.1, indel_rate=0.0)
    body = [e[i : i + 3] for i in range(0, len(e) - 3, 3)]
    assert not any(c in sim.STOPS for c in body)
