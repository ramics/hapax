import numpy as np

from python import swarm


def _translate(s):
    return "".join(swarm.CODON_AA[s[i : i + 3]] for i in range(0, len(s) - 2, 3))


def test_root_has_dual_region_open_in_both_frames():
    rng = np.random.default_rng(0)
    g = swarm.build_root(rng)
    seq = "".join(g.seq)
    dual = [p for p in range(len(seq)) if g.phase[1][p] >= 0]
    assert len(dual) > 200
    s = dual[0]
    assert "*" not in _translate(seq[s : dual[-1] + 1])
    assert sum(1 for p in range(len(seq)) if g.phase[1][p] >= 0 and g.phase[0][p] >= 0) == len(dual)


def test_selection_keeps_genes_open_and_mostly_synonymous():
    rng = np.random.default_rng(1)
    root = swarm.build_root(rng)
    child = swarm.mutate(rng, root, 300, omega=0.23, indel_per_sub=0.0)
    seq = "".join(child.seq)
    starts = [p for p in range(len(seq)) if child.phase[0][p] == 0]
    # no internal stops in track-0 codons except each gene's terminal stop
    stops = [p for p in starts if p + 3 <= len(seq) and seq[p : p + 3] in ("TAA", "TAG", "TGA")]
    assert len(stops) == 3
    syn = non = 0
    for p in starts:
        a, b = "".join(root.seq[p : p + 3]), seq[p : p + 3]
        if a != b:
            if swarm.CODON_AA[a] == swarm.CODON_AA[b]:
                syn += 1
            else:
                non += 1
    assert syn > non


def test_swarm_diversity_tracks_target():
    rng = np.random.default_rng(2)
    root = swarm.build_root(rng)
    vs = swarm.evolve_swarm(rng, root, 0.10, 8, 0.23)
    # indel-free comparison on coding track via orig coordinates
    d = []
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            a = {o: b for o, b in zip(vs[i].orig, vs[i].seq) if o >= 0}
            b = {o: c for o, c in zip(vs[j].orig, vs[j].seq) if o >= 0}
            common = set(a) & set(b)
            d.append(sum(a[o] != b[o] for o in common) / len(common))
    m = float(np.mean(d))
    assert 0.06 < m < 0.12, m


def test_world_truth_consistent():
    w = swarm.build_world(3, 0.05, coverage=5.0)
    assert len(w.reads) == len(w.truth)
    assert {t[1] for t in w.truth} == {0, 1, 2, 3}
