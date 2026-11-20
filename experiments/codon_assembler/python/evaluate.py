"""Score an assembly of a simulated swarm world against its ground truth.

Two evaluators, because they answer different questions:

* `membership` uses hapax's layout.tsv (which read went into which contig) and
  the simulator's truth.tsv. It is exact, but only hapax reports membership.
* `kmer` assigns every 200-nt window of every contig to the entity whose
  variants share the most 13-mers with it. It needs nothing but contigs, so it
  scores MEGAHIT, metaSPAdes, PLASS, PenguiN and hapax on equal terms; it is
  the evaluator used for any cross-assembler comparison.

Definitions (shared by both):

* **contigs per true entity**: contigs of at least `min_len` assigned to the
  entity. One population should give one contig, so 1.0 is the target; 0
  means the entity was missed.
* **chimaera**: a contig whose minority entity accounts for at least 10% of
  its reads (membership) or windows (kmer). Entities are separate populations
  that never co-occurred, so any such join is wrong. The sister-entity design
  of the simulator is what makes this number non-trivial.
* **recovered length**: fraction of the entity's root genome covered by the
  single largest contig assigned to it (membership: union of root intervals
  of its reads; kmer: contig length / root length, capped at 1).
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

K = 13
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    name, seq = None, []
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def load_truth(world):
    truth = {}
    with open(world / "truth.tsv") as f:
        next(f)
        for line in f:
            r, e, v, strand, s, t, rs, re_ = line.rstrip("\n").split("\t")
            truth[r] = (int(e), int(rs), int(re_))
    meta = json.loads((world / "world.json").read_text())
    root_len = {int(e): d["root_len"] for e, d in meta["entities"].items()}
    return truth, root_len


def _summary(per_contig, root_len, min_len):
    """per_contig: list of (length, entity_counts Counter, recovered fraction or None)."""
    entities = sorted(root_len)
    n_by_e = Counter()
    best_rec = defaultdict(float)
    chim = 0
    chim_bases = 0
    total = 0
    tot_bases = 0
    for length, cnt, rec in per_contig:
        if length < min_len or not cnt:
            continue
        total += 1
        tot_bases += length
        (e, n), *rest = cnt.most_common()
        n_all = sum(cnt.values())
        if rest and rest[0][1] / n_all >= 0.10:
            chim += 1
            chim_bases += length
        n_by_e[e] += 1
        r = rec if rec is not None else min(1.0, length / root_len[e])
        best_rec[e] = max(best_rec[e], r)
    cpe = [n_by_e[e] for e in entities]
    found = [c for c in cpe if c > 0]
    return {
        "contigs": total,
        "contigs_per_entity": cpe,
        "mean_contigs_per_found_entity": (sum(found) / len(found)) if found else 0.0,
        "entities_missed": sum(1 for c in cpe if c == 0),
        "chimeric_contigs": chim,
        "chimaera_rate": chim / total if total else 0.0,
        "chimeric_base_fraction": chim_bases / tot_bases if tot_bases else 0.0,
        "recovered": [round(best_rec[e], 4) for e in entities],
        "mean_recovered": sum(best_rec[e] for e in entities) / len(entities),
    }


def membership(world: Path, asm: Path, min_len=500):
    truth, root_len = load_truth(world)
    lens = {n.split()[0]: len(s) for n, s in read_fasta(asm / "contigs.fa")}
    members = defaultdict(list)
    with open(asm / "layout.tsv") as f:
        next(f)
        for line in f:
            c, r, *_ = line.split("\t")
            members[c].append(r)
    per = []
    for c, rs in members.items():
        cnt = Counter(truth[r][0] for r in rs)
        e = cnt.most_common(1)[0][0]
        iv = sorted((truth[r][1], truth[r][2]) for r in rs if truth[r][0] == e and truth[r][1] >= 0)
        cov, cur_s, cur_e = 0, None, None
        for s, t in iv:
            if cur_e is None or s > cur_e:
                if cur_e is not None:
                    cov += cur_e - cur_s
                cur_s, cur_e = s, t
            else:
                cur_e = max(cur_e, t)
        if cur_e is not None:
            cov += cur_e - cur_s
        per.append((lens.get(c, 0), cnt, cov / root_len[e]))
    return _summary(per, root_len, min_len)


def _kmers(s):
    s2 = s.translate(COMP)[::-1]
    out = set()
    for x in (s, s2):
        for i in range(len(x) - K + 1):
            out.add(x[i : i + K])
    return out


def kmer(world: Path, contigs: Path, min_len=500, win=200, step=100):
    _, root_len = load_truth(world)
    owner = defaultdict(set)  # kmer -> entities
    for name, s in read_fasta(world / "genomes.fa"):
        e = int(name.split("_")[0][1:])
        for k in _kmers(s):
            owner[k].add(e)
    per = []
    for name, s in read_fasta(contigs):
        cnt = Counter()
        for w0 in range(0, max(1, len(s) - win + 1), step):
            votes = Counter()
            w = s[w0 : w0 + win]
            for i in range(len(w) - K + 1):
                for e in owner.get(w[i : i + K], ()):
                    votes[e] += 1
            if votes:
                (e, n), *_ = votes.most_common()
                if n >= 5:
                    cnt[e] += 1
        per.append((len(s), cnt, None))
    return _summary(per, root_len, min_len)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("world", type=Path)
    ap.add_argument("asm", type=Path, help="hapax output dir, or a contigs FASTA for --kmer")
    ap.add_argument("--kmer", action="store_true")
    ap.add_argument("--min-len", type=int, default=500)
    a = ap.parse_args()
    if a.kmer:
        fa = a.asm / "contigs.fa" if a.asm.is_dir() else a.asm
        print(json.dumps(kmer(a.world, fa, a.min_len), indent=2))
    else:
        print(json.dumps(membership(a.world, a.asm, a.min_len), indent=2))


if __name__ == "__main__":
    main()
