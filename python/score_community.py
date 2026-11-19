"""Score an assembly of the mock metagenome, with novelty as the headline.

Three things are asked of each assembly, and the novelty question is the one
the other two exist to qualify:

1. **Per-genome recovery.** For every truth genome, the fraction covered by
   its best-matching contigs and the identity of that match. A genome counts
   as recovered at (cov, ident) when contigs assigned to it cover that much of
   it at that identity.

2. **Collapse.** For each pair (a known virus and its divergent relative),
   whether both members are recovered separately, or the relative's reads
   were absorbed into the known member's contigs. Collapse is measured by
   assignment: a relative is *collapsed* when it has no contig of its own
   whose best match is the relative rather than its neighbour, while the
   neighbour is recovered. This is how a novel virus disappears, and it is
   reported against the pair's divergence so the result is a threshold, not a
   verdict.

3. **False joins.** A contig whose 500 nt windows match different truth
   genomes is a chimaera; between a virus and a bacterium, or between two
   unrelated viruses, it is an invented sequence. Reported per assembly and
   with the genomes involved, because a chimaera between a known virus and
   its relative is the specific failure that produces a fake "novel" genome.

Contig-to-genome assignment uses minimap2 (via mappy) rather than edlib:
whole bacterial chromosomes make an all-against-all edit distance pointless,
and a metagenome assembly's contigs need a real aligner's placement anyway.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import mappy


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def parse_truth(truth_fa: Path):
    meta, seqs = {}, {}
    for header, s in read_fasta(truth_fa).items():
        parts = header.split()
        name = parts[0]
        seqs[name] = s
        d = {"length": len(s)}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                d[k] = None if v == "None" else v
        meta[name] = d
    return meta, seqs


def score(truth_fa: Path, contigs_fa: Path, min_len=500, win=500, cov_thr=0.9, ident_thr=0.95):
    meta, truth = parse_truth(truth_fa)
    aln = mappy.Aligner(str(truth_fa), preset="asm10", n_threads=4)
    if not aln:
        raise RuntimeError("failed to index truth")
    contigs = {k.split()[0]: v for k, v in read_fasta(contigs_fa).items() if len(v) >= min_len}

    covered = defaultdict(set)          # genome -> covered truth positions (block starts)
    ident_hits = defaultdict(list)      # genome -> (matching bases, aligned bases)
    assigned = {}                       # contig -> genome with the longest aligned block
    per_contig, chimaeras = [], []

    for cname, c in contigs.items():
        hits = [h for h in aln.map(c) if h.is_primary]
        if not hits:
            per_contig.append(dict(contig=cname, length=len(c), assigned=None))
            continue
        best = max(hits, key=lambda h: h.mlen)
        assigned[cname] = best.ctg
        for h in hits:
            g = h.ctg
            ident_hits[g].append((h.mlen, h.blen))
            for p in range(h.r_st, h.r_en, 100):
                covered[g].add(p // 100)

        # windowed genome assignment, for false joins
        wins = []
        for i in range(0, max(len(c) - win + 1, 1), win):
            w = c[i : i + win]
            if len(w) < win // 2:
                break
            wh = [x for x in aln.map(w)]
            wins.append(max(wh, key=lambda x: x.mlen).ctg if wh else None)
        genomes = [g for g in wins if g]
        switch = len(set(genomes)) > 1
        if switch:
            hosts = {meta.get(g, {}).get("host", g) for g in genomes}
            chimaeras.append(dict(contig=cname, length=len(c), genomes=sorted(set(genomes)),
                                  cross_host=len(hosts) > 1))
        per_contig.append(dict(contig=cname, length=len(c), assigned=best.ctg,
                               identity=round(best.mlen / max(best.blen, 1), 5),
                               aligned_bp=best.blen, windows=len(wins),
                               window_genomes=sorted(set(genomes)), chimaeric=switch))

    per_genome = {}
    for g, d in meta.items():
        tl = d["length"]
        frac = min(1.0, 100 * len(covered[g]) / tl) if tl else 0.0
        mlen = sum(x[0] for x in ident_hits[g])
        blen = sum(x[1] for x in ident_hits[g])
        own = [c for c, gg in assigned.items() if gg == g]
        per_genome[g] = dict(role=d.get("role"), relative=d.get("relative"), length=tl,
                             coverage_sim=float(d.get("coverage", 0) or 0),
                             covered_frac=round(frac, 4),
                             identity=round(mlen / blen, 5) if blen else None,
                             own_contigs=len(own),
                             recovered=bool(frac >= cov_thr and blen and mlen / blen >= ident_thr and own))

    pairs = []
    for g, d in meta.items():
        rel = d.get("relative")
        if d.get("role") != "novel" or not rel or rel not in meta:
            continue
        pairs.append(dict(novel=g, known=rel,
                          novel_recovered=per_genome[g]["recovered"],
                          novel_own_contigs=per_genome[g]["own_contigs"],
                          novel_covered=per_genome[g]["covered_frac"],
                          known_recovered=per_genome[rel]["recovered"],
                          collapsed=bool(per_genome[rel]["recovered"] and per_genome[g]["own_contigs"] == 0)))

    return dict(
        n_truth=len(meta), n_contigs=len(contigs),
        total_contig_bp=sum(len(c) for c in contigs.values()),
        unplaced_contigs=sum(1 for p in per_contig if p.get("assigned") is None),
        genomes_recovered=sum(1 for v in per_genome.values() if v["recovered"]),
        viral_genomes_recovered=sum(1 for k, v in per_genome.items()
                                    if v["recovered"] and v["role"] in ("known", "novel")),
        novel_recovered=sum(1 for p in pairs if p["novel_recovered"]),
        novel_collapsed=sum(1 for p in pairs if p["collapsed"]),
        chimaeras=len(chimaeras), cross_host_chimaeras=sum(1 for c in chimaeras if c["cross_host"]),
        pairs=pairs, per_genome=per_genome, chimaera_detail=chimaeras, per_contig=per_contig,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("truth", type=Path)
    ap.add_argument("contigs", type=Path)
    ap.add_argument("--min-len", type=int, default=500)
    ap.add_argument("--cov", type=float, default=0.9)
    ap.add_argument("--ident", type=float, default=0.95)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    r = score(a.truth, a.contigs, a.min_len, ident_thr=a.ident, cov_thr=a.cov)
    if a.json:
        a.json.write_text(json.dumps(r, indent=2))
    for k in ("per_genome", "chimaera_detail", "per_contig", "pairs"):
        r.pop(k)
    print(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()
