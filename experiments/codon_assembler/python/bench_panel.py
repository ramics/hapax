"""Build haplotype panels for the incumbent stress test from real sequence sets.

The point of this benchmark is to test published viral haplotype assemblers on
ground that is outside their own papers' benchmarks, so the truth set is real
sequence data rather than anything evolved in silico: within-host single
genome amplification (SGA) env sequence sets from the Zambia transmission
cohort (Derdeyn/Hunter, NCBI nuccore; accessions and dates in
data/MANIFEST.md). Each host's set is a sample of the haplotypes that actually
co-existed in one person, so the real length variation in the hypervariable
loops -- these sets span 2412-2604 nt for the same gene -- is in the truth set
rather than imposed by a model.

Panels:

* `within <host>` -- one host's SGA set. Mean pairwise divergence per host in
  this data runs 0.3% (ZM231F, ZM206F) to 3.0% with a 14% tail (ZM247F), which
  brackets what an intra-host swarm looks like.
* `mixed` -- sets from several hosts pooled, so between-host divergence (~10%+)
  sits on top of within-host swarms, as in a real metagenome. Joins between
  hosts are unambiguous chimaeras; joins between haplotypes of one host are
  the hard case.
* Abundance is log-normal (swarms are skewed) with a settable sigma; sigma 0
  gives the even mixtures the incumbents' own benchmarks use.

Nothing here simulates evolution or designs sequence variants. It samples
published sequences, assigns abundances, and emits reads via `sequencers`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    from . import sequencers as sq
except ImportError:  # pragma: no cover
    import sequencers as sq


def read_fasta(path) -> dict:
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        out[name] = "".join(buf).upper()
    return out


def build_panel(rng, refdir: Path, hosts, n_per_host=None, sigma=1.0):
    """Returns (names, seqs, abundances, host_of_name)."""
    names, seqs, host_of = [], [], {}
    for h in hosts:
        s = read_fasta(refdir / f"{h}.fa")
        # a handful of deposited sequences carry IUPAC ambiguity codes; drop them
        # rather than resolving, so the truth set is unaltered submitted sequence
        s = {k: v for k, v in s.items() if set(v) <= set("ACGT")}
        # identical sequences occur within a host's set (the same haplotype
        # amplified twice). No assembler can be asked to return both, and
        # keeping them would make recovery unreachable by construction, so a
        # panel holds distinct sequences only.
        seen, uniq = set(), {}
        for k in sorted(s):
            if s[k] not in seen:
                seen.add(s[k])
                uniq[k] = s[k]
        s = uniq
        ks = sorted(s)
        if n_per_host:
            ks = [ks[i] for i in rng.choice(len(ks), size=min(n_per_host, len(ks)), replace=False)]
        for k in ks:
            names.append(k)
            seqs.append(s[k])
            host_of[k] = h
    ab = rng.lognormal(0.0, sigma, size=len(seqs)) if sigma > 0 else np.ones(len(seqs))
    return names, seqs, ab / ab.sum(), host_of


def divergence_stats(seqs):
    """Pairwise divergence of the panel; the minimum is what sets the floor.

    An assembler cannot be asked to separate haplotypes closer than its own
    error model, so the minimum pairwise divergence is reported with every
    panel and quoted with every result.
    """
    import itertools

    import edlib

    d = [edlib.align(a, b, task="distance")["editDistance"] / max(len(a), len(b))
         for a, b in itertools.combinations(seqs, 2)]
    if not d:
        return None
    d.sort()
    return {"n_pairs": len(d), "min": round(d[0], 5), "median": round(d[len(d) // 2], 5),
            "mean": round(sum(d) / len(d), 5), "max": round(d[-1], 5)}


def write_panel(out: Path, names, seqs, ab, host_of, meta, reads_pe=None, reads_np=None):
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "truth_haplotypes.fa", "w") as f:
        for n, s, a in zip(names, seqs, ab):
            f.write(f">{n} host={host_of[n]} abundance={a:.6f}\n{s}\n")
    if reads_pe:
        r1, r2, t = reads_pe
        sq.write_fastq(out / "illumina_1.fq", r1, "I")
        sq.write_fastq(out / "illumina_2.fq", r2, "I")
        with open(out / "illumina_truth.tsv", "w") as f:
            f.write("pair\thaplotype\tstart\tend\tstrand\n")
            for rid, h, a, b, st in t:
                f.write(f"{rid}\t{names[h]}\t{a}\t{b}\t{st}\n")
    for chem, pack in (reads_np or {}).items():
        reads, t = pack
        sq.write_fastq(out / f"nanopore_{chem}.fq", reads, "5")
        with open(out / f"nanopore_{chem}.fa", "w") as f:
            for rid, s in reads:
                f.write(f">{rid}\n{s}\n")
        with open(out / f"nanopore_{chem}_truth.tsv", "w") as f:
            f.write("read\thaplotype\tstart\tend\tstrand\tsubs\tdels\tins\thp_errors\n")
            for row in t:
                f.write("\t".join([row[0], names[row[1]]] + [str(x) for x in row[2:]]) + "\n")
    (out / "panel.json").write_text(json.dumps(meta, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--refdir", type=Path, default=Path("data/env"))
    ap.add_argument("--hosts", nargs="+", required=True)
    ap.add_argument("--n-per-host", type=int, default=None)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--coverage", type=float, default=500.0, help="Illumina coverage")
    ap.add_argument("--read-len", type=int, default=250)
    ap.add_argument("--nanopore-reads", type=int, default=3000)
    ap.add_argument("--chem", nargs="+", default=["r9", "r10"])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    names, seqs, ab, host_of = build_panel(rng, a.refdir, a.hosts, a.n_per_host, a.sigma)
    pe = sq.paired_end(rng, seqs, ab, a.coverage, read_len=a.read_len)
    nps = {c: sq.nanopore(rng, seqs, ab, a.nanopore_reads, c) for c in a.chem}
    meta = vars(a) | {"out": str(a.out), "refdir": str(a.refdir), "n_haplotypes": len(seqs),
                      "lengths": [len(s) for s in seqs], "divergence": divergence_stats(seqs)}
    write_panel(a.out, names, seqs, ab, host_of, meta, pe, nps)
    print(f"{len(seqs)} haplotypes, {len(pe[0])} pairs, " + ", ".join(f"{c}:{len(nps[c][0])}" for c in nps))


if __name__ == "__main__":
    main()
