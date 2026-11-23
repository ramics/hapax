"""Score an assembler's output against a known haplotype panel.

The three numbers the stress test turns on, all from the same run, as
CLAUDE.md requires:

* **Recovery.** For each truth haplotype, the best contig by edit distance
  over the aligned span. A haplotype counts as *recovered* when some contig
  covers >= `cov` of it at >= `ident` identity (defaults 0.95 / 0.99, the
  thresholds SAVAGE and Haploflow report against). Reported both as a count
  and as the identity distribution, because "recovered at 99%" and
  "recovered at 99.9%" are different claims for a 0.3%-divergent swarm: a
  contig can hit 99% while being nobody's haplotype.
* **Chimaeras.** A contig is chimaeric when its two halves have different
  best-matching haplotypes *and* each half matches its own best better than
  the whole contig matches any single haplotype. Windowed: every 300 nt
  window is assigned a best haplotype, and a contig with switches between
  haplotypes of different hosts is an unambiguous chimaera; switches within
  one host are reported separately, since those may be real recombinants in
  the input or may be the assembler's invention.
* **Frameshifts.** Each contig is aligned to the best-matching truth
  haplotype; an indel in that alignment whose length is not a multiple of 3
  inside the coding span is a frameshift the assembler introduced, because
  the truth haplotypes are real deposited coding sequences and their own
  indels are codon-multiples. This is the number the codon-space argument
  lives or dies on, so it is counted per contig and per kb.

Alignment is edlib (infix mode, so a contig is placed inside a haplotype).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import edlib

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in reversed(s))


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


def best_match(contig, truth):
    """Best (identity, name, orientation, cigar, span) placing contig inside a haplotype."""
    best = None
    for q, orient in ((contig, "+"), (revcomp(contig), "-")):
        for name, t in truth.items():
            r = edlib.align(q, t, mode="HW", task="path")
            ident = 1.0 - r["editDistance"] / max(len(q), 1)
            if best is None or ident > best[0]:
                loc = r["locations"][0]
                best = (ident, name, orient, r["cigar"], loc)
    return best


def window_hosts(contig, truth, host_of, win=300, step=150):
    hits = []
    for i in range(0, max(len(contig) - win + 1, 1), step):
        w = contig[i : i + win]
        if len(w) < win // 2:
            break
        b = best_match(w, truth)
        hits.append((i, b[1], host_of.get(b[1], "?"), b[0]))
    return hits


def frameshifts(cigar):
    """Indel runs in an edlib path cigar whose length is not a multiple of 3."""
    n, fs, ind = "", 0, 0
    for ch in cigar:
        if ch.isdigit():
            n += ch
            continue
        k = int(n or 1)
        n = ""
        if ch in "ID":
            ind += 1
            if k % 3:
                fs += 1
    return fs, ind


def score(truth_fa: Path, contigs_fa: Path, min_len=500, cov=0.95, ident=0.99):
    truth_raw = read_fasta(truth_fa)
    truth, host_of, ab = {}, {}, {}
    for header, seq in truth_raw.items():
        parts = header.split()
        name = parts[0]
        truth[name] = seq
        for p in parts[1:]:
            if p.startswith("host="):
                host_of[name] = p[5:]
            if p.startswith("abundance="):
                ab[name] = float(p[10:])
    contigs = {k.split()[0]: v for k, v in read_fasta(contigs_fa).items() if len(v) >= min_len}

    per_contig, recovered = [], {}
    for cname, c in contigs.items():
        bident, bname, orient, cigar, loc = best_match(c, truth)
        fs, ind = frameshifts(cigar)
        hits = window_hosts(c, truth, host_of)
        hosts = [h for _, _, h, _ in hits]
        haps = [n for _, n, _, _ in hits]
        cross_host = len(set(hosts)) > 1
        within_host_switch = (not cross_host) and len(set(haps)) > 1
        per_contig.append(dict(contig=cname, length=len(c), best=bname, identity=round(bident, 5),
                               orientation=orient, frameshifts=fs, indels=ind,
                               fs_per_kb=round(1000 * fs / len(c), 3),
                               cross_host_chimaera=cross_host,
                               within_host_switch=within_host_switch,
                               window_haplotypes=len(set(haps))))
        tl = len(truth[bname])
        covered = (loc[1] - loc[0] + 1) / tl
        prev = recovered.get(bname)
        if prev is None or bident > prev[0]:
            recovered[bname] = (bident, covered, cname)

    full = {k: v for k, v in recovered.items() if v[0] >= ident and v[1] >= cov}
    n_fs_contigs = sum(1 for p in per_contig if p["frameshifts"])
    return dict(
        n_truth=len(truth), n_contigs=len(contigs),
        total_contig_bp=sum(len(c) for c in contigs.values()),
        n50=_n50([len(c) for c in contigs.values()]),
        recovered=len(full), recovered_frac=round(len(full) / max(len(truth), 1), 4),
        recovered_names=sorted(full),
        best_identity_median=round(_median([v[0] for v in recovered.values()]), 5) if recovered else None,
        cross_host_chimaeras=sum(1 for p in per_contig if p["cross_host_chimaera"]),
        within_host_switches=sum(1 for p in per_contig if p["within_host_switch"]),
        contigs_with_frameshift=n_fs_contigs,
        frameshifts_total=sum(p["frameshifts"] for p in per_contig),
        frameshifts_per_kb=round(1000 * sum(p["frameshifts"] for p in per_contig) / max(sum(len(c) for c in contigs.values()), 1), 4),
        per_contig=per_contig,
    )


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def _n50(lens):
    if not lens:
        return 0
    lens = sorted(lens, reverse=True)
    half = sum(lens) / 2
    t = 0
    for x in lens:
        t += x
        if t >= half:
            return x


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("truth", type=Path)
    ap.add_argument("contigs", type=Path)
    ap.add_argument("--min-len", type=int, default=500)
    ap.add_argument("--cov", type=float, default=0.95)
    ap.add_argument("--ident", type=float, default=0.99)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    r = score(a.truth, a.contigs, a.min_len, a.cov, a.ident)
    if a.json:
        a.json.write_text(json.dumps(r, indent=2))
    pc = r.pop("per_contig")
    print(json.dumps(r, indent=2))
    print(f"({len(pc)} contigs scored)")


if __name__ == "__main__":
    main()
