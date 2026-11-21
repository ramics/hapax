"""What mate separation buys: coding-region hits, as a function of insert length.

The hypothesis under test (Imo, 2026-11-20) is that a wide library template with
the mates far apart maximises the chance of hitting a coding region. The first
half of that is arithmetic on real annotation and can be settled without any new
sequencing: slide a template of length I across a real sequence, place a 107 nt
mate at its 5' end and a 93 nt mate at its 3' end (SRR23973193's geometry, mate 2
after its 14 nt non-template prefix is trimmed), and ask how often at least one
mate lands in coding sequence.

Why the answer is bounded before it is measured, and why it still has to be
measured. At insert = read length the two mates see the same region, so
P(at least one mate coding) is just the coding density d. As the gap grows the
mates decorrelate and the limit is 1 - (1-d)^2. The whole available gain is
therefore d(1-d) at most - large only when d is middling. Viral genomes run
d = 0.85-0.95, so the ceiling is a few points. But d is not the only thing that
matters: coding sequence is not scattered at random, it is in blocks with
untranslated ends, so the approach to the independent limit is governed by the
block structure of real sequence, which is exactly what a closed form cannot
tell you. Hence real annotation: GenBank CDS for reference genomes, and
Prodigal (pyrodigal, metagenomic mode) gene calls on the Berlin sample's own
assembled contigs, which is the composition a wastewater method actually meets.

The same sweep reports P(both mates coding), because for any score that needs a
coding read in *both* mates - the mate conjunction of
results/readlevel_features/REPORT.md §3.1 - width works against you, and the two
curves cross.

Each sequence is treated as a separate population of templates and results are
reported per sequence and pooled, with seeds differing in the random subsample of
start positions so the spread is a real spread rather than a rounding of the same
deterministic sweep.
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import subprocess
import time
import urllib.request
from pathlib import Path

EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text")

# Reference genomes with GenBank CDS annotation. Chosen to span what a wastewater
# RNA library actually contains plus the project's own spike and targets.
REFS = {
    "NC_001943.1": "Human astrovirus 1 (the Berlin probe's target family)",
    "NC_001542.1": "Rabies lyssavirus (the abundance sweep's held-out spike)",
    "K03455.1": "HIV-1 HXB2 (dual-coding stress case)",
    "NC_045512.2": "SARS-CoV-2 (dense, long ORF1ab)",
    "NC_001422.1": "phiX174 (overlapping genes, spike-in control)",
    "NC_000913.3": "E. coli K-12 MG1655 (high-density bacterial background)",
}


def fetch_gb(acc: str, cache: Path) -> Path:
    out = cache / f"{acc}.gb"
    if not out.exists() or out.stat().st_size < 1000:
        cache.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(EFETCH.format(acc=acc), timeout=300) as r:
            out.write_bytes(r.read())
    return out


def parse_gb(path: Path) -> tuple[str, list[tuple[int, int]]]:
    """Return (sequence, CDS intervals as 0-based half-open spans).

    Joined and complemented locations are flattened to their component spans:
    for a codingness track the strand and the join structure are irrelevant,
    only which bases are inside a CDS.
    """
    text = path.read_text()
    seq_part = text.split("ORIGIN", 1)[1] if "ORIGIN" in text else ""
    seq = "".join(c for c in seq_part if c.isalpha()).upper()
    spans: list[tuple[int, int]] = []
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("     CDS "):
            loc = line[21:].strip()
            while not loc.endswith(")") and ".." in loc and loc.count("(") > loc.count(")"):
                i += 1
                loc += lines[i][21:].strip()
            for part in loc.replace("join(", "").replace("complement(", "").replace(
                    "order(", "").replace(")", "").split(","):
                part = part.strip().replace("<", "").replace(">", "")
                if ".." in part:
                    a, b = part.split("..")[:2]
                    try:
                        spans.append((int(a) - 1, int(b)))
                    except ValueError:
                        pass
        i += 1
    return seq, spans


def coding_track(length: int, spans: list[tuple[int, int]]) -> bytearray:
    t = bytearray(length)
    for a, b in spans:
        a, b = max(0, a), min(length, b)
        for i in range(a, b):
            t[i] = 1
    return t


def prodigal_tracks(fasta: Path, min_len: int, max_seqs: int, seed: int):
    """Gene calls on real assembled contigs: the background's own coding density."""
    import pyrodigal

    orf = pyrodigal.GeneFinder(meta=True)
    seqs = []
    name, buf = None, []
    with fasta.open() as fh:
        for line in fh:
            if line[0] == ">":
                if name and len(s := "".join(buf)) >= min_len:
                    seqs.append((name, s))
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name and len(s := "".join(buf)) >= min_len:
        seqs.append((name, s))
    random.Random(seed).shuffle(seqs)
    seqs = seqs[:max_seqs]
    for nm, s in seqs:
        genes = orf.find_genes(s.encode())
        spans = [(g.begin - 1, g.end) for g in genes]
        yield nm, len(s), coding_track(len(s), spans)


def sweep(track: bytearray, inserts: list[int], len1: int, len2: int,
          min_coding: int, n_starts: int, rng: random.Random) -> dict:
    """P(mate coding) by insert length, from random template placements."""
    L = len(track)
    pre = [0] * (L + 1)
    for i, v in enumerate(track):
        pre[i + 1] = pre[i] + v
    out = {}
    for I in inserts:
        if I > L or I < max(len1, len2):
            continue
        starts = [rng.randrange(0, L - I + 1) for _ in range(n_starts)]
        n1 = n2 = both = either = 0
        for s in starts:
            a = pre[s + len1] - pre[s] >= min_coding
            e2 = s + I
            b = pre[e2] - pre[e2 - len2] >= min_coding
            n1 += a
            n2 += b
            both += a and b
            either += a or b
        k = len(starts)
        # correlation of codingness between the two mates at this separation
        pa, pb, pab = n1 / k, n2 / k, both / k
        den = (pa * (1 - pa) * pb * (1 - pb)) ** 0.5
        out[I] = {
            "p_mate1_coding": round(pa, 4),
            "p_mate2_coding": round(pb, 4),
            "p_either_coding": round(either / k, 4),
            "p_both_coding": round(pab, 4),
            "phi_codingness": round((pab - pa * pb) / den, 4) if den > 0 else None,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--len1", type=int, default=107)
    ap.add_argument("--len2", type=int, default=93)
    ap.add_argument("--min-coding", type=int, default=60,
                    help="bases of a mate inside a CDS for it to count as a coding read")
    ap.add_argument("--inserts", type=int, nargs="+",
                    default=[107, 150, 200, 250, 300, 400, 500, 700, 1000, 1500, 2000, 3000])
    ap.add_argument("--n-starts", type=int, default=20000)
    ap.add_argument("--contigs", type=Path, default=None,
                    help="assembled contigs for the Prodigal arm (real background)")
    ap.add_argument("--contig-min-len", type=int, default=3000)
    ap.add_argument("--max-contigs", type=int, default=300)
    ap.add_argument("--cache", type=Path, default=Path("data/mate_separation"))
    ap.add_argument("--out", type=Path, default=Path("results/mate_separation"))
    args = ap.parse_args()

    rng = random.Random(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    res = {"seed": args.seed, "len1": args.len1, "len2": args.len2,
           "min_coding": args.min_coding, "n_starts": args.n_starts,
           "references": {}, "contigs": None}

    for acc, label in REFS.items():
        seq, spans = parse_gb(fetch_gb(acc, args.cache))
        if not seq:
            continue
        track = coding_track(len(seq), spans)
        res["references"][acc] = {
            "label": label,
            "length": len(seq),
            "n_cds_spans": len(spans),
            "coding_density": round(sum(track) / len(track), 4),
            "sweep": sweep(track, args.inserts, args.len1, args.len2,
                           args.min_coding, args.n_starts, rng),
        }

    if args.contigs:
        per, dens, tot_b, tot_c = [], [], 0, 0
        agg: dict[int, list[dict]] = {}
        for nm, L, track in prodigal_tracks(args.contigs, args.contig_min_len,
                                            args.max_contigs, args.seed):
            d = sum(track) / L
            dens.append(d)
            tot_b += L
            tot_c += sum(track)
            s = sweep(track, args.inserts, args.len1, args.len2, args.min_coding,
                      max(400, args.n_starts // args.max_contigs), rng)
            for I, row in s.items():
                agg.setdefault(I, []).append(row)
            per.append({"contig": nm, "length": L, "coding_density": round(d, 4)})
        res["contigs"] = {
            "source": str(args.contigs),
            "n_contigs": len(per),
            "min_len": args.contig_min_len,
            "coding_density_pooled": round(tot_c / tot_b, 4) if tot_b else None,
            "coding_density_median_per_contig": round(statistics.median(dens), 4) if dens else None,
            "sweep_pooled": {
                I: {k: round(statistics.fmean([r[k] for r in rows if r[k] is not None]), 4)
                    for k in ("p_mate1_coding", "p_either_coding", "p_both_coding",
                              "phi_codingness")
                    if any(r[k] is not None for r in rows)}
                for I, rows in sorted(agg.items())
            },
            "per_contig": per[:50],
        }

    (args.out / f"mate_separation_s{args.seed}.json").write_text(json.dumps(res, indent=2))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/mate_separation.py", "git_sha": sha, "seed": args.seed,
        "params": vars(args) | {"contigs": str(args.contigs), "cache": str(args.cache),
                                "out": str(args.out)},
        "wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2, default=str))
    for acc, r in res["references"].items():
        print(acc, r["coding_density"], {I: v["p_either_coding"] for I, v in r["sweep"].items()})
    if res["contigs"]:
        print("contigs", res["contigs"]["coding_density_pooled"],
              {I: v["p_either_coding"] for I, v in res["contigs"]["sweep_pooled"].items()})


if __name__ == "__main__":
    main()
