"""How far apart must two positions be before their novelty decorrelates?

This is the question the wide-insert hypothesis turns on. Spreading the mates
apart only buys anything if the two mates become *independent observations*: a
second look at the same gene is nearly the same look. So the quantity to measure
is the correlation length of novelty along a genome, and then to compare it with
the separations a real library can deliver (measured in results/insert_profile/:
median insert 208 nt, p95 387, p99 ~515).

Novelty is defined here the way a detection method would have to see it: for a
genome whose own family has been removed from the database, which of its
positions still have a retained homolog, and which have none? That is a graded,
per-position track rather than a verdict about the genome, and it is built with
a sensitive profile-free protein search (pyhmmer's phmmer, HMMER3) of each
annotated protein against the family-removed RefSeq viral proteome. Per-residue
coverage comes from the alignment envelopes of the hits, which is why phmmer is
the right tool here and a bitscore per read is not: we need *where* on the
protein the homology lands.

The prior arm (results/abundance_sweep/REPORT.md) found that with Rhabdoviridae
removed only 2 of rabies' 5 proteins keep any homolog at all. If that pattern is
general, novelty is organised in gene-sized blocks of 0.7-6.4 kb, and a 200-500
nt insert cannot straddle a block boundary often enough to matter. This script
measures the block structure instead of assuming it, and then runs the same
template sweep as experiments/mate_separation.py over the novelty track:
P(at least one mate novel), P(both novel), and the correlation of the mates'
novelty, against insert length.

Seeds vary the random template placements. The phmmer search and the track it
produces are deterministic given the database, so the spread across seeds is the
placement spread only, which is what the sweep's error bars are about.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

from mate_separation import fetch_gb
from peptide_absence import lineage_accessions


def parse_cds(path: Path) -> tuple[str, list[dict]]:
    """Sequence plus CDS records: product, protein sequence, genomic positions.

    The per-residue to per-genome-position map is built by walking the location's
    spans in transcription order, so joins (ribosomal slippage in ORF1ab,
    spliced viral transcripts) and complement strands both land in the right
    place. A residue maps to the 3 nucleotides of its codon.
    """
    text = path.read_text()
    seq = "".join(c for c in text.split("ORIGIN", 1)[1] if c.isalpha()).upper() \
        if "ORIGIN" in text else ""
    lines = text.split("\n")
    out: list[dict] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("     CDS "):
            loc = lines[i][21:].strip()
            while loc.count("(") > loc.count(")"):
                i += 1
                loc += lines[i][21:].strip()
            rev = "complement" in loc
            spans = []
            for part in loc.replace("join(", "").replace("complement(", "").replace(
                    "order(", "").replace(")", "").split(","):
                part = part.strip().replace("<", "").replace(">", "")
                if ".." in part:
                    a, b = part.split("..")[:2]
                    try:
                        spans.append((int(a) - 1, int(b)))
                    except ValueError:
                        pass
            prod, trans = None, None
            j = i + 1
            while j < len(lines) and (lines[j].startswith("                     ")
                                      or lines[j].startswith("                    /")):
                if "/product=" in lines[j]:
                    prod = lines[j].split("/product=", 1)[1].strip().strip('"')
                if "/translation=" in lines[j]:
                    buf = [lines[j].split("/translation=", 1)[1].strip().strip('"')]
                    while not buf[-1].endswith('"') and '"' not in buf[-1][1:] and j + 1 < len(lines):
                        j += 1
                        buf.append(lines[j].strip())
                        if buf[-1].endswith('"'):
                            break
                    trans = "".join(buf).replace('"', "")
                j += 1
            if trans and spans:
                positions: list[int] = []
                order = sorted(spans, reverse=True) if rev else sorted(spans)
                for a, b in order:
                    rng = range(b - 1, a - 1, -1) if rev else range(a, b)
                    positions.extend(rng)
                out.append({"product": prod or "?", "protein": trans,
                            "spans": spans, "reverse": rev, "nt_positions": positions})
            i = j
            continue
        i += 1
    return seq, out


def write_db(faa: Path, exclude: set[str], out: Path) -> Path:
    """The family-removed proteome, as plain FASTA for phmmer."""
    if out.exists() and out.stat().st_size > 1000:
        return out
    keep = True
    with gzip.open(faa, "rt") as fh, out.open("w") as oh:
        for line in fh:
            if line[0] == ">":
                keep = line[1:].split()[0] not in exclude
            if keep:
                oh.write(line)
    return out


def _txt(v) -> str:
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)


def phmmer_coverage(proteins: list[dict], db: Path, threads: int,
                    evalue: float) -> dict[int, dict]:
    """Per-residue homology coverage for each protein, from phmmer envelopes."""
    import pyhmmer

    alpha = pyhmmer.easel.Alphabet.amino()
    ok = set("ACDEFGHIKLMNPQRSTVWYX")
    queries = []
    for n, p in enumerate(proteins):
        # GenBank translations carry U (selenocysteine), O, J and stray stop
        # characters that HMMER's amino alphabet rejects outright; they are mapped
        # to X so the query keeps its length and the residue-to-genome map stays
        # aligned, rather than dropping characters and shifting every coordinate.
        seq = "".join(c if c in ok else "X" for c in p["protein"].upper())
        queries.append(pyhmmer.easel.TextSequence(
            name=f"q{n}".encode(), sequence=seq).digitize(alpha))
    with pyhmmer.easel.SequenceFile(str(db), digital=True, alphabet=alpha) as sf:
        targets = sf.read_block()
    cov: dict[int, dict] = {}
    for hits in pyhmmer.hmmer.phmmer(queries, targets, cpus=threads, E=evalue):
        qi = int(_txt(hits.query.name)[1:])
        L = len(proteins[qi]["protein"])
        track = np.zeros(L, dtype=np.int32)
        best = None
        n_hits = 0
        for hit in hits:
            if hit.evalue > evalue:
                continue
            n_hits += 1
            for dom in hit.domains:
                if dom.i_evalue > evalue:
                    continue
                a = max(0, dom.alignment.hmm_from - 1)
                b = min(L, dom.alignment.hmm_to)
                track[a:b] += 1
            if best is None or hit.score > best[1]:
                best = (_txt(hit.name), float(hit.score), float(hit.evalue))
        cov[qi] = {"residues": L, "covered": track,
                   "n_hits": n_hits,
                   "best_hit": best[0] if best else None,
                   "best_bitscore": best[1] if best else None,
                   "best_evalue": best[2] if best else None}
    for n in range(len(proteins)):
        cov.setdefault(n, {"residues": len(proteins[n]["protein"]),
                           "covered": np.zeros(len(proteins[n]["protein"]), dtype=np.int32),
                           "n_hits": 0, "best_hit": None, "best_bitscore": None,
                           "best_evalue": None})
    return cov


def blocks(track: np.ndarray, value: int) -> list[int]:
    """Run lengths of a given value: the block structure of novelty."""
    out, run = [], 0
    for v in track:
        if v == value:
            run += 1
        elif run:
            out.append(run)
            run = 0
    if run:
        out.append(run)
    return out


def sweep_novelty(novel: np.ndarray, coding: np.ndarray, inserts: list[int],
                  len1: int, len2: int, min_nt: int, n_starts: int,
                  rng: random.Random) -> dict:
    """P(mate novel) and the mates' novelty correlation, by insert length."""
    L = novel.size
    pn = np.concatenate([[0], np.cumsum(novel)])
    pc = np.concatenate([[0], np.cumsum(coding)])
    out = {}
    for I in inserts:
        if I > L or I < max(len1, len2):
            continue
        a_n = np.zeros(n_starts)
        b_n = np.zeros(n_starts)
        for i in range(n_starts):
            s = rng.randrange(0, L - I + 1)
            e = s + I
            a_n[i] = (pn[s + len1] - pn[s]) >= min_nt
            b_n[i] = (pn[e] - pn[e - len2]) >= min_nt
        pa, pb = a_n.mean(), b_n.mean()
        both = float((a_n * b_n).mean())
        den = (pa * (1 - pa) * pb * (1 - pb)) ** 0.5
        out[I] = {
            "p_mate1_novel": round(float(pa), 4),
            "p_mate2_novel": round(float(pb), 4),
            "p_either_novel": round(float(((a_n + b_n) > 0).mean()), 4),
            "p_both_novel": round(both, 4),
            "p_both_if_independent": round(float(pa * pb), 4),
            "rho_novelty": round(float((both - pa * pb) / den), 4) if den > 0 else None,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--accession", default="NC_001542.1")
    ap.add_argument("--exclude-family", default="Rhabdoviridae")
    ap.add_argument("--faa", type=Path, default=Path("data/prot/viral.1.protein.faa.gz"))
    ap.add_argument("--gpff", type=Path, default=Path("data/prot/viral.1.protein.gpff.gz"))
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--len1", type=int, default=107)
    ap.add_argument("--len2", type=int, default=93)
    ap.add_argument("--min-nt", type=int, default=60)
    ap.add_argument("--inserts", type=int, nargs="+",
                    default=[107, 150, 200, 250, 300, 400, 500, 700, 1000, 1500,
                             2000, 3000, 5000])
    ap.add_argument("--n-starts", type=int, default=20000)
    ap.add_argument("--cache", type=Path, default=Path("data/novelty_blocks"))
    ap.add_argument("--out", type=Path, default=Path("results/novelty_blocks"))
    args = ap.parse_args()

    rng = random.Random(args.seed)
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    seq, cdss = parse_cds(fetch_gb(args.accession, Path("data/mate_separation")))
    if not cdss:
        raise SystemExit(f"no annotated CDS with translations in {args.accession}")
    excl = lineage_accessions(args.gpff, args.exclude_family)
    db = write_db(args.faa, excl,
                  args.cache / f"viral_no_{args.exclude_family.lower()}.faa")
    print(f"[db] {len(excl)} proteins of {args.exclude_family} removed -> {db}")

    # the search is deterministic, so cache its per-residue coverage: seeds differ
    # only in the template placements of the sweep, not in the homology track
    cov_cache = args.cache / f"cov_{args.accession}_no_{args.exclude_family.lower()}_e{args.evalue}.npz"
    meta_cache = cov_cache.with_suffix(".json")
    if cov_cache.exists() and meta_cache.exists():
        z = np.load(cov_cache)
        meta = json.loads(meta_cache.read_text())
        cov = {int(k): {"covered": z[f"c{k}"], **meta[k]} for k in meta}
        print(f"[cov] reusing {cov_cache}")
    else:
        cov = phmmer_coverage(cdss, db, args.threads, args.evalue)
        np.savez_compressed(cov_cache, **{f"c{n}": c["covered"] for n, c in cov.items()})
        meta_cache.write_text(json.dumps({
            str(n): {k: v for k, v in c.items() if k != "covered"}
            for n, c in cov.items()}))

    # genome-wide tracks: coding, and novel (coding but no retained homolog here)
    novel = np.zeros(len(seq), dtype=np.int8)
    coding = np.zeros(len(seq), dtype=np.int8)
    per_protein = []
    for n, p in enumerate(cdss):
        c = cov[n]
        covered = c["covered"]
        for res_i, nt_start in enumerate(range(0, len(p["nt_positions"]), 3)):
            if res_i >= covered.size:
                break
            for nt in p["nt_positions"][nt_start:nt_start + 3]:
                if nt < len(seq):
                    coding[nt] = 1
                    if covered[res_i] == 0:
                        novel[nt] = 1
        per_protein.append({
            "product": p["product"],
            "residues": int(c["residues"]),
            "fraction_residues_with_homolog": round(float((covered > 0).mean()), 4),
            "n_hits": int(c["n_hits"]),
            "best_hit": c["best_hit"],
            "best_bitscore": c["best_bitscore"],
            "best_evalue": c["best_evalue"],
        })
        print(f"  {p['product'][:40]:42s} {c['residues']:5d} aa  homolog-covered "
              f"{(covered > 0).mean():.3f}  best {c['best_bitscore']} {c['best_hit']}")

    nov_blocks = blocks(novel, 1)
    kn_blocks = [b for b in blocks(novel * 0 + (coding & (1 - novel)), 1)] \
        if False else blocks((coding.astype(int) - novel.astype(int) > 0).astype(int), 1)
    res = {
        "seed": args.seed,
        "accession": args.accession,
        "genome_length": len(seq),
        "exclude_family": args.exclude_family,
        "evalue": args.evalue,
        "proteins_removed_from_db": len(excl),
        "coding_fraction": round(float(coding.mean()), 4),
        "novel_fraction_of_genome": round(float(novel.mean()), 4),
        "novel_fraction_of_coding": round(float(novel.sum() / max(1, coding.sum())), 4),
        "per_protein": per_protein,
        "novel_block_lengths_nt": nov_blocks,
        "known_block_lengths_nt": kn_blocks,
        "novel_block_summary": {
            "n": len(nov_blocks),
            "median": statistics.median(nov_blocks) if nov_blocks else None,
            "max": max(nov_blocks) if nov_blocks else None,
            "mean": round(statistics.fmean(nov_blocks), 1) if nov_blocks else None,
        },
        "sweep": sweep_novelty(novel, coding, args.inserts, args.len1, args.len2,
                               args.min_nt, args.n_starts, rng),
    }
    (args.out / f"novelty_blocks_{args.accession}_s{args.seed}.json").write_text(
        json.dumps(res, indent=2))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (args.out / f"manifest_{args.accession}_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/novelty_blocks.py", "git_sha": sha, "seed": args.seed,
        "params": {k: str(v) for k, v in vars(args).items()},
        "wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2))
    print(json.dumps({k: v for k, v in res.items()
                      if k not in ("per_protein", "novel_block_lengths_nt",
                                   "known_block_lengths_nt", "sweep")}, indent=2))
    for I, v in res["sweep"].items():
        print(I, v)


if __name__ == "__main__":
    main()
