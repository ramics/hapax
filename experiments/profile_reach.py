#!/usr/bin/env python3
"""Profile reach, and the pair as one molecule rather than two reads.

WHY THIS EXISTS
---------------
SPEC.md §4.2 names "profile reach" as one of the five axes of the signal-versus-noise
measurement, and §7 says Pfam-A profiles are "the only axis still alive past 20%
divergence, where DIAMOND and every nucleotide index read zero". That claim is currently
an ASSERTION. It stands on the DIAMOND half of `results/spike_holdout/`'s divergence
ladder — which was measured — and on the absence of a measured profile half. This script
measures the profile half on the same ladder.

It answers two questions in one pass, because they share the ladder and the searches:

  §3  Where does profile search actually stop, relative to translated sequence search?
      The deliverable is NOT the profile hit rate. A Pfam domain hit means "carries a
      recognisable domain", and this project's background is itself mostly coding
      (`results/readlevel_features/`), so a hit rate is a familiarity signal wearing a
      novelty hat — the same way five single-read axes have already failed on direction
      rather than power. The reportable quantity is the GAP: queries with no significant
      translated-search hit but a significant profile hit. That gap is a distance
      measurement. Its polarity is measured here on known material (real library reads),
      not assumed.

  §4  Does treating the pair as one molecule extend that reach?
      A pair is one molecule of ~194 nt median insert, and 26-27% of pairs overlap their
      own mate (`results/quality_profile/`), so those merge into a single read of 60-199
      nt where a mate alone is 107 or 93. At 107 nt the true frame is the only stop-free
      frame in just 0.528 of coding reads (SPEC.md §4.1), so length should buy frame
      evidence as well as alignment. Measured as a PAIRED comparison on the same pairs —
      merged against mate 1 of the same molecule — because overlapping pairs are short
      inserts and comparing them to all mates would compare insert lengths instead.

CALIBRATION THAT TRAVELS WITH EVERY NUMBER (CLAUDE.md non-negotiables)
----------------------------------------------------------------------
  null            composition-matched mononucleotide shuffles of the same queries.
  known-positive  real reads from this library, which the cascade's indexes do identify.
                  This is the control that gives the gap its DIRECTION. If background
                  reads show the same no-DIAMOND-yes-Pfam gap, the gap describes the
                  database and not the query, and that is the finding.
  operating point every search threshold is reported as a curve, never as one number.

READ THE CONTROL FIRST. In all three failures of 2026-11-20 the control moved with the
test.

STATED BIASES, with direction
-----------------------------
1. Uniform substitution is a crude model of real divergence, exactly as in
   `results/spike_holdout/`. Rungs are ORDINAL, not family distances. Direction: uniform
   substitution spreads changes evenly instead of concentrating them in variable loops,
   so it is probably HARSHER on profiles than real divergence at the same percentage.
2. The merged read is simulated as the insert carrying errors at the same per-base rate,
   with NO error-correction credit for the overlap. We have a measured error rate but no
   quality model, and inventing one to award merging an accuracy benefit would be
   unfalsifiable. Direction: merged reach is a LOWER bound; a real merger with quality
   scores would do better.
3. The holdout is one genome, a picornavirus with a single long polyprotein ORF. Direction:
   favourable to profiles — a genome with little coding sequence would do worse. Said
   plainly rather than generalised from.
4. Background reads are the library's first N pairs in lane order, the same sampling
   `results/pair_ladder/`'s background uses. Direction: unknown, and it is a sample of
   one library on one day; drift makes every such figure a statement about SRR23973193
   (`results/dna_rna_pair/`).

MEASURED CONSTANTS, all from committed results in this repo
-----------------------------------------------------------
  read lengths          107 (mate 1) and 93 (mate 2), after the 14 nt non-template
                        prefix trim                     `results/nontemplate_prefix/`
  insert distribution   median 194, p95 363, p99 485     `results/insert_profile/`
  error rate            0.00068 along the body, 0.0019 at the final base, reference-free
                        from self-overlapping pairs      `results/quality_profile/`
                        (this supersedes `results/berlin_quality/`'s Q-derived rates,
                        which are upper bounds by ~2.6x)
  divergence rungs      0 / 5 / 10 / 20 / 30 %           `results/spike_holdout/`
  holdout genome        D00239, human rhinovirus A1B, the ONLY clean holdout of the
                        spike ladder — absent from this sample's own assembly, unlike
                        the other two genomes used there

USAGE
-----
  python3 experiments/profile_reach.py --seeds 1,2,3,4,5 --out results/profile_reach
Stages are cached in the output directory, so a search can be resumed without repeating
the simulation.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# measured constants (see docstring for provenance)
# ---------------------------------------------------------------------------
LEN_R1 = 107
LEN_R2 = 93
INSERT_MEDIAN = 194.0
INSERT_P95 = 363.0
INSERT_MIN = 60
INSERT_MAX = 1158
ERR_BODY = 0.00068
ERR_LAST = 0.0019
RUNGS = (0, 5, 10, 20, 30)

# lognormal fitted to the measured median and p95; p99 of the fit is 471 against a
# measured 485, so the fit is slightly tight in the tail. Stated, not corrected.
INSERT_MU = math.log(INSERT_MEDIAN)
INSERT_SIGMA = (math.log(INSERT_P95) - math.log(INSERT_MEDIAN)) / 1.6448536269514722

CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}
COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"
MIN_PEPTIDE = 10          # aa; shorter stop-free segments carry no profile information
OVERLAP_MIN = 20          # nt; minimum overlap to merge a real pair
OVERLAP_MAX_MM = 0.05     # mismatch fraction allowed in a real pair's overlap


def as_str(v) -> str:
    """pyhmmer returns names as bytes in some versions and str in others."""
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(nt: str) -> str:
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def six_frame_peptides(nt: str) -> list[str]:
    """Stop-free peptide segments of at least MIN_PEPTIDE aa, over all six frames.

    Splitting at stops rather than masking them is what a profile search needs, and the
    number of segments is itself a readout: at 107 nt the true frame is the only stop-free
    frame in barely half of coding reads, so a query normally offers several candidates
    and the instrument is asked to choose.
    """
    out = []
    for strand in (nt, revcomp(nt)):
        for frame in range(3):
            aa = translate(strand[frame:])
            for seg in aa.replace("X", "*").split("*"):
                if len(seg) >= MIN_PEPTIDE:
                    out.append(seg)
    return out


# ---------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------
def diverge(genome: str, pct: float, rng: random.Random) -> str:
    """Walk a genome away from itself by uniform substitution, as the spike ladder does."""
    if pct <= 0:
        return genome
    g = list(genome)
    n = len(g)
    k = int(round(n * pct / 100.0))
    for pos in rng.sample(range(n), k):
        cur = g[pos]
        if cur not in BASES:
            continue
        g[pos] = rng.choice([b for b in BASES if b != cur])
    return "".join(g)


def add_errors(seq: str, rng: random.Random) -> str:
    """Substitution errors at the library's reference-free measured rates."""
    s = list(seq)
    for i in range(len(s)):
        p = ERR_LAST if i == len(s) - 1 else ERR_BODY
        if rng.random() < p:
            cur = s[i]
            if cur in BASES:
                s[i] = rng.choice([b for b in BASES if b != cur])
    return "".join(s)


def draw_insert(rng: random.Random) -> int:
    while True:
        v = int(round(rng.lognormvariate(INSERT_MU, INSERT_SIGMA)))
        if INSERT_MIN <= v <= INSERT_MAX:
            return v


def simulate_pairs(genome: str, n: int, rng: random.Random) -> list[dict]:
    """Pairs from a genome, through this library's measured insert and error model.

    The merged condition carries NO error-correction credit (stated bias 2): it is the
    insert with one pass of errors, so it differs from a mate only in length.
    """
    pairs = []
    glen = len(genome)
    for i in range(n):
        ins = draw_insert(rng)
        if ins > glen:
            ins = glen
        start = rng.randrange(0, glen - ins + 1)
        frag = genome[start:start + ins]
        if rng.random() < 0.5:
            frag = revcomp(frag)
        m1 = add_errors(frag[:min(LEN_R1, ins)], rng)
        m2 = add_errors(revcomp(frag[max(0, ins - LEN_R2):]), rng)
        overlap = min(LEN_R1, ins) + min(LEN_R2, ins) - ins
        merged = add_errors(frag, rng) if overlap >= OVERLAP_MIN else None
        pairs.append({"id": i, "insert": ins, "m1": m1, "m2": m2, "merged": merged,
                      "overlap": overlap})
    return pairs


def shuffle_seq(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


# ---------------------------------------------------------------------------
# real reads: the known-positive control that gives every field its direction
# ---------------------------------------------------------------------------
def read_fastq_pairs(p1: Path, p2: Path, n: int, trim2: int = 14) -> list[dict]:
    pairs = []
    with gzip.open(p1, "rt") as f1, gzip.open(p2, "rt") as f2:
        while len(pairs) < n:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline().strip()
            f1.readline(); f1.readline()
            f2.readline()
            s2 = f2.readline().strip()[trim2:]
            f2.readline(); f2.readline()
            if "N" in s1 or "N" in s2:
                continue
            pairs.append({"id": len(pairs), "m1": s1, "m2": s2})
    return pairs


def merge_real_pair(m1: str, m2: str) -> tuple[str | None, int]:
    """Merge a real pair by suffix-prefix overlap, the way the error floor was measured.

    Returns (merged, overlap). The overlap is found by exhaustive offset scan against
    revcomp(m2); the acceptance rate of this function on the real library is itself a
    check, because `results/quality_profile/` measured 26-27% of pairs overlapping.
    """
    r2 = revcomp(m2)
    best = None
    for ov in range(min(len(m1), len(r2)), OVERLAP_MIN - 1, -1):
        a = m1[-ov:]
        b = r2[:ov]
        mm = sum(1 for x, y in zip(a, b) if x != y)
        if mm <= max(1, int(ov * OVERLAP_MAX_MM)):
            best = ov
            break
    if best is None:
        return None, 0
    return m1 + r2[best:], best


# ---------------------------------------------------------------------------
# query construction
# ---------------------------------------------------------------------------
def build_queries(seeds: list[int], genome: str, berlin: tuple[Path, Path],
                  n_sim: int, n_bg: int, log) -> list[tuple[str, str]]:
    """One flat list of (name, nucleotide sequence) for every cell of the design.

    Names encode seed | pile | rung | condition | pair id, so a single search pass can be
    split afterwards. Searching every cell in ONE pass is deliberate: a Pfam scan's cost
    is dominated by loading and setting up 21k models, and pricing anything here as a
    rate per query when most of it is a fixed scan is the error `results/diamond_scaling/`
    exists to stop.
    """
    queries: list[tuple[str, str]] = []
    overlap_stats = defaultdict(list)

    for seed in seeds:
        for rung in RUNGS:
            rng = random.Random((seed * 1000) + rung)
            div = diverge(genome, rung, rng)
            pairs = simulate_pairs(div, n_sim, rng)
            for pr in pairs:
                base = f"s{seed}|spike|d{rung}|"
                queries.append((base + f"m1|{pr['id']}", pr["m1"]))
                queries.append((base + f"m2|{pr['id']}", pr["m2"]))
                if pr["merged"] is not None:
                    queries.append((base + f"merged|{pr['id']}", pr["merged"]))
            overlap_stats[f"spike_d{rung}"].append(
                sum(1 for p in pairs if p["merged"] is not None) / len(pairs))

    # the known-positive control: real library reads, same conditions
    bg = read_fastq_pairs(berlin[0], berlin[1], n_bg * max(seeds))
    log(f"real pairs read: {len(bg)}")
    per = len(bg) // len(seeds)
    for si, seed in enumerate(seeds):
        chunk = bg[si * per:(si + 1) * per][:n_bg]
        merged_n = 0
        rng = random.Random(90000 + seed)
        for pr in chunk:
            base = f"s{seed}|background|dNA|"
            queries.append((base + f"m1|{pr['id']}", pr["m1"]))
            queries.append((base + f"m2|{pr['id']}", pr["m2"]))
            mg, ov = merge_real_pair(pr["m1"], pr["m2"])
            if mg is not None:
                merged_n += 1
                queries.append((base + f"merged|{pr['id']}", mg))
            # The null, composition-matched shuffles of the same reads. The MERGED shuffle
            # is the one that matters and it is why this condition exists: a longer query
            # can reach a fixed bitscore floor by chance more easily than a short one, so
            # a merged-against-mate gain measured at a fixed floor could be arithmetic
            # rather than recognition. The shuffle is searched at every length the real
            # comparison uses, so the chance rate is known at each of them instead of
            # assumed from the shortest.
            queries.append((f"s{seed}|shuffle|dNA|m1|{pr['id']}",
                            shuffle_seq(pr["m1"], rng)))
            queries.append((f"s{seed}|shuffle|dNA|m2|{pr['id']}",
                            shuffle_seq(pr["m2"], rng)))
            if mg is not None:
                queries.append((f"s{seed}|shuffle|dNA|merged|{pr['id']}",
                                shuffle_seq(mg, rng)))
        overlap_stats["background"].append(merged_n / max(1, len(chunk)))

    for k, v in sorted(overlap_stats.items()):
        log(f"overlap/merge rate {k}: {np.mean(v):.4f} (n={len(v)})")
    return queries


# ---------------------------------------------------------------------------
# searches
# ---------------------------------------------------------------------------
def write_fasta(path: Path, items) -> None:
    with open(path, "w") as fh:
        for name, seq in items:
            fh.write(f">{name}\n{seq}\n")


def run_diamond(queries, workdir: Path, dbfiles: list[Path], threads: int, log):
    """Translated search, as the reach instrument the profile half is measured against.

    --ultra-sensitive is an operating point and is reported as one; the per-hit numbers
    kept are bitscore, identity, alignment length AND E-value, because coverage and
    identity without the E-value beside them is what produced this project's withdrawn
    big-index reading.
    """
    dmnd = REPO / "bin_local" / "diamond"
    db = workdir / "prot.dmnd"
    if not db.exists():
        faa = workdir / "prot.faa"
        with open(faa, "w") as out:
            for f in dbfiles:
                op = gzip.open if f.suffix == ".gz" else open
                with op(f, "rt") as fh:
                    shutil.copyfileobj(fh, out)
        t0 = time.time()
        subprocess.run([str(dmnd), "makedb", "--in", str(faa), "-d", str(db),
                        "--threads", str(threads)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        log(f"diamond makedb: {time.time() - t0:.0f}s")
        faa.unlink()

    qfa = workdir / "queries.fa"
    write_fasta(qfa, queries)
    out = workdir / "diamond.tsv"
    # cache keyed on the exact query set, so a failure in a later stage does not pay for
    # this search twice and a changed query set can never silently reuse an old result
    key = hashlib.sha256("\n".join(n for n, _ in queries).encode()).hexdigest()
    keyfile = workdir / "diamond.queries.sha256"
    if out.exists() and keyfile.exists() and keyfile.read_text().strip() == key:
        log(f"diamond: reusing cached {out.name} for an identical query set")
        return parse_diamond(out), 0.0
    t0 = time.time()
    subprocess.run([str(dmnd), "blastx", "-d", str(db), "-q", str(qfa),
                    "-o", str(out), "--ultra-sensitive", "--max-target-seqs", "1",
                    "--threads", str(threads), "--quiet", "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"],
                   check=True)
    elapsed = time.time() - t0
    log(f"diamond blastx: {elapsed:.0f}s for {len(queries)} queries")
    keyfile.write_text(key)
    return parse_diamond(out), elapsed


def parse_diamond(out: Path) -> dict:
    best: dict[str, dict] = {}
    with open(out) as fh:
        for line in fh:
            q, subj, pid, ln, ev, bs = line.rstrip("\n").split("\t")
            score = float(bs)
            if q not in best or score > best[q]["bitscore"]:
                best[q] = {"subject": subj, "pident": float(pid), "length": int(ln),
                           "evalue": float(ev), "bitscore": score}
    return best


def run_pfam(queries, workdir: Path, pfam_hmm: Path, threads: int, log):
    """Pfam-A over six frames, reported as a curve with the curated threshold on it.

    The search runs WITHOUT --cut_ga and keeps the per-domain E-value, because the
    gathering thresholds are built for full-length proteins and are known to be
    conservative on fragments — a pilot on real 93 nt reads fired on 7 of 400 at GA, so
    reporting only GA would measure the threshold rather than the reach. Instead every
    hit down to E <= 10 is kept, and whether it clears its family's own GA cutoff is
    recorded as ONE operating point on that curve. That is the repo's rule about
    thresholds dominating everything else, applied to a profile search.

    Cost shape, measured before sizing anything (`results/diamond_scaling/`'s lesson):
    37.5 s of fixed model setup plus about 0.9 ms per peptide, so one pass over every cell
    of the design is far cheaper than a pass per cell.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    targets = []
    nseg = Counter()
    for name, nt in queries:
        peps = six_frame_peptides(nt)
        nseg[name] = len(peps)
        for i, pep in enumerate(peps):
            targets.append(pyhmmer.easel.TextSequence(
                name=f"{name}|f{i}".encode(), sequence=pep).digitize(alphabet))
    log(f"pfam targets: {len(targets)} peptide segments from {len(queries)} queries")

    best: dict[str, dict] = {}
    t0 = time.time()
    with pyhmmer.plan7.HMMFile(pfam_hmm) as hf:
        for hits in pyhmmer.hmmsearch(hf, targets, cpus=threads, E=10.0):
            fam = as_str(hits.query.name)
            ga = None
            try:
                cut = hits.query.cutoffs
                if cut.gathering_available():
                    ga = cut.gathering[1]        # per-domain gathering threshold
            except Exception:
                ga = None
            for hit in hits:
                qname = as_str(hit.name).rsplit("|f", 1)[0]
                for dom in hit.domains:
                    rec = {"family": fam, "evalue": float(dom.i_evalue),
                           "bitscore": float(dom.score),
                           "ga": bool(ga is not None and dom.score >= ga)}
                    if qname not in best or rec["bitscore"] > best[qname]["bitscore"]:
                        best[qname] = rec
    elapsed = time.time() - t0
    log(f"pfam hmmsearch: {elapsed:.0f}s, queries with a hit: {len(best)}")
    return best, dict(nseg), elapsed


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------
DIAMOND_CURVE = (30, 40, 50, 60, 80)      # bitscore operating points
PFAM_CURVE = (1e-2, 1e-3, 1e-5, 1e-10)    # per-domain E-value operating points


def summarise(queries, dia, pfam, nseg, seeds):
    """Per-cell summary. Every rate is reported as a curve over operating points.

    The headline field is `gap`: no significant translated hit, yes profile hit. Its
    polarity is only readable against the background and shuffle rows, which is why they
    are in the same table rather than in a separate control section.
    """
    cells = defaultdict(list)
    for name, nt in queries:
        seed, pile, rung, cond, _ = name.split("|")
        cells[(seed, pile, rung, cond)].append((name, nt))

    rows = []
    for (seed, pile, rung, cond), items in sorted(cells.items()):
        n = len(items)
        row = {"seed": seed, "pile": pile, "rung": rung, "condition": cond, "n": n,
               "mean_query_len": float(np.mean([len(s) for _, s in items])),
               "mean_peptide_segments": float(np.mean([nseg.get(q, 0) for q, _ in items]))}
        for bs in DIAMOND_CURVE:
            row[f"diamond_bs{bs}"] = sum(
                1 for q, _ in items if dia.get(q, {}).get("bitscore", 0) >= bs) / n
        for ev in PFAM_CURVE:
            row[f"pfam_e{ev:g}"] = sum(
                1 for q, _ in items if pfam.get(q, {}).get("evalue", 9e9) <= ev) / n
        row["pfam_ga"] = sum(1 for q, _ in items if pfam.get(q, {}).get("ga")) / n
        # the gap, at every pairing of the two curves, plus at the curated threshold
        for bs in DIAMOND_CURVE:
            for ev in PFAM_CURVE:
                row[f"gap_bs{bs}_e{ev:g}"] = sum(
                    1 for q, _ in items
                    if dia.get(q, {}).get("bitscore", 0) < bs
                    and pfam.get(q, {}).get("evalue", 9e9) <= ev) / n
            row[f"gap_bs{bs}_ga"] = sum(
                1 for q, _ in items
                if dia.get(q, {}).get("bitscore", 0) < bs
                and pfam.get(q, {}).get("ga")) / n
        hits = [dia[q] for q, _ in items if q in dia]
        if hits:
            row["diamond_median_evalue"] = float(np.median([h["evalue"] for h in hits]))
            row["diamond_median_pident"] = float(np.median([h["pident"] for h in hits]))
            row["diamond_median_alnlen"] = float(np.median([h["length"] for h in hits]))
            row["diamond_any_hit"] = len(hits) / n
        else:
            row["diamond_any_hit"] = 0.0
        rows.append(row)
    return rows


def paired_merged_vs_mate(queries, dia, pfam, seeds):
    """§4's paired comparison: merged against mate 1 OF THE SAME MOLECULE.

    Restricted to pairs that have a merged condition, because overlapping pairs are short
    inserts; comparing merged reads to all mate 1 reads would compare insert lengths and
    call it a result.
    """
    have = defaultdict(dict)
    for name, nt in queries:
        seed, pile, rung, cond, pid = name.split("|")
        have[(seed, pile, rung, pid)][cond] = (name, len(nt))

    rows = []
    groups = defaultdict(list)
    for key, conds in have.items():
        if "merged" in conds and "m1" in conds:
            groups[(key[0], key[1], key[2])].append(conds)
    for (seed, pile, rung), items in sorted(groups.items()):
        n = len(items)
        row = {"seed": seed, "pile": pile, "rung": rung, "n_pairs_overlapping": n,
               "mean_merged_len": float(np.mean([c["merged"][1] for c in items])),
               "mean_m1_len": float(np.mean([c["m1"][1] for c in items]))}
        for label, bs in (("bs40", 40), ("bs50", 50)):
            m = sum(1 for c in items
                    if dia.get(c["merged"][0], {}).get("bitscore", 0) >= bs) / n
            s = sum(1 for c in items
                    if dia.get(c["m1"][0], {}).get("bitscore", 0) >= bs) / n
            row[f"diamond_{label}_merged"] = m
            row[f"diamond_{label}_m1"] = s
        for ev in (1e-2, 1e-5):
            m = sum(1 for c in items
                    if pfam.get(c["merged"][0], {}).get("evalue", 9e9) <= ev) / n
            s = sum(1 for c in items
                    if pfam.get(c["m1"][0], {}).get("evalue", 9e9) <= ev) / n
            row[f"pfam_e{ev:g}_merged"] = m
            row[f"pfam_e{ev:g}_m1"] = s
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
def git_sha() -> str:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, check=True,
                             capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO, check=True,
                               capture_output=True, text=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--out", default="results/profile_reach")
    ap.add_argument("--n-sim", type=int, default=600,
                    help="simulated pairs per seed per rung")
    ap.add_argument("--n-bg", type=int, default=600,
                    help="real pairs per seed for the known-positive control")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--work", default=None)
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    outdir = REPO / args.out
    outdir.mkdir(parents=True, exist_ok=True)
    work = Path(args.work) if args.work else outdir / "work"
    work.mkdir(parents=True, exist_ok=True)

    logpath = outdir / "run.log"
    logfh = open(logpath, "a")

    def log(msg: str) -> None:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        logfh.write(line + "\n")
        logfh.flush()

    t_start = time.time()
    genome_path = REPO / "data" / "holdout" / "D00239.fa"
    genome = "".join(l.strip() for l in open(genome_path) if not l.startswith(">")).upper()
    log(f"holdout genome {genome_path.name}: {len(genome)} nt")

    berlin = (REPO / "data" / "berlin" / "SRR23973193_1.fastq.gz",
              REPO / "data" / "berlin" / "SRR23973193_2.fastq.gz")
    pfam_hmm = REPO / "data" / "pfam" / "Pfam-A.hmm"
    dbfiles = [REPO / "data" / "prot" / "viral.1.protein.faa.gz",
               REPO / "data" / "prot" / "uniprot_sprot.fasta.gz"]
    for p in (*berlin, pfam_hmm, *dbfiles):
        if not p.exists():
            log(f"MISSING INPUT: {p}")
            return 2

    queries = build_queries(seeds, genome, berlin, args.n_sim, args.n_bg, log)
    log(f"queries: {len(queries)}")

    dia, t_dia = run_diamond(queries, work, dbfiles, args.threads, log)
    pfam, nseg, t_pfam = run_pfam(queries, work, pfam_hmm, args.threads, log)

    rows = summarise(queries, dia, pfam, nseg, seeds)
    paired = paired_merged_vs_mate(queries, dia, pfam, seeds)

    with open(outdir / "cells.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    with open(outdir / "paired_merged.json", "w") as fh:
        json.dump(paired, fh, indent=1)

    script = Path(__file__).read_bytes()
    manifest = {
        "script": "experiments/profile_reach.py",
        "script_sha256": hashlib.sha256(script).hexdigest(),
        "git_sha": git_sha(),
        "seeds": seeds,
        "params": {
            "n_sim_pairs_per_seed_per_rung": args.n_sim,
            "n_background_pairs_per_seed": args.n_bg,
            "rungs_pct": list(RUNGS),
            "read_lengths": [LEN_R1, LEN_R2],
            "insert_lognormal": {"mu": INSERT_MU, "sigma": INSERT_SIGMA,
                                 "fitted_to": {"median": INSERT_MEDIAN, "p95": INSERT_P95},
                                 "truncated": [INSERT_MIN, INSERT_MAX]},
            "error_rates": {"body": ERR_BODY, "final_base": ERR_LAST},
            "min_peptide_aa": MIN_PEPTIDE,
            "diamond_mode": "--ultra-sensitive",
            "diamond_db": [f.name for f in dbfiles],
            "pfam_threshold": "--cut_ga (family gathering thresholds)",
            "overlap_min_nt": OVERLAP_MIN,
            "overlap_max_mismatch_frac": OVERLAP_MAX_MM,
        },
        "timings_s": {"diamond": t_dia, "pfam": t_pfam,
                      "total": time.time() - t_start},
        "inputs": {
            "holdout": "D00239 (human rhinovirus A1B), ENA fasta",
            "holdout_len_nt": len(genome),
            "library": "SRR23973193 (ENA fastq), mate 2 trimmed 14 nt",
            "pfam": str(pfam_hmm.name),
        },
    }
    with open(outdir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    log(f"done in {time.time() - t_start:.0f}s -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
