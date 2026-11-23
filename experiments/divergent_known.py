"""Every identified contig, with the distance at which it was identified.

THE GAP THIS CLOSES, and it is a design fault rather than a missing feature. The cascade's
subtrahend is one bit. A contig matching a reference at 99.9% and one matching it at 66.5%
both get the single word `identified`, are subtracted in tier 1, and never reach the
residue. `strains` only looks at the residue, so neither is ever reported. The holdout
thread demonstrated this with ground truth: a withheld genome assembles out of its library
as one contig covering 96% of it, hits the rolled-back reference at 66.5% protein identity
over 620 aa, is removed as known, and is never mentioned again. **The distance is computed
and then discarded.**

So this stage runs over the IDENTIFIED pile, which is the opposite of everything else in
this repo. The residue is where we look for sequence with no relative; this is the case of
sequence with a relative that is far away, and it has had no owner.

WHAT A ROW IS. A distance and a coverage, never a score. Identity on its own is
uninterpretable - this repo has withdrawn a claim for quoting one without its alignment
length - so every row carries the aligned length, the fraction of the subject covered, the
fraction of the contig covered, the read depth behind it and the evenness of that depth.
A reader orders the table by whichever of those they care about.

THE BAND QUESTION THIS ANSWERS FIRST. Before deciding where to draw a line, the shape of
the identity distribution over the whole identified pile has to be known: if the mass is
one continuum there is no line to draw and the output is a ranked table, and if there is a
gap between reference-matching sequence and a divergent tail then the gap is the line. The
histogram is computed over every contig with a significant hit and is printed whatever it
says. The holdout thread measured this on a different library; this measures our own,
because a band's shape is a property of the library and its assembly.

THREE CONTROL ARMS, drawn by the same code as the rows.

  1. Near-reference contigs from this same assembly. A field firing equally on them is
     describing assembly, not distance.
  2. A dinucleotide-preserving shuffle of the same contigs through the same search. A
     mononucleotide shuffle inflates ORFs, because stop-codon density is a dinucleotide
     property; that lesson is from `results/residue_contig_annotation/`.
  3. The holdout thread's withheld genome is the only real positive control the project
     has, and it is in a different library. Its operating point is recorded here, not its
     rows, and the old mutated spike ladder is a reference point rather than ground truth
     for a real distance.

WHAT THIS CANNOT SAY. Not whether a divergent contig is a new species: that is a taxonomic
judgement with demarcation criteria per family, and nothing here fetches them. Not how far
it really is from what is known: every identity here comes from RefSeq viral proteins plus
Swiss-Prot, which are narrow, so each is an UPPER BOUND on distance until the nr and
core_nt pass reports. Direction known, magnitude not.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import library_facts  # noqa: E402
from shuffle_control import shuffle_all  # noqa: E402  the one copy of the null's shuffle

from divergent_contig_workup import (BANDS, MIN_ALN_AA, SIGNIFICANT_EVALUE, band,
                                     call_orfs, diamond_blastp, gene_role, map_back,
                                     pfam_scan, read_fasta, sha256_of)

ROOT = Path(__file__).resolve().parents[1]

#: Histogram edges for the identity distribution. Five-point bins across the whole range,
#: because the question is whether the 60-70% region is sparse and a coarser bin would hide
#: exactly that.
EDGES = list(range(20, 101, 5))


def diamond_blastx(fa: Path, db: Path, out: Path, binpath: Path, threads: int,
                   mode: str) -> dict:
    """Whole-contig best hit, carrying the lengths that make a coverage computable.

    `slen` and `qlen` are in the output format on purpose: an identity without the
    fraction of the subject it covers cannot distinguish a diverged full-length protein
    from a conserved fragment of one, and those are different findings.

    25 targets rather than 1, and qstart/qend with them, because the best hit alone
    cannot answer whether ONE known protein explains the contig or whether several
    unrelated ones each explain a piece of it. The best hit is still whichever target
    wins on bitscore; the extra targets only feed the coverage gap below.
    """
    cmd = [str(binpath), "blastx", "-d", str(db), "-q", str(fa), "-o", str(out),
           "-p", str(threads), mode, "--quiet", "--evalue", "1e-3",
           "--max-target-seqs", "25", "--outfmt", "6",
           "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
           "qlen", "slen", "qstart", "qend", "stitle"]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    best, hsps, qlen = {}, {}, {}
    with out.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            row = {"subject": f[1], "pident": float(f[2]), "aln_aa": int(f[3]),
                   "evalue": float(f[4]), "bitscore": float(f[5]),
                   "contig_nt": int(f[6]), "subject_aa": int(f[7]),
                   "qstart": int(f[8]), "qend": int(f[9]),
                   "title": (f[10] if len(f) > 10 else "")[:200]}
            if row["evalue"] < SIGNIFICANT_EVALUE:
                hsps.setdefault(f[0], {}).setdefault(row["subject"], []).append(
                    (min(row["qstart"], row["qend"]), max(row["qstart"], row["qend"])))
                qlen[f[0]] = row["contig_nt"]
            row["subject_coverage"] = round(row["aln_aa"] / max(1, row["subject_aa"]), 3)
            row["contig_coverage"] = round(3 * row["aln_aa"] / max(1, row["contig_nt"]), 3)
            prev = best.get(f[0])
            if prev is None or row["bitscore"] > prev["bitscore"]:
                best[f[0]] = row
    return {"best": best, "hsps": hsps, "qlen": qlen,
            "seconds": round(time.time() - t0, 1)}


# The shuffle used to live here, reimplemented, and it was WRONG in the one way that
# manufactures results: the naive Eulerian walk dead-ends on most real sequence, and the
# fallback was `return s` - the input, unchanged, presented as a null. That has now fired
# three times in this repo on three separate copies (results/holdout_run/ discarded 80 jobs
# over it; results/divergent_known/ withdrew four published-looking numbers; and this copy
# produced the 0.79% null in the first holdout_cascade divergent run, which is withdrawn).
# So it is imported from ONE module and never written again: the second writer cannot know
# which non-obvious condition the first was fixing. shuffle_control raises rather than
# returning anything, and reports sequences it could not shuffle by name.
def histogram(values: list) -> list:
    out = []
    for lo, hi in zip(EDGES, EDGES[1:]):
        out.append({"from_percent": lo, "to_percent": hi,
                    "contigs": sum(1 for v in values if lo <= v < hi)})
    out.append({"from_percent": EDGES[-1], "to_percent": 100.1,
                "contigs": sum(1 for v in values if v >= EDGES[-1])})
    return out


FIELDS = ("pident", "aln_aa", "subject_coverage", "contig_coverage", "bitscore")


def quantiles(vals: list) -> dict:
    """p10/median/p90 rather than a mean, because every one of these fields is skewed
    and a mean of bitscores is not a readable number."""
    if not vals:
        return {}
    v = sorted(vals)
    def at(f):
        return v[min(len(v) - 1, int(f * (len(v) - 1) + 0.5))]
    return {"p10": round(at(0.10), 3), "median": round(at(0.50), 3),
            "p90": round(at(0.90), 3), "n": len(v)}


def discrimination(sig: dict, nullsig: dict, n_contigs: int,
                   n_null: int | None = None) -> dict:
    """What separates a real identified contig from a dinucleotide shuffle of itself,
    searched the same way, and what does not.

    This block exists because the smoke run killed the obvious reading. Over 2,855
    identified contigs the shuffled null's median whole-contig protein identity was
    HIGHER than the real one, 42.6% against 41.3%, and its p10, p90, alignment length,
    subject coverage and bitscore all sat inside the real spread. A chance alignment to a
    narrow protein index lands in the same band as real distant homology, so nothing about
    the alignment, once it exists, says which of the two it is.

    Scope, because the obvious over-reading of that is wrong. This calibrates whether a
    hit is REAL. It says nothing about whether identity measures the distance between a
    real contig and its real relative, which is the quantity this stage reports and which
    the holdout answers directly: a genome 30.28% diverged read 66.5% protein identity
    over 620 aa where a referenced comparator from the same library read 99.9%, ordered
    correctly and 30 points apart. The usable statement is the narrow one. A significant
    alignment's identity does not tell you the query is real, and on short alignments it
    is quantised and optimistic, so it ships with its coverage and its aligned length
    every time.

    What does separate is whether a significant alignment exists at all, which the null
    almost never produces. The direction of that bias matters: the null is a filter, not
    a shift. So this reports the existence rate first, the per-field comparison second
    with its honest reading, and a stricter operating point third, since the question a
    pipeline actually has to answer is whether any threshold on coverage and length
    empties the null while keeping the real hits.
    """
    base = {"readout": ("existence of a significant alignment, then every property of "
                        "that alignment, real against a dinucleotide-preserving shuffle "
                        "of the same contigs through the same search")}
    # The two arms can have different denominators, because a sequence invariant under a
    # dinucleotide shuffle is dropped from the null rather than passed through as itself.
    # Dividing both by the real arm's count would understate the null's rate.
    n_null = n_contigs if n_null is None else n_null
    rr = round(100.0 * len(sig) / max(1, n_contigs), 2)
    nr = round(100.0 * len(nullsig) / max(1, n_null), 2)
    base["existence"] = {
        "contigs_searched": n_contigs,
        "shuffles_searched": n_null,
        "real_with_a_significant_hit": len(sig),
        "real_percent": rr,
        "shuffled_with_a_significant_hit": len(nullsig),
        "shuffled_percent": nr,
        "bound_on_chance": (f"0 of {n_null} bounds chance at roughly 1 in {n_null}; it "
                            f"does not establish zero" if not nullsig else None),
        "ratio": round(rr / nr, 1) if nr else None,
        "reading": (f"{rr}% of real contigs against {nr}% of shuffles. The existence of "
                    f"the alignment is the calibrated quantity"),
    }
    fields = {}
    best = None
    for f in FIELDS:
        rv = [h[f] for h in sig.values() if h.get(f) is not None]
        nv = [h[f] for h in nullsig.values() if h.get(f) is not None]
        q, nq = quantiles(rv), quantiles(nv)
        share = None
        if rv and nv:
            share = round(100.0 * sum(1 for x in rv if x > nq["p90"]) / len(rv), 2)
        fields[f] = {"real": q, "shuffled_null": nq,
                     "percent_of_real_above_null_p90": share}
        if share is not None and (best is None or share > best[1]):
            best = (f, share)
    base["conditional_on_a_hit"] = {
        "readout": ("share of real hits above the shuffled null's p90, per field. A null "
                    "p90 puts 10% above it by construction, so 10% is the no-information "
                    "level and the whole table has to be read against that number"),
        "no_information_level_percent": 10.0,
        "fields": fields,
    }
    if best:
        base["conditional_on_a_hit"]["strongest_field"] = best[0]
        pid = fields["pident"]["percent_of_real_above_null_p90"]
        base["conditional_on_a_hit"]["reading"] = (
            f"cannot discriminate here. The strongest field is {best[0]} at {best[1]}% "
            f"against a no-information level of 10%, and identity itself reads {pid}%. "
            f"Once a shuffle clears the significance bar its alignment looks like a real "
            f"one on every field measured, so no property of a hit may be used to argue "
            f"that the hit is real")
    floors = []
    for cov, aln in ((0.5, 100), (0.8, 150), (0.9, 200)):
        keep = lambda d: sum(1 for h in d.values()
                             if h["subject_coverage"] >= cov and h["aln_aa"] >= aln)
        kr, kn = keep(sig), keep(nullsig)
        floors.append({
            "subject_coverage_at_least": cov, "aln_aa_at_least": aln,
            "real": kr, "shuffled_null": kn,
            "real_percent_of_searched": round(100.0 * kr / max(1, n_contigs), 2),
            "shuffled_percent_of_searched": round(100.0 * kn / max(1, n_contigs), 2),
        })
    base["at_a_stricter_operating_point"] = {
        "readout": ("the same two piles under coverage and length floors. This is the "
                    "only lever that moves the null, and it costs real hits, so the "
                    "trade is reported rather than chosen here"),
        "floors": floors,
    }
    return base


def _union(ivs) -> int:
    out, cur = 0, None
    for lo, hi in sorted(ivs):
        if cur is None or lo > cur[1]:
            if cur:
                out += cur[1] - cur[0] + 1
            cur = [lo, hi]
        else:
            cur[1] = max(cur[1], hi)
    if cur:
        out += cur[1] - cur[0] + 1
    return out


def coverage_gap(hsps: dict, qlen: dict, sig: dict,
                 nhsps: dict, nqlen: dict, nullsig: dict) -> dict:
    """How much more of a contig is explained by SEVERAL known proteins than by one.

    Two coverages of the same contig by the same significant HSPs. Per-subject is the
    union against whichever single subject covers most of it; any-subject is the union
    against all of them. The gap between them is the part of the contig that known
    proteins explain only collectively, which happens when a contig's genes land on
    unrelated database records.

    The holdout thread measured this on a library where the answer is known: the withheld
    genome reads +36.6 points, the 95.85th percentile among 193 length-matched identified
    contigs whose gap is median 0.1, and 34 contigs from already-referenced genomes in
    the same library read median, mean and MAXIMUM 0.0. They also killed the obvious
    confound: binned by number of significant subjects the median gap is 0.0 in every
    bin, the candidate's own bin of 5 subjects included.

    What that measurement lacked is a null, which is the arm this script has. Many weak
    HSPs scattered over a shuffled contig could manufacture a gap by chance, and if they
    do, the field is reading chance alignment density rather than gene architecture. So
    the shuffled arm is computed identically and reported beside the real one.
    """
    def gaps(hs, ql, keep):
        out = []
        for q, subj in hs.items():
            if q not in keep:
                continue
            n = max(1, ql.get(q, 1))
            per = max((_union(iv) for iv in subj.values()), default=0) / n
            allu = _union([iv for ivs in subj.values() for iv in ivs]) / n
            out.append({"contig": q, "subjects": len(subj),
                        "per_subject_coverage": round(per, 4),
                        "any_subject_coverage": round(allu, 4),
                        "gap_points": round(100.0 * (allu - per), 2)})
        return out

    real, null = gaps(hsps, qlen, sig), gaps(nhsps, nqlen, nullsig)

    def prof(rs):
        if not rs:
            return {}
        v = sorted(r["gap_points"] for r in rs)
        at = lambda f: v[min(len(v) - 1, int(f * (len(v) - 1) + 0.5))]
        by_n = collections.defaultdict(list)
        for r in rs:
            by_n[min(r["subjects"], 6)].append(r["gap_points"])
        return {"contigs": len(rs),
                "median_gap_points": round(at(0.50), 2),
                "p90_gap_points": round(at(0.90), 2),
                "p99_gap_points": round(at(0.99), 2),
                "max_gap_points": round(v[-1], 2),
                "share_above_10_points": round(sum(1 for x in v if x > 10) / len(v), 4),
                "share_above_30_points": round(sum(1 for x in v if x > 30) / len(v), 4),
                "median_by_significant_subject_count": {
                    str(k): round(statistics.median(x), 2)
                    for k, x in sorted(by_n.items())},
                }
    return {"readout": ("any-subject coverage minus per-subject coverage, in points of "
                        "contig length. High means several unrelated known proteins each "
                        "explain a piece and none explains the whole"),
            "real": prof(real), "shuffled_null": prof(null),
            "holdout_reference_point": (
                "the withheld genome reads +36.6 points on the holdout library, the "
                "95.85th percentile of its length-matched background, where 34 contigs "
                "from already-referenced genomes read a maximum of 0.0")}


COV_BANDS = [(0.0, 0.25), (0.25, 0.50), (0.50, 0.70), (0.70, 1.01)]
ID_BANDS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100.1)]


def joint_cells(sig: dict, nullsig: dict) -> dict:
    """Contigs per cell of contig coverage against identity, with the shuffled null.

    THE CELL IS THE OUTPUT, NOT A RANKING. Two independent measurements say so. The
    shuffled null here found coverage to be the only field reading above chance, identity
    the weakest. And the holdout library, where the answer is known, puts its withheld
    genome at the 93.46th percentile of identity among removed contigs - 7,104 contigs
    are MORE divergent than the one thing that mattered - while on coverage it lands in a
    cell holding 70 contigs, 0.92% of that pile. In the same library not one of 34 contigs
    from already-referenced genomes falls below 80% identity at either level, so
    recognisable material does not enter that cell: its occupants are other unnamed
    things. The withheld genome is still not the top row of its own cell, which is why
    this reports a cell with a count on it and no row ordering at all.
    """
    def cell(d):
        g = collections.defaultdict(int)
        for h in d.values():
            for ci, (clo, chi) in enumerate(COV_BANDS):
                if clo <= h["contig_coverage"] < chi:
                    break
            else:
                ci = len(COV_BANDS) - 1
            for ii, (ilo, ihi) in enumerate(ID_BANDS):
                if ilo <= h["pident"] < ihi:
                    break
            else:
                ii = len(ID_BANDS) - 1
            g[(ci, ii)] += 1
        return g
    rg, ng = cell(sig), cell(nullsig)
    tot = max(1, len(sig))
    rows = []
    for ci, (clo, chi) in enumerate(COV_BANDS):
        rows.append({
            "contig_coverage_from": clo, "contig_coverage_to": round(chi, 2),
            "cells": [{"identity_from": ilo, "identity_to": round(ihi, 1),
                       "contigs": rg.get((ci, ii), 0),
                       "percent_of_pile": round(100.0 * rg.get((ci, ii), 0) / tot, 2),
                       "shuffled_null": ng.get((ci, ii), 0)}
                      for ii, (ilo, ihi) in enumerate(ID_BANDS)],
        })
    return {"readout": ("contigs per cell of contig coverage against identity. A cell "
                        "with a count is the reportable object here; ordering rows "
                        "inside a cell by identity would order them by database "
                        "proximity"),
            "denominator": len(sig),
            "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=Path,
                    default=ROOT / "data/geo/sample_contigs_identified_norrna.fa",
                    help="the identified pile. rRNA contigs are excluded because a "
                         "ribosomal match says what the sequence is, not whose, and a "
                         "divergent rRNA is a different question with its own instrument")
    ap.add_argument("--min-contig", type=int, default=300,
                    help="below 300 nt a 60 aa alignment cannot be most of the contig")
    ap.add_argument("--divergent-below", type=float, default=90.0)
    ap.add_argument("--control-above", type=float, default=95.0)
    ap.add_argument("--min-orf-aa", type=int, default=50)
    ap.add_argument("--gene-level-max", type=int, default=600,
                    help="cap on contigs taken to the gene level, deepest first; the "
                         "whole-contig distance is reported for all of them regardless")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    library_facts.add_argument(ap)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="mate 2 bases trimmed when mapping reads back. This stage "
                         "trimmed 16 where the measured prefix is 14, so the value "
                         "comes from trim_r2_for_mapback in "
                         "experiments/libraries.json, which pins what these results "
                         "were produced with and states the discrepancy. Reconciling "
                         "it to 14 means re-running this stage")
    ap.add_argument("--map-limit", type=int, default=0)
    ap.add_argument("--dmnd", type=Path, default=ROOT / "data/prot/viral.dmnd")
    ap.add_argument("--dmnd-sprot", type=Path, default=ROOT / "data/prot/sprot.dmnd")
    ap.add_argument("--pfam", type=Path, default=ROOT / "data/pfam/Pfam-A.hmm")
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-mode", default="--sensitive")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    # The map-back trim was a bare default of 16 here while the measured prefix is 14.
    # It is resolved from experiments/libraries.json under its own key, so the value these
    # results were produced with is pinned and the discrepancy is recorded rather than
    # silently reconciled; changing it means re-running this stage.
    library_facts.resolve_args(a, {"trim_r2": "trim_r2_for_mapback"})

    # --dmnd-sprot was missing from this check and then skipped with `continue` wherever
    # it was absent, so a run without it reported Swiss-Prot-free distances under the
    # same name. Every index named is required.
    for p in (a.contigs, a.dmnd, a.dmnd_sprot, a.pfam, a.r1, a.r2):
        if not p.exists() or (p.is_file() and p.stat().st_size == 0):
            raise SystemExit(f"missing or empty input: {p}")
    if not os.access(a.diamond_bin, os.X_OK):
        raise SystemExit(f"--diamond-bin {a.diamond_bin} is not executable")
    a.work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rng = random.Random(a.seed)

    allseqs = read_fasta(a.contigs)
    seqs = {n: s for n, s in allseqs.items() if len(s) >= a.min_contig}
    print(f"{len(seqs)} identified contigs >= {a.min_contig} nt "
          f"of {len(allseqs)}", file=sys.stderr)
    fa = a.work / "identified.fa"
    fa.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))

    # The null goes through the identical search, so the two cannot differ by anything but
    # the sequence. Same names with a suffix, so a row can be traced back to its original.
    shuf = a.work / "identified_shuffled.fa"
    sh = shuffle_all(seqs, rng)
    shuffled = sh["shuffled"]
    # Counted into the output, not into a comment. An arm that silently lost members and an
    # arm that silently copied the test produce the same plausible-looking distribution, so
    # the denominator the null's percentage is computed against has to be visible.
    null_census = {"sequences_offered": len(seqs),
                   "sequences_shuffled": len(shuffled),
                   "unchanged_and_dropped": len(sh["unchanged_and_dropped"]),
                   "unchanged_names": sh["unchanged_and_dropped"][:20],
                   "returned_input": 0,
                   "note": ("shuffle_control raises rather than falling back, so "
                            "returned_input is 0 by construction and is reported anyway "
                            "because it is the field three broken copies of this null "
                            "would have had to show")}
    print(f"null: {len(shuffled)} of {len(seqs)} shuffled, "
          f"{len(sh['unchanged_and_dropped'])} invariant and dropped", file=sys.stderr)
    shuf.write_text("".join(f">{n}\n{s}\n" for n, s in shuffled.items()))

    # HSPs are pooled ACROSS the two indexes, because a contig whose genes land on a
    # RefSeq viral record and a Swiss-Prot record is exactly the case the coverage gap
    # is there to catch, and keeping the indexes apart would hide it.
    hits, secs = {}, {}
    hsps, qlen = {}, {}
    nullhits, nhsps, nqlen = {}, {}, {}
    for label, db in (("viral", a.dmnd), ("sprot", a.dmnd_sprot)):
        r = diamond_blastx(fa, db, a.work / f"whole_{label}.tsv", a.diamond_bin,
                           a.threads, a.diamond_mode)
        secs[f"blastx_{label}"] = r["seconds"]
        for q, row in r["best"].items():
            row["index"] = label
            if q not in hits or row["bitscore"] > hits[q]["bitscore"]:
                hits[q] = row
        for q, d in r["hsps"].items():
            for subj, iv in d.items():
                hsps.setdefault(q, {}).setdefault(f"{label}|{subj}", []).extend(iv)
        qlen.update(r["qlen"])
        rn = diamond_blastx(shuf, db, a.work / f"null_{label}.tsv", a.diamond_bin,
                            a.threads, a.diamond_mode)
        secs[f"blastx_null_{label}"] = rn["seconds"]
        for q, row in rn["best"].items():
            if q not in nullhits or row["bitscore"] > nullhits[q]["bitscore"]:
                nullhits[q] = row
        for q, d in rn["hsps"].items():
            for subj, iv in d.items():
                nhsps.setdefault(q, {}).setdefault(f"{label}|{subj}", []).extend(iv)
        nqlen.update(rn["qlen"])

    def significant(h):
        return h and h["evalue"] < SIGNIFICANT_EVALUE and h["aln_aa"] >= MIN_ALN_AA

    sig = {q: h for q, h in hits.items() if significant(h)}
    nullsig = {q: h for q, h in nullhits.items() if significant(h)}
    ids = [h["pident"] for h in sig.values()]
    nullids = [h["pident"] for h in nullsig.values()]

    divergent = {q: h for q, h in sig.items() if h["pident"] < a.divergent_below}
    control = {q: h for q, h in sig.items() if h["pident"] >= a.control_above}
    print(f"{len(sig)} with a significant protein hit, {len(divergent)} below "
          f"{a.divergent_below}%", file=sys.stderr)

    # Gene level on the divergent arm and a matched slice of the control. Capped and
    # ordered by contig length, because a longer contig carries more genes and the cap has
    # to be a stated rule rather than whatever the dictionary order was.
    take = sorted(divergent, key=lambda q: -len(seqs[q]))[:a.gene_level_max]
    ctrl_take = sorted(control, key=lambda q: -len(seqs[q]))[:a.gene_level_max]
    gl = {n: seqs[n] for n in take + ctrl_take}
    sub = a.work / "gene_level.fa"
    sub.write_text("".join(f">{n}\n{s}\n" for n, s in gl.items()))

    orfs = call_orfs(gl, a.min_orf_aa)
    print(f"{len(orfs)} ORFs on {len(gl)} contigs", file=sys.stderr)
    orfhits = {}
    for label, db in (("viral", a.dmnd), ("sprot", a.dmnd_sprot)):
        r = diamond_blastp(orfs, db, a.work, a.diamond_bin, a.threads,
                           a.diamond_mode, f"gl_{label}")
        secs[f"blastp_{label}"] = r["seconds"]
        for q, row in r["best"].items():
            if q not in orfhits or row["bitscore"] > orfhits[q]["bitscore"]:
                orfhits[q] = row
    pf = pfam_scan(orfs, a.pfam, a.work, a.threads)
    secs["pfam"] = pf["seconds"]
    depth = map_back(sub, a.r1, a.r2, a.trim_r2, a.work, a.map_limit)

    by_contig = collections.defaultdict(list)
    for o in orfs:
        by_contig[o["contig"]].append(o)

    def row_for(name, arm):
        h = sig[name]
        mine = sorted(by_contig.get(name, []), key=lambda o: o["start"])
        genes, pol = [], []
        for o in mine:
            g = orfhits.get(o["orf"])
            if not (g and g["evalue"] < SIGNIFICANT_EVALUE and g["aln_aa"] >= MIN_ALN_AA):
                continue
            role = gene_role(g.get("title"))
            genes.append({"role": role, "pident": g["pident"], "aln_aa": g["aln_aa"],
                          "aa_len": o["aa_len"], "subject": g["subject"],
                          "title": g.get("title", "")[:160],
                          "pfam": [d["pfam"] for d in pf["hits"].get(o["orf"], [])]})
            if role == "polymerase":
                pol.append(g["pident"])
        gids = [g["pident"] for g in genes]
        return {
            "arm": arm,
            "contig": name,
            "contig_nt": len(seqs[name]),
            # The distance, with everything needed to read it.
            "pident": h["pident"],
            "aln_aa": h["aln_aa"],
            "evalue": h["evalue"],
            "subject": h["subject"],
            "subject_title": h["title"],
            "subject_aa": h["subject_aa"],
            "subject_coverage": h["subject_coverage"],
            "contig_coverage": h["contig_coverage"],
            "index": h.get("index"),
            "divergence_class": band(h["pident"]),
            "polymerase_pident": max(pol) if pol else None,
            "genes": genes,
            "orf_identity_spread_points": (round(max(gids) - min(gids), 1)
                                           if len(gids) >= 2 else None),
            "read_support": depth["per_contig"].get(name, {}),
        }

    rows = ([row_for(n, "divergent") for n in take]
            + [row_for(n, "control_near_reference") for n in ctrl_take])
    rows.sort(key=lambda r: (r["pident"], -r["contig_nt"]))

    def arm_summary(rs):
        if not rs:
            return {}
        sp = [r["orf_identity_spread_points"] for r in rs
              if r["orf_identity_spread_points"] is not None]
        return {
            "contigs": len(rs),
            "median_pident": round(statistics.median([r["pident"] for r in rs]), 1),
            "median_aln_aa": statistics.median([r["aln_aa"] for r in rs]),
            "median_subject_coverage":
                round(statistics.median([r["subject_coverage"] for r in rs]), 3),
            "contigs_with_two_or_more_scorable_genes": len(sp),
            "median_orf_identity_spread_points": (round(statistics.median(sp), 2)
                                                  if sp else None),
            "median_mapped_depth": statistics.median(
                [r["read_support"].get("mapped_mean_depth", 0) for r in rs]),
        }

    # The two gene-level arms are NOT comparable as drawn, and the confound runs one way.
    # On the full pile the divergent arm's contigs have a median of 3 scorable genes over
    # 3,820 nt where the control arm has 1 over 552 nt, because a contig matching its
    # reference at 95% is usually a short fragment of a well-sampled organism. More genes
    # and more length give a spread more room, so the raw comparison is biased TOWARDS
    # the divergent arm. Matching on the number of scorable genes removes that, and the
    # gap survives it, so the matched figure is the one to quote.
    def matched_spread(rs_d, rs_c) -> dict:
        def by_gene_count(rs):
            out = collections.defaultdict(list)
            for r in rs:
                n = len([g for g in r.get("genes", []) if g.get("pident") is not None])
                if r.get("orf_identity_spread_points") is not None:
                    out[n].append(r["orf_identity_spread_points"])
            return out
        dd, cc = by_gene_count(rs_d), by_gene_count(rs_c)
        bands = []
        for n in sorted(set(dd) & set(cc)):
            if n < 2 or len(dd[n]) < 5 or len(cc[n]) < 5:
                continue
            bands.append({
                "scorable_genes": n,
                "divergent_contigs": len(dd[n]),
                "divergent_median_spread_points": round(statistics.median(dd[n]), 2),
                "control_contigs": len(cc[n]),
                "control_median_spread_points": round(statistics.median(cc[n]), 2),
            })
        return {"readout": ("median per-contig gene identity spread at equal numbers of "
                            "scorable genes. The unmatched comparison is biased towards "
                            "the divergent arm, whose contigs are longer and carry more "
                            "genes"),
                "bands": bands}

    out = {
        "unit": "one row per identified contig taken to the gene level",
        "question": ("how far was each identified contig from the reference that "
                     "identified it, which the cascade computes and then discards"),
        "a_row_is": "a distance and a coverage. There is no score on this table.",
        "scope": {
            "input": str(a.contigs),
            "contigs_in_file": len(allseqs),
            "contigs_at_least_min_length": len(seqs),
            "min_contig_nt": a.min_contig,
            "with_a_significant_protein_hit": len(sig),
            "below_divergent_threshold": len(divergent),
            "at_or_above_control_threshold": len(control),
            "taken_to_gene_level": len(gl),
            "gene_level_cap_per_arm": a.gene_level_max,
        },
        "identity_distribution": {
            "readout": ("contigs per five-point identity bin, over every identified "
                        "contig with a significant protein hit. This is the shape that "
                        "says whether a line can be drawn or only a ranking"),
            "real": histogram(ids),
            "dinucleotide_shuffled_null": histogram(nullids),
            "null_reading": (f"{len(nullsig)} of {len(seqs)} shuffled contigs clear the "
                             f"same significance bar, against {len(sig)} real"),
            "median_real": round(statistics.median(ids), 1) if ids else None,
            "median_null": round(statistics.median(nullids), 1) if nullids else None,
        },
        "joint_cells": joint_cells(sig, nullsig),
        "coverage_gap": coverage_gap(hsps, qlen, sig, nhsps, nqlen, nullsig),
        "gene_count_matched_spread": matched_spread(
            [r for r in rows if r["arm"] == "divergent"],
            [r for r in rows if r["arm"] == "control_near_reference"]),
        "discrimination": discrimination(sig, nullsig, len(seqs), len(shuffled)),
        "null_census": null_census,
        "operating_point": {
            "min_alignment_aa": MIN_ALN_AA,
            "significant_evalue": SIGNIFICANT_EVALUE,
            "diamond_mode": a.diamond_mode,
            "divergent_below_percent": a.divergent_below,
            "control_at_or_above_percent": a.control_above,
            "min_orf_aa": a.min_orf_aa,
            "pfam_threshold": "curators' gathering cutoffs",
            "mate_2_trimmed_nt": a.trim_r2,
            "mates_actually_streamed": depth["mates_streamed"],
            "null": "dinucleotide-preserving shuffle of the same contigs, same search",
        },
        "bands": [{"at_least_percent": lo, "reading": name} for lo, name in BANDS],
        "summary": {"divergent": arm_summary([r for r in rows if r["arm"] == "divergent"]),
                    "control_near_reference": arm_summary(
                        [r for r in rows if r["arm"] != "divergent"])},
        "external_positive_control": {
            "what": ("the holdout thread's withheld genome: it assembles as one contig "
                     "covering 96% of the withheld reference and hits it at 66.5% "
                     "protein identity over 620 aa, and the cascade removes it as "
                     "identified in tier 1"),
            "scope": ("a different library, so this is an operating point and not a row "
                      "in this table. The older mutated-spike ladder is a reference "
                      "point rather than ground truth for a real distance."),
        },
        "blind_spots": [
            "whole-contig protein identity does not separate real sequence from a "
        "dinucleotide shuffle of it: the null's median identity is as low as the "
        "real one, because a chance alignment to a narrow index is short and lands "
        "wherever it lands. Read identity only beside alignment length and subject "
        "coverage, which is what the discrimination block reports",
        "6.1% of dinucleotide shuffles clear E<1e-5 at 60 aa, which is high for a "
        "null. The untested candidate is low-complexity sequence, which a shuffle "
        "preserves and which aligns to anything; the control is the aligner's own "
        "compositional masking, not run here. Until that is measured the null rate is "
        "an upper bound on chance and the 16x existence ratio a lower bound",
        "every identity here comes from RefSeq viral proteins plus Swiss-Prot, which "
            "are narrow, so each is an UPPER BOUND on distance from what is known",
            "this says nothing about whether a divergent contig is a new species: "
            "demarcation criteria are per family and nothing here fetches them",
            "protein identity understates nucleotide divergence because synonymous "
            "changes are invisible to it",
            "a contig is already a consensus over reads and hides within-sample strain "
            "mixtures below about 3.8% divergence",
            "the gene-level arm is capped by contig length, so a short divergent contig "
            "may carry a whole-contig distance and no gene-level row",
        ],
        "rows": rows,
        "seed": a.seed,
        "timings_seconds": dict(secs, total=round(time.time() - t0, 1)),
    }
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                             capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, check=True,
                                    capture_output=True, text=True).stdout.strip())
    except Exception:
        sha, dirty = "unknown", True
    out["git_sha"] = sha + ("-dirty" if dirty else "")
    out["tree_dirty"] = dirty
    out["script_sha256"] = sha256_of(Path(__file__).resolve())
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"wrote {a.out} in {out['timings_seconds']['total']}s", file=sys.stderr)


if __name__ == "__main__":
    main()
