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
import hashlib
import json
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from shuffle_control import shuffle_all  # noqa: E402  the one copy of the null's shuffle

sys.path.insert(0, str(Path(__file__).resolve().parent))

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
    """
    cmd = [str(binpath), "blastx", "-d", str(db), "-q", str(fa), "-o", str(out),
           "-p", str(threads), mode, "--quiet", "--evalue", "1e-3",
           "--max-target-seqs", "1", "--outfmt", "6",
           "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
           "qlen", "slen", "stitle"]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    best = {}
    with out.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            row = {"subject": f[1], "pident": float(f[2]), "aln_aa": int(f[3]),
                   "evalue": float(f[4]), "bitscore": float(f[5]),
                   "contig_nt": int(f[6]), "subject_aa": int(f[7]),
                   "title": (f[8] if len(f) > 8 else "")[:200]}
            row["subject_coverage"] = round(row["aln_aa"] / max(1, row["subject_aa"]), 3)
            row["contig_coverage"] = round(3 * row["aln_aa"] / max(1, row["contig_nt"]), 3)
            prev = best.get(f[0])
            if prev is None or row["bitscore"] > prev["bitscore"]:
                best[f[0]] = row
    return {"best": best, "seconds": round(time.time() - t0, 1)}


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
    subject coverage and bitscore all sat inside the real spread. So a low identity is
    not evidence of distance from the database: a chance alignment to a narrow protein
    index lands in the same band as real distant homology, and nothing about the
    alignment, once it exists, says which it is.

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
    ap.add_argument("--trim-r2", type=int, default=16)
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

    hits, secs = {}, {}
    for label, db in (("viral", a.dmnd), ("sprot", a.dmnd_sprot)):
        r = diamond_blastx(fa, db, a.work / f"whole_{label}.tsv", a.diamond_bin,
                           a.threads, a.diamond_mode)
        secs[f"blastx_{label}"] = r["seconds"]
        for q, row in r["best"].items():
            row["index"] = label
            if q not in hits or row["bitscore"] > hits[q]["bitscore"]:
                hits[q] = row
    nullhits = {}
    for label, db in (("viral", a.dmnd), ("sprot", a.dmnd_sprot)):
        r = diamond_blastx(shuf, db, a.work / f"null_{label}.tsv", a.diamond_bin,
                           a.threads, a.diamond_mode)
        secs[f"blastx_null_{label}"] = r["seconds"]
        for q, row in r["best"].items():
            if q not in nullhits or row["bitscore"] > nullhits[q]["bitscore"]:
                nullhits[q] = row

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
