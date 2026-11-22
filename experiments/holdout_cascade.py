#!/usr/bin/env python3
"""What the committed cascade's own output says about a genome that was genuinely unknown.

THE QUESTION, AND WHY IT IS NOT THE ONE PR #29 ANSWERED
-------------------------------------------------------
`results/holdout_run/` measured the read-level search path on this library: a 145 nt read
from a 30.3%-divergent virus clears E < 1e-5 about one time in sixteen, and the identity it
reports flatters the truth by 13 points unless coverage is read beside it. That is the
instrument's resolution on a single read.

It is not what a person running the tool would see. `hapax run` assembles, subtracts
everything identifiable, routes what assembled-but-could-not-be-identified to tier 2, and
reports the residue in layers. A 4,147 nt genome at 230x does not reach a reader as 13,518
independent reads; it reaches them as a contig, in whichever pile the cascade put it, at
whatever rank the ordering gives it. So the question here is:

  **given the cascade's own output, with the clock turned back, is this genome anywhere a
  reader's eye would land - and what number beside it would say so?**

That is answerable from committed files and has a right answer nobody chose.

WHAT COUNTS AS AN ANSWER, AND WHAT WOULD NOT
--------------------------------------------
Not "did it get flagged". This project has measured, repeatedly, that an index can answer
*is this already known* and cannot answer *is this new*, because score ordering tracks
database proximity. A novelty score here would reintroduce the exact failure the work is
built around, and `hapax.py`'s own `strains` docstring refuses to build one. So what is
reported is the shape SPEC promises: an identity to the closest relative with its coverage
and its E-value, the tier that placed or failed to place the contig, and the contig's RANK
in each ordering the tool already produces. Whether that is enough to draw an eye is then a
statement about the output, not a threshold inside it.

THREE ARMS, SCORED BY THE SAME CODE IN THE SAME PASS
----------------------------------------------------
A row for the candidate alone says nothing. Every number below is computed by one function
over three populations:

- **candidate**: contigs of the genome deposited 2023-07-31, out of a library collected
  2022-05-10. Unknown to every database on earth when this water was sequenced.
- **leak**: contigs of a genome from the SAME library whose nearest relative is 99.9%
  identical and long predates the cutoff. This is the already-in-reference comparator, and
  it is also the alarm: if the candidate's row reads like this one, the holdout leaked; if
  this one does NOT read as identified, the holdout has removed more than it should and
  every other row is low for the wrong reason.
- **background**: every other contig in the assembly. Ordinary library material, which is
  what "already seen" looks like here.

THE STRATIFICATION THAT IS NOT OPTIONAL
---------------------------------------
Identity, coverage and identifiability all depend strongly on contig length - this project
measured protein identifying 26% of contigs against 6% of reads for that reason alone - and
the candidate's contigs are long because it is abundant. Comparing its row against the
background's pooled median would therefore compare lengths, not divergences. Background is
reported in length bins and the candidate's bin is named beside its row. The repo's rule is
to stratify rather than match when the matching variable is part of the question, and
length is: a novel thing that assembles well is long BECAUSE it is abundant.

THE DRIFT CONTROL
-----------------
Every distance is measured twice, against the date-held-out index and against the current
one, in the same pass. Two things make that necessary. The reference can only have grown,
so a distance against today's index is a floor and the pair of numbers shows how much of
the measured distance is an artefact of deposition date. And the rebuild of the held-out
index (`holdout_indexes.py`) could not reproduce PR #29's nucleotide record count exactly,
because upstream RefSeq dropped nine records in a day; the two columns are how that is
measured rather than argued about.

UNITS, STATED ONCE
------------------
A row is a CONTIG, not a read and not a pair. Coverage is the fraction of the CONTIG the
alignment covers, not of the subject, unless the column says subject. Identity is never
reported without its E-value and its coverage in the same row, because this repo has had to
retract an identity quoted without them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from search_hits import (  # noqa: E402  the one copy of the search-output parser
    DEFAULT_COLUMNS, SIGNIFICANCE_EVALUE, aligned_query_nt, best_per_query, coverage_gap,
    outfmt_args, query_coverage, read_hits, union_query_coverage)

ROOT = Path(__file__).resolve().parents[1]

#: Minimum fraction of a contig that must align to a reference genome for the contig to be
#: attributed to it. A contig is attributed by where its sequence comes from, not by its
#: best hit: a chimaeric or repeat-sharing contig would otherwise be assigned to whichever
#: genome won by a few bases. Reported per contig as well, so the threshold is visible.
ATTRIBUTE_MIN_COVERAGE = 0.50
#: What each of this script's two searches is, in `search_hits`' names. The name carries
#: the units: `diamond_blastx` has a nucleotide query and an amino-acid alignment length,
#: and every coverage computed from one of its rows knows that without being told again.
PROGRAM = {"blastn": "blastn", "diamond": "diamond_blastx"}
LENGTH_BINS = [(200, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 10 ** 9)]


def log(*a):
    print(*a, flush=True)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def diamond(binary: Path, query: Path, db: str, out: Path, threads: int,
            sensitivity: str = "--sensitive", evalue: str = "10") -> Path:
    """Translated search with DIAMOND, in the columns the rest of this script needs.

    DIAMOND rather than blastx, for two reasons and the second is the one that matters.
    Speed: blastx of 22 Mb of contigs against 700,000 proteins is hours, and this has to
    run over the whole assembly rather than a sample, because the arms must go through one
    invocation of one index. And identity: the CASCADE's protein tier is DIAMOND at
    `--sensitive` with a bitscore floor of 50, so a removal identity measured with blastx
    would be a second instrument's reading of the same question - and this repo's rule is
    never to compare numbers from two instruments. The point of these columns is to say
    what the cascade knew when it removed a contig, so they are computed with what the
    cascade used.

    The column list is `search_hits.DEFAULT_COLUMNS`, passed to the command and to the
    parser from the same place, so the two cannot disagree about which position holds the
    E-value. `qlen` and `slen` are in it rather than derived: a coverage computed against
    a length read from somewhere else is how a FASTA-parsing off-by-one becomes a claim.
    """
    cmd = [str(binary), "blastx", "--db", db, "--query", str(query), "--out", str(out),
           "--threads", str(threads), *outfmt_args(DEFAULT_COLUMNS, diamond=True),
           "--max-target-seqs", "5", "--evalue", evalue, sensitivity, "--quiet"]
    subprocess.run(cmd, check=True)
    return out


def blast(binary: Path, task: str, query: Path, db: str, out: Path, threads: int,
          evalue: str = "10") -> Path:
    """One BLAST call, tabular, with the columns every row in this script needs.

    The columns come from `search_hits.DEFAULT_COLUMNS`, which is also what the parser is
    given, so the command and the unpacking share one list. `qlen` and `slen` are fetched
    rather than derived so that coverage is computed from what BLAST itself believed the
    lengths were. A coverage computed against a length read from somewhere else is how an
    off-by-one in a FASTA parser becomes a biological claim.
    """
    cmd = [str(binary), "-query", str(query), "-db", db, "-out", str(out),
           *outfmt_args(DEFAULT_COLUMNS, diamond=False),
           "-evalue", evalue, "-num_threads", str(threads), "-max_target_seqs", "5"]
    if task == "blastn":
        cmd[0] = str(binary)
        cmd += ["-task", "blastn"]
    subprocess.run(cmd, check=True)
    return out


def parse_best(rows: list[dict]) -> dict[str, dict]:
    """Best HSP per query by bitscore, in the field names this script's rows use.

    The parsing, the tie rule and the units now come from `search_hits`; this is the
    rename into local field names and nothing else. It is a wrapper rather than a copy
    because of what the copy cost: for DIAMOND blastx the query is nucleotide and
    `length` is residues, so this function used to divide one by the other and every
    protein best-HSP coverage capped at exactly 0.3333, which reads as a plausible low
    coverage rather than as a bug. `search_hits.query_coverage` takes the span from
    `qstart`/`qend`, which are query nucleotides under every program here, so the unit
    cannot be got wrong and gaps are not counted as alignment. That makes the field
    slightly different from the `3 * aln_aa / qlen` the previous version computed, and
    the coordinate version is the one that is true; the union coverages, which are what
    the report quotes, were always built from the same coordinates and do not move.

    Best-by-bitscore rather than best-by-identity: identity on a short alignment is
    quantised and rises as the alignment shortens, so ordering by it systematically
    prefers the least informative hit. This repo has the retraction to prove it.
    """
    out: dict[str, dict] = {}
    for q, r in best_per_query(rows).items():
        out[q] = {"subject": r["sseqid"], "pident": r["pident"], "aln_len": r["length"],
                  "evalue": r["evalue"], "bitscore": r["bitscore"],
                  "qlen": r["qlen"], "slen": r["slen"],
                  "qstart": r["qstart"], "qend": r["qend"],
                  "aligned_query_nt": aligned_query_nt(r),
                  "query_coverage": round(query_coverage(r), 4)}
    return out


def per_subject_union(rows: list[dict]) -> dict[str, dict]:
    """Per query, the union of SIGNIFICANT HSP intervals on the best-covering subject.

    A single HSP cannot cover a multi-gene genome - this project already demoted exactly
    the right candidates with a coverage flag defined on one HSP, because a levivirus's
    best blastx HSP covers about 45% of it. Attribution therefore uses the union of
    intervals, so a contig re-derived in several pieces counts as covered.

    THE E-VALUE FLOOR IS NOT OPTIONAL AND IT IS NO LONGER WRITTEN HERE. Two versions of
    this union and the any-subject one lived side by side in this file, and both times
    they drifted on the filter rather than on the definition: first a missing floor on one
    side (394 of 7,734 protein rows and 508 of 2,269 nucleotide rows with a negative gap,
    to -0.50), then `>` against `>=` on the boundary (one row of 41,271, at -0.0453). The
    floor is now one `passes_floor` inside `search_hits`, passed one `evalue_below` by
    both modes, so a mismatch is no longer expressible. Attribution is unaffected either
    way: that pass runs at `evalue="1e-10"`, so every HSP it sees already clears the
    floor.
    """
    out = {}
    for q, u in union_query_coverage(rows, per_subject=True,
                                     evalue_below=SIGNIFICANCE_EVALUE).items():
        out[q] = {"subject": u["subject"], "covered_bases": u["covered_bases"],
                  "covered_fraction": round(u["covered_fraction"], 4),
                  "hsps": u["blocks"], "qlen": u["qlen"]}
    return out


def attribution_ties(rows: list[dict]) -> dict[str, dict]:
    """Per query, which reference genomes tie for the best union coverage, if any.

    WHY THIS IS A FIELD AND NOT A FOOTNOTE. `per_subject_union` keeps the single
    best-covering subject, and when two subjects cover exactly the same bases it keeps
    whichever the file happened to list first. Switching this script onto the shared union
    changed that order and relabelled 12 of 35,328 contigs from one reference to the
    other, at byte-identical covered fractions. Nothing moved - all 12 sit far below the
    0.50 attribution threshold and stay in `background` - but a label that flips when the
    iteration order changes is the array-order bias this repo has already withdrawn a
    result over, and the only safe version of it is the one that says so in the row.

    Built by splitting the rows by subject and unioning each subject's own intervals with
    the shared function, so there is no second copy of the union here and both subjects
    are measured under the same floor.
    """
    by_sub: dict[str, list[dict]] = {}
    for r in rows:
        by_sub.setdefault(r["sseqid"], []).append(r)
    per: dict[str, dict[str, int]] = {}
    for sub, rs in by_sub.items():
        for q, u in union_query_coverage(rs, per_subject=False,
                                         evalue_below=SIGNIFICANCE_EVALUE).items():
            per.setdefault(q, {})[sub] = u["covered_bases"]
    out = {}
    for q, d in per.items():
        top = max(d.values())
        tied = sorted(s for s, v in d.items() if v == top)
        out[q] = {"tied": len(tied) > 1, "tied_subjects": tied,
                  "covered_bases": top}
    return out


def gaps_of(rows: list[dict]) -> dict[str, dict]:
    """Both union definitions and the gap between them, from one call over one row set.

    WHY BOTH DEFINITIONS SHIP, and it is the ground-truth finding of this experiment. The
    per-subject union keeps whichever subject covers most of the query, which answers
    "does a single known protein explain this contig". That is the wrong question for a
    divergent multi-gene genome, whose genes have each drifted onto a *different* nearest
    relative: the withheld genome's RdRp and its second gene hit two different RefSeq
    records, so the per-subject figure caps at 0.4624 while known proteins between them
    explain 0.8286 of the contig.

    The gap between the two is a field in its own right, and unlike almost everything else
    this project has tried it has a direction fixed on known material: the 34 contigs from
    already-referenced genomes in this same library have a gap of 0.0 at median, mean AND
    max, while the withheld genome's is +36.6 points, the 95.85th percentile of its length
    band. It is not a gene-count artefact either - binned by number of significant
    subjects the median gap is 0.0 in every bin.

    `search_hits.coverage_gap` forms both unions from the same rows under the same floor
    and raises on a negative gap before any fraction is formed, which is the invariant this
    file used to assert at the row writer - too late to stop a committed table, and only
    because someone thought to look. What neither says is whether a shuffled contig
    produces a gap, which needs the dinucleotide shuffle arm that
    `experiments/divergent_known.py` owns.
    """
    return coverage_gap(rows, evalue_below=SIGNIFICANCE_EVALUE)


def length_bin(n: int) -> str:
    for lo, hi in LENGTH_BINS:
        if lo <= n < hi:
            return f"{lo}-{hi if hi < 10 ** 9 else 'inf'}"
    return "under_200"


def summarise(rows: list[dict], key: str) -> dict | None:
    """Median and range of a field over an arm, with n, or an explicit absence.

    A field that no row in the arm carries is reported as `not measured` rather than as a
    zero or a blank - the project has already shipped a table where five of six layers
    separated the residue from a control nobody had searched, and the fix was to write the
    truth into the file.
    """
    vals = [r[key] for r in rows if r.get(key) is not None]
    if not vals:
        return {"n": 0, "state": "not measured: no row in this arm carries this field"}
    return {"n": len(vals), "median": round(st.median(vals), 4),
            "min": round(min(vals), 4), "max": round(max(vals), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="carried for the repo convention. The contig attribution and the "
                         "searches here are deterministic given the assembly, so the seed "
                         "changes nothing and the manifest says so rather than implying a "
                         "spread this stage does not have")
    ap.add_argument("--contigs", type=Path, required=True)
    ap.add_argument("--refs", type=Path,
                    default=ROOT / "results/holdout_candidates/candidates.fa")
    ap.add_argument("--candidate", default="OQ722337.1")
    ap.add_argument("--leak", default="OQ722330.1")
    ap.add_argument("--holdout-nt-db", required=True)
    ap.add_argument("--holdout-prot-db", required=True)
    ap.add_argument("--full-nt-db", default="",
                    help="today's unfiltered index. Without it the drift control is "
                         "absent and the report must say so, not assume the holdout and "
                         "the current reference agree")
    ap.add_argument("--full-prot-db", default="")
    ap.add_argument("--blastbin", type=Path, required=True)
    ap.add_argument("--diamond-bin", type=Path, required=True)
    ap.add_argument("--sprot-db", default="",
                    help="the cascade's other protein tier. Not date-filtered; carried "
                         "as its own column because a contig it identifies is removed "
                         "regardless of what the viral index says")
    ap.add_argument("--classify", type=Path, default=None,
                    help="classify_contigs.py's contig_identification.json, for the "
                         "operating point and the sweeps it records")
    ap.add_argument("--identified-fa", type=Path, default=None,
                    help="the contigs classify_contigs.py called identified: tier 1's "
                         "subtrahend. Membership is read from the file the cascade "
                         "actually consumed rather than recomputed, so the tier column "
                         "cannot disagree with the run it describes")
    ap.add_argument("--unidentifiable-fa", type=Path, default=None,
                    help="the contigs it could not identify: tier 2's population, "
                         "observed and never subtracted")
    ap.add_argument("--run-dir", type=Path, default=None,
                    help="a hapax run directory, for the rank columns: where the "
                         "candidate's contigs sit in the orderings the tool actually "
                         "emits")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    timings = {}

    contigs = read_fasta(a.contigs)
    log(f"{len(contigs)} contigs, {sum(map(len, contigs.values())):,} nt")

    # ---- 1. attribute contigs to the two reference genomes -------------------
    refs = read_fasta(a.refs)
    for want in (a.candidate, a.leak):
        if want not in refs:
            log(f"FATAL: {want} is not in {a.refs}")
            return 2
    ref_fa = a.work / "arm_refs.fa"
    with open(ref_fa, "w") as w:
        for k in (a.candidate, a.leak):
            w.write(f">{k}\n{refs[k]}\n")
    subprocess.run([str(a.blastbin / "makeblastdb"), "-in", str(ref_fa),
                    "-dbtype", "nucl", "-out", str(a.work / "db_arm_refs")],
                   check=True, capture_output=True)
    t = time.perf_counter()
    att_tsv = blast(a.blastbin / "blastn", "blastn", a.contigs,
                    str(a.work / "db_arm_refs"), a.work / "attribute.tsv", a.threads,
                    evalue="1e-10")
    timings["attribute"] = round(time.perf_counter() - t, 1)
    att_rows = read_hits(att_tsv, program="blastn")
    att = per_subject_union(att_rows)
    att_best = parse_best(att_rows)
    att_tie = attribution_ties(att_rows)

    arms: dict[str, list[str]] = {"candidate": [], "leak": [], "background": []}
    for name in contigs:
        hit = att.get(name)
        if hit and hit["covered_fraction"] >= ATTRIBUTE_MIN_COVERAGE:
            arms["candidate" if hit["subject"] == a.candidate else "leak"].append(name)
        else:
            arms["background"].append(name)
    log(f"attribution at >={ATTRIBUTE_MIN_COVERAGE:.0%} of the contig: "
        + ", ".join(f"{k} {len(v)}" for k, v in arms.items()))

    # ---- 2. every contig scored against both indexes, same call --------------
    # Background is searched whole rather than sampled: the arms must go through one
    # invocation of one index, or the comparison between them is between two instruments.
    dbs = {"holdout_nt": (a.holdout_nt_db, "blastn"),
           "holdout_prot": (a.holdout_prot_db, "diamond")}
    if a.full_nt_db:
        dbs["full_nt"] = (a.full_nt_db, "blastn")
    if a.full_prot_db:
        dbs["full_prot"] = (a.full_prot_db, "diamond")
    if a.sprot_db:
        # Swiss-Prot is NOT date-filtered - UniProt does not publish a release that old
        # here - and it is a cascade tier, so a contig it identifies is removed whatever
        # the viral index says. Direction, one-sided: a post-cutoff Swiss-Prot can only
        # make a contig look MORE identifiable, never less, so including it cannot
        # manufacture novelty. It is carried as its own column rather than merged, so what
        # it alone contributes to a removal is readable instead of implied.
        dbs["sprot_undated"] = (a.sprot_db, "diamond")
    hits: dict[str, dict[str, dict]] = {}
    for label, (db, task) in dbs.items():
        t = time.perf_counter()
        if task == "diamond":
            tsv = diamond(a.diamond_bin, a.contigs, db,
                          a.work / f"hits_{label}.tsv", a.threads)
        else:
            tsv = blast(a.blastbin / task, task, a.contigs, db,
                        a.work / f"hits_{label}.tsv", a.threads)
        # One read of the file, one row set, and every derived number built from it, so
        # the best hit, the two unions and the gap cannot disagree about which HSPs the
        # search produced. `read_hits` raises on a missing file rather than returning no
        # hits, which is the distinction between "searched, zero hits" and "not searched"
        # that this repo has paid for twice.
        rows_l = read_hits(tsv, program=PROGRAM[task])
        hits[label] = parse_best(rows_l)
        hits[label + "__union"] = per_subject_union(rows_l)
        hits[label + "__gap"] = gaps_of(rows_l)
        timings[f"search_{label}"] = round(time.perf_counter() - t, 1)
        log(f"  {label}: {len(hits[label])} contigs with any hit "
            f"({timings[f'search_{label}']}s)")

    # ---- 3. the cascade's own tier for each contig ---------------------------
    # Tier membership comes from the two FASTAs the cascade was actually given, not from
    # a recomputation: a tier column derived independently could disagree with the run it
    # is describing, and this repo's rule is never to compare two instruments.
    tier_of: dict[str, str] = {}
    classify_meta = None
    if a.classify and a.classify.exists():
        d = json.loads(a.classify.read_text())
        classify_meta = {k: v for k, v in d.items()
                         if k in ("operating_point", "kraken2_confidence_sweep",
                                  "protein_bitscore_sweep", "identified_frac",
                                  "identified_bases_frac", "rrna_contig_frac",
                                  "tier2_min_len", "route_nucleotide_only",
                                  "route_protein_only", "files", "index_bytes")}
    for fa, tier in ((a.identified_fa, "identified"),
                     (a.unidentifiable_fa, "assembled_unidentifiable")):
        if fa and fa.exists():
            for nm in read_fasta(fa):
                tier_of[nm] = tier
    # A contig in neither file is one classify_contigs left out of tier 2 for being below
    # its minimum length. That is a third state, not a missing value, and it is named:
    # calling it unidentifiable would overstate the tier-2 population, which is the number
    # the whole three-way split turns on.
    below = ("below the tier-2 minimum length: not identified and not routed to tier 2, "
             "because no classifier in the 2024 benchmark operates at that length")
    if any(tier_of.values()):
        for nm in contigs:
            tier_of.setdefault(nm, below)

    # ---- 4. one row per contig ----------------------------------------------
    rows = []
    for arm, names in arms.items():
        for nm in names:
            seq = contigs[nm]
            row = {"contig": nm, "arm": arm, "length": len(seq),
                   "length_bin": length_bin(len(seq)),
                   "gc": round(sum(c in "GCgc" for c in seq) / max(1, len(seq)), 4),
                   "tier": tier_of.get(nm, "not measured: classify_contigs output absent"),
                   "attributed_to": (att.get(nm) or {}).get("subject"),
                   "attributed_covered_fraction":
                       (att.get(nm) or {}).get("covered_fraction"),
                   "attributed_identity_best_hsp":
                       (att_best.get(nm) or {}).get("pident"),
                   # True means both references cover exactly the same bases of this
                   # contig, so `attributed_to` is iteration order and not evidence.
                   "attributed_tie": (att_tie.get(nm) or {}).get("tied"),
                   }
            for label in dbs:
                b = hits[label].get(nm)
                u = hits[label + "__union"].get(nm)
                g = hits[label + "__gap"].get(nm)
                if b is None:
                    row[label] = {"state": "no hit at the reporting cutoff"}
                    row[f"{label}_identity"] = None
                    row[f"{label}_coverage"] = None
                    row[f"{label}_coverage_any"] = None
                    row[f"{label}_coverage_gap"] = None
                    row[f"{label}_significant_subjects"] = None
                    row[f"{label}_evalue"] = None
                    row[f"{label}_bitscore"] = None
                    row[f"{label}_significant"] = False
                    continue
                cu = round(g["per_subject"], 4) if g else None
                ca = round(g["any_subject"], 4) if g else None
                row[label] = {
                    "subject": b["subject"], "identity_percent": b["pident"],
                    "aligned": b["aln_len"], "query_coverage_best_hsp":
                        b["query_coverage"],
                    "query_coverage_union": cu,
                    "union_hsps": (u or {}).get("hsps"),
                    "query_coverage_union_any_subject": ca,
                    "significant_subjects": (g or {}).get("significant_subjects"),
                    "evalue": b["evalue"], "bitscore": b["bitscore"],
                    "significant": b["evalue"] < SIGNIFICANCE_EVALUE,
                }
                row[f"{label}_identity"] = b["pident"]
                row[f"{label}_coverage"] = cu
                row[f"{label}_coverage_any"] = ca
                # The gap, stored rather than left to be derived, so that a reader of the
                # TSV cannot accidentally compare the two definitions the wrong way round.
                # Taken from `search_hits.coverage_gap`, which formed both unions from the
                # same rows under one floor and already refused to return a negative one.
                # THE INVARIANT USED TO BE ASSERTED HERE AND THAT WAS THE WRONG PLACE: a
                # row writer sees the number after every chance to have got it wrong has
                # been taken, and both times this shipped broken the cause was upstream -
                # a missing E-value floor on one side (394 of 7,734 protein rows and 508
                # of 2,269 nucleotide rows negative, to -0.50), then `>` against `>=` on
                # the boundary (one row of 41,271, at -0.0453). An impossible sign is the
                # cheapest check there is, and it belongs where the subtraction happens.
                row[f"{label}_coverage_gap"] = round(g["gap"], 4) if g else None
                row[f"{label}_significant_subjects"] = (
                    g or {}).get("significant_subjects")
                row[f"{label}_evalue"] = b["evalue"]
                # Bitscore as well as E-value, because the two answer different questions
                # and only together say whether a removal was marginal. The cascade's
                # protein tier fires at bitscore 50; a contig removed at 840 is not near
                # its threshold, and no tightening of that threshold could have spared it
                # without also sparing almost every genuine removal. That is checkable
                # only if the bitscore is in the row.
                row[f"{label}_bitscore"] = b["bitscore"]
                row[f"{label}_significant"] = b["evalue"] < SIGNIFICANCE_EVALUE
            rows.append(row)

    # ---- 5. ranks in the orderings the tool actually emits -------------------
    # "Would a reader's eye land on it" is a question about ORDER, so it is answered with
    # ranks rather than with an adjective. Each ordering is one the tool already produces.
    by_len = sorted(rows, key=lambda r: -r["length"])
    rank_len = {r["contig"]: i + 1 for i, r in enumerate(by_len)}
    unid = [r for r in rows if r["tier"] == "assembled_unidentifiable"]
    rank_unid_len = {r["contig"]: i + 1
                     for i, r in enumerate(sorted(unid, key=lambda r: -r["length"]))}
    for r in rows:
        r["rank_by_length_all_contigs"] = rank_len[r["contig"]]
        r["rank_by_length_within_tier2"] = rank_unid_len.get(r["contig"])

    # ---- 6. arm summaries, stratified by length -----------------------------
    summary = {}
    for arm in arms:
        ar = [r for r in rows if r["arm"] == arm]
        s = {"contigs": len(ar), "total_nt": sum(r["length"] for r in ar),
             "length": summarise(ar, "length"),
             "tier_counts": {}}
        for r in ar:
            s["tier_counts"][r["tier"]] = s["tier_counts"].get(r["tier"], 0) + 1
        for label in dbs:
            sig = [r for r in ar if r[f"{label}_significant"]]
            s[label] = {
                "contigs_with_significant_hit": len(sig),
                "fraction_significant": round(len(sig) / len(ar), 4) if ar else None,
                "identity_of_significant": summarise(sig, f"{label}_identity"),
                "coverage_of_significant": summarise(sig, f"{label}_coverage"),
                "note": ("identity and coverage are over SIGNIFICANT hits only. Identity "
                         "on a sub-significant alignment describes the database's size, "
                         "not the contig"),
            }
        s["by_length_bin"] = {}
        for lo, hi in LENGTH_BINS:
            b = f"{lo}-{hi if hi < 10 ** 9 else 'inf'}"
            br = [r for r in ar if r["length_bin"] == b]
            if not br:
                continue
            s["by_length_bin"][b] = {
                "contigs": len(br),
                **{label: {"fraction_significant":
                           round(sum(r[f"{label}_significant"] for r in br) / len(br), 4),
                           "identity_of_significant":
                               summarise([r for r in br if r[f"{label}_significant"]],
                                         f"{label}_identity")}
                   for label in dbs},
            }
        summary[arm] = s

    out = {
        "unit": "one row per CONTIG. Not reads, not pairs",
        "attribute_min_coverage": ATTRIBUTE_MIN_COVERAGE,
        "significance_evalue": SIGNIFICANCE_EVALUE,
        "arms": {"candidate": a.candidate, "leak": a.leak,
                 "background": "every contig attributed to neither"},
        "indexes": {k: v[0] for k, v in dbs.items()},
        "drift_control": ("every distance is measured against the date-held-out index and "
                          "today's index in the same pass. The reference can only have "
                          "grown, so a distance against today's is a floor"
                          if a.full_nt_db else
                          "NOT MEASURED: no current index was given, so how much of each "
                          "distance is an artefact of deposition date is unknown here"),
        "classify_contigs": classify_meta or
                            "not measured: no classify_contigs output given, so the tier "
                            "column is absent rather than guessed",
        "summary": summary,
        "candidate_rows": [r for r in rows if r["arm"] == "candidate"],
        "leak_rows": [r for r in rows if r["arm"] == "leak"],
        "not_a_score": ("nothing here orders contigs by novelty. Identity is a distance "
                        "from a database and ordering by it orders by database proximity, "
                        "which this project has measured and refuses to do. The ranks "
                        "reported are by LENGTH, which is the ordering the tool already "
                        "emits, and they are a statement about what a reader would see"),
    }
    (a.out / f"holdout_cascade_s{a.seed}.json").write_text(json.dumps(out, indent=1))

    import csv
    with open(a.out / f"contig_rows_s{a.seed}.tsv", "w", newline="") as fh:
        cols = ["contig", "arm", "length", "length_bin", "gc", "tier",
                "attributed_to", "attributed_covered_fraction",
                "attributed_identity_best_hsp", "attributed_tie",
                "rank_by_length_all_contigs", "rank_by_length_within_tier2"]
        for label in dbs:
            cols += [f"{label}_identity", f"{label}_coverage",
                     f"{label}_coverage_any", f"{label}_coverage_gap",
                     f"{label}_significant_subjects", f"{label}_evalue",
                     f"{label}_bitscore", f"{label}_significant"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["arm"] != "candidate",
                                             r["arm"] != "leak", -r["length"])):
            w.writerow(r)

    timings["total"] = round(time.perf_counter() - t0, 1)
    man = {
        "script": "holdout_cascade.py", "script_sha256": sha256(Path(__file__)),
        "git_sha": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                  capture_output=True, text=True).stdout.strip(),
        "tree_dirty": bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                          capture_output=True, text=True).stdout.strip()),
        "seed": a.seed,
        "seed_note": ("deterministic given the assembly: the attribution and both "
                      "searches have no randomness, so five seeds here would be five "
                      "identical files and none is run"),
        "parameters": {k: str(v) for k, v in vars(a).items()},
        "inputs": {p.name: {"bytes": p.stat().st_size, "sha256_16": sha256(p)[:16]}
                   for p in (a.contigs, a.refs) if p.exists()},
        "contigs": len(contigs),
        "arm_sizes": {k: len(v) for k, v in arms.items()},
        # How many contigs the two references cover equally, where the attribution label
        # is iteration order. Counted rather than described: it is 12 of 35,328 here, all
        # below the 0.50 threshold, so none of them changes an arm.
        "attribution_ties": sum(1 for v in att_tie.values() if v["tied"]),
        "attribution_ties_above_threshold": sum(
            1 for nm, v in att_tie.items()
            if v["tied"] and (att.get(nm) or {}).get("covered_fraction", 0)
            >= ATTRIBUTE_MIN_COVERAGE),
        "timings_s": timings,
        "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python": platform.python_version(), "platform": platform.platform(),
    }
    (a.out / f"manifest_s{a.seed}.json").write_text(json.dumps(man, indent=1))
    log(f"\n{a.out}/holdout_cascade_s{a.seed}.json")
    for arm in ("candidate", "leak"):
        for r in sorted((x for x in rows if x["arm"] == arm),
                        key=lambda r: -r["length"])[:4]:
            log(f"  [{arm}] {r['contig']} {r['length']} nt tier={r['tier']} "
                f"rank_len={r['rank_by_length_all_contigs']}")
            for label in dbs:
                h = r[label]
                if "subject" in h:
                    log(f"      {label}: {h['identity_percent']}% over "
                        f"{h['aligned']} (cov {h['query_coverage_union']}) "
                        f"E={h['evalue']:.3g} -> {h['subject']}")
                else:
                    log(f"      {label}: {h['state']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
