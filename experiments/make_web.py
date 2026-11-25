"""Build the web view from a run directory. SPEC.md 8.3: a thin reader, nothing more.

WHY A GENERATOR AND NOT A PAGE. A number typed into a figure is a number without
provenance, and this project has had to withdraw claims whose prose drifted from its
outputs. So no number is written into the page by hand: this script reads the committed
files in a RUNDIR, bundles exactly what the panels draw, and writes one self-contained
HTML file with that bundle inlined. Change the run and the page changes; edit the page and
the next build overwrites you.

THE THREE RULES ARE ENFORCED HERE, not left to whoever reads the page:

  - **No figure without its operating point.** Every panel takes its operating point from
    the run directory and prints it. A panel whose operating point is missing from the files
    renders the words `operating point absent` where the chips would be, rather than
    rendering clean.
  - **No pile without its interval.** Fractions are drawn with their Wilson interval, which
    `hapax.py` attached from this run's own counts. A fraction with no interval in the
    bundle is drawn as a bare tick with `no interval` beside it.
  - **No score without its calibration run.** Each field's row carries its null and its
    known-positive, and `cannot discriminate here` is rendered as a value in the table -
    the one rendering decision in this project that is a research finding rather than a
    design choice.

WHAT THE PAGE DELIBERATELY DOES NOT HAVE. No novelty score, and no clustering. An index can
say whether something is already known and cannot say whether it is new, because score
ordering tracks database proximity; and the residue is one continuum with no clusters - no
cluster count beats 1 in 5 of 5 seeds against a null with its own mean and covariance. The
joint space is therefore drawn as a density with the two reference clouds underneath, which
is what signal against noise is actually read off.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import platform
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts
import page_descriptor
import run_rows
import web_controls
from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

#: Axes offered in the joint-space panel, as (bundle key, label, how to pull it from a row).
#: Kept short on purpose: 13 fields measure about 10 independent things (participation ratio
#: 10.5) and after abundance the largest residue-against-identified standardised difference
#: is 0.40, so more axes on the page would be more breadth and no more signal.
SCATTER_AXES = [
    ("gc", "GC fraction", lambda r: r.get("gc")),
    ("entropy3", "Triplet entropy", lambda r: r.get("entropy3")),
    ("log_abundance", "log2 library occurrences",
     lambda r: _log2((_sub(r, "library_abundance").get("median_31mer_count")))),
    ("longest_orf_aa", "Longest ORF (aa)",
     lambda r: _sub(r, "frame").get("longest_orf_aa")),
    ("CpG_oe", "CpG observed/expected",
     lambda r: _sub(r, "dinuc").get("CpG_oe")),
    ("UpA_oe", "UpA observed/expected",
     lambda r: _sub(r, "dinuc").get("UpA_oe")),
    ("fold_z", "Folding energy z",
     lambda r: _sub(r, "fold").get("z")),
    ("kmer16_duplicate_frac", "Repetitive fraction",
     lambda r: _sub(r, "repetitive").get("kmer16_duplicate_frac")),
    ("best_bitscore", "Best protein bitscore", lambda r: _best(r, "bitscore")),
    ("best_pident", "Protein identity %, significant hits only",
     lambda r: _best_significant(r, "pident")),
]

#: Identity is only plotted for hits at or below this E-value. Two reasons, and the second
#: is the one a reader notices. First, an identity without an E-value is the claim this
#: project has already retracted once: a 107 nt query against a large protein database
#: finds about 34 amino acids at chance, and the chance band reads the same on residue and
#: on already-identified reads, so plotting it draws the database rather than the sample.
#: Second, percent identity over a short alignment is a RATIO OF SMALL INTEGERS, so it can
#: only land on the values n/L for the alignment length L: at L=34 the nearest values to
#: 50% are 47.1% and 52.9% and nothing can fall between. Plotted unfiltered, the chance
#: band dominates, every short alignment length contributes its own comb, and the axis
#: shows dense vertical lines at 33.3%, 50% and 100% with empty space either side. That
#: banding is an artefact of integer arithmetic on short alignments and not a property of
#: the sample, which is exactly why it is not drawn.
SIGNIFICANT_EVALUE = 1e-5


def _sub(row: dict, field: str) -> dict:
    """A per-read sub-record, or an empty one where the run did not produce it.

    A field dropped on purpose with `--report-without` records WHY it is absent rather
    than a number, and a run that did not search is not a run that found nothing. Every
    axis here used `(r.get(f) or {}).get(k)`, which reads a mapping and throws on a
    reason, so one run with a field dropped by name would have taken the whole build
    down. The reason belongs in the census the page already prints, not in a traceback.
    """
    v = row.get(field)
    return v if isinstance(v, dict) else {}


def _log2(x):
    import math
    return None if x is None or x <= 0 else round(math.log2(x), 4)


def _best(row: dict, key: str):
    """The strongest protein hit on a row, whichever index found it.

    Reported with its E-value elsewhere and never without one: an identity quoted without
    its E-value is the retracted claim this project is most careful about, because a 107 nt
    query against a large database finds about 40 bases at about 88% by chance.
    """
    best = None
    for k, v in row.items():
        if k.startswith("closest_") and isinstance(v, dict) and v.get("evalue") is not None:
            if best is None or v["evalue"] < best["evalue"]:
                best = v
    return None if best is None else best.get(key)


def _best_significant(row: dict, key: str):
    """As `_best`, but only when the hit clears the significance threshold.

    A read with no significant hit returns None and is simply absent from that axis, which
    is the honest rendering: we do not know its identity to anything, and a point at the
    chance identity would claim we did.
    """
    best = None
    for k, v in row.items():
        if (k.startswith("closest_") and isinstance(v, dict)
                and v.get("evalue") is not None and v["evalue"] <= SIGNIFICANT_EVALUE):
            if best is None or v["evalue"] < best["evalue"]:
                best = v
    return None if best is None else best.get(key)


def _one(rundir: Path, pattern: str):
    """The single matching file in a run directory, or None.

    A run holds at most one of each of these per seed and the page reads one seed, so a
    glob with a first match is the whole lookup. Returning None rather than raising is
    deliberate: the panels that use these files print what is missing and which script
    writes it, which is more use to a reader than a traceback at build time.
    """
    hits = sorted(rundir.glob(pattern))
    if not hits:
        return None
    return json.loads(hits[0].read_text())


def _trim_ranking(R: dict, seed: int, cards: int = 400, per_facet: int = 60,
                  points: int = 6000) -> dict:
    """Keep the page readable when the ranking covers the whole residue.

    Searching every unplaced pair instead of a 6,000-pair sample took the subject table
    from 243 rows to 41,829 and ranking.json from 0.2 MB to 53 MB. Embedding that whole
    file would be a page nobody can open, so the page carries three different things and
    says which is which:

    - the AGGREGATES (the group mix) are computed here over every row, so the bar a
      reader sees is the whole ranking and not the part that fitted;
    - the SCATTER is a seeded sample, because 41,829 circles is megabytes of markup and
      the shape is the point. The sample is recorded with its seed and its size;
    - the CARDS are the highest-ranked rows, plus the highest-ranked of each host group
      and each kind, so that every filter the page offers has real rows behind it rather
      than whatever survived a global cut.

    The archive keeps the whole file. This only decides what travels into the HTML.
    """
    # An empty dict is not an absent file to a JS truth test, and the panels that read
    # this key check `if (!R)`. Returning {} for a run that never had the stage therefore
    # sailed past every absent guard and threw on the first field access, taking every
    # panel after it off the page. Absent is None.
    if not R:
        return None
    rows = R.get("subjects") or []
    if not rows:
        return R
    mix, readmix = {}, {}
    for r in rows:
        mix[r["host_group"]] = mix.get(r["host_group"], 0) + 1
        readmix[r["host_group"]] = (readmix.get(r["host_group"], 0)
                                    + (r.get("residue_reads") or 0))
    R["host_mix"] = mix
    # Reads, not subjects. A share of subjects is a share of the database; a share of
    # reads is a share of the water, which is the question a reader is actually asking,
    # and it is the occurrence-weighted unit this repo reports everything else in.
    R["host_mix_reads"] = readmix
    R["subjects_total"] = len(rows)
    R["reads_total"] = sum(readmix.values())

    rng = random.Random(seed)
    withx = [r for r in rows if r.get("median_identity_percent") is not None]
    sample = withx if len(withx) <= points else rng.sample(withx, points)
    R["points"] = [{"x": r["median_identity_percent"], "y": r["residue_reads"],
                    "g": r["host_group"], "w": r["weight_term"],
                    "a": r.get("median_aligned_aa"),
                    "n": r.get("organism") or r["subject"]} for r in sample]
    R["points_note"] = {"drawn": len(sample), "of": len(withx), "seed": seed,
                        "note": ("a seeded sample of the subjects, not the top of them: "
                                 "a cut by rank would empty the lower left corner, which "
                                 "is where most of the residue is")}

    keep, seen = [], set()
    def take(lst, n):
        for r in lst[:n]:
            if r["subject"] not in seen:
                seen.add(r["subject"])
                keep.append(r)
    take(rows, cards)
    for field in ("host_group", "sequence_kind"):
        vals = {r[field] for r in rows}
        for v in sorted(vals):
            take([r for r in rows if r[field] == v], per_facet)
    keep.sort(key=lambda r: r["rank"])
    R["subjects"] = keep
    R["subjects_embedded"] = len(keep)
    return R


def _names_for(rundir: Path, subjects, table_rows) -> dict:
    """The accession-to-organism map, cut down to what this page can display."""
    full = (load(rundir / "subject_names.json") or {}).get("names") or {}
    want = {s for s in subjects if s}
    for r in table_rows:
        for v in r.values():
            if isinstance(v, dict) and v.get("subject"):
                want.add(v["subject"])
    return {k: v for k, v in full.items() if k in want}


def _trim_workup(w):
    """Keep the summary and the rows the page draws, drop the per-ORF payload it does not.

    The whole workup is 625 KB because every ORF carries its Pfam list and its runner-up
    subjects. The page draws a ranked list and a per-contig detail, so each row keeps its
    own numbers and the genes that have a role, and the rest stays in the result file.
    """
    if not w:
        return None
    rows = []
    for r in w.get("rows", []):
        rows.append({k: r[k] for k in (
            "arm", "contig", "taxon", "on_the_panel", "contig_nt", "orfs",
            "whole_contig_pident", "whole_contig_aln_aa", "polymerase_pident",
            "polymerase_divergence_class", "orf_identity_spread_points",
            "orf_identity_min", "orf_identity_max", "rank_pident", "ranked_on",
            "divergence_class", "coding_density", "subject_disagreement_flag",
            "read_support") if k in r})
        rows[-1]["genes"] = [
            {"role": o["gene_role"], "pident": o["pident"], "aln_aa": o["aln_aa"],
             "aa_len": o["aa_len"], "title": o["best_title"],
             "pfam": [d["id"] for d in o["pfam"]]}
            for o in r.get("orf_rows", []) if o.get("significant")]
        # Coordinates, for the genome maps. A gene list without them can only be drawn
        # as a bar chart of identities, which is the thing the reader could not read:
        # with them the picture is the molecule, and the identity is a label on a part
        # of it rather than a number beside a name.
        rows[-1]["orfs_drawn"] = [
            {"s": o["start"], "e": o["end"], "d": o["strand"], "aa": o["aa_len"],
             "role": o["gene_role"], "pident": o["pident"],
             "title": o["best_title"], "sig": bool(o.get("significant")),
             "pfam": [x["id"] for x in o["pfam"]]}
            for o in r.get("orf_rows", [])]
    out = {k: w[k] for k in ("summary", "operating_point", "the_number_to_rank_by",
                             "bands", "blind_spots", "seed", "git_sha",
                             "script_sha256") if k in w} | {"rows": rows}
    out["summary"] = _flaggable(w.get("rows", []), out.get("summary") or {})
    return out


def _flaggable(rows, summary):
    """Fill in the chimaera flag's own denominator, which is not the arm's contig count.

    The flag fires only where a contig carries two or more ORFs with a subject, so a
    fraction over every contig in the arm is a fraction over a pile that includes
    contigs the test could not have fired on. Both readings are kept, each beside its
    count, because the page prints the conditional one and the archived result files
    quote the unconditional one. This is computed here as well as in the stage so that a
    run produced before the stage carried the field renders the same number rather than
    silently falling back to the wider denominator, which reads lower and in the
    direction of "the flag separates the two arms".
    """
    summary = {k: dict(v) if isinstance(v, dict) else v for k, v in summary.items()}
    for arm, s in summary.items():
        if not isinstance(s, dict) or "contigs_that_could_be_flagged" in s:
            continue
        rs = [r for r in rows if r.get("arm") == arm]
        flagged = sum(1 for r in rs if r.get("subject_disagreement_flag"))
        flaggable = sum(1 for r in rs
                        if (r.get("orfs_with_a_subject") or 0) >= 2
                        and (r.get("orfs") or 0) > 1)
        if flagged > flaggable:
            raise ValueError(
                f"arm {arm}: {flagged} contigs carry the subject-disagreement flag out "
                f"of {flaggable} that carry two or more named ORFs, which is the only "
                "way it can fire")
        s["contigs_that_could_be_flagged"] = flaggable
        s["fraction_flagged_of_the_contigs_that_could_be"] = (
            round(flagged / flaggable, 4) if flaggable else None)
    return summary


def load(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


ROOT = Path(__file__).resolve().parent.parent

#: The point on the spectrum below which a reported identity cannot be read. Both are
#: measured, not chosen. Below half the contig covered, per-read identity flatters towards
#: "familiar" by +17.2 points against +7.1 above 90% coverage (`results/holdout_run/`), and
#: below about 100 aligned residues identity is quantised in steps a reader will
#: over-interpret. A contig that fails either does not get a position on the axis; it goes
#: into a lane of its own, because a point drawn at a number that cannot be read is worse
#: than an admission.
SPECTRUM_MIN_COVERAGE = 0.5
SPECTRUM_MIN_ALN_AA = 100


def spectrum(inputs: dict) -> dict | None:
    """The distance spectrum: one continuous axis, no bucket boundary anywhere on it.

    WHY A SPECTRUM AND NOT TIERS. Identity against a reference orders contigs by database
    proximity, so any line drawn across it is a statement about the database rather than
    about novelty. The whole pile is therefore one ranked continuum, and the only
    categorical thing the panel draws is whether the measurement can be read at all.

    WHY THE POSITION COMES FROM ONE GENE AND NOT THE CONTIG. A whole-contig identity is
    whichever gene won on bitscore, and a divergent contig's own genes disagree by a median
    12.8 points against the near-reference arm's 2.45 (`results/divergent_known/`). So the
    axis takes the polymerase where one is called, the whole-contig hit otherwise, says
    which on every row, and draws the spread between a contig's genes as a whisker. A dot
    without its whisker is the error this panel exists to stop.

    WHY COVERAGE IS THE COLOUR. Identity and coverage are not two optional readouts; an
    identity without the fraction covered cannot tell a diverged full-length protein from a
    conserved fragment of one, and those are different findings. Putting coverage on the
    colour ramp means a thin alignment cannot be read as a confident number by accident.

    WHY THE BIG-INDEX RESULT IS A RING AND NOT AN ARROW. core_nt is a nucleotide index and
    this axis is protein identity. Drawing a contig as moving from one to the other would be
    comparing two instruments, which is the move this project refuses. So the 166 contigs
    that also went to core_nt carry a mark saying what the larger index did to them, in its
    own units, and the slot for a true displacement arrow waits on the protein pass against
    nr, which is the same instrument at a larger size.
    """
    # Every file this panel reads is one the page descriptor declared, by path, for THIS
    # sample. It used to read three repo-root result files whatever run directory it was
    # given, so the holdout page drew the Berlin contigs as its own spectrum and nothing
    # on the page said so.
    if inputs["spectrum"] is None:
        return None
    dk = load(inputs["spectrum"])
    big = {}
    cn = inputs["spectrum_big_index"]
    if cn is not None:
        for line in cn.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if not r.get("hit") or r["evalue"] >= 1e-5:
                big[r["read"]] = {"state": "no significant hit"}
                continue
            state = ("collapsed to a strain-level identification"
                     if r["pident"] >= 90 and r["aln_frac"] >= 0.5
                     else "still distant" if r["pident"] < 80 or r["aln_frac"] < 0.5
                     else "closer, still not an identification")
            big[r["read"]] = {"state": state, "pident": round(r["pident"], 1),
                              "frac": round(r["aln_frac"], 3)}

    on_axis, lane = [], []
    for r in dk["rows"]:
        glist = [g for g in r.get("genes") or [] if g.get("pident") is not None]
        genes = [g["pident"] for g in glist]
        # The coverage that belongs beside a gene-level identity is the fraction of THAT
        # gene aligned, not the fraction of the contig. A 4 kb contig with one scorable
        # gene covers a tenth of itself and the gene can still be aligned end to end; the
        # measured familiarity bias is a property of how much of the scored sequence the
        # alignment spans, so it is read on the same unit the identity came from.
        rank = next((g for g in glist if g.get("role") == "polymerase"), None)
        x = r.get("polymerase_pident")
        src = "polymerase"
        if x is None or rank is None:
            x, src, rank = r["pident"], "whole contig", None
        if rank is not None:
            cov = round(rank["aln_aa"] / max(1, rank.get("aa_len") or 1), 3)
            aln = rank["aln_aa"]
        else:
            cov, aln = r["subject_coverage"], r["aln_aa"]
        why = ("under half the scored gene covered" if cov < SPECTRUM_MIN_COVERAGE
               else "alignment under 100 aa" if aln < SPECTRUM_MIN_ALN_AA else None)
        row = {
            "name": r["contig"], "arm": r["arm"], "x": x, "src": src,
            # The whole-contig best hit travels on every row whether or not the row is
            # ranked on it, because SPEC 13.3 makes the choice of ranking a control: the
            # reader can put the figure back on the whole-contig identity and see what
            # that does. Both numbers are measured and neither is recomputed in the page.
            "x_whole": r["pident"], "cov_whole": r["subject_coverage"],
            "aln_whole": r["aln_aa"],
            "lo": min(genes) if len(genes) > 1 else None,
            "hi": max(genes) if len(genes) > 1 else None,
            "genes": len(genes), "cov": cov, "aln": aln,
            "nt": r["contig_nt"], "evalue": r["evalue"],
            "depth": (r.get("read_support") or {}).get("mapped_median_depth"),
            "subject": (r.get("subject_title") or r.get("subject") or "")[:90],
            "big": big.get(r["contig"]),
        }
        (lane if why else on_axis).append({**row, "why": why} if why else row)

    ex = dk.get("external_positive_control") or {}
    return {
        "on_axis": on_axis,
        "cannot_read": lane,
        "floors": {"min_coverage": SPECTRUM_MIN_COVERAGE,
                   "min_aln_aa": SPECTRUM_MIN_ALN_AA},
        "census": {
            "searched": (dk.get("discrimination") or {}).get("existence", {})
                        .get("contigs_searched"),
            "with_a_significant_hit": (dk.get("discrimination") or {})
                                      .get("existence", {}).get("real_with_a_significant_hit"),
            "shuffled_with_a_significant_hit": (dk.get("discrimination") or {})
                                               .get("existence", {})
                                               .get("shuffled_with_a_significant_hit"),
            "drawn": len(on_axis) + len(lane),
            "note": ("the panel draws the gene-level subset, which is every contig short "
                     "enough for per-ORF scoring; the census above it is the whole pile"),
        },
        "holdout_scope": ex or None,
        # Left empty on purpose rather than omitted, so the page says the axis is coming
        # instead of the reader discovering later that it was never there. Filled when the
        # structure pass reports, and only if it reaches past the protein search.
        "structure_axis": {"state": "not measured yet",
                           "what_would_fill_it": ("a per-gene structural distance - "
                                                  "ProstT5 plus Foldseek, or "
                                                  "profile-to-profile - taken only if it "
                                                  "reaches past DIAMOND on this material")},
    }


#: Every file outside the run directory that a panel may read, and what it is. A page
#: descriptor must state each one, either as a path or as `{"not_run": reason}`: a key it
#: leaves out is a refusal, because "nobody said" and "it does not apply" must not render
#: the same. Nothing here is ever found by globbing, because a sample that acquired an
#: input by filename would be reporting on somebody else's water.
PAGE_INPUTS = {
    "spectrum": "experiments/divergent_known.py output: the contigs the distance "
                "spectrum draws",
    "spectrum_big_index": "core_nt block marking what the larger index did to those "
                          "contigs",
}

#: The six fixed rows of a Methods entry, in order. Two spellings of the last two are in
#: use and both are accepted; anything else is a refusal.
METHOD_ROWS = [("Input",), ("Procedure",), ("Threshold",), ("Denominator",),
               ("Control", "Control arm"), ("Bias", "Bias direction")]
METHODS_FILE = ROOT / "results/web/methods.json"
NOT_MEASURED = "Not measured on this library."
#: Appended to a row that states its design but holds no number for this sample. A row
#: that printed the design alone would read as a measurement; NOT_MEASURED alone, which is
#: what every such row used to read, threw the design away with the number. On a library
#: nobody has written about, 28 of the 90 rows were blank, and nine of them were blank for
#: the wrong reason: a shuffled control arm and the direction of a bias are properties of
#: the procedure and hold whatever water is in the tube.
NO_VALUE = "No number for this row is recorded for this library."


def load_page(path: Path, rundir: Path | None) -> dict:
    """Read and check a page descriptor: the only per-sample input the generator takes.

    WHY A FILE AND NOT FLAGS. Two pages built from one generator with different flags are
    two generators: whatever the holdout page needed was a flag Berlin did not pass, and
    whatever Berlin said was typed into the renderer and printed on the holdout page too.
    So everything a page says about its sample that is not in the run directory - the
    files outside the run it reads, and the run itself - lives in one committed file per
    sample, and every sample goes through the same command with it.

    WHAT IS NO LONGER IN IT. The title and the one-line description of the water were
    typed here, which meant a sample nobody had touched could not have a page until
    somebody wrote a sentence for it, and the sentence was the only unattributable English
    on a page of numbers that all carry provenance. They are now assembled by
    `page_descriptor.derive()` from the library's own facts. A descriptor may still state
    either, and that is honoured - but as an override, recorded as one in the bundle
    beside what would have been said, because a silent substitution in either direction
    is the bug this repo keeps paying for.
    """
    d = json.loads(path.read_text())
    problems = []
    for k in ("sample", "rundir", "inputs"):
        if k not in d:
            problems.append(f"missing key {k!r}")
    smp = d.get("sample") or {}
    if not smp.get("accession"):
        problems.append("missing sample.accession")
    else:
        try:
            # `sample_facts`, where a descriptor carries it, is the water's own facts
            # supplied by the caller instead of read from the registry. That is what lets
            # a run on an input nobody has deposited from have a page at all: the
            # registry raises for a library it does not know, correctly, and a production
            # sample is by definition one it does not know. Supplied facts are recorded
            # as an override below, like a stated title.
            derived = page_descriptor.derive(smp["accession"], d.get("sample_facts"))
        except library_facts.UnknownLibrary:
            # A library the registry does not know is the ordinary case for this tool:
            # production water has no entry, and the registry raising is correct, since
            # the alternative is another library's facts. So the page is built stating
            # that it knows nothing about the water, and the statement is RECORDED as an
            # override rather than made quietly - the page then says it has no facts,
            # instead of carrying a sentence nobody can attribute.
            derived = page_descriptor.derive(smp["accession"], {})
            d["unknown_library"] = True
        # Derived unless overridden, and an override says so rather than looking derived.
        d["derived"] = derived
        d["overrides"] = sorted(
            [k for k in ("title",) if d.get(k) and d[k] != derived.get(k)] +
            [f"sample.{k}" for k in ("lede",)
             if smp.get(k) and smp[k] != derived.get(k)] +
            (["sample_facts"] if d.get("sample_facts") is not None else []) +
            (["unknown_library"] if d.get("unknown_library") else []))
        d.setdefault("title", derived.get("title"))
        if not smp.get("lede"):
            smp["lede"] = derived.get("lede")
            d["sample"] = smp
    inputs = d.get("inputs") or {}
    for k in sorted(set(PAGE_INPUTS) - set(inputs)):
        problems.append(f"inputs.{k} is not declared: give a path or "
                        f'{{"not_run": "why"}} ({PAGE_INPUTS[k]})')
    for k in sorted(set(inputs) - set(PAGE_INPUTS)):
        problems.append(f"inputs.{k} is not an input any panel reads")
    resolved = {}
    for k, v in inputs.items():
        if k not in PAGE_INPUTS:
            continue
        if not isinstance(v, dict) or (("path" in v) == ("not_run" in v)):
            problems.append(f"inputs.{k} must be exactly one of path or not_run")
            continue
        if "path" in v:
            f = ROOT / v["path"]
            if not f.exists():
                # Declared but absent is not "not run": somebody said it exists.
                problems.append(f"inputs.{k}: {v['path']} does not exist")
            resolved[k] = f
        else:
            resolved[k] = None
    if resolved.get("spectrum") is None and resolved.get("spectrum_big_index"):
        problems.append("spectrum inputs are declared without the spectrum itself")
    run = rundir if rundir is not None else ROOT / d.get("rundir", "")
    if not (run / "manifest.json").exists():
        problems.append(f"{run}/manifest.json is absent; this is not a hapax run "
                        f"directory (pass --rundir if the run lives elsewhere)")
    else:
        # The descriptor and the run must be the same sample. Without this a descriptor
        # could dress any run directory in another library's words.
        inp = (json.loads((run / "manifest.json").read_text()).get("inputs") or {})
        if smp.get("accession") and smp["accession"] not in str(inp.get("r1", "")):
            problems.append(f"sample.accession {smp['accession']} is not in the run's "
                            f"own input path {inp.get('r1')!r}")
    if problems:
        raise SystemExit(f"{path}: the page descriptor is refused:\n  - " +
                         "\n  - ".join(problems))
    d["_inputs"] = resolved
    d["_rundir"] = run
    return d


def _unmeasured(entry: dict, label: str) -> str:
    """What a per-library row says on a page whose library nobody has written about.

    The wording separates two things a blank row used to merge: a control that was never
    run, and a control that ran under a number nobody has written down. This row is the
    second, so it says the number is not recorded rather than that nothing was measured.

    The row's own number is not available and is never invented. What can still be printed
    is the design: which control arm this figure carries, and which way its bias points.
    Those are properties of the procedure, so they hold for whatever water is in the tube,
    and dropping them left a reader of a new library's page unable to see that a figure had
    a control at all - on SRR34217566, 28 of 90 rows read "Not measured on this library."
    with nine of them stating no number in the first place.
    """
    design = (entry.get("design") or {}).get(label)
    return f"{design} {NO_VALUE}" if design else NOT_MEASURED


def resolve_methods(accession: str) -> list[dict]:
    """The Methods entries as this sample's page prints them.

    Every page carries every entry, in one order, so a citation number means the same
    thing on every page. What differs is which rows are filled. A row that states a
    procedure or a threshold holds for every run and is written once. A row that states a
    number measured on one library is stored under that library's accession and printed
    only on that library's page; on any other page the row reads "Not measured on this
    library." That is the fix for the holdout page printing "SRR23973193, Berlin
    influent" as its own input row.

    An entry `about` one fixed sample (the withheld-genome fixture) is a claim about that
    sample wherever it is printed, so it must name the sample in its own Input row and can
    then never be read as a statement about the page it sits on.
    """
    doc = json.loads(METHODS_FILE.read_text())
    out, seen = [], set()
    for e in doc["entries"]:
        if e["id"] in seen:
            raise SystemExit(f"{METHODS_FILE}: entry {e['id']!r} appears twice")
        seen.add(e["id"])
        rows = e["rows"]
        labels = [r[0] for r in rows]
        if len(rows) != len(METHOD_ROWS) or any(
                lab not in ok for lab, ok in zip(labels, METHOD_ROWS)):
            raise SystemExit(f"{METHODS_FILE}: entry {e['id']!r} rows are {labels}; a "
                             f"Methods entry has exactly the six fixed rows in order")
        meas = (e.get("measured") or {}).get(accession) or {}
        if e["about"] != "tool":
            if e["about"] not in (rows[0][1] or ""):
                raise SystemExit(f"{METHODS_FILE}: entry {e['id']!r} is about "
                                 f"{e['about']} but its Input row does not say so")
            filled = rows
        else:
            blank = [lab for lab, txt in rows if txt is None]
            for lab in meas:
                if lab not in labels:
                    raise SystemExit(f"{METHODS_FILE}: {e['id']!r} measured row {lab!r} "
                                     f"for {accession} is not one of its rows")
                if lab not in blank:
                    raise SystemExit(
                        f"{METHODS_FILE}: {e['id']!r} measures row {lab!r} for "
                        f"{accession} where the row already holds text every page "
                        f"prints; one library's number must not shadow a statement "
                        f"made about the tool")
            for lab in (e.get("design") or {}):
                if lab not in blank:
                    raise SystemExit(
                        f"{METHODS_FILE}: {e['id']!r} states a design for row {lab!r}, "
                        f"which already holds text that every page prints; a design is "
                        f"what a row says when the number is this library's and the "
                        f"number is absent")
            filled = [[lab, txt if txt is not None
                       else (meas.get(lab) or _unmeasured(e, lab))]
                      for lab, txt in rows]
        out.append({"id": e["id"], "t": e["title"], "rows": filled})
    return out


#: One row per scored pair, per arm, as [e1, b1, e2, b2, depth, occupancy]: the strongest
#: E-value and bitscore of each mate, the pair's median 31-mer count in the library, and
#: how many other samples carry it. Six numbers a mate, so both arms of every pair the run
#: scored travel into the page rather than a sample of one of them.
SEARCH_FIELDS = ["evalue_mate1", "bitscore_mate1", "evalue_mate2", "bitscore_mate2",
                 "library_depth", "other_samples"]


def not_searched(piles: dict[str, list[dict]]) -> dict:
    """The fields a run dropped by name, with the run's own reason, from the rows.

    `hapax run --report-without X` writes the reason into the field where the value would
    be, as a plain string, so the rows themselves say which fields were not searched and
    why. That census is read once here and carried, rather than each reader inventing its
    own wording for an absence: *searched and found nothing* and *not searched* must not
    render the same, and this is the sentence that keeps them apart.
    """
    out = {}
    for rows in piles.values():
        for r in rows:
            for k, v in r.items():
                if isinstance(v, str) and k not in ("name", "set", "adapter") \
                        and v.lower().startswith("not searched"):
                    out.setdefault(k, v)
        if out:
            break
    return out


def search_table(piles: dict[str, list[dict]]) -> dict:
    """The per-read search numbers the significance controls filter, both arms together.

    WHY BOTH ARMS, ALWAYS, IN ONE PASS. A threshold quoted without its own control arm is
    the exact shape behind the readings this project has withdrawn: the hit rate of the
    residue at some floor is uninterpretable until the same floor is applied to material
    the cascade already identified. So the residue rows and the identified rows are
    carried together and filtered by the same predicate in the same pass, and a panel
    draws both or neither. The shuffled arm is NOT here, because this run has no shuffled
    search at read level - the composition fields have their shuffled null, the search
    does not - and the page says that rather than drawing two arms as though they were
    three.

    WHY PER MATE AND NOT PER PAIR. The two mates are not interchangeable: mate 2 carries a
    non-template prefix and is 14 nt shorter, which is 62 scorable 31-mers against 76, and
    it hits at a different rate. A per-pair number averages two biases pointing opposite
    ways, so the unit is a control and both are available.
    """
    out = {}
    for pile, rows in piles.items():
        acc = []
        for r in rows:
            per = {1: [None, None], 2: [None, None]}
            for k, v in r.items():
                if not (k.startswith("closest_") and isinstance(v, dict)):
                    continue
                m = 1 if str(v.get("mate", "1")) == "1" else 2
                e, b = v.get("evalue"), v.get("bitscore")
                if e is not None and (per[m][0] is None or e < per[m][0]):
                    per[m][0] = e
                if b is not None and (per[m][1] is None or b > per[m][1]):
                    per[m][1] = b
            acc.append([per[1][0], per[1][1], per[2][0], per[2][1],
                        _sub(r, "library_abundance").get("median_31mer_count"),
                        _sub(r, "cross_sample").get("distinct_samples")])
        out[pile] = acc
    return out


def _shown(path: Path) -> str:
    """The path as the page should print it: repo-relative where it is, absolute otherwise.

    A run directory is not required to sit inside the checkout - the production runs land
    under /opt/hapax/out and are fetched afterwards - and `relative_to` raises on one that
    does not, so this used to refuse any page whose inputs the run itself produced. The
    point of printing it at all is that a reader can find the file, which an absolute path
    serves and an exception does not.
    """
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def bundle(rundir: Path, scatter_per_pile: int, table_rows: int, seed: int,
           page: dict) -> dict:
    rng = random.Random(seed)
    manifest = load(rundir / "manifest.json")
    if manifest is None:
        raise SystemExit(f"{rundir}/manifest.json is absent; this is not a hapax run "
                         f"directory")
    piles: dict[str, list[dict]] = {"residue": [], "known": []}
    for r in run_rows.jsonl(rundir, "reads.jsonl") or []:
        piles.setdefault("residue" if r.get("set", "residue") == "residue"
                         else "known", []).append(r)

    scatter = {}
    for pile, rows in piles.items():
        take = rows if len(rows) <= scatter_per_pile else rng.sample(rows,
                                                                    scatter_per_pile)
        scatter[pile] = {"n_drawn": len(take), "n_total": len(rows),
                         "points": [[f(r) for _k, _l, f in SCATTER_AXES] for r in take]}

    # The table is a sample too, and says so. Drawing 20,000 rows into the page would
    # push it past the size limit and tell the reader nothing the sample does not.
    res = piles.get("residue", [])
    table = res if len(res) <= table_rows else rng.sample(res, table_rows)
    layers = {}
    for d in run_rows.jsonl(rundir, "layers.jsonl") or []:
        layers[d["name"]] = d["layer"]

    strains_subjects = [r.get("subject") for r in
                        ((load(rundir / "strains.json") or {}).get("subjects") or [])]

    search = search_table(piles)
    out = {
        "run": {k: manifest.get(k) for k in
                ("tool", "spec", "seed", "git_sha", "tree_dirty", "script_sha256",
                 "inputs", "complete", "incomplete_because", "wall_seconds",
                 "field_calibration", "known_positive_rows", "known_positive_note",
                 "interface_rules", "stages")},
        "parameters": manifest.get("parameters", {}),
        "spectrum": spectrum(page["_inputs"]),
        # The page's biology, as data. Every plain-English claim it makes about what
        # a gene does or what an organism is comes from here, so the renderer has no
        # vocabulary of its own and a wrong word is wrong in one place.
        "vocabulary": load(ROOT / "experiments/virus_vocabulary.json"),
        "architecture": _one(rundir, "architecture_s*.json"),
        # `derived` is what the library's own facts say this page should call itself and
        # the water; `overrides` names any of it the descriptor stated by hand, so a
        # hand-written sentence is visible as one instead of passing for a derived fact.
        "page": {"title": page["title"], "sample": page["sample"],
                 "derived": page.get("derived"), "overrides": page.get("overrides") or [],
                 "inputs": {k: ({"path": _shown(v)} if v is not None else
                                {"not_run": page["inputs"][k]["not_run"]})
                            for k, v in page["_inputs"].items()}},
        "methods": resolve_methods(page["sample"]["accession"]),
        "tier_table": load(rundir / "tier_table.json"),
        "layers": load(rundir / "layers.json"),
        "fields": load(rundir / "fields.json"),
        "linkage": load(rundir / "linkage.json"),
        "strains": load(rundir / "strains.json"),
        "census": load(rundir / "census.json"),
        "calibration_ladder": load(rundir / "calibration_ladder.json"),
        "blind_spots": (load(rundir / "blind_spots.json") or {}).get("blind_spots"),
        # The two public-health files. Absent is a legitimate state and the page says so
        # rather than drawing an empty panel: a run made before these stages existed is
        # still a valid run of the stages it had.
        "profile": _one(rundir, "public_health_profile_s*.json"),
        "strain_divergence": _one(rundir, "known_strain_divergence_s*.json"),
        # The gene-level workup of the divergent end. Absent is a legitimate state: the
        # divergence panel then ranks on the whole-contig hit and says so, which is what
        # it did before this stage existed.
        "contig_workup": _trim_workup(_one(rundir,
                                           "divergent_contig_workup_s*.json")),
        # Accession to organism, from the same local FASTA the DIAMOND database was built
        # from. Absent is a legitimate state: the page then shows the accession alone
        # rather than inventing a label.
        # Only the accessions the page can actually show. The full map is 14,743 records
        # and two megabytes of page nobody reads; the strains panel and the sampled read
        # table between them name a few hundred.
        "names": _names_for(rundir, strains_subjects, table),
        "significant_evalue": SIGNIFICANT_EVALUE,
        "ranking": _trim_ranking(load(rundir / "ranking.json"), seed),
        "axes": [{"key": k, "label": l} for k, l, _f in SCATTER_AXES],
        "scatter": scatter,
        "table": {"rows": [{"name": r.get("name"),
                            "layer": layers.get(r.get("name")),
                            # Carried so the table's own filters are the controls of
                            # SPEC 13.3 rather than a second set of numbers: the depth
                            # stratum and the occupancy cut read these, and the control
                            # arm reads the same fields on the identified rows.
                            "occupancy": _sub(r, "cross_sample")
                                         .get("distinct_samples"),
                            "depth": _sub(r, "library_abundance")
                                     .get("median_31mer_count"),
                            "values": [f(r) for _k, _l, f in SCATTER_AXES],
                            "closest": {k: v for k, v in r.items()
                                        if k.startswith("closest_")
                                        and isinstance(v, dict)},
                            # Same rule as _sub: a field dropped by name carries its
                            # reason, and the page reads a hit only where there is a
                            # list of them. The reason is a census entry, not a row.
                            "motif": (r["motif"][0] if isinstance(r.get("motif"), list)
                                      and r["motif"] else None),
                            "chimaera": (r.get("chimaera")
                                         if isinstance(r.get("chimaera"), dict)
                                         else None)}
                           for r in table],
                  "n_drawn": len(table), "n_total": len(res), "sample_seed": seed},
        "not_searched": not_searched(piles),
        "search": {"fields": SEARCH_FIELDS,
                   "piles": {k: v for k, v in search.items()},
                   "n": {k: len(v) for k, v in search.items()}},
        "built": {"script": "experiments/make_web.py", **provenance(__file__),
                  "rundir": (str(rundir.resolve().relative_to(ROOT))
                             if rundir.resolve().is_relative_to(ROOT) else str(rundir))},
    }
    # The controls come last because they are resolved against the bundle the panels will
    # draw: a slider's range is derived from the rows it filters, so it cannot run past
    # them. web_controls refuses the build where the declaration and the data disagree.
    out["controls"] = web_controls.resolve(
        rundir, out, {**search, "all": [r for v in search.values() for r in v]})
    return out


PAGE = """<title>Hapax run viewer</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: one column of instrument panels, each a readout with its operating point. */
:root {
  --bg: #f5f6f4; --panel: #ffffff; --edge: #d8dcd6; --fg: #15201c; --dim: #5c6b64;
  --identified: #41708a; --tier2: #8a7a41; --residue: #b4603a; --null: #9aa39e;
  --ok: #3f7d55; --warn: #a2622a;
  --g-phage: #78889a; --g-plant: #4f8a46; --g-fungal: #8168a6; --g-invert: #2f877f;
  --g-vert: #b07d22; --g-human: #c0392b; --g-unknown: #9aa39e;
  /* Chart slots, validated with the dataviz skill's checker in both modes under
     --pairs all: all checks pass, with one light-mode contrast WARN on --v3 that
     obliges the visible direct labels every ranked row and donut legend carries. */
  --v1: #2a78d6; --v2: #eb6834; --v3: #1baf7a;
  --cv1: #d6e3ea; --cv2: #a9c6d5; --cv3: #7ba6bf; --cv4: #4f85a5;
  --cv5: #2a5b76;
  --v2-soft: #f6c3ac; --v-mute: #9aa39e; --v-faint: #e4e8e3;
  --sans: "IBM Plex Sans", system-ui, -apple-system, Segoe UI, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
  --g-phage: #93a3b4; --g-plant: #74b069; --g-fungal: #a68cc8; --g-invert: #4fa89f;
  --g-vert: #d2a04a; --g-human: #e0695c; --g-unknown: #68756f;
  --v1: #3987e5; --v2: #d95926; --v3: #199e70;
  --cv1: #22333d; --cv2: #33566a; --cv3: #497e99; --cv4: #6ea2bf;
  --cv5: #9cc6dc;
  --v2-soft: #7d3f26; --v-mute: #68756f; --v-faint: #222b27;
} }
:root[data-theme="dark"] {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
  --g-phage: #93a3b4; --g-plant: #74b069; --g-fungal: #a68cc8; --g-invert: #4fa89f;
  --g-vert: #d2a04a; --g-human: #e0695c; --g-unknown: #68756f;
  --v1: #3987e5; --v2: #d95926; --v3: #199e70;
  --cv1: #22333d; --cv2: #33566a; --cv3: #497e99; --cv4: #6ea2bf;
  --cv5: #9cc6dc;
  --v2-soft: #7d3f26; --v-mute: #68756f; --v-faint: #222b27;
}
body { background: var(--bg); color: var(--fg); font-family: var(--sans);
  font-size: 15px; line-height: 1.5; }
.wrap { max-width: 1040px; margin: 0 auto;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 2.5rem 4rem; display: flex; flex-direction: column; gap: 1.75rem; }
h1 { font-size: 1.9rem; font-weight: 600; letter-spacing: -0.02em; margin: 0;
  text-wrap: balance; }
h2 { font-size: 1.05rem; font-weight: 600; margin: 0; letter-spacing: -0.01em; }
p { margin: 0; max-width: 68ch; }
.lede { color: var(--dim); max-width: 68ch; }
.panel { background: var(--panel); border: 1px solid var(--edge); border-radius: 4px;
  padding: 1.1rem 1.15rem; display: flex; flex-direction: column; gap: 0.8rem;
  min-width: 0; }
.op { display: flex; flex-wrap: wrap; gap: 0.35rem; }
.chip { font-family: var(--mono); font-size: 11.5px; color: var(--dim);
  border: 1px solid var(--edge); border-radius: 2px; padding: 0.1rem 0.4rem; }
.chip b { color: var(--fg); font-weight: 500; }
.absent { color: var(--warn); font-family: var(--mono); font-size: 12px; }
table { border-collapse: collapse; width: 100%; font-size: 13px;
  font-variant-numeric: tabular-nums; }
th, td { text-align: right; padding: 0.3rem 0.5rem; border-bottom: 1px solid var(--edge);
  white-space: nowrap; }
th:first-child, td:first-child { text-align: left; }
th { font-weight: 500; color: var(--dim); font-size: 11.5px; text-transform: uppercase;
  letter-spacing: 0.04em; }
.scroll { overflow-x: auto; min-width: 0; }
.mono { font-family: var(--mono); }
.dim { color: var(--dim); }
.bar { display: flex; height: 34px; border-radius: 3px; overflow: hidden;
  border: 1px solid var(--edge); }
.bar > span { display: block; }
.legend { display: flex; flex-wrap: wrap; gap: 0.9rem; font-size: 12.5px;
  color: var(--dim); }
.legend i { display: inline-block; width: 10px; height: 10px; border-radius: 2px;
  margin-right: 0.35rem; }
.controls { display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;
  font-size: 13px; }
select, input[type=range] { font: inherit; font-size: 13px; }
svg { display: block; max-width: 100%; height: auto; }
.grid2 { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 0.9rem; }
.note { font-size: 12.5px; color: var(--dim); max-width: 68ch; }
.rule { font-size: 12px; color: var(--dim); font-family: var(--mono); }
details > summary { cursor: pointer; font-size: 13px; color: var(--dim); }
.bs { border-left: 2px solid var(--residue); padding-left: 0.7rem; }
.bs h3 { font-size: 13.5px; margin: 0 0 0.2rem; font-weight: 600; }
.bs p { font-size: 12.5px; color: var(--dim); }
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(270px, 1fr));
  gap: 0.75rem; margin-top: 0.9rem; }
.card { border: 1px solid var(--edge); border-radius: 6px; padding: 0.8rem 0.85rem;
  min-width: 0; }
.ctitle { font-size: 14.5px; font-weight: 600; line-height: 1.25; overflow-wrap: anywhere; }
.csub { font-size: 12.5px; color: var(--dim); margin-top: 0.15rem;
  overflow-wrap: anywhere; }
.cbar { height: 5px; border-radius: 3px; background: var(--edge); margin: 0.6rem 0 0.5rem; }
.cbar > span { display: block; height: 100%; border-radius: 3px;
  background: var(--residue); }
.crow { font-size: 12.5px; margin-top: 0.3rem; overflow-wrap: anywhere; }
.idbar { height: 4px; border-radius: 2px; background: var(--edge); margin-top: 0.3rem; }
.idbar > span { display: block; height: 100%; border-radius: 2px;
  background: var(--identified); }
.tag { display: inline-block; margin-top: 0.55rem; font-size: 11.5px;
  padding: 0.1rem 0.45rem; border-radius: 10px; border: 1px solid var(--edge);
  color: var(--dim); }
.tag-only { border-color: var(--residue); color: var(--residue); }
.callout { border: 1px solid var(--warn); border-left-width: 3px; border-radius: 4px;
  padding: 0.7rem 0.8rem; font-size: 13px; line-height: 1.5; max-width: 76ch; }
.rankno { font-family: var(--mono); font-size: 11.5px; color: var(--dim); }
text.ax { font-size: 11px; fill: var(--dim); font-family: var(--sans); }
.axrow { display: grid; grid-template-columns: 4.2rem 1fr auto; gap: 0.5rem;
  align-items: center; font-size: 11.5px; color: var(--dim); margin-top: 0.3rem; }
.axlab { font-family: var(--mono); }
.axbar { background: var(--edge); border-radius: 2px; height: 6px; overflow: hidden; }
.axbar > i { display: block; height: 100%; border-radius: 2px; }
.axval { font-family: var(--mono); font-variant-numeric: tabular-nums; color: var(--fg); }
.hero { display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 0.9rem; }
.hero > div { border: 1px solid var(--edge); border-radius: 6px; padding: 0.85rem 0.9rem; }
.hero .bignum { font-size: 34px; }
.stack { display: flex; flex-direction: column; gap: 0.9rem; }
.barlab { font-size: 12px; color: var(--dim); font-family: var(--mono);
  margin-bottom: 0.25rem; }
/* Dynamism, kept to two behaviours so it reads as one page and not as effects:
   a bar grows from nothing the first time it is scrolled into view, and hovering any
   split desaturates the slices you are not pointing at. Both are suppressed under
   prefers-reduced-motion, and the grow is a no-op if the observer never fires, so the
   page is correct with JavaScript half-done. */
.bar > span, .axbar > i { transition: width 520ms cubic-bezier(.22,.61,.36,1); }
.bar:hover > span, .bites:hover .bbar > i, .flowwrap:hover .flowseg {
  filter: saturate(0.2); transition: filter 140ms linear; }
.bar > span:hover, .bite:hover .bbar > i, .flowseg:hover { filter: none; }
.bar > span, .flowseg, .bbar > i { cursor: default; }
@media (prefers-reduced-motion: reduce) {
  .bar > span, .axbar > i, .bbar > i { transition: none !important; }
}
.clicky { cursor: pointer; border-bottom: 1px dotted var(--edge); }
.clicky:hover { color: var(--fg); border-bottom-color: var(--dim); }
@keyframes fanin { from { opacity: 0 } }
.fanrib { animation: fanin 0.7s ease both; }
@media (prefers-reduced-motion: reduce) { .fanrib { animation: none } }
.panelgrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr));
  gap: 0.9rem 1.4rem; margin: 0.7rem 0; }
.pghead { font-size: 0.74rem; text-transform: uppercase; letter-spacing: 0.06em;
  color: var(--dim); margin: 0.2rem 0 0.35rem; }
.prow { padding: 0.28rem 0; border-bottom: 1px solid var(--edge); font-size: 0.84rem; }
.prow .pname { color: var(--dim); }
.prow.found .pname b { color: var(--fg); font-weight: 600; }
.prow .pval { font-size: 0.78rem; margin-top: 0.05rem; }
.prow:not(.found) .pval { font-size: 0.76rem; }
.chip { font-size: 0.66rem; text-transform: uppercase; letter-spacing: 0.05em;
  border: 1px solid var(--edge); border-radius: 3px; padding: 0.02rem 0.26rem;
  color: var(--dim); white-space: nowrap; }
.bar.big { height: 30px; }
.twocol { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  gap: 1.4rem; margin-top: 0.9rem; }
.trow { display: grid; grid-template-columns: 1fr 70px 62px 58px; gap: 0.5rem;
  align-items: center; font-size: 0.82rem; padding: 0.11rem 0; }
.tname { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.tbar { height: 9px; background: var(--edge); border-radius: 2px; overflow: hidden; }
.tbar span { display: block; height: 100%; }
.tnum, .tpct { text-align: right; font-family: var(--mono); font-size: 0.75rem; }
@media (max-width: 680px) {
  /* A taxon name clipped to "Tobamovirus..." identifies nothing, and every row in a
     genus reads the same. On a phone the name takes its own line and wraps. */
  .trow { grid-template-columns: 1fr auto; }
  .tname { white-space: normal; overflow: visible; text-overflow: clip;
    grid-column: 1 / -1; overflow-wrap: anywhere; }
}
.tpct { color: var(--dim); }
.drows { display: flex; flex-direction: column; gap: 0.25rem; margin: 0.6rem 0; }
.drow { display: grid; grid-template-columns: 1fr 160px 58px 1.1fr; gap: 0.6rem;
  align-items: center; font-size: 0.84rem; padding: 0.2rem 0;
  border-bottom: 1px solid var(--edge); }
.dbar { height: 10px; background: var(--edge); border-radius: 2px; overflow: hidden; }
.dbar span { display: block; height: 100%; }
.dnum { text-align: right; }
.dmeta { font-size: 0.76rem; }
.cards3 { display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 0.8rem; margin: 0.7rem 0; }
.cmp { margin: 1.1rem 0 0.6rem; }
.crow { display: grid; grid-template-columns: 210px 1fr 62px; gap: 0.6rem;
  align-items: center; padding: 0.3rem 0 0; font-size: 0.84rem; }
.crow .clab { text-align: right; }
.cbar { height: 14px; background: var(--edge); border-radius: 2px; overflow: hidden; }
.cbar span { display: block; height: 100%; }
.cnum { text-align: right; font-family: var(--mono); font-size: 0.8rem; }
.crow .csub { grid-column: 2 / 4; font-size: 0.74rem; padding-bottom: 0.3rem; }
.bites { display: flex; flex-direction: column; gap: 0.3rem; }
.bite { display: grid; grid-template-columns: 15rem 1fr 6rem 5rem; gap: 0.6rem;
  align-items: center; font-size: 12px; }
.bname { overflow-wrap: anywhere; }
.bbar { display: flex; height: 14px; border-radius: 3px; background: var(--edge);
  overflow: hidden; min-width: 0; }
.bbar > i { display: block; height: 100%; transition: width 420ms ease-out; }
.bbar > .before { background: var(--identified); opacity: 0.22; }
.bbar > .mine { background: var(--identified); }
.bnum, .bsize { text-align: right; font-variant-numeric: tabular-nums; }
@media (max-width: 680px) {
  .bite { grid-template-columns: 1fr 4.5rem; grid-template-areas: "n n" "b b" "v s"; }
  .bname { grid-area: n; } .bbar { grid-area: b; }
  .bnum { grid-area: v; text-align: left; } .bsize { grid-area: s; }
}
.ditem { border-top: 1px solid var(--edge); padding: 0.55rem 0; }
.gline { border-top: 1px solid var(--edge); padding: 0.5rem 0; line-height: 1.5; }
.maplegend { display: flex; flex-wrap: wrap; gap: 0.35rem 1.1rem; font-size: 12px;
  color: var(--dim); margin: 0.1rem 0 0.7rem; }
.maplegend span { display: inline-flex; align-items: center; gap: 0.35rem; }
.maplegend i { width: 11px; height: 11px; border-radius: 3px; display: inline-block; }
.gmap { border-top: 1px solid var(--edge); padding: 0.7rem 0 0.9rem; }
.gmaphead { font-size: 14.5px; margin-bottom: 0.1rem; }
.gintro { font-size: 13.5px; line-height: 1.55; margin: 0.1rem 0 0.6rem; }
.gsay { font-size: 13.5px; line-height: 1.55; color: var(--fg); margin-top: 0.15rem; }
.gline:first-of-type { border-top: 0; }
.ditem > summary { font-size: 13.5px; color: var(--fg); padding: 0.15rem 0; }
.ditem > summary b { font-weight: 600; }
.ditem > summary .dim { font-weight: 400; }
.ditem[open] > summary { margin-bottom: 0.5rem; }
.inner { display: flex; flex-direction: column; gap: 0.8rem; min-width: 0;
  padding: 0.2rem 0 0.9rem; }
.big { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 0.9rem; margin: 1rem 0 0.4rem; }
.bignum { font-size: 30px; font-weight: 650; line-height: 1.1;
  font-feature-settings: "tnum"; }
.biglab { font-size: 12.5px; color: var(--dim); margin-top: 0.25rem; }
.lede { font-size: 15px; max-width: 70ch; line-height: 1.55; }
/* ---- public-health front page: verdict heroes, ranked lists, donuts ---------
   Every ranked row and every donut slice carries a direct text label, because the
   validator's one light-mode contrast warning on --v3 is an obligation, not a note:
   identity here is never colour alone. */
.verdict { margin: 1rem 0 0.4rem; padding: 0.95rem 1.05rem; border-radius: 8px;
  background: var(--v-faint); border-left: 3px solid var(--v2); }
.verdict .vmain { font-size: 1.15rem; line-height: 1.35; text-wrap: balance; }
.verdict .vmain b { font-weight: 650; }
.verdict .vsub { font-size: 12.5px; color: var(--dim); margin-top: 0.4rem;
  max-width: 68ch; }
.ranked { display: flex; flex-direction: column; gap: 0.3rem; margin: 0.9rem 0 0.2rem; }
.rrow { display: grid; grid-template-columns: 1.6rem minmax(12rem, 24rem) 1fr auto;
  gap: 0.55rem; align-items: center; font-size: 13px; }
.rno { font-family: var(--mono); font-size: 0.72rem; color: var(--dim);
  text-align: right; }
.rno.wide { font-size: 0.78rem; color: var(--fg); }
.rrow:has(.rno.wide) { grid-template-columns: 3.2rem minmax(12rem, 24rem) 1fr auto; }
.rname { overflow-wrap: anywhere; }
.rname .rtag { font-size: 0.68rem; color: var(--dim); margin-left: 0.35rem; }
.rbar { height: 9px; min-width: 0; background: var(--v-faint); border-radius: 4px; }
.rbar > i { display: block; height: 100%; border-radius: 0 4px 4px 0;
  background: var(--v2); transition: width 500ms ease-out; }
.rval { font-family: var(--mono); font-size: 0.74rem; text-align: right;
  font-variant-numeric: tabular-nums; white-space: nowrap; }
.chip { display: inline-block; font-size: 0.66rem; padding: 0.05rem 0.35rem;
  border-radius: 3px; border: 1px solid var(--edge); color: var(--dim);
  margin-left: 0.35rem; vertical-align: 1px; }
.chip.strong { border-color: var(--v2); color: var(--v2); }
.more { margin: 0.5rem 0 0; }
.more > summary { font-size: 12.5px; }
.misses { display: flex; flex-wrap: wrap; gap: 0.3rem; margin-top: 0.5rem; }
.miss { font-size: 11.5px; color: var(--dim); border: 1px solid var(--edge);
  border-radius: 3px; padding: 0.1rem 0.4rem; }
.donuts { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
  gap: 1.2rem; margin: 1rem 0 0.3rem; }
.donut { text-align: center; }
.donut svg { width: 100%; max-width: 200px; height: auto; }
.donut .dhero { font-family: var(--mono); font-weight: 500; }
.dcap { font-size: 12.5px; color: var(--dim); margin-top: 0.35rem; }
.donutlg { display: flex; flex-direction: column; gap: 0.25rem; margin-top: 0.6rem;
  font-size: 12px; text-align: left; max-width: 15rem; margin-inline: auto; }
.donutlg > div { display: grid; grid-template-columns: 0.7rem 1fr auto; gap: 0.45rem;
  align-items: center; }
.donutlg i { display: block; width: 0.7rem; height: 0.7rem; border-radius: 2px; }
.donutlg span.v { font-family: var(--mono); font-size: 0.72rem;
  font-variant-numeric: tabular-nums; }
.split { display: flex; height: 26px; border-radius: 4px; overflow: hidden;
  margin: 1.1rem 0 0.8rem; gap: 2px; background: var(--v-faint); }
.split > i { display: block; min-width: 2px; transform-origin: left;
  animation: splitin 520ms ease-out both; }
@keyframes splitin { from { transform: scaleX(0); } to { transform: scaleX(1); } }
@media (prefers-reduced-motion: reduce) { .split > i { animation: none; } }
.splitlg { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
  gap: 0.6rem 1rem; font-size: 12.5px; }
.splitlg > div { display: grid; grid-template-columns: 0.7rem 1fr; gap: 0.5rem; }
.splitlg i { display: block; width: 0.7rem; height: 0.7rem; border-radius: 2px;
  margin-top: 0.25rem; }
.splitlg .dim { font-size: 11.5px; margin-top: 0.1rem; }
.arc { stroke-dasharray: var(--len) var(--gap);
  stroke-dashoffset: calc(var(--off) + var(--len));
  animation: arcin 700ms ease-out forwards; animation-delay: var(--d, 0ms); }
@keyframes arcin { to { stroke-dashoffset: var(--off); } }
@media (prefers-reduced-motion: reduce) {
  .arc { animation: none; stroke-dashoffset: var(--off); }
}
@media (max-width: 680px) {
  .rrow, .rrow:has(.rno.wide) { grid-template-columns: 2.6rem minmax(0, 1fr) auto;
    row-gap: 0.15rem; margin-bottom: 0.35rem; }
  .rrow .rname { grid-column: 2 / -1; }
  .rrow .rbar { grid-column: 2 / 3; }
  .donutlg { max-width: none; }
}
a { color: var(--identified); }
/* ---- citations ------------------------------------------------------------
   Every headline figure carries a link to the entry that states how it was
   computed. The link is an anchor and the highlight is :target, so a citation
   works with scripting disabled and cannot drift from its target. */
.cite { font-family: var(--mono); font-size: 0.68em; vertical-align: super;
  text-decoration: none; color: var(--v1); padding: 0 0.08em; }
.cite:hover { text-decoration: underline; }
.meth { border-top: 1px solid var(--edge); padding: 0.7rem 0;
  display: grid; grid-template-columns: 1.7rem minmax(0, 1fr); gap: 0 0.4rem;
  scroll-margin-top: 1rem; }
.meth .mno { font-family: var(--mono); font-size: 12px; color: var(--dim); }
.meth h3 { font-size: 13.5px; font-weight: 600; margin: 0 0 0.35rem; }
.mgrid { display: grid; grid-template-columns: 6.6rem minmax(0, 1fr);
  gap: 0.2rem 0.7rem; margin: 0; font-size: 12.5px; }
.mgrid dt { color: var(--dim); }
.mgrid dd { margin: 0; }
.meth:target { background: var(--v2-soft); border-radius: 6px;
  box-shadow: 0 0 0 0.5rem var(--v2-soft); }
@media (max-width: 680px) {
  .meth { grid-template-columns: 1.3rem minmax(0, 1fr); }
  .mgrid { grid-template-columns: minmax(0, 1fr); gap: 0 0; }
  .mgrid dt { margin-top: 0.3rem; }
}
/* ---- THE PLAYGROUND ------------------------------------------------------
   Controls sit with the figure they move, so a reader never has to work out
   which panel a slider belongs to, and the ones nobody moves sit behind one
   disclosure. The strip at the top of the page carries no control of its own:
   it says what the current reading is set to, hands back a link to it, and
   returns everything to the shipped defaults in one click. */
.play { border: 1px solid var(--edge); border-radius: 10px; padding: 0.6rem 0.8rem;
  margin: 0 0 1rem; display: flex; flex-wrap: wrap; gap: 0.5rem 1rem;
  align-items: center; justify-content: space-between; background: var(--panel); }
.playline { font-size: 12.5px; min-width: 14rem; flex: 1 1 22rem; }
.playline b { font-weight: 600; }
.playbtns { display: flex; gap: 0.4rem; flex-wrap: wrap; }
.play button, .ctls button { font: inherit; font-size: 12px; padding: 0.22rem 0.5rem;
  border: 1px solid var(--edge); border-radius: 6px; background: transparent;
  color: var(--fg); cursor: pointer; }
.play button:hover, .ctls button:hover { border-color: var(--identified); }
.ctls { display: flex; flex-wrap: wrap; gap: 0.35rem 0.9rem; align-items: center;
  margin: 0.2rem 0 0.5rem; font-size: 12.5px; }
.ctl { display: flex; flex-wrap: wrap; gap: 0.3rem; align-items: center; }
.ctl input[type=range] { width: 11rem; max-width: 46vw; accent-color: var(--identified); }
.clab { color: var(--dim); }
.cval { font-family: var(--mono); font-size: 11.5px; }
.cship { color: var(--dim); }
.cmoved { color: var(--warn); }
.cbtns { display: flex; gap: 0.25rem; flex-wrap: wrap; }
.cbtns button.on { border-color: var(--identified); background: var(--v2-soft); }
.ctl-absent { display: block; font-size: 12px; color: var(--dim);
  border-left: 2px solid var(--warn); padding: 0.15rem 0 0.15rem 0.5rem;
  flex-basis: 100%; }
.cmore { flex-basis: 100%; }
.cmore > summary { cursor: pointer; color: var(--dim); font-size: 12px; }
.cmore[open] { border-left: 1px solid var(--edge); padding-left: 0.6rem; }
.cmore .ctl { margin: 0.3rem 0; }
.setting { font-family: var(--mono); font-size: 11.5px; color: var(--dim);
  margin: 0.3rem 0 0; }
.setting b { font-family: var(--sans); color: var(--fg); }
.relabel { color: var(--warn); }
:focus-visible { outline: 2px solid var(--identified); outline-offset: 2px; }
@media (prefers-reduced-motion: reduce) { * { transition: none !important; } }
</style>
<div class="wrap" id="app"></div>
<script id="bundle" type="application/json">__BUNDLE__</script>
<script>
const D = JSON.parse(document.getElementById('bundle').textContent);
const app = document.getElementById('app');
const esc = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const pct = x => (x == null ? '-' : (100 * x).toFixed(2) + '%');
const num = (x, d = 3) => (x == null ? '-' : (+x).toFixed(d));
const gb = b => b == null ? '-' : b >= 1e9 ? (b/1e9).toFixed(2)+' GB' : b >= 1e6 ? (b/1e6).toFixed(1)+' MB' : b+' B';

function el(tag, cls, html) { const e = document.createElement(tag);
  if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; }

function panel(title, sub) { const p = el('section', 'panel');
  p.append(el('h2', null, esc(title)));
  if (sub) p.append(el('p', 'note', sub));
  return p; }

/* ---- THE PLAYGROUND -------------------------------------------------------
   SPEC 13.3. Every control here re-filters and re-sorts numbers the run already
   measured; nothing on this page recomputes a measurement, and nothing can be read
   at a setting the run did not reach. Four things the store guarantees, because
   each one exists because of a reading this project has had to withdraw:

     - A control's range and its measured points come from the bundle, which
       experiments/web_controls.py derived from the run's own files. The page never
       invents a point, so a slider cannot run past its measurement.
     - A discrete control snaps to those points and never interpolates. Kraken2
       confidence was run at 0, 0.05, 0.1 and 0.2, and a factor of twelve in removal
       lives inside that gap.
     - Every reading carries its setting, the shipped default is marked, and the
       state is in the URL, so a reading taken here is citable and one click returns.
     - A control changes what is shown and never what was computed. Where moving one
       re-labels material, the label re-renders with it and says both: at your
       setting, and as shipped. */
const CTLDROPPED = [];
const CTLDEF = D.controls || {controls: [], fixed: []};
const CTLBY = {};
(CTLDEF.controls || []).forEach(c => { CTLBY[c.id] = c; });
const CTLV = {};
const CTLSUBS = [];

/* Two controls are a pair: the threshold's measured points are the selected tier's,
   so they are read off that tier rather than held twice. Two copies of a measured
   grid is two instruments. */
function ctlGrid(id) {
  const c = CTLBY[id];
  if (!c || c.absent) return null;
  if (c.points_from_control) {
    const g = ((CTLBY[c.points_from_control] || {}).grid || {})[ctlGet(c.points_from_control)];
    return g ? g.points : null;
  }
  return c.points || null;
}
function ctlShipped(id) {
  const c = CTLBY[id];
  if (!c) return null;
  if (c.points_from_control) {
    const g = ((CTLBY[c.points_from_control] || {}).grid || {})[ctlGet(c.points_from_control)];
    return g ? g.shipped : null;
  }
  return c.shipped;
}
function ctlGet(id) {
  return Object.prototype.hasOwnProperty.call(CTLV, id) ? CTLV[id] : ctlShipped(id);
}
function ctlMoved(id) {
  const c = CTLBY[id];
  if (!c || c.absent) return false;
  return String(ctlGet(id)) !== String(ctlShipped(id));
}
function ctlValid(id, v) {
  const c = CTLBY[id];
  if (!c || c.absent) return false;
  const g = ctlGrid(id);
  if (g) return g.some(p => String(p.value) === String(v));
  const b = c.bounds;
  return !!b && typeof v === 'number' && v >= b.min - 1e-9 && v <= b.max + 1e-9;
}
function ctlFmt(id, v) {
  const c = CTLBY[id] || {};
  const g = ctlGrid(id);
  if (g) {
    const p = g.filter(q => String(q.value) === String(v))[0];
    return p ? String(p.label) : String(v);
  }
  if (v == null) return 'not set';
  if (c.unit === 'fraction') return (100 * v).toFixed(0) + '%';
  if (c.unit === 'x') return num(v, 0) + 'x';
  return num(v, c.unit === 'aa' || c.unit === 'nt' || c.unit === 'bits' ||
             c.unit === 'samples' ? 0 : 2) + (c.unit ? ' ' + c.unit : '');
}
/* The setting string. It goes under every controlled figure, so a screenshot of one
   carries the threshold it was taken at: a fraction without its operating point is
   unreadable, including the ones this project has quoted at itself. */
function ctlSetting(ids) {
  return ids.filter(id => CTLBY[id] && !CTLBY[id].absent).map(id =>
    CTLBY[id].label + ' ' + ctlFmt(id, ctlGet(id)) +
    (ctlMoved(id) ? ' (moved, shipped ' + ctlFmt(id, ctlShipped(id)) + ')'
                  : ' (shipped)')).join('   |   ');
}
function ctlAnyMoved(ids) { return ids.some(ctlMoved); }
function settingLine(ids) {
  const d = el('div', 'setting');
  // `Read at:` AND NOTHING ELSE IS NOT AN OPERATING POINT. Where every control on a
  // figure is absent on this run, the join is empty and the page printed a bold label
  // with nothing after it. The absences are the answer, so they are what it prints.
  const paint = () => {
    const said = ctlSetting(ids);
    if (said) { d.innerHTML = '<b>Read at:</b> ' + esc(said); return; }
    const why = ids.filter(id => CTLBY[id] && CTLBY[id].absent)
                   .map(id => CTLBY[id].label + ': ' + CTLBY[id].absent);
    d.innerHTML = '<b>Read at:</b> the shipped setting. ' +
      (why.length ? esc('No control on this figure can move on this run. ' +
                        why.join('. ') + '.')
                  : esc('This figure has no control.'));
  };
  paint(); ctlOn(paint);
  return d;
}
function ctlURL() {
  const on = Object.keys(CTLBY).filter(ctlMoved).map(id => id + '~' + ctlGet(id));
  return on.length ? '#c=' + encodeURIComponent(on.join(';')) : '';
}
function ctlWriteURL() {
  try {
    const h = ctlURL();
    history.replaceState(null, '', h || (location.pathname + location.search));
  } catch (e) { /* file:// in some browsers refuses replaceState; the page still works */ }
}
function ctlReadURL() {
  const m = /[#&]c=([^&]*)/.exec(location.hash || '');
  if (!m) return;
  decodeURIComponent(m[1]).split(';').forEach(kv => {
    const i = kv.indexOf('~');
    if (i < 0) return;
    const id = kv.slice(0, i), raw = kv.slice(i + 1);
    const v = /^-?[0-9]*[.]?[0-9]+(e[-+]?[0-9]+)?$/i.test(raw) ? +raw : raw;
    // A setting the run cannot reach is dropped rather than clamped: a link that
    // silently moves to the nearest measured point is a citation to a reading nobody
    // took. The strip says when a link was dropped.
    if (ctlValid(id, v)) CTLV[id] = v; else CTLDROPPED.push(id + '~' + raw);
  });
}
function ctlOn(fn) { CTLSUBS.push(fn); }
function ctlSet(id, v) {
  if (!ctlValid(id, v)) return;
  CTLV[id] = v;
  // The paired threshold follows its tier: the point that was selected on one tier's
  // grid is meaningless on another's, so it returns to that tier's shipped point.
  Object.keys(CTLBY).forEach(k => {
    if ((CTLBY[k] || {}).points_from_control === id) delete CTLV[k];
  });
  ctlWriteURL();
  CTLSUBS.forEach(f => f(id));
}
function ctlReset() {
  Object.keys(CTLV).forEach(k => { delete CTLV[k]; });
  ctlWriteURL();
  CTLSUBS.forEach(f => f(null));
}

function ctlWidget(id) {
  const c = CTLBY[id];
  const box = el('div', 'ctl');
  if (!c) {
    box.className = 'ctl ctl-absent';
    box.setAttribute('data-ctl', id);
    box.textContent = 'the page asks for a control called ' + id + ' that nothing ' +
      'declared, so there is nothing to move';
    return box;
  }
  box.setAttribute('data-ctl', id);
  if (c.absent) {
    box.className = 'ctl ctl-absent';
    box.innerHTML = '<b>' + esc(c.label) + '</b> cannot move on this run: ' +
      esc(c.absent);
    return box;
  }
  box.append(el('span', 'clab', esc(c.label)));
  const g = ctlGrid(id);
  const out = el('span', 'cval');
  let range = null, btnRow = null, builtFor = null;
  if (g && g.length <= 6) {
    btnRow = el('span', 'cbtns');
    // Rebuilt whenever the grid underneath changes, because a paired control's points
    // are the selected tier's: buttons left over from another tier's sweep would offer
    // settings that tier never measured, which is the one thing a snapping control is
    // for. The grid is read from the bundle each time, never cached here.
    const build = () => {
      const gg = ctlGrid(id) || [];
      const sig = gg.map(p => p.value).join('|');
      if (sig === builtFor) return;
      builtFor = sig;
      btnRow.innerHTML = '';
      gg.forEach(p => {
        const b = el('button', null, esc(p.label));
        b.type = 'button';
        b.setAttribute('data-v', p.value);
        b.onclick = () => ctlSet(id, p.value);
        btnRow.append(b);
      });
    };
    build();
    ctlOn(build);
    box.append(btnRow);
  } else {
    range = document.createElement('input');
    range.type = 'range';
    range.setAttribute('aria-label', c.label);
    if (g) { range.min = 0; range.max = g.length - 1; range.step = 1;
             range.oninput = () => { const gg = ctlGrid(id);
               ctlSet(id, gg[Math.min(+range.value, gg.length - 1)].value); }; }
    else { range.min = c.bounds.min; range.max = c.bounds.max;
           range.step = c.bounds.step; range.oninput = () => ctlSet(id, +range.value); }
    box.append(range);
  }
  box.append(out);
  const paint = () => {
    const v = ctlGet(id), gg = ctlGrid(id);
    out.innerHTML = esc(ctlFmt(id, v)) + (ctlMoved(id)
      ? ' <span class="cmoved">moved, shipped ' + esc(ctlFmt(id, ctlShipped(id))) +
        '</span>'
      : ' <span class="cship">shipped</span>');
    Array.prototype.forEach.call(box.querySelectorAll('.cbtns button'), b => {
      b.classList.toggle('on', String(b.getAttribute('data-v')) === String(v));
      // The shipped point is named in the readout rather than marked on the button: a
      // mark a reader has to decode is one more thing between them and the number.
      if (String(b.getAttribute('data-v')) === String(ctlShipped(id))) {
        b.title = 'the setting this run shipped';
      }
      b.setAttribute('aria-pressed', String(b.getAttribute('data-v')) === String(v));
    });
    if (range) {
      if (gg) {
        range.max = gg.length - 1;
        const i = gg.map(p => String(p.value)).indexOf(String(v));
        if (i >= 0) range.value = i;
      } else { range.value = v; }
    }
  };
  paint(); ctlOn(paint);
  return box;
}

/* A handful visible, the rest behind one disclosure, in the same summary-then-detail
   shape as the panels. A control the reader never moves should be invisible; a control
   that cannot move honestly should be visible as a stated constraint, which is why the
   absent ones are rendered rather than dropped. */
function ctlRow(ids) {
  const d = el('div', 'ctls');
  const head = ids.filter(id => CTLBY[id] && !CTLBY[id].absent &&
                                CTLBY[id].placement === 'headline');
  const rest = ids.filter(id => head.indexOf(id) < 0);
  head.forEach(id => d.append(ctlWidget(id)));
  if (rest.length) {
    const det = el('details', 'cmore');
    const s = document.createElement('summary');
    s.textContent = 'More controls (' + rest.length + ')';
    det.append(s);
    rest.forEach(id => det.append(ctlWidget(id)));
    d.append(det);
  }
  return d;
}
// The URL is read once, before any panel draws, so a cited link opens at the reading
// it names rather than at the defaults and then jumping.
ctlReadURL();

/* A control whose null exists only at the shipped setting. Moving it renders the words
   rather than a clean-looking number: the withdrawn readings of this project are what a
   threshold quoted without its own null looks like. */
function nullLine(ids, atShipped) {
  const d = el('div', 'note');
  const paint = () => {
    d.innerHTML = ctlAnyMoved(ids)
      ? '<b>Null not measured at this setting.</b> ' + esc(atShipped.moved)
      : '<b>Null, at the shipped setting.</b> ' + esc(atShipped.shipped);
  };
  paint(); ctlOn(paint);
  return d;
}

/* ---- GENE MAPS ------------------------------------------------------------
   A genome drawn as the molecule, shared by every panel that draws one, so the
   colour of a polymerase cannot mean one thing in the summary and another in the
   holdout panel. Roles carry a plain-English gloss because a reader who does not
   know what a maturation protein is cannot read a legend that only names it. */
/* Built from D.vocabulary, never written here. A gloss the renderer invents is a
   biological claim with no provenance sitting beside numbers that have one. */
const VOC = D.vocabulary || {roles: {}, architectures: {}};
const GENE_ROLE = {};
Object.keys(VOC.roles || {}).forEach(k => {
  GENE_ROLE[k] = ['var(' + VOC.roles[k].colour_token + ')', VOC.roles[k].plain];
});
if (!GENE_ROLE.other) GENE_ROLE.other = ['var(--null)', 'role not called'];
const geneRoleOf = o => (o && GENE_ROLE[o.role] ? o.role : 'other');

/* The architecture call for a contig, from the file experiments/architecture_call.py
   wrote. Returns null where the rule declined, which is a census entry and not a verdict
   of "not a phage". */
const ARCH = {};
((D.architecture || {}).rows || []).forEach(r => { ARCH[r.contig] = r; });
const archOf = c => ARCH[c] || null;
function archWords(a) {
  if (!a || !a.architecture) return null;
  return (VOC.architectures || {})[a.architecture] || null;
}

function geneLegend(keys) {
  const lg = el('div', 'maplegend');
  lg.innerHTML = (keys || ['polymerase', 'coat', 'maturation', 'other']).map(k =>
    '<span><i style="background:' + GENE_ROLE[k][0] + '"></i>' + k + ' &mdash; ' +
    GENE_ROLE[k][1] + '</span>').join('');
  return lg;
}

/* nt: the contig length. orfs: {s, e, d, aa, role, pident, title, pfam, sig}.
   Returns an element, scrollable, so a long genome keeps its scale on a phone
   rather than being squeezed until the short genes vanish. */
function geneMap(nt, orfs, rightLabel) {
  const MW = 1000, ML = 8, MR = 8, MH = 74, TRK = 46;
  const mx = v => ML + (v / nt) * (MW - ML - MR);
  let g = '<svg viewBox="0 0 ' + MW + ' ' + MH + '" width="100%" ' +
    'style="min-width:520px" role="img" aria-label="gene map"><line x1="' + mx(0) +
    '" y1="' + TRK + '" x2="' + mx(nt) + '" y2="' + TRK +
    '" stroke="var(--edge)" stroke-width="3"/>';
  orfs.slice().sort((a, b) => a.s - b.s).forEach(o => {
    const k = geneRoleOf(o), x0 = mx(Math.min(o.s, o.e)), x1 = mx(Math.max(o.s, o.e));
    const w = Math.max(3, x1 - x0), tip = Math.min(10, w / 3);
    // an arrow, not a rectangle: strand is a fact about the molecule
    const pts = o.d < 0
      ? [[x0, TRK], [x0 + tip, TRK - 11], [x0 + w, TRK - 11], [x0 + w, TRK + 11],
         [x0 + tip, TRK + 11]]
      : [[x0, TRK - 11], [x0 + w - tip, TRK - 11], [x0 + w, TRK],
         [x0 + w - tip, TRK + 11], [x0, TRK + 11]];
    g += '<polygon points="' + pts.map(q => q[0].toFixed(1) + ',' + q[1]).join(' ') +
      '" fill="' + GENE_ROLE[k][0] + '" opacity="' + (o.sig === false ? 0.45 : 1) +
      '"><title>' + esc(o.role || 'role not called') + ', ' + o.aa + ' aa' +
      (o.pident != null ? ', ' + num(o.pident, 1) + '% to ' + esc(o.title || '')
        : ', no significant protein hit') +
      ((o.pfam || []).length ? '\\\ndomain ' + esc(o.pfam.join(', ')) : '') +
      '</title></polygon>';
    if (o.pident != null && w > 74) {
      g += '<text x="' + (x0 + w / 2) + '" y="' + (TRK + 4) + '" text-anchor="middle" ' +
        'style="font-size:12px;font-weight:600;fill:var(--panel);' +
        'font-family:var(--sans)">' + num(o.pident, 0) + '%</text>';
    }
    if (w > 54) {
      g += '<text x="' + (x0 + w / 2) + '" y="' + (TRK - 18) +
        '" text-anchor="middle" class="ax">' + esc(k === 'other'
          ? (o.pident != null ? 'other' : 'no hit') : k) + '</text>';
    }
  });
  g += '<text x="' + mx(0) + '" y="' + (TRK + 26) + '" class="ax">0</text>';
  g += '<text x="' + mx(nt) + '" y="' + (TRK + 26) + '" text-anchor="end" class="ax">' +
    esc(rightLabel || ((nt / 1000).toFixed(1) + ' kb')) + '</text></svg>';
  const hostd = el('div', 'scroll'); hostd.innerHTML = g;
  return hostd;
}



/* ---- METHODS: ONE ENTRY PER HEADLINE FIGURE -------------------------------
   The page states numbers at the top and how they were obtained at the bottom, and
   the two are joined by a link rather than by the reader remembering. Each entry
   carries the six things that make a figure readable here: its input, the procedure,
   the threshold, the denominator, the control arm, and the direction of the bias that
   the control does not remove. The entries are data, from results/web/methods.json,
   resolved for this page's sample by make_web.py: a row measured on another library
   reads "Not measured on this library" here rather than borrowing that library's number.
   The numbering is positional and the same on every page. */
const METHODS = D.methods || [];
const MIDX = {};
METHODS.forEach((m, i) => { MIDX[m.id] = i + 1; });
const cite = id => '<a class="cite" href="#m-' + id + '" title="' +
  (MIDX[id] ? 'method ' + MIDX[id] : 'method') + '">[' + (MIDX[id] || '?') + ']</a>';


/* A PILE FANNING INTO ITS PARTS, WITH ITS CONTROL DRAWN THE SAME WAY UNDERNEATH.
   Two of this page's measurements are splits of one pile into named parts, and both of
   them only mean anything against the same split of material we can already name - this
   repo has had three results turn out to be the control moving with the test. Drawing
   both rows from one function is the cheapest way to keep them honestly comparable:
   neither row can quietly get a different scale, a different order or a different
   colour from the other. The interval is drawn where the file carries one, because a
   share without its interval is the thing the spec forbids. */
function fanFig(rows, opts) {
  const W = 920, RH = (opts && opts.rowHeight) || 190, GAP = 42;
  const srcX = 18, srcW = 64, dstX = 452, dstW = 64, labX = dstX + dstW + 10;
  const H = rows.length * RH + (rows.length - 1) * GAP + 30;
  const keys = rows[0].slices.map(s => s.k);
  let g = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' +
    esc((opts && opts.alt) || 'a pile fanning into its parts') + '">';
  let delay = 0;
  rows.forEach((row, ri) => {
    const y0 = 20 + ri * (RH + GAP);
    const tot = row.slices.reduce((a, c) => a + c.v, 0) || 1;
    // A destination block never falls below MIN, or a 0.17% pile has no bar to label and
    // the label lands on the pile above it. The floor is paid for out of the blocks that
    // can afford it rather than out of the figure's height, so the stack still ends where
    // the source bar ends and the two rows stay on one scale.
    const MIN = 12, nz = row.slices.filter(s => s.v > 0).length || 1;
    let pool = RH - 10 * (nz - 1), rest = tot, hd = {};
    row.slices.forEach(sl => {
      if (sl.v > 0 && pool * sl.v / (rest || 1) < MIN) {
        hd[sl.k] = MIN; pool -= MIN; rest -= sl.v;
      }
    });
    row.slices.forEach(sl => {
      if (sl.v > 0 && hd[sl.k] == null) hd[sl.k] = Math.max(MIN, pool * sl.v / (rest || 1));
    });
    let ys = y0, yd = y0;
    const dim = row.control ? 0.42 : 1;
    g += '<text x="' + srcX + '" y="' + (y0 - 7) + '" font-size="11.5" ' +
      'fill="var(--fg)" font-weight="500">' + esc(row.title) +
      ' <tspan font-family="var(--mono)" font-weight="400" fill="var(--dim)">' +
      row.total.toLocaleString() + '</tspan></text>';
    g += '<rect x="' + srcX + '" y="' + y0 + '" width="' + srcW + '" height="' + RH +
      '" rx="3" fill="var(--edge)" opacity="' + (row.control ? 0.5 : 0.9) + '"/>';
    row.slices.forEach(sl => {
      const hs = RH * sl.v / tot;
      const h = hd[sl.k] || 0;
      if (sl.v > 0) {
        const col = sl.col || 'var(--null)';
        const x1 = srcX + srcW, x2 = dstX;
        const cx = (x1 + x2) / 2;
        g += '<path class="fanrib" style="animation-delay:' + delay + 'ms" d="M' + x1 +
          ' ' + ys.toFixed(1) + ' C' + cx + ' ' + ys.toFixed(1) + ' ' + cx + ' ' +
          yd.toFixed(1) + ' ' + x2 + ' ' + yd.toFixed(1) + ' L' + x2 + ' ' +
          (yd + h).toFixed(1) + ' C' + cx + ' ' + (yd + h).toFixed(1) + ' ' + cx +
          ' ' + (ys + hs).toFixed(1) + ' ' + x1 + ' ' + (ys + hs).toFixed(1) +
          ' Z" fill="' + col + '" opacity="' + (0.22 * dim + 0.04) + '"/>';
        g += '<rect x="' + dstX + '" y="' + yd.toFixed(1) + '" width="' + dstW +
          '" height="' + h.toFixed(1) + '" rx="2" fill="' + col + '" opacity="' +
          (0.95 * dim) + '"><title>' + esc(sl.lab + ': ' + sl.v.toLocaleString() +
          ' of ' + tot.toLocaleString() + ' (' + (100 * sl.v / tot).toFixed(2) + '%)') +
          '</title></rect>';
        const pct = (100 * sl.v / tot);
        let txt = esc(sl.lab) + ' <tspan font-family="var(--mono)" fill="var(--fg)">' +
          pct.toFixed(pct < 10 ? 2 : 1) + '%</tspan>';
        if (sl.ci) txt += ' <tspan font-size="10" fill="var(--dim)">' +
          (100 * sl.ci[0]).toFixed(1) + ' to ' + (100 * sl.ci[1]).toFixed(1) + '</tspan>';
        g += '<text x="' + labX + '" y="' + (yd + h / 2 + 4).toFixed(1) +
          '" font-size="11" fill="var(--dim)">' + txt + '</text>';
        delay += 70;
      }
      ys += hs;
      yd += h + (sl.v > 0 ? 10 : 0);
    });
  });
  g += '</svg>';
  const holder = el('div', 'scroll flowwrap');
  holder.innerHTML = g;
  return holder;
}

/* WHICH PILE THESE ROWS CAME FROM, on the panel's face. The rows of a search over the
   whole residue and a search over a 6,000-pair sample of it look identical, and they
   differ by more than two orders of magnitude in what they cover - so a reader who is
   not told assumes the wider one. The denominator is printed with it, because that is
   the number the shares on these rows are taken over. */
function scopeLine(R) {
  const sc = (R && R.scope) || {};
  const res = ((R && R.piles) || {}).residue || {};
  const bits = [];
  if (sc.residue_pairs) {
    bits.push('<b>Every unplaced pair in the library</b>: ' + num(sc.residue_pairs, 0) +
      ' pairs, ' + num(sc.residue_mates, 0) + ' mates searched, against a control of ' +
      num(sc.control_pairs, 0) + ' already-placed pairs through the same search in the ' +
      'same pass.');
  } else {
    bits.push('<b>The per-read report\\'s block sample of the residue</b>: ' +
      num(res.searched, 0) + ' mates, not the whole pile.');
  }
  if (sc.bias) bits.push(esc(sc.bias));
  const e = el('div', 'callout');
  e.innerHTML = bits.join(' ');
  return e;
}

/* No figure without its operating point: this is the only way a panel gets its chips,
   and a missing value renders as the words rather than as a clean blank. */
function opChips(obj, keys) {
  const d = el('div', 'op');
  let any = false;
  (keys || Object.keys(obj || {})).forEach(k => {
    const v = (obj || {})[k];
    if (v == null) return;
    any = true;
    d.append(el('span', 'chip', esc(k.replace(/_/g, ' ')) + ' <b>' + esc(v) + '</b>'));
  });
  if (!any) d.append(el('span', 'absent', 'operating point absent from the run directory'));
  return d;
}

/* No pile without its interval. */
function interval(iv) {
  if (!iv) return '<span class="absent">no interval</span>';
  return '<span class="mono">' + pct(iv.point) + '</span> <span class="dim mono">[' +
    pct(iv.lo) + ', ' + pct(iv.hi) + ']</span>';
}

/* WHETHER A TIER RETURNS A NAME IS A PROPERTY OF THE INDEX, NOT OF ITS DIRECTORY NAME.
   Two panels each carried a literal map of tier names that can name, keyed `k2_viral`,
   `k2_standard08` and so on. Production pulls its Kraken2 databases through
   `hapax fetch-db`, which names the directory after the release it downloaded, so the
   tiers arrive as `k2_viral_20261022.d` and matched nothing: every Kraken2 read moved out
   of "named" and into "recognised, unnamed", the named band read 2.24% where the cascade
   had named 4.07%, and the reconciliation panel went on asserting that its two named rows
   agreed while they stood 15 points apart. Classification is by what the tier IS, and a
   tier this function does not recognise is returned as unknown so that a new index lands
   in a stated category rather than quietly in the larger one. */
function tierKind(name) {
  const n = String(name || '');
  if (/^k2_/.test(n) || /^prot_/.test(n)) return 'names';
  if (n === 'rrna' || n === 'artefact' || /^selfasm_/.test(n)) return 'recognises';
  if (/^observe_/.test(n)) return 'observes';
  return 'unknown';
}

function tierCanName(name) { return tierKind(name) === 'names'; }

// ---------- header
(function () {
  const r = D.run || {};
  const h = el('header', null, '');
  const inp = (r.inputs || {});
  // Name the sample, not the tool. Whoever opens this is looking at one water on one
  // day, and the accession is the only thing that lets them find it again or compare it
  // with another; "run viewer" told them which program wrote the page.
  const acc = ((inp.r1 || '').match(/([A-Z]{3}\\d{6,})/) || [])[1];
  // THE WATER IS NAMED BY ITS OWN FACTS, NEVER BY THIS FILE. The heading read
  // `<accession>, untreated wastewater` and the line under it `Untargeted RNA sequencing
  // of untreated wastewater`, both literals earned by the discriminating dataset. On
  // the third production library the archive does not say influent or effluent - that
  // library's own basis says so in as many words, and the ribosomal fraction differs by
  // two orders of magnitude between the two - so the page's first sentence made a claim
  // its own facts refused to make. `page.derived.facts` is the block, each fact with how it was
  // established, and a fact that is absent is left out rather than guessed.
  const facts = ((D.page || {}).derived || {}).facts || {};
  const matrix = facts.matrix || null;
  h.append(el('h1', null, acc ? esc(acc) + (matrix ? ', ' + esc(matrix) : '')
                              : 'Hapax run viewer'));
  h.append(el('p', 'lede',
    esc(((D.page || {}).sample || {}).lede || '') +
    ' Screening, composition, divergence of assembled viruses, and the unnamed ' +
    'remainder. Every figure carries its operating point and control. No figure ' +
    'scores novelty.'));
  // Two scopes on one page, said once at the top. The classified profile and the
  // residue search read every pair; the cascade reads a block, because it cannot be run
  // over the library in this container. A reader who assumes one scope for all of it
  // will misread the panels against each other.
  const lib = (D.profile || {}).library_pairs;
  h.append(opChips({sample: acc, library: lib ? lib.toLocaleString() + ' read pairs' : null,
    // The molecule is a fact about the library, and `untargeted shotgun` was a claim
    // about its preparation that no fact here carries: PRJNA948850 alone holds 99
    // hybrid-capture re-runs of the same RNA, so the word is not free.
    'library type': facts.molecule || null,
    'cascade block': ((D.parameters || {}).pairs || 0).toLocaleString() + ' pairs',
    seed: r.seed, 'mate 2 trim': inp.trim_r2,
    'git sha': (r.git_sha || '').slice(0, 10),
    complete: String(!!r.complete)}));
  if (r.incomplete_because && r.incomplete_because.length) {
    const d = el('div', 'note');
    d.innerHTML = '<b>Incomplete run.</b> ' +
      r.incomplete_because.map(esc).join('; ') +
      '. A stage that could not run is absent here rather than estimated.';
    h.append(d);
  }
  h.append(el('p', 'rule', (r.interface_rules || []).join(' / ')));
  app.append(h);
})();

// ---------- hero: the library as 100 squares
// One square is one percent of the read pairs in this run. The point of drawing it this
// way rather than as another stacked bar is that a reader can count the dark squares,
// which is the whole argument of the project in one picture. Every value comes from
// tier_table.json; nothing here is typed in.
(function () {
  const t = D.tier_table, b = (t && t.breakdown) || {};
  if (!b.identified || !b.residual) return;
  const p = panel('Library partition, one square per percent',
    'One square is one percent of this run\\'s read pairs. Pale squares were named by ' +
    'searching; dark squares were not.');
  const frac = {identified: b.identified.point,
    assembled_unidentified: (b.assembled_unidentified || {}).point || 0,
    residual: b.residual.point};
  const COL = {identified: 'var(--identified)', assembled_unidentified: 'var(--tier2)',
    residual: 'var(--residue)'};
  const LAB = {identified: 'identified outright',
    assembled_unidentified: 'assembled, nothing could name it',
    residual: 'no name and no contig'};
  // Largest-remainder apportionment, so the 100 squares sum to 100 and the rounding
  // error does not quietly land on the smallest pile - which is the pile of interest.
  const keys = Object.keys(frac);
  const exact = keys.map(k => 100 * frac[k]);
  const base = exact.map(Math.floor);
  let left = 100 - base.reduce((a, c) => a + c, 0);
  keys.map((k, i) => [i, exact[i] - base[i]]).sort((a, c) => c[1] - a[1])
    .forEach(([i]) => { if (left > 0) { base[i]++; left--; } });
  const squares = [];
  keys.forEach((k, i) => { for (let j = 0; j < base[i]; j++) squares.push(k); });
  const W = 920, gap = 6, cols = 20, rows = 5;
  const side = (W - gap * (cols - 1)) / cols;
  const H = (side + gap) * rows - gap;
  let svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="xMinYMin meet" ' +
    'role="img" aria-label="one hundred squares, ' + base.map((n, i) =>
      n + ' ' + LAB[keys[i]]).join(', ') + '">';
  squares.forEach((k, i) => {
    const x = (i % cols) * (side + gap), y = Math.floor(i / cols) * (side + gap);
    svg += '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' +
      side.toFixed(1) + '" height="' + side.toFixed(1) + '" rx="4" fill="' + COL[k] +
      '" opacity="' + (k === 'identified' ? '0.30' : '0.95') + '"><title>' +
      esc(LAB[k]) + '</title></rect>';
  });
  svg += '</svg>';
  const holder = el('div'); holder.innerHTML = svg; p.append(holder);
  const lg = el('div', 'legend');
  lg.innerHTML = keys.map((k, i) => '<span><i style="background:' + COL[k] +
    ';opacity:' + (k === 'identified' ? '0.30' : '0.95') + '"></i>' +
    base[i] + ' &middot; ' + esc(LAB[k]) + ' (' + pct(frac[k]) + ')</span>').join('');
  p.append(lg);

  const g = el('div', 'hero');
  function card(v, lab) { const d = el('div');
    d.append(el('div', 'bignum', v)); d.append(el('div', 'biglab', lab)); g.append(d); }
  card(pct(frac.assembled_unidentified + frac.residual),
    'of the library that the largest databases available could not name, at the operating ' +
    'point in the chips below');
  const R = D.ranking, hc = (R && R.host_confidence_summary) || {};
  if (hc.subjects != null) card(String(hc.subjects),
    'known organisms the leftovers resemble from a distance, none of them closely');
  /* THE SAME PAIRING OF TWO MAXIMA, IN THE HEADLINE CARD. The warning panel further
     down takes the closest row and prints that row's own alignment length; this card
     printed the maximum identity beside the maximum aligned length, which are maxima
     over the rows taken separately and usually different subjects. It read "96.3%, and
     over only 60 amino acids" where the 96.3% row is a 31 aa alignment. The producer
     ships `closest_row_aligned_aa` for exactly this, and a field saying in words that
     the two maxima must not be printed as one match; the page was ignoring both. The
     length is taken from the rows where they are available, so a ranking written before
     that field renders the same, and the card is dropped rather than printed without a
     length, because an identity over an unstated alignment length is the half of the
     pair that reads as the whole. */
  const closestRow = (D.ranking && D.ranking.subjects || [])
    .filter(r => r.median_identity_percent != null)
    .sort((a, b) => (b.median_identity_percent - a.median_identity_percent) ||
                    ((b.median_aligned_aa || 0) - (a.median_aligned_aa || 0)))[0];
  const closestAa = (closestRow && closestRow.median_aligned_aa) != null
    ? closestRow.median_aligned_aa : hc.closest_row_aligned_aa;
  if (hc.best_identity_percent != null && closestAa != null) {
    card(num(hc.best_identity_percent, 1) + '%',
      'the closest any of those resemblances gets, over ' + closestAa +
      ' amino acids of that one match');
  }
  // The wall, read off the only positive control the tool has. A spiked genome whose
  // reads SURVIVE tier 1 is a genome the cascade failed to recognise, so the first rung
  // where the median genome is more than half unrecognised is where recognition breaks.
  const cal = D.calibration_ladder && (D.calibration_ladder.result || D.calibration_ladder);
  const rungs = (cal && cal.rungs) || [];
  if (rungs.length) {
    const med = r => {
      const v = Object.values(r.genomes || {})
        .map(g => (g.survival || {}).survives_tier1).filter(x => x != null).sort((a, c) => a - c);
      return v.length ? v[Math.floor(v.length / 2)] : null;
    };
    const ord = rungs.slice().sort((a, c) => a.divergence - c.divergence);
    const die = ord.find(r => { const v = med(r); return v != null && v >= 0.5; });
    if (die) card(num(100 * die.divergence, 0) + '%',
      'divergence at which the cascade stops recognising more than half of a known virus ' +
      'planted in this very library');
  }
  p.append(g);
  p.append(el('p', 'note',
    'The dark squares are the subject of this page. Unnameable sequence is the normal ' +
    'state of a wastewater metagenome. What follows measures its composition, its ' +
    'distance from anything known, and where the instrument stops reading.'));
  p.append(opChips((t || {}).operating_point || {}));
  app.append(p);
})();

// ---------- in plain words
// A reader who is not a bioinformatician should be able to leave this page knowing what
// was done and what was found. Every figure in it is pulled from the same files the rest
// of the page reads, so it cannot drift from the panels below it.
(function () {
  const t = D.tier_table, S = D.strains, L = D.linkage, C = D.census;
  const p = panel('Run summary',
    'Inputs, thresholds and outputs of this run.');
  const b = (t && t.breakdown) || {};
  const lede = el('p', 'lede',
    esc(((D.page || {}).sample || {}).lede || '') + ' The cascade removes ' +
    'identifiable sequence and the remainder is characterised. The reported quantity ' +
    'is distance to the nearest known relative, with the confidence in that ' +
    'distance.' + cite('cascade'));
  p.append(lede);
  const g = el('div', 'big');
  function big(v, lab) { const d = el('div');
    d.append(el('div', 'bignum', v)); d.append(el('div', 'biglab', lab)); g.append(d); }
  if (b.identified) big(pct(b.identified.point), 'of the library identified');
  if (b.assembled_unidentified) big(pct(b.assembled_unidentified.point),
    'assembled, unidentified');
  if (b.residual) big(pct(b.residual.point),
    'residual: no name, no contig');
  if (S && S.subjects_total != null) big(String(S.subjects_total),
    'distinct subjects matched by residual reads');
  p.append(g);
  const bits = [];
  if (S && (S.subjects || []).length) {
    const top = (S.subjects || []).slice().sort((a, c) => c.residue_reads - a.residue_reads)[0];
    const nm = (D.names || {})[top.subject];
    bits.push('The leftover reads most often resemble ' +
      '<b>' + esc((nm && (nm.organism || nm.name)) || top.subject) + '</b>, at ' +
      num(top.median_identity_percent, 0) + '% identity over ' + top.median_aligned_aa +
      ' amino acids. That is a relative, not an identification.');
  }
  const dj = L && L.result && L.result.mates_disjoint && L.result.mates_disjoint['10'];
  if (dj) {
    bits.push('Two reads that came off the same fragment of genetic material find each ' +
      'other <b>' + num(dj.enrichment, 0) + ' times</b> more often than chance here, ' +
      'counting only the pairs that do not overlap each other. That is a real signal and ' +
      'still far too weak to gather up the reads of one unknown organism: it finds the ' +
      'partner itself in ' + pct(dj.recovery) + ' of attempts.');
  }
  const cu = C && C.result && C.result.residue && C.result.residue.piles &&
    C.result.residue.piles.unknown_neighbourhood;
  if (cu) {
    bits.push('About <b>' + pct(cu.point) + '</b> of the leftovers sit next to sequence ' +
      'that is itself unidentified, so they are not simply the unlabelled edges of ' +
      'things we did recognise.');
  }
  bits.push('The honest limit: nothing on this page can see an organism present at less ' +
    'than about one read pair in a million, and recognition of a relative fails ' +
    'somewhere between 10% and 20% divergence. Both limits are measured, and both are ' +
    'in the blind spots panel at the foot of the page with the experiment that found them.');
  const ul = el('ul');
  bits.forEach(t2 => { const li = document.createElement('li');
    li.className = 'lede'; li.style.marginBottom = '0.5rem'; li.innerHTML = t2;
    ul.append(li); });
  p.append(ul);
  app.append(p);
})();

// ---------- how the library splits
// The whole instrument as one picture: three bars at three scales, joined by wedges that
// state their own magnification. A single-scale flow would be a lie twice over - the
// residue is a ninth of the library and the named part a seventeenth of that, so both
// would vanish - and the second wedge also changes the UNIT, from read pairs of the whole
// block to the mates the per-read report actually searched. The wedge says so on its face
// rather than in a footnote, because that is the join a reader would otherwise assume.
(function () {
  const t = D.tier_table, b = (t && t.breakdown) || {};
  const R = D.ranking, K = R && R.kind_summary;
  if (!b.identified || !K) return;
  const p = panel('Cascade partition, nested',
    'Every read pair in three steps. Each panel is the highlighted slice of the one ' +
    'before it, at full height.');

  const COL = {identified: 'var(--identified)', assembled_unidentified: 'var(--tier2)',
    residual: 'var(--residue)', named: 'var(--identified)', none: 'var(--null)',
    putative_phage: 'var(--g-phage)', virus_no_phage_signal: 'var(--residue)',
    putative_bacterial: 'var(--g-vert)', putative_eukaryotic: 'var(--g-plant)',
    cellular_kind_not_placed: 'var(--g-unknown)'};

  const W = 920, H = 430, top = 34, bot = 46, bw = 96;
  const X = [26, 392, 800];            // left edge of each bar
  const inner = H - top - bot;

  // Each stage: a list of [key, value, label], the key of the slice the next stage zooms
  // into, and the caption printed on the wedge that leaves it.
  const st1 = [["identified", b.identified.point, "identified outright"],
               ["assembled_unidentified", (b.assembled_unidentified || {}).point || 0,
                "assembled, nothing could name it"],
               ["residual", b.residual.point, "no name, no contig"]];
  const kinds = Object.keys(K.by_kind).filter(k => K.by_kind[k].reads > 0)
    .sort((a, c) => K.by_kind[c].reads - K.by_kind[a].reads);
  const st2 = [["named", K.named, "has a named relative"],
               ["none", K.searched - K.named, "points at nothing we can name"]];
  const st3 = kinds.map(k => [k, K.by_kind[k].reads, K.by_kind[k].label]);

  function layout(items) {
    const tot = items.reduce((a, c) => a + c[1], 0) || 1;
    let y = top;
    return items.map(([k, v, lab]) => {
      const h = inner * v / tot;
      const o = {k: k, v: v, lab: lab, y: y, h: h, frac: v / tot};
      y += h;
      return o;
    });
  }
  const L = [layout(st1), layout(st2), layout(st3)];
  const seg = (i, k) => L[i].filter(o => o.k === k)[0];

  let g = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="how the ' +
    'library splits: ' + st1.map(x => x[2]).join(', ') + ', then ' +
    st2.map(x => x[2]).join(' against ') + ', then ' +
    st3.map(x => x[2]).join(', ') + '">';

  // The wedges go down first so the bars sit on top of them.
  function wedge(i, k, caption, zoom) {
    const a = seg(i, k), x1 = X[i] + bw, x2 = X[i + 1];
    g += '<path d="M' + x1 + ' ' + a.y.toFixed(1) + ' L' + x2 + ' ' + top +
      ' L' + x2 + ' ' + (top + inner) + ' L' + x1 + ' ' + (a.y + a.h).toFixed(1) +
      ' Z" fill="' + (COL[k] || 'var(--null)') + '" opacity="0.10"/>';
    g += '<text x="' + ((x1 + x2) / 2) + '" y="' + (top - 16) + '" text-anchor="middle" ' +
      'font-size="11" fill="var(--dim)">' + esc(caption) + '</text>';
    g += '<text x="' + ((x1 + x2) / 2) + '" y="' + (H - bot + 28) + '" ' +
      'text-anchor="middle" font-size="11" fill="var(--dim)" ' +
      'font-family="var(--mono)">blown up ' + zoom + '</text>';
  }
  wedge(0, 'residual', 'the leftovers, on their own scale',
        '&times;' + (1 / seg(0, 'residual').frac).toFixed(1));
  wedge(1, 'named', 'the few with a relative, on their own scale',
        '&times;' + (1 / seg(1, 'named').frac).toFixed(1));

  const TITLE = ['every read pair in the run', 'the leftover mates that were searched',
                 'the ones with a named relative'];
  const TOT = [t.pairs, K.searched, K.named];
  L.forEach((col, i) => {
    g += '<text x="' + (X[i] + bw / 2) + '" y="' + (top - 16) + '" text-anchor="middle" ' +
      'font-size="11.5" fill="var(--fg)" font-weight="500">' +
      (TOT[i] || 0).toLocaleString() + '</text>';
    g += '<text x="' + (X[i] + bw / 2) + '" y="' + (H - bot + 28) + '" ' +
      'text-anchor="middle" font-size="11" fill="var(--dim)">' + esc(TITLE[i]) + '</text>';
    col.forEach(o => {
      if (o.h <= 0) return;
      g += '<rect class="flowseg" data-k="' + o.k + '" x="' + X[i] + '" y="' +
        o.y.toFixed(1) + '" width="' + bw + '" height="' + Math.max(1, o.h - 1.5)
        .toFixed(1) + '" rx="2" fill="' + (COL[o.k] || 'var(--null)') + '" opacity="' +
        (o.k === 'identified' || o.k === 'none' ? '0.34' : '0.95') + '"><title>' +
        esc(o.lab + ': ' + o.v.toLocaleString() + ' (' + (100 * o.frac).toFixed(2) +
            '% of this panel)') + '</title></rect>';
      // Only label a slice tall enough to carry one, and never over a wedge.
      if (o.h >= 17) {
        // The last column's labels run leftwards: there is no room to its right and a
        // clipped label is worse than none.
        const lastCol = i === L.length - 1;
        g += '<text x="' + (lastCol ? X[i] - 8 : X[i] + bw + 8) + '" y="' +
          (o.y + o.h / 2 + 4).toFixed(1) + '" font-size="11" fill="var(--dim)" ' +
          'text-anchor="' + (lastCol ? 'end' : 'start') + '">' + esc(o.lab) +
          ' <tspan font-family="var(--mono)" fill="var(--fg)">' +
          (100 * o.frac).toFixed(o.frac < 0.1 ? 2 : 1) + '%</tspan></text>';
      }
    });
  });
  g += '</svg>';
  const holder = el('div', 'scroll flowwrap'); holder.innerHTML = g; p.append(holder);
  // Click a slice of the last panel to see the records behind it. The diagram and the
  // list are then showing one pile, which is the question a reader has at exactly this
  // point and otherwise has to carry down the page by hand.
  Array.prototype.forEach.call(holder.querySelectorAll('.flowseg'), n => {
    const k = n.getAttribute('data-k');
    if (!K.by_kind[k]) return;
    n.style.cursor = 'pointer';
    n.addEventListener('click', () => {
      if (window.__hapaxShowKind) window.__hapaxShowKind(k);
    });
  });
  p.append(el('p', 'note',
    '<b>The second join changes unit and sample.</b> Panels one and two are read ' +
    'PAIRS over the whole block. Panels three and four are MATES, and only the ' +
    K.searched.toLocaleString() + ' the per-read report searched, so the wedge between ' +
    'them is a re-basing. No percentage runs straight across this figure.' +
    cite('library')));
  app.append(p);
})();

// ---------- each index takes a bite
// The cascade as a column of splits rather than a table. One row per index, in the order
// they actually ran: pale is what earlier indexes had already removed, saturated is this
// index's OWN marginal bite, and the rest is what was still left after it. Marginal and
// not standalone, because the standalone figures overlap and a reader who adds them up
// gets more than the library.
(function () {
  const t = D.tier_table, curve = (t && t.curve) || [];
  if (!curve.length) return;
  const p = panel('Marginal removal per index',
    'One row per index, in run order. The dark slice is what that index removed and ' +
    'nothing before it had.');
  const rows = el('div', 'bites');
  curve.forEach((c, i) => {
    const before = c.cumulative_removed - c.marginal_removed;
    const r = el('div', 'bite');
    r.innerHTML =
      '<div class="bname mono">' + esc(c.tier) + '</div>' +
      '<div class="bbar">' +
        '<i class="before" style="width:' + (100 * before).toFixed(3) + '%"></i>' +
        '<i class="mine" style="width:' + (100 * c.marginal_removed).toFixed(3) + '%">' +
        '</i></div>' +
      '<div class="bnum mono">' +
        (c.marginal_removed >= 0.0005
          ? '+' + (100 * c.marginal_removed).toFixed(2) + ' pts'
          : '<span class="dim">nothing new</span>') +
      '</div>' +
      '<div class="bsize mono dim">' + gb(c.index_bytes) + '</div>';
    r.title = c.tier + ': removed ' + (100 * c.tier_alone_removed).toFixed(2) +
      '% on its own, ' + (100 * c.marginal_removed).toFixed(2) +
      ' points that nothing before it had taken, leaving ' +
      (100 * c.residue).toFixed(2) + '% of the library';
    rows.append(r);
  });
  p.append(rows);
  const last = curve[curve.length - 1];
  /* THE CAPTION NAMES THE ROWS THE TABLE DREW, AND IT NAMED A DIFFERENT ONE. It read
     "the assembly of this sample is the biggest single bite in the cascade, at 2.3 points
     from 7.1 MB, while the largest index here is tens of gigabytes and buys fractions of
     a point". Three things were wrong on this run and none of them is visible beside a
     table of ten rows: the biggest bite is the ribosomal index at 75.82 points, the
     self-assembly row the sentence picked is the smaller of the two self-assembly rows
     and left the 9.73-point one out, and the largest index bought 1.44 points rather than
     a fraction of one. The point the sentence exists to make - that the index which buys
     the most per byte is the sample's own assembly - is true and is now computed, with
     every number read off the same rows the bars are drawn from. */
  /* A RATIO NEEDS A DENOMINATOR THAT IS THE SAME KIND OF THING ON EVERY ROW. Per byte
     over every row makes the artefact filter the winner: it is 120 bytes of patterns,
     not a reference, so the ratio is a small number over a smaller one and says nothing
     about consulting an index. The comparison is between references, and the floor that
     separates a reference from a rule is stated in the sentence rather than assumed. */
  const REFERENCE_MIN_BYTES = 1e6;
  const sized = curve.filter(c => (c.index_bytes || 0) >= REFERENCE_MIN_BYTES);
  const byMarginal = curve.slice()
    .sort((a, b) => (b.marginal_removed || 0) - (a.marginal_removed || 0));
  const biggest = byMarginal[0] || {};
  const byBytes = sized.slice().sort((a, b) => b.index_bytes - a.index_bytes);
  const largest = byBytes[0] || {};
  const perByte = sized.slice().sort((a, b) =>
    ((b.marginal_removed || 0) / b.index_bytes) -
    ((a.marginal_removed || 0) / a.index_bytes));
  const bestValue = perByte[0] || {};
  const nothing = curve.filter(c => (c.marginal_removed || 0) < 0.0005).length;
  p.append(el('p', 'note',
    '<b>Which index you consult beats how big it is.</b> ' +
    'Of the references over ' + gb(REFERENCE_MIN_BYTES) + ', the most per byte is ' +
    esc(bestValue.tier) + ', ' +
    num(100 * (bestValue.marginal_removed || 0), 2) + ' points from ' +
    gb(bestValue.index_bytes) + ', while the largest index here, ' +
    esc(largest.tier) + ' at ' + gb(largest.index_bytes) + ', buys ' +
    num(100 * (largest.marginal_removed || 0), 2) + ' points. ' +
    'The biggest single bite is ' + esc(biggest.tier) + ' at ' +
    num(100 * (biggest.marginal_removed || 0), 2) + ' points. ' +
    (nothing
      ? '<b>' + nothing + ' row' + (nothing === 1 ? '' : 's') + ' take' +
        (nothing === 1 ? 's' : '') + ' nothing at all</b>: by the time ' +
        (nothing === 1 ? 'it runs, everything it could' : 'they run, everything they ' +
         'could') + ' name has already been named by something else. '
      : '') +
    'Everything left after the last row is ' + pct(last.residue) +
    ' of the library, which is the subject of this page.'));
  p.append(opChips((t || {}).operating_point || {}));
  app.append(p);
})();

// ---------- phage, virus, or something living
// The split Imo asked for, in READS and with its control beside it. Subject counts and
// read counts answer different questions - 204 of 243 subjects being phage is not the
// claim that most of the residue's named reads are phage - so the picture is drawn in
// reads and the subject count rides along in the tooltip.
(function () {
  const R = D.ranking, K = R && R.kind_summary;
  const p = panel('Kind assignment of matched subjects',
    'Leftover reads with a named relative, by what kind of thing that relative is.');
  if (!K) { p.append(el('div', 'absent', 'no kind summary in ranking.json'));
    app.append(p); return; }
  /* NOTHING SEARCHED IS A CENSUS, AND THE INTERVAL IS CORRECTLY NULL. The producer
     writes `share_named: null` where the denominator is zero, which is right: a share of
     nothing has no value and inventing 0% would read as "no leftover has a relative"
     rather than "no leftover was searched". This panel then read `.point` off that null,
     threw, and took every panel after it off the page - the whole report, from one
     legitimate zero. So the zero is rendered, and it says which zero it is. */
  if (!K.searched) {
    p.append(el('div', 'absent', 'no leftover mate in this run was searched for a ' +
      'named relative, so there is no kind to assign. That is a count of this run and ' +
      'not a finding: zero searched is not zero named.'));
    app.append(p); return; }
  const KCOL = {putative_phage: 'var(--g-phage)',
    virus_no_phage_signal: 'var(--residue)',
    putative_bacterial: 'var(--g-vert)',
    putative_eukaryotic: 'var(--g-plant)',
    cellular_kind_not_placed: 'var(--g-unknown)'};

  // Two bars at ONE scale: the residue, and the already-named pile it is read against.
  // Drawn first and together, because 5.8% means nothing until you see 23.8% beside it.
  function scaleBar(label, named, total, sub) {
    const row = el('div');
    row.append(el('div', 'barlab', esc(label)));
    const bar = el('div', 'bar');
    const a = document.createElement('span');
    a.style.width = (100 * named / total) + '%';
    a.style.background = 'var(--identified)';
    a.title = named + ' of ' + total;
    const b = document.createElement('span');
    b.style.flex = '1';
    b.style.background = 'var(--null)';
    b.style.opacity = '0.28';
    bar.append(a, b);
    row.append(bar);
    row.append(el('div', 'note', sub));
    return row;
  }
  const g = el('div', 'stack');
  g.append(scaleBar('the leftovers', K.named, K.searched,
    '<b>' + pct(K.share_named.point) + '</b> of ' + K.searched.toLocaleString() +
    ' searched mates have a significant relative. ' +
    (K.searched - K.named).toLocaleString() + ' have none.'));
  const C = K.known_pile_control;
  if (C && C.searched) g.append(scaleBar('reads the cascade already named, the control',
    C.significant, C.searched,
    '<b>' + pct(C.share.point) + '</b> of ' + C.searched.toLocaleString() +
    ' searched mates. These have a known relative by construction, and even they do ' +
    'not all clear the same threshold.'));
  p.append(g);

  // Then the sliver, magnified. Labelled as magnified, on its own scale, because a
  // breakdown drawn at the scale above would be four invisible slices.
  const kinds = Object.keys(K.by_kind).filter(k => K.by_kind[k].reads > 0)
    .sort((a, b) => K.by_kind[b].reads - K.by_kind[a].reads);
  p.append(el('div', 'barlab',
    'those ' + K.named + ' reads, magnified: what kind of record they point at'));
  const bar = el('div', 'bar');
  kinds.forEach(k => {
    const v = K.by_kind[k];
    const sp = document.createElement('span');
    sp.style.width = (100 * v.reads / K.named) + '%';
    sp.style.background = KCOL[k] || 'var(--null)';
    sp.title = v.label + ': ' + v.reads + ' reads across ' + v.subjects + ' records';
    bar.append(sp);
  });
  p.append(bar);
  const lg = el('div', 'legend');
  lg.innerHTML = kinds.map(k => {
    const v = K.by_kind[k];
    return '<span data-k="' + k + '" class="clicky"><i style="background:' +
      (KCOL[k] || 'var(--null)') + '"></i>' +
      esc(v.label) + ' <b>' + v.reads + '</b> reads, ' + v.subjects + ' records ' +
      '<span class="mono">(' + pct(v.share_of_searched.point) + ' of all leftovers)' +
      '</span></span>';
  }).join('');
  Array.prototype.forEach.call(lg.children, n => {
    const k = n.getAttribute('data-k');
    if (!k) return;
    n.addEventListener('click', () => {
      if (window.__hapaxShowKind) window.__hapaxShowKind(k);
    });
  });
  p.append(lg);

  const big = el('div', 'hero');
  function card(v, lab, col) { const d = el('div');
    const n = el('div', 'bignum', v); if (col) n.style.color = col;
    d.append(n); d.append(el('div', 'biglab', lab)); big.append(d); }
  card(pct(1 - K.share_named.point), 'of the leftovers point at nothing we can name',
    'var(--null)');
  kinds.slice(0, 3).forEach(k => card(pct(K.by_kind[k].share_of_searched.point),
    K.by_kind[k].label, KCOL[k]));
  p.append(big);

  const det = el('details');
  det.innerHTML = '<summary>what these three buckets can and cannot tell you</summary>';
  const dd = el('div', 'note');
  dd.innerHTML = '<p><b>' + esc(K.caveat) + '</b></p>' +
    kinds.map(k => '<p><b>' + esc(K.by_kind[k].label) + '</b> &mdash; ' +
      esc(K.by_kind[k].reading) + '</p>').join('') +
    '<p>' + esc(K.unit) + '. Denominator: ' + esc(R.denominator_note || '') + '.</p>';
  det.append(dd); p.append(det);
  app.append(p);
})();

// ---------- worth a second look
(function () {
  const R = D.ranking;
  const p = panel('Subject ranking, four axes',
    'Named relatives of the leftover reads, ordered on four measured axes. All four ' +
    'are shown on every row.');
  if (!R) { p.append(el('div', 'absent',
    'ranking.json absent: run experiments/concern_ranking.py over this directory'));
    app.append(p); return; }
  p.append(scopeLine(R));

  // The honest headline of the whole panel, and it is a negative. Put it above the
  // leaderboard rather than under it: a reader who scrolls the cards without seeing this
  // will read the host colours as claims about the sample.
  const hc = R.host_confidence_summary || {};
  // ONE ROW'S IDENTITY WITH THAT ROW'S OWN ALIGNMENT LENGTH. The summary carries the
  // maximum identity and the maximum aligned length as two independent maxima, and this
  // sentence used to put them together: on the second production library that read
  // "100.0% identity over 53 amino acids" when the 100% was a 21 aa alignment and the
  // 53 aa belonged to a different subject. A pairing of two maxima is a subject that
  // does not exist, and identity over a short alignment is a ratio of small integers,
  // which is exactly why the length has to be the one that identity was measured over.
  // The producer now carries that length itself (closest_row_aligned_aa, with a field
  // saying the two maxima are separate maxima), and the two agree on both measured
  // libraries - 21 aa here, 35 on the other - which is the cross-check worth having.
  // The page keeps deriving it from the rows so a ranking written before that field
  // renders the same rather than silently losing the length.
  const best = (R.subjects || []).filter(r => r.median_identity_percent != null)
    .sort((a, b) => (b.median_identity_percent - a.median_identity_percent) ||
                    ((b.median_aligned_aa || 0) - (a.median_aligned_aa || 0)))[0];
  const warn = el('div', 'callout');
  warn.innerHTML = '<b>Read this before the colours below.</b> No residue subject is ' +
    'close enough for its host to transfer. ' +
    (best ? 'The closest of all ' + hc.subjects + ' reads ' +
      num(best.median_identity_percent, 1) + '% identity over ' +
      best.median_aligned_aa + ' amino acids, '
      : 'Not one of the ' + hc.subjects + ' carries an identity, ') +
    '<b>' + hc['at_or_above_0.5'] + ' of ' + hc.subjects +
    '</b> clear half confidence, and the median is ' + num(hc.median, 2) +
    '. The host groups below label a database record. The surveillance question has ' +
    'run past the edge of what the references answer, which is itself the ' +
    'finding.' + cite('residue');
  p.append(warn);

  const GCOL = {bacteriophage: 'var(--g-phage)', plant: 'var(--g-plant)',
    fungal_or_protist: 'var(--g-fungal)', invertebrate: 'var(--g-invert)',
    vertebrate_non_human: 'var(--g-vert)', human_associated: 'var(--g-human)',
    unknown_host: 'var(--g-unknown)'};
  // The table's own display label, never the key with its underscores stripped:
  // "vertebrate non human" read as "human" at a glance on a narrow card, which is the
  // one misreading this panel cannot afford.
  const GLAB = g => ((R.host_table.groups[g] || {}).label) || g.replace(/_/g, ' ');
  const KCOL2 = {putative_phage: 'var(--g-phage)',
    virus_no_phage_signal: 'var(--residue)', putative_bacterial: 'var(--g-vert)',
    putative_eukaryotic: 'var(--g-plant)',
    cellular_kind_not_placed: 'var(--g-unknown)'};

  const ctl = el('div', 'controls');
  const sel = el('select');
  [['second_look', 'overall, worth a second look'],
   ['distance', 'furthest from anything known'],
   ['residue_reads', 'most reads pointing at it'],
   ['host', 'host group of the relative']].forEach(([v, t]) => {
    const o = document.createElement('option'); o.value = v; o.textContent = t;
    sel.append(o); });
  const gsel = el('select');
  const gnames = ['every group'].concat(Object.keys(R.host_table.groups));
  gnames.forEach(g => { const o = document.createElement('option'); o.value = g;
    o.textContent = g === 'every group' ? g : GLAB(g); gsel.append(o); });
  const ksel = el('select');
  const kinds = (R.kind_summary && R.kind_summary.by_kind) || {};
  [['every kind', 'phage, virus and cellular together']].concat(
    Object.keys(kinds).filter(k => kinds[k].subjects > 0)
      .map(k => [k, kinds[k].label])).forEach(([v, t]) => {
    const o = document.createElement('option'); o.value = v; o.textContent = t;
    ksel.append(o); });
  ctl.append(el('span', 'dim', 'order by'), sel, el('span', 'dim', 'show'), ksel, gsel);
  p.append(ctl);

  // The group mix as one bar, because 201 of 243 being phage is the single most
  // informative picture on this page and a table buries it.
  // Over EVERY subject in the ranking, not only the rows embedded below: the bar is the
  // shape of the whole thing and the cards are a window onto the top of it.
  const mix = R.host_mix || {};
  const mixn = R.subjects_total || (R.subjects || []).length;
  const mixbar = el('div', 'bar');
  Object.keys(R.host_table.groups).forEach(g => {
    if (!mix[g]) return;
    const sp = document.createElement('span');
    sp.style.width = (100 * mix[g] / mixn) + '%';
    sp.style.background = GCOL[g] || 'var(--null)';
    sp.title = g + ': ' + mix[g];
    mixbar.append(sp);
  });
  p.append(mixbar);
  const mlg = el('div', 'legend');
  mlg.innerHTML = Object.keys(R.host_table.groups).filter(g => mix[g]).map(g =>
    '<span><i style="background:' + (GCOL[g] || 'var(--null)') + '"></i>' +
    esc(GLAB(g)) + ' ' + mix[g] + '</span>').join('');
  p.append(mlg);

  if (R.subjects_embedded && R.subjects_total &&
      R.subjects_embedded < R.subjects_total) {
    p.append(el('p', 'note',
      'The bar above and the map at the end of this page are over all ' +
      R.subjects_total.toLocaleString() + ' subjects. The cards below are the ' +
      R.subjects_embedded.toLocaleString() + ' highest ranked, including the highest ' +
      'ranked of every group and every kind, so each filter has real rows behind it. ' +
      'The whole table is in ranking.json in the run directory.'));
  }

  const grid = el('div', 'cards');
  function draw() {
    const v = sel.value, g = gsel.value, k = ksel.value;
    const list = (R.subjects || [])
      .filter(r => k === 'every kind' || r.sequence_kind === k)
      .filter(r => g === 'every group' || r.host_group === g)
      .slice().sort((a, b) =>
        v === 'distance' ? b.distance - a.distance :
        v === 'residue_reads' ? b.residue_reads - a.residue_reads :
        v === 'host' ? (b.host_term - a.host_term) || (b.second_look - a.second_look) :
        b.second_look - a.second_look).slice(0, 48);
    grid.innerHTML = '';
    list.forEach(r => {
      const c = el('div', 'card');
      c.style.borderLeft = '3px solid ' + (GCOL[r.host_group] || 'var(--null)');
      const axis = (lab, val, col, extra) =>
        '<div class="axrow"><span class="axlab">' + lab + '</span>' +
        '<span class="axbar"><i style="width:' + (100 * Math.max(0, Math.min(1, val)))
        .toFixed(1) + '%;background:' + col + '"></i></span>' +
        '<span class="axval">' + extra + '</span></div>';
      c.innerHTML =
        '<div class="rankno">#' + r.rank + '</div>' +
        '<div class="ctitle">' + esc(r.organism || r.subject) + '</div>' +
        (r.protein ? '<div class="csub">' + esc(r.protein) + '</div>' : '') +
        '<div class="tag" style="border-color:' + (KCOL2[r.sequence_kind] ||
            'var(--edge)') + ';color:' + (KCOL2[r.sequence_kind] || 'var(--dim)') + '">' +
          esc(r.sequence_kind_label || '') +
          (r.sequence_kind_matched_keyword ? ' &middot; "' +
            esc(r.sequence_kind_matched_keyword) + '"' : '') +
        '</div> ' +
        '<div class="tag" style="border-color:' + (GCOL[r.host_group] || 'var(--edge)') +
          ';color:' + (GCOL[r.host_group] || 'var(--dim)') + '">' +
          (r.host_group_not_applicable
            ? 'host: not applicable, it is the organism'
            : esc(r.host_group_label || GLAB(r.host_group)) +
              (r.host_group_matched_keyword ? ' &middot; "' +
                esc(r.host_group_matched_keyword) + '"' : ' &middot; no rule matched')) +
        '</div>' +
        axis('how far', r.distance, 'var(--residue)',
             num(r.median_identity_percent, 1) + '% id') +
        axis('how much', r.weight_term, 'var(--identified)', r.residue_reads + ' reads') +
        axis('host', r.host_term, GCOL[r.host_group] || 'var(--null)',
             'conf ' + num(r.host_confidence, 2)) +
        axis('alone', r.only_in_residue ? 1 : 0, 'var(--tier2)',
             r.only_in_residue ? 'residue only' : 'also identified') +
        '<div class="crow dim mono">E ' + (+r.median_evalue).toExponential(1) +
        ' &middot; ' + r.median_aligned_aa + ' aa &middot; ' + esc(r.subject) + '</div>';
      grid.append(c);
    });
    if (!list.length) grid.innerHTML =
      '<div class="absent">no subject in that combination</div>';
  }
  sel.onchange = draw; gsel.onchange = draw; ksel.onchange = draw; draw();
  // Published so the split diagrams above can drive this list. It only ever sets a
  // filter this panel already offers, so a click cannot show anything the select
  // could not; if the kind is not one of the options the call does nothing.
  window.__hapaxShowKind = function (kind) {
    const ok = Array.prototype.some.call(ksel.options, o => o.value === kind);
    if (!ok) return false;
    ksel.value = kind; draw();
    p.scrollIntoView({behavior: 'smooth', block: 'start'});
    return true;
  };
  p.append(grid);

  const det = el('details');
  det.innerHTML = '<summary>the formula, in full, and what it is blind to</summary>';
  const dd = el('div', 'note');
  dd.innerHTML = '<p class="mono">' + esc(R.formula) + '</p>' +
    Object.entries(R.formula_terms).map(([k, v]) =>
      '<p><b>' + esc(k) + '</b> &mdash; ' + esc(v) + '</p>').join('') +
    '<p><b>' + esc(R.not_a_risk_assessment) + '</b></p>' +
    '<p>' + esc(R.host_table.subjects_no_rule_placed) + ' of ' +
      esc(R.host_table.subjects_total) + ' names no rule in <span class="mono">' +
      esc(R.host_table.source) + '</span> could place. ' +
      esc(R.host_table.coverage_note) + '.</p>' +
    '<p>' + esc(R.interval_note) + '.</p>' +
    '<ul>' + (R.blind_to || []).map(b => '<li>' + esc(b) + '</li>').join('') + '</ul>';
  det.append(dd); p.append(det);
  app.append(p);
})();

// ---------- where the relatives sit
(function () {
  const R = D.ranking;
  if (!R || !(R.subjects || []).length) return;
  const pn = R.points_note;
  const p = panel('Identity against read support, per subject',
    'Distance from each named relative against the number of reads pointing at it.' + (pn ? ' Drawn from a seeded sample of ' + pn.drawn.toLocaleString() +
    ' of the ' + pn.of.toLocaleString() + ' subjects, because a cut by rank would empty ' +
    'the lower left corner and that is where most of the residue is.' : ''));
  const GCOL = {bacteriophage: 'var(--g-phage)', plant: 'var(--g-plant)',
    fungal_or_protist: 'var(--g-fungal)', invertebrate: 'var(--g-invert)',
    vertebrate_non_human: 'var(--g-vert)', human_associated: 'var(--g-human)',
    unknown_host: 'var(--g-unknown)'};
  const W = 920, H = 420, m = {l: 62, r: 20, t: 18, b: 52};
  // The sampled points when the ranking is bigger than the page can carry, mapped onto
  // the same field names the full rows use so the drawing code is unchanged.
  const rows = (R.points || []).length
    ? R.points.map(q => ({median_identity_percent: q.x, residue_reads: q.y,
                          host_group: q.g, weight_term: q.w, median_aligned_aa: q.a,
                          organism: q.n}))
    : R.subjects.filter(r => r.median_identity_percent != null);
  // The axis is fixed from 40 to 100 rather than fitted to the data, because the empty
  // right-hand third IS the result: a fitted axis would silently rescale the gap away
  // and make a cloud of 50% matches look like it reaches identification.
  const obslo = Math.min(...rows.map(r => r.median_identity_percent));
  const obshi = Math.max(...rows.map(r => r.median_identity_percent));
  const xmin = 40, xmax = 100;
  const ymax = Math.max(...rows.map(r => r.residue_reads));
  const px = v => m.l + (v - xmin) / ((xmax - xmin) || 1) * (W - m.l - m.r);
  const py = v => H - m.b - (Math.log1p(v) / Math.log1p(ymax)) * (H - m.t - m.b);
  let s = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="relatives">';
  // Where an identification would sit, drawn as a region rather than asserted in prose.
  s += '<rect x="' + px(90).toFixed(1) + '" y="' + m.t + '" width="' +
    (px(100) - px(90)).toFixed(1) + '" height="' + (H - m.t - m.b) +
    '" fill="var(--identified)" opacity="0.08"/>';
  s += '<text x="' + ((px(90) + px(100)) / 2).toFixed(1) + '" y="' + (m.t + 16) +
    '" font-size="11" fill="var(--dim)" text-anchor="middle">an identification ' +
    'would be here</text>';
  s += '<rect x="' + m.l + '" y="' + m.t + '" width="' + (W - m.l - m.r) + '" height="' +
    (H - m.t - m.b) + '" fill="none" stroke="var(--edge)"/>';
  [1, 2, 5, 10, 25].filter(v => v <= ymax).forEach(v => {
    s += '<line x1="' + m.l + '" x2="' + (W - m.r) + '" y1="' + py(v).toFixed(1) +
      '" y2="' + py(v).toFixed(1) + '" stroke="var(--edge)" stroke-dasharray="2 3"/>' +
      '<text x="' + (m.l - 8) + '" y="' + (py(v) + 4).toFixed(1) + '" font-size="11" ' +
      'fill="var(--dim)" text-anchor="end" font-family="var(--mono)">' + v + '</text>';
  });
  rows.forEach(r => {
    s += '<circle cx="' + px(r.median_identity_percent).toFixed(1) + '" cy="' +
      py(r.residue_reads).toFixed(1) + '" r="' +
      (3 + 2.5 * r.weight_term).toFixed(1) + '" fill="' +
      (GCOL[r.host_group] || 'var(--null)') + '" opacity="0.72"><title>' +
      esc((r.organism || r.subject) + ' - ' + r.residue_reads + ' reads at ' +
          r.median_identity_percent + '%') + '</title></circle>';
  });
  [40, 50, 60, 70, 80, 90, 100].forEach(v => {
    s += '<text x="' + px(v).toFixed(1) + '" y="' + (H - m.b + 17) +
      '" font-size="11" fill="var(--dim)" font-family="var(--mono)" ' +
      'text-anchor="middle">' + v + '%</text>';
  });
  s += '<text x="' + ((W + m.l - m.r) / 2) + '" y="' + (H - 12) + '" text-anchor="middle"' +
    ' font-size="11.5" fill="var(--dim)">identity to the nearest named relative</text>';
  s += '<text transform="translate(16 ' + ((H - m.b + m.t) / 2) + ') rotate(-90)" ' +
    'text-anchor="middle" font-size="11.5" fill="var(--dim)">residue reads pointing at ' +
    'it (log)</text></svg>';
  const holder = el('div', 'scroll'); holder.innerHTML = s; p.append(holder);
  p.append(el('p', 'note',
    (function () {
      // The aligned length is the whole reason the shaded band is not an identification,
      // so it is computed from the drawn points rather than described. This sentence
      // used to list the band's lengths one by one, which was right when the band held
      // two points and unreadable once the search covered the whole residue.
      const mid = a => { const v = a.filter(x => x != null).sort((x, y) => x - y);
        return v.length ? v[Math.floor(v.length / 2)] : null; };
      const inb = rows.filter(r => r.median_identity_percent >= 90);
      let t = 'The cloud sits between ' + num(obslo, 1) + '% and ' + num(obshi, 1) +
        '% identity, over a median of ' + mid(rows.map(r => r.median_aligned_aa)) +
        ' amino acids. ';
      t += inb.length === 0 ? 'The shaded band is empty. '
        : inb.length.toLocaleString() + ' of the ' + rows.length.toLocaleString() +
          ' drawn here reach the shaded band, at a median alignment of ' +
          mid(inb.map(r => r.median_aligned_aa)) + ' amino acids. Over an alignment ' +
          'that short, high identity is ordinary rather than an identification, which ' +
          'is why every hit on this page is filtered to E &lt; ' + D.significant_evalue +
          ' and still reported with its aligned length. ';
      return t; })() +
    'Colour is ' +
    'the matched record\\'s host group and carries the same warning as the panel above: ' +
    'at these identities it describes the database entry, not the sample. ' +
    (function () { const n = rows.filter(r => r.host_group === 'human_associated').length;
      return n === 0
        ? 'Note what is not here: not one of these ' + rows.length + ' records falls in ' +
          'the human-associated group at all.'
        : n + ' of these ' + rows.length + ' records are in the human-associated group, ' +
          'at host confidences the panel above says are too low to transfer.'; })()));
  app.append(p);
})();

// ---------- the breakdown
(function () {
  const t = D.tier_table;
  const p = panel('Class fractions with intervals',
    'Identified, assembled without identification, and left over. Each segment ' +
    'carries its interval.');
  if (!t) { p.append(el('div', 'absent', 'tier_table.json absent: the cascade did not run'));
    app.append(p); return; }
  const last = (t.curve || [])[(t.curve || []).length - 1];
  const b = t.breakdown || {};
  const ident = b.identified || (last && last.cumulative_removed_interval);
  const tier2 = b.assembled_unidentified;
  const resid = b.residual || (last && last.residue_interval);
  const bar = el('div', 'bar');
  const seg = (iv, colour) => { const s = document.createElement('span');
    s.style.width = (100 * ((iv && iv.point) || 0)) + '%'; s.style.background = colour;
    return s; };
  bar.append(seg(ident, 'var(--identified)'));
  if (tier2) bar.append(seg(tier2, 'var(--tier2)'));
  bar.append(seg(resid, 'var(--residue)'));
  p.append(bar);
  const lg = el('div', 'legend');
  lg.innerHTML =
    '<span><i style="background:var(--identified)"></i>identified ' + interval(ident) +
      '</span>' +
    (tier2 ? '<span><i style="background:var(--tier2)"></i>assembled, not identified ' +
      interval(tier2) + '</span>' : '') +
    '<span><i style="background:var(--residue)"></i>residual reads ' + interval(resid) +
      '</span>';
  p.append(lg);
  p.append(opChips(t.operating_point));
  p.append(el('p', 'note', tier2 ? 'The middle segment is observed and never subtracted: ' +
    'these reads map to contigs that assembled without being identified, from ' +
    esc(b.tier2_from || '') + '.' :
    'Tier 2 has no population in this run, so the bar has two segments rather than three.'));
  // The reads scored are a contiguous block of the FASTQ, not a random sample of the
  // library, so the Wilson interval above is narrower than the sampling supports. Both
  // are shown: hiding either would be a choice presented as a measurement.
  const de = b.design_effect;
  if (de && ident && ident.lo_block_corrected != null) {
    const w = iv => '<span class="mono">' + pct(iv.lo_block_corrected) + ', ' +
      pct(iv.hi_block_corrected) + '</span>';
    const d2 = el('p', 'note');
    d2.innerHTML = '<b>Widened for block sampling.</b> These reads are a contiguous ' +
      'block of the file rather than a random sample of the library, and a Wilson ' +
      'interval assumes the latter. Measured over all 8.7 million pairs, the real ' +
      'between-block spread exceeds the binomial one by ' + num(de.identified, 2) +
      'x on identified and ' + num(de.residual, 2) + 'x on the residue, with random ' +
      'samples of the same size reading 1.0 as the control. So the honest intervals are ' +
      'identified [' + w(ident) + ']' + (tier2 ? ', assembled-not-identified [' +
      w(tier2) + ']' : '') + ', residual [' + w(resid) + ']. ' + esc(de._source) + '.';
    p.append(d2);
  }
  app.append(p);
})();

// ---------- the cascade, with the sweep live
(function () {
  const t = D.tier_table;
  const p = panel('Tier removal and threshold sweep',
    'Per index: removal alone, marginal removal, and residue. The slider sweeps the ' +
    'threshold of that tier, which moves a tier further than database size does.');
  if (!t) { p.append(el('div', 'absent', 'tier_table.json absent')); app.append(p); return; }
  const swept = (t.curve || []).filter(c => c.threshold_sweep);
  /* THE CONTROL THAT MAKES THE CASE FOR SNAPPING. Kraken2 confidence was run at 0, 0.05,
     0.1 and 0.2 and nowhere between, and one database's removal of this library moves
     from 42.6% to 3.5% across that gap. So the threshold snaps to the points the run
     measured, the point the run shipped is marked, and nothing is interpolated: a smooth
     slider here would invent a factor of twelve. */
  const SWIDS = ['sweep_tier', 'sweep_point', 'union_tiers'];
  const out = el('div', 'mono');
  if (swept.length) { p.append(ctlRow(SWIDS)); p.append(out); }
  const draw = () => {
    const tier = ctlGet('sweep_tier');
    const g = ((CTLBY.sweep_tier || {}).grid || {})[tier];
    const c = swept.filter(x => x.tier === tier)[0];
    if (!c || !g) { out.textContent = ''; return; }
    const k = String(ctlGet('sweep_point'));
    const pt = (g.points || []).filter(x => String(x.value) === k)[0];
    out.innerHTML = esc(g.quantity + ' ' + k) + ' removes ' + pct(pt ? pt.removed : null) +
      ' of the library on its own' +
      (ctlMoved('sweep_point')
        ? ' <span class="relabel">at your setting; this run shipped ' +
          esc(g.quantity + ' ' + g.shipped) + ', which removes ' +
          pct((g.points.filter(x => String(x.value) === String(g.shipped))[0] || {})
              .removed) + ' and is the number every total on this page is built from' +
          '</span>'
        : ' <span class="dim">at the setting this run shipped</span>') +
      '  (span ' + num(g.span_points, 2) + ' points across the measured points, and the ' +
      'marginal contribution of this tier to the shipped union is ' + pct(g.marginal) +
      ', which is not recomputed here: marginal contribution depends on the order the ' +
      'tiers were applied in)';
  };
  if (swept.length) { ctlOn(draw); draw(); p.append(settingLine(SWIDS)); }

  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>tier</th><th>index</th><th>alone</th><th>marginal</th>' +
    '<th>cumulative</th><th>residue</th><th>sweep span</th></tr></thead>';
  const body = el('tbody');
  (t.curve || []).forEach(c => {
    const tr = el('tr');
    const span = c.sweep_span_points ? Object.values(c.sweep_span_points)[0] : null;
    tr.innerHTML = '<td>' + esc(c.tier) + '</td><td class="mono">' + gb(c.index_bytes) +
      '</td><td class="mono">' + pct(c.tier_alone_removed) + '</td><td class="mono">' +
      pct(c.marginal_removed) + '</td><td>' + interval(c.cumulative_removed_interval) +
      '</td><td>' + interval(c.residue_interval) + '</td><td class="mono">' +
      (span == null ? '<span class="dim">-</span>' : num(span, 2) + ' pts') + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  p.append(el('p', 'note', 'Intervals are Wilson on this run\\'s ' + (t.pairs || 0) +
    ' pairs. A multi-seed run reports across-seed spread instead, which is the stronger ' +
    'statement and the one the headline numbers use.'));
  app.append(p);
})();

// ---------- cost curve
(function () {
  const t = D.tier_table;
  const p = panel('Removal against index size',
    'Removal against bytes consulted, measured points only. The far end is not ' +
    'extrapolated: where the curve stops is where the measurement stops.');
  if (!t) { p.append(el('div', 'absent', 'tier_table.json absent')); app.append(p); return; }
  const pts = (t.curve || []).map(c => ({x: c.cumulative_index_bytes || 1,
    y: c.cumulative_removed, tier: c.tier})).filter(d => d.y != null);
  const W = 920, H = 300, m = {l: 54, r: 14, t: 14, b: 54};
  const xs = pts.map(d => Math.log10(Math.max(d.x, 1)));
  const x0 = Math.min(...xs), x1 = Math.max(...xs, x0 + 1);
  const sx = v => m.l + (Math.log10(Math.max(v, 1)) - x0) / (x1 - x0) * (W - m.l - m.r);
  const sy = v => H - m.b - (v || 0) * (H - m.t - m.b);
  let svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="removal against index bytes">';
  for (let g = 0; g <= 1.0001; g += 0.25) {
    svg += '<line x1="' + m.l + '" x2="' + (W - m.r) + '" y1="' + sy(g) + '" y2="' + sy(g) +
      '" stroke="var(--edge)" fill="none"/><text x="' + (m.l - 8) + '" y="' + (sy(g) + 4) +
      '" text-anchor="end" font-size="11" fill="var(--dim)" font-family="var(--mono)">' +
      (100 * g).toFixed(0) + '%</text>';
  }
  svg += '<path d="' + pts.map((d, i) => (i ? 'L' : 'M') + sx(d.x) + ' ' + sy(d.y)).join(' ') +
    '" fill="none" stroke="var(--identified)" stroke-width="2"/>';
  pts.forEach(d => {
    svg += '<circle cx="' + sx(d.x) + '" cy="' + sy(d.y) + '" r="3.5" fill="var(--identified)"/>';
  });
  // Labels are clamped inside the drawing and staggered, because several tiers sit within
  // a few bytes of each other at the right-hand end and a label that leaves the viewBox is
  // simply missing from the figure.
  pts.forEach((d, i) => {
    const x = sx(d.x), near = x > W - m.r - 80;
    svg += '<text x="' + (near ? W - m.r - 4 : x) + '" y="' +
      (sy(d.y) - 10 - (i % 3) * 12) + '" font-size="10.5" fill="var(--dim)" ' +
      'text-anchor="' + (near ? 'end' : 'middle') + '" font-family="var(--mono)">' +
      esc(d.tier) + '</text>';
  });
  const lo = pts.reduce((a, d) => Math.min(a, d.x), Infinity);
  const hi = pts.reduce((a, d) => Math.max(a, d.x), 0);
  svg += '<text x="' + ((W + m.l) / 2) + '" y="' + (H - 16) + '" text-anchor="middle" ' +
    'font-size="11.5" fill="var(--dim)">cumulative bytes consulted (log scale, ' +
    gb(lo) + ' to ' + gb(hi) + ')</text></svg>';
  p.append(el('div', 'scroll', svg));
  p.append(opChips(t.operating_point));
  app.append(p);
})();

// ---------- mate linkage
(function () {
  const L = D.linkage;
  const p = panel('Mate recovery in the joint space',
    'Each mate is a point in the per-mate joint space, and the question is whether a ' +
    'read\\'s own mate is among its nearest neighbours. It is the only label-free ground ' +
    'truth for "same molecule" a real library carries, and it is an upper bound on ' +
    'grouping rather than an estimate of it: two mates share an insert.');
  if (!L) { p.append(el('div', 'absent',
    'linkage.json absent: the grouping step has no calibration in this run, so no ' +
    'grouping may be read off this page'));
    app.append(p); return; }
  const r = L.result || {};
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>population</th><th>k=1</th><th>k=10</th><th>k=50</th>' +
    '<th>chance at k=10</th></tr></thead>';
  const body = el('tbody');
  const show = {all: 'all mates', mates_disjoint: 'mates sharing no 31-mer (the honest one)',
    mates_overlap: 'mates that overlap', from_mate1: 'from mate 1', from_mate2: 'from mate 2',
    coverage_only: 'coverage axes only', composition_only: 'composition axes only',
    permuted_partner_null: 'permuted-partner null (read this first)'};
  Object.keys(show).forEach(k => {
    const v = r[k]; if (!v) return;
    const tr = el('tr');
    const f = q => v[q] ? pct(v[q].recovery) + ' <span class="dim mono">' +
      (v[q].enrichment ? v[q].enrichment + 'x' : '') + '</span>' : '-';
    tr.innerHTML = '<td>' + esc(show[k]) + '</td><td>' + f('1') + '</td><td>' + f('10') +
      '</td><td>' + f('50') + '</td><td class="mono dim">' +
      (v['10'] ? pct(v['10'].chance) : '-') + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  if (r.scale_ladder) {
    const s2 = el('div', 'scroll');
    const t2 = el('table');
    t2.innerHTML = '<thead><tr><th>pairs in the space</th><th>recovery at 10</th>' +
      '<th>enrichment</th><th>disjoint recovery</th><th>disjoint enrichment</th></tr></thead>';
    const b2 = el('tbody');
    r.scale_ladder.forEach(d => {
      const tr = el('tr');
      tr.innerHTML = '<td class="mono">' + d.pairs + '</td><td class="mono">' +
        pct(d.recovery_at_10) + '</td><td class="mono">' + d.enrichment_at_10 +
        'x</td><td class="mono">' + pct(d.disjoint_recovery_at_10) + '</td><td class="mono">' +
        d.disjoint_enrichment_at_10 + 'x</td>';
      b2.append(tr);
    });
    t2.append(b2); s2.append(t2);
    p.append(el('p', 'note', 'How it behaves as the space fills. This, rather than the ' +
      'level, is what decides whether linkage is usable at library scale.'));
    p.append(s2);
  }
  p.append(el('p', 'note', esc(L.bias || '')));
  app.append(p);
})();

// ---------- layers
(function () {
  const L = D.layers;
  const p = panel('Residue layers, ordered, with controls',
    'A layer is a sort order with its evidence attached. Each sits beside the same ' +
    'rule applied to identified reads; where that control was never computed the page ' +
    'says so.');
  if (!L) { p.append(el('div', 'absent', 'layers.json absent')); app.append(p); return; }
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>layer</th><th>residue</th><th>already identified</th>' +
    '<th>reading</th></tr></thead>';
  const body = el('tbody');
  (L.layer_order_most_to_least_interesting || []).forEach(k => {
    const v = (L.against_known_positive || {})[k] || {};
    const tr = el('tr');
    const kp = v.known_positive && (v.field_coverage || {}).known &&
      Object.values(v.field_coverage.known).every(c => c > 0)
      ? interval(v.known_positive)
      : '<span class="absent">not measured</span>';
    tr.innerHTML = '<td>' + esc(k.replace(/_/g, ' ')) + '</td><td>' +
      interval(v.residue) + '</td><td>' + kp + '</td><td class="dim" ' +
      'style="white-space:normal;text-align:left;max-width:34ch">' +
      esc(v.reading || '') + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  p.append(el('p', 'rule', esc(L.not_a_score || '')));
  const d = el('details');
  d.append(el('summary', null, 'the rule behind each layer'));
  d.append(el('pre', 'mono note', esc(JSON.stringify(L.rule || {}, null, 1))));
  p.append(d);
  app.append(p);
})();

// ---------- field calibration
(function () {
  const F = D.fields || {};
  const fields = F.fields || {};
  const p = panel('Field calibration: null and known-positive',
    'A null says what the field reads on sequence carrying no information. A ' +
    'known-positive says whether it reads differently on sequence the cascade identified. ' +
    '"Cannot discriminate here" is a value in this table, not a blank.');
  const state = (D.run || {}).field_calibration || {};
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>field</th><th>residue</th><th>shuffled null</th>' +
    '<th>identified reads</th><th>state</th></tr></thead>';
  const body = el('tbody');
  const mean = o => o && o.mean != null ? num(o.mean, 4) + ' <span class="dim mono">&plusmn;' +
    num(o.sd, 3) + '</span>' : (o && o.any_hit_frac != null ? pct(o.any_hit_frac) : '-');
  Object.keys(fields).forEach(k => {
    const v = fields[k] || {};
    const st = state[k];
    const txt = typeof st === 'string' ? st : (st && st.state) || 'not recorded';
    const tr = el('tr');
    tr.innerHTML = '<td>' + esc(k) + '</td><td class="mono">' + mean(v.residue || v) +
      '</td><td class="mono">' + mean(v.shuffled_null) + '</td><td class="mono">' +
      (v.subtracted_reads_any_hit_frac != null ? pct(v.subtracted_reads_any_hit_frac)
        : (v.subtracted_reads_bitscore ? num(v.subtracted_reads_bitscore.mean, 2)
          : '<span class="absent">not measured</span>')) +
      '</td><td class="' + (/cannot/.test(txt) ? 'absent' : 'dim') +
      '" style="white-space:normal;text-align:left;max-width:30ch">' + esc(txt) + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  // A string here is the per-read report saying the test did not run; it used to be 0.0,
  // which this sentence rendered as "fires on 0%".
  if (typeof F.chimaera_split_rate === 'string') {
    p.append(el('p', 'note', 'Chimaera test ' + esc(F.chimaera_split_rate) + '.'));
  } else if (F.chimaera_split_rate != null) {
    p.append(el('p', 'note', 'Chimaera fires on ' + pct(F.chimaera_split_rate) +
      ' of residue reads and ' + pct(F.chimaera_split_rate_on_subtracted_reads) +
      ' of reads already identified, which is the test\\'s own false-positive rate. ' +
      'That is why a read-level chimaera call is not made on this page.'));
  }
  p.append(el('p', 'note', esc((D.run || {}).known_positive_note || '')));
  app.append(p);
})();

// ---------- the linkage census
(function () {
  const C = D.census;
  const p = panel('Insert-adjacency labels of residue reads',
    'A paired library answers adjacency for free. For each residue read, the other reads ' +
    'covering the same sequence are found and THEIR other mates are labelled: each is ' +
    'sequence one insert away from something the read also covers. That separates reads ' +
    'that are unidentified because our references are a gene catalogue from reads whose ' +
    'neighbourhood is as unknown as they are.');
  if (!C) { p.append(el('div', 'absent', 'census.json absent: the census did not run'));
    app.append(p); return; }
  const R = C.result || {};
  const order = ['beside_known', 'beside_assembled_unidentified', 'unknown_neighbourhood',
                 'mixed', 'no_witnesses'];
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>the query sits...</th><th>residue</th>' +
    '<th>already identified</th><th>permuted-label null</th></tr></thead>';
  const body = el('tbody');
  order.forEach(k => {
    const a = (R.residue && R.residue.piles || {})[k];
    const b = (R.identified && R.identified.piles || {})[k];
    const nl = (R.residue && R.residue.permuted_label_null || {})[k];
    if (!a && !b) return;
    const tr = el('tr');
    tr.innerHTML = '<td>' + esc(k.replace(/_/g, ' ')) + '</td><td>' + interval(a) +
      '</td><td>' + interval(b) + '</td><td class="dim">' + interval(nl) + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  p.append(el('p', 'note', 'Read the two right-hand columns first. The null says what the ' +
    'piles would read from the neighbour label mix alone, so a pile that matches it is ' +
    'reporting the composition of the library and nothing about adjacency.'));
  p.append(el('p', 'note', esc(C.bias || '')));
  app.append(p);
})();

// ---------- new strains of known viruses
(function () {
  const S = D.strains;
  const p = panel('Candidate divergent strains, grouped by subject',
    'Grouped by nearest relative: a candidate strain is a subject several residue ' +
    'reads point at from a distance. Subjects are ordered by read count and never by ' +
    'identity: ' +
    'identity is a distance from a database, and ordering by it would order by database ' +
    'proximity.');
  if (!S) { p.append(el('div', 'absent',
    'strains.json absent: the per-read report did not run in this directory'));
    app.append(p); return; }
  const NM = D.names || {};
  const rows = (S.subjects || []).slice();
  const maxr = Math.max(1, ...rows.map(r => r.residue_reads || 0));
  const ctl = el('div', 'controls');
  const sel = el('select');
  [['reads', 'most reads pointing at it'],
   ['identity', 'furthest from its relative'],
   ['residue_only', 'seen only in the residue']].forEach(([v, t]) => {
    const o = document.createElement('option'); o.value = v; o.textContent = t;
    sel.append(o); });
  ctl.append(el('span', 'dim', 'order by'), sel);
  p.append(ctl);
  const grid = el('div', 'cards');
  function label(r) {
    const n = NM[r.subject];
    if (!n) return {title: r.subject, sub: null, acc: null};
    const org = n.organism || null, nm = n.name || null;
    return {title: org || nm || r.subject, sub: org && nm ? nm : null, acc: r.subject};
  }
  function draw() {
    const v = sel.value;
    const list = rows.slice().sort((a, b) =>
      v === 'identity' ? (a.median_identity_percent - b.median_identity_percent) :
      v === 'residue_only' ? ((b.known_reads === 0) - (a.known_reads === 0)) ||
                             (b.residue_reads - a.residue_reads) :
      (b.residue_reads - a.residue_reads)).slice(0, 60);
    grid.innerHTML = '';
    list.forEach(r => {
      const L = label(r);
      const c = el('div', 'card');
      const only = S.known_pile_searched_reads && r.known_reads === 0;
      const id = r.median_identity_percent;
      // The identity bar runs 0 to 100 and is filled to the median, so a short bar is a
      // distant relative. It is never used to rank: the ranking control above is explicit
      // about what it sorts by, because ordering by identity orders by database proximity.
      c.innerHTML =
        '<div class="ctitle">' + esc(L.title) + '</div>' +
        (L.sub ? '<div class="csub">' + esc(L.sub) + '</div>' : '') +
        '<div class="cbar"><span style="width:' + Math.max(2, 100 * (r.residue_reads /
          maxr)).toFixed(1) + '%"></span></div>' +
        '<div class="crow"><b>' + r.residue_reads + '</b> residue reads point here' +
        (S.known_pile_searched_reads ? ', <b>' + r.known_reads + '</b> identified reads do'
          : '') + '</div>' +
        '<div class="crow">closest at <b>' + num(id, 1) + '%</b> identity over ' +
        r.median_aligned_aa + ' amino acids' +
        '<div class="idbar"><span style="width:' + Math.max(1, Math.min(100, id)).toFixed(1)
        + '%"></span></div></div>' +
        '<div class="crow dim mono">E ' + (+r.median_evalue).toExponential(1) +
        ' median, ' + (+r.worst_evalue).toExponential(1) + ' worst' +
        (L.acc ? ' &middot; ' + esc(L.acc) : '') + ' &middot; ' + esc(r.index) + '</div>' +
        '<div class="tag ' + (only ? 'tag-only' : 'tag-also') + '">' +
        (!S.known_pile_searched_reads ? 'control not searched' :
         only ? 'only in the residue' : 'also in the identified pile') + '</div>';
      grid.append(c);
    });
  }
  sel.onchange = draw; draw();
  p.append(grid);
  p.append(el('p', 'note', S.subjects_total + ' subjects at E < ' +
    S.significance_evalue + ', the 60 highest shown under the chosen order. Names come ' +
    'from the header of the matched record in the same local database that was searched, ' +
    'so a name is what that record calls itself and never a classification of the read: ' +
    'a read matching a dengue polyprotein over 34 amino acids is not dengue, which is ' +
    'why the identity and the E-value sit beside every name. ' + esc(S.caveat) +
    (S.known_pile_searched_reads ? '' : ' In this run the already-subtracted pile was ' +
      'not searched, so every zero in that column is absence rather than a measurement ' +
      'and no row may be read as residue-only.')));
  app.append(p);
})();

// ---------- the joint space
(function () {
  const S = D.scatter || {};
  const p = panel('Joint-space density, residue against identified',
    'Each surviving read a point, with the already-identified reads drawn underneath as ' +
    'the reference cloud. This is where signal against noise is read off, and it replaces ' +
    'the panel that would have shown a novelty score. It is drawn as a density and never ' +
    'as a clustering: no cluster count beats one against a null with the same mean and ' +
    'covariance.');
  const axes = D.axes || [];
  if (!axes.length || !(S.residue && S.residue.points.length)) {
    p.append(el('div', 'absent', 'reads.jsonl absent: no per-read rows in this run'));
    app.append(p); return;
  }
  const ctl = el('div', 'controls');
  const mk = (id, init) => { const s = el('select'); s.id = id;
    axes.forEach((a, i) => { const o = document.createElement('option');
      o.value = i; o.textContent = a.label; if (i === init) o.selected = true; s.append(o); });
    return s; };
  const sx = mk('axis-x', 0), sy = mk('axis-y', 2);
  ctl.append(el('span', 'dim', 'x'), sx, el('span', 'dim', 'y'), sy);
  p.append(ctl);
  const holder = el('div', 'scroll');
  p.append(holder);
  const W = 920, H = 420, m = {l: 58, r: 14, t: 14, b: 46};
  function draw() {
    const i = +sx.value, j = +sy.value;
    const sets = [['known', 'var(--null)', 0.5], ['residue', 'var(--residue)', 0.75]];
    let vals = [];
    sets.forEach(([k]) => (S[k] ? S[k].points : []).forEach(pt => {
      if (pt[i] != null && pt[j] != null) vals.push([pt[i], pt[j]]); }));
    if (!vals.length) { holder.innerHTML = '<div class="absent">neither axis is ' +
      'populated on these rows</div>'; return; }
    const xr = [Math.min(...vals.map(v => v[0])), Math.max(...vals.map(v => v[0]))];
    const yr = [Math.min(...vals.map(v => v[1])), Math.max(...vals.map(v => v[1]))];
    const px = v => m.l + (v - xr[0]) / ((xr[1] - xr[0]) || 1) * (W - m.l - m.r);
    const py = v => H - m.b - (v - yr[0]) / ((yr[1] - yr[0]) || 1) * (H - m.t - m.b);
    let s = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="joint space">';
    s += '<rect x="' + m.l + '" y="' + m.t + '" width="' + (W - m.l - m.r) + '" height="' +
      (H - m.t - m.b) + '" fill="none" stroke="var(--edge)"/>';
    sets.forEach(([k, col, op]) => {
      (S[k] ? S[k].points : []).forEach(pt => {
        if (pt[i] == null || pt[j] == null) return;
        s += '<circle cx="' + px(pt[i]).toFixed(1) + '" cy="' + py(pt[j]).toFixed(1) +
          '" r="1.7" fill="' + col + '" opacity="' + op + '"/>';
      });
    });
    [[xr[0], 0], [xr[1], 1]].forEach(([v, side]) => {
      s += '<text x="' + (side ? W - m.r : m.l) + '" y="' + (H - m.b + 16) +
        '" font-size="11" fill="var(--dim)" font-family="var(--mono)" text-anchor="' +
        (side ? 'end' : 'start') + '">' + num(v, 2) + '</text>';
    });
    [[yr[0], H - m.b], [yr[1], m.t + 10]].forEach(([v, y]) => {
      s += '<text x="' + (m.l - 8) + '" y="' + y + '" font-size="11" fill="var(--dim)" ' +
        'font-family="var(--mono)" text-anchor="end">' + num(v, 2) + '</text>';
    });
    s += '<text x="' + ((W + m.l) / 2) + '" y="' + (H - 10) + '" text-anchor="middle" ' +
      'font-size="11.5" fill="var(--dim)">' + esc(axes[i].label) + '</text>';
    s += '<text transform="translate(14 ' + ((H - m.b + m.t) / 2) + ') rotate(-90)" ' +
      'text-anchor="middle" font-size="11.5" fill="var(--dim)">' + esc(axes[j].label) +
      '</text></svg>';
    holder.innerHTML = s;
  }
  // The identity axis needs its own warning, because the shape it draws is arithmetic
  // rather than biology and it looked like a defect in the measurement the first time it
  // was seen on this page.
  const idnote = el('p', 'note');
  function notes() {
    const picked = [axes[+sx.value].key, axes[+sy.value].key];
    idnote.hidden = !picked.includes('best_pident');
    idnote.innerHTML = 'Identity is plotted only for hits at E &lt; ' +
      (D.significant_evalue || 1e-5) + ', so a read with no significant relative is ' +
      'absent from this axis rather than shown at its chance identity. Even so the ' +
      'points fall in rows: percent identity over a ' +
      'short alignment is a ratio of small whole numbers, so at 34 aligned amino acids ' +
      'the only values near half are 47.1% and 52.9% and nothing can land between them. ' +
      'The banding is integer arithmetic, not structure in the sample. Unfiltered it is ' +
      'far worse, with dense lines at 33.3%, 50% and 100% drawn almost entirely by ' +
      'chance-level matches.';
  }
  p.append(idnote);
  sx.onchange = () => { draw(); notes(); };
  sy.onchange = () => { draw(); notes(); };
  draw(); notes();
  const lg = el('div', 'legend');
  lg.innerHTML = '<span><i style="background:var(--residue)"></i>residue (' +
    (S.residue ? S.residue.n_drawn : 0) + ' of ' + (S.residue ? S.residue.n_total : 0) +
    ' drawn)</span><span><i style="background:var(--null)"></i>already identified (' +
    (S.known ? S.known.n_drawn : 0) + ' of ' + (S.known ? S.known.n_total : 0) +
    ' drawn)</span>';
  p.append(lg);
  app.append(p);
})();

// ---------- divergence ladder
(function () {
  const C = D.calibration_ladder;
  const p = panel('Holdout spike, divergence ladder',
    'A holdout genome spiked at known divergence through the same cascade. The only ' +
    'positive control here; the rungs are ordinal.');
  if (!C) {
    p.append(el('div', 'absent',
      'no calibration ladder in this run: it was not requested, so no ' +
      'sensitivity-at-divergence figure may be read off this page'));
    app.append(p); return;
  }
  const R = C.result || C;
  const rungs = R.rungs || [];
  if (!rungs.length) {
    p.append(el('pre', 'mono note scroll', esc(JSON.stringify(R, null, 1).slice(0, 3000))));
    app.append(p); return;
  }
  // One line per spiked genome: how much of it the cascade still removes as the genome is
  // walked away from its own reference. The rungs are ordinal, so the x axis is drawn at
  // equal spacing and labelled with the divergence rather than scaled by it.
  const gids = Object.keys(rungs[0].genomes || {});
  const W = 920, H = 380, m = {l: 54, r: 150, t: 16, b: 50};
  const xs = i => m.l + (rungs.length === 1 ? 0 : i / (rungs.length - 1)) * (W - m.l - m.r);
  const ys = v => H - m.b - v * (H - m.t - m.b);
  const cols = ['var(--residue)', 'var(--identified)', 'var(--null)'];
  let s = '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="divergence ladder">';
  [0, 0.25, 0.5, 0.75, 1].forEach(v => {
    s += '<line x1="' + m.l + '" x2="' + (W - m.r) + '" y1="' + ys(v).toFixed(1) +
      '" y2="' + ys(v).toFixed(1) + '" stroke="var(--edge)"/>';
    s += '<text x="' + (m.l - 8) + '" y="' + (ys(v) + 4).toFixed(1) + '" font-size="11" ' +
      'fill="var(--dim)" text-anchor="end" font-family="var(--mono)">' +
      (100 * v).toFixed(0) + '%</text>';
  });
  rungs.forEach((r, i) => {
    s += '<text x="' + xs(i).toFixed(1) + '" y="' + (H - m.b + 18) + '" font-size="11.5" ' +
      'fill="var(--dim)" text-anchor="middle">' + (100 * r.divergence).toFixed(0) +
      '%</text>';
  });
  gids.forEach((gid, gi) => {
    const col = cols[gi % cols.length];
    let d = '';
    rungs.forEach((r, i) => {
      const v = ((r.genomes[gid] || {}).survival || {}).survives_tier1;
      if (v == null) return;
      d += (d ? ' L' : 'M') + xs(i).toFixed(1) + ' ' + ys(v).toFixed(1);
      s += '<circle cx="' + xs(i).toFixed(1) + '" cy="' + ys(v).toFixed(1) +
        '" r="3.5" fill="' + col + '"/>';
    });
    s += '<path d="' + d + '" fill="none" stroke="' + col + '" stroke-width="2"/>';
    // The lines converge at the right, so the end labels are stacked at a fixed pitch
    // instead of sitting on their own line: three labels on one pixel row is three
    // labels nobody can read.
    s += '<text x="' + (W - m.r + 8) + '" y="' + (m.t + 14 + gi * 17) +
      '" font-size="11.5" fill="' + col + '" font-family="var(--mono)">' + esc(gid) +
      '</text>';
  });
  s += '<text x="' + ((W + m.l - m.r) / 2) + '" y="' + (H - 10) + '" text-anchor="middle" ' +
    'font-size="11.5" fill="var(--dim)">divergence from the reference the spike was ' +
    'built on</text>';
  s += '<text transform="translate(14 ' + ((H - m.b + m.t) / 2) + ') rotate(-90)" ' +
    'text-anchor="middle" font-size="11.5" fill="var(--dim)">spike reads surviving the ' +
    'cascade</text></svg>';
  const holder = el('div', 'scroll'); holder.innerHTML = s; p.append(holder);
  p.append(el('p', 'note',
    'A line that stays on the floor is a genome the cascade still catches. Where a line ' +
    'leaves the floor is where this tool stops recognising a relative, and it is the ' +
    'single most important number on the page: a negative above that divergence means ' +
    'nothing. ' + R.spike_pairs_per_rung + ' spiked pairs per rung, inserts drawn from ' +
    'this library\\'s own distribution, sequencing error applied at ' +
    (R.error_rates || []).map(e => num(e, 5)).join(' and ') + ' per base.'));
  const det = el('details');
  det.innerHTML = '<summary>which index removes the spike, rung by rung</summary>';
  const sc2 = el('div', 'scroll');
  const tiers = Object.keys((((rungs[0].genomes[gids[0]] || {}).survival || {})
    .removed_by_tier) || {});
  let th = '<thead><tr><th>divergence</th><th>genome</th>' +
    tiers.map(t2 => '<th>' + esc(t2) + '</th>').join('') + '</tr></thead><tbody>';
  rungs.forEach(r => gids.forEach(gid => {
    const rb = ((r.genomes[gid] || {}).survival || {}).removed_by_tier || {};
    th += '<tr><td class="mono">' + (100 * r.divergence).toFixed(0) + '%</td>' +
      '<td class="mono dim">' + esc(gid) + '</td>' +
      tiers.map(t2 => '<td class="mono">' + (rb[t2] == null ? '-' :
        (100 * rb[t2]).toFixed(1) + '%') + '</td>').join('') + '</tr>';
  }));
  const tb2 = el('table'); tb2.innerHTML = th + '</tbody>';
  sc2.append(tb2); det.append(sc2); p.append(det);
  app.append(p);
})();

// ---------- the distance spectrum
/* ONE AXIS, NO BUCKET EDGE ON IT. Identity against a reference orders contigs by how
   close the database happens to sit, so a line drawn across this axis is a statement
   about the database and not about novelty. The only categorical call the panel makes is
   whether the measurement can be read at all, and those rows leave the axis entirely
   rather than being placed at a number nobody should read. */
(function () {
  const S = D.spectrum;
  const p = panel('The distance spectrum, one contig per row',
    'Every identified contig scored gene by gene, laid out by how far its conserved gene ' +
    'sits from the closest protein two reference sets hold. Far on the left, near on the ' +
    'right, and no boundary anywhere between them. ' + cite('spectrum'));
  if (!S || !S.on_axis || !S.on_axis.length) {
    const why = (((D.page || {}).inputs || {}).spectrum || {}).not_run;
    p.append(el('div', 'absent', 'not run: ' + (why || 'no divergent-stage result in ' +
      'this build')));
    app.append(p); return;
  }
  const RAMP = ['var(--cv1)', 'var(--cv2)', 'var(--cv3)', 'var(--cv4)', 'var(--cv5)'];
  function covColour(c) {
    if (c == null) return 'var(--null)';
    return RAMP[Math.min(4, Math.max(0, Math.floor(c * 5 - 0.0001)))];
  }
  const BIG = {
    'collapsed to a strain-level identification': 'var(--v3)',
    'closer, still not an identification': 'var(--v1)',
    'still distant': 'var(--v2)',
    'no significant hit': 'var(--null)'
  };

  const ctl = el('div', 'controls');
  const armSel = el('select');
  ['both arms', 'divergent only', 'near-reference only'].forEach(v => {
    const o = document.createElement('option'); o.textContent = v; armSel.append(o); });
  const bands = el('input'); bands.type = 'checkbox'; bands.id = 'spec-bands';
  const bl = el('label'); bl.htmlFor = 'spec-bands';
  bl.textContent = 'draw the old three buckets';
  const asTable = el('input'); asTable.type = 'checkbox'; asTable.id = 'spec-table';
  const tl = el('label'); tl.htmlFor = 'spec-table'; tl.textContent = 'as a table';
  ctl.append(el('span', 'dim', 'show'), armSel, bands, bl, asTable, tl);
  p.append(ctl);
  // The floors this figure is read at are controls, because they decide whether a
  // number can be read at all rather than what it says: a row under the coverage floor
  // is not a far-left row, it is a row whose identity cannot be read. Moving them moves
  // rows between the axis and the lane, and every row that moves says both labels.
  const SPECC = ['spec_rank', 'spec_min_coverage', 'spec_min_aln_aa', 'min_contig_nt',
                 'min_depth'];
  p.append(ctlRow(SPECC));
  p.append(settingLine(SPECC));

  // The chart keeps a minimum width and scrolls inside its own box on a narrow screen,
  // rather than shrinking its labels to four pixels. The page itself still never
  // scrolls sideways.
  const host = el('div', 'scroll'); p.append(host);
  const foot = el('div', 'note'); p.append(foot);

  /* Every row re-read at the reader's floors, from the numbers already on it. The axis
     and the lane are one pile re-partitioned, never two piles: a row that crosses
     carries both labels, because a contig that quietly changed pile when a slider moved
     is the failure this is written against. */
  function rowsAt() {
    const whole = ctlGet('spec_rank') === 'whole';
    const minc = ctlGet('spec_min_coverage'), mina = ctlGet('spec_min_aln_aa');
    const minnt = ctlGet('min_contig_nt'), mind = ctlGet('min_depth');
    const out = {axis: [], lane: [], below: 0, crossed: 0};
    (S.on_axis || []).concat(S.cannot_read || []).forEach(r0 => {
      const x = whole ? r0.x_whole : r0.x;
      const cov = whole ? r0.cov_whole : r0.cov;
      const aln = whole ? r0.aln_whole : r0.aln;
      if ((r0.nt != null && r0.nt < minnt) || (r0.depth != null && r0.depth < mind)) {
        out.below++; return;
      }
      const why = x == null ? 'no identity to read at this ranking'
        : (cov == null || cov < minc)
          ? 'under ' + Math.round(minc * 100) + '% of the scored gene covered'
          : (aln == null || aln < mina) ? 'alignment under ' + mina + ' aa' : null;
      const r = Object.assign({}, r0, {x: x, cov: cov, aln: aln, why: why,
        src: whole ? 'whole contig' : r0.src,
        shipped_axis: r0.why == null});
      r.crossed = (why == null) !== r.shipped_axis;
      if (r.crossed) out.crossed++;
      (why ? out.lane : out.axis).push(r);
    });
    return out;
  }
  function armed(rows) {
    const want = armSel.value;
    return rows.filter(r => want === 'both arms' ||
      (want === 'divergent only' ? r.arm === 'divergent'
                                 : r.arm !== 'divergent'));
  }
  function relabel(r) {
    if (!r.crossed) return '';
    return '\\n' + (r.why ? 'in the lane at your setting, on the axis as shipped'
                          : 'on the axis at your setting, in the lane as shipped');
  }
  function rowsNow() { return armed(rowsAt().axis); }

  function drawTable(rows) {
    const sc = el('div', 'scroll');
    const tb = el('table');
    tb.innerHTML = '<thead><tr><th>contig</th><th>arm</th><th>ranked on</th>' +
      '<th>identity</th><th>gene spread</th><th>gene covered</th><th>aligned</th>' +
      '<th>depth</th><th>closest protein</th></tr></thead>';
    const body = el('tbody');
    rows.slice().sort((a, b) => a.x - b.x).slice(0, 400).forEach(r => {
      const tr = el('tr');
      tr.innerHTML = '<td class="mono">' + esc(short(r.name)) + '</td><td class="dim">' +
        (r.arm === 'divergent' ? 'divergent' : 'near-reference') + '</td><td class="dim">' +
        esc(r.src) + '</td><td class="mono">' + num(r.x, 1) + '%</td><td class="mono">' +
        (r.hi == null ? '<span class="dim">one gene</span>'
                      : num(r.hi - r.lo, 1) + ' pts over ' + r.genes) +
        '</td><td class="mono">' + num(r.cov * 100, 0) + '%</td><td class="mono">' +
        r.aln + ' aa</td><td class="mono">' + (r.depth == null ? '-' : num(r.depth, 0)) +
        '</td><td>' + esc(r.subject) + '</td>';
      body.append(tr);
    });
    tb.append(body); sc.append(tb);
    return sc;
  }

  function short(n) {
    const m = String(n).match(/NODE_(\\d+)_length_(\\d+)/);
    return m ? 'NODE ' + m[1] + ' / ' + m[2] + ' nt' : String(n).slice(0, 28);
  }

  function draw() {
    host.innerHTML = '';
    const R = rowsAt();
    const rows = armed(R.axis), lane = armed(R.lane);
    if (asTable.checked) { host.append(drawTable(rows)); setFoot(rows, R, lane); return; }
    const W = 1000, PAD_L = 54, PAD_R = 18, PAD_T = 48;
    const H_MAIN = 300, H_LANE = 92, GAP = 30;
    const H = PAD_T + H_MAIN + GAP + H_LANE + 64;  // room for two legend rows
    const x = v => PAD_L + (Math.min(100, Math.max(20, v)) - 20) / 80 * (W - PAD_L - PAD_R);
    const dmax = 600;
    const y = d => {
      const v = Math.log10(Math.max(1, Math.min(dmax, d == null ? 1 : d)));
      return PAD_T + H_MAIN - v / Math.log10(dmax) * (H_MAIN - 14);
    };
    let s = '<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" ' +
      'style="min-width:720px" role="img" ' +
      'aria-label="contigs by the identity of their conserved gene against mapped depth">';

    if (bands.checked) {
      [[20, 40, 'no name survives'], [40, 90, 'intermediate'], [90, 100, 'strain']]
        .forEach((b, i) => {
          s += '<rect x="' + x(b[0]) + '" y="' + PAD_T + '" width="' +
            (x(b[1]) - x(b[0])) + '" height="' + H_MAIN + '" fill="var(--v-faint)" ' +
            'opacity="' + (i === 1 ? 0.9 : 0.5) + '"/>';
          s += '<text x="' + (x(b[0]) + 6) + '" y="' + (PAD_T + 13) +
            '" class="ax">' + b[2] + '</text>';
        });
    }

    for (let v = 20; v <= 100; v += 10) {
      s += '<line x1="' + x(v) + '" y1="' + PAD_T + '" x2="' + x(v) + '" y2="' +
        (PAD_T + H_MAIN) + '" stroke="var(--edge)" stroke-width="1"/>';
      s += '<text x="' + x(v) + '" y="' + (PAD_T + H_MAIN + 16) +
        '" text-anchor="middle" class="ax">' + v + '%</text>';
    }
    [1, 10, 100].forEach(d => {
      s += '<line x1="' + PAD_L + '" y1="' + y(d) + '" x2="' + (W - PAD_R) + '" y2="' +
        y(d) + '" stroke="var(--edge)" stroke-width="1" stroke-dasharray="2 4"/>';
      s += '<text x="' + (PAD_L - 8) + '" y="' + (y(d) + 4) +
        '" text-anchor="end" class="ax">' + d + 'x</text>';
    });
    s += '<text x="' + (PAD_L - 8) + '" y="' + (PAD_T - 10) +
      '" text-anchor="end" class="ax">depth</text>';

    rows.forEach(r => {
      const yy = y(r.depth);
      if (r.hi != null && r.hi - r.lo > 0.05) {
        s += '<line x1="' + x(r.lo) + '" y1="' + yy + '" x2="' + x(r.hi) + '" y2="' +
          yy + '" stroke="' + covColour(r.cov) + '" stroke-width="2" ' +
          'stroke-linecap="round" opacity="0.55"/>';
      }
      const ring = r.big && BIG[r.big.state];
      s += '<circle cx="' + x(r.x) + '" cy="' + yy + '" r="' + (ring ? 5 : 4) +
        '" fill="' + covColour(r.cov) + '" stroke="' +
        (ring ? ring : 'var(--panel)') + '" stroke-width="' + (ring ? 2 : 1) + '">' +
        '<title>' + esc(short(r.name)) + ' - ' +
        (r.arm === 'divergent' ? 'divergent arm' : 'near-reference arm') + '\\n' +
        num(r.x, 1) + '% identity, ranked on the ' + r.src + ', E=' + r.evalue + '\\n' +
        (r.hi == null ? 'one scorable gene'
          : r.genes + ' scorable genes spanning ' + num(r.lo, 1) + ' to ' +
            num(r.hi, 1) + '%, a spread of ' + num(r.hi - r.lo, 1) + ' points') + '\\n' +
        num(r.cov * 100, 0) + '% of that gene aligned over ' + r.aln + ' aa\\n' +
        'mapped depth ' + (r.depth == null ? 'not measured' : num(r.depth, 0) + 'x') +
        '\\n' + esc(r.subject) +
        (r.big ? '\\nagainst core_nt, a 319 GB nucleotide index: ' + r.big.state +
          (r.big.pident ? ' (' + num(r.big.pident, 1) + '% over ' +
            num(r.big.frac * 100, 0) + '% of the contig)' : '') : '') +
        relabel(r) + '</title></circle>';
    });

    const LY = PAD_T + H_MAIN + GAP + 34;
    s += '<rect x="' + PAD_L + '" y="' + (LY - 26) + '" width="' + (W - PAD_L - PAD_R) +
      '" height="62" fill="var(--v-faint)" rx="3"/>';
    s += '<text x="' + (PAD_L + 8) + '" y="' + (LY - 10) + '" class="ax">' +
      'off the axis: ' + lane.length + ' contigs whose identity cannot be read at ' +
      'this setting' + '</text>';
    const perRow = Math.ceil(Math.max(1, lane.length) / 3);
    const lw = (W - PAD_L - PAD_R - 16) / Math.max(1, perRow);
    lane.forEach((r, i) => {
      s += '<circle cx="' + (PAD_L + 8 + (i % perRow) * lw + lw / 2) + '" cy="' +
        (LY + 2 + Math.floor(i / perRow) * 9) +
        '" r="2.6" fill="var(--null)" opacity="0.85"><title>' +
        esc(short(r.name)) + ' - ' + esc(r.why) +
        '\\nreported identity ' + num(r.x, 1) + '%, which is why it is not drawn at it' +
        relabel(r) + '</title></circle>';
    });
    s += '<text x="' + ((PAD_L + W - PAD_R) / 2) + '" y="' + (PAD_T + H_MAIN + 34) +
      '" text-anchor="middle" class="ax">identity of the conserved gene, against the ' +
      'closest protein in RefSeq viral plus Swiss-Prot</text>';
    const LGY = H - 28;
    let lx = PAD_L;
    s += '<text x="' + lx + '" y="' + (LGY - 1) + '" class="ax">gene covered</text>';
    lx += 78;
    RAMP.forEach((c, i) => {
      s += '<rect x="' + (lx + i * 15) + '" y="' + (LGY - 11) + '" width="13" ' +
        'height="10" rx="2" fill="' + c + '"/>';
    });
    s += '<text x="' + (lx - 4) + '" y="' + (LGY - 1) +
      '" text-anchor="end" class="ax"></text>';
    s += '<text x="' + (lx + 78) + '" y="' + (LGY - 1) + '" class="ax">' +
      '20% to 100%</text>';
    lx += 150;
    s += '<line x1="' + lx + '" y1="' + (LGY - 6) + '" x2="' + (lx + 26) + '" y2="' +
      (LGY - 6) + '" stroke="var(--cv4)" stroke-width="2" stroke-linecap="round"/>';
    s += '<text x="' + (lx + 32) + '" y="' + (LGY - 1) + '" class="ax">' +
      'the spread between the genes of one contig</text>';
    lx = PAD_L;
    const LGY2 = H - 8;
    [['var(--v3)', 'core_nt identifies it'], ['var(--v1)', 'core_nt gets closer'],
     ['var(--v2)', 'core_nt leaves it distant'],
     ['var(--warn)', 'the withheld genome']].forEach(g => {
      s += '<circle cx="' + (lx + 6) + '" cy="' + (LGY2 - 6) + '" r="4.5" ' +
        'fill="var(--panel)" stroke="' + g[0] + '" stroke-width="2"/>';
      s += '<text x="' + (lx + 15) + '" y="' + (LGY2 - 1) + '" class="ax">' + g[1] +
        '</text>';
      lx += 14 + g[1].length * 5.6;
    });
    s += '</svg>';
    host.innerHTML = s;
    setFoot(rows, R, lane);
  }

  function setFoot(rows, R, lane) {
    const c = S.census || {};
    const thin = lane.filter(r => /covered/.test(r.why)).length;
    const shortAln = lane.length - thin;
    foot.innerHTML =
      '<b>How to read it.</b> Position is the identity of the conserved gene - the ' +
      'polymerase on ' + rows.filter(r => r.src === 'polymerase').length +
      ' of ' + rows.length + ' contigs drawn, the whole-contig hit on the rest, stated ' +
      'per row - because a contig that carries one identity carries whichever gene won ' +
      'on bitscore, and the genes of a divergent contig disagree by a median 12.8 ' +
      'points. The ' +
      'bar through each dot is that disagreement. Colour is how much of the scored gene ' +
      'the alignment spans, pale to deep, because identity flatters towards familiar by ' +
      '+17.2 points when coverage is under half against +7.1 above 90%. Height is ' +
      'mapped depth. ' +
      '<b>No line is drawn across the axis</b>, because identity against a reference ' +
      'orders contigs by database proximity and any threshold on it would be a ' +
      'statement about the database. ' +
      '<b>The lane underneath is not a middle bucket</b>: it holds the ' +
      lane.length + ' contigs where the measurement failed - ' + thin +
      ' with under ' + Math.round(ctlGet('spec_min_coverage') * 100) +
      '% of the scored gene aligned and ' + shortAln + ' aligned over under ' +
      ctlGet('spec_min_aln_aa') + ' aa - and it is a census, not a verdict. ' +
      '<b>A ring</b> marks the ' + S.on_axis.filter(r => r.big).length + ' contigs that ' +
      'also went to core_nt, 319 GB: green collapsed to a strain-level identification ' +
      'there, blue got closer without reaching one, orange stayed distant. The ring ' +
      'is a label and not a move along this axis, ' +
      'because core_nt is a nucleotide index and this axis is protein identity; a ' +
      'displacement arrow waits on the protein pass against nr, which is the same ' +
      'instrument at a larger size. ' +
      (S.structure_axis ? '<b>Structure.</b> ' + esc(S.structure_axis.state) + ': ' +
        esc(S.structure_axis.what_would_fill_it) + '. ' : '') +
      '<b>Scope.</b> ' + (c.searched || '?') + ' identified contigs searched, ' +
      (c.with_a_significant_hit || '?') + ' with a significant protein hit and ' +
      (c.shuffled_with_a_significant_hit === 0 ? 'none at all from ' + c.searched +
        ' dinucleotide shuffles of the same contigs'
        : c.shuffled_with_a_significant_hit + ' from the shuffled null') + '. This panel ' +
      'draws the ' + c.drawn + ' scored gene by gene; ' + rows.length + ' shown at the ' +
      'current filter. ' +
      '<b>At your setting.</b> ' + R.axis.length + ' on the axis and ' + R.lane.length +
      ' in the lane, against ' + S.on_axis.length + ' and ' + S.cannot_read.length +
      ' as shipped' +
      (R.crossed ? ', with ' + R.crossed + ' contigs crossing between them - each one ' +
        'says both labels where it is drawn' : ', with none crossing between them') +
      (R.below ? ', and ' + R.below + ' held back by the length and depth floors' : '') +
      '.';
  }

  armSel.onchange = draw; bands.onchange = draw; asTable.onchange = draw;
  ctlOn(draw);
  // The shuffled contigs were searched once, at the operating point the run used, so the
  // null travels with the figure only at the shipped floors. Past that edge the panel
  // says so rather than leaving a clean-looking count of what cannot be read.
  p.append(nullLine(SPECC, {
    shipped: 'none of the ' + ((S.census || {}).searched || '?') + ' dinucleotide ' +
      'shuffles of these same contigs produced a significant protein hit, at the floors ' +
      'this figure ships at.',
    moved: 'the shuffled contigs were searched once, at the shipped floors. Nothing was ' +
      'measured at this one, so the lane and the axis counts here have no null beside ' +
      'them; the arms drawn in the same pass are the two real ones.'}));
  draw();
  app.append(p);
})();

// ---------- the residue table
(function () {
  const T = D.table;
  const p = panel('Per-read residue table, un-thresholded',
    'One row per surviving pair, un-thresholded. The closest relative is printed with its ' +
    'E-value every time: an identity without one is the claim this project has already had ' +
    'to withdraw, because a 107 nt query against a large database finds about 40 bases at ' +
    'about 88% by chance.');
  if (!T || !T.rows.length) { p.append(el('div', 'absent', 'no rows in this run'));
    app.append(p); return; }
  const ctl = el('div', 'controls');
  const f = el('select'); f.id = 'layer-filter';
  ['all layers', ...new Set(T.rows.map(r => r.layer).filter(Boolean))]
    .forEach(v => { const o = document.createElement('option'); o.textContent = v;
      f.append(o); });
  const hide = el('input'); hide.type = 'checkbox'; hide.id = 'hide-undet';
  const hl = el('label'); hl.htmlFor = 'hide-undet';
  hl.textContent = 'hide the rows the fields cannot discriminate';
  ctl.append(el('span', 'dim', 'show'), f, hide, hl);
  p.append(ctl);
  /* THE FLOOR DOES NOT DELETE A ROW HERE, IT RE-LABELS IT. A read whose relative clears
     the floor the reader set would have been subtracted as identified at that setting,
     and the honest rendering is to say both: at your setting, and as shipped. A row that
     quietly changed pile when a slider moved is the failure this is written against.
     `cannot discriminate here` is a value, so those rows are never filtered away by
     default; the switch that hides them ships off. */
  const TIDS = ['evalue_floor', 'read_depth_stratum', 'min_occupancy'];
  p.append(ctlRow(TIDS));
  const sc = el('div', 'scroll'); p.append(sc);
  const note = el('p', 'note'); 
  function inStratum(d) {
    const st = ctlGet('read_depth_stratum');
    if (st === 'all' || d == null) return st === 'all';
    if (st === '1') return d <= 1;
    if (st === '2-9') return d >= 2 && d <= 9;
    if (st === '10-99') return d >= 10 && d <= 99;
    return d >= 100;
  }
  function draw() {
    const want = f.value;
    const e = ctlGet('evalue_floor');
    const minocc = ctlGet('min_occupancy');
    let relabelled = 0, undet = 0, held = 0;
    const tb = el('table');
    tb.innerHTML = '<thead><tr><th>read</th><th>layer</th><th>closest relative</th>' +
      '<th>identity</th><th>E</th><th>occurrences</th><th>profile</th></tr></thead>';
    const body = el('tbody');
    const rows = T.rows.filter(r => {
      if (!(want === 'all layers' || r.layer === want)) return false;
      if (hide.checked && (r.layer === 'undetermined' || !r.layer)) return false;
      if (!inStratum(r.depth)) { held++; return false; }
      if (minocc && (r.occupancy == null || r.occupancy < minocc)) { held++; return false; }
      return true;
    });
    rows.slice(0, 400).forEach(r => {
      let best = null;
      Object.values(r.closest || {}).forEach(v => {
        if (v && v.evalue != null && (!best || v.evalue < best.evalue)) best = v; });
      const tr = el('tr');
      const occ = r.values[2];
      const nm = best && (D.names || {})[best.subject];
      const shown = !best ? '<span class="dim">none</span>'
        : nm && (nm.organism || nm.name)
          ? esc(nm.organism || nm.name) +
            ' <span class="dim mono">' + esc(best.subject) + '</span>'
          : '<span class="mono">' + esc(best.subject) + '</span>';
      const clears = best && best.evalue != null && best.evalue <= e;
      if (clears) relabelled++;
      if (r.layer === 'undetermined' || !r.layer) undet++;
      tr.innerHTML = '<td class="mono">' + esc(r.name) + '</td><td class="dim">' +
        esc((r.layer || 'cannot discriminate here').replace(/_/g, ' ')) +
        (clears ? '<br><span class="relabel">tier 1 at your setting, residue as ' +
          'shipped</span>' : '') + '</td><td>' + shown +
        '</td><td class="mono">' + (best ? num(best.pident, 1) + '%' : '-') +
        '</td><td class="mono">' + (best ? best.evalue.toExponential(1) : '-') +
        '</td><td class="mono">' + (occ == null ? '-' : Math.round(Math.pow(2, occ))) +
        '</td><td class="mono">' + (r.motif ? esc(r.motif.profile) + ' <span class="dim">' +
          r.motif.evalue.toExponential(1) + '</span>' : '<span class="dim">-</span>') +
        '</td>';
      body.append(tr);
    });
    tb.append(body);
    sc.innerHTML = ''; sc.append(tb);
    note.innerHTML = T.n_drawn + ' of ' + T.n_total + ' residue pairs sampled for ' +
      'this page with seed ' + T.sample_seed + '; ' + rows.length + ' pass the filters ' +
      'and the first 400 are shown. ' +
      (held ? held + ' rows are held back by the depth stratum and the occupancy cut. '
            : '') +
      (relabelled ? '<span class="relabel">' + relabelled + ' of the rows shown would ' +
        'have been subtracted as identified at your floor and are labelled as such, not ' +
        'moved.</span> ' : '') +
      (Object.keys(D.not_searched || {}).length
        ? '<b>Not searched on this run:</b> ' +
          Object.keys(D.not_searched).map(k => esc(k) + ' (' +
            esc(D.not_searched[k]) + ')').join(', ') + '. A column that is empty for ' +
          'that reason is a field nobody looked at, which is not the same as a field ' +
          'that found nothing. '
        : '') +
      undet + ' of the rows shown are <i>undetermined</i> - the fields cannot ' +
      'discriminate on them - which is a value and not a missing one, so it is never ' +
      'filtered away unless the switch above is turned on. The full rows are in ' +
      'reads.jsonl in the run directory.';
  }
  f.onchange = draw; hide.onchange = draw; ctlOn(draw); draw();
  p.append(note);
  p.append(settingLine(TIDS));
  app.append(p);
})();

// ---------- blind spots
(function () {
  const B = D.blind_spots;
  const p = panel('Blind spots',
    'Each entry names the measurement that established it.' + cite('nulls'));
  if (!B) { p.append(el('div', 'absent', 'blind_spots.json absent from the run directory'));
    app.append(p); return; }
  const g = el('div', 'grid2');
  B.forEach(b => {
    const d = el('div', 'bs');
    d.append(el('h3', null, esc(b.statement)));
    d.append(el('p', null, '<b>Evidence.</b> ' + esc(b.evidence) + '<br><b>So.</b> ' +
      esc(b.consequence)));
    g.append(d);
  });
  p.append(g);
  app.append(p);
})();

// ---------- provenance
(function () {
  const b = D.built || {}, r = D.run || {};
  const p = panel('Provenance', null);
  const chips = opChips({'built by': b.script, 'script sha256': (b.script_sha256 || '').slice(0, 12),
    'git sha': (b.git_sha || '').slice(0, 10), 'tree dirty': String(!!b.tree_dirty),
    'run directory': b.rundir, 'run wall seconds': r.wall_seconds});
  p.append(chips);
  const d = el('details');
  d.append(el('summary', null, 'every stage, with the code that produced it'));
  d.append(el('pre', 'mono note scroll', esc(JSON.stringify(r.stages || [], null, 1))));
  p.append(d);
  // What this page read from outside the run directory, as its descriptor declared it.
  // A "not run" here is the reason a panel says so, stated once where it can be checked.
  const pin = (D.page || {}).inputs || {};
  const d3 = el('details');
  d3.append(el('summary', null, 'files outside the run directory, as this page declares them'));
  d3.append(el('pre', 'mono note scroll', esc(Object.keys(pin).map(k =>
    k + ': ' + (pin[k].path || 'not run - ' + pin[k].not_run)).join('\\n'))));
  p.append(d3);
  app.append(p);
})();

// ---------- what a leftover read sits next to (a split, with its control)
(function () {
  const C = D.census, r = C && C.result;
  if (!r || !r.residue || !r.identified) return;
  const p = panel('Coverage-neighbour labels, residue against identified',
    'Take a leftover pair, find every other read in the library that covers the same ' +
    'sequence, and ask what those neighbours were called. The same question asked of ' +
    'reads the cascade already named is the row underneath, and it is the row to read ' +
    'first: it is what this measurement looks like on material that is not a mystery.');
  const ORDER = [
    ['unknown_neighbourhood', 'neighbours, none of them namable', 'var(--residue)'],
    ['beside_assembled_unidentified', 'beside a contig nothing could name', 'var(--tier2)'],
    ['no_witnesses', 'no neighbours at all', 'var(--null)'],
    ['beside_known', 'beside named sequence', 'var(--identified)'],
    ['mixed', 'both kinds at once', 'var(--g-unknown)']];
  const row = (src, title, control) => {
    const piles = src.piles || {};
    return {title: title, control: control,
      total: src.queries || 0,
      slices: ORDER.map(([k, lab, col]) => {
        const v = piles[k] || {};
        return {k: k, lab: lab, col: col, v: v.k || 0,
                ci: (v.lo != null && v.hi != null) ? [v.lo, v.hi] : null};
      })};
  };
  p.append(fanFig([row(r.residue, 'leftover pairs', false),
                   row(r.identified, 'pairs the cascade already named (the control)', true)],
                  {alt: 'what a leftover read sits next to, against the same split on ' +
                        'reads that were already named'}));
  const res = r.residue.piles || {}, idp = r.identified.piles || {};
  const note = el('div', 'callout');
  note.innerHTML = '<b>The split itself is the finding.</b> A leftover read sits in a ' +
    'neighbourhood of sequence that is also unnamable ' +
    num(100 * (res.unknown_neighbourhood || {}).point, 1) + '% of the time, against ' +
    num(100 * (idp.unknown_neighbourhood || {}).point, 1) + '% for a named read &mdash; ' +
    'so the unnamable material keeps company with itself rather than being scattered ' +
    'through ordinary reads. The other end matters as much: ' +
    num(100 * (res.no_witnesses || {}).point, 1) + '% of leftover pairs have no ' +
    'neighbours at all, and nothing in this figure or anywhere else on this page can ' +
    'say anything about those beyond what one molecule carries.';
  p.append(note);
  app.append(p);
})();
// ---------- the layers as a split, against the known pile
(function () {
  const L = D.layers;
  if (!L || !L.counts || !L.counts.residue) return;
  const p = panel('Layer assignment against the identified control',
    'Each leftover read put in the most interesting layer it qualifies for, most ' +
    'interesting first. The same rule run over reads the cascade already named is ' +
    'underneath: a layer that fills up equally on both is a layer that is describing ' +
    'the library rather than the mystery.');
  const COL = {recurs_and_unplaceable: 'var(--residue)', distant_relative: 'var(--tier2)',
    profile_only: 'var(--g-unknown)', placed: 'var(--identified)',
    artefact_shaped: 'var(--g-vert)', undetermined: 'var(--null)'};
  const LAB = {recurs_and_unplaceable: 'recurs, and nothing can place it',
    distant_relative: 'a distant relative exists', profile_only: 'only a profile sees it',
    placed: 'placed after all', artefact_shaped: 'artefact shaped',
    undetermined: 'nothing in this run decides it'};
  const order = L.layer_order_most_to_least_interesting || Object.keys(COL);
  const row = (counts, title, control) => ({
    title: title, control: control,
    total: Object.keys(counts).reduce((a, k) => a + counts[k], 0),
    slices: order.map(k => ({k: k, lab: LAB[k] || k, col: COL[k] || 'var(--null)',
                             v: counts[k] || 0}))});
  p.append(fanFig([row(L.counts.residue, 'leftover reads', false),
                   row(L.counts.known || {}, 'already-named reads (the control)', true)],
                  {alt: 'the leftover reads laid out in layers, against the same layers ' +
                        'on reads that were already named'}));
  const rc = L.counts.residue, kc = L.counts.known || {};
  const tr = Object.keys(rc).reduce((a, k) => a + rc[k], 0) || 1;
  const tk = Object.keys(kc).reduce((a, k) => a + kc[k], 0) || 1;
  const note = el('div', 'callout');
  note.innerHTML = '<b>Read the control first.</b> The top layer takes ' +
    num(100 * (rc.recurs_and_unplaceable || 0) / tr, 1) + '% of the leftovers and ' +
    num(100 * (kc.recurs_and_unplaceable || 0) / tk, 1) + '% of the already-named reads, ' +
    'so recurrence on its own is not what separates them. What does separate them is the ' +
    'bottom layer: ' + num(100 * (rc.undetermined || 0) / tr, 1) + '% of leftovers ' +
    'against ' + num(100 * (kc.undetermined || 0) / tk, 1) + '% of named reads fall ' +
    'through every rule this run has. <b>This ordering is a layout, not a score.</b>';
  p.append(note);
  app.append(p);
})();

// ---------- public health: what is here that we watch for
// THE FIRST QUESTION, and the only one on this page a reader might act on. It is
// answered with a ranked list and a level, not a table of everything: 27 of the 44
// watchlist taxa were not seen, and 27 rows of "not found" at the top of a page is a
// list nobody reads. They are all still here, one toggle away, because a taxon looked
// for and not found is a result.
(function () {
  const P = D.profile;
  const p = panel('Watchlist screen',
    '44 watchlist taxa against the whole library. Levels are read pairs, not cases.');
  if (!P) {
    p.append(el('div', 'absent', 'public_health_profile_s*.json absent: run ' +
      'experiments/public_health_profile.py over this library'));
    app.append(p); return;
  }
  const CONF = '0.1';
  const at = P.at_confidence[CONF] || {};
  const rows = at.panel || [];
  const strictRows = {};
  ((P.at_confidence['0.2'] || {}).panel || []).forEach(r => { strictRows[r.name] = r; });
  const hit = rows.filter(r => r.detected).sort((a, b) => b.pairs - a.pairs);
  const miss = rows.filter(r => !r.detected);
  const nul = P.shuffled_null || {};
  const top = hit[0];

  // The headline is a LEVEL, because that is what separates a finding from a trace and
  // it is the thing a reader has to carry away. Written as one in N: every rate here is
  // around one pair in a million, where a percentage prints 0.00% and says nothing.
  const oneIn = x => {
    if (!x) return '-';
    const n = 1 / x;
    return n >= 1e6 ? '1 in ' + (n / 1e6).toFixed(1) + ' million'
         : n >= 1e3 ? '1 in ' + Math.round(n / 1e3) + ',000'
         : '1 in ' + Math.round(n);
  };
  // Rows at or above one pair in a hundred thousand, which the sentence below counts.
  const above = hit.filter(r => r.fraction && r.fraction >= 1e-5);
  const hero = el('div', 'verdict');
  hero.innerHTML =
    '<div class="vmain">' + (top
      ? 'Highest watchlist detection: <b>' + esc(top.name) + '</b>, at <b>' +
        oneIn(top.fraction) + ' read pairs</b>.' + cite('watchlist')
      : 'No watchlist taxon detected.' + cite('watchlist')) + '</div>' +
    // COUNTED, NOT ASSERTED. This read `No row reaches 1 pair in 100,000`, which was
    // true of the library the sentence was written on and false on the page it was read
    // from: the top row directly above it said 1 in 44,000. A claim about the rows is
    // computed from the rows. `Untargeted` went the same way as the heading's
    // `untreated`: the preparation of an unseen library is not a fact this code holds.
    '<div class="vsub">' + hit.length + ' of ' + rows.length + ' detected at ' +
    'confidence ' + CONF + '. ' + (above.length
      ? above.length + (above.length === 1 ? ' row reaches' : ' rows reach') +
        ' 1 pair in 100,000 or more.'
      : 'No row reaches 1 pair in 100,000.') +
    ' RNA sequencing of wastewater has low sensitivity for human pathogens at these ' +
    'concentrations.' + cite('library') + '</div>';
  p.append(hero);

  const rank = el('div', 'ranked');
  const maxv = (top && top.fraction) || 1;
  hit.forEach((r, i) => {
    const holds = ((strictRows[r.name] || {}).pairs || 0) > 0;
    const row = el('div', 'rrow');
    row.innerHTML =
      '<div class="rno mono">' + (i + 1) + '</div>' +
      '<div class="rname"><b>' + esc(r.name) + '</b> <span class="dim">' +
        esc(r.group) + '</span>' + (holds
          ? '<span class="chip strong">holds at confidence 0.2</span>'
          : '<span class="chip">trace</span>') + '</div>' +
      // One hue for every bar. Length is the level and the chip is the operating
      // point; colouring the bar by robustness as well put the two encodings in
      // conflict, with the longest bars drawn palest.
      '<div class="rbar"><i style="width:' +
        Math.max(2, 100 * r.fraction / maxv).toFixed(1) + '%"></i></div>' +
      '<div class="rval mono">' + oneIn(r.fraction) + '</div>';
    row.title = r.pairs.toLocaleString() + ' read pairs. ' + r.why_on_the_panel;
    rank.append(row);
  });
  p.append(rank);

  const more = el('details', 'more');
  const sum = document.createElement('summary');
  sum.innerHTML = miss.length + ' not detected';
  more.append(sum);
  const mgrid = el('div', 'misses');
  miss.forEach(r => {
    const m = el('div', 'miss');
    m.innerHTML = '<b>' + esc(r.name) + '</b><span class="dim"> under ' +
      oneIn(r.upper_bound_if_absent) + '</span>';
    m.title = r.why_on_the_panel;
    mgrid.append(m);
  });
  more.append(mgrid);
  more.append(el('p', 'note',
    'Bounds are rule-of-three upper limits at this sequencing depth.' +
    cite('watchlist') + ' Absence requires a targeted assay.'));
  p.append(more);

  const cal = el('div', 'callout');
  cal.innerHTML = '<b>Shuffled control.</b>' + cite('nulls') +
    ' The same panel over composition-matched ' +
    'shuffled reads detects <b>' + ((nul['0.0'] || {}).panel_rows_detected) + ' of ' +
    rows.length + ' rows at confidence 0</b> and ' +
    ((nul[CONF] || {}).panel_rows_detected) + ' at confidence ' + CONF + '. ' +
    'Misassignment to a close relative is uncontrolled. Treat a one-pair row as ' +
    'grounds for a targeted test.';
  p.append(cal);
  app.append(p);
})();

// ---------- public health: what this water is
(function () {
  const P = D.profile;
  const p = panel('What this water is',
    'Kraken2 PlusPF-08 over all read pairs.' + cite('composition') +
    ' First circle: the fraction carrying any name. Second: the composition of that ' +
    'fraction.');
  if (!P) { p.append(el('div', 'absent', 'public_health_profile_s*.json absent'));
    app.append(p); return; }

  const confs = Object.keys(P.at_confidence).sort();
  const ctl = el('div', 'controls');
  const sel = document.createElement('select');
  confs.forEach(c => { const o = document.createElement('option');
    o.value = c;
    o.textContent = c === '0.0' ? 'loose (0)' : c === '0.2' ? 'strict (0.2)'
                  : 'working (' + c + ')';
    sel.append(o); });
  sel.value = confs.indexOf('0.1') >= 0 ? '0.1' : confs[0];
  ctl.append(el('span', 'dim', 'how sure the classifier has to be'), sel);
  p.append(ctl);

  // Two donuts, because there are two part-to-whole questions here and putting them on
  // one ring would mean a 93% slice that hides everything else. The left is the whole
  // library; the right renormalises over the part that got a name, which is the only
  // honest way to show the bacteria-to-virus split and is labelled as such.
  const GROUPS = [
    ['Bacteria', 'bacteria', 'var(--v1)'],
    ['Viruses', 'viruses', 'var(--v2)'],
    ['other living things', 'other living things', 'var(--v3)'],
    ['shallow', 'named, but only to a broad group', 'var(--v-mute)'],
  ];
  const wrap = el('div', 'donuts');
  p.append(wrap);
  const taxaBox = el('div', 'twocol');
  p.append(taxaBox);

  function donut(slices, centre, centreSub, caption) {
    const R = 54, SW = 22, C = 2 * Math.PI * R;
    const box = el('div', 'donut');
    let off = 0, g = '';
    const tot = slices.reduce((a, s) => a + s.v, 0) || 1;
    slices.forEach((s, i) => {
      const len = C * s.v / tot;
      // A 2px surface gap between fills, so adjacent slices never touch.
      // The dash geometry is handed to CSS as custom properties, not as presentation
      // attributes: the reveal animation sets stroke-dashoffset, and a keyframe beats
      // an attribute, so an attribute here would be overwritten and every slice would
      // start at twelve o'clock on top of the others. Drawn once that way.
      const drawn = Math.max(0, len - 2);   // 2px surface gap between fills
      g += '<circle class="arc" r="' + R + '" cx="70" cy="70" fill="none" stroke="' +
        s.col + '" stroke-width="' + SW + '" style="--len:' + drawn.toFixed(2) +
        ';--gap:' + (C - drawn).toFixed(2) + ';--off:' + (-off).toFixed(2) +
        ';--d:' + (90 * i) + 'ms"><title>' + esc(s.lab + ': ' + s.pct) +
        '</title></circle>';
      off += len;
    });
    box.innerHTML = '<svg viewBox="0 0 140 140" role="img" aria-label="' + esc(caption) +
      '"><g transform="rotate(-90 70 70)">' + g + '</g>' +
      '<text x="70" y="68" text-anchor="middle" font-size="21" font-weight="600" ' +
      'fill="var(--fg)">' + esc(centre) + '</text>' +
      '<text x="70" y="84" text-anchor="middle" font-size="9.5" fill="var(--dim)">' +
      esc(centreSub) + '</text></svg>' +
      '<div class="dcap">' + esc(caption) + '</div>';
    return box;
  }

  function taxaList(title, list, total, col) {
    const box = el('div', 'tcol');
    box.append(el('div', 'pghead', title));
    if (!list || !list.length) {
      box.append(el('div', 'dim', 'nothing placed in this group'));
      return box;
    }
    const top = list[0].pairs || 1;
    const draw = (t, hidden) => {
      const r = el('div', 'trow' + (hidden ? ' hid' : ''));
      const pc = 100 * t.pairs / total;
      r.innerHTML = '<div class="tname">' + esc(t.name) + '</div>' +
        '<div class="tbar"><span style="width:' +
        Math.max(2, 100 * t.pairs / top).toFixed(1) + '%;background:' + col +
        '"></span></div>' +
        '<div class="tnum mono">' + t.pairs.toLocaleString() + '</div>';
      r.title = t.pairs.toLocaleString() + ' pairs of ' + total.toLocaleString() +
        ' (' + (pc >= 0.01 ? pc.toFixed(3) : pc.toExponential(1)) + '%)';
      return r;
    };
    list.slice(0, 8).forEach(t => box.append(draw(t, false)));
    if (list.length > 8) {
      const d = el('details', 'more');
      const s2 = document.createElement('summary');
      s2.textContent = (list.length - 8) + ' more';
      d.append(s2);
      list.slice(8).forEach(t => d.append(draw(t, false)));
      box.append(d);
    }
    return box;
  }

  function draw() {
    const a = P.at_confidence[sel.value] || {};
    const total = a.pairs || 1;
    const by = {};
    (a.by_superkingdom || []).forEach(g => { by[g.group] = g.pairs; });
    const unnamed = by.unclassified || 0;
    const named = total - unnamed;
    const shallow = by['classified, but above any of these groups'] || 0;
    const other = (by.Eukaryota || 0) + (by.Archaea || 0);
    wrap.innerHTML = '';
    wrap.append(donut(
      [{lab: 'put on a named organism', v: named, col: 'var(--v1)',
        pct: (100 * named / total).toFixed(2) + '% of the library'},
       {lab: 'nothing in this reference could name it', v: unnamed,
        col: 'var(--v-faint)', pct: (100 * unnamed / total).toFixed(2) + '%'}],
      (100 * named / total).toFixed(1) + '%', 'can be named',
      'of all ' + total.toLocaleString() + ' read pairs, how much any reference names'));
    const inner = [
      {lab: 'bacteria', v: by.Bacteria || 0, col: 'var(--v1)'},
      {lab: 'viruses', v: by.Viruses || 0, col: 'var(--v2)'},
      {lab: 'other living things', v: other, col: 'var(--v3)'},
      {lab: 'named only to a broad group', v: shallow, col: 'var(--v-mute)'},
    ].filter(s => s.v > 0);
    // A SHARE THAT IS NOT ZERO IS NEVER PRINTED AS ZERO. The viral share of the named
    // part rounded to `0%` in the middle of the ring on a library whose own virus table,
    // directly below, lists eight named taxa and hundreds of reads. A reader takes the
    // big number in the ring, and `0%` beside a list of viruses is the page arguing with
    // itself. Anything present but below the printed precision is shown as a bound.
    const small = (v, whole, dp) => {
      const share = 100 * v / (whole || 1);
      const floor = Math.pow(10, -dp);
      return v > 0 && share < floor ? '<' + floor.toFixed(dp) + '%'
                                    : share.toFixed(dp) + '%';
    };
    inner.forEach(s => { s.pct = small(s.v, named, 1) + ' of the named part'; });
    const d2 = donut(inner,
      small(by.Viruses || 0, named, 1), 'of that is viral',
      'and what that named part is');
    // The legend sits under its own donut, never off to the side, and every slice is
    // named in text: one of these three hues carries a contrast warning against the
    // light surface, which obliges a visible label rather than colour alone.
    const lg = el('div', 'donutlg');
    lg.innerHTML = inner.map(s => '<div><i style="background:' + s.col + '"></i>' +
      '<span>' + esc(s.lab) + '</span><span class="v">' + s.pct.split(' ')[0] +
      '</span></div>').join('');
    d2.append(lg);
    wrap.append(d2);

    taxaBox.innerHTML = '';
    taxaBox.append(
      taxaList('the viruses, most reads first', (a.top_taxa || {}).Viruses, total,
               'var(--v2)'),
      taxaList('the bacteria, most reads first', (a.top_taxa || {}).Bacteria, total,
               'var(--v1)'));
  }
  sel.onchange = draw; draw();

  const note = el('div', 'callout');
  note.innerHTML = '<b>Caveats.</b>' + cite('composition') +
    ' This is RNA abundance rather than cell counts: one ' +
    'bacterium with a thousand ribosomes contributes a thousand times its cell count. ' +
    'The named fraction runs from ' +
    pct(1 - (((P.at_confidence['0.0'] || {}).by_superkingdom || [])
      .filter(g => g.group === 'unclassified')[0] || {}).fraction) +
    ' at the loosest setting to ' +
    pct(1 - (((P.at_confidence['0.2'] || {}).by_superkingdom || [])
      .filter(g => g.group === 'unclassified')[0] || {}).fraction) +
    ' at the strictest. Unnamed sequence is not evidence of novelty; see below.';
  p.append(note);
  app.append(p);
})();

// ---------- public health: known viruses, reference or drifted
// RANKED ON THE POLYMERASE, not on a whole-contig best hit. That hit is whichever gene
// won on bitscore, and the control arm says how much that choice moves a row: on contigs
// matching their reference a contig's own genes agree to 2.35 points, on the divergent
// arm they disagree by 17.4. So the page shows the gene it ranked on, the spread, and the
// genes themselves.
(function () {
  const S = D.strain_divergence, W = D.contig_workup;
  const p = panel('How far the assembled viruses sit from anything published',
    'One line. Left is far from everything known, right is a match to a published ' +
    'genome. Every assembled virus has a place on it.' + cite('contigid'));
  if (!S) {
    p.append(el('div', 'absent', 'known_strain_divergence_s*.json absent: run ' +
      'experiments/known_strain_divergence.py'));
    app.append(p); return;
  }
  const taxa = (S.taxa || []).filter(t => t.median_identity_percent != null)
    .sort((a, b) => a.median_identity_percent - b.median_identity_percent);
  const belowNinety = taxa.filter(t => t.median_identity_percent < 90);
  const rows = W ? (W.rows || []).filter(r => r.arm === 'divergent')
    .slice().sort((a, b) => a.rank_pident - b.rank_pident) : [];
  const ctrl = W ? (W.summary || {}).control_near_reference : null;
  const dsum = W ? (W.summary || {}).divergent : null;
  const deep = rows.filter(r => (r.read_support || {}).mapped_mean_depth >= 100
                                && r.contig_nt >= 3000);

  /* A length claim with no chimaera number beside it is the thing CLAUDE.md refuses,
     and this page was making one: the flag was computed by the stage, carried in the
     payload and read by nothing, so the Methods entry describing it sat uncited while
     the section above it counted complete genomes. Both arms are printed from the same
     summary, flagged over the contigs that could be flagged, and the sentence says what
     the comparison means rather than leaving two percentages side by side. A missing
     denominator is reported as not testable, never as 0%. */
  /* A difference of two fractions is points, not a percent: "within 15.38%" invites a
     reader to take it as a ratio of the two rates. */
  function points(g) { return (g * 100).toFixed(1) + ' percentage points'; }

  function chimaeraLine() {
    const pair = [['divergent', dsum], ['control_near_reference', ctrl]];
    const got = pair.filter(([, s]) => s && s.contigs_that_could_be_flagged);
    if (!got.length) {
      return '<div class="vsub">No contig in this run carries two or more named genes, ' +
        'so the chimaera check could not fire on either arm' + cite('chimaera') +
        '. Nothing here is dropped on it.</div>';
    }
    const say = got.map(([k, s]) =>
      pct(s.fraction_flagged_of_the_contigs_that_could_be) + ' of the ' +
      s.contigs_that_could_be_flagged +
      (k === 'divergent' ? ' divergent' : ' near-reference control') +
      ' contigs that could be').join(' and ');
    /* Whether the two arms read alike is a property of the run, not a sentence to
       type: on a library where the flag did separate them, the clause below would be
       asserting the opposite of its own numbers. It is read off the gap, and either way
       nothing is dropped on the flag, which is the part that does not depend on it. */
    const gap = (dsum && ctrl
                 && dsum.fraction_flagged_of_the_contigs_that_could_be != null
                 && ctrl.fraction_flagged_of_the_contigs_that_could_be != null)
      ? Math.abs(dsum.fraction_flagged_of_the_contigs_that_could_be
                 - ctrl.fraction_flagged_of_the_contigs_that_could_be)
      : null;
    const reading = gap == null
      ? ' Only one arm could be tested here, so the flag has no comparison on this run.'
      : (gap <= 0.2
         ? ' The two arms read within ' + points(gap) + ' of each other, so the flag ' +
           'tracks how the assembler joined genes rather than how far they sit from a ' +
           'reference.'
         : ' The two arms differ by ' + points(gap) + ', which this page reports without ' +
           'reading it as divergence: the flag has not been calibrated against an ' +
           'assembly whose joins are known.');
    return '<div class="vsub">Genes on one contig disagree about which organism they ' +
      'come from on ' + say + cite('chimaera') + '.' + reading +
      ' Nothing is dropped on it.</div>';
  }

  const hero = el('div', 'verdict');
  hero.innerHTML = '<div class="vmain"><b>' + belowNinety.length + ' of ' + taxa.length +
    '</b> assembled viruses sit <b>below 90% identity</b> to their closest published ' +
    'protein' + cite('contigid') +
    '.</div>' +
    // The two counts in this section come from two piles, and the page used to print the
    // first as a sentence and the second as a bare "15 of them" under the figure, leaving
    // a reader to guess at the relation. They are stated side by side and NOT as a
    // subset: the gene-level workup selects its own contigs, so 15 of 34 would be a claim
    // about an overlap this page has not measured.
    '<div class="vmain"><b>' + rows.length + '</b> contigs are scored gene by gene, ' +
    'which is what the figure below lays out' +
    (deep.length ? ', and <b>' + deep.length + '</b> of those are complete genomes ' +
      'carried at over 100x read depth' + cite('depth') : '') + '.</div>' +
    '<div class="vsub">Identity here means one gene against the single closest protein ' +
    'in the reference sets, at E&lt;1e-5 over at least 60 aa. Below about 90% that ' +
    'protein is a nearest neighbour rather than an identification, and the name it ' +
    'lends is borrowed.</div>' + chimaeraLine();
  p.append(hero);

  /* THE LINE. One quantity, one chart, no legend: everything that would need one is
     written on the drawing itself. The three bar blocks this replaced used the same
     orange bar for read depth, then distance, then identity, on one screen. */
  /* The headline's own controls, sitting with the headline: which gene the row is ranked
     on, and the two floors that decide whether a contig is in the characterisable set at
     all. Nothing is recomputed - both identities and the depth were measured - and the
     line under the rug says what moved. */
  const HIDS = ['spec_rank', 'min_contig_nt', 'min_depth'];
  // The pile the figure draws FROM, which is not the pile the headline counts: the
  // sentence above says how many assembled viruses there are, and the figure marks the
  // ones scored gene by gene that clear the shipped length and depth floors. Those were
  // 34 and 15 on the third library through this path, printed as "34" and "15 of them"
  // with nothing saying what the 15 were 15 of.
  const rugTotal = rows.length;
  const rugHost = el('div');
  const rugNote = el('div', 'note');
  // `shipped` reads the defaults rather than the current state, because the comparison
  // the note makes is against what this page ships and not against where the reader
  // happened to leave a slider, or against whatever a cited link set.
  function rugRows(shipped) {
    const g = shipped ? ctlShipped : ctlGet;
    const whole = g('spec_rank') === 'whole';
    const minnt = g('min_contig_nt'), mind = g('min_depth');
    return rows.map(r => Object.assign({}, r, {
      rank_pident: whole ? r.whole_contig_pident : r.rank_pident,
      ranked_on: whole ? 'the whole-contig best hit' : r.ranked_on
    })).filter(r => r.rank_pident != null &&
      (r.contig_nt == null || r.contig_nt >= minnt) &&
      ((r.read_support || {}).mapped_mean_depth == null ||
       (r.read_support || {}).mapped_mean_depth >= mind));
  }
  if (rows.length) {
    p.append(ctlRow(HIDS));
    p.append(rugHost);
    p.append(rugNote);
    p.append(settingLine(HIDS));
  }
  function drawRug() {
    const rows = rugRows();
    const deep = rows.filter(r => (r.read_support || {}).mapped_mean_depth >= 100
                                  && r.contig_nt >= 3000);
    if (!rows.length) {
      rugHost.innerHTML = '';
      rugNote.innerHTML = 'No contig is left at this setting, so there is no figure to ' +
        'draw rather than an empty axis to read.';
      return;
    }
    const W2 = 1000, L = 20, R = 18, TOP = 54, BOT = 46, H = TOP + 96 + BOT;
    const x = v => L + 36 + (Math.min(100, Math.max(20, v)) - 20) / 80 *
      (W2 - L - R - 36);
    const BASE = TOP + 58;
    let s = '<svg viewBox="0 0 ' + W2 + ' ' + H + '" width="100%" ' +
      'style="min-width:660px" role="img" aria-label="assembled viruses by identity ' +
      'to the closest published protein">';
    // the rug: one tick per contig, stacked where they pile up
    const bins = {};
    rows.forEach(r => {
      const b = Math.round(r.rank_pident * 2) / 2;
      bins[b] = (bins[b] || 0) + 1;
    });
    Object.keys(bins).forEach(k => {
      const n = bins[k], h = Math.min(46, 6 + n * 4);
      s += '<rect x="' + (x(+k) - 2.5) + '" y="' + (BASE - h) + '" width="5" height="' +
        h + '" rx="1.5" fill="var(--identified)" opacity="0.8"/>';
    });
    s += '<line x1="' + x(20) + '" y1="' + BASE + '" x2="' + x(100) + '" y2="' + BASE +
      '" stroke="var(--fg)" stroke-width="1.5"/>';
    for (let v = 20; v <= 100; v += 10) {
      s += '<line x1="' + x(v) + '" y1="' + BASE + '" x2="' + x(v) + '" y2="' +
        (BASE + 5) + '" stroke="var(--dim)" stroke-width="1"/>';
      s += '<text x="' + x(v) + '" y="' + (BASE + 19) +
        '" text-anchor="middle" class="ax">' + v + '%</text>';
    }
    s += '<text x="' + x(20) + '" y="' + (BASE + 36) + '" class="ax">' +
      'nothing published is close enough to name it</text>';
    s += '<text x="' + x(100) + '" y="' + (BASE + 36) +
      '" text-anchor="end" class="ax">the same virus, already sequenced</text>';
    // the handful worth naming
    deep.slice().sort((a, b) => a.rank_pident - b.rank_pident).forEach((r, i) => {
      const xx = x(r.rank_pident), yy = BASE - 54 - (i % 2) * 16;
      s += '<line x1="' + xx + '" y1="' + yy + '" x2="' + xx + '" y2="' + (BASE - 4) +
        '" stroke="var(--residue)" stroke-width="1" stroke-dasharray="2 3"/>';
      s += '<circle cx="' + xx + '" cy="' + yy + '" r="4" fill="var(--residue)"><title>' +
        esc(r.taxon) + '</title></circle>';
    });
    if (deep.length) {
      s += '<text x="' + (x(Math.min.apply(null, deep.map(r => r.rank_pident))) - 10) +
        '" y="' + (BASE - 50) + '" text-anchor="end" class="ax">' + deep.length +
        ' finished genomes, deep enough to work on</text>';
    }
    s += '</svg>';
    rugHost.innerHTML = '';
    const host = el('div', 'scroll'); host.innerHTML = s; rugHost.append(host);
    rugNote.innerHTML = 'One mark per contig, ' + rows.length +
      (rows.length === rugTotal
        ? ' of them'
        : ' of the ' + rugTotal + ' scored gene by gene; the rest sit below the length ' +
          'or depth floor on the control row above') +
      ', stacked where several land on the same identity. Position is the gene the row ' +
      'was ranked on, the polymerase wherever one is called.' + cite('genespread') +
      (ctlAnyMoved(HIDS)
        ? ' <span class="relabel">At your setting: ' + rows.length + ' contigs and ' +
          deep.length + ' finished genomes, against ' + rugRowsShipped.n + ' and ' +
          rugRowsShipped.deep + ' as shipped. The sentence above the figure is the ' +
          'shipped reading and does not move, because it is what the rest of this page ' +
          'is built from.</span>'
        : '');
  }
  const rugRowsShipped = (function () {
    const r = rugRows(true);
    return {n: r.length, deep: r.filter(x => (x.read_support || {}).mapped_mean_depth
      >= 100 && x.contig_nt >= 3000).length};
  })();
  if (rows.length) { drawRug(); ctlOn(drawRug); }

  /* THE FEW THAT MATTER, DRAWN AS WHAT THEY ARE. A bar chart of identities beside a Latin
     name tells a reader neither what the thing is nor what the number is a property of.
     The molecule itself is the picture: genes at their real coordinates, coloured by the
     job they do, with the identity written on the gene it belongs to. */
  const ROLE = GENE_ROLE, roleOf = geneRoleOf;

  if (deep.length) {
    const d = el('div', 'bites');
    d.append(el('div', 'pghead', 'the finished genomes, drawn to scale'));
    // What these are is read from the architecture file, which applied the committed
    // vocabulary's rule and recorded the evidence. The renderer counts rows.
    const calls = deep.map(r => archOf(r.contig)).filter(Boolean);
    const oneArch = calls.length === deep.length && calls.length > 0 &&
      calls.every(c => c.architecture && c.architecture === calls[0].architecture);
    const AW = oneArch ? archWords(calls[0]) : null;
    if (AW) {
      const full = calls.filter(c => c.has_every_expected_role).length;
      const intro = el('p', 'gintro');
      intro.innerHTML = 'All ' + deep.length + ' are called <b>' +
        esc(deep.length === 1 ? AW.label : (AW.label_plural || AW.label)) +
        '</b>: ' + esc(AW.plain) + ', and ' + esc(AW.does_not) + '. ' +
        (full ? full + ' of them carry ' + esc(AW.expected_roles_plain) + '. ' : '') +
        'Each is drawn below at its real gene positions; the arrow is the direction the ' +
        'gene is read and the number on it is identity to the closest published ' +
        'protein.' + cite('vocab');
      d.append(intro);
    }
    d.append(geneLegend());

    deep.slice().sort((a, b) => (b.read_support || {}).mapped_mean_depth -
                                (a.read_support || {}).mapped_mean_depth).forEach(r => {
      const rs = r.read_support || {}, os = (r.orfs_drawn || []).slice()
        .sort((a, b) => a.s - b.s);
      const card = el('div', 'gmap');
      const kb = (r.contig_nt / 1000).toFixed(1);
      const hd = el('div', 'gmaphead');
      hd.innerHTML = '<b>' + esc(r.taxon) + '</b> <span class="dim">&mdash; the name of ' +
        'its nearest published relative, not an identification</span>';
      card.append(hd);

      card.append(geneMap(r.contig_nt, os.map(o => Object.assign({}, o,
        {pident: o.sig ? o.pident : null})), kb + ' kb'));

      // what it is, in words a reader can check against the picture above. Every
      // claim here is a field of the architecture row or the vocabulary it applied.
      const A = archOf(r.contig), AX = archWords(A);
      const say = el('div', 'gsay');
      say.innerHTML =
        (AX && !AW ? '<b>' + esc(AX.label) + '.</b> ' : '') +
        (A && !A.architecture
          ? '<span class="dim">No architecture called: ' + esc(A.not_called_because) +
            '</span> ' : '') +
        '<b>' + kb + ' kb</b> at <b>' + Math.round(rs.mapped_mean_depth) +
        'x</b>, covered ' + (rs.breadth_covered > 0.98 ? 'end to end'
          : pct(rs.breadth_covered) + ' of its length') +
        (rs.evenness_within_2x_of_median < 0.6
          ? ' but unevenly, so part of it may be a different molecule' : '') + '. ' +
        (A && A.most_diverged_role && A.identity_spread_points != null
          ? 'Its <b>' + esc(A.most_diverged_role) + '</b> is the gene furthest from ' +
            'anything published and its <b>' + esc(A.least_diverged_role) +
            '</b> the closest, ' + num(A.identity_spread_points, 1) + ' points apart.'
          : (far ? 'Furthest from anything published: the <b>' +
             esc(far.role || 'unnamed gene') + '</b> at ' + num(far.pident, 1) + '%.'
             : '')) + ' ' +
        (A && A.genes_called > A.genes_with_a_role
          ? (A.genes_called - A.genes_with_a_role) + ' of its ' + A.genes_called +
            ' genes carry no role. ' : '') +
        (A && A.nearest_subject_provenance === 'environmental_record'
          ? '<span class="dim">Its nearest published relative is ' +
            esc(((VOC.subject_provenance || {}).environmental_record || {}).plain ||
                'an environmental record') + '.</span>' : '');
      card.append(say);
      d.append(card);
    });
    p.append(d);
  }

  /* "CAN WE RANK BY SCARINESS." No, and the reason is measured rather than cautious: a
     rank over sequence is a rank over database proximity, and the host attached to a
     matched record is a property of that record. What can be shown is what the water's
     unplaceable reads point at, in reads rather than in database entries, with the
     confidence that the host label transfers at all. */
  const R = D.ranking;
  if (R && R.host_mix_reads && R.reads_total) {
    const GL = ((R.host_table || {}).groups) || {};
    const ORD = ['bacteriophage', 'unknown_host', 'plant', 'vertebrate_non_human',
                 'human_associated', 'invertebrate', 'fungal_or_protist'];
    const hc = R.host_confidence_summary || {};
    const seen = Object.keys(R.host_mix_reads).sort(
      (a, b) => R.host_mix_reads[b] - R.host_mix_reads[a]);
    const tot = R.reads_total;
    const hum = R.host_mix_reads.human_associated || 0;
    const COL = {bacteriophage: 'var(--v1)', unknown_host: 'var(--null)',
                 plant: 'var(--v3)', vertebrate_non_human: 'var(--tier2)',
                 human_associated: 'var(--v2)', invertebrate: 'var(--cv3)',
                 fungal_or_protist: 'var(--cv2)'};

    const blk = el('div', 'bites');
    blk.append(el('div', 'pghead', 'could any of it infect a person'));
    const lead = el('div', 'callout');
    lead.innerHTML =
      'Of the <b>' + tot.toLocaleString() + '</b> leftover reads that point at a named ' +
      'relative at all, <b>' + pct((R.host_mix_reads.bacteriophage || 0) / tot) +
      '</b> point at viruses of bacteria and <b>' + pct(hum / tot) + '</b> at anything ' +
      'human-associated' + cite('concern') + '. ' +
      // WHAT THOSE READS ARE IS NOT SOMETHING THIS PAGE KNOWS. The sentence read
      // "Those N reads are the enteric viruses sewage is expected to carry, matched
      // over about 34 amino acids": an identification the host group does not support,
      // an expectation about sewage this run did not measure, and an alignment length
      // measured on another library. At N = 0 it also read "Those 0 reads are the
      // enteric viruses", which is the sentence arguing with its own number.
      (hum === 0
        ? 'No leftover read points at a record whose host is human-associated. '
        : 'A host group is the host of the matched RECORD, so those ' + hum +
          ' reads are reads whose nearest named relative carries that label. ') +
      '<b>Host confidence never exceeds ' + num(hc.max, 2) + '</b> anywhere in the ' +
      Number(hc.subjects || 0).toLocaleString() + ' records, so no relative here is ' +
      'close enough for its host to transfer: these are labels on the database, not ' +
      'claims about the water.';
    blk.append(lead);

    // one bar, in reads, direct-labelled. A share of records would be a share of the
    // database; a share of reads is a share of the sample.
    const BW = 1000, BH = 100, BY = 28, BAR = 26;
    let bs = '<svg viewBox="0 0 ' + BW + ' ' + BH + '" width="100%" ' +
      'style="min-width:560px" role="img" aria-label="what the named leftovers point ' +
      'at, by read">';
    let cx = 0, lab = [];
    seen.forEach(k => {
      const w = (R.host_mix_reads[k] / tot) * BW;
      if (w <= 0) return;
      bs += '<rect x="' + cx + '" y="' + BY + '" width="' + Math.max(0.6, w - 2) +
        '" height="' + BAR + '" rx="2" fill="' + (COL[k] || 'var(--null)') +
        '"><title>' + esc((GL[k] || {}).label || k) + ': ' +
        R.host_mix_reads[k].toLocaleString() + ' reads, ' +
        (R.host_mix[k] || 0).toLocaleString() + ' records</title></rect>';
      lab.push([cx + w / 2, k, w]);
      cx += w;
    });
    // only the shares with room get an inline label; the rest are called out beneath
    lab.forEach(l => {
      if (l[2] < 110) return;
      bs += '<text x="' + l[0] + '" y="' + (BY + BAR + 16) +
        '" text-anchor="middle" class="ax">' + esc((GL[l[1]] || {}).label || l[1]) +
        ' ' + pct(R.host_mix_reads[l[1]] / tot) + '</text>';
    });
    const hx = lab.filter(l => l[1] === 'human_associated')[0];
    if (hx) {
      bs += '<line x1="' + hx[0] + '" y1="' + (BY - 2) + '" x2="' + hx[0] + '" y2="' +
        (BY - 14) + '" stroke="var(--v2)" stroke-width="1.5"/>';
      bs += '<text x="' + Math.min(BW - 4, hx[0] + 6) + '" y="' + (BY - 18) +
        '" text-anchor="end" class="ax">human-associated, ' + pct(hum / tot) +
        ' of named reads</text>';
    }
    bs += '</svg>';
    const bh = el('div', 'scroll'); bh.innerHTML = bs; blk.append(bh);

    const ha = (R.subjects || []).filter(x => x.host_group === 'human_associated')
      .sort((a, b) => (b.residue_reads || 0) - (a.residue_reads || 0)).slice(0, 6);
    if (ha.length) {
      blk.append(el('p', 'note', 'The human-associated records, largest first: ' +
        ha.map(x => esc(x.organism || x.subject) + ' (' + x.residue_reads + ' reads, ' +
          num(x.median_identity_percent, 1) + '% over ' + x.median_aligned_aa + ' aa)')
          .join('; ') + '. At that alignment length a phage protein and a bacterial one ' +
        'cannot be told apart either, so read these as the direction of a match rather ' +
        'than as a finding.'));
    }

    const no = el('div', 'callout');
    no.innerHTML = '<b>There is no scariness number here, and that is a result rather ' +
      'than a gap.</b> Ordering by distance from the database orders by the database' +
      cite('nulls') + ', and the host attached to a match is a property of the matched ' +
      'record. What can be ranked is how much of something is present, which is a ' +
      'property of the water, and the table below is ordered that way. What would earn ' +
      'the word is recurrence: the same sequence appearing at one plant, then climbing ' +
      'across dates. That is measured at read level across 512 libraries' +
      cite('crossplant') + ' and has not yet been run over these contigs.';
    blk.append(no);
    p.append(blk);
  }

  /* EVERYTHING ELSE, as one table behind a disclosure. */
  if (rows.length) {
    const det = el('details', 'more');
    const sm = document.createElement('summary');
    sm.textContent = 'all ' + rows.length + ' divergent contigs, most abundant first';
    det.append(sm);
    const sc = el('div', 'scroll');
    const tb = el('table');
    tb.innerHTML = '<thead><tr><th>nearest named relative</th><th>identity</th>' +
      '<th>ranked on</th><th>genes disagree by</th><th>length</th><th>depth</th>' +
      '</tr></thead>';
    const body = el('tbody');
    rows.slice().sort((a, b) => ((b.read_support || {}).mapped_mean_depth || 0) -
                                ((a.read_support || {}).mapped_mean_depth || 0))
      .forEach(r => {
      const tr = el('tr');
      tr.innerHTML = '<td>' + esc(r.taxon) + '</td><td class="mono">' +
        num(r.rank_pident, 1) + '%</td><td class="dim">' + esc(r.ranked_on) +
        '</td><td class="mono">' + (r.orf_identity_spread_points == null
          ? '<span class="dim">one gene</span>'
          : num(r.orf_identity_spread_points, 1) + ' pts') +
        '</td><td class="mono">' + (r.contig_nt / 1000).toFixed(1) + ' kb</td>' +
        '<td class="mono">' + Math.round((r.read_support || {}).mapped_mean_depth || 0) +
        'x</td>';
      tr.title = r.contig;
      body.append(tr);
    });
    tb.append(body); sc.append(tb); det.append(sc); p.append(det);
  }

  /* The gene-level calibration, as one sentence instead of a fourth bar chart. */
  if (dsum && ctrl && dsum.by_gene_role && ctrl.by_gene_role) {
    const g = dsum.by_gene_role, c = ctrl.by_gene_role;
    const cal2 = el('div', 'callout');
    const say = k => g[k] ? k + ' ' + num(g[k].median_pident, 1) + '% against ' +
      (c[k] ? num(c[k].median_pident, 1) + '%' : 'no control') : null;
    cal2.innerHTML =
      'Which gene you measure decides the answer. On these contigs the genes of one ' +
      'contig disagree by a median of <b>' +
      num(dsum.median_orf_identity_spread_points, 1) + ' points</b>, against <b>' +
      num(ctrl.median_orf_identity_spread_points, 2) + '</b> on contigs that match ' +
      'their reference, so a single identity per contig reports whichever gene won on ' +
      'bitscore.' + cite('genespread') + ' By role, divergent against near-reference: ' +
      ['polymerase', 'maturation', 'coat'].map(say).filter(Boolean).join(', ') + '.';
    p.append(cal2);
  }

  const bs = el('details', 'more');
  const bsum = document.createElement('summary');
  bsum.textContent = 'limitations of this measurement';
  bs.append(bsum);
  const bl = document.createElement('ul');
  bl.className = 'note';
  ((S.blind_spots || []).concat((W && W.blind_spots) || [])).forEach(x => {
    const li = document.createElement('li');
    li.textContent = String(x).replace(/\\.$/, '') + '.';
    bl.append(li);
  });
  bs.append(bl);
  p.append(bs);
  app.append(p);
})();

// ---------- public health: the part nobody can name
(function () {
  const R = D.ranking, K = R && R.kind_summary, b = (D.tier_table || {}).breakdown;
  const cur = (D.tier_table || {}).curve || [];
  const tier = k => (cur.filter(c => c.tier === k)[0] || {});
  /* ONE FORMAT, EVERY SAMPLE: this section used to disappear on a library whose
     ranking or cascade breakdown was missing, which made a page with no remainder
     section and a page with an empty one look the same. It stays, and says which
     file it needs. */
  const p = panel('The unnamed remainder',
    'The fraction no reference names.' + cite('cascade'));
  if (!K || !b) {
    p.append(el('div', 'absent', 'not measured in this run: ' +
      (!K ? 'ranking.json carries no kind summary' : 'tier_table.json carries no ' +
       'breakdown') + ', so the share of the library no reference names cannot be ' +
      'quoted here'));
    app.append(p); return;
  }
  const own = tier('selfasm_sample_contigs_identified_norrna').marginal_removed;
  const rr = tier('rrna').tier_alone_removed;

  /* Each tier in the recognised-but-unnamed band with its own marginal share. The
     ribosomal index itself has its own band above, so it is excluded here rather than
     counted twice. A tier `tierKind` does not recognise is listed under its own name, so
     a new index appears as itself instead of swelling a label written for another one. */
  function unnamedBandParts(curve, rrAlone) {
    const NAMES = {artefact: 'adapter and low complexity',
                   selfasm_sample_contigs_rrna: 'own-sample ribosomal contigs',
                   selfasm_sample_contigs_identified_norrna:
                     'own-sample contigs carrying no name'};
    const parts = (curve || [])
      .filter(c => tierKind(c.tier) === 'recognises' && c.tier !== 'rrna'
                   && (c.marginal_removed || 0) > 0)
      .map(c => [NAMES[c.tier] || c.tier, c.marginal_removed])
      .sort((a, b) => b[1] - a[1]);
    if (!parts.length) {
      return 'no tier in this run removed reads without naming them';
    }
    return parts.map(x => x[0] + ' ' + pct(x[1])).join(', ');
  }

  const hero = el('div', 'verdict');
  hero.innerHTML = '<div class="vmain"><b>' +
    pct(1 - (K.named / (K.searched || 1))) + '</b> of residual reads have no ' +
    'significant protein match in any index available here.' + cite('residue') +
    '</div><div class="vsub">' + pct(own) + ' of the library matches contigs assembled ' +
    'from this sample that carry no name; ' + pct(rr) + ' is ribosomal RNA.' +
    cite('cascade') + '</div>';
  p.append(hero);

  // One picture of the whole library, split by WHAT WE CAN SAY about each part rather
  // than by which tier fired. The two middle bands are the honest answer to "removed
  // because there is an unknown contig": that sequence is recognised and still
  // nameless, which is a different thing from both named and left over.
  (function () {
    const brk = (D.tier_table || {}).breakdown, curve = (D.tier_table || {}).curve || [];
    if (!brk || !curve.length) return;
    let named = 0;
    curve.forEach(c => { if (tierCanName(c.tier)) named += c.marginal_removed || 0; });
    const ident = brk.identified.point;
    const bands = [
      {lab: 'named', v: named, col: 'var(--v1)',
       sub: 'a taxon or an accession'},
      /* THE BANDS MUST PARTITION THE IDENTIFIED FRACTION, SO EVERY ONE IS A MARGINAL.
         This band used the ribosomal index's own reach (76.06%) while the band beside it
         was computed as identified minus named minus this, so the stacked parts summed
         to 0.24 points more than the band they sat under. The index's own reach is a
         different quantity, it is the one the sentence above the figure quotes, and it
         belongs there and not in a stack. */
      {lab: 'ribosomal RNA', v: tier('rrna').marginal_removed, col: 'var(--v1)',
       op: 0.45, sub: 'identifies a ribosome, not an organism'},
      /* THE SUB-LABEL IS A DECOMPOSITION, NOT A GUESS AT ONE. It read "mostly
         own-sample contigs (2.26%), rest adapter and low complexity" over a band of
         14.18%: the 2.26% named only the non-ribosomal self-assembly tier and left the
         9.73% of own-sample RIBOSOMAL contigs out, while the artefact filters it
         attributed the remainder to removed 0.61%. Every tier in the band is listed with
         its own marginal share, largest first, in the shape the reconciliation panel
         below already used. */
      {lab: 'recognised, unnamed',
       v: Math.max(0, ident - named - tier('rrna').marginal_removed),
       col: 'var(--v-mute)', sub: unnamedBandParts(curve, rr)},
      /* THE SUB-LABELS SAY WHAT THE PILE IS, NOT WHAT IT MEASURED SOMEWHERE ELSE.
         These two read "median depth 6" and "median depth 1", which are one library's
         measurements written into the renderer. No stage in a run produces a median
         depth per pile, so on any other sample they are a number with no provenance -
         and on a small run the page contradicted itself on one screen, printing depth 6
         two lines above a sentence that said depth 1. What each pile IS can be said from
         its definition, and needs no measurement. */
      {lab: 'assembles, unidentified', v: brk.assembled_unidentified.point,
       col: 'var(--v2-soft)', sub: 'assembled here and nothing could name it'},
      {lab: 'residual', v: brk.residual.point, col: 'var(--v2)',
       sub: 'no name, and no contig here either'},
    ].filter(b => b.v > 0);
    const tot = bands.reduce((a, b) => a + b.v, 0) || 1;
    /* THE INVARIANT, ASSERTED WHERE THE BANDS ARE BUILT. These five are the whole
       library split by what can be said about each part, so they sum to one. They did
       not: the named band was short by every Kraken2 tier because a literal map of
       undated directory names had stopped matching, and the ribosomal band carried an
       index's own reach where the band beside it had been computed from a marginal. Both
       faults are invisible in a stacked bar, which normalises by its own total and so
       draws any set of numbers as a full width. A throw here fails the build, which is
       the only tell a figure drawn from a wrong total can have. */
    if (Math.abs(tot - 1) > 0.005) {
      throw new Error('the library bands sum to ' + tot.toFixed(4) +
        ' rather than 1: ' + bands.map(b => b.lab + ' ' + b.v.toFixed(4)).join(', '));
    }
    const bar = el('div', 'split');
    bar.innerHTML = bands.map((b, i) =>
      '<i title="' + esc(b.lab + ': ' + pct(b.v)) + '" style="flex:' +
      (b.v / tot).toFixed(5) + ';background:' + b.col +
      (b.op ? ';opacity:' + b.op : '') + ';animation-delay:' + (i * 110) +
      'ms"></i>').join('');
    p.append(bar);
    const lg = el('div', 'splitlg');
    lg.innerHTML = bands.map(b =>
      '<div><i style="background:' + b.col + (b.op ? ';opacity:' + b.op : '') +
      '"></i><div><b>' + esc(b.lab) + '</b> <span class="mono">' + pct(b.v) +
      '</span><div class="dim">' + esc(b.sub) + '</div></div></div>').join('');
    p.append(lg);
  })();

  const note = el('div', 'callout');
  /* The two sentences that remain are about what a SEARCH can establish, which is a
     property of the method and holds on every sample. The three that were here with them
     - artefact under 1%, no abundant unknown genome, median depth 1 - are measurements of
     one library, and this run produces none of them: no stage writes a per-pile median
     depth, and the artefact fraction of the residue is not computed anywhere. Printed
     anyway they were a claim about somebody's water taken from another sample, which is
     the fault the vocabulary file exists to prevent, applied to the numbers and not yet
     to the claims beside them. */
  note.innerHTML = '<b>Unnamed is not novel.</b> A search establishes that sequence ' +
    'is already known. It cannot establish novelty, because any search-derived score ' +
    'tracks distance to the database.' + cite('residue');
  p.append(note);
  app.append(p);
})();

// ---------- detail: naming against recognising
// NOT on the front page. A public-health reader does not care which tier fired, and
// "removed because it matched an unknown contig" is a fact about our plumbing. It stays
// in the page because it is the reconciliation between two figures that otherwise look
// like a contradiction, and because somebody asked the question once and will again.
(function () {
  const P = D.profile, curve = (D.tier_table || {}).curve || [];
  const brk = (D.tier_table || {}).breakdown;
  if (!P || !curve.length || !brk) return;
  // The subtitle quotes the two figures this panel reconciles, so they are computed
  // before it is written rather than typed into it. They were "7%" and "77%", which is
  // Berlin's answer printed on whatever library the page was built from.
  let named = 0, unnamed = 0; const unnamedRows = [];
  // Which tiers can name is a property of the tier, declared here where it can be read.
  // A Kraken2 or DIAMOND hit carries a taxon or an accession; an artefact filter, a
  // SILVA match and a match to a contig from this very sample tell you the sequence has
  // been seen before and nothing about what it is.
  const WHY = {artefact: 'adapter and low complexity',
               rrna: 'ribosomal RNA',
               selfasm_sample_contigs_rrna: 'own-sample ribosomal contigs',
               selfasm_sample_contigs_identified_norrna: 'own-sample unnamed contigs'};
  curve.forEach(c => {
    const m = c.marginal_removed || 0;
    if (tierCanName(c.tier)) { named += m; }
    else { unnamed += m; if (m > 0.002) unnamedRows.push([c.tier, m]); }
  });
  const k2here = 1 - (((P.at_confidence['0.1'] || {}).by_superkingdom || [])
    .filter(g => g.group === 'unclassified')[0] || {}).fraction;
  const p = panel('Named against recognised: cascade reconciliation',
    'Reconciling ' + pct(named) + ' named against ' + pct(brk.identified.point) +
    ' identifiable.');
  const box = el('div', 'cmp');
  const rowbar = (lab, v, col, sub) => {
    const r = el('div', 'crow');
    r.innerHTML = '<div class="clab">' + lab + '</div>' +
      '<div class="cbar"><span style="width:' + (100 * v).toFixed(2) + '%;background:' +
      col + '"></span></div><div class="cnum mono">' + pct(v) + '</div>' +
      '<div class="csub dim">' + sub + '</div>';
    box.append(r);
  };
  rowbar('<b>named</b>, one index', k2here, 'var(--v1)',
    'Kraken2 PlusPF-08, confidence 0.1, whole library.');
  /* Whether the two named rows agree is a comparison, not a sentence: the claim
     survived the tier names changing underneath it and went on asserting agreement
     across a fifteen-point gap. The two differ in scope as well as index, and the note
     under the panel says so, so the row states the gap and leaves the reading there. */
  rowbar('<b>named</b>, whole cascade', named, 'var(--v1)',
    'sum of every cascade tier returning a taxon or accession, ' +
    ((Math.abs(named - k2here) <= 0.02)
      ? 'within ' + ((Math.abs(named - k2here)) * 100).toFixed(1) +
        ' percentage points of the row above.'
      : ((named - k2here) * 100).toFixed(1) +
        ' percentage points from the row above, which is one index over the whole ' +
        'library against the cascade over its block.'));
  rowbar('<b>recognised, unnamed</b>', unnamed, 'var(--v-mute)',
    unnamedRows.map(r => (WHY[r[0]] || r[0]) + ', ' + pct(r[1])).join('; ') + '.');
  rowbar('residual', brk.residual.point, 'var(--v2)',
    'plus ' + pct((brk.assembled_unidentified || {}).point) +
    ' assembling into unidentified contigs.');
  p.append(box);
  p.append(el('p', 'note',
    'Most removal is by match to unnamed sequence. Scopes differ: the index row is ' +
    'the whole library, the cascade rows its ' +
    ((D.parameters || {}).pairs || 0).toLocaleString() + '-pair block.'));
  app.append(p);
})();
// ---------- methods
// Built after every panel, so a citation can only point at an entry that exists. The
// entries are static text: what varies between runs is in the result files and on the
// panels above, and a method that changed with the run would not be a method.
(function () {
  const p = panel('Methods',
    'One entry per headline figure, cited from it by number. Each states the input, ' +
    'the procedure, the threshold, the denominator, the control arm and the direction ' +
    'of the bias the control leaves behind.');
  METHODS.forEach((m, i) => {
    const d = el('div', 'meth');
    d.id = 'm-' + m.id;
    d.innerHTML = '<div class="mno">' + (i + 1) + '</div><div><h3>' + esc(m.t) +
      '</h3><dl class="mgrid">' +
      m.rows.map(r => '<dt>' + esc(r[0]) + '</dt><dd>' + r[1] + '</dd>').join('') +
      '</dl></div>';
    p.append(d);
  });
  app.append(p);
})();


/* NO GROUND-TRUTH SECTION, ON ANY SAMPLE. There used to be one here: the withheld
   genome of results/holdout_run/, first as its own panel on the one library that had
   one, then - worse - as a section on every page that other samples would render
   empty. Both are the fixture leaking into the tool. Imo, 2026-11-23: "WHY WOULD WE,
   WITH UNKNOWN WASTEWATER, HAVE A GENOME WITH A DISTANCE FIXED BEFORE THE RUN? BY
   DEFAULT WE KNOW NOTHING. WE ALWAYS KNOW NOTHING." A genome whose distance was
   settled before the run validates this tool; it is not something the tool found in
   anybody's water, so it is not in the report. That measurement lives in the
   validation suite - tests/ and results/holdout_run/ - and SPEC 13.4 says so. This
   report contains what was measured on this water and nothing about what we happened
   to know in advance. */

// ---------- summary and detail
// Runs last, after every panel has rendered. It reorganises rather than regenerates: the
// detail panels keep the markup they were built with and are only moved inside a
// disclosure, so a figure can never disagree with the summary that cites it.
// ---------- the significance floor, both arms in one pass
/* THE CONTROL THIS PAGE EXISTS TO MAKE VISIBLE. The reading this project has most had
   to withdraw was a hit rate taken with no E-value floor at all: 92.3% of shuffled reads
   "hit" something, and with a floor it is 0.00%. So the floor is a control, it is always
   stated, and the identified arm of the same run moves with it through the same
   predicate in the same pass. The shuffled arm is absent at read level in this run and
   the panel says so rather than drawing two arms as though they were three. */
(function () {
  const SR = D.search;
  const p = panel('Significance floor, both arms at your setting',
    'The fraction of reads with a relative at the floor you set, for the residue and ' +
    'for material the cascade already identified, counted in the same pass. ' +
    cite('controls'));
  if (!SR || !SR.piles || !(SR.piles.residue || []).length) {
    p.append(el('div', 'absent', 'reads.jsonl holds no scored pairs in this run, so ' +
      'there is nothing to re-filter'));
    app.append(p); return;
  }
  const IDS = ['evalue_floor', 'search_unit'];
  p.append(ctlRow(IDS));
  const fig = el('div'); p.append(fig);
  const foot = el('div', 'note'); p.append(foot);
  p.append(settingLine(IDS));
  p.append(nullLine(IDS, {
    shipped: 'absent at every setting, including this one: this run has no shuffled ' +
      'search at read level. The composition fields carry their shuffled null; the ' +
      'search does not, and a panel that drew one anyway would be inventing it.',
    moved: 'absent at every setting, this one included, for the same reason: there is ' +
      'no shuffled search arm at read level in this run.'}));

  const ARMS = [['residue', 'the residue', 'var(--residual)'],
                ['known', 'already identified', 'var(--identified)']];
  function rate(rows, e, unit) {
    let hit = 0, n = 0;
    rows.forEach(r => {
      const m1 = r[0] != null && r[0] <= e;
      const m2 = r[2] != null && r[2] <= e;
      if (unit === 'mate') {
        // Per mate, and only mates that exist: mate 2 is 14 nt shorter after the
        // non-template prefix comes off, so pooling the two averages two biases that
        // point opposite ways.
        n += 2; hit += (m1 ? 1 : 0) + (m2 ? 1 : 0);
      } else { n += 1; hit += (m1 || m2) ? 1 : 0; }
    });
    return {hit: hit, n: n, frac: n ? hit / n : null};
  }
  function draw() {
    const e = ctlGet('evalue_floor');
    const unit = ctlGet('search_unit');
    fig.innerHTML = '';
    const W = 940, BARH = 30, H = 20 + ARMS.length * (BARH + 26);
    let svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" role="img" ' +
      'aria-label="fraction with a significant relative, both arms">';
    const out = {};
    ARMS.forEach((a, i) => {
      const r = rate(SR.piles[a[0]] || [], e, unit);
      out[a[0]] = r;
      const y = 14 + i * (BARH + 26);
      svg += '<text x="0" y="' + (y - 2) + '" class="ax">' + a[1] + '</text>';
      svg += '<rect x="0" y="' + y + '" width="' + W + '" height="' + BARH +
        '" fill="var(--v-faint)" rx="4"/>';
      svg += '<rect x="0" y="' + y + '" width="' + ((r.frac || 0) * W) + '" height="' +
        BARH + '" fill="' + a[2] + '" rx="4"/>';
      svg += '<text x="' + Math.min(W - 6, (r.frac || 0) * W + 8) + '" y="' +
        (y + BARH - 9) + '" class="ax" text-anchor="' +
        ((r.frac || 0) > 0.8 ? 'end' : 'start') + '">' + pct(r.frac) + ' of ' + r.n +
        ' ' + (unit === 'mate' ? 'mates' : 'pairs') + '</text>';
    });
    svg += '</svg>';
    fig.innerHTML = svg;
    const d = out.residue, k = out.known;
    const pts = (d.frac != null && k.frac != null)
      ? (100 * (k.frac - d.frac)).toFixed(2) : null;
    foot.innerHTML = '<b>How to read it.</b> The gap, not either bar: at this floor the ' +
      'identified arm is ' + (pts == null ? 'not comparable' : pts +
      ' points above the residue') + '. Both arms are the same run, the same floor and ' +
      'the same pass, which is the only way either number means anything. ' +
      '<b>The floor is the measurement.</b> Dropping it admits chance alignments - a ' +
      '107 nt query against a large protein database finds about 34 amino acids at ' +
      'chance - and those are drawn from whatever the database holds most of, so a ' +
      'loose floor shows the database rather than the water.' +
      (ctlMoved('search_unit') ? ' <span class="relabel">Counted per mate at your ' +
        'setting, per pair as shipped: every percentage of this library elsewhere on ' +
        'the page is of pairs.</span>' : '');
  }
  draw(); ctlOn(draw);
  app.append(p);
})();

// ---------- what cannot move, and why
/* A control that cannot move honestly is visible as a stated constraint rather than
   absent from the page. Left off, the page would read as though these were simply not
   thought about; each row is a measurement that says why the number cannot be re-taken
   at output time. */
(function () {
  const p = panel('Controls that cannot move, and why',
    'Everything a reader might reasonably expect to adjust here that cannot be adjusted ' +
    'honestly, with the reason in each case. ' + cite('controls'));
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>control</th><th>why it cannot move at output time</th>' +
    '</tr></thead>';
  const body = el('tbody');
  (CTLDEF.fixed || []).forEach(f => {
    const tr = el('tr');
    tr.innerHTML = '<td>' + esc(f.control) + '</td><td>' + esc(f.why) + '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  /* The controls that would move on a run that measured them, stated here rather than
     left off. This panel is also where they are guaranteed to be visible: the panel a
     control belongs to may itself be a census of a stage that was not run, and a page
     that simply omitted the control would read as though the setting did not exist. */
  const gone = (CTLDEF.controls || []).filter(c => c.absent);
  if (gone.length) {
    p.append(el('p', 'note', 'And ' + gone.length +
      (gone.length === 1 ? ' control that moves' : ' controls that move') +
      ' on a run which measured them, but not on this one. ' +
      'That is not the same thing as a control that cannot move honestly, so it is said ' +
      'separately:'));
    const box = el('div', 'ctls');
    gone.forEach(c => box.append(ctlWidget(c.id)));
    p.append(box);
  }
  app.append(p);
})();

(function () {
  // The front page answers a public-health reader's questions in the order they ask
  // them. Everything else on this page is how those answers were arrived at, and goes
  // below the fold: it was all on the front page before, and the front page was a
  // method report wearing a sample report's clothes.
  /* ONE REPORT FORMAT, EVERY SAMPLE. These sections are the report: the same titles in
     the same order on every library, whether or not the run measured what they hold. A
     section with nothing in it says so in its own place and stays there, because a page
     that grows and drops sections by sample is a page a reader cannot compare with
     another one - and because a missing section and a measured absence must not look
     the same. The rule this replaces moved an empty summary panel into Detail, which
     kept any one page tidy at the cost of every page having a different shape. */
  const SUMMARY = ['Watchlist screen', 'What this water is',
                   'How far the assembled viruses sit from anything published',
                   'The unnamed remainder', 'Methods'];
  const panels = Array.prototype.filter.call(app.children,
    n => n.tagName === 'SECTION' && n.classList.contains('panel'));
  const detail = panels.filter(n => {
    const h = n.querySelector('h2');
    if (!h) return false;
    return SUMMARY.indexOf(h.textContent) < 0;
  });
  if (!detail.length) return;
  const firstSentence = t => {
    const whole = String(t || '').replace(/\\s+/g, ' ').trim();
    const m = whole.match(/^.*?\\.(?=\\s|$)/);
    let out = m ? m[0] : whole;
    if (out.length > 130) out = out.slice(0, 130).replace(/\\s+\\S*$/, '') + '...';
    return out;
  };
  // The summary panels answer a reader's questions in a fixed order, and a panel that
  // happens to be built last should not therefore be read last. Reorder to SUMMARY.
  const keep = panels.filter(n => detail.indexOf(n) < 0);
  keep.slice().sort((a, b) => SUMMARY.indexOf(a.querySelector('h2').textContent) -
                              SUMMARY.indexOf(b.querySelector('h2').textContent))
    .forEach(n => app.append(n));

  const wrap = el('section', 'panel');
  wrap.append(el('h2', null, 'Detail'));
  wrap.append(el('p', 'note',
    'Every measurement the summary is read from, one section each, closed by default. ' +
    'The summary is drawn from these same files.'));
  detail.forEach(n => {
    const d = el('details', 'ditem');
    const h = n.querySelector('h2');
    const note = n.querySelector('p.note');
    const sum = document.createElement('summary');
    sum.innerHTML = '<b>' + esc(h.textContent) + '</b> <span class="dim">' +
      esc(note ? firstSentence(note.textContent) : '') + '</span>';
    h.remove();
    n.classList.remove('panel');
    n.classList.add('inner');
    d.append(sum);
    d.append(n);
    wrap.append(d);
  });
  app.append(wrap);
})();

// ---------- the playground strip
/* Carries no control of its own: the controls sit with the figures they move. This says
   what the page is currently set to, hands back a link to that reading, and returns
   everything to the shipped defaults in one click. */
(function () {
  const bar = el('div', 'play');
  const line = el('div', 'playline');
  const btns = el('div', 'playbtns');
  const reset = el('button', null, 'Return to the shipped defaults');
  reset.type = 'button'; reset.onclick = ctlReset;
  const copy = el('button', null, 'Copy a link to this reading');
  copy.type = 'button';
  copy.onclick = () => {
    const url = location.href.split('#')[0] + ctlURL();
    const done = () => { copy.textContent = 'Link copied';
      setTimeout(() => { copy.textContent = 'Copy a link to this reading'; }, 1600); };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(url).then(done, () => { line.append(
        el('span', 'dim', ' ' + esc(url))); });
    } else { line.append(el('span', 'dim', ' ' + esc(url))); }
  };
  btns.append(reset, copy);
  const paint = () => {
    const moved = Object.keys(CTLBY).filter(ctlMoved);
    const absent = (CTLDEF.controls || []).filter(c => c.absent).length;
    line.innerHTML = moved.length
      ? '<b>Moved from the shipped defaults:</b> ' + esc(ctlSetting(moved)) +
        '. Every figure these touch says so, and nothing has been recomputed.'
      : '<b>Shipped defaults.</b> Every control sits beside the figure it moves, and ' +
        'each one re-filters numbers this run already measured. ' +
        (absent ? absent + (absent === 1 ? ' control cannot move on this run, and says'
                                          : ' controls cannot move on this run, and say') +
          ' why. ' : '') +
        'Nothing here produces a novelty call, because there is no score to threshold.';
    if (CTLDROPPED.length) {
      line.append(el('div', 'dim', 'A setting in the link you opened is not one this ' +
        'run measured, so it was dropped rather than moved to the nearest point: ' +
        esc(CTLDROPPED.join(', '))));
    }
    reset.disabled = !moved.length;
  };
  paint(); ctlOn(paint);
  bar.append(line, btns);
  app.prepend(bar);
})();

// ---------- reveal
// Grow every split from nothing the first time it is scrolled into view. The width that
// is animated TO is the one the panel already wrote, so this cannot change a figure: if
// the observer never runs, or motion is reduced, the bars are simply already correct.
(function () {
  const reduce = window.matchMedia &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  if (reduce || !window.IntersectionObserver) return;
  const bars = Array.prototype.slice.call(
    document.querySelectorAll('.bar, .bbar, .axbar'));
  const stash = [];
  bars.forEach(bar => {
    const parts = Array.prototype.slice.call(bar.children).map(c => [c, c.style.width]);
    if (!parts.length) return;
    parts.forEach(([c]) => { c.style.width = '0%'; });
    stash.push([bar, parts]);
  });
  const io = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if (!e.isIntersecting) return;
      io.unobserve(e.target);
      const hit = stash.filter(x => x[0] === e.target)[0];
      if (!hit) return;
      hit[1].forEach(([c, w], i) => {
        window.setTimeout(() => { c.style.width = w || ''; }, 40 * i);
      });
    });
  }, {rootMargin: '0px 0px -8% 0px'});
  stash.forEach(([bar]) => io.observe(bar));
  // Belt and braces: whatever happens to the observer - a panel inside a closed
  // disclosure, a browser that never fires it, a thumbnail capture - every bar is back
  // at its true width a few seconds in. The animation is decoration; the width is data.
  window.setTimeout(() => {
    stash.forEach(([bar, parts]) => parts.forEach(([c, w]) => {
      if (!c.style.width || c.style.width === '0%') c.style.width = w || '';
    }));
  }, 2500);
})();
</script>
"""


#: Run in headless Chromium against the page this build is about to write. The bundle is
#: per-sample, so which panels render, and which citations they print, is only known once
#: the page runs - a static check over the source can say every cite() names an entry, and
#: cannot say a panel threw and silently took every later panel with it, which is how the
#: holdout page once rendered four panels out of twenty.
RENDER_CHECK = r"""
const {chromium} = require('playwright');
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage();
  const errors = [];
  p.on('pageerror', e => errors.push(String(e)));
  await p.goto('file://' + process.argv[2]);
  await p.waitForTimeout(500);
  const r = await p.evaluate(() => {
    window.D = D;
    const app = document.getElementById('app');
    const h = n => ((n.querySelector('h2') || {}).textContent || '');
    const top = [...app.children].filter(n => n.tagName === 'SECTION').map(h)
      .filter(t => t !== 'Detail');
    const detail = [...document.querySelectorAll('.ditem > summary b')].map(x => x.textContent);
    const cited = [...new Set([...document.querySelectorAll('a.cite')]
      .map(a => a.getAttribute('href').slice(3)))];
    const dangling = cited.filter(id => !document.getElementById('m-' + id));
    // Every control the generator resolved has to be ON the page: one that resolved and
    // then rendered nowhere is a setting the reader cannot see or return from, and one
    // that is absent has to be visible as a stated constraint rather than dropped. The
    // build fails on either, which is the rule that a control the generator cannot
    // verify fails the build rather than rendering.
    const shown = new Set([...document.querySelectorAll('[data-ctl]')]
      .map(n => n.getAttribute('data-ctl')));
    const absentShown = new Set([...document.querySelectorAll('.ctl-absent[data-ctl]')]
      .map(n => n.getAttribute('data-ctl')));
    const declared = ((window.D || {}).controls || {}).controls || [];
    const missing = declared.filter(c => !shown.has(c.id)).map(c => c.id);
    const absent_not_stated = declared.filter(c => c.absent && !absentShown.has(c.id))
      .map(c => c.id);
    // A reading has to carry its setting, so a controlled panel carries a setting line.
    const settings = document.querySelectorAll('.setting').length;
    return {top, detail, cited, dangling, missing, absent_not_stated, settings,
            controls: declared.length};
  });
  console.log(JSON.stringify({errors, ...r}));
  await b.close();
})().catch(e => { console.log(JSON.stringify({fatal: String(e)})); });
"""


def _render_check_dir() -> Path:
    """Where to put the check script, which decides whether it can find playwright.

    Node resolves a bare `require` from the SCRIPT's directory upward, not from the
    working directory. The script used to be written beside the page, which on any real
    run is a run directory under /tmp with no node_modules above it: `npm install
    playwright` in the checkout, which is where anyone would run it, then had no effect
    and the check reported itself as unavailable. So the script goes where the module is
    installed, and beside the page only when nothing in the checkout can hold it.
    """
    for candidate in (ROOT, Path.cwd()):
        if (candidate / "node_modules" / "playwright").exists():
            return candidate
    return ROOT if os.access(ROOT, os.W_OK) else Path.cwd()


def render_check(page_path: Path) -> dict | None:
    """Run the written page; None where no browser is available, and the build says so."""
    tmp = _render_check_dir() / ".render_check.js"
    try:
        tmp.write_text(RENDER_CHECK)
        r = subprocess.run(["node", str(tmp), str(page_path.resolve())],
                           capture_output=True, text=True, timeout=120)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return None
    finally:
        tmp.unlink(missing_ok=True)
    try:
        out = json.loads(r.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError):
        return None
    if "fatal" in out:
        if "Cannot find module" in out["fatal"]:
            return None
        raise SystemExit(f"the render check could not run the page: {out['fatal']}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build the run viewer for one sample from its page descriptor. Every "
                    "sample goes through this same command; what differs is the "
                    "descriptor, never a flag.")
    ap.add_argument("page", type=Path,
                    help="a page descriptor, results/web/pages/<sample>.json")
    ap.add_argument("--rundir", type=Path, default=None,
                    help="where the run directory is, if not at the descriptor's rundir "
                         "(a location, not a choice of content: the run must be the "
                         "descriptor's sample)")
    ap.add_argument("--out", type=Path, default=None,
                    help="the HTML file to write; default <rundir>/web/index.html")
    ap.add_argument("--scatter-per-pile", type=int, default=2500,
                    help="points drawn per pile in the joint space. The page states the "
                         "number drawn and the number available, because a sample that "
                         "does not say it is one reads as the whole pile")
    ap.add_argument("--table-rows", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1, help="draws those samples")
    args = ap.parse_args()

    t0 = time.perf_counter()
    page_desc = load_page(args.page, args.rundir)
    rundir = page_desc["_rundir"]
    b = bundle(rundir, args.scatter_per_pile, args.table_rows, args.seed, page_desc)
    b["built"]["page_descriptor"] = str(args.page)
    out = args.out or (rundir / "web" / "index.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(b, separators=(",", ":"))
    if "</script" in payload:
        raise SystemExit("the bundle contains a closing script tag and would break the "
                         "page; a value in the run directory needs escaping")
    # A citation that points at no method entry is worse than no citation: it reads as
    # a reference and goes nowhere. Two checks, because the entries are now data. This
    # one is static: every cite() in the source names an entry in the methods file. The
    # render check below is the one that knows which citations this sample's page prints.
    cited = set(re.findall(r"cite\('([a-z_]+)'\)", PAGE))
    defined = {m["id"] for m in b["methods"]}
    missing = sorted(cited - defined)
    if missing:
        raise SystemExit(f"citations with no method entry: {', '.join(missing)}")
    unused = sorted(defined - cited)

    # Two run directories make two pages, and an artifact gallery lists both by their
    # <title>. One title for both is two rows nobody can tell apart.
    page = PAGE.replace("<title>Hapax run viewer</title>",
                        f"<title>{page_desc['title']}</title>")
    page = page.replace("__BUNDLE__", payload)
    # A page whose script does not parse renders as a blank panel-less document and says
    # nothing about why, which is exactly how one stray apostrophe cost a build here. Node
    # is present in this container; where it is not, the check is skipped rather than
    # faked, and the build says which it did.
    script = page[page.rindex("<script>") + len("<script>"):page.rindex("</script>")]
    tmp = out.parent / ".page_check.js"
    checked = "not checked: node is not on PATH"
    try:
        tmp.write_text(script)
        r = subprocess.run(["node", "--check", str(tmp)], capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"the generated page does not parse:\n{r.stderr.strip()}")
        checked = "node --check passed"
    except FileNotFoundError:
        pass
    finally:
        tmp.unlink(missing_ok=True)
    out.write_text(page)
    size = out.stat().st_size

    rc = render_check(out)
    if rc is None:
        rendered = "render not checked: node or playwright is not available"
    else:
        bad = []
        if rc["errors"]:
            bad.append("the page threw while rendering, so every panel after the throw "
                       "is missing: " + "; ".join(rc["errors"]))
        if rc["dangling"]:
            bad.append("rendered citations with no Methods entry on this page: " +
                       ", ".join(rc["dangling"]))
        if rc.get("missing"):
            bad.append("controls the generator resolved that the page never renders, so "
                       "the reader cannot see or return from them: " +
                       ", ".join(rc["missing"]))
        if rc.get("absent_not_stated"):
            bad.append("controls that cannot move on this run and are not visible as a "
                       "stated constraint: " + ", ".join(rc["absent_not_stated"]))
        if rc.get("controls") and not rc.get("settings"):
            bad.append("the page has controls and no setting line, so a reading taken "
                       "here would not carry the setting it was taken at")
        if bad:
            out.unlink()
            raise SystemExit("render check failed, page not written:\n  - " +
                             "\n  - ".join(bad))
        rendered = (f"render check passed: {len(rc['top'])} summary panels, "
                    f"{len(rc['detail'])} in detail, {len(rc['cited'])} methods cited, "
                    f"{rc.get('controls', 0)} controls over "
                    f"{rc.get('settings', 0)} setting lines")
    print(f"{out}  {size / 1e6:.2f} MB  {checked}; {rendered}  "
          f"residue rows {b['table']['n_total']}, drawn {b['table']['n_drawn']}; "
          f"scatter {b['scatter'].get('residue', {}).get('n_drawn', 0)}+"
          f"{b['scatter'].get('known', {}).get('n_drawn', 0)}; "
          f"{len(cited)} of {len(defined)} methods citable"
          + (f", never cited: {', '.join(unused)}" if unused else ""))
    if rc is not None:
        print("summary: " + " | ".join(rc["top"]))
        print("not run: " + (", ".join(k for k, v in page_desc["_inputs"].items()
                                       if v is None) or "nothing"))
    if size > 16_000_000:
        raise SystemExit(f"{size} bytes is over the 16 MB page limit; lower "
                         f"--scatter-per-pile or --table-rows")
    print(f"built in {time.perf_counter() - t0:.1f}s "
          f"({platform.python_version()})", flush=True)


if __name__ == "__main__":
    main()
