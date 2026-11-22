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
import re
import platform
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

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
     lambda r: _log2(((r.get("library_abundance") or {}).get("median_31mer_count")))),
    ("longest_orf_aa", "Longest ORF (aa)",
     lambda r: (r.get("frame") or {}).get("longest_orf_aa")),
    ("CpG_oe", "CpG observed/expected",
     lambda r: (r.get("dinuc") or {}).get("CpG_oe")),
    ("UpA_oe", "UpA observed/expected",
     lambda r: (r.get("dinuc") or {}).get("UpA_oe")),
    ("fold_z", "Folding energy z",
     lambda r: (r.get("fold") or {}).get("z")),
    ("kmer16_duplicate_frac", "Repetitive fraction",
     lambda r: (r.get("repetitive") or {}).get("kmer16_duplicate_frac")),
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
    return {k: w[k] for k in ("summary", "operating_point", "the_number_to_rank_by",
                              "bands", "blind_spots", "seed", "git_sha",
                              "script_sha256") if k in w} | {"rows": rows}


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


def spectrum() -> dict | None:
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
    dk = load(ROOT / "results/divergent_known/divergent_known_s1.json")
    if dk is None:
        return None
    big = {}
    cn = ROOT / "results/divergent_contigs/core_nt_block_001.jsonl"
    if cn.exists():
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
            "lo": min(genes) if len(genes) > 1 else None,
            "hi": max(genes) if len(genes) > 1 else None,
            "genes": len(genes), "cov": cov, "aln": aln,
            "nt": r["contig_nt"], "evalue": r["evalue"],
            "depth": (r.get("read_support") or {}).get("mapped_median_depth"),
            "subject": (r.get("subject_title") or r.get("subject") or "")[:90],
            "big": big.get(r["contig"]),
        }
        (lane if why else on_axis).append({**row, "why": why} if why else row)

    # THE ONE POINT ON THE PLOT WHOSE TRUE DISTANCE IS KNOWN. Everything else on this
    # axis is a distance to the nearest thing a database happens to hold, which is an
    # upper bound on identity. The holdout genome was deposited after the reference cutoff
    # and withheld from every index, so its real divergence is a measured 30.28% at
    # nucleotide level rather than an inference, and it belongs on the plot labelled as
    # such. Its three genes also show the gene spread doing its job: 24.3 points.
    hv = load(ROOT / "results/conserved_positions/conserved_positions_s1.json") or {}
    hold = None
    for row in hv.get("per_contig_roles") or []:
        if row.get("arm") == "extra_candidate_withheld":
            hold = {"by_role": row["by_role"], "spread": row["spread_points"],
                    "true_nucleotide_divergence_percent": 30.28,
                    "note": ("a genome withheld from every index; its position here is a "
                             "measurement of this axis, not an estimate of the genome")}
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
        "holdout": hold,
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


def bundle(rundir: Path, scatter_per_pile: int, table_rows: int, seed: int) -> dict:
    rng = random.Random(seed)
    manifest = load(rundir / "manifest.json")
    if manifest is None:
        raise SystemExit(f"{rundir}/manifest.json is absent; this is not a hapax run "
                         f"directory")
    rows_path = rundir / "reads.jsonl"
    piles: dict[str, list[dict]] = {"residue": [], "known": []}
    if rows_path.exists():
        with open(rows_path) as fh:
            for line in fh:
                r = json.loads(line)
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
    lp = rundir / "layers.jsonl"
    if lp.exists():
        with open(lp) as fh:
            for line in fh:
                d = json.loads(line)
                layers[d["name"]] = d["layer"]

    strains_subjects = [r.get("subject") for r in
                        ((load(rundir / "strains.json") or {}).get("subjects") or [])]

    return {
        "run": {k: manifest.get(k) for k in
                ("tool", "spec", "seed", "git_sha", "tree_dirty", "script_sha256",
                 "inputs", "complete", "incomplete_because", "wall_seconds",
                 "field_calibration", "known_positive_rows", "known_positive_note",
                 "interface_rules", "stages")},
        "parameters": manifest.get("parameters", {}),
        "spectrum": spectrum(),
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
        "ranking": _trim_ranking(load(rundir / "ranking.json") or {}, seed),
        "axes": [{"key": k, "label": l} for k, l, _f in SCATTER_AXES],
        "scatter": scatter,
        "table": {"rows": [{"name": r.get("name"),
                            "layer": layers.get(r.get("name")),
                            "values": [f(r) for _k, _l, f in SCATTER_AXES],
                            "closest": {k: v for k, v in r.items()
                                        if k.startswith("closest_")
                                        and isinstance(v, dict)},
                            "motif": (r.get("motif") or [None])[0],
                            "chimaera": r.get("chimaera")}
                           for r in table],
                  "n_drawn": len(table), "n_total": len(res), "sample_seed": seed},
        "built": {"script": "experiments/make_web.py", **provenance(__file__),
                  "rundir": str(rundir)},
    }


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


/* ---- METHODS: ONE ENTRY PER HEADLINE FIGURE -------------------------------
   The page states numbers at the top and how they were obtained at the bottom, and
   the two are joined by a link rather than by the reader remembering. Each entry
   carries the six things that make a figure readable here: its input, the procedure,
   the threshold, the denominator, the control arm, and the direction of the bias that
   the control does not remove. An entry with no control arm says so in that row.
   The numbering is positional, so an entry cannot be cited by a number that moved. */
const METHODS = [
  {id: 'library', t: 'Library and read handling', rows: [
    ['Input', 'SRR23973193, Berlin influent, total RNA, 8,731,224 read pairs.'],
    ['Procedure', 'Mate 2 carries about 14 nt of non-template 5-prime sequence; 16 nt ' +
      'are trimmed, leaving mates of 107 and 93 nt. Every percentage on this page is ' +
      'of read pairs.'],
    ['Threshold', 'None.'],
    ['Denominator', 'Read pairs in the analysed block.'],
    ['Control', 'Per-base error rate 0.00068, measured on self-overlapping pairs, so ' +
      'quality-derived error rates run 2.6x conservative.'],
    ['Bias', 'One plant, one date. In this plant the ribosomal fraction moves from ' +
      '38.0% to 67.1% and the identifiable fraction from 22.7% to 40.7% over eight ' +
      'months. Every fraction here is specific to this run.'],
  ]},
  {id: 'watchlist', t: 'Watchlist screen', rows: [
    ['Input', '44 taxa, each with a recorded reason for inclusion.'],
    ['Procedure', 'Kraken2 against PlusPF-08. A pair counts for a taxon when it is ' +
      'assigned to that taxon or below it.'],
    ['Threshold', 'Confidence 0.1 for the ranking; 0.2 for the per-row chip. ' +
      'Threshold dominates database size by an order of magnitude: one database ' +
      'removes 42.8%, 21.7% or 3.5% of this library on its setting alone.'],
    ['Denominator', 'All read pairs.'],
    ['Control', 'The same panel over composition-matched shuffled reads, counted on ' +
      'the panel at both settings.'],
    ['Bias', 'Misassignment to a close relative inflates a detection and is ' +
      'uncontrolled. Non-detections carry rule-of-three upper bounds at this depth, ' +
      'which are statements about sequencing depth and not about the catchment.'],
  ]},
  {id: 'composition', t: 'Library composition', rows: [
    ['Input', 'All read pairs.'],
    ['Procedure', 'Kraken2 PlusPF-08, assignments rolled up to superkingdom. The ' +
      'second ring renormalises over the named fraction.'],
    ['Threshold', 'Confidence selectable on the panel; the named fraction spans the ' +
      'range quoted in the caveat.'],
    ['Denominator', 'First ring all pairs, second ring the named pairs.'],
    ['Control', 'None applicable: the figure is one classifier output.'],
    ['Bias', 'RNA abundance, not cell counts: one bacterium with a thousand ribosomes ' +
      'contributes a thousand times its cell count. The unnamed share is this index ' +
      'reach at this setting; other tiers recognise most of it.'],
  ]},
  {id: 'cascade', t: 'Cascade removal', rows: [
    ['Input', 'The read pairs. Each subtraction index runs independently over all of ' +
      'them, so per-tier removal is order-free.'],
    ['Procedure', 'Union across tiers. The ribosomal tier precedes the protein tiers, ' +
      'because against a protein index ribosomal reads stop looking novel first. ' +
      'Marginal removal is the increment a tier adds to the union.'],
    ['Threshold', 'Kraken2 confidence as above; DIAMOND bitscore 50 for the protein ' +
      'tier, swept across 30 to 100 in the output file.'],
    ['Denominator', 'Read pairs.'],
    ['Control', 'Shuffled reads, which score zero at every tier.'],
    ['Bias', 'A tier can match sequence that nothing names, so identifiable exceeds ' +
      'named: 77.10% identifiable against 7.49% nameable. Five seeds are contiguous ' +
      'blocks of one flow cell with overdispersion about 1.5x, so the reported spread ' +
      'is observed and is not an estimate of dispersion.'],
  ]},
  {id: 'contigid', t: 'Contig identification and divergence', rows: [
    ['Input', 'Contigs assembled from this library.'],
    ['Procedure', 'DIAMOND blastx against RefSeq viral proteins and Swiss-Prot. Every ' +
      'identity is reported with its aligned length and the fraction of the subject ' +
      'covered.'],
    ['Threshold', 'E below 1e-5 and at least 60 aa aligned. Identity is quantised at ' +
      'shorter alignment lengths and is not plotted there.'],
    ['Denominator', 'Contigs with a significant hit.'],
    ['Control', 'A dinucleotide-preserving shuffle of the same contigs through the ' +
      'same search. 89.5% of real contigs carry a significant alignment and 0 of ' +
      '11,361 shuffles do; the whole shuffled library yields one 47 aa alignment at ' +
      'E=8.4e-4, which clears neither floor. Existence of the alignment is therefore ' +
      'the calibrated quantity, and no field conditional on a hit has a null to be ' +
      'compared against. This control asks whether a hit is real. It does not ask ' +
      'whether identity measures the distance between a real contig and its real ' +
      'relative, which the holdout answers: a genome 30.28% diverged was reported at ' +
      '66.5% protein identity over 620 aa against a referenced comparator in the same ' +
      'library at 99.9%, ordered correctly by 30 points.'],
    ['Bias', 'Both reference sets are narrow, so each identity is an upper bound on ' +
      'identity and each divergence a lower bound on divergence. That bound is ' +
      'measured at nucleotide level against NCBI core_nt, 319 GB: of 149 divergent ' +
      'contigs, 25 collapse to 90% identity or above over at least half their length ' +
      'and 83 stay under 80% identity or under half covered, while all 135 ' +
      'near-reference controls return at a median 97.55% over 98.8% of their length. ' +
      'The protein equivalent against nr has not landed.'],
  ]},
  {id: 'genespread', t: 'Per-gene identity and ranking', rows: [
    ['Input', 'ORFs called by Prodigal in metagenome mode, minimum 50 aa.'],
    ['Procedure', 'DIAMOND blastp per ORF, with the gene role read from the subject ' +
      'title. Rows rank on the polymerase where one is called; the panel states when ' +
      'a row ranks on the whole contig instead.'],
    ['Threshold', 'E below 1e-5, at least 60 aa.'],
    ['Denominator', 'For the spread, contigs with two or more scorable genes.'],
    ['Control', 'Contigs matching their reference at 95% or above, same procedure. ' +
      'Their genes agree to a median 2.35 points where the divergent arm disagrees by ' +
      '17.4. Over the whole identified pile, 600 divergent against 302 near-reference ' +
      'contigs, the same comparison reads 12.8 against 2.45. Which gene carries the ' +
      'distance is read paired, role against role on contigs carrying both as separate ' +
      'ORFs: the maturation protein sits below the polymerase on 122 of 131 such ' +
      'contigs, a median 6.9 points, and below the capsid on 28 of 31. The ' +
      'near-reference arm is a coin toss there, which is what makes the divergent arm ' +
      'readable.'],
    ['Bias', 'The arms are not matched as drawn: divergent contigs carry a median of 3 ' +
      'scorable genes over 3,820 nt against 1 over 552 nt, and more genes give a ' +
      'spread more room, which biases the gap upward. At equal gene counts the gap ' +
      'holds: 5.6 against 1.3 points at two genes, 16.9 against 2.65 at three. One ' +
      'identity per contig is whichever gene won on bitscore, usually the more ' +
      'conserved gene. Protein identity also understates nucleotide divergence, since ' +
      'synonymous changes are invisible to it. Two confounds on the role ordering were ' +
      'measured and neither explains it: distinct reference proteins matched per ORF ' +
      'are 0.604 for the polymerase against 0.598 for the maturation protein, and ' +
      'inside matched alignment-length bands the gap is 7 to 8 points. The capsid is ' +
      'not placed in the ordering, because half its alignments run under 150 aa where ' +
      'identity is quantised. Two roles called on one ORF are one number compared with ' +
      'itself and are excluded.'],
  ]},
  {id: 'depth', t: 'Mapped depth and evenness', rows: [
    ['Input', 'Every read pair in the library, mate 2 trimmed as above.'],
    ['Procedure', 'Reads mapped back onto the contigs with minimap2. Depth is the ' +
      'mean over contig positions; evenness is the fraction of positions within 2x ' +
      'of the median.'],
    ['Threshold', 'Primary alignments only.'],
    ['Denominator', 'Contig positions.'],
    ['Control', 'The assembler coverage header on the same contigs, which understates ' +
      'mapped depth about fourfold on this set.'],
    ['Bias', 'A region shared with another molecule collects that molecule reads, ' +
      'which raises depth without raising evidence. Evenness is the field that ' +
      'exposes it.'],
  ]},
  {id: 'chimaera', t: 'Chimaera check', rows: [
    ['Input', 'The same contigs and the same per-gene hits.'],
    ['Procedure', 'Genes on one contig are tested for whether they point at the same ' +
      'subject organism.'],
    ['Threshold', 'As the gene level above.'],
    ['Denominator', 'Contigs with two or more scorable genes.'],
    ['Control', 'The near-reference arm. The flag fires on 46.31% of divergent ' +
      'contigs and 39.26% of controls.'],
    ['Bias', 'The flag tracks how the assembler joins genes, so it carries almost no ' +
      'information about divergence and nothing is dropped on it. Split-mapping, the ' +
      'read-level version, has a 26% false-positive rate on known background.'],
  ]},
  {id: 'residue', t: 'Residue census and recurrence', rows: [
    ['Input', 'Read pairs surviving every subtraction tier.'],
    ['Procedure', 'Artefact signatures, mapped depth, and occupancy across the other ' +
      'dates from the same plant.'],
    ['Threshold', 'Significance as above. Occupancy is counted over 41 dates, which ' +
      'are 85 distinct waters across 215 runs.'],
    ['Denominator', 'Residue pairs seen once, against identified pairs at the same ' +
      'depth.'],
    ['Control', 'Identified reads at matched depth: 21.5% of residue reads appear on ' +
      'none of the other dates against 0.9%, a factor of 24. Artefact signatures ' +
      'together account for under 1% of the pile.'],
    ['Bias', 'A query drawn from the library being searched is abundant there by ' +
      'construction, which made 54% of background reads look like they were growing ' +
      'against 6% under the swap control, so queries come from outside the library ' +
      'searched. The median depth of this pile is 1, which is the floor of what any ' +
      'read-level field can read.'],
  ]},
  {id: 'spectrum', t: 'The distance spectrum', rows: [
    ['Input', 'Every identified contig of at least 300 nt that is short enough for ' +
      'per-ORF scoring, from the divergent stage, both arms.'],
    ['Procedure', 'ORFs called by Prodigal in metagenome mode, DIAMOND blastx and blastp ' +
      'against RefSeq viral proteins and Swiss-Prot, roles read from Pfam at gathering ' +
      'cutoffs. Position on the axis is the identity of the polymerase where one is ' +
      'called and the whole-contig hit otherwise, stated per row. The bar is the range ' +
      'across that contig own scorable genes. Colour is the fraction of the scored gene ' +
      'the alignment spans. Height is median mapped depth from minimap2.'],
    ['Threshold', 'E below 1e-5 and at least 60 aa aligned for a hit to exist. A contig ' +
      'is drawn at an identity only if at least 50% of the scored gene is aligned over ' +
      'at least 100 aa; the rest go to a lane off the axis.'],
    ['Denominator', 'Contigs with a significant protein hit and a scorable gene. The ' +
      'panel prints that count against the whole searched pile.'],
    ['Control', 'Two, drawn by the same code. The near-reference arm, contigs matching ' +
      'their reference at 95% or above, lands at the right-hand end by construction and ' +
      'is what makes the left-hand end readable. And a dinucleotide-preserving shuffle ' +
      'of the same contigs, which produces 0 significant alignments from 11,361, so the ' +
      'existence of every point on this axis is calibrated even though no identity on ' +
      'it is.'],
    ['Bias', 'Identity against a reference orders contigs by database proximity, so no ' +
      'threshold is drawn on the axis and the three-bucket overlay is off by default. ' +
      'Both reference sets are narrow, so each position is an upper bound on identity; ' +
      'against core_nt, 319 GB, 25 of 149 divergent contigs collapsed to a strain-level ' +
      'identification and 83 stayed distant, which is the measured size of that bound at ' +
      'nucleotide level. The ring recording that is a label, not a move along this axis, ' +
      'because comparing a nucleotide identity with a protein one compares two ' +
      'instruments. Identity also flatters towards familiar as coverage falls, by +17.2 ' +
      'points below half coverage against +7.1 above 90%, which is why coverage is on ' +
      'the colour ramp rather than in a column.'],
  ]},
  {id: 'concern', t: 'What the named leftovers point at', rows: [
    ['Input', 'Every unplaced mate in the library, 2,744,042 of them, and the named ' +
      'relative each one points at where it has a significant one. 214,502 reads over ' +
      '41,829 distinct records.'],
    ['Procedure', 'DIAMOND --sensitive against RefSeq viral proteins plus Swiss-Prot, ' +
      'then the matched record is placed in a host group by an auditable keyword table ' +
      'rather than a classifier, so a reader can disagree with one line. Two separate ' +
      'axes are kept: what the record IS, and what it INFECTS. A cellular record is ' +
      'not given a host, because it is the organism.'],
    ['Threshold', 'E < 1e-5. Host confidence ramps with identity from 40% to 90% and ' +
      'with aligned length from 20 to 60 aa, so a distant or short match is damped ' +
      'towards unknown rather than counted at face value.'],
    ['Denominator', 'Reads, not records: a share of records is a share of the ' +
      'database. Both are shown on hover.'],
    ['Control arm', 'The same search over 400,000 reads the cascade had already named, ' +
      'and a dinucleotide-shuffled null that returns no significant hit at all.'],
    ['Bias direction', 'Towards the well-sequenced. A group is large here when its ' +
      'members are in the reference sets, not when it is abundant in the water, and ' +
      'phage records dominate both. The 22.7% whose host could not be placed is counted ' +
      'as unknown and never as safe. Nothing in this is a measure of harm.']]},
  {id: 'crossplant', t: 'Recurrence across plants', rows: [
    ['Input', '512 panel libraries from 32 wastewater plants on three continents, ' +
      'against reads held out of the library searched.'],
    ['Procedure', 'K-mer witnessing of each query read against the union of the panel, ' +
      'paired within plant so that plants are compared with themselves.'],
    ['Threshold', 'Exact 31-mer witness; depth of the query read recorded beside it, ' +
      'because depth decides the outcome far more than the pile a read came from.'],
    ['Denominator', 'Read pairs of the held-out library.'],
    ['Control arm', 'Identified reads at matched depth, and a swap control against ' +
      'reads from the library itself, which grow 54% by construction.'],
    ['Bias direction', 'Blind to anything present once rather than to anything novel: ' +
      'at depth 1 residue survives 79.0% and identified reads 70.5%, a gap of 8.5 ' +
      'points. Contigs have not been put through it.']]},
  {id: 'nulls', t: 'Shuffled nulls', rows: [
    ['Input', 'The same reads or contigs, shuffled with composition preserved.'],
    ['Procedure', 'The identical search at the identical threshold, in the same pass ' +
      'as the test it calibrates. The shuffle is the Altschul and Erikson Eulerian ' +
      'walk, preserving length and dinucleotide composition exactly, from one shared ' +
      'module that raises on failure. A shuffle that cannot be computed is reported, ' +
      'never replaced by its input: an earlier version returned the input on a ' +
      'dead-ended walk for 324 of 11,361 contigs, and those 324 supplied every hit ' +
      'the null appeared to produce.'],
    ['Threshold', 'As the test.'],
    ['Denominator', 'As the test. A null of zero is reported as a census over the ' +
      'number shuffled, because an empty table and a search that never ran look the ' +
      'same.'],
    ['Control', 'This row is the control.'],
    ['Bias', 'A null states what a field reads on meaningless input. It does not ' +
      'state what a field reads on real input for the wrong reason, which only a ' +
      'positive control can remove. Where none exists the direction of the bias is ' +
      'given alone. Operating point decides a null: the weak-hit band reads 92.3% on ' +
      'shuffled reads with no E-value floor and 0.00% with one.'],
  ]},
];
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

// ---------- header
(function () {
  const r = D.run || {};
  const h = el('header', null, '');
  const inp = (r.inputs || {});
  // Name the sample, not the tool. Whoever opens this is looking at one water on one
  // day, and the accession is the only thing that lets them find it again or compare it
  // with another; "run viewer" told them which program wrote the page.
  const acc = ((inp.r1 || '').match(/([A-Z]{3}\d{6,})/) || [])[1];
  h.append(el('h1', null, acc ? esc(acc) + ', untreated wastewater'
                              : 'Hapax run viewer'));
  h.append(el('p', 'lede',
    'Untargeted RNA sequencing of untreated wastewater. Screening, composition, ' +
    'divergence of assembled viruses, and the unnamed remainder. Every figure carries ' +
    'its operating point and control. No figure scores novelty.'));
  // Two scopes on one page, said once at the top. The classified profile and the
  // residue search read every pair; the cascade reads a block, because it cannot be run
  // over the library in this container. A reader who assumes one scope for all of it
  // will misread the panels against each other.
  const lib = (D.profile || {}).library_pairs;
  h.append(opChips({sample: acc, library: lib ? lib.toLocaleString() + ' read pairs' : null,
    'library type': 'RNA, untargeted shotgun',
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
    'One square is one percent of this run read pairs. Pale squares were named by ' +
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
  if (hc.best_identity_percent != null) card(num(hc.best_identity_percent, 1) + '%',
    'the closest any of those resemblances gets, and over only ' + hc.best_aligned_aa +
    ' amino acids');
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
    'Sewage from one treatment plant in Berlin, sequenced as RNA. The cascade removes ' +
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
  p.append(el('p', 'note',
    '<b>The assembly of this sample is the biggest single bite in the cascade</b>, ' +
    'at ' +
    num(100 * (curve.filter(c => /selfasm_sample_contigs_identified/.test(c.tier))[0]
      || {marginal_removed: 0}).marginal_removed, 1) +
    ' points from ' + gb((curve.filter(c => /selfasm_sample_contigs_identified/
      .test(c.tier))[0] || {}).index_bytes) +
    ', while the largest index here is tens of gigabytes and buys fractions of a point. ' +
    'And <b>several rows take nothing at all</b>: by the time they run, everything they ' +
    'could name has already been named by something smaller. ' +
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
  const warn = el('div', 'callout');
  warn.innerHTML = '<b>Read this before the colours below.</b> No residue subject is ' +
    'close enough for its host to transfer. The best of all ' + hc.subjects + ' is ' +
    num(hc.best_identity_percent, 1) + '% identity over ' + hc.best_aligned_aa +
    ' amino acids, <b>' + hc['at_or_above_0.5'] + ' of ' + hc.subjects +
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
  const ctl = el('div', 'controls');
  const sel = el('select');
  swept.forEach((c, i) => { const o = document.createElement('option');
    o.value = i; o.textContent = c.tier; sel.append(o); });
  const range = document.createElement('input');
  range.type = 'range'; range.id = 'sweep-range'; range.min = 0; range.step = 1;
  const out = el('span', 'mono');
  ctl.append(el('span', 'dim', 'tier'), sel, el('span', 'dim', 'threshold'), range, out);
  if (swept.length) p.append(ctl);
  const draw = () => {
    const c = swept[+sel.value]; if (!c) return;
    const name = Object.keys(c.threshold_sweep)[0];
    const sweep = c.threshold_sweep[name];
    const ks = Object.keys(sweep);
    range.max = ks.length - 1;
    const k = ks[Math.min(+range.value, ks.length - 1)];
    out.textContent = name.replace('_sweep', '') + ' ' + k + ' removes ' + pct(sweep[k]) +
      '  (span ' + num((c.sweep_span_points || {})[name], 2) + ' points across the sweep)';
  };
  sel.onchange = () => { range.value = 0; draw(); };
  range.oninput = draw;
  if (swept.length) draw();

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
    p.append(el('div', 'absent', 'no divergent-stage result in this build'));
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

  // The chart keeps a minimum width and scrolls inside its own box on a narrow screen,
  // rather than shrinking its labels to four pixels. The page itself still never
  // scrolls sideways.
  const host = el('div', 'scroll'); p.append(host);
  const foot = el('div', 'note'); p.append(foot);

  function rowsNow() {
    const want = armSel.value;
    return S.on_axis.filter(r => want === 'both arms' ||
      (want === 'divergent only' ? r.arm === 'divergent'
                                 : r.arm !== 'divergent'));
  }

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
    const rows = rowsNow();
    if (asTable.checked) { host.append(drawTable(rows)); setFoot(rows); return; }
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
        '</title></circle>';
    });

    const hv = S.holdout;
    if (hv && hv.by_role) {
      const ids = Object.keys(hv.by_role).map(k => hv.by_role[k]);
      const lo = Math.min.apply(null, ids), hi = Math.max.apply(null, ids);
      const yy = PAD_T - 20;
      s += '<line x1="' + x(lo) + '" y1="' + yy + '" x2="' + x(hi) + '" y2="' + yy +
        '" stroke="var(--warn)" stroke-width="2" stroke-linecap="round"/>';
      Object.keys(hv.by_role).forEach(k => {
        s += '<circle cx="' + x(hv.by_role[k]) + '" cy="' + yy + '" r="4.5" ' +
          'fill="var(--panel)" stroke="var(--warn)" stroke-width="2"><title>' +
          'the withheld genome, ' + esc(k) + ' at ' + num(hv.by_role[k], 1) +
          '% - true divergence 30.28% at nucleotide level</title></circle>';
      });
      s += '<text x="' + (x(lo) - 10) + '" y="' + (yy + 4) +
        '" text-anchor="end" class="ax">the one genome whose true distance is known, ' +
        num(hv.spread, 1) + ' points across three genes</text>';
    }

    const LY = PAD_T + H_MAIN + GAP + 34;
    s += '<rect x="' + PAD_L + '" y="' + (LY - 26) + '" width="' + (W - PAD_L - PAD_R) +
      '" height="62" fill="var(--v-faint)" rx="3"/>';
    s += '<text x="' + (PAD_L + 8) + '" y="' + (LY - 10) + '" class="ax">' +
      'off the axis: ' + S.cannot_read.length + ' contigs whose identity cannot be read' +
      '</text>';
    const perRow = Math.ceil(S.cannot_read.length / 3);
    const lw = (W - PAD_L - PAD_R - 16) / Math.max(1, perRow);
    S.cannot_read.forEach((r, i) => {
      s += '<circle cx="' + (PAD_L + 8 + (i % perRow) * lw + lw / 2) + '" cy="' +
        (LY + 2 + Math.floor(i / perRow) * 9) +
        '" r="2.6" fill="var(--null)" opacity="0.85"><title>' +
        esc(short(r.name)) + ' - ' + esc(r.why) +
        '\\nreported identity ' + num(r.x, 1) + '%, which is why it is not drawn at it' +
        '</title></circle>';
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
    setFoot(rows);
  }

  function setFoot(rows) {
    const c = S.census || {};
    const thin = S.cannot_read.filter(r => /covered/.test(r.why)).length;
    const shortAln = S.cannot_read.length - thin;
    foot.innerHTML =
      '<b>How to read it.</b> Position is the identity of the conserved gene - the ' +
      'polymerase on ' + S.on_axis.filter(r => r.src === 'polymerase').length +
      ' of ' + S.on_axis.length + ' contigs, the whole-contig hit on the rest, stated ' +
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
      S.cannot_read.length + ' contigs where the measurement failed - ' + thin +
      ' with under ' + Math.round(S.floors.min_coverage * 100) +
      '% of the scored gene aligned and ' + shortAln + ' aligned over under ' +
      S.floors.min_aln_aa + ' aa - and it is a census, not a verdict. ' +
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
      'current filter.';
  }

  armSel.onchange = draw; bands.onchange = draw; asTable.onchange = draw;
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
  ctl.append(el('span', 'dim', 'show'), f);
  p.append(ctl);
  const sc = el('div', 'scroll'); p.append(sc);
  function draw() {
    const want = f.value;
    const tb = el('table');
    tb.innerHTML = '<thead><tr><th>read</th><th>layer</th><th>closest relative</th>' +
      '<th>identity</th><th>E</th><th>occurrences</th><th>profile</th></tr></thead>';
    const body = el('tbody');
    T.rows.filter(r => want === 'all layers' || r.layer === want)
      .slice(0, 400).forEach(r => {
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
      tr.innerHTML = '<td class="mono">' + esc(r.name) + '</td><td class="dim">' +
        esc((r.layer || '-').replace(/_/g, ' ')) + '</td><td>' + shown +
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
  }
  f.onchange = draw; draw();
  p.append(el('p', 'note', T.n_drawn + ' of ' + T.n_total + ' residue pairs sampled for ' +
    'this page with seed ' + T.sample_seed + ', first 400 of the filter shown. The full ' +
    'rows are in reads.jsonl in the run directory.'));
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
  const hero = el('div', 'verdict');
  hero.innerHTML =
    '<div class="vmain">' + (top
      ? 'Highest watchlist detection: <b>' + esc(top.name) + '</b>, at <b>' +
        oneIn(top.fraction) + ' read pairs</b>.' + cite('watchlist')
      : 'No watchlist taxon detected.' + cite('watchlist')) + '</div>' +
    '<div class="vsub">' + hit.length + ' of ' + rows.length + ' detected at ' +
    'confidence ' + CONF + '. No row reaches 1 pair in 100,000. Untargeted RNA ' +
    'sequencing has low sensitivity for human pathogens at these ' +
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
    inner.forEach(s => { s.pct = (100 * s.v / (named || 1)).toFixed(1) + '% of the named part'; });
    const d2 = donut(inner,
      (100 * (by.Viruses || 0) / (named || 1)).toFixed(0) + '%', 'of that is viral',
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

  const hero = el('div', 'verdict');
  hero.innerHTML = '<div class="vmain"><b>' + belowNinety.length + ' of ' + taxa.length +
    '</b> assembled viruses sit <b>below 90% identity</b> to their closest published ' +
    'protein' + cite('contigid') +
    '.</div>' + (deep.length ? '<div class="vmain">Of the ' + rows.length +
      ' worked up gene by gene, <b>' + deep.length + '</b> are complete genomes carried ' +
      'at over 100x read depth' + cite('depth') + '.</div>' : '') +
    '<div class="vsub">Identity here means one gene against the single closest protein ' +
    'in the reference sets, at E&lt;1e-5 over at least 60 aa. Below about 90% that ' +
    'protein is a nearest neighbour rather than an identification, and the name it ' +
    'lends is borrowed.</div>';
  p.append(hero);

  /* THE LINE. One quantity, one chart, no legend: everything that would need one is
     written on the drawing itself. The three bar blocks this replaced used the same
     orange bar for read depth, then distance, then identity, on one screen. */
  if (rows.length) {
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
    const host = el('div', 'scroll'); host.innerHTML = s; p.append(host);
    p.append(el('p', 'note', 'One mark per contig, ' + rows.length + ' of them, stacked ' +
      'where several land on the same identity. Position is the gene the row was ranked on, the polymerase ' +
      'wherever one is called.' + cite('genespread')));
  }

  /* THE FEW THAT MATTER, DRAWN AS WHAT THEY ARE. A bar chart of identities beside a Latin
     name tells a reader neither what the thing is nor what the number is a property of.
     The molecule itself is the picture: genes at their real coordinates, coloured by the
     job they do, with the identity written on the gene it belongs to. */
  const ROLE = {polymerase: ['var(--v1)', 'copies the genome'],
                coat: ['var(--v3)', 'builds the shell'],
                maturation: ['var(--v2)', 'attaches to the host cell'],
                other: ['var(--null)', 'another gene, role not called']};
  const roleOf = o => (o && ROLE[o.role] ? o.role : 'other');

  if (deep.length) {
    const d = el('div', 'bites');
    d.append(el('div', 'pghead', 'the finished genomes, drawn to scale'));
    const allPhage = deep.every(r => (r.orfs_drawn || []).some(o =>
      /ssRNA phage|Levi_coat|Phage_mat-A|RNA_replicase_B/i
        .test((o.title || '') + ' ' + o.pfam.join(' '))));
    if (allPhage) {
      const intro = el('p', 'gintro');
      intro.innerHTML = 'All ' + deep.length + ' are <b>single-stranded RNA ' +
        'bacteriophages</b>: viruses that infect bacteria, not people. Each is a strand ' +
        'of RNA a few kilobases long carrying the three genes such a phage needs, drawn ' +
        'below at their real positions. The arrow is the direction the gene is read and ' +
        'the number on it is identity to the closest published protein.';
      d.append(intro);
    }
    const lg = el('div', 'maplegend');
    lg.innerHTML = ['polymerase', 'coat', 'maturation', 'other'].map(k =>
      '<span><i style="background:' + ROLE[k][0] + '"></i>' + k + ' &mdash; ' +
      ROLE[k][1] + '</span>').join('');
    d.append(lg);

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

      // the molecule
      const MW = 1000, ML = 8, MR = 8, MH = 74, TRK = 46;
      const mx = v => ML + (v / r.contig_nt) * (MW - ML - MR);
      let g = '<svg viewBox="0 0 ' + MW + ' ' + MH + '" width="100%" ' +
        'style="min-width:520px" role="img" aria-label="gene map of ' + esc(r.taxon) +
        '"><line x1="' + mx(0) + '" y1="' + TRK + '" x2="' + mx(r.contig_nt) + '" y2="' +
        TRK + '" stroke="var(--edge)" stroke-width="3"/>';
      os.forEach(o => {
        const k = roleOf(o), x0 = mx(Math.min(o.s, o.e)), x1 = mx(Math.max(o.s, o.e));
        const w = Math.max(3, x1 - x0), tip = Math.min(10, w / 3);
        // drawn as an arrow, because strand is a fact about the molecule and a plain
        // rectangle throws it away
        const pts = o.d < 0
          ? [[x0, TRK], [x0 + tip, TRK - 11], [x0 + w, TRK - 11], [x0 + w, TRK + 11],
             [x0 + tip, TRK + 11]]
          : [[x0, TRK - 11], [x0 + w - tip, TRK - 11], [x0 + w, TRK],
             [x0 + w - tip, TRK + 11], [x0, TRK + 11]];
        g += '<polygon points="' + pts.map(q => q[0].toFixed(1) + ',' + q[1]).join(' ') +
          '" fill="' + ROLE[k][0] + '" opacity="' + (o.sig ? 1 : 0.45) + '"><title>' +
          esc(o.role || 'role not called') + ', ' + o.aa + ' aa' +
          (o.sig ? ', ' + num(o.pident, 1) + '% to ' + esc(o.title || '') : ', no ' +
            'significant protein hit') + (o.pfam.length ? '\\\ndomain ' +
            esc(o.pfam.join(', ')) : '') + '</title></polygon>';
        if (o.sig && w > 74) {
          g += '<text x="' + (x0 + w / 2) + '" y="' + (TRK + 4) +
            '" text-anchor="middle" style="font-size:12px;font-weight:600;' +
            'fill:var(--panel);font-family:var(--sans)">' + num(o.pident, 0) +
            '%</text>';
        }
        if (w > 54) {
          g += '<text x="' + (x0 + w / 2) + '" y="' + (TRK - 18) +
            '" text-anchor="middle" class="ax">' + esc(k === 'other'
              ? (o.sig ? 'other' : 'no hit') : k) + '</text>';
        }
      });
      g += '<text x="' + mx(0) + '" y="' + (TRK + 26) + '" class="ax">0</text>';
      g += '<text x="' + mx(r.contig_nt) + '" y="' + (TRK + 26) +
        '" text-anchor="end" class="ax">' + kb + ' kb</text>';
      g += '</svg>';
      const hostd = el('div', 'scroll'); hostd.innerHTML = g; card.append(hostd);

      // what it is, in words a reader can check against the picture above
      const roles = os.filter(o => o.sig).map(roleOf);
      const phage = os.some(o => /ssRNA phage|Levi_coat|Phage_mat-A|RNA_replicase_B/i
        .test((o.title || '') + ' ' + o.pfam.join(' ')));
      const envOnly = os.filter(o => o.sig).length &&
        os.filter(o => o.sig).every(o => /ssRNA phage SRR|MAG:/.test(o.title || ''));
      const far = os.filter(o => o.sig).sort((a, b) => a.pident - b.pident)[0];
      const say = el('div', 'gsay');
      say.innerHTML =
        (allPhage ? '' : (phage ? '<b>A single-stranded RNA bacteriophage.</b> ' :
          '<b>A virus</b>, architecture not called from its genes. ')) +
        '<b>' + kb + ' kb</b> at <b>' + Math.round(rs.mapped_mean_depth) +
        'x</b>, covered ' + (rs.breadth_covered > 0.98 ? 'end to end'
          : pct(rs.breadth_covered) + ' of its length') +
        (rs.evenness_within_2x_of_median < 0.6
          ? ' but unevenly, so part of it may be a different molecule' : '') + '. ' +
        (far ? 'Furthest from anything published: the <b>' +
          esc(far.role || 'unnamed gene') + '</b> at ' + num(far.pident, 1) + '%. ' : '') +
        (roles.length < os.length
          ? (os.length - roles.length === 1
             ? 'One of its ' + os.length + ' genes has no significant match at all. '
             : (os.length - roles.length) + ' of its ' + os.length + ' genes have no ' +
               'significant match at all. ') : '') +
        (envOnly ? '<span class="dim">Its nearest published relative is itself an ' +
          'unnamed sequence from another wastewater sample.</span>' : '');
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
      'human-associated' + cite('concern') + '. Those ' + hum + ' reads are the enteric ' +
      'viruses sewage is expected to carry, matched over about 34 amino acids. ' +
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
  if (!K || !b) return;
  const p = panel('The unnamed remainder',
    'The fraction no reference names.' + cite('cascade'));
  const own = tier('selfasm_sample_contigs_identified_norrna').marginal_removed;
  const rr = tier('rrna').tier_alone_removed;

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
    const CAN_NAME = {k2_viral: 1, k2_standard08: 1, k2_pluspf08: 1, k2_pluspf16: 1,
                      prot_viral: 1, prot_sprot: 1};
    let named = 0;
    curve.forEach(c => { if (CAN_NAME[c.tier]) named += c.marginal_removed || 0; });
    const ident = brk.identified.point;
    const bands = [
      {lab: 'named', v: named, col: 'var(--v1)',
       sub: 'a taxon or an accession'},
      {lab: 'ribosomal RNA', v: rr, col: 'var(--v1)', op: 0.45,
       sub: 'identifies a ribosome, not an organism'},
      {lab: 'recognised, unnamed', v: Math.max(0, ident - named - rr),
       col: 'var(--v-mute)', sub: 'mostly own-sample contigs (' + pct(own) +
         '), rest adapter and low complexity'},
      {lab: 'assembles, unidentified', v: brk.assembled_unidentified.point,
       col: 'var(--v2-soft)', sub: 'median depth 6'},
      {lab: 'residual', v: brk.residual.point, col: 'var(--v2)',
       sub: 'median depth 1'},
    ].filter(b => b.v > 0);
    const tot = bands.reduce((a, b) => a + b.v, 0) || 1;
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
  note.innerHTML = '<b>Unnamed is not novel.</b> A search establishes that sequence ' +
    'is already known. It cannot establish novelty, because any search-derived score ' +
    'tracks distance to the database. Artefact signatures account for under 1% of this ' +
    'pile, it contains no abundant unknown genome, and its median depth is 1.' +
    cite('residue');
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
  const p = panel('Named against recognised: cascade reconciliation',
    'Reconciling 7% named against 77% identifiable.');
  // Which tiers can name is a property of the tier, declared here where it can be read.
  // A Kraken2 or DIAMOND hit carries a taxon or an accession; an artefact filter, a
  // SILVA match and a match to a contig from this very sample tell you the sequence has
  // been seen before and nothing about what it is.
  const CAN_NAME = {k2_viral: 1, k2_standard08: 1, k2_pluspf08: 1, k2_pluspf16: 1,
                    prot_viral: 1, prot_sprot: 1};
  const WHY = {artefact: 'adapter and low complexity',
               rrna: 'ribosomal RNA',
               selfasm_sample_contigs_rrna: 'own-sample ribosomal contigs',
               selfasm_sample_contigs_identified_norrna: 'own-sample unnamed contigs'};
  let named = 0, unnamed = 0; const unnamedRows = [];
  curve.forEach(c => {
    const m = c.marginal_removed || 0;
    if (CAN_NAME[c.tier]) { named += m; }
    else { unnamed += m; if (m > 0.002) unnamedRows.push([c.tier, m]); }
  });
  const k2here = 1 - (((P.at_confidence['0.1'] || {}).by_superkingdom || [])
    .filter(g => g.group === 'unclassified')[0] || {}).fraction;
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
  rowbar('<b>named</b>, whole cascade', named, 'var(--v1)',
    'sum of every cascade tier returning a taxon or accession. Agrees with the row ' +
    'above.');
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

// ---------- summary and detail
// Runs last, after every panel has rendered. It reorganises rather than regenerates: the
// detail panels keep the markup they were built with and are only moved inside a
// disclosure, so a figure can never disagree with the summary that cites it.
(function () {
  // The front page answers a public-health reader's questions in the order they ask
  // them. Everything else on this page is how those answers were arrived at, and goes
  // below the fold: it was all on the front page before, and the front page was a
  // method report wearing a sample report's clothes.
  const SUMMARY = ['Watchlist screen', 'What this water is',
                   'How far the assembled viruses sit from anything published',
                   'The unnamed remainder', 'Methods'];
  const panels = Array.prototype.filter.call(app.children,
    n => n.tagName === 'SECTION' && n.classList.contains('panel'));
  const detail = panels.filter(n => {
    const h = n.querySelector('h2');
    return h && SUMMARY.indexOf(h.textContent) < 0;
  });
  if (!detail.length) return;
  const firstSentence = t => {
    const whole = String(t || '').replace(/\\s+/g, ' ').trim();
    const m = whole.match(/^.*?\\.(?=\\s|$)/);
    let out = m ? m[0] : whole;
    if (out.length > 130) out = out.slice(0, 130).replace(/\\s+\\S*$/, '') + '...';
    return out;
  };
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", type=Path)
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
    b = bundle(args.rundir, args.scatter_per_pile, args.table_rows, args.seed)
    out = args.out or (args.rundir / "web" / "index.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(b, separators=(",", ":"))
    if "</script" in payload:
        raise SystemExit("the bundle contains a closing script tag and would break the "
                         "page; a value in the run directory needs escaping")
    # A citation that points at no method entry is worse than no citation: it reads as
    # a reference and goes nowhere. The ids are both in the page source, so the check is
    # static and runs on every build rather than on whoever happens to click one.
    cited = set(re.findall(r"cite\('([a-z_]+)'\)", PAGE))
    defined = set(re.findall(r"\{id: '([a-z_]+)', t:", PAGE))
    missing = sorted(cited - defined)
    if missing:
        raise SystemExit(f"citations with no method entry: {', '.join(missing)}")
    unused = sorted(defined - cited)

    page = PAGE.replace("__BUNDLE__", payload)
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
    print(f"{out}  {size / 1e6:.2f} MB  {checked}  "
          f"residue rows {b['table']['n_total']}, drawn {b['table']['n_drawn']}; "
          f"scatter {b['scatter'].get('residue', {}).get('n_drawn', 0)}+"
          f"{b['scatter'].get('known', {}).get('n_drawn', 0)}; "
          f"{len(cited)} of {len(defined)} methods cited"
          + (f", never cited: {', '.join(unused)}" if unused else ""))
    if size > 16_000_000:
        raise SystemExit(f"{size} bytes is over the 16 MB page limit; lower "
                         f"--scatter-per-pile or --table-rows")
    print(f"built in {time.perf_counter() - t0:.1f}s "
          f"({platform.python_version()})", flush=True)


if __name__ == "__main__":
    main()
