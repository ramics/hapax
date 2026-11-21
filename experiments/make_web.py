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
    mix = {}
    for r in rows:
        mix[r["host_group"]] = mix.get(r["host_group"], 0) + 1
    R["host_mix"] = mix
    R["subjects_total"] = len(rows)

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


def load(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


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
  --sans: "IBM Plex Sans", system-ui, -apple-system, Segoe UI, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
  --g-phage: #93a3b4; --g-plant: #74b069; --g-fungal: #a68cc8; --g-invert: #4fa89f;
  --g-vert: #d2a04a; --g-human: #e0695c; --g-unknown: #68756f;
} }
:root[data-theme="dark"] {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
  --g-phage: #93a3b4; --g-plant: #74b069; --g-fungal: #a68cc8; --g-invert: #4fa89f;
  --g-vert: #d2a04a; --g-human: #e0695c; --g-unknown: #68756f;
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
a { color: var(--identified); }
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
    'What is in this water, what we screened for and did not find, and how far the ' +
    'things we recognise have drifted from the references. Then the larger part: the ' +
    'sequence nothing could name, and what is and is not known about it. Every number ' +
    'is printed beside the operating point and the control that make it readable, and ' +
    'nothing here scores novelty: an index can say whether something is already known ' +
    'and cannot say whether it is new.'));
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
  const p = panel('What is in the water',
    'Each square is one percent of the read pairs in this run. Pale squares were named ' +
    'by searching; the dark ones are what no search could name.');
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
    'The dark squares are the subject of this page. They are not a discovery: ' +
    'unnameable sequence is the normal state of a wastewater metagenome, and the ' +
    'measurements below are an attempt to say what is actually in that pile, how far ' +
    'from anything known it sits, and where the instrument stops being able to tell.'));
  p.append(opChips((t || {}).operating_point || {}));
  app.append(p);
})();

// ---------- in plain words
// A reader who is not a bioinformatician should be able to leave this page knowing what
// was done and what was found. Every figure in it is pulled from the same files the rest
// of the page reads, so it cannot drift from the panels below it.
(function () {
  const t = D.tier_table, S = D.strains, L = D.linkage, C = D.census;
  const p = panel('In plain words',
    'What this run did, and what came out of it.');
  const b = (t && t.breakdown) || {};
  const lede = el('p', 'lede',
    'This is sewage from one treatment plant in Berlin, sequenced as RNA. Most of what ' +
    'is in it is already known to science, so the job is to take the known away and look ' +
    'hard at what is left. Nothing here guesses that a sequence is new. The page says ' +
    'how far a piece of sequence is from its nearest known relative, and how confident ' +
    'that distance is, and stops there.');
  p.append(lede);
  const g = el('div', 'big');
  function big(v, lab) { const d = el('div');
    d.append(el('div', 'bignum', v)); d.append(el('div', 'biglab', lab)); g.append(d); }
  if (b.identified) big(pct(b.identified.point), 'of the library was identified outright');
  if (b.assembled_unidentified) big(pct(b.assembled_unidentified.point),
    'assembled into something, but nothing could name it');
  if (b.residual) big(pct(b.residual.point),
    'left over: no name, and no contig either');
  if (S && S.subjects_total != null) big(String(S.subjects_total),
    'known organisms the leftovers resemble from a distance');
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
  const p = panel('How the library splits',
    'Every read pair in this run, taken apart in three steps. Each panel is the ' +
    'highlighted slice of the one before it, blown up to full height.');

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
    '<b>The second join changes the unit and the sample, and that is marked on it.</b> ' +
    'The first two panels are read PAIRS over the whole block this run read. The third ' +
    'and fourth are MATES, and only the ' + K.searched.toLocaleString() + ' of them the ' +
    'per-read report searched, so the wedge between them is a re-basing and not a ' +
    'subdivision. Reading a percentage straight across the whole figure would be wrong, ' +
    'which is why none is drawn that way.'));
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
  const p = panel('Each index takes a bite',
    'The cascade, one row per index, in the order they ran. The dark slice is what that ' +
    'index removed that nothing before it had already taken.');
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
    'Two things this shape says that the numbers alone do not. <b>The sample\\'s own ' +
    'assembly is the biggest single bite in the cascade</b> at ' +
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
  const p = panel('Phage, virus, or something living',
    'Of the leftover reads, the few that have a named relative at all, and what kind of ' +
    'thing that relative is.');
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
  const p = panel('Worth a second look',
    'The named relatives the leftover reads point at, ordered by how much of a person\\'s ' +
    'attention each one deserves. Four measured things go in and all four are shown on ' +
    'every row, so the order can be argued with one row at a time.');
  if (!R) { p.append(el('div', 'absent',
    'ranking.json absent: run experiments/concern_ranking.py over this directory'));
    app.append(p); return; }
  p.append(scopeLine(R));

  // The honest headline of the whole panel, and it is a negative. Put it above the
  // leaderboard rather than under it: a reader who scrolls the cards without seeing this
  // will read the host colours as claims about the sample.
  const hc = R.host_confidence_summary || {};
  const warn = el('div', 'callout');
  warn.innerHTML = '<b>Read this before the colours below.</b> Nothing in this residue is ' +
    'close enough to a known organism for that organism\\'s host to transfer. The best ' +
    'match of all ' + hc.subjects + ' is ' + num(hc.best_identity_percent, 1) +
    '% identity over ' + hc.best_aligned_aa + ' amino acids, <b>' + hc['at_or_above_0.5'] +
    ' of ' + hc.subjects + '</b> clear even half confidence, and the median is ' +
    num(hc.median, 2) + '. So the host groups here are labels on the matched database ' +
    'record, not claims about what is in the water, and the ranking below is really ' +
    'distance and quantity wearing a colour. <b>That is the finding, not a caveat:</b> ' +
    'this is what it looks like when a surveillance question runs past the edge of what ' +
    'the references can answer.';
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
  const p = panel('Where the relatives sit',
    'The named relatives at once: how far the reads are from each one, against how many ' +
    'reads point at it. One picture of what the leftovers are related to, and how ' +
    'distantly.' + (pn ? ' Drawn from a seeded sample of ' + pn.drawn.toLocaleString() +
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
  const p = panel('The breakdown',
    'One bar over the library: what was identified, what assembled without being ' +
    'identified, and what is left. Each segment carries the interval from this run\\'s ' +
    'own counts.');
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
  const p = panel('The cascade',
    'Each index in turn: what it removes alone, what it adds, and what is left. The ' +
    'slider is the tier\\'s own threshold sweep, because the threshold moves a tier by ' +
    'more than a database\\'s size does.');
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
  const p = panel('What each index costs',
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
  const p = panel('Can a read find a read it is related to?',
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
  const p = panel('The residue in layers, most to least worth a look',
    'A layer is a sort order with its evidence attached, not a score. Each is laid ' +
    'beside the same rule applied to reads the cascade did identify, and where that ' +
    'control was never computed the page says so instead of showing a separation.');
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
  const p = panel('Every field beside its calibration',
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
  if (F.chimaera_split_rate != null) {
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
  const p = panel('What the residue sits next to',
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
  const p = panel('New strains of known viruses',
    'Grouped by the relative the reads are nearest to, because a new strain is a subject ' +
    'several residue reads point at from a distance, not a read with a low score. ' +
    'Subjects are ordered by how many reads point at them and never by identity: ' +
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
  const p = panel('The joint space',
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
  const p = panel('This run\\'s divergence ladder',
    'A holdout genome spiked at known divergence and pushed through the same cascade. It ' +
    'is the only positive control the tool has, and the rungs are ordinal.');
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

// ---------- the residue table
(function () {
  const T = D.table;
  const p = panel('The residue, read by read',
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
    'A fixed panel, not a footnote. Each names the measurement that established it, so it ' +
    'can be checked rather than taken.');
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
  const p = panel('What a leftover read sits next to',
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
  const p = panel('Layering the leftovers',
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

// ---------- public health: what we screened for
// The first thing a public-health reader asks of a sewage sample is not what is novel in
// it. It is whether anything they watch for is present, and what was looked for and not
// found. An empty row is a result here, so every panel entry is drawn whether or not it
// was seen.
(function () {
  const P = D.profile;
  const p = panel('What we screened for',
    'A committed list of organisms a public-health reader would want asked about, every ' +
    'one reported whether it was found or not. Found means reads in this library share ' +
    'enough k-mers with a reference for the classifier to place them there. It is not a ' +
    'diagnosis, not a viable organism, and not a case.');
  if (!P) {
    p.append(el('div', 'absent', 'public_health_profile_s*.json absent: run ' +
      'experiments/public_health_profile.py over this library'));
    app.append(p); return;
  }
  const CONF = '0.1';
  const at = P.at_confidence[CONF] || {};
  const rows = at.panel || [];
  const hit = rows.filter(r => r.detected);
  const miss = rows.filter(r => !r.detected);
  const lo = (P.at_confidence['0.2'] || {}).panel || [];
  const hi = (P.at_confidence['0.0'] || {}).panel || [];
  const byName = l => { const m = {}; l.forEach(r => { m[r.name] = r; }); return m; };
  const LO = byName(lo), HI = byName(hi);

  const nul = P.shuffled_null || {};
  const head = el('div', 'callout');
  head.innerHTML = '<b>' + hit.length + ' of ' + rows.length + ' found, ' + miss.length +
    ' not found</b>, over ' + (at.pairs || 0).toLocaleString() + ' read pairs at ' +
    'classifier confidence ' + CONF + '. <b>That threshold is not a preference, it is ' +
    'the calibration:</b> the same panel run over this library with every read shuffled ' +
    'lights up ' + ((nul['0.0'] || {}).panel_rows_detected) + ' of ' + rows.length +
    ' rows at confidence 0 and ' + ((nul[CONF] || {}).panel_rows_detected) + ' at ' +
    CONF + '. Counts at all three thresholds are on every found row, smallest last; a ' +
    'detection that survives the strictest one is the strongest claim here.';
  p.append(head);

  // One in N, not a percentage. Every rate on this panel is of the order of one pair in
  // a million, where a percentage to two places prints 0.00% for everything and a reader
  // learns nothing. One in N is also the unit this project talks in.
  const oneIn = x => {
    if (!x) return '-';
    const n = 1 / x;
    return n >= 1e6 ? '1 in ' + (n / 1e6).toFixed(1) + ' million'
         : n >= 1e3 ? '1 in ' + Math.round(n / 1e3) + ',000'
         : '1 in ' + Math.round(n);
  };
  const grid = el('div', 'panelgrid');
  // The panel file's own order, which groups the enteric viruses wastewater
  // surveillance was built on before the organisms it was later extended to. Sorting by
  // how many were found would put the sample in charge of the running order.
  const groups = [];
  rows.slice().sort((a, b) => rows.indexOf(a) - rows.indexOf(b))
    .forEach(r => { if (groups.indexOf(r.group) < 0) groups.push(r.group); });
  groups.forEach(g => {
    const mine = rows.filter(r => r.group === g);
    const sec = el('div', 'pgroup');
    sec.append(el('div', 'pghead', esc(g) + ' <span class="dim">' +
      mine.filter(r => r.detected).length + ' of ' + mine.length + '</span>'));
    mine.forEach(r => {
      const row = el('div', 'prow' + (r.detected ? ' found' : ''));
      const strict = (LO[r.name] || {}).pairs || 0;
      if (r.detected) {
        const sweep = [(HI[r.name] || {}).pairs || 0, r.pairs, strict];
        row.innerHTML = '<div class="pname"><b>' + esc(r.name) + '</b>' +
          (strict ? ' <span class="chip">holds at 0.2</span>' : '') + '</div>' +
          '<div class="pval"><span class="mono">' + r.pairs.toLocaleString() +
          '</span> pair' + (r.pairs === 1 ? '' : 's') + ', ' + oneIn(r.fraction) +
          '<span class="dim mono"> &nbsp;' + sweep.join(' / ') + '</span></div>';
      } else {
        row.innerHTML = '<div class="pname">' + esc(r.name) + '</div>' +
          '<div class="pval dim">not found, under ' +
          oneIn(r.upper_bound_if_absent) + '</div>';
      }
      row.title = r.why_on_the_panel;
      sec.append(row);
    });
    grid.append(sec);
  });
  p.append(grid);
  const foot = el('div', 'note');
  foot.innerHTML = '<b>What a not-found row does not mean.</b> ' + esc(P.panel_note);
  p.append(foot);
  app.append(p);
})();

// ---------- public health: composition, with the threshold as a control
(function () {
  const P = D.profile;
  const p = panel('What one reference index can name',
    'This is the REACH OF ONE INDEX, not the composition of the water. Every read pair ' +
    'put in the broadest group Kraken2 can place it in, against one database at one ' +
    'threshold. Move the threshold and watch the picture move: that movement is the ' +
    'single most important thing to understand about any number like this one.');
  if (!P) {
    p.append(el('div', 'absent', 'public_health_profile_s*.json absent'));
    app.append(p); return;
  }
  const COL = {Viruses: 'var(--g-phage)', Bacteria: 'var(--g-vert)',
    Archaea: 'var(--g-fungal)', Eukaryota: 'var(--g-plant)',
    unclassified: 'var(--null)', 'classified, above superkingdom': 'var(--g-unknown)'};
  // "Unclassified" is the classifier's word and it reads, to anyone who is not holding
  // the classifier in their head, as "unknown sequence". It is not: it is sequence this
  // one database could not place at this one threshold. The page renames it on the way
  // in, because the mislabel is what made a reader take the 93% for the water.
  const GLABEL = g => g === 'unclassified' ? 'not named by this index at this threshold'
                    : g === 'classified, but above any of these groups'
                      ? 'named, but no deeper than these groups is defined' : g;
  const confs = Object.keys(P.at_confidence).sort();
  const ctl = el('div', 'controls');
  const sel = document.createElement('select');
  confs.forEach(c => { const o = document.createElement('option');
    o.value = c; o.textContent = 'confidence ' + c; sel.append(o); });
  sel.value = confs.indexOf('0.1') >= 0 ? '0.1' : confs[0];
  ctl.append(el('span', 'dim', 'classifier threshold'), sel);
  p.append(ctl);

  const bar = el('div', 'bar big');
  const lg = el('div', 'legend');
  // The whole library is mostly unclassified at any usable threshold, so a single bar
  // is one grey block and tells a reader nothing about the part that was placed. The
  // second bar renormalises over the classified part alone and is labelled as such:
  // it answers "of what we could name, how much is viral", which is the question, and
  // it must never be read as a share of the library.
  const bar2 = el('div', 'bar big');
  const lg2 = el('div', 'legend');
  const cols = el('div', 'twocol');
  p.append(el('div', 'pghead', 'every read pair in the library'), bar, lg,
           el('div', 'pghead', 'the classified part alone, renormalised'), bar2, lg2,
           cols);

  function taxaList(title, list, total, col) {
    const box = el('div', 'tcol');
    box.append(el('div', 'pghead', title));
    if (!list || !list.length) {
      box.append(el('div', 'dim', 'nothing placed in this group'));
      return box;
    }
    const top = list[0].pairs || 1;
    list.slice(0, 12).forEach(t => {
      const r = el('div', 'trow');
      const b = el('div', 'tbar');
      const f = document.createElement('span');
      f.style.width = (100 * t.pairs / top) + '%';
      f.style.background = col;
      b.append(f);
      // Counts, with a percentage whose precision follows the value. A column of
      // "0.00%" is what a fixed two places does to a profile whose tail is most of it.
      const pc = 100 * t.pairs / total;
      const sig = pc >= 1 ? pc.toFixed(2)
                : pc >= 0.01 ? pc.toFixed(3) : pc.toExponential(1);
      r.append(el('div', 'tname', esc(t.name)), b,
        el('div', 'tnum', t.pairs.toLocaleString()),
        el('div', 'tpct dim', sig + '%'));
      r.title = t.pairs.toLocaleString() + ' pairs of ' + total.toLocaleString();
      box.append(r);
    });
    return box;
  }

  function draw() {
    const a = P.at_confidence[sel.value] || {};
    const total = a.pairs || 1;
    bar.innerHTML = '';
    (a.by_superkingdom || []).forEach(g => {
      const sp = document.createElement('span');
      sp.style.width = (100 * g.fraction) + '%';
      sp.style.background = COL[g.group] || 'var(--null)';
      sp.title = GLABEL(g.group) + ': ' + g.pairs.toLocaleString() + ' pairs (' +
        (100 * g.fraction).toFixed(2) + '%)';
      bar.append(sp);
    });
    const cl = (a.by_superkingdom || []).filter(g => g.group !== 'unclassified');
    const clTot = cl.reduce((x, g) => x + g.pairs, 0) || 1;
    bar2.innerHTML = '';
    cl.forEach(g => {
      const sp = document.createElement('span');
      sp.style.width = (100 * g.pairs / clTot) + '%';
      sp.style.background = COL[g.group] || 'var(--null)';
      sp.title = g.group + ': ' + (100 * g.pairs / clTot).toFixed(2) +
        '% of the classified part';
      bar2.append(sp);
    });
    lg2.innerHTML = cl.map(g =>
      '<span><i style="background:' + (COL[g.group] || 'var(--null)') + '"></i>' +
      esc(GLABEL(g.group)) + ' <b>' + (100 * g.pairs / clTot).toFixed(1) + '%</b></span>')
      .join('') + '<span class="dim">of the ' + clTot.toLocaleString() +
      ' pairs placed at all</span>';
    lg.innerHTML = (a.by_superkingdom || []).map(g =>
      '<span><i style="background:' + (COL[g.group] || 'var(--null)') + '"></i>' +
      esc(GLABEL(g.group)) + ' <b>' + (100 * g.fraction).toFixed(2) + '%</b> <span class="dim ' +
      'mono">' + g.pairs.toLocaleString() + '</span></span>').join('');
    cols.innerHTML = '';
    cols.append(taxaList('viruses, most reads first', (a.top_taxa || {}).Viruses, total,
                         COL.Viruses),
                taxaList('bacteria, most reads first', (a.top_taxa || {}).Bacteria, total,
                         COL.Bacteria));
  }
  sel.onchange = draw; draw();

  // THE COMPARISON A READER NEEDS AT EXACTLY THIS POINT. Left alone, the 93% above
  // reads as "93% of this water is unknown", which is what happened, and it sits next
  // to a residue figure of 11% and looks like a contradiction. It is not one: the two
  // numbers come from instruments with different jobs. Drawing the cascade's own
  // removal beside this index, on one scale, with the tiers split by whether they can
  // produce a NAME, is the only version of this that cannot be misread.
  //
  // Which tiers can name is a property of the tier, declared here and visible: a
  // Kraken2 or DIAMOND hit carries a taxon or an accession, while an artefact filter, a
  // SILVA match and a match to a contig assembled from this very sample tell you the
  // sequence has been seen before and nothing about what it is.
  const CAN_NAME = {k2_viral: 1, k2_standard08: 1, k2_pluspf08: 1, k2_pluspf16: 1,
                    prot_viral: 1, prot_sprot: 1};
  const WHY = {artefact: 'adapter and low complexity, removed and not named',
               rrna: 'matches SILVA ribosomal RNA, which says it is a ribosome and not ' +
                     'whose',
               selfasm_sample_contigs_rrna: 'matches this sample own ribosomal contigs',
               selfasm_sample_contigs_identified_norrna:
                 'matches a contig assembled from this very sample'};
  const curve = (D.tier_table || {}).curve || [];
  const brk = (D.tier_table || {}).breakdown;
  if (curve.length && brk) {
    let named = 0, unnamed = 0;
    const unnamedRows = [];
    curve.forEach(c => {
      const m = c.marginal_removed || 0;
      if (CAN_NAME[c.tier]) { named += m; } else {
        unnamed += m;
        if (m > 0.002) unnamedRows.push([c.tier, m]);
      }
    });
    const a01 = P.at_confidence['0.1'] || {};
    const k2here = 1 - ((a01.by_superkingdom || [])
      .filter(g => g.group === 'unclassified')[0] || {}).fraction;
    const box = el('div', 'cmp');
    box.append(el('div', 'pghead',
                  'the same library, by the two instruments on this page'));
    const rowbar = (lab, v, col, sub) => {
      const r = el('div', 'crow');
      r.innerHTML = '<div class="clab">' + lab + '</div>' +
        '<div class="cbar"><span style="width:' + (100 * v).toFixed(2) + '%;background:' +
        col + '"></span></div>' +
        '<div class="cnum mono">' + pct(v) + '</div>' +
        '<div class="csub dim">' + sub + '</div>';
      box.append(r);
    };
    rowbar('<b>named</b> by this index', k2here, 'var(--g-phage)',
      'Kraken2, PlusPF-08, confidence 0.1, whole library. A taxon you can act on.');
    rowbar('<b>named</b> by the whole cascade', named, 'var(--g-phage)',
      'every cascade tier that returns a taxon or an accession, added up. That it lands ' +
      'beside the figure above is the cross-check.');
    rowbar('<b>recognised but not named</b>', unnamed, 'var(--identified)',
      unnamedRows.map(r => (WHY[r[0]] || r[0]) + ', ' + pct(r[1])).join('; ') + '.');
    rowbar('left over after everything', brk.residual.point, 'var(--residue)',
      'the pile the rest of this page is about, plus ' +
      pct((brk.assembled_unidentified || {}).point) + ' that assembles into contigs ' +
      'nothing can name.');
    const cnote = el('div', 'callout');
    cnote.innerHTML = '<b>93% unnamed and 11% left over are both true, and they are not ' +
      'the same measurement.</b> Naming a read means putting a taxon on it, which one ' +
      'reference index does for ' + pct(k2here) + ' of pairs and the whole cascade for ' +
      pct(named) + '. Recognising a read means matching it to sequence we have seen ' +
      'before, which is a far weaker requirement and covers ' + pct(unnamed) + ' more: ' +
      'mostly ribosomal RNA and contigs assembled from this very sample, neither of ' +
      'which tells you what the organism is. <b>Most of what this pipeline removes as ' +
      'known, it removes by matching sequence nobody has named.</b> The scopes differ ' +
      'too: the index figure is the whole library, the cascade figures are its ' +
      ((D.parameters || {}).pairs || 0).toLocaleString() + '-pair block.';
    p.append(box, cnote);
  }

  const note = el('div', 'callout');
  // The JSON fields are sentences without terminators, so the page supplies them
  // rather than the files carrying page punctuation.
  note.innerHTML = '<b>This is RNA, not organisms.</b> ' + esc(P.reading) +
    '. <b>And the threshold is the measurement.</b> ' + esc(P.bias) + '.';
  p.append(note);
  app.append(p);
})();

// ---------- public health: known viruses, reference or drifted
(function () {
  const S = D.strain_divergence;
  const p = panel('Known viruses: the reference, or a drifted strain',
    'For each virus this sample assembled, how far its own sequence has drifted from the ' +
    'closest reference protein. Low identity over a LONG alignment is a divergent strain ' +
    'and is the row worth attention; low identity over a short one says nothing and is ' +
    'marked unreadable rather than scored.');
  if (!S) {
    p.append(el('div', 'absent', 'known_strain_divergence_s*.json absent: run ' +
      'experiments/known_strain_divergence.py'));
    app.append(p); return;
  }
  const taxa = (S.taxa || []).filter(t => t.median_identity_percent != null);
  const head = el('div', 'callout');
  const div = taxa.filter(t => t.median_identity_percent < 90);
  head.innerHTML = '<b>' + taxa.length + ' taxa have an alignment long enough to read an ' +
    'identity from, and ' + div.length + ' of them sit below 90% identity to their ' +
    'closest reference protein.</b> This measures drift from the DATABASE. It is blind ' +
    'to strain mixtures inside the sample, because an assembler fuses lineages below ' +
    'about 3.8% divergence at any depth.';
  p.append(head);
  // Most drifted first. The file is ordered by assembled mass, which is the right
  // order for an archive and the wrong one here: a reader opening this panel is asking
  // which of the things we recognise has moved furthest from the reference.
  const list = el('div', 'drows');
  taxa.slice().sort((a, b) => a.median_identity_percent - b.median_identity_percent)
    .slice(0, 24).forEach(t => {
    const r = el('div', 'drow');
    const far = t.median_identity_percent < 90;
    r.innerHTML = '<div class="dname"><b>' + esc(t.name) + '</b>' +
      (t.on_the_panel ? ' <span class="chip">on the panel</span>' : '') + '</div>' +
      '<div class="dbar"><span style="width:' + t.median_identity_percent + '%;' +
      'background:' + (far ? 'var(--residue)' : 'var(--identified)') + '"></span></div>' +
      '<div class="dnum mono">' + num(t.median_identity_percent, 1) + '%</div>' +
      '<div class="dmeta dim">over ' + t.median_aligned_aa + ' aa, ' + t.contigs +
      ' contig' + (t.contigs === 1 ? '' : 's') + ', longest ' +
      t.longest_nt.toLocaleString() + ' nt</div>';
    r.title = t.reading;
    list.append(r);
  });
  p.append(list);
  const bs = el('div', 'note');
  bs.innerHTML = '<b>What this cannot see.</b> ' +
    (S.blind_spots || []).map(b => esc(b)).join('. ') + '.';
  p.append(bs);
  app.append(p);
})();

// ---------- public health: the part nothing could name
// The whole of the rest of this page is about this pile. On the front page it gets one
// panel and a number, because a public-health reader needs to know it exists and how big
// it is before they need to know how it was characterised.
(function () {
  const R = D.ranking, K = R && R.kind_summary, b = (D.tier_table || {}).breakdown;
  if (!K || !b) return;
  const p = panel('What nothing could name',
    'The reads no index places, which on this library is most of them. This is the ' +
    'research question the rest of this page is about; here it is as one number and ' +
    'what is known about it.');
  const big = el('div', 'cards3');
  const card = (v, l, c) => {
    const d = el('div');
    const n = el('div', 'bignum', v); if (c) n.style.color = c;
    d.append(n); d.append(el('div', 'biglab', l));
    return d;
  };
  big.append(
    card(pct(b.residual.point), 'of the library survives every index the cascade runs',
         'var(--residue)'),
    card(pct(1 - (K.named / (K.searched || 1))),
         'of those reads point at nothing nameable even in protein', 'var(--residue)'),
    card(pct((b.identified || {}).point), 'identifiable, which is what gets subtracted',
         'var(--identified)'));
  p.append(big);
  // The reconciliation itself lives in the index panel above, where a reader first
  // meets the number that needs it. Repeating it here would be the same paragraph
  // twice; a pointer is enough, and it carries the one figure this panel needs.
  const cur = (D.tier_table || {}).curve || [];
  const own = (cur.filter(c => c.tier === 'selfasm_sample_contigs_identified_norrna')[0]
               || {}).marginal_removed;
  if (own != null) {
    const rec = el('div', 'note');
    rec.innerHTML = 'Identifiable is not the same as named: ' + pct(own) + ' of this ' +
      'library is removed purely by matching contigs assembled from the sample itself, ' +
      'which carry no name at all. The panel above puts the two side by side.';
    p.append(rec);
  }

  const note = el('div', 'callout');
  note.innerHTML = '<b>Unnamed is not novel.</b> An index can say whether something is ' +
    'already known; it cannot say whether it is new, because any score built from a ' +
    'search tracks distance to the database. What is measured instead, in the detail ' +
    'below: this pile is not artefact, holds no abundant hidden genome, and sits at a ' +
    'median depth of one read, which is below every instrument that needs a second copy.';
  p.append(note);
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
  const SUMMARY = ['What we screened for', 'What one reference index can name',
                   'Known viruses: the reference, or a drifted strain',
                   'What nothing could name'];
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
    'Every measurement the summary above is read from, one section each, closed by ' +
    'default. Nothing here is a different number from the summary: the summary is ' +
    'drawn from these same files.'));
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
          f"{b['scatter'].get('known', {}).get('n_drawn', 0)}")
    if size > 16_000_000:
        raise SystemExit(f"{size} bytes is over the 16 MB page limit; lower "
                         f"--scatter-per-pile or --table-rows")
    print(f"built in {time.perf_counter() - t0:.1f}s "
          f"({platform.python_version()})", flush=True)


if __name__ == "__main__":
    main()
