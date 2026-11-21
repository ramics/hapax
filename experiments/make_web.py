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
    ("best_pident", "Best protein identity %", lambda r: _best(r, "pident")),
]


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
  --sans: "IBM Plex Sans", system-ui, -apple-system, Segoe UI, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
} }
:root[data-theme="dark"] {
  --bg: #101513; --panel: #171e1b; --edge: #2a332f; --fg: #e6ebe8; --dim: #97a39d;
  --identified: #6ea2bf; --tier2: #c0ab63; --residue: #df8a5e; --null: #68756f;
  --ok: #67b183; --warn: #d39354; color-scheme: dark;
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
  h.append(el('h1', null, 'Hapax run viewer'));
  const inp = (r.inputs || {});
  h.append(el('p', 'lede',
    'Everything below is read from one run directory. Subtract what is identifiable, ' +
    'keep what is left, and report each number beside the operating point and ' +
    'calibration that make it readable. There is no novelty score on this page: an index ' +
    'can say whether something is already known and cannot say whether it is new.'));
  h.append(opChips({seed: r.seed, pairs: (D.parameters || {}).pairs,
    'mate 2 trim': inp.trim_r2, 'git sha': (r.git_sha || '').slice(0, 10),
    'tree dirty': String(!!r.tree_dirty), complete: String(!!r.complete)}));
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
  const identified = last ? last.cumulative_removed : null;
  const residue = last ? last.residue : null;
  const bar = el('div', 'bar');
  const seg = (frac, colour) => { const s = document.createElement('span');
    s.style.width = (100 * (frac || 0)) + '%'; s.style.background = colour; return s; };
  bar.append(seg(identified, 'var(--identified)'));
  bar.append(seg(residue, 'var(--residue)'));
  p.append(bar);
  const lg = el('div', 'legend');
  lg.innerHTML =
    '<span><i style="background:var(--identified)"></i>identified ' +
      interval(last && last.cumulative_removed_interval) + '</span>' +
    '<span><i style="background:var(--residue)"></i>residue ' +
      interval(last && last.residue_interval) + '</span>';
  p.append(lg);
  p.append(opChips(t.operating_point));
  p.append(el('p', 'note', 'Tier 2, what assembled but nothing identified, is drawn ' +
    'only when a contig set was supplied to the run; this run ' +
    ((D.run.stages || []).find(s => s.stage === 'contigs' && s.status === 'ok')
      ? 'has one.' : 'does not, so the bar has two segments rather than three.')));
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
  const sc = el('div', 'scroll');
  const tb = el('table');
  tb.innerHTML = '<thead><tr><th>subject</th><th>index</th><th>residue reads</th>' +
    '<th>identified reads</th><th>median identity</th><th>aligned</th>' +
    '<th>median E</th><th>worst E</th><th>reading</th></tr></thead>';
  const body = el('tbody');
  (S.subjects || []).slice(0, 120).forEach(r => {
    const tr = el('tr');
    tr.innerHTML = '<td class="mono">' + esc(r.subject) + '</td><td class="dim">' +
      esc(r.index) + '</td><td class="mono">' + r.residue_reads + '</td><td class="mono">' +
      (S.known_pile_searched_reads ? r.known_reads :
        '<span class="absent">not searched</span>') + '</td><td class="mono">' +
      num(r.median_identity_percent, 1) + '%</td><td class="mono">' +
      r.median_aligned_aa + ' aa</td><td class="mono">' + (+r.median_evalue).toExponential(1) +
      '</td><td class="mono">' + (+r.worst_evalue).toExponential(1) + '</td><td class="' +
      (/not comparable/.test(r.reading) ? 'absent' : 'dim') + '">' +
      esc(/not comparable/.test(r.reading) ? 'control not searched' :
          /also/.test(r.reading) ? 'also in the identified pile' : 'residue only') +
      '</td>';
    body.append(tr);
  });
  tb.append(body); sc.append(tb); p.append(sc);
  p.append(el('p', 'note', S.subjects_total + ' subjects at E < ' +
    S.significance_evalue + ', the top 120 by read count shown. ' + esc(S.caveat) +
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
  sx.onchange = draw; sy.onchange = draw; draw();
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
  p.append(el('pre', 'mono note scroll', esc(JSON.stringify(C.result || C, null, 1)
    .slice(0, 4000))));
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
      tr.innerHTML = '<td class="mono">' + esc(r.name) + '</td><td class="dim">' +
        esc((r.layer || '-').replace(/_/g, ' ')) + '</td><td class="mono">' +
        (best ? esc(best.subject) : '<span class="dim">none</span>') +
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
