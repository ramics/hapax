#!/usr/bin/env python3
"""Turn a production run directory into the tables its report quotes.

Separate from the run for the reason every `report_*.py` here is separate: the report's
numbers come from the committed JSON rather than from prose, and prose is what this project
has had to withdraw. Every table carries the spread across seeds, because five seeds with
no spread hides whether the five agree - and the spread is an observed range, never a
dispersion estimate (CLAUDE.md).

It is deliberately generic over the sample. The one per-sample thing it reads is the
library's own entry in `experiments/libraries.json` and `experiments/stage_references.json`,
so a second sample needs no edit here; and it prints what each stage did NOT do beside what
it did, because for a library with no truth panel the absences are most of the result.
"""
from __future__ import annotations

import argparse
import json
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_facts
import stage_refs

ROOT = Path(__file__).resolve().parents[1]


def spread(xs, pct=True, dp=2):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "-"
    k = 100 if pct else 1
    lo, hi = f"{min(xs) * k:.{dp}f}", f"{max(xs) * k:.{dp}f}"
    # A range printed over identical values is noise that reads like a measurement. It
    # fired on every index size and on the per-seed pair count, which are fixed by
    # construction across seeds: "0.000 (0.000-0.000)" invites the reader to look for a
    # spread in a quantity that cannot have one. Collapse to the single figure, which also
    # keeps a one-seed run and a five-seed run of a constant looking the same.
    if len(xs) == 1 or lo == hi:
        return lo
    return (f"{st.mean(xs) * k:.{dp}f} ({lo}-{hi})")


def load(path: Path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _ran(entry: dict, ran: dict, key: str) -> str:
    """One sample fact as the run used it, saying which of the two sources answered.

    A fact absent from the registry is not an unknown: a trim and a block are measured
    off the reads, and a library whose registry entry deliberately records neither still
    ran on a definite value. Printing the registry alone said "unmeasured" of a trim the
    pass had measured, which is the one direction of error this project treats as
    expensive, since it reads as a gap in the measurement rather than in the lookup.
    """
    stated, used = entry.get(key), ran.get(key)
    fmt = (lambda v: f"{v:,}" if isinstance(v, int) else str(v))
    if stated is not None:
        return f"{fmt(stated)}, stated" + (
            f" (the run recorded {fmt(used)})" if used is not None and used != stated
            else "")
    if used is not None:
        return f"{fmt(used)}, measured on this run"
    return "neither stated nor recorded"


def _shown(path: Path) -> str:
    """The path as the repo names it, or in full where it sits outside the repo.

    A run directory is not always under the checkout: it can be a work directory on an
    instance or a copy under a test's tmp_path, and `relative_to` raises on those rather
    than falling back, which turned one provenance line into a crash. The same fix is in
    `make_web._shown` for the same reason.
    """
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="recorded; inputs are fixed files")
    ap.add_argument("--library", required=True)
    ap.add_argument("--results", type=Path, default=None,
                    help="the run directory's parent; default results/prod_<library>")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    res = args.results or ROOT / f"results/prod_{args.library}"
    if not res.exists():
        raise SystemExit(f"{res} does not exist")
    found = sorted(res.glob("run_s*"))
    if not found:
        raise SystemExit(f"{res} holds no run_s* directory")
    # A DIRECTORY IS NOT A SEED. A run still in flight, or one fetched while in flight,
    # leaves a `run_sN` directory holding part of its output and no manifest. Folding
    # those into the seed count prints a census over seeds it has no data for: this
    # table read "0 failures across 5 seeds" and "0 of 5 complete" off four finished
    # seeds and one directory containing a single file, while the cascade table beside
    # it correctly said four. The two are separated and the partial ones are named, so
    # that a count is over the seeds that produced a manifest and nothing is averaged
    # over an absence.
    runs = [p for p in found if (p / "manifest.json").exists()]
    started = [p.name.removeprefix("run_") for p in found if p not in runs]
    if not runs:
        raise SystemExit(
            f"{res} holds {len(found)} run_s* directories and not one manifest.json, "
            "so no seed of this run has finished")

    entry = library_facts.library(args.library)
    # The parameters of the first seed's own manifest: what the run actually ran on.
    ran = ((load(runs[0] / "manifest.json") or {}).get("parameters") or {})
    truth = stage_refs.resolve(args.library, "truth_panel")
    cross = stage_refs.resolve(args.library, "crossplant_panel")

    L = []
    w = L.append
    w(f"# Tables — production run over {args.library}")
    w("")
    w(f"Generated by `experiments/report_prod.py` from `{_shown(res)}`. "
      f"Every number here is read out of a committed file in that directory; nothing is "
      f"typed.")
    w("")
    w("## The sample")
    w("")
    w("| field | value |")
    w("|---|---|")
    w(f"| library | {args.library} |")
    w(f"| label | {entry.get('label', '-')} |")
    w(f"| project | {entry.get('bioproject', '-')} |")
    w(f"| layout | {entry.get('layout', '-')} |")
    w(f"| pairs in the library | {entry.get('pairs', '-'):,} |"
      if isinstance(entry.get("pairs"), int) else "| pairs in the library | - |")
    w(f"| read length | {entry.get('read_length', '-')} |")
    # THE VALUE THE RUN USED, not the value the registry holds. For a library whose trim
    # and block are deliberately absent from the registry - because both are measured on
    # the instance, which is the point of the one-command path - these two rows read
    # "unmeasured" and "-" on a run that had measured them both and printed them in its
    # own operating point two sections down. The registry answers what was supplied; only
    # the run answers what ran, and where they differ the stated one wins, so both are
    # shown with the source named.
    w(f"| mate-2 trim used | {_ran(entry, ran, 'trim_r2')} |")
    w(f"| block stride | {_ran(entry, ran, 'block')} |")
    w(f"| truth panel | {truth['state']}"
      + (f": {truth.get('panel')}" if truth["state"] == "resolved" else "") + " |")
    w(f"| cross-plant panel | {cross['state']}"
      + (f": {cross.get('panel')}" if cross["state"] == "resolved" else "") + " |")
    w("")
    if truth["state"] == "absent":
        w(f"**No ground truth for this water.** {truth['reason']}")
        w("")

    # ---- what was measured off the FASTQ, beside what the registry states. Imo, on being
    # shown the hand-supplied ledger: "why are we not discovering the facts from the sample
    # itself?" Four of these come off the file; where the registry states a different value
    # the stated one is what ran, because stating it is the override, so both are printed
    # and the disagreement is visible rather than reconciled.
    meas = load(res / "measured_library.json")
    if meas:
        w("## Measured off the file, beside what the registry states")
        w("")
        w("| field | measured from the reads | registry | which one ran |")
        w("|---|---|---|---|")
        for k in ("pairs", "read_length", "layout", "block"):
            mv, ev = meas.get(k), entry.get(k)
            ran = "registry" if ev is not None else "measured"
            w(f"| {k} | {mv} | {ev if ev is not None else '-'} | {ran} |")
        w("")
        w("Not measurable from the file, so supplied with a basis: "
          + ", ".join(f"`{k}`" for k in (meas.get("not_measurable_from_the_file") or {}))
          + ".")
        w("")

    # ---- the trim, which for a new library is a measurement and not a lookup
    # Either filename: `prod_run.sh` states `--out .../trim_from_prefix.json` while the
    # script's own default keys the name by accession, so a run driven by the launcher and
    # a run driven by hand leave different names. The section this builds is the one that
    # records a REFUSAL, and a missing file would render as no section at all.
    dec = (load(res / "nontemplate_prefix/trim_from_prefix.json")
           or load(res / f"nontemplate_prefix/trim_from_prefix_{args.library}.json"))
    if dec:
        w("## The mate-2 non-template prefix, measured on this library")
        w("")
        w("| field | value |")
        w("|---|---|")
        w(f"| seeds | {dec.get('seeds')} |")
        w(f"| mate-2 5' mode per seed | {dec.get('modes')} |")
        w(f"| mode share per seed | {dec.get('mode_shares')} |")
        w(f"| mate-1 5' mode per seed | {dec.get('r1_modes')} |")
        w(f"| anchored mate-2 reads per seed | {dec.get('anchored_r2')} |")
        named = dec.get("trim_r2")
        w(f"| trim named by the measurement | "
          f"{'refused' if named is None else named} |")
        w(f"| trim the run used | {entry.get('trim_r2', 'unmeasured')}, stated |")
        # One refusal per seed, five times the same sentence with one number changed, is a
        # wall of prose that hides how many distinct reasons fired. The count is what a
        # reader needs and the sentence is what they need once.
        probs = dec.get("problems") or []
        if probs:
            seen = {}
            for q in probs:
                key = "".join(c for c in q if not c.isdigit())
                seen[key] = seen.get(key, 0) + 1
            for key, n in seen.items():
                first = next(q for q in probs
                             if "".join(c for c in q if not c.isdigit()) == key)
                w(f"| refused, {n} of {len(probs)} seeds | {first} |")
        else:
            w("| refusals | none |")
        w("")
        w(f"Rule applied: {dec.get('rule')}")
        w("")
        w(f"Bias: {dec.get('bias')}")
        w("")

    # ---- the assembly, with its own argv and version
    asm = load(res / "assembly/assembly_baseline.json")
    if asm:
        w("## The assembly (whole library, committed baseline arm)")
        w("")
        w("| field | value |")
        w("|---|---|")
        for k in ("arm", "flags", "threads", "seconds", "contigs", "bases", "n50", "max",
                  "ge500", "ge1000", "ge3000", "peak_rss_bytes"):
            if k in asm:
                v = asm[k]
                w(f"| {k} | {v:,} |" if isinstance(v, int) else f"| {k} | {v} |")
        w("")
    man = load(res / "assembly/manifest_s1.json")
    if man:
        w(f"MEGAHIT `{man.get('megahit_version', '?')}`, "
          f"minimap2 `{man.get('minimap2_version', '?')}`, "
          f"reads sha256 `{str(man.get('reads_sha256'))[:16]}`.")
        w("")

    # ---- the contig split: what identification could and could not name
    cid = load(res / "contig_identification.json")
    if cid:
        w("## The contig split (identification, never assembly)")
        w("")
        w("| field | value |")
        w("|---|---|")
        for k, v in cid.items():
            if isinstance(v, bool) or not isinstance(v, (int, float, str)):
                continue
            if k.startswith("_"):
                continue
            # A fraction is printed to four places. The unrounded float carries sixteen
            # digits of which about three are meaningful at n = 323,430 contigs, and a
            # table that shows all sixteen invites a reader to quote them.
            if isinstance(v, float):
                w(f"| {k} | {v:.4f} |")
            else:
                w(f"| {k} | {v:,} |" if isinstance(v, int) else f"| {k} | {v} |")
        w("")

    # ---- the cascade, across seeds, with the operating point it was read at
    tiers = [(p.name, load(p / "tiers.json"), load(p / "tier_table.json")) for p in runs]
    have = [(n, t, tt) for n, t, tt in tiers if t]
    if have:
        op = next((tt.get("operating_point") for _, _, tt in have if tt), {})
        w("## The cascade")
        w("")
        w(f"Operating point: `{json.dumps(op)}`. Percentages are of read **pairs**, the "
          f"occurrence-weighted unit. Spread is the observed range over "
          f"{len(have)} seed{'' if len(have) == 1 else 's'} and is not a dispersion "
          f"estimate.")
        w("")
        w("| quantity | mean (min-max) |")
        w("|---|---|")
        w(f"| pairs read per seed | {spread([t['pairs'] for _, t, _ in have], False, 0)} |")
        w(f"| union removed | {spread([t['union_removed'] for _, t, _ in have])} |")
        w(f"| residue | {spread([1 - t['union_removed'] for _, t, _ in have])} |")
        if all(t.get("tier3_pairs") is not None and t.get("pairs") for _, t, _ in have):
            w(f"| tier 3 (unplaceable reads) | "
              f"{spread([t['tier3_pairs'] / t['pairs'] for _, t, _ in have])} |")
        w("")
        # THE CASCADE CURVE LIVES IN `curve`, NOT IN `tiers`, and reading it from the wrong
        # key is why this table rendered as a column of dashes on both libraries it has ever
        # been run on. `tiers[name]` is each tier's own internal breakdown - the artefact
        # tier's entropy sweep, a database's confidence sweep - while the per-tier removal
        # a reader wants (alone, marginal, cumulative) is one row per tier in `curve`, and
        # the index sizes are a top-level dict keyed by tier. A table of dashes renders as
        # cleanly as a table of numbers, which is the same failure this repo keeps paying
        # for: an empty output and a step that never ran look identical.
        curves = [{r["tier"]: r for r in (t.get("curve") or [])} for _, t, _ in have]
        sizes = [t.get("index_bytes") or {} for _, t, _ in have]
        order = list(dict.fromkeys(k for c in curves for k in c))
        if order:
            w("| tier | alone | marginal | cumulative | index GB |")
            w("|---|---|---|---|---|")
            for tier in order:
                rows = [c.get(tier) or {} for c in curves]
                gb = [(s.get(tier) or (r.get("index_bytes") or 0)) / 1e9
                      for s, r in zip(sizes, rows)]
                gb_cell = spread(gb, False, 3) if any(gb) else "-"
                w(f"| {tier} | {spread([r.get('tier_alone_removed') for r in rows])} | "
                  f"{spread([r.get('marginal_removed') for r in rows])} | "
                  f"{spread([r.get('cumulative_removed') for r in rows])} | {gb_cell} |")
            w("")
            # A tier that is measured and never subtracted has no row in the curve, and its
            # absence from the table above has to be stated or it reads as a tier that was
            # dropped. The unidentifiable-contig tier is the one this applies to: SPEC's
            # subtrahend is identification, so matching an unnamed contig is recorded and
            # does not remove the read. Subtracting it inflated removal to 85.9% against
            # the true 77.0% once already.
            observed = [k for _, t, _ in have for k in (t.get("tiers") or {})
                        if k not in order]
            if observed:
                w(f"Measured and not subtracted, so carrying no row above: "
                  f"`{'`, `'.join(dict.fromkeys(observed))}`. Matching sequence nobody has "
                  f"named is recorded against the read and does not remove it.")
                w("")

    # ---- what every stage did, and what it refused: for this sample that IS the result
    # THE STAGE CENSUS IS COLLAPSED ACROSS SEEDS, AND WHY. Every stage has to appear, because
    # "ran and found nothing" and "never ran" look identical in an absence - that is the census
    # rule. But five seeds of a thirteen-stage run is 65 rows of which 65 were identical here,
    # and a table nobody reads to the end hides a disagreement as effectively as omitting it.
    # So the row is keyed by (stage, status, reason) with the seeds that produced it named:
    # agreement collapses to one line, and a seed that differs becomes its own line, which is
    # the only case worth looking at. The stage order is the run's own, never sorted.
    w("## What every stage did, with the seeds that agree collapsed")
    w("")
    w("| seeds | stage | status | reason |")
    w("|---|---|---|---|")
    census: dict[tuple, list[str]] = {}
    for p in runs:
        m = load(p / "manifest.json") or {}
        for s in m.get("stages", []):
            reason = (s.get("reason") or "").replace("|", "/").replace("\n", " ")
            if len(reason) > 220:
                reason = reason[:217] + "..."
            census.setdefault((s.get("stage"), s.get("status"), reason), []).append(p.name)
    n = len(runs)
    for (stage, status, reason), seeds in census.items():
        cell = (f"all {n}" if len(seeds) == n
                else ", ".join(s.removeprefix("run_") for s in seeds))
        w(f"| {cell} | {stage} | {status} | {reason} |")
    w("")
    if all(len(v) == n for v in census.values()):
        w(f"Every stage returned the same status and reason on all {n} seeds, so no row here "
          f"names a seed.")
        w("")
    if started:
        w(f"Seeds with a directory and no manifest, so counted nowhere above: "
          f"**{', '.join(started)}**. Those seeds had not finished when this table was "
          f"written.")
        w("")
    done = [load(p / "manifest.json") or {} for p in runs]
    # THREE COUNTS, NOT ONE, BECAUSE THEY MEAN DIFFERENT THINGS AND ONE HID THE OTHERS
    # TWICE, IN BOTH DIRECTIONS.
    #
    # `complete` is false whenever any stage produced no output, and a stage deliberately
    # left out of a production pass (SPEC 13.4) is one of those, so a correct run reported
    # "complete runs: 0 of 5" and read as five failures. The fix for that printed the
    # failure count alone and dropped completeness entirely, which is the opposite error:
    # a run whose stages all exited cleanly while writing nothing then reads as fine, and
    # PLAN.md 1.5 is that fault. A stage exiting zero is not a filled report section.
    #
    # So the line carries all three: how many stages FAILED, how many produced NO OUTPUT
    # without failing, and the plain completeness figure with the deliberate omissions
    # named beside it. A reader can tell a refusal from an omission from an empty success.
    # Judging a run on its outputs rather than its exit codes is PLAN.md step 4, where the
    # report build itself starts refusing; this line is the visible version in the interim.
    failed = sum(len([s for s in (m.get("stages") or []) if s.get("status") == "failed"])
                 for m in done)
    unsel = sorted({s.get("stage") for m in done for s in (m.get("stages") or [])
                    if s.get("status") == "skipped"})
    empty = sum(len([s for s in (m.get("stages") or [])
                     if s.get("status") not in ("failed", "skipped") and not s.get("outputs")])
                for m in done)
    complete = sum(1 for m in done if m.get("complete"))
    w(f"Stages that failed: **{failed}** across {len(done)} seeds; stages that exited "
      f"without failing and wrote no output: **{empty}**; runs marked `complete`: "
      f"**{complete} of {len(done)}**. Completeness is low by construction here because "
      f"{len(unsel)} stages are not part of a production pass "
      f"(`{'`, `'.join(str(u) for u in unsel)}`), which is the manifest recording a "
      f"deliberate omission and not a failure, and the figure is printed anyway so that "
      f"an empty stage cannot hide behind a clean exit code. Wall clock per seed: "
      f"{[round(m.get('wall_seconds', 0) / 60, 1) for m in done]} minutes.")
    w("")

    out = args.out or res / "TABLES.md"
    out.write_text("\n".join(L) + "\n")
    print(f"written {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
