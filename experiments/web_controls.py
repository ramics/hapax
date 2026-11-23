"""Resolve the page's sensitivity controls against one run's own files.

SPEC 13.3. The page is a playground because the report is un-thresholded by
construction: every per-read and per-contig row already carries identity, alignment
length, E-value, bitscore, coverage, depth and occupancy, so moving a control
re-filters and re-sorts what was computed and never recomputes it.

WHY THIS IS A SEPARATE PASS AND NOT INLINE IN THE PAGE. A slider is a way to take a
reading, so it needs the same provenance as any other reading this project ships. Three
things therefore happen here rather than in the renderer:

  - **The range comes from the data, never from the declaration.** `results/web/controls.json`
    says which quantity a control moves and what travels with it; every bound and every
    measured point is read from the run directory. A range typed into a file is a range
    with no provenance, and a slider that runs past where the measurement exists is how a
    clean-looking number gets taken at a setting nobody measured.
  - **The shipped default is checked against the measured grid.** A discrete control whose
    default is not one of the points the run actually measured means the page would mark
    the wrong setting as shipped, which is worse than having no marker: the reader would
    cite a threshold the run never used. That is a refusal, not a warning.
  - **Absent is absent, and it is not zero.** A control whose inputs this run does not
    hold resolves with `absent` set to the reason, and the page renders it as a stated
    constraint. A control whose inputs are present but inconsistent with the declaration
    fails the build. The shape this rules out is the one that has cost this project twice:
    a default where a refusal belongs.

The three sweeps the cascade already writes are the reason rule 2 exists. Kraken2
confidence was run at 0, 0.05, 0.1 and 0.2 and nowhere between, and one database's removal
of the same library moves from 42.6% to 3.5% across that gap, so the control snaps to the
measured points. Interpolating it would invent a factor of twelve.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTROLS_FILE = ROOT / "results/web/controls.json"

#: Which operating-point field each sweep in `tier_table.json` is a sweep of. This is a
#: property of the run format, so it lives beside the code that reads it rather than in
#: the declaration file. A sweep this does not know is a refusal: without the mapping the
#: page cannot say which of the measured points is the one the run shipped, and a sweep
#: drawn with no shipped marker reads as though every point were equally the result.
SWEEP_OPERATING_FIELD = {
    "entropy_sweep": "entropy_min",
    "mlen_sweep": "min_match_length_nt",
    "confidence_sweep": "k2_confidence",
    "bitscore_sweep": "diamond_bitscore",
}

REQUIRED_KEYS = ("id", "label", "kind", "unit", "placement", "panels", "arms")
KINDS = ("discrete", "continuous")
PLACEMENTS = ("headline", "detail")


def _fail(problems: list[str]) -> None:
    if problems:
        raise SystemExit(f"{CONTROLS_FILE}: the control declaration is refused:\n  - " +
                         "\n  - ".join(problems))


def _finite(xs):
    return [x for x in xs if isinstance(x, (int, float)) and math.isfinite(x)]


def _bounds(values, step, shipped):
    """A continuous control's range, from the values the run actually holds.

    The range is widened to include the shipped default where the data does not reach it,
    and says that it was: a default outside its own slider would be unreachable, and
    silently moving the default instead would mean the page marked a setting as shipped
    that the run never used.
    """
    vs = _finite(values)
    if not vs:
        return None
    lo, hi = min(vs), max(vs)
    widened = False
    if shipped is not None and math.isfinite(shipped):
        if shipped < lo:
            lo, widened = shipped, True
        if shipped > hi:
            hi, widened = shipped, True
    if hi <= lo:
        hi = lo + step
    return {"min": round(lo, 6), "max": round(hi, 6), "step": step,
            "n_values": len(vs), "widened_to_include_shipped": widened}


def _decades(evalues):
    """The E-value decades this run's own hits span, strongest first.

    Decades and not a smooth axis: an E-value floor is quoted, cited and argued about in
    decades, and the interesting thing about this control is not where between 1e-5 and
    1e-6 it sits but that the withdrawn reading of a 92.3% hit rate was taken with no
    floor at all.
    """
    vs = [v for v in _finite(evalues) if v > 0]
    if not vs:
        return None
    lo = max(-300, math.floor(math.log10(min(vs))))
    hi = min(2, math.ceil(math.log10(max(vs))))
    return [{"value": float(f"1e{e}"), "label": f"1e{e}" if e else "1 (no floor)"}
            for e in range(hi, lo - 1, -1)]


def _spec_rows(bundle: dict):
    S = bundle.get("spectrum")
    if not S:
        return None
    return list(S.get("on_axis") or []) + list(S.get("cannot_read") or [])


def _resolve_one(c: dict, rundir: Path, bundle: dict, search: dict,
                 problems: list[str]) -> dict:
    """One control, with its range, its measured points and its arms, or its absence."""
    out = {k: c[k] for k in ("id", "label", "kind", "unit", "placement", "panels",
                             "arms")}
    for k in ("note", "why_discrete", "shipped_basis", "absent_means"):
        if c.get(k):
            out[k] = c[k]

    missing = [f for f in (c.get("requires") or []) if not (rundir / f).exists()]
    if missing:
        out["absent"] = (f"this run holds no {', '.join(missing)}" +
                         (f". {c['absent_means']}" if c.get("absent_means") else ""))
        return out
    for key in c.get("requires_bundle") or []:
        if bundle.get(key) is None:
            why = (((bundle.get("page") or {}).get("inputs") or {})
                   .get(key) or {}).get("not_run")
            out["absent"] = (f"the {key} stage was not run on this library" +
                             (f": {why}" if why else ""))
            return out

    cid = c["id"]
    tt = bundle.get("tier_table") or {}
    op = tt.get("operating_point") or {}

    # A quantity this run's files do not carry AT ALL is an absence, and it renders as a
    # stated constraint: the run did not measure it, which is a different thing from a
    # declaration that disagrees with the data. What is a refusal is a measured quantity
    # whose SHIPPED point cannot be resolved - then the page would mark a setting as
    # shipped that the run never used, and a reader would cite it.
    if cid == "evalue_floor":
        pts = _decades([r[0] for r in search["all"]] + [r[2] for r in search["all"]])
        if pts is None:
            out["absent"] = ("reads.jsonl is present but no read in it has a hit with an "
                             "E-value, so there is no floor to move: that is a census of "
                             "this run and not a threshold")
            return out
        out["points"] = pts
        out["shipped"] = bundle["significant_evalue"]
    elif cid == "min_bitscore":
        bs = [r[1] for r in search["all"]] + [r[3] for r in search["all"]]
        out["bounds"] = _bounds(bs, 1, op.get("diamond_bitscore"))
        out["shipped"] = op.get("diamond_bitscore")
        if out["bounds"] is None:
            out.pop("bounds")
            out["absent"] = ("no read row in this run carries a bitscore, so there is "
                             "nothing for this floor to move over")
            return out
        if out["shipped"] is None:
            problems.append("min_bitscore: tier_table.json is present but its operating "
                            "point has no diamond_bitscore, so the page cannot say which "
                            "setting this run shipped")
    elif cid in ("search_unit", "spec_rank", "read_depth_stratum"):
        out["points"] = c["points"]
        out["shipped"] = c["shipped"]
    elif cid in ("spec_min_coverage", "spec_min_aln_aa", "min_contig_nt", "min_depth"):
        rows = _spec_rows(bundle)
        field, step = {"spec_min_coverage": ("cov", 0.01),
                       "spec_min_aln_aa": ("aln", 5),
                       "min_contig_nt": ("nt", 50),
                       "min_depth": ("depth", 1)}[cid]
        floors = (bundle["spectrum"].get("floors") or {})
        shipped = (floors.get("min_coverage") if cid == "spec_min_coverage" else
                   floors.get("min_aln_aa") if cid == "spec_min_aln_aa" else c["shipped"])
        if shipped is None:
            problems.append(f"{cid}: the spectrum is present but states no floor for it, "
                            f"so the page cannot mark the shipped setting")
        out["shipped"] = shipped
        out["bounds"] = _bounds([r.get(field) for r in rows], step, shipped)
        if out["bounds"] is None:
            out.pop("bounds")
            out["absent"] = (f"no contig on this run's spectrum carries {field!r}, so "
                             f"there is nothing for this floor to move over")
            return out
    elif cid == "min_occupancy":
        out["shipped"] = c["shipped"]
        out["bounds"] = _bounds([r[5] for r in search["all"]], 1, c["shipped"])
        if out["bounds"] is None:
            out.pop("bounds")
            out["absent"] = ("this run measured no cross-sample occupancy, so there is "
                             "no count to cut on. Occupancy comes from the recurrence "
                             "pass over the other dates of the same plant, and a run "
                             "without it has not searched them rather than found nothing")
            return out
    elif cid == "sweep_tier":
        swept = [c2 for c2 in (tt.get("curve") or []) if c2.get("threshold_sweep")]
        if not swept:
            out["absent"] = ("tier_table.json is present but no tier in it carries a "
                             "threshold sweep, so there is no measured grid to move on")
            return out
        out["points"] = [{"value": s["tier"], "label": s["tier"]} for s in swept]
        out["shipped"] = swept[0]["tier"]
        # The sweeps each tier measured, and which of their points the run shipped,
        # travel with the tier rather than being re-derived in the page.
        grid = {}
        for s in swept:
            name, sweep = next(iter(s["threshold_sweep"].items()))
            if name not in SWEEP_OPERATING_FIELD:
                problems.append(f"sweep_point: {s['tier']} carries a {name!r} this build "
                                f"cannot tie to an operating point; add it to "
                                f"SWEEP_OPERATING_FIELD rather than drawing a sweep with "
                                f"no shipped marker")
                continue
            opv = op.get(SWEEP_OPERATING_FIELD[name])
            keys = list(sweep.keys())
            at = next((k for k in keys if opv is not None and
                       abs(float(k) - float(opv)) < 1e-9), None)
            if at is None:
                problems.append(
                    f"sweep_point: {s['tier']}'s {name} was measured at "
                    f"{', '.join(keys)} but this run's operating point is "
                    f"{SWEEP_OPERATING_FIELD[name]}={opv!r}, which is not one of them. "
                    f"A default off the measured grid would mark a setting as shipped "
                    f"that the run never used")
            grid[s["tier"]] = {
                "quantity": name.replace("_sweep", ""),
                "operating_field": SWEEP_OPERATING_FIELD[name],
                "points": [{"value": k, "label": k, "removed": sweep[k]} for k in keys],
                "shipped": at,
                "span_points": (s.get("sweep_span_points") or {}).get(name),
                "alone": s.get("tier_alone_removed"),
                "marginal": s.get("marginal_removed"),
                "cumulative": (s.get("cumulative_removed_interval") or {}).get("point"),
            }
        out["grid"] = grid
    elif cid == "sweep_point":
        # Its points are the selected tier's, carried on sweep_tier's grid: the page reads
        # them from there so the two controls cannot disagree about what was measured.
        out["points_from_control"] = "sweep_tier"
    else:
        problems.append(f"{cid}: declared but this build has no resolver for it, so its "
                        f"range would be whatever the page invented")
        return out

    if out.get("kind") == "discrete" and not out.get("absent"):
        pts = out.get("points")
        if pts is not None:
            vals = [p["value"] for p in pts]
            if out.get("shipped") not in vals:
                problems.append(f"{cid}: the shipped default {out.get('shipped')!r} is "
                                f"not one of the measured points {vals}")
    return out


def resolve(rundir: Path, bundle: dict, search: dict) -> dict:
    """Every control, resolved against this run, plus the ones that cannot move.

    `search` is the compact per-read search table the page filters, which is also where
    the continuous bounds come from: the control and the figure read the same numbers, so
    a slider cannot run past the rows it filters.
    """
    doc = json.loads(CONTROLS_FILE.read_text())
    problems: list[str] = []
    seen = set()
    resolved = []
    for c in doc["controls"]:
        for k in REQUIRED_KEYS:
            if k not in c:
                problems.append(f"{c.get('id', '?')}: missing {k!r}")
        cid = c.get("id")
        if cid in seen:
            problems.append(f"{cid!r} is declared twice; a control id is what the URL "
                            f"state and every quoted setting are keyed on")
        seen.add(cid)
        if c.get("kind") not in KINDS:
            problems.append(f"{cid}: kind {c.get('kind')!r} is not one of {KINDS}")
        if c.get("placement") not in PLACEMENTS:
            problems.append(f"{cid}: placement {c.get('placement')!r} is not one of "
                            f"{PLACEMENTS}")
        if not c.get("panels"):
            problems.append(f"{cid}: no panel reads it, so moving it would do nothing")
        for arm in ("control", "null"):
            if not (c.get("arms") or {}).get(arm):
                problems.append(f"{cid}: arms.{arm} is not stated. A control that does "
                                f"not say where its {arm} arm comes from is the shape "
                                f"behind the readings this project has withdrawn")
        resolved.append(_resolve_one(c, rundir, bundle, search, problems))

    if not doc.get("fixed"):
        problems.append("the cannot-move table is empty; a control that cannot move "
                        "honestly has to be visible as a stated constraint rather than "
                        "absent from the page")
    for f in doc.get("fixed") or []:
        if not f.get("control") or not f.get("why"):
            problems.append(f"a cannot-move row is missing control or why: {f}")
    _fail(problems)
    return {"controls": resolved, "fixed": doc["fixed"]}
