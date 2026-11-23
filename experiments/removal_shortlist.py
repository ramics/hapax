#!/usr/bin/env python3
"""Rank the pile tier 1 removed by how much of each contig its best hit EXPLAINS.

WHAT THIS PROMOTES, AND WHY IT IS A STAGE AND NOT A ONE-OFF. `results/holdout_cascade/`
measured the one mechanism that catches a genuinely novel genome being subtracted as known,
on the one positive this project has: a virus deposited after the library that contains it,
withheld from every index, which assembles out of its own library as a single 3,990 nt contig
covering 96% of itself, hits the held-out reference at 66.5% protein identity over 620 aa
with a bitscore of 840 against the protein tier's threshold of 50, and is removed in tier 1.
**The pipeline would have deleted it.** That measurement existed as a report over one
directory. This script is the same arithmetic over any run's `removal_evidence.tsv`, so a
new sample gets the shortlist without anyone re-deriving it.

THE AXIS IS COVERAGE, AND IDENTITY IS THE AXIS THAT FAILS. Of 22,260 removed background
contigs, 7,601 have a significant protein hit, and their median identity to the nearest
known protein is 40.0%. The withheld genome at 66.5% is the **93.46th percentile**: it is
MORE similar to a known protein than 93% of what the cascade removed, so a rule of "flag
the low-identity removals" reaches 7,104 other contigs before it. That rule is dead and
measured. Adding the fraction of the contig the best hit explains turns the same pile into a
cell of **69 contigs, 0.91% of the removals**, which `tests/test_removal_shortlist.py`
reproduces from the committed rows.

AND THE DIRECTION IS PINNED ON KNOWN MATERIAL FROM THE SAME WATER, which is the half that
makes the cell readable rather than arbitrary: of the leak arm's 34 contigs - sequence from
this same library whose relatives are in the index - **not one** sits below 80% identity at
either level, and none enters the candidate's cell. So this stage always draws its control
arm in the same pass by the same code (SPEC 6.1): near-reference removals from this run's
own pile. A band with no control arm beside it is reported as uncontrolled, by name.

WHAT IT REFUSES TO BE. Not a flag and not a score. Inside its own cell the candidate ranks
fourth by bitscore at nucleotide level and second at protein, so a box that fired on the top
row would have missed it: what is defensible is the CELL WITH A COUNT ON IT, and this script
writes the whole joint table beside the band so the band can be moved at output time
(SPEC 13.3) rather than being the only thing a reader is shown. Ordering within the cell is
an open question and the output says so instead of implying the top row matters.

WHAT IT DOES NOT COMPUTE. It runs no search. It reads the evidence the cascade already
wrote and then discarded - `classify_contigs.py` emits `removal_evidence.tsv` in input order
with no score on it - so it cannot produce its own shuffled null, and it will not pretend to:
pass `--shuffled-evidence` to supply one, or `--without-null` to record its absence by name.
The citations for the null (0 of 22,110 shuffled contigs significant on the holdout library,
0 of 11,361 on Berlin) travel in the output as citations, labelled as such.

ROUTES ARE NOT POOLED. Only the protein route carries an identity AND a coverage on every
row, so the band is a protein-level band. The nucleotide side of the cascade records
matching bases for the ribosomal route and a taxid with a confidence for Kraken2: three
kinds of number, and averaging them would be comparing two instruments. The nucleotide
equivalent of this cell is measurable (137 contigs, 6.64%, in `results/holdout_cascade/`)
but not from this file, and the output says which.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import stage_refs
from loud import Preflight
from provenance import provenance

ROOT = Path(__file__).resolve().parents[1]

#: The joint grid. Fixed here rather than taken from the command line because these are the
#: cells `results/holdout_cascade/` published and the band resolves to one of them; a reader
#: moving the band at output time moves it over the same grid.
COVERAGE_BANDS = [(0.70, 1.01), (0.50, 0.70), (0.25, 0.50), (0.00, 0.25)]
IDENTITY_BANDS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100.01)]


def read_evidence(path: Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = {"contig", "route", "pident", "evalue", "bitscore", "query_coverage"}
    if rows and not need <= set(rows[0]):
        raise SystemExit(f"{path} is missing columns {sorted(need - set(rows[0]))}: this "
                         f"is not a removal_evidence.tsv written by classify_contigs.py")
    return rows


def best_protein_rows(rows: list[dict], evalue: float) -> tuple[dict, dict]:
    """One row per contig: its best protein hit by bitscore, and a census of what was cut.

    THE CENSUS IS NOT DECORATION. "No significant protein hit" and "no protein row at all"
    are different facts about a contig - the first says the protein tier looked and found
    nothing worth reporting, the second that this contig was removed by a nucleotide tier
    and has no identity to place on this axis at all. A stage that silently dropped both
    into the same denominator would be quoting a share of a pile it had redefined.
    """
    census = {"rows": len(rows), "contigs": len({r["contig"] for r in rows}),
              "protein_rows": 0, "unrecorded_route_rows": 0,
              "protein_rows_without_coverage": 0,
              "contigs_with_a_protein_row": 0,
              "contigs_with_a_significant_protein_row": 0,
              "contigs_removed_by_other_routes_only": 0,
              "routes": {}}
    best: dict[str, dict] = {}
    for r in rows:
        census["routes"][r["route"]] = census["routes"].get(r["route"], 0) + 1
        if r["route"] == "unrecorded":
            census["unrecorded_route_rows"] += 1
        if r["route"] != "protein":
            continue
        census["protein_rows"] += 1
        if not r.get("query_coverage"):
            # Coverage IS the axis. A protein row with a blank coverage cannot be scored
            # and must not be read as 0.0, which would put it in the lowest cell.
            census["protein_rows_without_coverage"] += 1
            continue
        try:
            e, b = float(r["evalue"]), float(r["bitscore"])
            cov, pid = float(r["query_coverage"]), float(r["pident"])
        except (TypeError, ValueError):
            census["protein_rows_without_coverage"] += 1
            continue
        if e >= evalue:
            continue
        prev = best.get(r["contig"])
        if prev is None or b > prev["bitscore"]:
            best[r["contig"]] = {
                "contig": r["contig"],
                "contig_nt": int(r["contig_nt"]) if r.get("contig_nt") else None,
                "index": r.get("index", ""), "subject": r.get("subject_or_taxon", ""),
                "identity_percent": pid, "aln_aa": r.get("aln_aa", ""),
                "evalue": e, "bitscore": b, "explains_fraction": cov,
                "tier": r.get("tier", "")}
    prot_contigs = {r["contig"] for r in rows if r["route"] == "protein"}
    census["contigs_with_a_protein_row"] = len(prot_contigs)
    census["contigs_with_a_significant_protein_row"] = len(best)
    census["contigs_removed_by_other_routes_only"] = census["contigs"] - len(prot_contigs)
    return best, census


def cell_of(row: dict) -> tuple[tuple[float, float], tuple[int, int]] | None:
    cb = next((p for p in COVERAGE_BANDS if p[0] <= row["explains_fraction"] < p[1]), None)
    ib = next((p for p in IDENTITY_BANDS if p[0] <= row["identity_percent"] < p[1]), None)
    return (cb, ib) if cb and ib else None


def joint_table(best: dict) -> list[dict]:
    """Every cell with its count. The whole table ships, not only the band's cell."""
    out = []
    n = len(best) or 1
    for cb in COVERAGE_BANDS:
        for ib in IDENTITY_BANDS:
            k = sum(1 for r in best.values()
                    if cb[0] <= r["explains_fraction"] < cb[1]
                    and ib[0] <= r["identity_percent"] < ib[1])
            out.append({"coverage_lo": cb[0], "coverage_hi": min(cb[1], 1.0),
                        "identity_lo": ib[0], "identity_hi": min(ib[1], 100),
                        "contigs": k, "percent_of_pile": round(100 * k / n, 4)})
    return out


def in_band(best: dict, cov: tuple[float, float], ident: tuple[float, float]) -> list[dict]:
    return [r for r in best.values()
            if cov[0] <= r["explains_fraction"] < cov[1]
            and ident[0] <= r["identity_percent"] < ident[1]]


def main() -> int:
    t0 = time.perf_counter()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="unused by the arithmetic - the input is a committed table - but "
                         "recorded, because every result in this repo carries its seed")
    ap.add_argument("--evidence", type=Path, required=True,
                    help="removal_evidence.tsv from classify_contigs.py: one row per route "
                         "per removed contig. REQUIRED and never defaulted: a default here "
                         "would be one sample's removals")
    ap.add_argument("--shuffled-evidence", type=Path, default=None,
                    help="the same table computed over shuffled contigs, which is this "
                         "stage's null. Absent, --without-null must say so by name")
    ap.add_argument("--without-null", action="store_true",
                    help="record ON PURPOSE that no null was supplied in this run, with "
                         "the citations for the nulls measured elsewhere. Without this, a "
                         "missing --shuffled-evidence refuses to start, because a band "
                         "quoted with a silently absent null is the shape this repo has "
                         "withdrawn figures over")
    ap.add_argument("--band-coverage", type=float, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="override the shipped band's coverage edges; the shipped values "
                         "come from experiments/stage_references.json with their provenance")
    ap.add_argument("--band-identity", type=float, nargs=2, default=None,
                    metavar=("LO", "HI"), help="override the shipped band's identity edges")
    ap.add_argument("--significance-evalue", type=float, default=None,
                    help="the significance floor; shipped at 1e-5, because the weaker band "
                         "describes the database rather than the contig")
    ap.add_argument("--near-reference-identity", type=float, default=None,
                    help="identity at or above which a removal counts as near-reference "
                         "material from this water: the control arm. Shipped at 97%%, the "
                         "top band of SPEC 6.1")
    ap.add_argument("--references", type=Path, default=None,
                    help="stage_references.json, for tests")
    ap.add_argument("--out", type=Path, required=True, help="output directory")
    ap.add_argument("--preflight", action="store_true")
    args = ap.parse_args()

    pre = Preflight("removal_shortlist.py")
    pre.file(args.evidence, "removal evidence")
    if args.shuffled_evidence is not None:
        pre.file(args.shuffled_evidence, "shuffled-contig evidence (the null)")
    else:
        pre.require(args.without_null,
                    "no --shuffled-evidence and no --without-null: this stage will not "
                    "report a band with an unstated null. Supply the null, or drop it by "
                    "name")
    try:
        band = stage_refs.band(stage_refs.load(args.references))
    except (stage_refs.MissingReference, KeyError) as e:
        band = None
        pre.require(False, f"the shipped band could not be resolved: {e}")
    pre.done(args.preflight)

    cov = tuple(args.band_coverage or band["coverage"])
    ident = tuple(args.band_identity or band["identity_percent"])
    evalue = args.significance_evalue or band["significance_evalue"]
    near = args.near_reference_identity or band["near_reference_identity_percent"]
    args.out.mkdir(parents=True, exist_ok=True)

    rows = read_evidence(args.evidence)
    if not rows:
        raise SystemExit(f"{args.evidence} has no rows. An empty removal table is not an "
                         f"empty shortlist: it means the cascade removed nothing or wrote "
                         f"nothing, and this stage cannot tell which")
    best, census = best_protein_rows(rows, evalue)
    if census["protein_rows"] == 0:
        # Refusal, not a zero. With no protein row anywhere, either the protein tier did
        # not run or it removed nothing, and this file cannot distinguish those. The
        # difference decides whether an empty shortlist means anything.
        raise SystemExit(
            f"{args.evidence} carries no protein-route rows ({census['rows']} rows over "
            f"{census['contigs']} contigs, routes {census['routes']}). The explained "
            f"fraction is a protein-level axis, and this file cannot say whether the "
            f"protein tier ran and found nothing or never ran. Re-run the cascade with its "
            f"protein indexes, or leave this stage out of --stages.")
    if not best:
        shortlist: list[dict] = []
    else:
        shortlist = sorted(in_band(best, cov, ident), key=lambda r: r["explains_fraction"])

    # The control arm, drawn in the same pass by the same code: near-reference removals from
    # this run's own pile. The holdout's leak arm is the measured version of this - 34 known
    # contigs from the same water, none of them in the candidate's cell - and a run with no
    # near-reference removals has no such arm, which is reported rather than left blank.
    control = [r for r in best.values() if r["identity_percent"] >= near]
    control_in_band = in_band({r["contig"]: r for r in control}, cov, ident)
    identities = sorted(r["identity_percent"] for r in best.values())

    out = {
        "what_this_is": ("the pile tier 1 removed, ranked by how much of each contig its "
                         "best protein hit explains. A cell with a count on it, never a "
                         "flag and never a score"),
        "spec": "SPEC.md 5.1",
        "evidence_file": str(args.evidence),
        "operating_point": {
            "significance_evalue": evalue,
            "band_coverage": list(cov), "band_identity_percent": list(ident),
            "near_reference_identity_percent": near,
            "band_is_shipped_default": (args.band_coverage is None
                                        and args.band_identity is None),
            "band_provenance": band,
        },
        "pile": {
            "contigs_scored": len(best),
            "a_contig_is_scored_when": (f"it has a protein-route row with an E-value below "
                                        f"{evalue} and a recorded query coverage"),
            "identity_median_percent": round(st.median(identities), 3) if identities else None,
            "identity_quartiles_percent": (
                [round(q, 3) for q in st.quantiles(identities, n=4)]
                if len(identities) > 3 else None),
            "explains_fraction_median": round(
                st.median([r["explains_fraction"] for r in best.values()]), 4
            ) if best else None,
        },
        "census": census,
        "joint_table": joint_table(best),
        "band": {
            "contigs": len(shortlist),
            "percent_of_scored_pile": round(100 * len(shortlist) / len(best), 4) if best else None,
            "ordering": "ascending explained fraction",
            "ordering_is_validated": False,
            "ordering_note": band["ordering_within_the_cell_is_open"],
        },
        "control_arm": {
            "state": "measured" if control else "not measured",
            "what_it_is": (f"removals at or above {near}% identity: near-reference material "
                           f"from this same library, the arm that says which cells known "
                           f"sequence occupies"),
            "contigs": len(control),
            "contigs_inside_the_band": len(control_in_band),
            "direction_measured_on_known_material": (
                "On the holdout library's 34 leak contigs - known sequence from the same "
                "water - not one sits below 80% identity at either level and none enters "
                "the band. So crowding inside the band is other background, not "
                "recognisable material (results/holdout_cascade/)."),
            "why_not_measured": (None if control else
                                 f"no removal in this run reached {near}% identity, so this "
                                 f"run has no near-reference arm and the band below is "
                                 f"uncontrolled here"),
        },
        "null": ({"state": "measured", "file": str(args.shuffled_evidence)}
                 if args.shuffled_evidence else
                 {"state": "not supplied: --without-null",
                  "why": band["null_citations"][
                      "why_it_is_a_citation_and_not_this_run_s_null"],
                  "citations": band["null_citations"][
                      "shuffled_contigs_with_a_significant_hit"]}),
        "what_this_refuses": [
            "Ranking by identity: measured wrong on the one positive, which sits at the "
            "93.46th percentile of removal identity.",
            "A per-contig flag: the candidate is not the top row of its own cell.",
            "A nucleotide-level band from this file: only the protein route carries an "
            "identity and a coverage on every row. The nucleotide cell is 137 contigs, "
            "6.64%, measured in results/holdout_cascade/ and not reproducible from here.",
        ],
    }
    if args.shuffled_evidence:
        sbest, scensus = best_protein_rows(read_evidence(args.shuffled_evidence), evalue)
        out["null"].update({
            "contigs_scored": len(sbest), "census": scensus,
            "contigs_inside_the_band": len(in_band(sbest, cov, ident))})

    (args.out / "removal_shortlist.json").write_text(json.dumps(out, indent=2))
    cols = ["contig", "contig_nt", "explains_fraction", "identity_percent", "aln_aa",
            "evalue", "bitscore", "index", "subject", "tier"]
    with open(args.out / "removal_shortlist.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in shortlist:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    (args.out / "removal_shortlist_manifest.json").write_text(json.dumps({
        "script": "experiments/removal_shortlist.py",
        **provenance(__file__),
        "seed": args.seed,
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "consulted": [{"role": "removal evidence", "path": str(args.evidence),
                       "state": "searched", "queries": census["contigs"],
                       "hits": len(best)}] +
                     ([{"role": "null", "path": str(args.shuffled_evidence),
                        "state": "searched", "queries": out["null"]["census"]["contigs"],
                        "hits": out["null"]["contigs_scored"]}]
                      if args.shuffled_evidence else
                      [{"role": "null", "path": None, "state": "not searched",
                        "why": "dropped by name with --without-null"}]),
        "wall_seconds": round(time.perf_counter() - t0, 2),
    }, indent=2))
    print(f"{len(best):,} removals scored; band {cov[0]}-{cov[1]} coverage x "
          f"{ident[0]}-{ident[1]}% identity holds {len(shortlist):,} contigs "
          f"({out['band']['percent_of_scored_pile']}% of the scored pile); "
          f"control arm {len(control):,} near-reference removals, "
          f"{len(control_in_band)} of them inside the band")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
