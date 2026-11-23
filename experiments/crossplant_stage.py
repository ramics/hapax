#!/usr/bin/env python3
"""Subtract other plants' water from this run, with the panel's rules carried as checks.

WHAT THIS PROMOTES. `results/crossplant/` measured, over 530 other wastewater libraries from
32 plants on four continents, what a union of other plants' backgrounds removes from a
held-out library and from this sample's own piles. Both arms are merged. The measurement
exists; what did not exist is a way for a NEW sample to get it without someone re-deriving
the rules, and the rules are the part that is easy to lose:

  1. **Query reads must come from OUTSIDE the library searched.** Own-library reads are
     abundant by construction - measured at 54% false growth against a 6% swap control - so
     the query library is never a panel member and the queries are drawn from spots outside
     each library's own subject subset. This script REFUSES to run a library that is in the
     panel it was asked to subtract against.
  2. **Report paired per plant, never pooled.** Within-plant rates run 0.1530 to 0.9593, a
     six-fold spread, so a pooled within-minus-between difference compares plants instead of
     comparing within against between. The paired gap is +0.1432 (`most1`), positive in 32
     of 32 plants. This script refuses to emit the pooled reading and fails if the paired
     strata are missing from the report it reads.
  3. **Contig queries go in 90 nt windows, the same width as the read queries**, so the
     instrument does not change between the two arms. A width that disagrees with the
     panel's own query length is a refusal, not a warning.
  4. **Both operating points travel.** `any` (one exact 31-mer present) and `most` (half or
     more) disagree by a factor of four about the thing the measurement is for: at depth 1
     the identified-versus-residue gap is 8.5 points at `most` and 42.0 at `any`. Neither
     ships alone, so the output carries both or fails.
  5. **Depth decides and the pile label almost nothing.** Every saturation figure is read
     down the depth column, not across the pile row: at depth 100+ the union takes
     everything including our own abundant residue, and at depth 1 identified reads survive
     at 0.7050 against the residue's 0.7900. A union is blind to anything present once, not
     to novelty, and that sentence travels in the output beside the number.

WHAT IT ORCHESTRATES RATHER THAN REIMPLEMENTS. `crossplant_contig_queries.py` cuts the
windows, `crossplant_report.py` and `crossplant_contig_report.py` produce the per-plant
matrices, the saturation curves and the controls - including the paired strata, which were
already computed there. This script resolves the panel, enforces the five rules above,
drives those scripts and writes one record saying which panel was used, whether it was
reused from its staged copy or re-staged, and what the panel cannot say.

THE PANEL IS STAGED FOR REUSE, AND WHAT IS NOT STAGED IS STATED. The reduced query x subject
matrices of the passes already made live in S3 and are what a report reads. The 531 subject
FASTQ subsets are NOT staged: each pass re-downloads them on the instance, which is the bulk
of its wall time. That is a cost, not a failure, and this stage records which of the two it
did rather than letting a silent re-download look like reuse.

WHAT THE PANEL CANNOT SAY, carried from the committed reference file: "absent" here means no
31 exact bases shared with any panel library, and past roughly 10-15% divergence that is not
the same as having no relative there - one substitution destroys up to 31 of a 90 nt mate's
60 31-mers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import stage_refs
from loud import Preflight, run_tool
from provenance import provenance

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent

#: The four per-mate operating points the reports write. Both sensitivities for both mates,
#: because the mates are not interchangeable and the two points disagree by a factor of four.
REQUIRED_OPS = ("most1", "most2", "any1", "any2")
#: The strata that make the per-plant reading paired. A report without these is a report
#: that can only be read pooled, and the pooled reading is the one this measurement refuses.
REQUIRED_STRATA = ("within_plant_other_water", "other_plant_other_project",
                   "paired_within_minus_other_project", "paired_plants_positive")


def panel_members(panel_tsv: Path) -> dict[str, dict]:
    with open(panel_tsv) as fh:
        return {r["run"]: r for r in csv.DictReader(fh, delimiter="\t")}


def own_library_refusal(library: str, members: dict[str, dict]) -> str | None:
    """RULE 1, as a value: why this library may not be subtracted against this panel.

    A library inside the panel it is being subtracted against is searched partly against
    itself, and own-library reads are abundant by construction - 54% false growth against
    a 6% swap control - so the rate it returns is a property of the draw and not of the
    water. Returned as a string rather than raised so the stage can list it beside its
    other preflight problems, and so a test can assert the rule without a fixture library.
    """
    if library not in members:
        return None
    return (f"{library} is a member of this panel, so its query reads would be searched "
            f"against its own library. Query reads must come from outside the library "
            f"searched; own-library reads are abundant by construction (54% false growth "
            f"against a 6% swap control). Hold it out of the panel or pick another panel")


def eligible_subjects(members: dict[str, dict], depth: int,
                      headroom: int = 200_000) -> list[str]:
    """Which panel rows a pass can actually scan at the fixed subject depth.

    NOT THE SAME NUMBER AS THE ROWS IN THE TABLE, and the difference is the point: the
    committed panel holds 536 rows and the pass scans 530 of them, because a library below
    the fixed depth is LEFT OUT rather than scanned shallower - a subject at a lower depth
    would make a between-plant difference partly a difference in library size, which is
    what the fixed depth exists to remove. The headroom is the query pool drawn from spots
    outside the subject, which is rule 1 in sequence form.
    """
    return sorted(run for run, r in members.items()
                  if r.get("paired") == "1"
                  and int(r.get("read_count") or 0) >= depth + headroom)


def sha16(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def missing_paired_rows(rows: list[dict]) -> list[str]:
    """Which required (operating point, stratum) rows a per-plant report does not carry.

    Separate from the stage so it can be asserted on a table directly: the committed
    `results/crossplant/plant_summary.tsv` must come back with nothing missing, which is
    the check that this stage's idea of a readable report is the one the published report
    satisfies rather than a stricter one invented here.
    """
    have = {(r.get("op"), r.get("stratum")) for r in rows}
    return [f"{op}/{s}" for op in REQUIRED_OPS for s in REQUIRED_STRATA
            if (op, s) not in have]


def main() -> int:
    t0 = time.perf_counter()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded; the reports read seeds 1-5 and quote the spread, which "
                         "is the only honest dispersion this repo reports")
    ap.add_argument("--library", required=True,
                    help="the query run accession, as experiments/stage_references.json "
                         "spells it. REQUIRED: the panel a library may be subtracted "
                         "against is a property of the library, and a default would "
                         "subtract one sample's plants from another's water")
    ap.add_argument("--stages", nargs="+", default=["queries", "report"],
                    choices=["queries", "report"],
                    help="`queries` cuts this run's contigs into panel-width windows "
                         "(free, in-container); `report` reads a completed pass's run "
                         "directory and writes the per-plant and saturation tables")
    ap.add_argument("--contig-dir", type=Path, default=None,
                    help="directory holding sample_contigs_unidentified.fa and "
                         "sample_contigs_identified.fa: tier 2 and its known-background "
                         "comparator from the same assembly. Required for `queries`")
    ap.add_argument("--run-dir", type=Path, default=None,
                    help="a completed cross-plant pass downloaded from S3 (packed/, "
                         "sample/). Required for `report`")
    ap.add_argument("--arm", default="reads", choices=["reads", "contigs"],
                    help="which arm's report to write from --run-dir. The two arms use the "
                         "same panel at the same depth and the same 90 nt query width")
    ap.add_argument("--window", type=int, default=None,
                    help="contig window width; defaults to the panel's own query length so "
                         "the instrument does not change between arms. A disagreement is a "
                         "refusal")
    ap.add_argument("--references", type=Path, default=None)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--preflight", action="store_true")
    args = ap.parse_args()

    pre = Preflight("crossplant_stage.py")
    resolved = None
    try:
        resolved = stage_refs.resolve(args.library, "crossplant_panel",
                                      stage_refs.load(args.references))
    except (stage_refs.UnknownLibrary, stage_refs.MissingReference) as e:
        pre.require(False, str(e))
    members: dict[str, dict] = {}
    if resolved is not None and resolved["state"] == "resolved":
        for problem in stage_refs.missing_files(resolved):
            pre.require(False, f"cross-plant panel {resolved['panel']}: {problem}")
        panel_tsv = resolved["files"].get("panel_table")
        if panel_tsv and panel_tsv.exists():
            members = panel_members(panel_tsv)
            own = own_library_refusal(args.library, members)
            pre.require(own is None, f"panel {resolved['panel']}: {own}")
        want = args.window or resolved["record"]["query_length_nt"]
        # RULE 3, as a refusal: one width across both arms or the arms are two instruments.
        pre.require(want == resolved["record"]["contig_window_nt"] ==
                    resolved["record"]["query_length_nt"],
                    f"window {want} nt disagrees with the panel's query length "
                    f"{resolved['record']['query_length_nt']} nt. The contig arm and the "
                    f"read arm must offer the same number of 31-mers per query or they are "
                    f"two instruments and their rates are not comparable")
        if "queries" in args.stages:
            pre.require(args.contig_dir is not None, "`queries` needs --contig-dir")
            if args.contig_dir is not None:
                for name in ("sample_contigs_unidentified.fa",
                             "sample_contigs_identified.fa"):
                    pre.file(args.contig_dir / name, f"contig set {name}")
        if "report" in args.stages:
            pre.require(args.run_dir is not None, "`report` needs --run-dir")
            if args.run_dir is not None:
                pre.require((args.run_dir / "packed").is_dir(),
                            f"{args.run_dir}/packed is not a directory: a cross-plant "
                            f"report reads the reduced matrices, and an absent packed/ is "
                            f"a pass that did not finish rather than a pass with no hits")
    elif resolved is not None and resolved["state"] == "absent":
        pass  # handled below as a measured absence
    pre.done(args.preflight)

    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    out_path = args.out / "crossplant_stage.json"
    consulted: list[dict] = []

    def write(payload: dict) -> None:
        payload.setdefault("library", args.library)
        payload.setdefault("spec", "SPEC.md 7.1")
        out_path.write_text(json.dumps(payload, indent=2))
        (args.out / "crossplant_stage_manifest.json").write_text(json.dumps({
            "script": "experiments/crossplant_stage.py", **provenance(__file__),
            "seed": args.seed,
            "parameters": {k: (str(v) if isinstance(v, Path) else
                               [str(x) for x in v] if isinstance(v, list) else v)
                           for k, v in vars(args).items()},
            "consulted": consulted,
            "wall_seconds": round(time.perf_counter() - t0, 2)}, indent=2))

    if resolved["state"] == "absent":
        write({"status": "skipped", "reason": resolved["reason"],
               "not_the_same_as": ("a panel that was searched and removed nothing. Nothing "
                                   "was subtracted and no cross-plant figure exists for "
                                   "this library")})
        print(f"crossplant: skipped for {args.library}: {resolved['reason']}")
        return 0

    record = resolved["record"]
    panel_tsv = resolved["files"]["panel_table"]
    window = args.window or record["query_length_nt"]
    plants = sorted({f"{r['project']}:{r['plant'] or 'untreated'}"
                     for r in members.values()})
    eligible = eligible_subjects(members, record["subject_depth_pairs"])
    payload = {
        "status": "ok",
        "panel": {
            "name": resolved["panel"],
            "table": str(panel_tsv), "table_sha256_16": sha16(panel_tsv),
            "rows": len(members),
            "eligible_subjects_at_depth": len(eligible),
            "why_fewer_than_rows": ("a library below the fixed subject depth is left out "
                                    "rather than scanned shallower; the committed pass "
                                    "scanned 530 of this table's 536 rows"),
            "plants": len(plants),
            "projects": sorted({r["project"] for r in members.values()}),
            "subject_depth_pairs": record["subject_depth_pairs"],
            "query_length_nt": record["query_length_nt"],
            "contig_window_nt": window,
            "witness_instrument": record["witness_instrument"],
            "staged_for_reuse": record["staged_for_reuse"],
            "query_library_is_not_a_member": args.library not in members,
        },
        "rules_enforced": record["rules_that_travel_with_this_panel"],
        "operating_points": {"required": list(REQUIRED_OPS),
                             "both_travel": record["both_points_travel"]},
        "what_this_panel_cannot_say": record["what_this_panel_cannot_say"],
        "stages": {},
    }
    consulted.append({"role": "cross-plant panel", "path": str(panel_tsv),
                      "state": "searched", "queries": len(members),
                      "hits": len(eligible)})

    if "queries" in args.stages:
        qdir = args.work / "queries"
        qdir.mkdir(parents=True, exist_ok=True)
        # --out is the query FILE, not a directory: the script writes its manifest beside
        # it. Named explicitly here so the stage's own output can point at both.
        qfile = qdir / "contig_queries.tsv"
        run_tool([sys.executable, "-I", str(HERE / "crossplant_contig_queries.py"),
                  "--contig-dir", str(args.contig_dir), "--out", str(qfile),
                  "--window", str(window)],
                 "crossplant_contig_queries.py", cwd=ROOT, stdout=subprocess.DEVNULL)
        if not qfile.exists() or qfile.stat().st_size == 0:
            raise SystemExit(f"the contig query set {qfile} was not written or is empty")
        rows = sum(1 for _ in open(qfile))
        payload["stages"]["queries"] = {
            "state": "built", "file": str(qfile), "rows": rows, "window_nt": window,
            "manifest": str(qdir / "contig_queries_manifest.json"),
            "why_windows": ("a single 90 nt read at depth one is a weak probe of whether a "
                            "molecule is background elsewhere; a contig offers many "
                            "independent windows of the same molecule, and the two windows "
                            "of a query pair are a replicate rather than a second "
                            "instrument because they are the same length and the same kit"),
            "next_step": ("upload beside the panel under the pass's query-set prefix and "
                          "launch with experiments/ec2_crossplant.py --query-mode contigs; "
                          "subject subsets are not staged, so that pass pays a fresh "
                          "download"),
        }
        consulted.append({"role": "this run's contigs", "path": str(args.contig_dir),
                          "state": "searched", "queries": rows, "hits": rows})

    if "report" in args.stages:
        script = ("crossplant_report.py" if args.arm == "reads"
                  else "crossplant_contig_report.py")
        rdir = args.out / f"{args.arm}_arm"
        run_tool([sys.executable, "-I", str(HERE / script),
                  "--run-dir", str(args.run_dir), "--panel", str(panel_tsv),
                  "--out", str(rdir)], script, cwd=ROOT)
        stage = {"state": "reported", "script": f"experiments/{script}",
                 "run_dir": str(args.run_dir), "out": str(rdir)}
        if args.arm == "reads":
            summary = rdir / "plant_summary.tsv"
            controls = rdir / "controls.json"
            if not summary.exists():
                raise SystemExit(f"{script} wrote no plant_summary.tsv, so there is no "
                                 f"per-plant reading to check")
            rows = read_tsv(summary)
            # RULE 2 and RULE 4, as refusals after the fact: a report missing the paired
            # strata or an operating point can only be read in the shape this measurement
            # refuses, so the stage fails rather than passing it on.
            missing = missing_paired_rows(rows)
            if missing:
                raise SystemExit(
                    f"the per-plant report is missing {len(missing)} required rows "
                    f"({missing[:6]}...). Paired per-plant strata at both operating points "
                    f"for both mates are what makes this readable: pooled, the comparison "
                    f"is between plants rather than within against between, and one "
                    f"operating point alone hides a factor of four.")
            paired = {}
            for r in rows:
                if r["stratum"] == "paired_within_minus_other_project":
                    paired[r["op"]] = {
                        "plants": int(r["plants"]), "median": r["witnessed_median"],
                        "plant_min": r["plant_min"], "plant_max": r["plant_max"]}
                if r["stratum"] == "paired_plants_positive":
                    paired.setdefault(r["op"], {})["plants_positive"] = \
                        int(r["witnessed_median"])
            stage["paired_per_plant"] = paired
            stage["pooled_reading"] = (
                "not emitted. Within-plant rates run 0.1530 to 0.9593, so a pooled "
                "within-minus-between difference compares plants rather than comparing "
                "within against between")
            if controls.exists():
                c = json.loads(controls.read_text())
                depth = c.get("subject_pairs_min_max")
                # The fixed depth is what makes a difference between plants a difference
                # between waters rather than between library sizes. An unequal depth is
                # recorded as a confound below rather than quietly accepted.
                stage["controls"] = {
                    "subject_pairs_min_max": depth,
                    "subject_depth_is_fixed": bool(depth and depth[0] == depth[1]),
                    "why_fixed_depth": ("subjects are all cut to the same depth, so a "
                                        "difference between plants is not a difference in "
                                        "library size"),
                    "shuffled_floor": c.get("shuffled_in_other_libraries"),
                    "own_library_ceiling": c.get(
                        "own_library_other_spots_min_over_libraries"),
                }
                if depth and depth[0] != depth[1]:
                    stage["controls"]["depth_warning"] = (
                        f"subject depths range {depth}, so a between-plant difference in "
                        f"this run is confounded with library size and the plant matrix "
                        f"must not be read as a plant difference")
            stage["read_the_depth_column"] = (
                "Depth decides and the pile label almost nothing: at depth 100+ the union "
                "takes everything including our own abundant residue, and at depth 1 "
                "identified reads survive at 0.7050 against the residue's 0.7900. A union "
                "is blind to anything present once, not to novelty.")
            consulted.append({"role": "cross-plant pass (reads arm)",
                              "path": str(args.run_dir), "state": "searched",
                              "queries": len(rows), "hits": len(rows)})
        else:
            stage["tables"] = sorted(p.name for p in rdir.glob("contig_*"))
            stage["unit"] = ("a contig, not a query row: a contig's windows are aggregated "
                             "into one witnessed-window fraction per subject, with "
                             "identified contigs from the same assembly as the comparator")
            consulted.append({"role": "cross-plant pass (contig arm)",
                              "path": str(args.run_dir), "state": "searched",
                              "queries": 1, "hits": 1})
        payload["stages"][args.arm] = stage

    for name in ("queries", "report"):
        if name not in args.stages:
            payload["stages"].setdefault(name, {
                "state": "not selected",
                "why": f"{name} was not in --stages, so nothing of it ran"})
    write(payload)
    print(f"crossplant: panel {resolved['panel']}, {len(eligible)} of {len(members)} "
          f"rows scannable at depth {record['subject_depth_pairs']:,} over "
          f"{len(plants)} plants; " +
          ", ".join(f"{k}={v.get('state')}" for k, v in payload["stages"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
