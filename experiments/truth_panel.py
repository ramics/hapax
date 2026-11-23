#!/usr/bin/env python3
"""Check a run's assembly against the depositors' own genomes from the same water.

WHAT THIS PROMOTES. `results/path_settings/` scored ten assembly arms against a ground truth
nobody chose: the depositors assembled 25 virus genomes out of the same library and deposited
them a year later, and `results/holdout_candidates/depth.tsv` records each one's depth in the
run by mapping the reads back - 0.4x to 700x, at 0.08% to 30.4% divergence from public
reference. 53 further genomes in the same committed panel have exactly zero reads mapped in
that run, which makes a false-join control free and drawn from the same pool. That was run by
hand, once, over an assembly ladder. This script makes it a check any run gets.

WHY A RUN NEEDS IT RATHER THAN A SWEEP. The number this project reports about its own
pipeline - how much of the library is residue - moves 8.79 points on assembler flags that do
not improve recovery of a single real genome. A residue fraction quoted without a recovery
figure from the same pass is therefore unreadable: a setting that shrinks the residue by
fusing novel sequence into an abundant known contig reads as an improvement on every summary
number in this repo. So recovery and residue size ship as a PAIR, and this stage is the half
that cannot be computed from the residue.

THREE THINGS IT IMPORTS RATHER THAN RETYPES, each of which was a defect once:

  - `path_settings.paf` aligns with SECONDARY ALIGNMENTS ON. With `--secondary=no` each
    contig is credited to one genome only, so one member of every near-identical pair reads
    0.00 breadth: recovery read 9/25 instead of 20/25 and a 499x genome looked like an
    assembly failure.
  - `path_settings.lineage_groups` groups accessions that are the same thing to within
    sequencing error, and purity is measured against the GROUP. OQ722323 is a 5,048 nt
    partial deposit 99.84% identical to the 9,467 nt OQ722327; per-accession purity reports a
    correct recovery of the longer relative as absorption into foreign sequence.
  - `path_settings.build_panel` derives present and absent from the committed depth table
    rather than from a hand-written list, so the control arm cannot drift from the panel.

It imports them; it does not copy them and does not edit them. Editing `path_settings.py`
would change the `script_sha256` under a committed five-seed directory nobody is re-running.

WHAT THE PANEL CANNOT TEST, AND WHY THAT TRAVELS IN THE OUTPUT. The panel is 25 separate
molecules, not a lineage mixture, so it CANNOT test the sub-3.8% fusion that
`results/pair_ladder/` measures - and it cannot be made to: it does contain at least six
pairs above 99.4% identity, but at 99.9% identity a 200 nt alignment at >=95% cannot say
which member a contig came from, so breadth is the wrong instrument for that question
whatever the panel holds. Nor does anything in it sit at the median depth 1 where most of the
residue lives. Those sentences are resolved from the committed reference file and written
into the stage's own output, because a recovery figure read without them is read as a
completeness claim about the whole library.
"""

from __future__ import annotations

import argparse
import json
import statistics as st
import sys
import time
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import stage_refs
from loud import Preflight
from path_settings import build_panel, lineage_groups, read_fasta, score
from provenance import provenance

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    t0 = time.perf_counter()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded; the panel and the assembly are both fixed inputs, so "
                         "this stage is deterministic and the seed identifies the run")
    ap.add_argument("--library", required=True,
                    help="the run accession, as experiments/stage_references.json spells "
                         "it. REQUIRED and never inferred: the present/absent split is per "
                         "library, and scoring one library's assembly against another's "
                         "panel would report a recovery for genomes that are not there")
    ap.add_argument("--contigs", type=Path, default=None,
                    help="this run's assembly. Without it the stage composes and describes "
                         "the panel and scores nothing, which is the free half")
    ap.add_argument("--residue-fraction", type=float, default=None,
                    help="the residue fraction from the SAME pass, so recovery and residue "
                         "size travel together. Recorded as absent rather than omitted")
    ap.add_argument("--minimap2", type=Path, default=None,
                    help="required when --contigs is given")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-identity", type=float, default=None)
    ap.add_argument("--min-aln", type=int, default=None)
    ap.add_argument("--min-mosaic", type=int, default=None)
    ap.add_argument("--references", type=Path, default=None)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--preflight", action="store_true")
    args = ap.parse_args()

    pre = Preflight("truth_panel.py")
    resolved = None
    try:
        resolved = stage_refs.resolve(args.library, "truth_panel",
                                      stage_refs.load(args.references))
    except (stage_refs.UnknownLibrary, stage_refs.MissingReference) as e:
        # An unknown library is a refusal. A KNOWN library with no panel is not: it is
        # handled below as a measured absence with its committed reason.
        pre.require(False, str(e))
    if resolved is not None and resolved["state"] == "resolved":
        for problem in stage_refs.missing_files(resolved):
            pre.require(False, f"truth panel {resolved['panel']}: {problem}")
        if args.contigs is not None:
            pre.file(args.contigs, "assembly")
            pre.exe(args.minimap2, "minimap2")
    pre.done(args.preflight)

    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    out_path = args.out / "truth_panel.json"
    manifest_path = args.out / "truth_panel_manifest.json"

    def write(payload: dict, consulted: list[dict]) -> None:
        payload.setdefault("library", args.library)
        payload.setdefault("spec", "SPEC.md 4.3")
        out_path.write_text(json.dumps(payload, indent=2))
        manifest_path.write_text(json.dumps({
            "script": "experiments/truth_panel.py", **provenance(__file__),
            "seed": args.seed,
            "parameters": {k: (str(v) if isinstance(v, Path) else v)
                           for k, v in vars(args).items()},
            "consulted": consulted,
            "wall_seconds": round(time.perf_counter() - t0, 2)}, indent=2))

    if resolved["state"] == "absent":
        # Berlin's case, and it is a result rather than a gap: no depositor genome set
        # exists for that water, so this run's positive control is the spiked ladder. The
        # file says which, so a reader cannot mistake a skipped panel for a passed one.
        write({"status": "skipped", "reason": resolved["reason"],
               "not_the_same_as": ("a panel that was run and found nothing. Nothing was "
                                   "scored here and no recovery figure exists for this "
                                   "library")},
              [{"role": "truth panel", "path": None, "state": "not searched",
                "why": resolved["reason"]}])
        print(f"truth_panel: skipped for {args.library}: {resolved['reason']}")
        return 0

    panel_record = resolved["record"]
    op = panel_record["operating_point"]
    files = resolved["files"]
    ns = types.SimpleNamespace(
        run=args.library,
        depth=files["depth_table"], candidates=files["genomes"],
        shortlist=files["divergence_table"], work=args.work,
        minimap2=args.minimap2, threads=args.threads,
        min_identity=args.min_identity if args.min_identity is not None
        else op["min_identity"],
        min_aln=args.min_aln if args.min_aln is not None else op["min_aln"],
        min_mosaic=args.min_mosaic if args.min_mosaic is not None else op["min_mosaic"])

    panel_fa, meta = build_panel(ns)
    present, absent = meta["present"], meta["absent"]
    if not present:
        raise SystemExit(
            f"{args.library} resolves to truth panel {resolved['panel']} but the depth "
            f"table records no genome with a mapped read for it. An empty present arm is "
            f"not a recovery of zero: it means this library is not the one the panel was "
            f"measured in.")
    depths = sorted(v["mean_depth"] for v in present.values())
    divs = [d for d in (meta["divergence_pct"].get(a) for a in present) if d is not None]
    payload = {
        "status": "ok",
        "panel": resolved["panel"],
        "what_this_is": panel_record["_what"],
        "operating_point": {**op, "min_identity": ns.min_identity,
                            "min_aln": ns.min_aln, "min_mosaic": ns.min_mosaic},
        "composition": {
            "present": len(present), "absent_control": len(absent),
            "present_is": ("mapped_reads > 0 in the committed depth table for this run, "
                           "measured by mapping the reads back rather than asserted"),
            "absent_is": ("exactly 0 mapped reads in the same table, same panel, same "
                          "code, so the false-join control costs nothing and cannot be "
                          "accused of coming from a different distribution of sequence"),
            "mean_depth_range": [depths[0], depths[-1]],
            "mean_depth_median": st.median(depths),
            "divergence_percent_range": [min(divs), max(divs)] if divs else None,
            "panel_bases": sum(len(s) for s in read_fasta(panel_fa).values()),
        },
        "what_this_panel_cannot_test": panel_record["what_this_panel_cannot_test"],
        "purity_is_against_the_group": panel_record["purity_is_against_the_group"],
        "residue_fraction_from_the_same_pass": (
            args.residue_fraction if args.residue_fraction is not None else
            "not measured in this pass: a recovery figure here does not license a residue "
            "claim, and residue size moves 8.79 points on assembler flags that move "
            "recovery by at most one genome (results/path_settings/)"),
    }
    consulted = [{"role": "truth panel genomes", "path": str(files["genomes"]),
                  "state": "searched", "queries": len(present) + len(absent),
                  "hits": len(present)},
                 {"role": "depth table", "path": str(files["depth_table"]),
                  "state": "searched", "queries": 1, "hits": 1}]

    if args.contigs is None:
        payload["scoring"] = {
            "state": "not scored",
            "why": "--contigs was not given, so the panel was composed and described only"}
        write(payload, consulted)
        print(f"truth_panel: {args.library} panel composed: {len(present)} present, "
              f"{len(absent)} absent controls; nothing scored (--contigs absent)")
        return 0

    groups = lineage_groups(ns, panel_fa)
    scored = score(ns, "run", args.contigs, panel_fa, meta, groups)
    rows = [r for r in scored["genomes"].values() if r["present"]]
    purities = [r["best_contig_purity"] for r in rows if r["best_contig_purity"] is not None]
    payload["scoring"] = {
        "state": "scored",
        "contigs_file": str(args.contigs),
        "near_identical_groups": scored["panel_present"] and len(set(groups.values())),
        "recovered_90pct": scored["recovered_90pct"],
        "recovered_50pct": scored["recovered_50pct"],
        "of_present": len(rows),
        "mean_covered_fraction": round(scored["mean_covered_frac"], 4),
        "median_best_contig_purity": scored["median_purity"],
        "purity_floor_note": ("a genome recovered at >=90% breadth sitting below 0.9 purity "
                              "is the absorption failure mode; zero arms showed one in "
                              "results/path_settings/"),
        "absorbed_recoveries": sum(1 for r in rows if r["covered_frac"] >= 0.9
                                   and (r["best_contig_purity"] or 1.0) < 0.9),
        "inter_group_mosaics": scored["inter_group_mosaics"],
        "mosaic_detail": scored["mosaic_detail"],
        "false_joins_on_the_absent_arm": scored["false_joins_absent_arm"],
        "false_join_is": (f"any contig covering at least {ns.min_aln} nt of a genome with "
                          f"zero mapped reads in this run"),
        "purity_median_excludes": len(rows) - len(purities),
    }
    payload["per_genome"] = scored["genomes"]
    consulted.append({"role": "assembly", "path": str(args.contigs), "state": "searched",
                      "queries": len(read_fasta(args.contigs)),
                      "hits": sum(1 for r in rows if r["contigs_touching"])})
    write(payload, consulted)
    print(f"truth_panel: {args.library}: {scored['recovered_90pct']}/{len(rows)} genomes at "
          f">=90% breadth, {scored['inter_group_mosaics']} inter-group mosaics, "
          f"{len(scored['false_joins_absent_arm'])} false joins on {len(absent)} zero-read "
          f"controls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
