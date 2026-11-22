"""Which committed results could the silent-failure bugs fixed with `loud.py` have touched?

WHY. The fail-loudly pass fixed six latent bugs of one family - a stage doing less than
asked and reporting success. Fixing code does not say whether a committed number already
rests on one of them, so this reads every committed run manifest and result file and
checks each bug's precondition against what the files record. It reads results only; it
re-runs nothing.

LIMIT, stated: these checks see what manifests recorded. A reference regenerated IN PLACE
(same path, new content) between two runs sharing a work directory is invisible to them,
which is exactly why the cache is now keyed to size and mtime.

Usage: python3 experiments/fail_loudly_audit.py --seed 1   (the seed is unused; the
convention requires it)
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def manifests() -> list[tuple[str, dict]]:
    out = []
    for f in sorted(glob.glob(str(ROOT / "results/**/*.json"), recursive=True)):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        if isinstance(d, dict) and "script" in d and isinstance(d.get("parameters"), dict):
            out.append((os.path.relpath(f, ROOT), d))
    return out


def stem(p: str) -> str:
    n = Path(p).name
    for suf in (".fa.gz", ".fasta.gz", ".fa", ".fasta", ".dmnd"):
        if n.endswith(suf):
            return n[: -len(suf)]
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=ROOT / "results/fail_loudly_audit")
    args = ap.parse_args()
    ms = manifests()
    by = defaultdict(list)
    for f, d in ms:
        by[Path(d["script"]).name].append((f, d))
    verdict = {}

    # 1. hapax.py runs (relative paths; the two index-passing bugs)
    hx = [(f, d) for f, d in ms if d.get("tool") == "hapax run" or d["script"] == "hapax.py"]
    rel = [f for f, d in hx
           if any(isinstance(v, str) and v.endswith((".gz", ".fa", ".fasta", ".dmnd"))
                  and not v.startswith("/") for v in d["parameters"].values())]
    verdict["relative_paths_in_hapax"] = {
        "hapax_runs": [f for f, _ in hx], "affected": rel,
        "checked": "every committed hapax run manifest, for a relative input path",
        "verdict": ("no committed output came from hapax.py" if not hx else
                    "every committed hapax run was given absolute paths" if not rel
                    else "see affected")}

    # 2. cascade: a requested tier missing from the curve, or Kraken2 operating
    # confidence outside the sweep (which made every Kraken2 tier remove nothing)
    casc = []
    for f, d in by["subtract_known.py"]:
        res = (Path(ROOT / f).parents[1] / "tiers.json"
               if f.endswith("stages/cascade_manifest.json")
               else Path(ROOT / f.replace("manifest_s", "subtract_s")))
        if not res.exists():
            continue
        r = json.loads(res.read_text())
        p = d["parameters"]
        tiers = {c["tier"]: c["tier_alone_removed"] for c in r["curve"]}
        want = {Path(x).name for x in p["k2_dbs"]}
        want |= {"selfasm_" + Path(x).stem for x in p.get("selfasm", [])}
        want |= {"prot_" + Path(x).stem for x in p.get("dmnd", [])}
        if p.get("silva"):
            want.add("rrna")
        sweep = [float(x) for x in p["k2_confidence"]]
        casc.append({
            "manifest": f, "requested_tiers_missing_from_curve": sorted(want - set(tiers)),
            "operating_confidence_in_sweep": float(p["k2_operating_confidence"]) in sweep,
            "kraken2_tiers_removing_nothing": sorted(
                t for t in (Path(x).name for x in p["k2_dbs"]) if tiers.get(t) == 0),
            "has_script_sha256": "script_sha256" in d,
            "git_sha": d.get("git_sha")})
    bad = [c for c in casc if c["requested_tiers_missing_from_curve"]
           or not c["operating_confidence_in_sweep"] or c["kraken2_tiers_removing_nothing"]]
    verdict["cascade_silent_tier_skip_and_k2_operating_point"] = {
        "runs": casc, "affected": [c["manifest"] for c in bad],
        "verdict": ("every committed cascade run has every requested tier in its curve, its "
                    "operating confidence in its sweep, and every Kraken2 tier removing "
                    "something" if not bad else "see affected")}
    shas = defaultdict(set)
    for c in casc:
        shas[str(Path(c["manifest"]).parent)].add(c["git_sha"])
    verdict["cascade_manifest_without_script_hash"] = {
        "affected": [c["manifest"] for c in casc if not c["has_script_sha256"]],
        "git_shas_per_directory": {k: sorted(v) for k, v in shas.items()},
        "verdict": "numbers unaffected; provenance gap. A directory whose seeds carry more "
                   "than one SHA, or a -dirty SHA, cannot be pinned to code without a "
                   "content hash"}

    # 3. census: a class with no reference labels its neighbours `none`; --silva unread
    cen = [{"manifest": f, "references": sorted(d.get("references", {})),
            "silva_in_parameters_but_unread": bool(d["parameters"].get("silva"))}
           for f, d in by["linkage_census.py"]]
    need = {"rrna", "identified_contig", "unidentified_contig"}
    verdict["census_class_without_reference"] = {
        "runs": cen, "affected": [c["manifest"] for c in cen
                                  if set(c["references"]) != need],
        "verdict": ("every committed census run labelled all three classes"
                    if all(set(c["references"]) == need for c in cen) else "see affected")}
    verdict["census_silva_accepted_never_read"] = {
        "affected": [c["manifest"] for c in cen if c["silva_in_parameters_but_unread"]],
        "verdict": "numbers unaffected (the census never read SILVA); these manifests name "
                   "an index in `parameters` that was not consulted"}

    # 4. per-read report: an optional field absent, and motif_rate 0.0 when not run
    rr = []
    for f, d in by["residue_report.py"]:
        cal = ROOT / f.replace("manifest_s", "calibration_s")
        if not cal.exists():
            continue
        c = json.loads(cal.read_text())
        p = d["parameters"]
        pfam_ran = not p.get("skip_fields") and Path(p["pfam"]).name in {
            Path(k).name for k in d.get("index_bytes", {})}
        rr.append({"manifest": f, "skip_fields": p.get("skip_fields"),
                   "motif_rate": c.get("motif_rate"),
                   "motif_rate_is_unmeasured_zero": (c.get("motif_rate") == 0.0
                                                     and p.get("skip_fields")),
                   "closest_fields": sorted(k for k in c["fields"]
                                            if k.startswith("closest_")),
                   "dmnd_requested": sorted("closest_" + stem(x) for x in p["dmnd"]),
                   "pfam_index_listed_in_index_bytes": not pfam_ran and Path(
                       p["pfam"]).name in {Path(k).name for k in d.get("index_bytes", {})}})
    # The same unmeasured zeros inside hapax run directories, whose fields.json is the
    # per-read report's calibration file copied in.
    for fj in sorted(glob.glob(str(ROOT / "results/**/fields.json"), recursive=True)):
        c = json.loads(Path(fj).read_text())
        rr.append({"manifest": os.path.relpath(fj, ROOT), "skip_fields": None,
                   "motif_rate": c.get("motif_rate"),
                   "motif_rate_is_unmeasured_zero": False,
                   "closest_fields": sorted(k for k in c.get("fields", {})
                                            if k.startswith("closest_")),
                   "dmnd_requested": [], "pfam_index_listed_in_index_bytes": False})
    for r in rr:
        cal = ROOT / (r["manifest"] if r["manifest"].endswith("fields.json") else
                      r["manifest"].replace("manifest_s", "calibration_s"))
        c = json.loads(cal.read_text())
        # A chimaera rate of exactly 0.0 beside a control that was never computed means
        # the test did not run: whenever it runs on a run with a known pile, the control
        # is computed too.
        r["chimaera_rate_is_unmeasured_zero"] = (
            c.get("chimaera_split_rate") == 0.0
            and c.get("chimaera_split_rate_on_subtracted_reads") is None)
        r["fields_silently_absent"] = sorted(
            k for k, key in (("cross_sample", "cross_sample_recurrence"), ("rfam", "rfam"),
                             ("crispr", "crispr"), ("abundance", "library_abundance"))
            if c.get(key) is None and not r["skip_fields"])
    verdict["residue_report_unmeasured_zero"] = {
        "runs": rr,
        "affected": [r["manifest"] for r in rr if r["motif_rate_is_unmeasured_zero"]
                     or r["chimaera_rate_is_unmeasured_zero"]],
        "fields_silently_absent": {r["manifest"]: r["fields_silently_absent"]
                                   for r in rr if r["fields_silently_absent"]},
        "dmnd_dropped": [r["manifest"] for r in rr
                         if set(r["dmnd_requested"]) - set(r["closest_fields"])],
        "verdict": "see affected: a rate of 0.0 written where the test never ran (motif: "
                   "Pfam absent or --skip-fields; chimaera: reference absent)"
                   if any(r["motif_rate_is_unmeasured_zero"]
                          or r["chimaera_rate_is_unmeasured_zero"] for r in rr) else "none"}

    # 5. divergence ladder: an index skipped because absent
    sp = []
    for f, d in by["spike_holdout.py"]:
        p = d["parameters"]
        if p.get("abundance_only"):
            continue
        ib = set(d.get("index_bytes", {}))
        want = [p["silva"], p["identified"], p["unidentified"], *p["dmnd"]]
        sp.append({"manifest": f, "missing": [w for w in want if w not in ib]})
    verdict["spike_ladder_index_skipped"] = {
        "runs": sp, "affected": [s["manifest"] for s in sp if s["missing"]],
        "verdict": ("every committed ladder run consulted every requested index it "
                    "recorded" if not any(s["missing"] for s in sp) else "see affected"),
        "limit": "Kraken2 databases are not recorded in this script's index_bytes, so "
                 "their presence cannot be checked from committed output"}

    # 6. index caches keyed by stem: same work dir + stem, different reference
    cache = defaultdict(set)
    for f, d in ms:
        p = d["parameters"]
        w = p.get("work")
        if not w:
            continue
        for k in ("selfasm", "observe", "chimaera_ref", "cross_sample_ref", "identified",
                  "unidentified", "rrna_contigs", "identified_contigs",
                  "unidentified_contigs"):
            v = p.get(k)
            for x in (v if isinstance(v, list) else [v] if v else []):
                cache[(w, Path(x).stem)].add(x)
        if p.get("silva"):
            cache[(w, "silva")].add(tuple(sorted(p["silva"])) if isinstance(
                p["silva"], list) else p["silva"])
    coll = {f"{w} :: {s}": sorted(map(str, v)) for (w, s), v in cache.items() if len(v) > 1}
    verdict["stale_index_cache"] = {
        "affected": coll,
        "verdict": ("no two committed runs sharing a work directory used different "
                    "references under one cache name" if not coll else "see affected"),
        "limit": "a reference regenerated in place under the same path is invisible here"}

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "audit.json").write_text(json.dumps(
        {"script": "experiments/fail_loudly_audit.py", "seed": args.seed,
         "manifests_read": len(ms), "verdict": verdict}, indent=2))
    for k, v in verdict.items():
        print(f"{k}: {v['verdict']}")
        if v.get("affected"):
            print(f"    affected: {v['affected']}")


if __name__ == "__main__":
    main()
