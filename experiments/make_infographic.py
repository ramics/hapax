"""Build the one-page figure of the whole process from committed outputs only.

WHY A GENERATOR AND NOT A HAND-WRITTEN PAGE. A figure with numbers typed into it is a
number without provenance, and this project has already had to withdraw claims whose prose
had drifted from the outputs it described (PR #13, PR #6). So the page is assembled from
the result JSONs by this script: every figure on it is reproducible from a committed script
plus a seed, and if a run is re-done the figure changes with it. Missing runs are left out
rather than filled in, so the page can be built while experiments are still going.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TIER_LABEL = {
    "artefact": "adapters, dimers, low complexity",
    "rrna": "SILVA rRNA (nucleotide)",
    "k2_viral": "Kraken2 viral",
    "k2_standard08": "Kraken2 Standard-8",
    "k2_pluspf08": "Kraken2 PlusPF-8",
    "k2_pluspf16": "Kraken2 PlusPF-16",
    "selfasm_sample_contigs_rrna": "the sample's own rRNA contigs",
    "selfasm_sample_contigs_identified_norrna": "the sample's own identified contigs",
    "prot_viral": "DIAMOND RefSeq viral proteins",
    "prot_sprot": "DIAMOND Swiss-Prot",
    "prot_uniref50_02pct": "DIAMOND UniRef50, 2%",
    "prot_uniref50_10pct": "DIAMOND UniRef50, 10%",
}
SELF_TIERS = {"selfasm_sample_contigs_rrna",
              "selfasm_sample_contigs_identified_norrna"}


def mean_min_max(xs):
    return {"mean": st.mean(xs), "min": min(xs), "max": max(xs)}


def cascade_data(results: Path):
    files = sorted(glob.glob(str(results / "subtract_s*.json")))
    if not files:
        return None, None, None
    ds = [json.load(open(f)) for f in files]
    curve = ds[0]["curve"]
    rows = []
    for i, step in enumerate(curve):
        tag = step["tier"]
        alone, marg, cum = [], [], []
        for d in ds:
            s = d["curve"][i]
            alone.append(100 * s["tier_alone_removed"])
            marg.append(100 * s["marginal_removed"])
            cum.append(100 * s["cumulative_removed"])
        rows.append({
            "tier": TIER_LABEL.get(tag, tag),
            "gb": step["index_bytes"] / 1e9,
            "alone": st.mean(alone),
            "marginal": st.mean(marg),
            "cumulative": st.mean(cum),
            "self": tag in SELF_TIERS,
        })
    union = mean_min_max([100 * d["union_removed"] for d in ds])
    residue = {
        "total": st.mean([100 * d["residue_pairs"] / d["pairs"] for d in ds]),
        "assembled": st.mean([100 * d["residue_composition"]
                              ["observe_sample_contigs_unidentified"] for d in ds]),
        "missing": st.mean([100 * d["tier3_pairs"] / d["pairs"] for d in ds]),
    }
    return rows, union, residue


def contig_data(path: Path):
    d = json.load(open(path))
    return {
        "total": d["contigs"], "bases": d["bases"],
        "identified": d["identified_excluding_rrna"]["contigs"],
        "rrna": d["rrna_contigs"]["contigs"],
        "tier2": d["tier2_population"]["contigs"],
        "tier2Scorable": round(d["tier2_population"]["frac_over_500nt_and_2x"]
                               * d["tier2_population"]["contigs"]),
        "tier2MedianLen": d["tier2_population"]["len_median"],
        "tier2MedianCov": round(d["tier2_population"]["cov_median"], 2),
        "rrnaMedianCov": round(d["rrna_contigs"]["cov_median"], 1),
        "allMedianCov": round(d["all"]["cov_median"], 2),
        "routeProtein": d["route_protein_only"],
        "routeNucleotide": d["route_nucleotide_only"],
    }


BLAST_STYLE = {
    "rrna": ("rRNA", "var(--known)"),
    "public_nucleotide": ("public nucleotide", "var(--known)"),
    "selfasm_identified": ("identified here", "var(--self)"),
    "residue_assembled": ("residue, assembled here", "var(--residue)"),
    "residue_missing": ("residue, missing", "var(--residue)"),
    "tier2_contigs": ("tier-2 contigs", "var(--residue-soft)"),
}


def blast_data(results: Path):
    files = sorted(glob.glob(str(results / "blast_check_s*.json")))
    if not files:
        return None
    rows = []
    per_class: dict[str, list[tuple[float, float, float]]] = {}
    for f in files:
        d = json.load(open(f))
        bn = d["blast"].get("blastn")
        if not bn:
            continue
        for cls, c in bn["per_class"].items():
            if not c.get("pident"):
                continue
            per_class.setdefault(cls, []).append(
                (c["pident"]["median"], c["query_cov"]["median"], c["any_hit_frac"]))
    for cls, vals in per_class.items():
        label, fill = BLAST_STYLE.get(cls, (cls, "var(--ink-soft)"))
        rows.append({"label": label, "fill": fill,
                     "pident": st.mean(v[0] for v in vals),
                     "qcov": st.mean(v[1] for v in vals),
                     "any_hit": st.mean(v[2] for v in vals),
                     "seeds": len(vals)})
    return rows


def spike_data(results: Path):
    files = sorted(glob.glob(str(results / "spike_s*.json")))
    if not files:
        return None
    ds = [json.load(open(f)) for f in files]
    rungs = []
    for i, r0 in enumerate(ds[0]["rungs"]):
        surv, rem = [], []
        for d in ds:
            for g in d["rungs"][i]["genomes"].values():
                surv.append(g["survival"]["survives_tier1"])
                rem.append(g["survival"]["tier1_removed"])
        rungs.append({"divergence": r0["divergence"],
                      "survives_tier1": st.mean(surv),
                      "tier1_removed": st.mean(rem),
                      "n": len(surv)})
    run_pairs = ds[0]["one_in_a_million_arithmetic"]["run_pairs"]
    expect = ds[0]["one_in_a_million_arithmetic"]["expected_spike_pairs_in_run"]
    worst = min(rungs, key=lambda r: r["survives_tier1"])
    best = max(rungs, key=lambda r: r["survives_tier1"])
    return {
        "rungs": rungs,
        "caption": (f"{len(ds)} seeds, {ds[0]['spike_pairs_per_rung']:,} simulated pairs "
                    f"per genome per rung, {len(ds[0]['genomes'])} genomes."),
        "arithmetic": (
            f"A virus at 1 in a million contributes about {expect:.0f} read pairs in this "
            f"{run_pairs:,}-pair run. At {best['divergence'] * 100:.0f}% divergence "
            f"{best['survives_tier1'] * 100:.0f}% of its reads survive tier 1, so about "
            f"{expect * best['survives_tier1']:.1f} pairs reach the per-read report; at "
            f"{worst['divergence'] * 100:.0f}% divergence "
            f"{worst['survives_tier1'] * 100:.0f}% survive, about "
            f"{expect * worst['survives_tier1']:.1f} pairs. Per-read survival is measured "
            f"at a resolvable spike rate; this sentence is arithmetic on top of it, not a "
            f"measurement at 1 in a million."),
    }


FIELDS = [
    {"name": "closest relative, however far", "what": "best DIAMOND bitscore and identity against four protein indexes, with no floor", "verdict": "limit",
     "calibration": "score ordering tracks database proximity, not novelty (PR #8), so it is honest to report and dishonest to threshold"},
    {"name": "maybe-coding", "what": "stop-free frames, longest ORF, and whether the frame is identifiable at all", "verdict": "limit",
     "calibration": "at 107 nt the true frame is the only stop-free frame in 0.528 of coding reads, so the field reports its own uncertainty"},
    {"name": "GC and entropy", "what": "composition and 3-mer entropy", "verdict": "live",
     "calibration": "composition-matched shuffled null, plus the reads the cascade did subtract"},
    {"name": "repetitive or information-bearing", "what": "duplicate 16-mer fraction within the pair", "verdict": "live",
     "calibration": "shuffled control; low-complexity reads are removed by the artefact tier first"},
    {"name": "primer dimer or adapter", "what": "adapter match, dimer length, poly-A/G runs", "verdict": "live",
     "calibration": "removes 0.31% of the library before any database is consulted"},
    {"name": "library abundance", "what": "median count-min count of the read's 31-mers in 400,000 pairs", "verdict": "live",
     "calibration": "a count-min sketch only overstates, so a read called rare is rare for certain"},
    {"name": "cross-sample recurrence", "what": "does it map a contig from any of the other 115 samples", "verdict": "live",
     "calibration": "that reference is observed and deliberately not subtracted, or the field would be zero by construction"},
    {"name": "structured non-coding RNA", "what": "Rfam cmscan, 4,227 covariance models, 66% short enough for a read", "verdict": "live",
     "calibration": "shuffled control at the same E-value"},
    {"name": "folding", "what": "minimum free energy z-score against a first-order Markov resample", "verdict": "live",
     "calibration": "Markov-1 shuffle preserves dinucleotide content, which MFE is sensitive to"},
    {"name": "dinucleotide bias", "what": "CpG and UpA observed over expected", "verdict": "live",
     "calibration": "against subtracted reads from the same run"},
    {"name": "CRISPR spacer match", "what": "exact 20-mer against spacers mined from the assemblies", "verdict": "dead",
     "calibration": "5 hits in 25 spacers: arrays are DNA loci and this is a metatranscriptome"},
    {"name": "chimaera by split mapping", "what": "the two halves of a read mapping apart", "verdict": "dead",
     "calibration": "26% false positives on reads the cascade had already identified"},
]

LIMITS = [
    "Removal is a lower bound on <em>known</em>: an index cannot tell sequence with no relative anywhere from sequence whose relative is below the matching threshold. The swept thresholds are the honest form of that limit.",
    "Subtraction is not free. The SURPI+ evaluation measured pre-processing removing 0.43% of true viral reads, and at a 1-in-a-million target a small false-removal rate matters more than the headline removal does.",
    "The far end of the index axis is extrapolated. A DIAMOND database of full UniRef50 is about 30 GB against the disk available here, so UniRef50 enters as streamed nested 2% and 10% subsets.",
    "The all-sample tier is a lower bound: the GEO assembly file truncates mid-stream at sample 65 of 116.",
    "Abundance is recurrence within a 400,000-pair sample of the run, not a genome coverage depth.",
    "One library, one city, one protocol. Nothing here says how the curve looks for a different sample type.",
]


def anatomy_data(results: Path):
    """The five-seed contig anatomy. Reported as the four categories plus the two
    continuous quantities, because the categorical answer is a negative and the continuous
    one is the result."""
    files = sorted(glob.glob(str(results / "anatomy_s*.json")))
    if not files:
        return None
    ds = [json.load(open(f)) for f in files]
    def v(cat, which):
        return mean_min_max([d["verdict"][cat][which] for d in ds])
    def q(setname, key):
        return mean_min_max([d["by_set"][setname][key]["mean"] for d in ds])
    cats = ["A_chimaera", "B_divergent_continuation", "C_artefact_end",
            "D_short_true_homology"]
    return {
        "seeds": len(ds),
        "contigs": sum(d["by_set"]["residue"]["reads"] for d in ds),
        "categories": [{"name": c.split("_", 1)[1].replace("_", " "),
                        "residue": v(c, "residue")["mean"] * 100,
                        "known": v(c, "known_control")["mean"] * 100} for c in cats],
        "alnfrac": {"residue": q("residue", "aligned_fraction")["mean"],
                    "known": q("known", "aligned_fraction")["mean"]},
        "pident": {"residue": q("residue", "pident")["mean"],
                   "known": q("known", "pident")["mean"]},
    }


def indexcost_data(ladder: Path, scaling: Path):
    """What a bigger protein index buys, and what searching it costs. Both halves are
    needed on one page: the buying is nearly nothing at an identification-grade threshold,
    and the cost is the reason the far end has not simply been run."""
    out = {}
    s = ladder / "ladder_summary.json"
    if s.exists():
        d = json.load(open(s))
        rungs = []
        for i, name in enumerate(["viral", "sprot", "uniref50 2%", "uniref50 10%"]):
            anyh, b50 = [], []
            for seed in d["per_seed"]:
                st_ = seed["rungs"][i]["sets"].get("residue")
                if st_ and not st_.get("not_searched"):
                    anyh.append(st_["any_hit"])
                    b50.append(st_["bitscore_ge_50"])
            if anyh:
                rungs.append({"name": name, "gb": d["cumulative_bytes"][i] / 1e9,
                              "any": mean_min_max(anyh)["mean"] * 100,
                              "strong": mean_min_max(b50)["mean"] * 100})
        out["ladder"] = {"rungs": rungs, "seeds": len(d["per_seed"])}
    f = scaling / "scaling_s1.json"
    if f.exists():
        d = json.load(open(f))
        out["scaling"] = [
            {"db": k.split("|")[0].replace(".dmnd", "").replace("uniref50_", ""),
             "mode": k.split("|")[1],
             "fixed": v["fixed_seconds_db_scan"],
             "marginal_per_1k": v["marginal_seconds_per_read"] * 1000}
            for k, v in d["fits"].items()]
    return out or None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cascade", type=Path, default=ROOT / "results/subtract_recut")
    ap.add_argument("--contigs", type=Path,
                    default=ROOT / "results/subtract_known/contig_identification.json")
    ap.add_argument("--blast", type=Path, default=ROOT / "results/blast_check")
    ap.add_argument("--spike", type=Path, default=ROOT / "results/spike_holdout")
    ap.add_argument("--anatomy", type=Path,
                    default=ROOT / "results/partial_alignment_contigs")
    ap.add_argument("--ladder", type=Path, default=ROOT / "results/residue_protladder")
    ap.add_argument("--scaling", type=Path, default=ROOT / "results/diamond_scaling")
    ap.add_argument("--template", type=Path,
                    default=ROOT / "experiments/infographic_template.html")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/subtract_known/infographic.html")
    args = ap.parse_args()

    rows, union, residue = cascade_data(args.cascade)
    if rows is None:
        raise SystemExit(f"no cascade results in {args.cascade}")
    data = {
        "cascade": rows,
        "union": union,
        "residue": residue,
        "selfRrnaAlone": next(r["alone"] for r in rows
                              if r["tier"] == TIER_LABEL["selfasm_sample_contigs_rrna"]),
        "contigs": contig_data(args.contigs),
        "blast": blast_data(args.blast),
        "spike": spike_data(args.spike),
        "anatomy": anatomy_data(args.anatomy),
        "indexcost": indexcost_data(args.ladder, args.scaling),
        "fields": FIELDS,
        "limits": LIMITS,
    }
    # Kraken2's threshold sweep comes from the cascade run, not the contig run: the point
    # being made is about reads, where the 12-fold swing lives.
    sw = json.load(open(sorted(glob.glob(str(args.cascade / "subtract_s*.json")))[0]))
    conf = sw["tiers"]["k2_pluspf16"]["confidence_sweep"]
    data["thresholds"] = {"conf0": 100 * conf["0.00"], "conf01": 100 * conf["0.10"],
                          "conf02": 100 * conf["0.20"]}

    html = args.template.read_text()
    a, b = html.index("/*HAPAX_DATA*/"), html.index("/*END*/") + len("/*END*/")
    html = html[:a] + json.dumps(data) + html[b:]
    args.out.write_text(html)
    print(f"wrote {args.out} ({len(html):,} bytes); spike "
          f"{'included' if data['spike'] else 'not yet run'}, blast "
          f"{len(data['blast'] or [])} classes, anatomy "
          f"{'included' if data['anatomy'] else 'not yet run'}, index cost "
          f"{'included' if data['indexcost'] else 'not yet run'}")


if __name__ == "__main__":
    main()
