"""The reach table, generated from the instance's rows. Nothing here is typed by hand.

WHAT A "MOST DISTANT RELATIVE RECOVERED" IS, AND WHY IT NEEDS SAYING. Among the alignments
an arm reports for one query that clear a significance floor, the most distant relative is
the one with the lowest identity. That ordering is only meaningful because the floor is
applied first: identity on a 30 aa alignment is quantised and optimistic, so the lowest
identity in an unfiltered pile is whichever chance match happened to be worst, which is a
property of the database's size and not of the method's reach. The floor is therefore the
same for every sequence arm - `results/divergent_known/`'s operating point, E < 1e-5 over at
least 60 aligned residues - and the shuffled null is what says it is the right floor.

WHY THE PROFILE ARMS HAVE NO IDENTITY COLUMN. hmmsearch reports a profile, not a subject
sequence, so there is no pairwise identity to report and inventing one by aligning to the
profile's consensus would be a number with no provenance. Those arms are ranked on
significance and their reach is read as "did it place a domain where the sequence arms placed
nothing", which is a different and weaker claim. Stated rather than papered over.

WHY "FOUND SOMETHING THE BASELINE MISSED" IS COMPUTED ONLY INSIDE THE SHARED UNIVERSE. An
arm searching Pfam or PDB finds subjects the baseline could not have found, because the
baseline never saw them. Counting those as new reach would measure the subject set, which is
this project's standing confound - 1.4 MB of the sample's own contigs removed more than 272
MB of SILVA. So the difference is taken only over Swiss-Prot plus RefSeq viral, where every
arm had the same chance, and the other universes are reported as their own rows with no
difference computed.

READ THE CONTROL FIRST. Every table carries the near-reference arm beside the divergent one
and the shuffled null beside both, drawn by this same code from the same run. An arm that
reaches further on the test arm and no further on the control has extended reach; one that
moves on both is reading its own threshold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The operating point of `results/divergent_known/`, reused so this run's numbers sit on
#: the same scale as the committed baseline. Both halves matter: the E-value keeps chance
#: matches out and the length floor keeps the identity column readable.
EVALUE_FLOOR = 1e-5
MIN_ALN_AA = 60

#: The subject universe every sequence arm shares. A difference against the baseline is
#: computed here and nowhere else.
SHARED_DB = "swissprot+viral"

#: Which arm is the baseline to beat. DIAMOND's most sensitive mode, because the question is
#: what the production search can reach when it is allowed to be slow, not what it reaches
#: at its default.
BASELINE_ARMS = ("diamond_ultra_sensitive", "diamond_very_sensitive", "diamond_sensitive")

ARM_ORDER = ("diamond_sensitive", "diamond_very_sensitive", "diamond_ultra_sensitive",
             "blastp", "sw_exhaustive", "hmmsearch_ga", "hmmsearch_open", "jackhmmer",
             "hhsearch", "foldseek_afdb", "foldseek_pdb", "codon")

ARM_LABEL = {
    "diamond_sensitive": "DIAMOND --sensitive",
    "diamond_very_sensitive": "DIAMOND --very-sensitive",
    "diamond_ultra_sensitive": "DIAMOND --ultra-sensitive",
    "blastp": "BLASTP",
    "sw_exhaustive": "Smith-Waterman, exhaustive",
    "hmmsearch_ga": "hmmsearch, Pfam-A at GA",
    "hmmsearch_open": "hmmsearch, Pfam-A at E<=10",
    "jackhmmer": "jackhmmer, 5 iterations",
    "hhsearch": "HHsearch, profile vs Pfam profiles",
    "foldseek_afdb": "Foldseek vs AlphaFold/Swiss-Prot",
    "foldseek_pdb": "Foldseek vs PDB",
    "codon": "tblastn vs RefSeq viral genomes",
}

#: Arms with no pairwise identity to report. Kept explicit so a None in the identity column
#: is a declared property of the instrument rather than a missing value.
NO_IDENTITY = {"hmmsearch_ga", "hmmsearch_open", "jackhmmer", "hhsearch"}


def significant(row: dict) -> bool:
    """One floor for every arm, so the arms can be compared at all.

    Profile arms are floored on the E-value alone: `aln` for them is an envelope on the
    query, not an aligned length against a subject sequence, and applying a residue floor
    to it would be comparing two different quantities. hmmsearch at GA has already passed
    its curated cutoff, which is a stricter statement than any E-value here, and it is kept
    as its own arm for that reason.
    """
    if row.get("evalue") is None:
        return False
    if row["arm"] in NO_IDENTITY:
        return row["evalue"] <= EVALUE_FLOOR
    return row["evalue"] <= EVALUE_FLOOR and (row.get("aln") or 0) >= MIN_ALN_AA


def fmt(value, spec: str = ".1f") -> str:
    if value is None:
        return "-"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "-"
    return format(value, spec)


def ev(value) -> str:
    if value is None:
        return "-"
    return "0" if value == 0 else f"{value:.1e}"


def load(run_dir: Path) -> tuple[list[dict], dict, dict, dict]:
    rows: list[dict] = []
    for path in sorted(run_dir.glob("rows_*.jsonl")):
        with open(path) as handle:
            for line in handle:
                if line.strip():
                    rows.append(json.loads(line))
    timings = {}
    for path in sorted(run_dir.glob("timing_*.json")):
        t = json.loads(path.read_text())
        timings[t["arm"]] = t
    orfs = {r["orf"]: r for r in json.loads((run_dir / "orfs.json").read_text())} \
        if (run_dir / "orfs.json").exists() else {}
    stage = json.loads((run_dir / "stage_manifest.json").read_text()) \
        if (run_dir / "stage_manifest.json").exists() else {}
    return rows, timings, orfs, stage


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="the instance's output, as downloaded from its results prefix")
    ap.add_argument("--panel", type=Path, required=True, help="panel.json from the run")
    ap.add_argument("--out", type=Path, default=ROOT / "results/reach")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    rows, timings, orfs, stage = load(args.run_dir)
    panel = json.loads(args.panel.read_text())
    arm_of_query = {r["qid"]: r["arm"] for r in panel["rows"]}
    taxon_of_query = {r["qid"]: r["taxon"] for r in panel["rows"]}
    queries_by_pile: dict[str, list[str]] = {}
    for row in panel["rows"]:
        queries_by_pile.setdefault(row["arm"], []).append(row["qid"])

    # Every alignment is on an ORF; the panel's unit is a contig. A contig's reach is the
    # most distant relative any of its ORFs recovered, because a complete genome's genes sit
    # at different distances and `results/divergent_contigs/` measured a 17.4-point median
    # spread between them. Reporting one number per contig without saying which gene it came
    # from is the error that report named, so the gene is carried through.
    for row in rows:
        orf = orfs.get(row["orf"], {})
        row["query"] = orf.get("query") or row["orf"].split("__orf")[0]
        row["pile"] = arm_of_query.get(row["query"], "unknown")
        row["aa_len"] = orf.get("aa_len")

    arms_present = [a for a in ARM_ORDER if any(r["arm"] == a for r in rows)]
    sig = [r for r in rows if significant(r)]

    def most_distant(arm: str, query: str) -> dict | None:
        """The lowest-identity significant alignment for one contig under one arm.

        Ties and the no-identity arms both fall back to the weakest surviving E-value, which
        is the only distance proxy those instruments offer.
        """
        pool = [r for r in sig if r["arm"] == arm and r["query"] == query]
        if not pool:
            return None
        if arm in NO_IDENTITY:
            return max(pool, key=lambda r: r["evalue"])
        return min(pool, key=lambda r: (r["pident"], -r["evalue"]))

    def best(arm: str, query: str) -> dict | None:
        pool = [r for r in sig if r["arm"] == arm and r["query"] == query]
        return min(pool, key=lambda r: r["evalue"]) if pool else None

    # ---------------------------------------------------------------- the tables
    lines: list[str] = []
    out: dict = dict(
        unit="one contig", question=panel["question"], seed=args.seed,
        operating_point=dict(evalue_floor=EVALUE_FLOOR, min_alignment_aa=MIN_ALN_AA,
                             shared_subject_set=SHARED_DB,
                             baseline=next((a for a in BASELINE_ARMS
                                            if a in arms_present), None),
                             profile_arms_floored_on_evalue_only=sorted(NO_IDENTITY)),
        panel=panel["arms"], shuffle=panel["shuffle"], stage=stage, arms=arms_present,
    )
    baseline = out["operating_point"]["baseline"]

    lines.append("# Reach by method: the generated tables\n")
    lines.append(f"Generated by `experiments/report_reach.py` from "
                 f"`{args.run_dir.name}`. Seed {args.seed}. "
                 f"Floor for every sequence arm: E < {EVALUE_FLOOR:g} over at least "
                 f"{MIN_ALN_AA} aligned residues.\n")

    # §1 reach per pile per arm
    lines.append("\n## 1. How far each method reaches, by pile\n")
    lines.append("`queries with a hit` is out of the pile's size. `most distant` is the "
                 "median over that pile of each contig's lowest-identity significant "
                 "alignment; `furthest single` is the lowest of them.\n")
    lines.append("| arm | pile | queries with a hit | most distant (median %id) | "
                 "furthest single | median aligned aa | median query coverage |")
    lines.append("|---|---|---|---|---|---|---|")
    per_arm_pile: dict[tuple[str, str], dict] = {}
    for arm in arms_present:
        for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                     "shuffled"):
            qs = queries_by_pile.get(pile, [])
            if not qs:
                continue
            hits = [most_distant(arm, q) for q in qs]
            hits = [h for h in hits if h]
            idents = [h["pident"] for h in hits if h.get("pident") is not None]
            cell = dict(
                arm=arm, pile=pile, queries=len(qs), with_a_hit=len(hits),
                median_most_distant_pident=(round(statistics.median(idents), 1)
                                            if idents else None),
                furthest_single_pident=(round(min(idents), 1) if idents else None),
                median_aligned_aa=(round(statistics.median([h["aln"] for h in hits]))
                                   if hits else None),
                median_query_coverage=(round(statistics.median([h["qcov"] for h in hits]), 3)
                                       if hits else None))
            per_arm_pile[(arm, pile)] = cell
            lines.append(
                f"| {ARM_LABEL.get(arm, arm)} | {pile} | {len(hits)}/{len(qs)} | "
                f"{fmt(cell['median_most_distant_pident'])} | "
                f"{fmt(cell['furthest_single_pident'])} | "
                f"{fmt(cell['median_aligned_aa'], '.0f')} | "
                f"{fmt(cell['median_query_coverage'], '.2f')} |")
    out["reach_by_pile"] = list(per_arm_pile.values())

    # §2 the holdout, which is the only row with a known answer
    lines.append("\n## 2. The holdout: the one query whose distance nobody chose\n")
    lines.append("True nucleotide divergence 30.28%. The already-referenced companion from "
                 "the same library is the zero-distance end of the same instrument, so the "
                 "two columns together say whether an arm orders distance correctly at all."
                 "\n")
    lines.append("| arm | most distant relative | %id | aligned aa | query coverage | "
                 "E-value | best hit %id | companion best %id |")
    lines.append("|---|---|---|---|---|---|---|---|")
    holdout_rows = []
    for arm in arms_present:
        distant = most_distant(arm, "hold")
        top = best(arm, "hold")
        companion = best(arm, "leak")
        row = dict(arm=arm,
                   subject=distant["subject"] if distant else None,
                   pident=distant.get("pident") if distant else None,
                   aln=distant["aln"] if distant else None,
                   qcov=round(distant["qcov"], 3) if distant else None,
                   evalue=distant["evalue"] if distant else None,
                   best_pident=top.get("pident") if top else None,
                   companion_best_pident=companion.get("pident") if companion else None)
        holdout_rows.append(row)
        lines.append(
            f"| {ARM_LABEL.get(arm, arm)} | {row['subject'] or '-'} | "
            f"{fmt(row['pident'])} | {fmt(row['aln'], '.0f')} | {fmt(row['qcov'], '.2f')} | "
            f"{ev(row['evalue'])} | {fmt(row['best_pident'])} | "
            f"{fmt(row['companion_best_pident'])} |")
    out["holdout"] = holdout_rows

    # §3 what each arm found that the baseline did not, inside the shared universe only
    lines.append("\n## 3. What each arm found that the baseline did not\n")
    if baseline is None:
        lines.append("No baseline arm ran, so no difference can be taken.\n")
    else:
        lines.append(f"Baseline: **{ARM_LABEL[baseline]}**. Computed over `{SHARED_DB}` "
                     f"only, where every sequence arm had the same chance; an arm searching "
                     f"Pfam or PDB is searching subjects the baseline never saw, so no "
                     f"difference is taken for those and their rows say so.\n")
        lines.append("| arm | subjects the baseline missed | on test-arm contigs | "
                     "on near-reference contigs | on the null | most distant of them |")
        lines.append("|---|---|---|---|---|---|")
        base_pairs = {(r["query"], r["subject"]) for r in sig
                      if r["arm"] == baseline and r.get("db") == SHARED_DB}
        novel_rows = []
        for arm in arms_present:
            if arm == baseline:
                continue
            arm_rows = [r for r in sig if r["arm"] == arm]
            if not arm_rows or arm_rows[0].get("db") != SHARED_DB:
                novel_rows.append(dict(arm=arm, comparable=False,
                                       db=(arm_rows[0].get("db") if arm_rows else None)))
                lines.append(f"| {ARM_LABEL.get(arm, arm)} | not comparable: searches "
                             f"`{arm_rows[0].get('db') if arm_rows else 'nothing'}` | - | "
                             f"- | - | - |")
                continue
            extra = [r for r in arm_rows if (r["query"], r["subject"]) not in base_pairs]
            by_pile: dict[str, int] = {}
            for r in extra:
                by_pile[r["pile"]] = by_pile.get(r["pile"], 0) + 1
            with_id = [r for r in extra if r.get("pident") is not None]
            far = min(with_id, key=lambda r: r["pident"]) if with_id else None
            novel_rows.append(dict(
                arm=arm, comparable=True, extra_pairs=len(extra), by_pile=by_pile,
                most_distant_extra=(dict(query=far["query"], subject=far["subject"],
                                         pident=round(far["pident"], 1), aln=far["aln"],
                                         evalue=far["evalue"]) if far else None)))
            if far is None:
                furthest = "-"
            else:
                furthest = (f"{far['subject']} at {far['pident']:.1f}% "
                            f"over {far['aln']} aa")
            lines.append(
                f"| {ARM_LABEL.get(arm, arm)} | {len(extra)} | "
                f"{by_pile.get('divergent', 0)} | {by_pile.get('near_reference', 0)} | "
                f"{by_pile.get('shuffled', 0)} | {furthest} |")
        out["against_baseline"] = novel_rows

    # §4 the two controls, which decide every row above
    lines.append("\n## 4. The controls\n")
    lines.append("Five null piles, one per shuffle seed. The null column is the per-seed "
                 "range, not a pooled total: a null that reads zero once is one draw, and "
                 "the range is what makes it a result.\n")
    lines.append("| arm | significant alignments on the null, per seed (min-max) | "
                 "null contigs with one | near-reference contigs with one | "
                 "test contigs with one |")
    lines.append("|---|---|---|---|---|")
    null_seeds = panel.get("null_seeds") or [1]
    seed_of_query = {r["qid"]: r.get("null_seed") for r in panel["rows"]}
    control_rows = []
    for arm in arms_present:
        null_hits = [r for r in sig if r["arm"] == arm and r["pile"] == "shuffled"]
        per_seed = {sd: sum(1 for r in null_hits if seed_of_query.get(r["query"]) == sd)
                    for sd in null_seeds}
        row = dict(arm=arm, null_alignments=len(null_hits),
                   null_alignments_per_seed=per_seed,
                   null_alignments_min=min(per_seed.values()) if per_seed else 0,
                   null_alignments_max=max(per_seed.values()) if per_seed else 0,
                   null_queries_with_one=len({r["query"] for r in null_hits}),
                   null_queries=len(queries_by_pile.get("shuffled", [])),
                   near_reference_with_one=per_arm_pile.get(
                       (arm, "near_reference"), {}).get("with_a_hit"),
                   divergent_with_one=per_arm_pile.get(
                       (arm, "divergent"), {}).get("with_a_hit"))
        control_rows.append(row)
        span = (f"{row['null_alignments_min']}" if
                row['null_alignments_min'] == row['null_alignments_max']
                else f"{row['null_alignments_min']}-{row['null_alignments_max']}")
        lines.append(
            f"| {ARM_LABEL.get(arm, arm)} | {span} | "
            f"{row['null_queries_with_one']}/{row['null_queries']} | "
            f"{fmt(row['near_reference_with_one'], '.0f')}/"
            f"{len(queries_by_pile.get('near_reference', []))} | "
            f"{fmt(row['divergent_with_one'], '.0f')}/"
            f"{len(queries_by_pile.get('divergent', []))} |")
    out["controls"] = control_rows
    lines.append(f"\nThe null is {panel['shuffle']['implementation']}, fallback per seed "
                 f"{panel['shuffle']['fallback_records']} out of "
                 f"{panel['shuffle'].get('records_shuffled_per_seed')} records each. "
                 f"{panel['shuffle']['note']}.\n")

    # §5 cost
    lines.append("\n## 5. Cost\n")
    lines.append("| arm | wall clock (s) | alignments reported | seconds per query ORF |")
    lines.append("|---|---|---|---|")
    n_orfs = len({r["orf"] for r in rows}) or 1
    cost_rows = []
    for arm in arms_present:
        # DIAMOND's three modes share one timing entry, because they run in one call of the
        # arm. Attributing a third of it to each would be an invented number.
        key = "diamond" if arm.startswith("diamond") else (
            "hmmsearch" if arm.startswith("hmmsearch") else
            ("foldseek_afdb" if arm == "foldseek_afdb" else
             ("foldseek_pdb" if arm == "foldseek_pdb" else
              ("exhaustive" if arm == "sw_exhaustive" else arm))))
        t = timings.get(key, {})
        row = dict(arm=arm, timing_key=key, seconds=t.get("seconds"),
                   rows=sum(1 for r in rows if r["arm"] == arm),
                   seconds_per_orf=(round(t["seconds"] / n_orfs, 2)
                                    if t.get("seconds") else None),
                   ok=t.get("ok"))
        cost_rows.append(row)
        lines.append(f"| {ARM_LABEL.get(arm, arm)} | {fmt(row['seconds'], '.0f')} | "
                     f"{row['rows']} | {fmt(row['seconds_per_orf'], '.2f')} |")
    out["cost"] = cost_rows
    failed = {a: t for a, t in timings.items() if t.get("ok") is False}
    if failed:
        lines.append("\nArms that did not run, and why:\n")
        for arm, t in sorted(failed.items()):
            first = (t.get("traceback", "").strip().splitlines() or ["no traceback"])[-1]
            lines.append(f"- `{arm}`: {first}")
        out["failed_arms"] = {a: t.get("traceback") for a, t in failed.items()}

    # §6 profiles: does a second profile in the same clan reach further
    pfam = [r for r in sig if r["arm"] == "hmmsearch_ga"]
    if pfam:
        lines.append("\n## 6. More than one profile per role\n")
        lines.append("The literature on this family reports that several profiles per "
                     "protein are needed to find its divergent members. Counted here as "
                     "distinct Pfam models placed on one ORF, which is the observable "
                     "version of that claim.\n")
        lines.append("| pile | ORFs with a domain | distinct models | "
                     "ORFs with 2+ models | models per ORF (max) |")
        lines.append("|---|---|---|---|---|")
        prof_rows = []
        for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                     "shuffled"):
            pool = [r for r in pfam if r["pile"] == pile]
            if not pool and pile not in ("divergent", "near_reference"):
                continue
            by_orf: dict[str, set] = {}
            for r in pool:
                by_orf.setdefault(r["orf"], set()).add(r.get("profile_name") or r["subject"])
            multi = sum(1 for s in by_orf.values() if len(s) > 1)
            row = dict(pile=pile, orfs_with_a_domain=len(by_orf),
                       distinct_models=len({m for s in by_orf.values() for m in s}),
                       orfs_with_two_or_more=multi,
                       max_models_on_one_orf=max((len(s) for s in by_orf.values()),
                                                 default=0))
            prof_rows.append(row)
            lines.append(f"| {pile} | {row['orfs_with_a_domain']} | "
                         f"{row['distinct_models']} | {multi} | "
                         f"{row['max_models_on_one_orf']} |")
        out["profiles"] = prof_rows

    # §7 codon space and the saturation question
    codon = [r for r in rows if r["arm"] == "codon" and significant(r)]
    if codon:
        lines.append("\n## 7. Codon space: is the synonymous signal saturated where these "
                     "candidates sit\n")
        lines.append("One row per pile, over significant tblastn alignments. `dS undefined` "
                     "is the Jukes-Cantor correction's argument going non-positive, which "
                     "is what saturation looks like when it is not clipped to a finite "
                     "number.\n")
        lines.append("| pile | alignments | median %id | dS undefined | median dS | "
                     "median dN | median dN/dS | spread of dN/dS (p10-p90) |")
        lines.append("|---|---|---|---|---|---|---|---|")
        codon_rows = []
        for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                     "shuffled"):
            pool = [r for r in codon if r["pile"] == pile]
            if not pool:
                continue
            ds = [r["d_syn"] for r in pool if r.get("d_syn") is not None]
            dn = [r["d_nonsyn"] for r in pool if r.get("d_nonsyn") is not None]
            ratios = sorted(r["dnds"] for r in pool if r.get("dnds") is not None)
            undefined = sum(1 for r in pool if r.get("synonymous_saturated"))
            row = dict(
                pile=pile, alignments=len(pool),
                median_pident=round(statistics.median([r["pident"] for r in pool]), 1),
                ds_undefined=undefined,
                ds_undefined_share=round(undefined / len(pool), 3),
                median_ds=round(statistics.median(ds), 3) if ds else None,
                median_dn=round(statistics.median(dn), 3) if dn else None,
                median_dnds=round(statistics.median(ratios), 3) if ratios else None,
                dnds_p10=round(ratios[int(0.1 * (len(ratios) - 1))], 3) if ratios else None,
                dnds_p90=round(ratios[int(0.9 * (len(ratios) - 1))], 3) if ratios else None)
            codon_rows.append(row)
            lines.append(
                f"| {pile} | {len(pool)} | {fmt(row['median_pident'])} | "
                f"{undefined}/{len(pool)} | {fmt(row['median_ds'], '.3f')} | "
                f"{fmt(row['median_dn'], '.3f')} | {fmt(row['median_dnds'], '.3f')} | "
                f"{fmt(row['dnds_p10'], '.2f')}-{fmt(row['dnds_p90'], '.2f')} |")
        out["codon"] = codon_rows

    # §8 one row per contig, so a reader can see which gene carried the number
    lines.append("\n## 8. Per contig, under the baseline and under the arm that reached "
                 "furthest\n")
    lines.append("| contig | pile | nearest named relative | baseline most distant %id | "
                 "best arm | its most distant %id | gene |")
    lines.append("|---|---|---|---|---|---|---|")
    contig_rows = []
    for row in panel["rows"]:
        if row["arm"] == "shuffled":
            continue
        qid = row["qid"]
        base = most_distant(baseline, qid) if baseline else None
        cands = [(a, most_distant(a, qid)) for a in arms_present if a not in NO_IDENTITY]
        cands = [(a, h) for a, h in cands if h and h.get("pident") is not None]
        winner = min(cands, key=lambda t: t[1]["pident"]) if cands else None
        entry = dict(
            qid=qid, pile=row["arm"], taxon=taxon_of_query.get(qid), nt=row["nt"],
            mapped_depth=row.get("mapped_depth"),
            polymerase_pident=row.get("polymerase_pident"),
            baseline_most_distant=(round(base["pident"], 1)
                                   if base and base.get("pident") is not None else None),
            best_arm=winner[0] if winner else None,
            best_arm_most_distant=(round(winner[1]["pident"], 1) if winner else None),
            best_arm_gene=winner[1]["orf"] if winner else None,
            best_arm_subject=winner[1]["subject"] if winner else None)
        contig_rows.append(entry)
        lines.append(
            f"| {qid} | {row['arm']} | {taxon_of_query.get(qid, '-')[:34]} | "
            f"{fmt(entry['baseline_most_distant'])} | "
            f"{ARM_LABEL.get(entry['best_arm'], '-')} | "
            f"{fmt(entry['best_arm_most_distant'])} | "
            f"{(entry['best_arm_gene'] or '-').split('__')[-1]} |")
    out["per_contig"] = contig_rows

    out["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out["git_sha"] = stage.get("commit", "unknown")
    (args.out / "TABLES.md").write_text("\n".join(lines) + "\n")
    (args.out / f"reach_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {args.out / 'TABLES.md'} and reach_s{args.seed}.json")
    print(f"arms: {', '.join(arms_present)}")
    if failed:
        print(f"arms that failed: {', '.join(sorted(failed))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
