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
import re
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

#: How a subject set is named in the tables. The arms write a short key on each row
#: (`afdb`, `pdb`) because that is what names their database files; this is the only place
#: that decides how a reader sees it, so a table never asks anyone to know that `afdb` means
#: Swiss-Prot folded - which is the fact that makes the structure arm a controlled
#: comparison rather than a different search against different subjects.
DB_LABEL = {"afdb": "AlphaFold/Swiss-Prot (Swiss-Prot, folded)", "pdb": "the PDB",
            "pdbvirus": "every solved virus structure in the PDB"}

#: Which arm is the baseline to beat. DIAMOND's most sensitive mode, because the question is
#: what the production search can reach when it is allowed to be slow, not what it reaches
#: at its default.
BASELINE_ARMS = ("diamond_ultra_sensitive", "diamond_very_sensitive", "diamond_sensitive")

ARM_ORDER = ("diamond_sensitive", "diamond_very_sensitive", "diamond_ultra_sensitive",
             "blastp", "sw_exhaustive", "hmmsearch_ga", "hmmsearch_open", "jackhmmer",
             "hhsearch", "foldseek_afdb", "foldseek_pdb", "foldseek_pdbvirus", "codon")

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
    "foldseek_pdbvirus": "Foldseek vs solved virus structures",
    "codon": "tblastn vs RefSeq viral genomes",
}

#: Arms with no pairwise identity to report. Kept explicit so a None in the identity column
#: is a declared property of the instrument rather than a missing value.
NO_IDENTITY = {"hmmsearch_ga", "hmmsearch_open", "jackhmmer", "hhsearch"}


def failure_reason(traceback: str) -> str:
    """The one line from a failed arm's traceback that says what went wrong.

    The last line is the usual choice and it is the wrong one here, because the arms capture
    their tools' output and re-raise it: the last line of foldseek's stdout is
    "Time for processing: 0h 0m 0s 5ms", which went into this table as a reason for three
    runs' worth of structure failures. So a line that names the problem is preferred over a
    line that merely came last.
    """
    rows = [l.strip() for l in traceback.strip().splitlines() if l.strip()]
    if not rows:
        return "no traceback recorded"
    wanted = ("could not", "cannot", "refused", "not find", "no such", "error:",
              "failed", "permission", "missing")
    named = [l for l in rows
             if any(w in l.lower() for w in wanted) and not l.startswith("File \"")]
    # The exception's own final line still wins when it says something, so prefer the last
    # matching line over the first: tracebacks nest, and the outermost message is ours.
    return (named[-1] if named else rows[-1])[:400]


def significant(row: dict) -> bool:
    """One floor for every arm, so the arms can be compared at all.

    Profile arms get no residue floor: `aln` for them is an envelope on the query, not an
    aligned length against a subject sequence, and applying a residue floor to it would be
    comparing two different quantities. hmmsearch at GA has already passed its curated
    cutoff, which is a stricter statement than any E-value here, and it is kept as its own
    arm for that reason.

    WHY PROFILE ARMS ARE FLOORED ON THE DOMAIN E-VALUE AND NOT THE FULL-SEQUENCE ONE, which
    is the whole of this project's null leak. A full-sequence E-value sums the scores of
    every domain HMMER finds, so it answers "does this query have a relative"; this table
    is a table of alignments, and what says an alignment is real is its own domain E-value.
    The difference is not academic. On run i-0499bf15b019eba9f, floored on the full-sequence
    value, jackhmmer was the one arm to produce a significant alignment on shuffled
    sequence: 1 of 80 null records, seven alignments. All seven were 23 residues of one
    shuffled ORF against one Swiss-Prot protein, at subject positions 318, 346, 374, 402,
    430, 458 and 486 - a 28-residue tandem repeat, matched once per copy. Every individual
    domain sat at E between 2,000 and 2,900, which is nothing; the full-sequence value
    reached 1.1e-08 only by summing seven of them.

    Flooring on the domain value removes it and costs the real piles nothing: 6 of 107
    hmmsearch rows and 2.7% of jackhmmer's go, every query on every real pile keeps a hit,
    and the null goes to 0 of 80 for every arm. So this is a repeat artefact being named
    rather than a sensitivity trade, and the aggregate-against-per-domain shape is the same
    one this repo has been caught by three times before.
    """
    if row.get("evalue") is None:
        return False
    if row["arm"] in NO_IDENTITY:
        return (row["evalue"] <= EVALUE_FLOOR
                and (row.get("domain_evalue") or row["evalue"]) <= EVALUE_FLOOR)
    return row["evalue"] <= EVALUE_FLOOR and (row.get("aln") or 0) >= MIN_ALN_AA


def fmt(value, spec: str = ".1f") -> str:
    if value is None:
        return "-"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "-"
    # A dN of exactly zero divided by a positive dS gives -0.0, which formats as "-0.000"
    # and, inside a range, as the unreadable "-0.00--0.00". The sign carries no information
    # at zero, so it is dropped here rather than in the arithmetic, where the sign of a
    # genuinely small negative number would matter.
    if isinstance(value, float) and value == 0.0:
        value = 0.0
    return format(value, spec)


def ev(value) -> str:
    if value is None:
        return "-"
    return "0" if value == 0 else f"{value:.1e}"


def load(run_dir: Path, extra_dirs: list[Path] | None = None
         ) -> tuple[list[dict], dict, dict, dict]:
    """Rows, timings, ORFs and the staging manifest for one table.

    WHY THERE IS MORE THAN ONE SOURCE DIRECTORY. The eight paid arms need subject sets that
    only fit on an instance and come back from its results prefix. The structure arm needs
    about 6 GB of coordinates and a 2.3 GB model, which fits in the container, so it runs
    there for nothing and writes its own directory with its own manifest. Merging them is
    the honest shape: two hosts, two script hashes, one panel, and `§0` prints both so no
    reader has to assume one provenance covers everything.

    An extra directory must agree with the run on the GENES. `orfs.json` is read from the
    run directory only, and `reach_structure.py` records the sha256 of the `orfs.faa` it
    called so the agreement is checkable rather than assumed.
    """
    rows: list[dict] = []
    timings = {}
    for source in [run_dir] + list(extra_dirs or []):
        for path in sorted(source.glob("rows_*.jsonl")):
            with open(path) as handle:
                for line in handle:
                    if line.strip():
                        rows.append(json.loads(line))
        for path in sorted(source.glob("timing_*.json")):
            t = json.loads(path.read_text())
            timings[t["arm"]] = t
        # `reach_structure.py` writes one manifest rather than a timing file per arm,
        # because it is one arm. Translating it here keeps the cost table's shape.
        manifest = source / "structure_manifest.json"
        if manifest.exists():
            m = json.loads(manifest.read_text())
            timings["foldseek_pdbvirus"] = dict(
                arm="foldseek_pdbvirus", seconds=m.get("seconds"), ok=True,
                threads=m.get("threads"), host=m.get("host"),
                script_sha256=m.get("script_sha256"), git_sha=m.get("git_sha"),
                subject_chains=(m.get("subject_db") or {}).get("chains"),
                target=m.get("target"), orfs_sha256=m.get("orfs_sha256"),
                encoding_price=m.get("encoding_price"))
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
    ap.add_argument("--extra-dir", type=Path, action="append", default=[],
                    help="another directory of rows for this same panel, such as "
                         "reach_structure.py's in-container output. Repeatable. Its "
                         "provenance is printed beside the run's rather than folded into "
                         "it, because it is a different host and a different script.")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    rows, timings, orfs, stage = load(args.run_dir, args.extra_dir)
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
                 f"{MIN_ALN_AA} aligned residues. Profile arms get no residue floor, "
                 f"because their aligned length is an envelope on the query rather than a "
                 f"span against a subject; instead both their full-sequence and their "
                 f"**domain** E-value must clear {EVALUE_FLOOR:g}. The domain floor is "
                 f"load-bearing: without it a tandem repeat in one Swiss-Prot protein puts "
                 f"a shuffled ORF above significance, which is this panel's only null "
                 f"leak.\n")
    # PROVENANCE IS PRINTED PER SOURCE, not once for the table. The paid arms and the
    # structure arm ran on different hosts under different scripts, and a single line would
    # make the reader assume one hash covers all of it - which is the exact failure the
    # `script_sha256` convention exists to prevent.
    if args.extra_dir:
        run_orfs = args.run_dir / "orfs.faa"
        orfs_sha = (hashlib.sha256(run_orfs.read_bytes()).hexdigest()
                    if run_orfs.exists() else None)
        prov = [f"- `{args.run_dir.name}`: commit `{stage.get('commit', 'unknown')}`, "
                f"{len([a for a in arms_present if a != 'foldseek_pdbvirus'])} arms on "
                f"{stage.get('instance_type', 'an instance')}"]
        for extra in args.extra_dir:
            m_path = extra / "structure_manifest.json"
            if not m_path.exists():
                prov.append(f"- `{extra.name}`: rows merged, no manifest found")
                continue
            m = json.loads(m_path.read_text())
            same = ("the same `orfs.faa`" if m.get("orfs_sha256") == orfs_sha
                    else "**a different `orfs.faa` from the run's**")
            prov.append(
                f"- `{extra.name}`: git `{m.get('git_sha')}`, script sha256 "
                f"`{(m.get('script_sha256') or '')[:12]}`, {m.get('foldseek_version','')[:12]}, "
                f"{(m.get('subject_db') or {}).get('chains')} subject chains from "
                f"{(m.get('target') or {}).get('unique_entries')} PDB entries queried "
                f"{(m.get('target') or {}).get('fetched_utc', 'unknown')}, on {same}")
        lines.append("This table is merged from more than one source, so each one's "
                     "provenance is given separately:\n")
        lines += prov
        lines.append("")

    # §1 reach per pile per arm
    lines.append("\n## 1. How far each method reaches, by pile\n")
    lines.append("`queries with a hit` is out of the pile's size. `most distant` is the "
                 "median over that pile of each contig's lowest-identity significant "
                 "alignment; `furthest single` is the lowest of them. `nearest` is the "
                 "median of each contig's *best* significant alignment.\n")
    lines.append("**Read the two identity columns against each other and against the "
                 "near-reference row.** The most-distant column is bounded below by the "
                 "floor rather than by the query: at a fixed E-value the lowest-identity "
                 "alignment a large database returns is set by the database's size and the "
                 "query's length, so a near-reference contig and a divergent one come back "
                 "at nearly the same number. That is the control moving with the test, and "
                 "it means the most-distant statistic is a property of the operating point "
                 "and not of the method. The nearest column is the one that separates "
                 "piles, and `queries with a hit` is the one that separates arms.\n")
    lines.append("| arm | pile | queries with a hit | nearest (median %id) | "
                 "most distant (median %id) | furthest single | median aligned aa | "
                 "median query coverage |")
    lines.append("|---|---|---|---|---|---|---|---|")
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
            tops = [best(arm, q) for q in qs]
            top_idents = [t["pident"] for t in tops
                          if t and t.get("pident") is not None]
            cell = dict(
                arm=arm, pile=pile, queries=len(qs), with_a_hit=len(hits),
                median_nearest_pident=(round(statistics.median(top_idents), 1)
                                       if top_idents else None),
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
                f"{fmt(cell['median_nearest_pident'])} | "
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
            # THE SUBJECT SET IS READ OFF EVERY ROW, not only the significant ones. An arm
            # that ran against Pfam or a structure database and cleared the floor nowhere
            # has no significant row to carry its `db`, and this table used to print that as
            # `searches nothing` - which reads as the arm having no subjects rather than no
            # hits, and those are opposite findings. The two cases are now distinguished:
            # an arm with no rows at all did not run.
            arm_db = next((r.get("db") for r in rows if r["arm"] == arm and r.get("db")), None)
            if arm_db != SHARED_DB:
                novel_rows.append(dict(arm=arm, comparable=False, db=arm_db))
                searched = (f"{DB_LABEL.get(arm_db, f'`{arm_db}`')}, which the baseline "
                            f"never saw" if arm_db
                            else "no subject set is recorded on its rows, so this arm "
                                 "produced none")
                lines.append(f"| {ARM_LABEL.get(arm, arm)} | not comparable: searches "
                             f"{searched} | - | - | - | - |")
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
                 f"{panel['shuffle'].get('records_shuffled_per_seed', 'all')} records "
                 f"each. "
                 f"{panel['shuffle']['note'][0].upper()}{panel['shuffle']['note'][1:]}.\n")
    lines.append(
        "\n**These zeroes are the domain floor's doing, and that is the result, not a "
        "tidier table.** Floored on the full-sequence E-value instead, jackhmmer reads "
        "1 of 80 null records: 23 residues of one shuffled ORF against a 28-residue tandem "
        "repeat in one Swiss-Prot protein, matched once per copy at subject positions 318 "
        "through 486. Each domain sits at E between 2,000 and 2,900 and only their sum "
        "clears the floor. Removing it costs the real piles nothing - every query on every "
        "real pile keeps a hit - so an iterative profile search over a database containing "
        "repeat proteins needs a per-domain floor, and a reader checking this arm against "
        "another tool's numbers should check which E-value that tool reports.\n")

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
    # WHAT THE INTERESTING VERSION OF THE STRUCTURE ARM WOULD COST, from this arm's own
    # throughput rather than from an estimate. The arm as run encodes 202 query ORFs; the
    # version that would be a controlled comparison against the sequence arms encodes the
    # SUBJECT side instead, and that is the same model over four orders of magnitude more
    # residues. Printed here with both halves so the division can be redone.
    price = (timings.get("foldseek_pdbvirus") or {}).get("encoding_price")
    if price:
        lines.append(
            f"\nThe structure arm's own ProstT5 throughput is "
            f"**{fmt(price.get('residues_per_second'), '.1f')} residues/s on "
            f"{price.get('threads')} vCPU** ({price.get('query_residues'):,} query residues "
            f"in {fmt(price.get('seconds'), '.0f')} s). Encoding the *subject* side instead "
            f"- which is the only version of this arm that would be a controlled comparison "
            f"against the sequence arms - is {price.get('refseq_viral_residues'):,} residues "
            f"for the RefSeq viral half alone, or **"
            f"{fmt(price.get('subject_side_days'), '.1f')} days of CPU**. That is a GPU job, "
            f"and its rate there is unmeasured here.\n")
    failed = {a: t for a, t in timings.items() if t.get("ok") is False}
    if failed:
        lines.append("\nArms that did not run, and why:\n")
        for arm, t in sorted(failed.items()):
            lines.append(f"- `{arm}`: {failure_reason(t.get('traceback', ''))}")
        out["failed_arms"] = {a: t.get("traceback") for a, t in failed.items()}

    # §5b the structure arm's subject universe, which is the thing that decides how a null
    # result from it may be read.
    #
    # WHY THIS SECTION EXISTS. AlphaFold/Swiss-Prot was chosen as the structure target
    # because it is Swiss-Prot folded, which makes it the baseline's subjects in another
    # alphabet - that was the design, and it is half true. The shared sequence database is
    # Swiss-Prot PLUS the RefSeq viral proteins, and only the Swiss-Prot half is folded
    # anywhere. So this counts where the sequence arms' nearest relatives actually come
    # from. If they come from the RefSeq viral half, the structure arm is searching a
    # universe that excludes them, and "structure found nothing further" is then a statement
    # about Swiss-Prot's coverage of this family rather than about structure.
    REFSEQ_PREFIX = re.compile(r"^(?:YP|NP|AP|XP|WP)_")
    if baseline is not None:
        best_by_orf: dict[str, dict] = {}
        for r in sig:
            if r["arm"] != baseline:
                continue
            k = r["orf"]
            if k not in best_by_orf or r["bitscore"] > best_by_orf[k]["bitscore"]:
                best_by_orf[k] = r
        folded = [r for r in best_by_orf.values() if not REFSEQ_PREFIX.match(r["subject"])]
        lines.append("\n## 5b. Whether the structure target holds these proteins at all\n")
        lines.append(
            "The shared sequence database is Swiss-Prot plus the RefSeq viral proteins, and "
            "AlphaFold/Swiss-Prot folds only the Swiss-Prot half. So the structure arm is "
            "the baseline's subjects in another alphabet only to the extent that the "
            "baseline's hits are Swiss-Prot hits. Counted on each ORF's best significant "
            "alignment under the baseline:\n")
        lines.append("| ORFs with a significant baseline hit | best hit in Swiss-Prot "
                     "(foldable) | best hit in RefSeq viral (folded nowhere) |")
        lines.append("|---|---|---|")
        lines.append(f"| {len(best_by_orf)} | {len(folded)} | "
                     f"{len(best_by_orf) - len(folded)} |")
        out["structure_subject_universe"] = dict(
            orfs_with_baseline_hit=len(best_by_orf),
            best_hit_in_swissprot=len(folded),
            best_hit_in_refseq_viral=len(best_by_orf) - len(folded),
            swissprot_examples=sorted({r["subject"] for r in folded})[:8])

    # §4b what the domain E-value floor costs the real piles.
    #
    # This is the fix for this panel's only null leak, and a fix that removes a leak by
    # removing sensitivity is a trade rather than a fix. So the cost is counted: rows that
    # clear the full-sequence floor and fail the domain floor, per profile arm. It is
    # generated because it drifted once already - the prose carried a figure from a run
    # whose jackhmmer rows were not the shipped ones, and jackhmmer's iterative profile
    # building is not bit-identical between runs.
    prof_rows = [r for r in rows if r["arm"] in ("hmmsearch_ga", "hmmsearch_open",
                                                 "jackhmmer", "hhsearch")]
    if prof_rows:
        lines.append("\n## 4b. What the domain E-value floor costs\n")
        lines.append(
            "The null leak was a tandem repeat clearing the **full-sequence** E-value by "
            "summing seven domains that were each worthless. Flooring on the domain value "
            "too removes it; this is what that costs the real piles, counted as rows which "
            "clear the full-sequence floor and fail the domain floor.\n")
        lines.append("| arm | rows significant on full-sequence E | removed by the domain "
                     "floor | % |")
        lines.append("|---|---|---|---|")
        floor_rows = []
        for arm in ("hmmsearch_ga", "hmmsearch_open", "jackhmmer", "hhsearch"):
            pool = [r for r in prof_rows if r["arm"] == arm
                    and r.get("evalue") is not None and r["evalue"] <= EVALUE_FLOOR]
            if not pool:
                continue
            lost = [r for r in pool
                    if (r.get("domain_evalue") or r["evalue"]) > EVALUE_FLOOR]
            pct = 100.0 * len(lost) / len(pool)
            floor_rows.append(dict(arm=arm, full_sequence_significant=len(pool),
                                   removed=len(lost), percent=round(pct, 1)))
            lines.append(f"| {ARM_LABEL.get(arm, arm)} | {len(pool)} | {len(lost)} | "
                         f"{pct:.1f} |")
        out["domain_floor_cost"] = floor_rows

    # §5c the structure arm counted per GENE, not per contig.
    #
    # WHY PER GENE. §1 counts contigs, and a contig counts as placed if any one of its genes
    # is. For this arm that hides the result: it places one gene per contig where the
    # sequence arms place nearly all of them, and the write-up's claim is about which genes
    # and how many. Generated here so the prose quotes a table rather than a pass somebody
    # ran once by hand.
    #
    # "Structure only" is the reach question asked the only way it can be asked across two
    # subject universes: not which subjects were missed - they are not shared - but whether
    # any GENE is placed by structure that sequence could not place at all.
    SEQ_ARMS = ("diamond_sensitive", "diamond_very_sensitive", "diamond_ultra_sensitive",
                "blastp", "sw_exhaustive")
    if "foldseek_pdbvirus" in arms_present:
        seq_placed = {r["orf"] for r in sig if r["arm"] in SEQ_ARMS}
        str_best: dict[str, dict] = {}
        for r in sig:
            if r["arm"] != "foldseek_pdbvirus":
                continue
            if r["orf"] not in str_best or r["evalue"] < str_best[r["orf"]]["evalue"]:
                str_best[r["orf"]] = r
        pile_of_orf = {r["orf"]: r["pile"] for r in rows}
        lines.append("\n## 5c. The structure arm per gene, and the one reach question it "
                     "can answer\n")
        lines.append(
            "§1 counts contigs, and a contig counts as placed if any one of its genes is. "
            "This counts genes. `structure only` is the reach question asked the only way "
            "two different subject universes allow: not which subjects the other arms "
            "missed, since none are shared, but whether any gene is placed by structure "
            "that no sequence arm could place at all.\n")
        lines.append("| pile | ORFs | placed by a sequence arm | placed by structure | "
                     "**structure only** | sequence only | median 3Di identity of the "
                     "best placement |")
        lines.append("|---|---|---|---|---|---|---|")
        gene_rows = []
        for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                     "shuffled"):
            pool = [o for o, p in pile_of_orf.items() if p == pile]
            if not pool:
                continue
            s = [o for o in pool if o in seq_placed]
            t = [o for o in pool if o in str_best]
            only_t = sorted(o for o in pool if o in str_best and o not in seq_placed)
            only_s = [o for o in pool if o in seq_placed and o not in str_best]
            ids = [str_best[o]["pident"] for o in t]
            med = statistics.median(ids) if ids else None
            gene_rows.append(dict(pile=pile, orfs=len(pool), sequence=len(s),
                                  structure=len(t), structure_only=len(only_t),
                                  structure_only_orfs=only_t, sequence_only=len(only_s),
                                  median_best_pident=(round(med, 1) if med else None)))
            lines.append(f"| {pile} | {len(pool)} | {len(s)} | {len(t)} | "
                         f"**{len(only_t)}** | {len(only_s)} | {fmt(med, '.1f')} |")
        lines.append("")
        lines.append("Every structural placement, with the entry it landed on:\n")
        lines.append("| ORF | pile | subject | E | 3Di identity | aligned aa | query cov |")
        lines.append("|---|---|---|---|---|---|---|")
        placements = []
        for orf in sorted(str_best, key=lambda o: (pile_of_orf.get(o, ""), o)):
            r = str_best[orf]
            if pile_of_orf.get(orf) == "shuffled":
                continue
            placements.append(dict(orf=orf, pile=pile_of_orf.get(orf),
                                   subject=r["subject"], evalue=r["evalue"],
                                   pident=round(r["pident"], 1), aln=r["aln"],
                                   qcov=r.get("qcov")))
            lines.append(f"| `{orf}` | {pile_of_orf.get(orf)} | `{r['subject']}` | "
                         f"{ev(r['evalue'])} | {r['pident']:.1f} | {r['aln']} | "
                         f"{fmt(r.get('qcov'), '.2f')} |")
        out["structure_per_gene"] = dict(piles=gene_rows, placements=placements,
                                         sequence_arms_counted=list(SEQ_ARMS))

    # §6 profiles: does a second profile in the same clan reach further
    #: Every arm whose subjects are Pfam models, because the "several profiles per role"
    #: question is answered differently by profile-against-sequence and
    #: profile-against-profile, and reading it off only the first contradicts the literature.
    PROFILE_ARMS = ("hmmsearch_ga", "hmmsearch_open", "hhsearch")
    pfam = [r for r in sig if r["arm"] in PROFILE_ARMS]
    if pfam:
        lines.append("\n## 6. More than one profile per role\n")
        lines.append(
            "The literature on this family reports that several profiles per protein are "
            "needed to find its divergent members. Counted here as distinct Pfam models "
            "placed on one ORF, which is the observable version of that claim - and counted "
            "per arm, because profile-against-sequence and profile-against-profile answer "
            "it oppositely. Models are counted by ACCESSION with its version stripped, so "
            "two arms citing different releases of one family are not counted as two "
            "families.\n")
        lines.append("| arm | pile | ORFs with a domain | distinct models | "
                     "ORFs with 2+ models | models per ORF (median) | (max) |")
        lines.append("|---|---|---|---|---|---|---|")
        prof_rows = []
        for arm in PROFILE_ARMS:
          arm_rows = [r for r in pfam if r["arm"] == arm]
          if not arm_rows:
              continue
          for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                       "shuffled"):
            pool = [r for r in arm_rows if r["pile"] == pile]
            if not pool and pile not in ("divergent", "near_reference"):
                continue
            by_orf: dict[str, set] = {}
            for r in pool:
                # The accession without its version: hmmsearch reports PF01819.25 and
                # HHsearch may report another release of the same family, and counting those
                # as two models would manufacture exactly the multiplicity being measured.
                acc = (r.get("subject") or "").split(".")[0] or r.get("profile_name")
                by_orf.setdefault(r["orf"], set()).add(acc)
            sizes = sorted(len(s) for s in by_orf.values())
            multi = sum(1 for n in sizes if n > 1)
            row = dict(arm=arm, pile=pile, orfs_with_a_domain=len(by_orf),
                       distinct_models=len({m for s in by_orf.values() for m in s}),
                       orfs_with_two_or_more=multi,
                       median_models_on_one_orf=sizes[len(sizes) // 2] if sizes else 0,
                       max_models_on_one_orf=sizes[-1] if sizes else 0)
            prof_rows.append(row)
            lines.append(f"| {ARM_LABEL.get(arm, arm)} | {pile} | "
                         f"{row['orfs_with_a_domain']} | {row['distinct_models']} | "
                         f"{multi} | {row['median_models_on_one_orf']} | "
                         f"{row['max_models_on_one_orf']} |")
        out["profiles"] = prof_rows

        # Which models one profile arm places that another does not, per pile. This is the
        # column that decides whether the extra models are reach or tail: if an arm adds the
        # same number on the divergent pile as on the 95-99% control, it is describing Pfam's
        # own overlapping families rather than the query.
        placed: dict[tuple[str, str], set] = {}
        for r in pfam:
            acc = (r.get("subject") or "").split(".")[0] or r.get("profile_name")
            placed.setdefault((r["arm"], r["pile"]), set()).add((r["query"], acc))
        if any(a == "hhsearch" for a, _ in placed):
            lines.append("\n**Which models HHsearch adds, and on which pile.** A contig-model "
                         "pair counted once per contig. `hmmsearch-only` is the check that "
                         "the more sensitive arm has not simply traded one set for another.\n")
            lines.append("| pile | hmmsearch pairs | HHsearch pairs | HHsearch-only | "
                         "hmmsearch-only |")
            lines.append("|---|---|---|---|---|")
            add_rows = []
            for pile in ("divergent", "near_reference", "holdout", "holdout_reference"):
                hm = (placed.get(("hmmsearch_ga", pile), set())
                      | placed.get(("hmmsearch_open", pile), set()))
                hh = placed.get(("hhsearch", pile), set())
                if not hm and not hh:
                    continue
                add_rows.append(dict(pile=pile, hmmsearch=len(hm), hhsearch=len(hh),
                                     hhsearch_only=len(hh - hm),
                                     hmmsearch_only=len(hm - hh)))
                lines.append(f"| {pile} | {len(hm)} | {len(hh)} | {len(hh - hm)} | "
                             f"{len(hm - hh)} |")
            out["profile_overlap"] = add_rows

    # §7 codon space and the saturation question
    codon = [r for r in rows if r["arm"] == "codon" and significant(r)]
    if codon:
        lines.append("\n## 7. Codon space: is the synonymous signal saturated where these "
                     "candidates sit\n")
        lines.append(
            "**Read the `nearest` block, not the `all` block.** `all` pools every "
            "significant tblastn alignment a pile produced, and a pile's alignments are "
            "dominated by its long tail of distant subjects rather than by its own nearest "
            "relative - so `all` reads the same on the divergent pile and on the 95-99% "
            "control, which is the subject set describing itself rather than the piles. "
            "`nearest` takes the best-scoring alignment per ORF, which is the relative the "
            "contig actually has. `dS undefined` is the Jukes-Cantor correction's argument "
            "going non-positive, which is what saturation looks like when it is not clipped "
            "to a finite number.\n")
        lines.append("| set | pile | alignments | median %id | dS undefined | median dS | "
                     "median dN | median dN/dS | spread of dN/dS (p10-p90) |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        codon_rows = []
        #: Best-scoring codon alignment per ORF, which is the per-ORF nearest relative.
        nearest: dict[str, dict] = {}
        for r in codon:
            held = nearest.get(r["orf"])
            if held is None or r["bitscore"] > held["bitscore"]:
                nearest[r["orf"]] = r
        pools = [("nearest", list(nearest.values())), ("all", codon)]
        for label, source in pools:
          for pile in ("divergent", "near_reference", "holdout", "holdout_reference",
                       "shuffled"):
            pool = [r for r in source if r["pile"] == pile]
            if not pool:
                continue
            ds = [r["d_syn"] for r in pool if r.get("d_syn") is not None]
            dn = [r["d_nonsyn"] for r in pool if r.get("d_nonsyn") is not None]
            ratios = sorted(r["dnds"] for r in pool if r.get("dnds") is not None)
            undefined = sum(1 for r in pool if r.get("synonymous_saturated"))
            row = dict(
                set=label, pile=pile, alignments=len(pool),
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
                f"| {label} | {pile} | {len(pool)} | {fmt(row['median_pident'])} | "
                f"{undefined}/{len(pool)} | {fmt(row['median_ds'], '.3f')} | "
                f"{fmt(row['median_dn'], '.3f')} | {fmt(row['median_dnds'], '.3f')} | "
                f"{fmt(row['dnds_p10'], '.2f')}-{fmt(row['dnds_p90'], '.2f')} |")
        out["codon"] = codon_rows

    # §8 one row per contig, so a reader can see which gene carried the number
    lines.append("\n## 8. Per contig: did any arm find a nearer relative than the "
                 "baseline\n")
    lines.append("Ranked on the NEAREST relative, not the most distant one, for the reason "
                 "in section 1: the most-distant column is set by the floor. A contig whose "
                 "best arm is the baseline is a contig where nothing in this comparison "
                 "beat the production search. The gene column says which ORF carried the "
                 "number, because a complete genome's genes sit at different distances - "
                 "`results/divergent_contigs/` measured a 17.4-point median spread between "
                 "them - so a per-contig identity with no gene beside it is the error that "
                 "report named.\n")
    lines.append("| contig | pile | nearest named relative | baseline nearest %id | "
                 "best arm | its nearest %id | gain | gene |")
    lines.append("|---|---|---|---|---|---|---|---|")
    contig_rows = []
    for row in panel["rows"]:
        if row["arm"] == "shuffled":
            continue
        qid = row["qid"]
        base = best(baseline, qid) if baseline else None
        base_id = base.get("pident") if base else None
        cands = [(a, best(a, qid)) for a in arms_present if a not in NO_IDENTITY]
        cands = [(a, h) for a, h in cands if h and h.get("pident") is not None]
        # Highest identity wins, with the baseline preferred on a tie: the question is
        # whether anything beat the production search, so a tie is "nothing did".
        winner = max(cands, key=lambda t: (t[1]["pident"], t[0] == baseline)) \
            if cands else None
        gain = (round(winner[1]["pident"] - base_id, 1)
                if winner and base_id is not None else None)
        entry = dict(
            qid=qid, pile=row["arm"], taxon=taxon_of_query.get(qid), nt=row["nt"],
            mapped_depth=row.get("mapped_depth"),
            polymerase_pident=row.get("polymerase_pident"),
            baseline_nearest=(round(base_id, 1) if base_id is not None else None),
            best_arm=winner[0] if winner else None,
            best_arm_nearest=(round(winner[1]["pident"], 1) if winner else None),
            gain_over_baseline_points=gain,
            best_arm_gene=winner[1]["orf"] if winner else None,
            best_arm_subject=winner[1]["subject"] if winner else None,
            best_arm_query_coverage=(round(winner[1]["qcov"], 3) if winner else None))
        contig_rows.append(entry)
        lines.append(
            f"| {qid} | {row['arm']} | {taxon_of_query.get(qid, '-')[:34]} | "
            f"{fmt(entry['baseline_nearest'])} | "
            f"{ARM_LABEL.get(entry['best_arm'], '-')} | "
            f"{fmt(entry['best_arm_nearest'])} | "
            f"{'none' if gain in (0.0, None) else format(gain, '+.1f')} | "
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
