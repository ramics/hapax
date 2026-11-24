"""Order what a person should look at first, from quantities that were already measured.

WHAT THIS IS, AND WHAT IT IS CAREFUL NOT TO BE. Imo asked for a ranking of how much
attention a finding deserves. This project's founding constraint is that a search cannot
score novelty, because score ordering tracks database proximity - measured, repeatedly -
so nothing here scores novelty. What it does instead is order NAMED RELATIVES by how much
a person should look at them, out of four things that were each measured elsewhere and
each ship with their own calibration:

  how far      the median identity and aligned length of the match. A distance from a
               database, reported as a distance and never converted into a probability of
               being new.
  how much     how many residue read pairs point at the subject, with a Wilson interval
               widened by the measured block design effect.
  what it is   the coarse host group of the MATCHED RECORD, from the committed keyword
               table in host_groups.json. This is a statement about the known relative.
  how alone    whether anything the cascade already identified points at the same subject.
               A subject that only the residue points at is the interesting case.

IT IS A TRIAGE ORDER, NOT A RISK ASSESSMENT. It says which rows a human should read first
in an afternoon. It is not a probability of harm, not a clinical finding, and not a
statement that anything here infects anybody. Every card shows its own four components so
the order can be disagreed with row by row, and the formula is printed on the page rather
than buried here.

THE ONE PIECE OF REAL MODELLING, AND IT ENCODES A TENSION RATHER THAN HIDING IT.
Distance does two opposite things at once: the further a read is from its nearest known
relative, the more interesting it is, and the less that relative's host group can be
trusted to transfer. A score that multiplied "far" by "scary host" would be reporting the
product of a signal and its own unreliability. So the host term is damped towards the
unknown-host weight by a confidence that falls with distance and with short alignments:

    host_term = concern(group) * c + concern(unknown) * (1 - c)

with c built from identity and aligned length. At 95% identity over 100 aa a human-
associated relative keeps nearly all its weight; at 45% over 30 aa it keeps almost none
and the row is ranked as an unknown, which is what it is.

WHY IT IS A PIPELINE STAGE, HAVING ARGUED THAT IT WAS NOT. This docstring used to say that
keeping it out of `hapax.py` was the honest signal that no number in the run changes because
of it. That was right about the script and wrong about the report. It measures nothing and
changes no number - it reads a finished run directory and reorders it - but the run viewer's
unnamed-remainder section is built from its output, and a section with no stage behind it is
a hole in the tool rather than in the water: the first production pass over a library nobody
had touched produced a page whose summary sections were a census of absent files. So it is
stage `ranking`, last in STAGE_ORDER because it reads the run's own outputs, and it still
measures nothing: no number anywhere else in the run moves if it is skipped, which is why it
is also the one stage that costs seconds and consults no index.

WHAT ITS HOST COLUMN IS. A label on a database record, never a claim about the water. Over
41,829 residue subjects of the discriminating library host confidence maxes at 0.40, and on
the second library through production, SRR23702504, every one of its 58-83 ranked subjects
per seed is unknown_host with a maximum of 0.256 and a closest relative at 78.6% identity
over 52 aa. Nothing there is close enough to a named record for that record's host to
transfer, so the ceiling travels with the colours wherever they are drawn.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]
TABLE = Path(__file__).resolve().parent / "host_groups.json"

#: Identity at which the matched record's host group is believed in full, and the identity
#: at which it is believed not at all. Between them the confidence ramps linearly. The
#: lower end is set at 40% because the repo's own BLAST ground truth puts the chance band
#: for a short query at about that identity, so below it the name itself is not reliable.
IDENTITY_FULL = 90.0
IDENTITY_NONE = 40.0
#: Aligned length at which a match is long enough for the name to mean anything. 34 amino
#: acids is this library's median significant alignment, which is short, so the length
#: term is deliberately harsh below it.
ALIGN_FULL = 60.0
ALIGN_NONE = 20.0

#: The significance threshold every count here is taken at. Never quote an identity
#: without it: the chance band reads identically on residue and already-identified reads.
SIGNIFICANT_EVALUE = 1e-5

#: Design effect on the residue's counts, from results/block_position_bias/. The reads are
#: a contiguous block of the FASTQ and a Wilson interval assumes a random sample.
RESIDUE_DESIGN_EFFECT = 1.09


def ramp(x: float, lo: float, hi: float) -> float:
    if x is None:
        return 0.0
    if hi == lo:
        return 1.0
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def wilson(k: int, n: int, z: float = 1.96, design: float = 1.0) -> dict:
    if not n:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d * design
    return {"point": p, "lo": max(0.0, c - h), "hi": min(1.0, c + h),
            "k": k, "n": n, "design_effect": design}


def classify(text: str, table: dict, index: str | None = None) -> tuple[str, str | None]:
    """First matching rule wins, and the matched keyword is returned so the page can show
    why a row was grouped as it was. A grouping a reader cannot audit is a grouping a
    reader has to take on trust, which is not what this project does.

    Two kinds of keyword, because the table needs both and conflating them is how it goes
    wrong. `any` entries are PREFIXES anchored at the left only, so `Picornavir` catches
    `Picornaviridae`; most virus-family keywords have to work that way. `word` entries are
    anchored at both ends, for host words that are substrings of ordinary English -
    `bat` must not fire on `battery`, and a prefix rule cannot express that.

    A plain substring test put `Bastrovirus` into the human-associated group because
    `Astrovir` occurs inside it: the exact failure mode of a keyword table, caught on the
    first run by reading the matched keyword printed beside the row, which is why it is
    printed.
    """
    # Descending priority, then file order. Python's sort is stable, so a table with no
    # priorities anywhere behaves exactly as it did before this existed.
    for rule in sorted(table["rules"], key=lambda r: -r.get("priority", 0)):
        if rule.get("only_if_index") and rule["only_if_index"] != index:
            continue
        for kw in rule.get("any", ()):
            if re.search(r"(?<![A-Za-z])" + re.escape(kw), text, re.IGNORECASE):
                return rule["group"], kw
        for kw in rule.get("word", ()):
            if re.search(r"(?<![A-Za-z])" + re.escape(kw) + r"(?![A-Za-z])",
                         text, re.IGNORECASE):
                return rule["group"], kw
    return "unknown_host", None


def _hit(kw: str, text: str, word: bool = False) -> bool:
    pat = r"(?<![A-Za-z])" + re.escape(kw) + (r"(?![A-Za-z])" if word else "")
    return re.search(pat, text, re.IGNORECASE) is not None


def sequence_kind(text: str, index: str | None, table: dict) -> tuple[str, str | None]:
    """What KIND of thing the matched record is: a phage, some other virus, or something
    living. This is a different question from the host group, which says what the record
    INFECTS, and a row needs both - a phage record has a bacterial host and is not itself
    a bacterium.

    The split phage-against-bacterial is the hardest of the three and the honest reason is
    biological, not a limit of this table: prophages live inside bacterial genomes, so one
    protein can truthfully be both. At the identities this residue reaches - a few dozen
    amino acids at about half identity - nothing here settles that, so these buckets
    describe the RECORD that was matched and never the read.
    """
    kr = table["kind_rules"]
    if index in kr["viral_indexes"]:
        for kw in kr["phage_name_words"]:
            if _hit(kw, text):
                return "putative_phage", kw
        return "virus_no_phage_signal", None
    # A cellular index can still hold a viral entry, so check before assuming.
    if _hit("virus", text) or _hit("phage", text):
        for kw in kr["phage_name_words"]:
            if _hit(kw, text):
                return "putative_phage", kw
        return "virus_no_phage_signal", None
    for kw in kr["prokaryote_genera"]:
        if _hit(kw, text, word=True):
            return "putative_bacterial", kw
    for kw in kr.get("eukaryote_genera", ()):
        if _hit(kw, text, word=True):
            return "putative_eukaryotic", kw
    return "cellular_kind_not_placed", None


CELLULAR_KINDS = ("putative_bacterial", "putative_eukaryotic", "cellular_kind_not_placed")


def pile_counts(rundir) -> dict:
    """How many of the scored mates have a significant named relative at all, per pile.

    The denominator matters more than the numerator here and had been wrong: shares were
    being taken over the residual PAIRS of the whole cascade block while the numerator
    counted MATES among the sample the per-read report actually searched. Two different
    units and a factor of about 1.8. Counted here from reads.jsonl so both come from one
    file, and the known pile is counted in the same pass because a share of the residue
    without its control says nothing.
    """
    path = rundir / "reads.jsonl"
    if not path.exists():
        return {}
    out: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            d = json.loads(line)
            pile = out.setdefault(d.get("set") or "unknown",
                                  {"scored": 0, "searched": 0, "significant": 0})
            pile["scored"] += 1
            best = None
            for h in (d.get("closest_viral"), d.get("closest_sprot")):
                if h and h.get("evalue") is not None:
                    best = h["evalue"] if best is None else min(best, h["evalue"])
            if best is not None:
                pile["searched"] += 1
                if best < SIGNIFICANT_EVALUE:
                    pile["significant"] += 1
    return out


def kind_summary(rows: list[dict], piles: dict, table: dict) -> dict:
    """The three piles Imo asked for, in READS rather than subjects, with the pile that
    matters most beside them: the mates that point at nothing significant at all.

    Subject counts and read counts answer different questions and the read count is the
    one a reader wants - 201 of 243 subjects being phage is not the same claim as most of
    the residue's named reads being phage. Both ship.

    The known pile is carried in the same structure, because a share of the residue with
    no control beside it is the mistake this repo keeps catching itself making.
    """
    res = piles.get("residue") or {}
    searched = res.get("searched") or 0
    kinds = table["kind_rules"]["kinds"]
    by_kind = {k: {"subjects": 0, "reads": 0, "label": v["label"],
                   "reading": v["reading"]} for k, v in kinds.items()}
    for r in rows:
        b = by_kind[r["sequence_kind"]]
        b["subjects"] += 1
        b["reads"] += r["residue_reads"] or 0
    named = sum(b["reads"] for b in by_kind.values())
    for b in by_kind.values():
        b["share_of_searched"] = (wilson(b["reads"], searched,
                                         design=RESIDUE_DESIGN_EFFECT)
                                  if searched else None)
    known = piles.get("known") or {}
    return {
        "unit": "residue MATES with a significant named relative, at E < "
                f"{SIGNIFICANT_EVALUE:g}",
        "searched": searched,
        "named": named,
        "no_significant_relative": searched - named,
        "share_named": wilson(named, searched, design=RESIDUE_DESIGN_EFFECT)
        if searched else None,
        "by_kind": by_kind,
        "known_pile_control": {
            "searched": known.get("searched"),
            "significant": known.get("significant"),
            "share": wilson(known.get("significant") or 0, known.get("searched") or 0,
                            design=RESIDUE_DESIGN_EFFECT) if known.get("searched") else None,
            "reading": "the same count on reads the cascade already identified. It is the "
                       "scale the residue's figure is read against: these are reads with "
                       "a known relative by construction, and even they do not all have a "
                       "significant one at this threshold",
        },
        "caveat": table["kind_rules"]["_the_hard_split"],
    }


def from_full_search(path: Path) -> tuple[dict, dict, dict]:
    """A whole-residue search, reshaped into what the block-sample report produces.

    WHY RESHAPE RATHER THAN BRANCH. Everything downstream of here - the keyword axes, the
    damping, the composite, the intervals - is the same arithmetic whichever pile the
    subjects came from, and it is arithmetic this repo has already had to debug once. So
    the only thing that changes with the source is the rows and the denominators, and
    they are made to look identical here.

    WHAT THE DENOMINATOR MEANS HERE, AND IT IS NOT THE SAME WORD. In the block sample,
    `searched` meant the mates the per-read report actually put through a search, which
    was a fraction of the pile. In this file every unplaced mate went through the search,
    so the denominator is the whole pile and a share is a share of the residue rather than
    a share of a sample of it. That makes the number larger in meaning and smaller in
    value than the block figure, which is why the two are reported as different
    instruments and never merged.

    The control's per-subject counts bound ALREADY SEEN from below: it is a sample of the
    placed pile, so a subject missing from it may still be placed material this sample did
    not happen to draw. The row says which pile it was compared against.
    """
    d = json.loads(path.read_text())
    res = d["piles"]["residue"]
    ctl = d["piles"]["control"]
    ctl_reads = ctl.get("control_subject_reads") or {}
    rows = []
    for r in res.get("subjects") or []:
        known = int(ctl_reads.get(r["subject"], 0))
        rows.append({**r, "known_reads": known,
                     "reading": ("present and partly identified: the already-placed "
                                 "sample points at this subject too"
                                 if known else
                                 "pointed at only by unplaced reads in this sample of "
                                 "the placed pile")})
    S = {"subjects": rows,
         "significance_evalue": d["operating_point"]["significance_evalue"],
         "subjects_total": len(rows),
         "known_pile_searched_reads": ctl["mates"]}
    piles = {"residue": {"scored": res["mates"], "searched": res["mates"],
                         "significant": res["mates_with_significant_hit"],
                         "with_any_reported_hit": res["searched_mates_with_any_hit"]},
             "known": {"scored": ctl["mates"], "searched": ctl["mates"],
                       "significant": ctl["mates_with_significant_hit"],
                       "with_any_reported_hit": ctl["searched_mates_with_any_hit"]}}
    meta = {"source": path.name,
            "residue_pairs": res["pairs"], "residue_mates": res["mates"],
            "control_pairs": ctl["pairs"],
            "library_pairs": d["library_pairs"],
            "unplaced_fraction": d["unplaced_fraction"],
            "operating_point": d["operating_point"],
            "bias": d["bias"],
            "wall_seconds": d.get("wall_seconds"),
            "script_sha256": d.get("script_sha256"),
            "git_sha": d.get("git_sha"),
            "scope": ("every unplaced pair in the library, not a block sample. The pile "
                      "comes from a nucleotide-only label table, so it is a little more "
                      "generous than the cascade's residue and the named-relative rate "
                      "is if anything a small overestimate")}
    return S, piles, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--table", type=Path, default=TABLE)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--full-residue", type=Path, default=None,
                    help="a full_residue_search_s*.json to rank instead of the block "
                         "sample's strains.json. Default: use one if the run dir has it")
    ap.add_argument("--block-sample", action="store_true",
                    help="rank the block sample even when a full-residue search is there")
    args = ap.parse_args()

    full = args.full_residue
    if full is None and not args.block_sample:
        found = sorted(args.rundir.glob("full_residue_search_s*.json"))
        full = found[0] if found else None
    scope = None
    if full is not None:
        S, full_piles, scope = from_full_search(full)
    else:
        strains = args.rundir / "strains.json"
        if not strains.exists():
            raise SystemExit(f"{strains} absent: the per-read report did not run here")
        S = json.loads(strains.read_text())
        full_piles = None
    names = (json.loads((args.rundir / "subject_names.json").read_text())
             if (args.rundir / "subject_names.json").exists() else {}).get("names", {})
    table = json.loads(args.table.read_text())
    groups = table["groups"]

    tt = (json.loads((args.rundir / "tier_table.json").read_text())
          if (args.rundir / "tier_table.json").exists() else {})
    # The share's denominator is the residue mates the per-read report actually SEARCHED,
    # counted from the same file as the numerator. It used to be the residual pairs of the
    # whole cascade block, which is both a different unit and about 1.8x larger, so every
    # share was understated by that factor. Numerator and denominator now come from one
    # pass over one file, which is the only way this stays true when either changes.
    piles = full_piles if full_piles is not None else pile_counts(args.rundir)
    n_denom = (piles.get("residue") or {}).get("searched") or 0
    denom_note = (("every unplaced mate in the library, all of which went through this "
                   "search, so a share here is a share of the residue rather than of a "
                   "sample of it")
                  if full_piles is not None else
                  ("residue mates searched by the per-read report (reads.jsonl), not the "
                   "residual pairs of the cascade block: the numerator counts mates in "
                   "that same sample"))

    rows, unplaced = [], 0
    for r in S.get("subjects", []):
        subj = r["subject"]
        nm = names.get(subj) or {}
        text = " ".join(x for x in [nm.get("organism"), nm.get("name"), subj] if x)
        kind, kindkw = sequence_kind(text, r.get("index"), table)
        group, kw = classify(text, table, r.get("index"))
        # A cellular record has no host: it IS the organism. Letting the host rules speak
        # here is a category error, and a live one - it had turned a Bacillus subtilis
        # transporter into a bacteriophage and an Arabidopsis kinase into a plant virus.
        # The row keeps the unknown weight rather than a guess, and says why on its face.
        host_na = kind in CELLULAR_KINDS
        if host_na:
            group, kw = "unknown_host", None
        if group == "unknown_host":
            unplaced += 1
        ident = r.get("median_identity_percent")
        aln = r.get("median_aligned_aa")

        # Confidence that the matched record's host group transfers to these reads.
        conf = ramp(ident, IDENTITY_NONE, IDENTITY_FULL) * ramp(aln, ALIGN_NONE, ALIGN_FULL)
        host_term = (groups[group]["concern"] * conf
                     + groups["unknown_host"]["concern"] * (1 - conf))
        # Distance: 0 at an identical match, 1 at the chance floor. Clamped rather than
        # extrapolated, because below the chance band the identity is not a distance.
        distance = 1.0 - ramp(ident, IDENTITY_NONE, 100.0)
        # Weight: log of the reads pointing here, normalised by the largest subject in this
        # run, so it is an ordering within the run and never an abundance claim.
        reads = r.get("residue_reads") or 0
        alone = 1.0 if (S.get("known_pile_searched_reads")
                        and (r.get("known_reads") or 0) == 0) else 0.0

        rows.append({
            "subject": subj,
            "organism": nm.get("organism"),
            "protein": nm.get("name"),
            "index": r.get("index"),
            "host_group": group,
            "host_group_matched_keyword": kw,
            "host_group_label": groups[group].get("label", group.replace("_", " ")),
            "host_group_reading": (
                "not applicable: the matched record is not a virus, so it has no host. "
                "Weighted as unknown rather than as safe."
                if host_na else groups[group]["reading"]),
            "host_group_not_applicable": host_na,
            "sequence_kind": kind,
            "sequence_kind_label": table["kind_rules"]["kinds"][kind]["label"],
            "sequence_kind_reading": table["kind_rules"]["kinds"][kind]["reading"],
            "sequence_kind_matched_keyword": kindkw,
            "host_confidence": round(conf, 4),
            "host_term": round(host_term, 4),
            "distance": round(distance, 4),
            "median_identity_percent": ident,
            "identity_range": r.get("identity_range"),
            "median_aligned_aa": aln,
            "median_evalue": r.get("median_evalue"),
            "worst_evalue": r.get("worst_evalue"),
            "residue_reads": reads,
            "known_reads": r.get("known_reads"),
            "only_in_residue": bool(alone),
            "share_of_residue": wilson(reads, int(n_denom), design=RESIDUE_DESIGN_EFFECT)
            if n_denom else None,
        })

    maxreads = max((r["residue_reads"] for r in rows), default=1) or 1
    for r in rows:
        r["weight_term"] = round(math.log1p(r["residue_reads"]) / math.log1p(maxreads), 4)
        # The composite, written out so the page can print it verbatim. Multiplicative,
        # because a row needs all four to deserve a person's morning: plenty of reads, a
        # real distance, a host group that survives damping, and nobody already looking
        # at it. The +0.25 floors keep a single weak component from zeroing a row out.
        r["second_look"] = round(
            (0.25 + 0.75 * r["weight_term"]) *
            (0.25 + 0.75 * r["distance"]) *
            r["host_term"] *
            (1.0 if r["only_in_residue"] else 0.6), 5)
    rows.sort(key=lambda r: -r["second_look"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    out = {
        "unit": "a named relative that residue reads point at, at E < "
                f"{S.get('significance_evalue')}",
        "not_a_risk_assessment":
            "this is a triage order for a person's attention, built from four measured "
            "quantities. It is not a probability of harm, not a clinical finding, and not "
            "a statement that anything here infects anyone. Nothing in it scores novelty: "
            "an index can say whether something is already known and cannot say whether "
            "it is new",
        "formula": "(0.25 + 0.75*weight) * (0.25 + 0.75*distance) * host_term * "
                   "(1.0 if only_in_residue else 0.6)",
        "formula_terms": {
            "weight": "log1p(residue reads) over log1p(the largest subject in this run); "
                      "an ordering within the run, never an abundance claim",
            "distance": "1 - identity ramped from 40% to 100%. 40% is where this repo's "
                        "own BLAST ground truth puts the chance band for a short query, "
                        "so below it an identity is not a distance",
            "host_term": "the matched record's group weight, damped towards unknown by a "
                         "confidence that falls with distance and with short alignments. "
                         "Distance raises interest and lowers host confidence at once, "
                         "and this is where that tension is paid for rather than hidden",
            "only_in_residue": "0.6 if the already-identified pile points at the same "
                               "subject, because then it is something we have, seen again",
        },
        "host_confidence_ramp": {"identity_full": IDENTITY_FULL,
                                 "identity_none": IDENTITY_NONE,
                                 "aligned_aa_full": ALIGN_FULL,
                                 "aligned_aa_none": ALIGN_NONE},
        "host_table": {"groups": groups,
                       "source": str(args.table.relative_to(ROOT)),
                       "subjects_no_rule_placed": unplaced,
                       "subjects_total": len(rows),
                       "coverage_note": "an unplaced name is grouped unknown_host and "
                                        "weighted in the middle, never at zero"},
        # The whole split in one block, with its control beside it. This is what the page
        # draws as its summary, so it is computed here and never recomputed in the page.
        "kind_summary": kind_summary(rows, piles, table),
        "piles": piles,
        "denominator_note": denom_note,
        # Which pile this ranking was built from, so the page can say it on its face
        # rather than a reader assuming the wider one.
        "scope": scope or {"source": "strains.json",
                           "scope": "the per-read report's block sample of the residue"},
        "interval_note": f"shares of the residue carry a Wilson interval widened by "
                         f"{RESIDUE_DESIGN_EFFECT}x, the measured design effect of "
                         f"reading a contiguous block of the FASTQ "
                         f"(results/block_position_bias/)",
        "blind_to": [
            "anything below about 1 read pair in a million, which is most of what a "
            "surveillance tool would want to catch early",
            "anything whose nearest relative is further than roughly 10-20% divergence, "
            "because then there is no named relative to rank at all",
            "any subject the keyword table could not place, which is counted above and "
            "ranked as unknown rather than as safe",
        ],
        # The headline of this panel, and it is a negative. If no relative is close
        # enough for its host group to be believed, then the host axis carries nothing on
        # this library and the ranking is distance and weight wearing a colour. Printing
        # the distribution rather than the best case is the only way a reader can tell.
        "host_confidence_summary": {
            "subjects": len(rows),
            "max": round(max((r["host_confidence"] for r in rows), default=0.0), 4),
            "median": round(sorted(r["host_confidence"] for r in rows)[len(rows) // 2], 4)
            if rows else None,
            "at_or_above_0.5": sum(1 for r in rows if r["host_confidence"] >= 0.5),
            "at_or_above_0.2": sum(1 for r in rows if r["host_confidence"] >= 0.2),
            "best_identity_percent": max((r["median_identity_percent"] or 0)
                                         for r in rows) if rows else None,
            "best_aligned_aa": max((r["median_aligned_aa"] or 0) for r in rows)
            if rows else None,
            "reading": "how much of the matched record's host group transfers to the "
                       "reads. Built from identity and aligned length. If the maximum is "
                       "low, no relative in this residue is close enough for its host to "
                       "mean anything, and the host colours on this page are labels on "
                       "the database rather than claims about the sample",
        },
        "subjects": rows,
        **provenance(__file__),
    }
    dest = args.out or (args.rundir / "ranking.json")
    dest.write_text(json.dumps(out, indent=1))
    print(f"{len(rows)} subjects ranked, {unplaced} unplaced by the table -> {dest}")
    for r in rows[:10]:
        print(f"  {r['rank']:3d} {r['second_look']:.3f} {r['host_group']:22s} "
              f"{(r['organism'] or r['subject'])[:44]:44s} "
              f"{r['residue_reads']:4d} reads  {r['median_identity_percent']}%")


if __name__ == "__main__":
    main()
