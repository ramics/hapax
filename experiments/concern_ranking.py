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

WHY IT IS NOT A PIPELINE STAGE. It measures nothing. It reads a finished run directory and
reorders it, so it lives beside `subject_names.py` as a presentation step. Keeping it out
of `hapax.py` is the honest signal that no number in the run changes because of it.
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


def classify(text: str, table: dict) -> tuple[str, str | None]:
    """First matching rule wins, and the matched keyword is returned so the page can show
    why a row was grouped as it was. A grouping a reader cannot audit is a grouping a
    reader has to take on trust, which is not what this project does.

    Keywords match only at a word boundary. A plain substring test put `Bastrovirus` into
    the human-associated group because `Astrovir` occurs inside it, which is the exact
    failure mode of a keyword table and was caught on the first run by reading the matched
    keyword printed beside the row - which is why it is printed.
    """
    for rule in table["rules"]:
        for kw in rule["any"]:
            if re.search(r"(?<![A-Za-z])" + re.escape(kw), text, re.IGNORECASE):
                return rule["group"], kw
    return "unknown_host", None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--table", type=Path, default=TABLE)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    strains = args.rundir / "strains.json"
    if not strains.exists():
        raise SystemExit(f"{strains} absent: the per-read report did not run here")
    S = json.loads(strains.read_text())
    names = (json.loads((args.rundir / "subject_names.json").read_text())
             if (args.rundir / "subject_names.json").exists() else {}).get("names", {})
    table = json.loads(args.table.read_text())
    groups = table["groups"]

    tt = (json.loads((args.rundir / "tier_table.json").read_text())
          if (args.rundir / "tier_table.json").exists() else {})
    residue_pairs = None
    br = (tt.get("breakdown") or {}).get("residual")
    if br:
        residue_pairs = br.get("k")
    n_denom = residue_pairs or (tt.get("pairs") or 0)

    rows, unplaced = [], 0
    for r in S.get("subjects", []):
        subj = r["subject"]
        nm = names.get(subj) or {}
        text = " ".join(x for x in [nm.get("organism"), nm.get("name"), subj] if x)
        group, kw = classify(text, table)
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
            "host_group_reading": groups[group]["reading"],
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
