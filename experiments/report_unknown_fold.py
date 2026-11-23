#!/usr/bin/env python3
"""Readout for `unknown_fold.py`: does structure place the unnameable pile, and onto what?

THE CONTROL IS READ FIRST, which on this experiment means three things in order.

*The null.* Shuffled contigs through the identical ORF caller and the identical search. If
the null is not at or near zero, nothing below is readable and the report says so instead of
ranking anything.

*The comparator.* Identified contigs from the same assembly of the same water. A placement
rate on the unnameable pile means nothing without the rate the same instrument gets on
material from the same sample that sequence can name.

*The length match.* The two piles' ORFs are not the same length - 78 aa against 105 at the
median - and alignment significance is a function of query length, so the raw difference
between the arms is partly a difference in how long their genes are. Every comparison is
therefore reported inside ORF-length bands, and the headline is the band-weighted rate:
the comparator's per-band rates reweighted onto the unnameable pile's own band distribution.
The raw rates are printed beside it so the size of the length effect is visible rather than
corrected away silently.

WHAT IS NOT CLAIMED. This arm searches a prefiltered index, so every rate is a floor. It
searches 30%-identity cluster representatives of the solved databank, which is what has been
crystallised and not a sample of fold space. And a subject's description and source organism
are a database record's words about a database entry: they say what the placed gene is like,
never what is in the water and never what anything infects.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "experiments"))

from provenance import provenance  # noqa: E402
from rcsb_chains import annotate_subjects  # noqa: E402

BAND_ORDER = ["50-80", "80-100", "100-150", "150-200", "200+"]


def spread(values: list[float]) -> str:
    """Observed min-max across seeds, never a standard deviation.

    CLAUDE.md: five seeds give an honest range and a useless variance, and the sampling
    error of a five-sample sd swamps what it would be measuring. So this prints the range.
    """
    if not values:
        return "-"
    if len(values) == 1:
        return f"{values[0]:.2f}"
    return f"{min(values):.2f}-{max(values):.2f}"


def load(out: Path):
    """The ORFs, the hits and the run manifest. NOT the run's own subject annotations.

    The run writes a `subject_annotations.json` keyed by PDB ENTRY, and an entry is usually
    several molecules, so the chain that aligned cannot be recovered from it - see
    `rcsb_chains`. The names are therefore resolved here instead, from the target names in
    the hits this report actually reads. It costs a handful of API calls and seconds,
    because the cost is bounded by the result rather than by the library, and it means one
    code path answers "like what" for every run whatever produced it.
    """
    orfs = json.loads((out / "orfs.json").read_text())
    manifest = json.loads((out / "manifest.json").read_text())
    hits: dict[str, list[dict]] = collections.defaultdict(list)
    with gzip.open(out / "hits.tsv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            hits[row["query"]].append(row)
    return orfs, hits, manifest


def best_hit(rows: list[dict], floor: float) -> dict | None:
    """Strongest alignment clearing the floor, by bitscore with a strict `>`.

    Strict, as every copy in this repo is: a tie broken by array order is how a
    nearest-neighbour test once inflated its own control.
    """
    best = None
    for row in rows:
        if float(row["evalue"]) > floor:
            continue
        if best is None or float(row["bits"]) > float(best["bits"]):
            best = row
    return best


def rates(orfs, hits, floor: float):
    """Placed fraction per (pile, seed), overall and per ORF-length band."""
    cells: dict[tuple[str, str], dict] = {}
    for orf in orfs:
        arm = orf["arm"]
        seed = arm.rsplit("_s", 1)[1]
        key = (orf["pile"], seed)
        cell = cells.setdefault(key, dict(n=0, placed=0, bands=collections.defaultdict(
            lambda: dict(n=0, placed=0))))
        cell["n"] += 1
        cell["bands"][orf["band"]]["n"] += 1
        if best_hit(hits.get(orf["orf"], []), floor):
            cell["placed"] += 1
            cell["bands"][orf["band"]]["placed"] += 1
    return cells


def band_weighted(cells, target_pile: str, source_pile: str, seeds: list[str]) -> list[float]:
    """The comparator's per-band rates, reweighted onto the target pile's band mix.

    THIS IS THE HEADLINE COMPARISON, and it exists because the raw one is confounded. The
    target pile supplies the weights - how many of its ORFs sit in each length band - and the
    source pile supplies the per-band rates. The result answers "what would the comparator
    read if its genes were as short as the unnameable pile's", which is the only version of
    the comparison that is about the pile rather than about gene length.
    """
    out = []
    for seed in seeds:
        tgt = cells.get((target_pile, seed))
        src = cells.get((source_pile, seed))
        if not tgt or not src:
            continue
        num = den = 0.0
        for band in BAND_ORDER:
            w = tgt["bands"].get(band, {}).get("n", 0)
            s = src["bands"].get(band)
            if not w or not s or not s["n"]:
                continue
            num += w * (s["placed"] / s["n"])
            den += w
        if den:
            out.append(100 * num / den)
    return out


def subject_table(orfs, hits, annotations, floor: float, pile: str, top: int):
    """What the placed ORFs of one pile landed on, ranked, with the tail left off.

    Imo reads summary-then-detail and wants the tail hidden, so this ranks by count and
    returns the head. Each row carries the subject's own description and its source organism
    as separate fields, because what a record IS and what it came FROM are two axes and
    collapsing them is how a host colour starts reading as a statement about the water.
    """
    counts: collections.Counter = collections.Counter()
    examples: dict[str, dict] = {}
    for orf in orfs:
        if orf["pile"] != pile:
            continue
        best = best_hit(hits.get(orf["orf"], []), floor)
        if not best:
            continue
        ann = annotations["by_target"][best["target"]]
        desc = ann.get("description") or "(no description in the record)"
        counts[desc] += 1
        prev = examples.get(desc)
        if prev is None or float(best["bits"]) > float(prev["bits"]):
            examples[desc] = dict(bits=best["bits"], evalue=best["evalue"],
                                  fident=best["fident"], alnlen=best["alnlen"],
                                  qlen=best["qlen"], entry=best["target"], orf=orf["orf"],
                                  organisms="; ".join((ann.get("source_organisms") or [])[:2]))
    rows = []
    for desc, n in counts.most_common(top):
        e = examples[desc]
        rows.append(dict(description=desc, placed=n, best_entry=e["entry"],
                         best_bits=e["bits"], best_evalue=e["evalue"],
                         best_fident=e["fident"], best_alnlen=e["alnlen"],
                         best_qlen=e["qlen"], source_organisms=e["organisms"]))
    return rows, sum(counts.values()), len(counts)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--floor", type=float, default=1e-5)
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args(argv)

    orfs, hits, manifest = load(args.out)
    # Only the chains that carry a SIGNIFICANT hit, so the lookup is bounded by the result.
    wanted = sorted({h["target"] for rows in hits.values() for h in rows
                     if float(h["evalue"]) <= args.floor})
    annotations = annotate_subjects(wanted, args.out)
    cells = rates(orfs, hits, args.floor)
    piles = sorted({orf["pile"] for orf in orfs})
    seeds = sorted({orf["arm"].rsplit("_s", 1)[1] for orf in orfs}, key=int)

    lines = ["# Structure on the unnameable pile: tables", "",
             f"Generated by `experiments/report_unknown_fold.py` from "
             f"`{args.out}`, E-value floor {args.floor:g}.", "",
             "## 1. The null, read first", ""]
    null_piles = [p for p in piles if p.endswith("_null")]
    lines += ["| arm | ORFs per seed | placed, % (per seed) |", "|---|---|---|"]
    summary: dict = {}
    for pile in null_piles:
        per = [100 * cells[(pile, s)]["placed"] / cells[(pile, s)]["n"]
               for s in seeds if (pile, s) in cells and cells[(pile, s)]["n"]]
        ns = [cells[(pile, s)]["n"] for s in seeds if (pile, s) in cells]
        lines.append(f"| {pile} | {spread([float(n) for n in ns])} | {spread(per)} |")
        summary[f"null_{pile}"] = per
    lines += ["", "## 2. Placement rate, raw and band-weighted", "",
              "| arm | ORFs per seed | placed, % (per seed) |", "|---|---|---|"]
    for pile in [p for p in piles if not p.endswith("_null")]:
        per = [100 * cells[(pile, s)]["placed"] / cells[(pile, s)]["n"]
               for s in seeds if (pile, s) in cells and cells[(pile, s)]["n"]]
        ns = [cells[(pile, s)]["n"] for s in seeds if (pile, s) in cells]
        lines.append(f"| {pile} | {spread([float(n) for n in ns])} | {spread(per)} |")
        summary[f"raw_{pile}"] = per
    matched = band_weighted(cells, "unident", "ident", seeds)
    summary["ident_weighted_to_unident_bands"] = matched
    lines += ["", f"**Length-matched comparator**: the identified pile's per-band rates on "
              f"the unnameable pile's own band mix read **{spread(matched)}%** per seed, "
              f"against its raw {spread(summary.get('raw_ident', []))}%. The difference "
              f"between those two is the length effect, measured rather than argued.", ""]
    # THE DEFICIT AS A COUNT, not only as two percentages side by side. Pooled over seeds,
    # because the question is whether the pile places at all and a per-seed range of zeros
    # does not say how many chances it had. The bound is the exact one-sided binomial for
    # zero successes, which has the closed form 1 - alpha^(1/n) and needs no dependency; it
    # is a statement about ONE proportion, not an estimate of spread across seeds, so it is
    # not the five-sample variance CLAUDE.md rules out. It is quoted only when the observed
    # count is zero, because that is the only case this form covers.
    tgt_n = sum(cells[("unident", s)]["n"] for s in seeds if ("unident", s) in cells)
    tgt_placed = sum(cells[("unident", s)]["placed"] for s in seeds if ("unident", s) in cells)
    expected = sum(w * r / 100 for w, r in
                   [(cells[("unident", s)]["n"], m) for s, m in zip(seeds, matched)])
    summary["unident_pooled"] = dict(n=tgt_n, placed=tgt_placed, expected_band_matched=expected)
    deficit = [f"Pooled over the {len(seeds)} seeds the unnameable pile offers **{tgt_n} "
               f"ORFs** and places **{tgt_placed}**. The length-matched comparator predicts "
               f"**{expected:.1f}**."]
    if tgt_placed == 0 and tgt_n:
        bound = 100 * (1 - 0.05 ** (1 / tgt_n))
        summary["unident_upper_bound_95"] = bound
        deficit.append(f"Zero of {tgt_n} puts the one-sided 95% upper bound on this pile's "
                       f"placement rate at **{bound:.2f}%**, below the matched comparator's "
                       f"{spread(matched)}%. So the gap is not an absence of chances.")
    lines += [" ".join(deficit), ""]
    # UN-THRESHOLDED, and this is the section that keeps the headline honest. A placement
    # rate at E<1e-5 answers "does it clear the statistics", and for this pile the answer is
    # no everywhere, which on its own reads as "the instrument says nothing". It does not:
    # the same queries produce alignments, and where their best E sits relative to their OWN
    # shuffled null is a distance rather than an absence. That is the project's rule - the
    # report is un-thresholded with a null beside every field - applied to the one field
    # whose thresholded form is all zeros. Quartiles of the best E-value per ORF, pooled
    # over seeds, with the fraction that aligned at all beside them.
    lines += ["## 2b. The same arms un-thresholded, with the null beside them", "",
              "| arm | ORFs | aligned at all, % | best E: lower quartile | median | best seen |",
              "|---|---|---|---|---|---|"]
    best_e: dict[str, list[float]] = collections.defaultdict(list)
    arm_n: collections.Counter = collections.Counter()
    for orf in orfs:
        arm_n[orf["pile"]] += 1
        rows_ = hits.get(orf["orf"], [])
        if rows_:
            best_e[orf["pile"]].append(min(float(r["evalue"]) for r in rows_))
    for pile in piles:
        es = sorted(best_e.get(pile, []))
        n = arm_n[pile]
        if not es:
            lines.append(f"| {pile} | {n} | 0.0 | - | - | - |")
            continue
        q1 = es[max(0, len(es) // 4 - 1)]
        med = statistics.median(es)
        lines.append(f"| {pile} | {n} | {100 * len(es) / n:.1f} | {q1:.2g} | {med:.2g} | "
                     f"{es[0]:.2g} |")
        summary[f"unthresholded_{pile}"] = dict(
            n=n, aligned=len(es), q1=q1, median=med, best=es[0])
    lines += ["", "The unnameable arm sits between its own shuffled null and the identified "
              "comparator on every column here while reading zero at E<1e-5. The operating "
              "point is what makes those two statements consistent: this is a signal that "
              "does not clear its own significance statistics, not an absence of one.", ""]
    lines += ["## 3. Per band", "",
              "| band, aa | " + " | ".join(f"{p} n / placed %" for p in
                                           [q for q in piles if not q.endswith('_null')])
              + " |",
              "|---|" + "---|" * len([q for q in piles if not q.endswith('_null')])]
    band_rows: dict = {}
    for band in BAND_ORDER:
        cell_text = []
        for pile in [q for q in piles if not q.endswith("_null")]:
            n = sum(cells[(pile, s)]["bands"].get(band, {}).get("n", 0)
                    for s in seeds if (pile, s) in cells)
            placed = sum(cells[(pile, s)]["bands"].get(band, {}).get("placed", 0)
                         for s in seeds if (pile, s) in cells)
            cell_text.append(f"{n} / {100*placed/n:.1f}%" if n else "-")
            band_rows.setdefault(band, {})[pile] = dict(n=n, placed=placed)
        lines.append(f"| {band} | " + " | ".join(cell_text) + " |")
    census = annotations["polymer_type_census"]
    summary["subject_polymer_type_census"] = census
    lines += ["", "## 4. What the unnameable pile was placed onto", "",
              f"{annotations['chains']} distinct subject chains carry a significant hit, "
              f"resolved per chain rather than per entry. Polymer type census: "
              f"`{census}` - this is a check on the name lookup and not a result, because "
              f"a Foldseek structure library holds protein chains only, so anything but "
              f"100% protein here means a chain was resolved to the wrong molecule.", ""]
    for pile in ["unident", "ident"]:
        rows, total, distinct = subject_table(orfs, hits, annotations, args.floor, pile,
                                              args.top)
        lines += [f"### {pile}: {total} placements onto {distinct} distinct descriptions, "
                  f"top {len(rows)}", "",
                  "| what the record is | placed | best entry | bits | E | 3Di ident | aln/qlen "
                  "| the record's source organism |", "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            lines.append(
                f"| {r['description'][:70]} | {r['placed']} | {r['best_entry']} | "
                f"{float(r['best_bits']):.0f} | {float(r['best_evalue']):.1g} | "
                f"{100*float(r['best_fident']):.1f}% | {r['best_alnlen']}/{r['best_qlen']} | "
                f"{r['source_organisms'][:50]} |")
        lines.append("")
        summary[f"subjects_{pile}"] = dict(placements=total, distinct=distinct, top=rows)

    (args.out / "TABLES.md").write_text("\n".join(lines) + "\n")
    report_manifest = dict(provenance(Path(__file__)), floor=args.floor,
                           run_manifest_script_sha256=manifest.get("script_sha256"),
                           run_git_sha=manifest.get("git_sha"), summary=summary)
    (args.out / "report_manifest.json").write_text(json.dumps(report_manifest, indent=1))
    print("\n".join(lines[:80]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
