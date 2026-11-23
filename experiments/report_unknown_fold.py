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
    orfs = json.loads((out / "orfs.json").read_text())
    annotations = json.loads((out / "subject_annotations.json").read_text())
    manifest = json.loads((out / "manifest.json").read_text())
    hits: dict[str, list[dict]] = collections.defaultdict(list)
    with gzip.open(out / "hits.tsv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            hits[row["query"]].append(row)
    return orfs, hits, annotations, manifest


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
        entry = best["target"].split("_")[0].split(".")[0].upper()
        ann = annotations.get(entry, {})
        desc = (ann.get("descriptions") or ["(no description in the record)"])[0]
        counts[desc] += 1
        prev = examples.get(desc)
        if prev is None or float(best["bits"]) > float(prev["bits"]):
            examples[desc] = dict(bits=best["bits"], evalue=best["evalue"],
                                  fident=best["fident"], alnlen=best["alnlen"],
                                  qlen=best["qlen"], entry=entry, orf=orf["orf"],
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

    orfs, hits, annotations, manifest = load(args.out)
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
    lines += ["", "## 4. What the unnameable pile was placed onto", ""]
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
