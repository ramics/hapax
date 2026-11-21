"""Aggregate the BLAST ground truth across seeds, per class, with spread.

WHY IT NEEDS AGGREGATING AT ALL. This is the only check in the project that asks the real
public databases what the cascade's labels are worth, and it was reported from one seed
while the repo rule is five. One seed of 60 sequences per class is 60 sequences; five is
300, and more to the point it shows whether the per-class separation is stable or whether a
class's number moves by twenty points between blocks.

WHAT TO READ AND WHAT TO IGNORE. `any_hit` is nearly useless and is printed anyway, because
its uselessness is the result: almost everything hits something, so a hit-or-no-hit rule
cannot work and nobody should build one. **The strong-hit fraction is the column that
matters** - high identity over most of the query - because identity alone counts a third of
a read at 95% as a 95% match, which it is not.

THE ERROR COLUMN IS PART OF THE TABLE, AND ONE KEY DECIDES WHETHER A SEED CAN BE USED AT
ALL. The public service occasionally returns malformed JSON, and the first version of this
check scored those batches as MISSES - a dropped batch and a batch of non-matches are
indistinguishable in a hit rate, which is how a 50% any-hit figure appeared for a class that
really sits at 100%. The fix excludes those queries from the denominator and records how
many per class in `dropped_failed_batch`. A seed file written before that fix has no such
key, and its losses cannot be repaired afterwards: the rows carry `hit: false` for a query
nobody ever got an answer for. So the presence of the key is the test, a seed without it is
excluded and named, and the aggregate never silently mixes the two conventions. This is the
same rule the protein ladder needed for a seed truncated by a full disk.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CLASSES = ("rrna", "public_nucleotide", "selfasm_identified", "residue_assembled",
           "residue_missing", "tier2_contigs")

#: What each class is, in one line, because a table of six slugs is unreadable without it
#: and the distinction between the two residue classes is the architecture's own.
GLOSS = {
    "rrna": "removed by SILVA",
    "public_nucleotide": "removed by Kraken2, not by SILVA",
    "selfasm_identified": "removed only by this sample's own identified contigs",
    "residue_assembled": "survived, maps to an unidentifiable contig of this sample (tier 2)",
    "residue_missing": "survived, maps to no contig at all (tier 3)",
    "tier2_contigs": "unidentifiable contigs >=500 nt, >=2x coverage",
}


def spread(xs) -> str:
    if not xs:
        return "n/a"
    if len(xs) == 1:
        return f"{xs[0] * 100:.1f} (1 seed)"
    return f"{st.mean(xs) * 100:.1f} ({min(xs) * 100:.1f}-{max(xs) * 100:.1f})"


def med(xs) -> str:
    return f"{st.mean(xs):.1f}" if xs else "n/a"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=ROOT / "results/blast_check")
    ap.add_argument("--out", type=Path, default=ROOT / "results/blast_check/TABLES.md")
    args = ap.parse_args()

    paths = sorted(glob.glob(str(args.dir / "blast_check_s*.json")))
    if not paths:
        raise SystemExit(f"nothing in {args.dir}")
    loaded = [(Path(p).name, json.loads(Path(p).read_text())) for p in paths]
    seeds, prefix = [], []
    for name, s in loaded:
        per = next((s["blast"][pr]["per_class"] for pr in s["blast"] if pr != "errors"),
                   {})
        if any("dropped_failed_batch" not in v for v in per.values()):
            prefix.append((name, s["seed"]))
            continue
        seeds.append(s)
    if not seeds:
        raise SystemExit(f"every seed file predates the dropped-batch fix: {prefix}")

    programs = sorted({p for s in seeds for p in s["blast"] if p != "errors"})
    lines = ["# BLAST ground truth against the public databases", "",
             f"{len(seeds)} seed(s): {', '.join(str(s['seed']) for s in seeds)}. "
             f"60 sequences per class per seed, labelled by the same references and "
             f"thresholds as the cascade.", ""]

    if prefix:
        lines += ["**" + str(len(prefix)) + " seed(s) excluded as pre-fix** ("
                  + ", ".join(f"{n} (seed {sd})" for n, sd in prefix)
                  + "): written before failed batches were excluded from the denominator, "
                    "so their misses and their lost queries are the same number and cannot "
                    "be separated now. Queued to re-run rather than averaged in.", ""]
    errors = {s["seed"]: len(s["blast"].get("errors", [])) for s in seeds}
    if any(errors.values()):
        lines += [f"**Batches lost to malformed replies and correctly dropped**, per seed: "
                  + ", ".join(f"seed {k}: {v}" for k, v in errors.items() if v)
                  + ". Those queries are dropped rather than scored as misses - scoring "
                    "them as misses is how an earlier version of this table reported 50% "
                    "any-hit for a class that sits at 100%.", ""]

    for program in programs:
        db = next((s["blast"][program].get("database") for s in seeds
                   if program in s["blast"]), "?")
        lines += [f"## {program} against {db}", "",
                  "| class | what it is | any hit | strong hit | median identity | "
                  "median query coverage |", "|---|---|---|---|---|---|"]
        for cls in CLASSES:
            per = [s["blast"][program]["per_class"].get(cls) for s in seeds
                   if program in s["blast"]]
            per = [p for p in per if p]
            if not per:
                continue
            lines.append(
                f"| `{cls}` | {GLOSS[cls]} | {spread([p['any_hit_frac'] for p in per])} "
                f"| **{spread([p['strong_hit_frac'] for p in per])}** "
                f"| {med([p['pident']['median'] for p in per if p['pident']['n']])} "
                f"| {med([p['query_cov']['median'] for p in per if p['query_cov']['n']])} |")
        lines.append("")

    sizes = {cls: [s["class_fracs_in_block"][cls] for s in seeds
                   if cls in s.get("class_fracs_in_block", {})] for cls in CLASSES}
    lines += ["## How big each class is in the library", "",
              "The sampled 60 per class says nothing about abundance, so the proportions "
              "come from the labelled block itself.", "",
              "| class | share of pairs |", "|---|---|"]
    for cls, xs in sizes.items():
        if xs:
            lines.append(f"| `{cls}` | {spread(xs)} |")

    lines += ["", "## What the table says", "",
              "- **Any hit is not a test.** Read the strong-hit column instead; the two "
              "disagree by tens of points on exactly the classes the architecture cares "
              "about.",
              "- The two residue classes separate from the removed ones continuously, in "
              "identity and coverage together, which is the same shape the contig anatomy "
              "and the coverage-shape work both arrived at independently.",
              "- 60 per class per seed gives percentage-level resolution and no better. "
              "Every figure here is a proportion with an interval in the per-seed files; "
              "the range across seeds is in brackets above."]
    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
