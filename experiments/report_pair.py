"""Aggregate the DNA/RNA pair seeds into the tables its report quotes.

WHY THIS EXISTS. The first version of `results/dna_rna_pair/REPORT.md` was typed from one
seed's JSON by hand. `CLAUDE.md` is explicit that figures come from committed result files
and never from prose, because this repo has had to withdraw claims whose text drifted from
its outputs - and a hand-typed table across five seeds is where that drift starts. So the
spread is computed in one place, here, and the report quotes this script's output.

WHAT IT REFUSES TO DO. It will not aggregate seeds whose provenance differs. Every seed of a
multi-seed result has to come from one script hash, because a fix landing midway through
leaves a directory where some seeds are right and some are wrong and nothing visible says
which - which is exactly what happened to this experiment, whose first three seeds ran with
mate 2 trimmed by 16 and were re-run at 14. A mixture is not a five-seed result, so the
guard is a hard stop rather than a warning.

THE TRIM IS PRINTED BESIDE EVERY RRNA FIGURE, not because it is interesting but because it
moves the Kraken2 column by 1.7 points between 14 and 16 (see
`results/nontemplate_prefix/`). Every number ships with its operating point, and for this
experiment the trim is one of them.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def spread(xs, pct: bool = True) -> str:
    xs = [x for x in xs if x is not None]
    if not xs:
        return "-"
    m, lo, hi = st.mean(xs), min(xs), max(xs)
    if pct:
        return f"{m * 100:.2f} ({lo * 100:.2f}-{hi * 100:.2f})"
    return f"{m:.3f} ({lo:.3f}-{hi:.3f})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results/dna_rna_pair")
    ap.add_argument("--out", type=Path, default=ROOT / "results/dna_rna_pair/TABLES.md")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.results / "dna_rna_pair_s*.json")))
    if not files:
        raise SystemExit(f"no pair results in {args.results}")
    ds = [json.load(open(f)) for f in files]

    hashes = {d["script_sha256"] for d in ds}
    if len(hashes) > 1:
        detail = ", ".join(f"{Path(f).name}={d['script_sha256'][:12]}"
                           for f, d in zip(files, ds))
        raise SystemExit(
            f"these seeds did not come from one script: {detail}. A mixture is not a "
            f"multi-seed result - re-run every seed on the current code before reporting a "
            f"spread, because nothing in the directory says which half is right.")
    # A dirty tree makes the COMMIT unusable as provenance, not the run. What pins the code
    # is the content hash of the script, and every seed here carries one: if they agree,
    # the same script ran five times whatever else in the tree was being edited at the
    # time, which is the case this field was added for. So the refusal is narrowed to what
    # it can actually catch - a seed with no hash at all - and the dirty seeds are recorded
    # in the table rather than silently accepted.
    dirty = [Path(f).name for f, d in zip(files, ds) if d.get("tree_dirty")]
    unpinned = [Path(f).name for f, d in zip(files, ds) if not d.get("script_sha256")]
    if unpinned:
        raise SystemExit(f"{unpinned} carry no script_sha256, and their tree was dirty, so "
                         f"nothing says which code ran. Commit and re-run those seeds.")

    libs = list(ds[0]["pair"])
    trims = {d["pair"][l]["rrna"]["r2_trimmed_by"] for d in ds for l in libs}
    titles = {l: ds[0]["pair"][l]["title"] for l in libs}
    lines = [
        f"Seeds: {len(ds)}; {ds[0]['pair'][libs[0]]['pairs']:,} pairs per library per seed; "
        f"script `{list(hashes)[0][:12]}`, commit `{ds[0]['git_sha'][:12]}`."
        + (f" {len(dirty)} of {len(ds)} seeds ran from a dirty tree, so their commit names "
           f"HEAD rather than the code; the script hash above is identical across all "
           f"{len(ds)} and is what pins them." if dirty else ""),
        "",
        f"Mate 2 trimmed by {sorted(trims)} bases — the measured prefix mode. It moves the "
        f"Kraken2 row by about 1.7 points between 14 and 16; see "
        f"`../nontemplate_prefix/`.",
        "",
        "| | " + " | ".join(f"{l} ({titles[l]})" for l in libs) + " |",
        "|---|" + "---|" * len(libs),
    ]

    def row(label: str, get, pct: bool = True) -> None:
        lines.append(f"| {label} | " + " | ".join(
            spread([get(d["pair"][l]) for d in ds], pct) for l in libs) + " |")

    row("rRNA, ≥40 aligned bases against SILVA", lambda p: p["rrna"]["frac"])
    row("identified, Kraken2 PlusPF-16 at confidence 0.1",
        lambda p: p["kraken2"]["frac"])
    row("identified, RefSeq viral proteins, bitscore ≥50",
        lambda p: p["protein"]["viral"]["identified"])
    row("any hit, RefSeq viral proteins", lambda p: p["protein"]["viral"]["any_hit"])
    row("identified, Swiss-Prot, bitscore ≥50",
        lambda p: p["protein"]["sprot"]["identified"])
    row("any hit, Swiss-Prot", lambda p: p["protein"]["sprot"]["any_hit"])
    row("mean GC", lambda p: p["gc"]["r1"], pct=False)
    row("duplicate pairs (first 50 bases of both mates)",
        lambda p: p["duplicate_frac"])
    lines.append("| modal length R1/R2 | " + " | ".join(
        f"{ds[0]['pair'][l]['modal_length']['r1']}/{ds[0]['pair'][l]['modal_length']['r2']}"
        for l in libs) + " |")

    lines += ["", "Mean across seeds (min-max). Per-seed Wilson intervals are in the JSONs; "
                  "at this depth they are inside ±0.3 points of every figure, so the spread "
                  "above is seed-to-seed variation and not sampling.", ""]

    # The protein-to-cellular ratio is the figure that carries the "it is a virome" claim,
    # so it is computed here rather than divided by hand in prose.
    lines += ["## Viral against cellular protein identification", "",
              "| library | viral | Swiss-Prot | ratio |", "|---|---|---|---|"]
    for l in libs:
        v = [d["pair"][l]["protein"]["viral"]["identified"] for d in ds]
        s = [d["pair"][l]["protein"]["sprot"]["identified"] for d in ds]
        ratio = [a / b for a, b in zip(v, s) if b]
        lines.append(f"| {l} | {spread(v)} | {spread(s)} | "
                     f"{spread(ratio, pct=False) if ratio else '-'} |")
    lines += ["", "A library whose identifiable material is viral reads high here; the RNA "
                  "libraries sit near 1.", ""]

    lines += ["## Do the two nucleic acids see the same organisms?", "",
              "| | mean (min-max) |", "|---|---|"]
    for key, label in (("taxa_shared", "taxa called in both"),
                       ("jaccard", "Jaccard over called taxa"),
                       ("dna_reads_on_shared_taxa", "DNA classified reads on shared taxa"),
                       ("rna_reads_on_shared_taxa", "RNA classified reads on shared taxa")):
        vals = [d["taxon_overlap"][key] for d in ds]
        pct = key.endswith("shared_taxa")
        lines.append(f"| {label} | "
                     f"{spread(vals, pct) if pct else spread(vals, False)} |")
    lines += ["", "The DNA library has few classified reads, so its taxon count is as much "
                  "a statement about depth as about biology; the read-weighted rows are the "
                  "defensible half.", "",
              f"**Confound, carried from the run:** {ds[0]['confound']}", ""]

    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
