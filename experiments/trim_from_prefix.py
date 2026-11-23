#!/usr/bin/env python3
"""Read `nontemplate_prefix_map.py`'s seeds and print the trim, or refuse to name one.

WHY THIS IS A SCRIPT. The trim of a new library is the first number a production run needs
and the one nothing downstream can check. For SRR23973193 a human read five JSON files and
wrote `trim_r2: 14` into `experiments/libraries.json` with a sentence of prose; that is
exactly the step that cannot be re-run, and on the next library it is the step most likely
to be done by eye. The rule it applied is small enough to write down, so it is written down
here and applied by code, which means a library where the rule does NOT hold comes back as
a refusal with a reason instead of a number somebody chose.

THE RULE, AND WHERE EACH PART COMES FROM (`results/nontemplate_prefix/`):

1. **The mode of mate 2's 5' soft clip is the candidate.** It is the length at which the
   clip spikes, measured on untrimmed reads against SILVA with the body of the same read as
   the anchor.
2. **It must be the same mode in every seed.** Berlin's was 14 in 5 of 5. Seeds are
   disjoint blocks of one flow cell, so disagreement between them is not sampling noise
   about a library-wide constant: it means there is no single length to cut.
3. **Mate 1's mode must be 0.** Mate 1 has ~2 bases of 5' excess and no mode, and that
   difference is what separates a non-template prefix from Illumina's priming bias. A
   library where both mates spike at the same length has an instrument property, not a kit
   prefix, and cutting one mate would be arbitrary.
4. **The mode must carry a real share of the anchored reads.** Berlin's mode share is
   0.507-0.518: the prefix is a distribution, not a constant, and a "mode" holding 3% of
   reads is the top of a flat histogram. The floor is a stated operating point, not a
   measured threshold - the only calibration behind it is that the one library known to
   have a prefix sits at 0.51 and the one known not to (the DNA library) has no spike at
   all.
5. **Enough reads must anchor to say anything.** A library with almost no ribosomal content
   gives SILVA nothing to anchor to, and its clip distribution is then a statement about
   the reference rather than about the library - which is the committed reading of the DNA
   library's row. Below the floor this refuses rather than reporting its mode.

A refusal is the useful answer. It says: measure this library another way, or state the
trim on the command line and own it, and SPEC 4.1 allows exactly that ("stating the value
on the command line always wins and is never checked, because stating it is the override").

The output is one line on stdout, the integer, so a shell can read it; everything else goes
to stderr and to the JSON written beside the seeds.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I drops the script dir

from provenance import provenance

#: Operating points, stated because thresholds dominate: both are the shape of Berlin's
#: measurement rather than a measured cut, and both are flags on the command line so a run
#: that moves them says so in its manifest.
MIN_MODE_SHARE = 0.25
MIN_ANCHORED = 2_000


def decide(rows: list[dict], library: str, min_share: float, min_anchored: int) -> dict:
    """The five rules of the docstring, each with the number that failed it."""
    problems, modes, shares, anchored = [], [], [], []
    for r in rows:
        row = r["libraries"].get(library)
        if row is None:
            problems.append(f"seed file has no row for {library}: it holds "
                            f"{', '.join(sorted(r['libraries']))}")
            continue
        r2 = row.get("r2_five_prime_clip") or {}
        r1 = row.get("r1_five_prime_clip") or {}
        n = (row.get("anchored") or {}).get("r2")
        anchored.append(n)
        if not r2 or n is None:
            problems.append("a seed anchored no mate-2 reads at all")
            continue
        modes.append(r2.get("mode"))
        shares.append(r2.get("mode_share"))
        if n < min_anchored:
            problems.append(
                f"only {n} mate-2 reads anchored to SILVA in one seed, under the floor of "
                f"{min_anchored}: the clip distribution of a library with little ribosomal "
                f"content describes the reference, not the library")
        if r1.get("mode") not in (0, None) and r1.get("mode") == r2.get("mode"):
            problems.append(
                f"both mates spike at {r1['mode']}: that is an instrument property rather "
                f"than one mate's kit prefix, and cutting one mate would be arbitrary")
        if r2.get("mode_share") is not None and r2["mode_share"] < min_share:
            problems.append(
                f"the mate-2 mode {r2.get('mode')} holds only "
                f"{r2['mode_share']:.1%} of anchored reads, under the floor of "
                f"{min_share:.0%}: that is the top of a flat histogram, not a prefix")
    distinct = sorted({m for m in modes if m is not None})
    if len(distinct) > 1:
        problems.append(
            f"the mate-2 5' mode is not the same in every seed ({distinct}). Seeds are "
            f"disjoint blocks of one flow cell, so that is not sampling noise about a "
            f"library-wide constant: there is no single length to cut")
    if not distinct:
        problems.append("no seed produced a mate-2 5' mode")
    return {"library": library, "seeds": len(rows), "modes": modes,
            "mode_shares": shares, "anchored_r2": anchored,
            "r1_modes": [(r["libraries"].get(library, {}).get("r1_five_prime_clip") or {})
                         .get("mode") for r in rows],
            "trim_r2": distinct[0] if (distinct and not problems) else None,
            "problems": problems,
            "min_mode_share": min_share, "min_anchored": min_anchored}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--prefix-json", type=Path, nargs="+", required=True,
                    help="nontemplate_prefix_s*.json, one per seed")
    ap.add_argument("--library", required=True,
                    help="the key inside those files, i.e. the accession")
    ap.add_argument("--min-mode-share", type=float, default=MIN_MODE_SHARE)
    ap.add_argument("--min-anchored", type=int, default=MIN_ANCHORED)
    ap.add_argument("--out", type=Path, default=None,
                    help="where to write the decision; defaults beside the first seed")
    args = ap.parse_args()

    rows = []
    for p in args.prefix_json:
        if not p.exists():
            raise SystemExit(f"{p} is missing")
        rows.append(json.loads(p.read_text()))
    if len(rows) < 2:
        raise SystemExit("one seed cannot show that the mode is the same in every seed; "
                         "pass every seed's file")

    d = decide(rows, args.library, args.min_mode_share, args.min_anchored)
    out = {"script": "experiments/trim_from_prefix.py", **provenance(__file__),
           "parameters": {"library": args.library,
                          "prefix_json": [str(p) for p in args.prefix_json]},
           **d,
           "rule": "the mate-2 5' soft-clip mode, required to be identical in every seed, "
                   "to hold at least min_mode_share of anchored reads, to sit beside a "
                   "mate-1 mode that is not the same number, and to rest on at least "
                   "min_anchored anchored reads per seed",
           "bias": "adapter read-through clips 3' ends and inflates the within-read "
                   "control, which makes a real 5' excess harder to see rather than "
                   "easier. Every failure here therefore errs towards refusing to trim, "
                   "and not trimming leaves kit sequence in the reads"}
    dst = args.out or args.prefix_json[0].parent / f"trim_from_prefix_{args.library}.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")

    if d["problems"]:
        print(f"trim_from_prefix: REFUSING to name a trim for {args.library}:",
              file=sys.stderr)
        for p in dict.fromkeys(d["problems"]):
            print(f"  - {p}", file=sys.stderr)
        print(f"  modes {d['modes']}, shares {d['mode_shares']}, "
              f"anchored {d['anchored_r2']}; written {dst}", file=sys.stderr)
        return 2
    lo, hi = min(d["mode_shares"]), max(d["mode_shares"])
    print(f"trim_from_prefix: {args.library} mode {d['trim_r2']} in {d['seeds']}/"
          f"{d['seeds']} seeds, mode share {lo:.3f}-{hi:.3f}, "
          f"mate-1 modes {d['r1_modes']}; written {dst}", file=sys.stderr)
    print(d["trim_r2"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
