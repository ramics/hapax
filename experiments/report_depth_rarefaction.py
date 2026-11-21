#!/usr/bin/env python3
"""Tables for results/depth_rarefaction/ from its committed per-level TSVs.

Separate from the generator on purpose. `depth_rarefaction.py` produces the measurement
and its manifest pins the code that produced it; this script only reads those committed
files, so the analysis can be extended without invalidating the manifest's script hash.
That is the repo's rule about a provenance change forcing a re-run, applied by not
needing one.

Two readouts:

1. THE CURVE. The fully-unwitnessed fraction of a fixed set of real molecules against the
   size of the library around them, over a measured 8x range.

2. THE TRANSITION, which is the better number and is pure measurement. Of the reads that
   are fully unwitnessed at one level, what fraction has gained a witness one doubling
   later? No curve fit, no extrapolation, no estimator: the same named reads are scored
   against both libraries and the join is on the name. It answers the actual question —
   "if we sequenced this water twice as deep, how much of what we now see once would we
   then see twice" — rather than a proxy for it.

   Reported per transition rather than pooled, because the rate is NOT constant: each
   doubling adds twice the absolute number of reads of the one before it, and witnessing
   depends on absolute depth, so the rate rises along the ladder. That direction makes the
   last measured transition a floor for the next one rather than a ceiling, and it is why
   this script will still not print a 10x figure.

Usage: python3 experiments/report_depth_rarefaction.py [--dir results/depth_rarefaction]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
LEVELS = (0.125, 0.25, 0.5, 1.0)


def load_level(path: Path, mate: str) -> dict[str, tuple[int, float]]:
    """name -> (median depth, fraction of k-mers with no witness) for QUERY_ONLY rows.

    Columns by name from the file's own header; block 9 only, because the tool reports
    raw counts for block 9 and own-subtracted counts elsewhere, which are two different
    definitions of depth.
    """
    out = {}
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["block"] != "9":
                continue
            if int(r[f"n{mate}"]) > 0:
                out[r["name"]] = (int(r[f"med{mate}"]), float(r[f"zero{mate}"]))
    return out


def rng_fmt(vals, pct=True, dp=1):
    a = np.array(vals, dtype=float)
    if pct:
        return f"{100 * a.mean():.{dp}f} ({100 * a.min():.{dp}f}–{100 * a.max():.{dp}f})"
    return f"{a.mean():.{dp}f} ({a.min():.{dp}f}–{a.max():.{dp}f})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/depth_rarefaction")
    args = ap.parse_args()
    d = REPO / args.dir
    man = json.load(open(d / "manifest.json"))
    seeds = man["seeds"]

    data = {}       # (seed, level, mate) -> dict
    for seed in seeds:
        for f in LEVELS:
            p = d / f"depth_s{seed}_f{f}.tsv.gz"
            if not p.exists():
                p = d / f"depth_s{seed}_f{f}.tsv"
            if not p.exists():
                raise SystemExit(f"missing depth_s{seed}_f{f}.tsv[.gz] in {d}")
            for mate in ("1", "2"):
                data[(seed, f, mate)] = load_level(p, mate)

    libpairs = {}
    for f in LEVELS:
        q = d / f"depth_s{seeds[0]}_f{f}.tsv.gz"
        if not q.exists():
            q = d / f"depth_s{seeds[0]}_f{f}.tsv"
        with (gzip.open if q.suffix == ".gz" else open)(q, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                libpairs[f] = int(r["library_pairs"])
                break

    out = []
    w = out.append
    w("# What deeper sequencing of this library would buy")
    w("")
    w(f"Generated from the committed per-level TSVs by "
      f"`experiments/report_depth_rarefaction.py`. {len(seeds)} seeds, each drawing its "
      f"own query block of {man['params']['query_pairs']:,} real pairs AND its own nested "
      f"subsample, so the spread covers both. Query pairs are excluded from every library "
      f"at every level, so a depth is always occurrences elsewhere.")
    w("")

    w("## 1. The curve")
    w("")
    w("| effort | library pairs | mate | median depth | fully unwitnessed | depth < 5 |")
    w("|---|---|---|---|---|---|")
    for f in LEVELS:
        for mate in ("1", "2"):
            med, z0, lt5 = [], [], []
            for seed in seeds:
                v = np.array([x[0] for x in data[(seed, f, mate)].values()], dtype=float)
                zs = np.array([x[1] for x in data[(seed, f, mate)].values()], dtype=float)
                med.append(np.median(v))
                z0.append(float((zs >= 1.0).mean()))
                lt5.append(float((v < 5).mean()))
            w(f"| {f:g}x | {libpairs[f]:,} | {mate} | {rng_fmt(med, pct=False, dp=0)} | "
              f"{rng_fmt(z0)} | {rng_fmt(lt5)} |")
    w("")
    w("'Fully unwitnessed' is every scorable 31-mer of the read having no witness anywhere "
      "else in the library — the strict version, not 'median depth zero'. A read with an "
      "internal zero run of about 31 in otherwise covered sequence is a substitution "
      "error, not a rare molecule, and the strict definition excludes it by construction.")
    w("")

    w("## 2. The transition: what one doubling actually rescues")
    w("")
    w("Of the reads fully unwitnessed at one level, the fraction with any witness one "
      "doubling later. The same named reads scored against both libraries; no fit and no "
      "estimator.")
    w("")
    w("| doubling | mate | n unwitnessed before | gained a witness | reached depth >= 5 |")
    w("|---|---|---|---|---|")
    trans = {}
    for a, b in zip(LEVELS, LEVELS[1:]):
        for mate in ("1", "2"):
            n_before, gained, deep = [], [], []
            for seed in seeds:
                before = data[(seed, a, mate)]
                after = data[(seed, b, mate)]
                names = [n for n, (m, z) in before.items() if z >= 1.0]
                n_before.append(len(names))
                g = sum(1 for n in names if after.get(n, (0, 1.0))[1] < 1.0)
                dp = sum(1 for n in names if after.get(n, (0, 1.0))[0] >= 5)
                gained.append(g / max(1, len(names)))
                deep.append(dp / max(1, len(names)))
            trans[(a, b, mate)] = gained
            w(f"| {a:g}x -> {b:g}x | {mate} | {int(np.mean(n_before)):,} | "
              f"{rng_fmt(gained)} | {rng_fmt(deep)} |")
    w("")
    last = trans[(LEVELS[-2], LEVELS[-1], "1")]
    w(f"The rate RISES along the ladder, because each doubling adds twice the absolute "
      f"number of reads of the one before it while witnessing depends on absolute depth. "
      f"So the last measured transition — {rng_fmt(last)}% of mate 1's fully unwitnessed "
      f"reads gaining a witness when the library goes from half to full — is a FLOOR for "
      f"what the next doubling would do, not a ceiling.")
    w("")
    w("No 10x figure is reported. Interpolation inside the measured 8x range is "
      "measurement; one doubling is a short extrapolation anchored on three measured "
      "transitions; ten times is neither, and this repo has already withdrawn one claim "
      "for reading a curve past where it was measured.")
    w("")

    w("## Stated biases, with direction")
    w("")
    w("- **Subsampling a finished library models re-sequencing the SAME library deeper, "
      "not sequencing new water.** The same plant's rRNA fraction went 38.0% to 67.1% in "
      "eight months, so new water is a different population. Direction: optimistic about "
      "new water, roughly right about deeper sequencing of this library.")
    w("- **Duplicates read as depth.** A rescued read whose witness is a PCR duplicate of "
      "its own molecule has gained confirmation that it is not an error, but no new "
      "sequence. Direction: inflates the rescue rate as a measure of assembly gain, while "
      "remaining honest as a measure of 'seen more than once'.")
    w("- **A random query block is not the residue.** These are real library reads of every "
      "kind, so the curve describes the library. The unplaced residue sits at the bottom "
      "of this distribution by definition and its own transition rate will be lower than "
      "the library average; the ladder here bounds the shape, not that pile's value.")
    w("")
    w("## Provenance")
    w("")
    w(f"- generator `{man['script']}`, content hash `{man['script_sha256'][:16]}`")
    w(f"- this report `experiments/report_depth_rarefaction.py`, content hash "
      f"`{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}`")
    w(f"- git SHA `{man['git_sha']}`")
    w(f"- k = {man['params']['k']}, levels {man['params']['levels']}, "
      f"{man['params']['query_pairs']:,} query pairs per seed, mate 2 trimmed "
      f"{man['params']['trim_r2']} nt")
    w(f"- total runtime {man['total_s']:.0f}s; it shared the machine with another job, so "
      f"the timings in `manifest.json` are operational, not benchmarks")
    w("")

    (d / "TABLES.md").write_text("\n".join(out) + "\n")
    print(f"wrote {d / 'TABLES.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
