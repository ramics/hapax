#!/usr/bin/env python3
"""Per-pile tables for results/depth_rarefaction_piles/ from its committed TSVs.

`results/depth_rarefaction/` answered "what would a doubling of depth buy" on a random
block of the library, and flagged its own limit: a random block is not the residue, so the
ladder bounded the shape rather than giving that pile's value. This reads the same
instrument run with the query block drawn per pair-level PILE, so each pile has its own
number — above all the unplaced residue, the ~15.7% of pairs where neither mate matches
SILVA or any of this sample's own contigs at 40 bases.

The piles come from the residue thread's shared per-mate label table and the pair rule is
stated in `experiments/depth_rarefaction.py`. They partition the library, and `unplaced`
reproduces `results/residue_census/` §2 to the pair before anything is read off it.

Reported per mate, never pooled: mate 2 is 93 nt against mate 1's 107, 62 scorable 31-mers
against 76.

Usage: python3 experiments/report_depth_piles.py [--dir results/depth_rarefaction_piles]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
LEVELS = (0.125, 0.25, 0.5, 1.0)
PILES = ("rrna", "identified", "other_ctg", "tier2", "unplaced")
LABEL = {
    "rrna": "ribosomal",
    "identified": "identified (non-rRNA contig)",
    "other_ctg": "contig outside the three-way split",
    "tier2": "residue on an unidentifiable contig",
    "unplaced": "**residue, unplaced by anything**",
}


def open_maybe_gz(p: Path):
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def find(d: Path, seed: int, f: float) -> Path:
    for name in (f"depth_s{seed}_f{f}.tsv.gz", f"depth_s{seed}_f{f}.tsv"):
        if (d / name).exists():
            return d / name
    raise SystemExit(f"missing depth_s{seed}_f{f}.tsv[.gz] in {d}")


def load(path: Path, mate: str) -> dict[str, dict[str, tuple[int, float]]]:
    """pile -> name -> (median depth, fraction of k-mers with no witness elsewhere).

    Block 9 only: the tool reports raw counts there and own-subtracted counts elsewhere,
    two different definitions of depth. Columns by header name, never by position.
    """
    out: dict[str, dict[str, tuple[int, float]]] = defaultdict(dict)
    with open_maybe_gz(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["block"] != "9":
                continue
            if int(r[f"n{mate}"]) <= 0:
                continue
            pile = r["name"].split("|", 1)[0]
            out[pile][r["name"]] = (int(r[f"med{mate}"]), float(r[f"zero{mate}"]))
    return out


def rng_fmt(vals, pct=True, dp=1):
    a = np.array(vals, dtype=float)
    if pct:
        return f"{100 * a.mean():.{dp}f} ({100 * a.min():.{dp}f}–{100 * a.max():.{dp}f})"
    return f"{a.mean():.{dp}f} ({a.min():.{dp}f}–{a.max():.{dp}f})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/depth_rarefaction_piles")
    args = ap.parse_args()
    d = REPO / args.dir
    man = json.load(open(d / "manifest.json"))
    seeds = man["seeds"]

    data = {}
    for seed in seeds:
        for f in LEVELS:
            p = find(d, seed, f)
            for mate in ("1", "2"):
                data[(seed, f, mate)] = load(p, mate)

    libpairs = {}
    for f in LEVELS:
        with open_maybe_gz(find(d, seeds[0], f)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                libpairs[f] = int(r["library_pairs"])
                break

    out = []
    w = out.append
    w("# The depth ladder, per pile")
    w("")
    w(f"Generated from the committed per-level TSVs by "
      f"`experiments/report_depth_piles.py`. {len(seeds)} seeds, "
      f"{man['params']['per_pile']:,} query pairs per pile per seed, drawn at random "
      f"across the whole library and excluded from every library at every level.")
    w("")
    w("Piles are pair-level, from the residue thread's shared per-mate label table; the "
      "pair rule is in `experiments/depth_rarefaction.py` and `unplaced` reproduces "
      "`results/residue_census/` §2 to the pair (1,372,021 pairs, 15.71%).")
    w("")

    w("## 1. Fully unwitnessed fraction, by pile and effort (mate 1)")
    w("")
    w("| pile | 0.125x | 0.25x | 0.5x | 1x |")
    w("|---|---|---|---|---|")
    for pile in PILES:
        cells = []
        for f in LEVELS:
            vals = []
            for seed in seeds:
                v = data[(seed, f, "1")].get(pile, {})
                zs = np.array([x[1] for x in v.values()], dtype=float)
                vals.append(float((zs >= 1.0).mean()) if zs.size else float("nan"))
            cells.append(rng_fmt(vals))
        w(f"| {LABEL[pile]} | " + " | ".join(cells) + " |")
    w("")

    w("## 2. Median depth, by pile and effort (mate 1)")
    w("")
    w("| pile | 0.125x | 0.25x | 0.5x | 1x |")
    w("|---|---|---|---|---|")
    for pile in PILES:
        cells = []
        for f in LEVELS:
            vals = []
            for seed in seeds:
                v = data[(seed, f, "1")].get(pile, {})
                a = np.array([x[0] for x in v.values()], dtype=float)
                vals.append(float(np.median(a)) if a.size else float("nan"))
            cells.append(rng_fmt(vals, pct=False, dp=0))
        w(f"| {LABEL[pile]} | " + " | ".join(cells) + " |")
    w("")

    w("## 3. The transition: what one doubling rescues, by pile")
    w("")
    w("Of the reads fully unwitnessed at one level, the fraction with any witness one "
      "doubling later. The same named reads scored against both libraries, joined on the "
      "name; no fit and no estimator.")
    w("")
    for mate in ("1", "2"):
        w(f"**Mate {mate}.**")
        w("")
        w("| pile | 0.125x -> 0.25x | 0.25x -> 0.5x | 0.5x -> 1x | n unwitnessed at 0.5x |")
        w("|---|---|---|---|---|")
        for pile in PILES:
            cells, nlast = [], []
            for a, b in zip(LEVELS, LEVELS[1:]):
                vals = []
                for seed in seeds:
                    before = data[(seed, a, mate)].get(pile, {})
                    after = data[(seed, b, mate)].get(pile, {})
                    names = [n for n, (m, z) in before.items() if z >= 1.0]
                    if not names:
                        vals.append(float("nan"))
                        continue
                    g = sum(1 for n in names if after.get(n, (0, 1.0))[1] < 1.0)
                    vals.append(g / len(names))
                    if a == LEVELS[-2]:
                        nlast.append(len(names))
                cells.append(rng_fmt(vals))
            n_txt = f"{int(np.mean(nlast)):,}" if nlast else "0"
            w(f"| {LABEL[pile]} | " + " | ".join(cells) + f" | {n_txt} |")
        w("")

    # the headline comparison, computed rather than typed
    lib = REPO / "results" / "depth_rarefaction" / "TABLES.md"
    unpl = []
    for seed in seeds:
        before = data[(seed, LEVELS[-2], "1")].get("unplaced", {})
        after = data[(seed, LEVELS[-1], "1")].get("unplaced", {})
        names = [n for n, (m, z) in before.items() if z >= 1.0]
        unpl.append(sum(1 for n in names if after.get(n, (0, 1.0))[1] < 1.0) / len(names))
    # the qualifier that decides what the rescue is worth: a witness is not a depth
    deep = {}
    for pile in PILES:
        vals = []
        for seed in seeds:
            before = data[(seed, LEVELS[-2], "1")].get(pile, {})
            after = data[(seed, LEVELS[-1], "1")].get(pile, {})
            names = [n for n, (m, z) in before.items() if z >= 1.0]
            if not names:
                continue
            vals.append(sum(1 for n in names
                            if after.get(n, (0, 1.0))[0] >= 5) / len(names))
        deep[pile] = vals

    w("## 4. A witness is not a depth")
    w("")
    w("Of the same fully-unwitnessed reads, the fraction reaching median depth 5 or more "
      "one doubling later — the depth where assembly and classification begin to work at "
      "all (`results/spike_abundance/`). Mate 1.")
    w("")
    w("| pile | gained a witness | reached depth >= 5 |")
    w("|---|---|---|")
    for pile in PILES:
        if not deep[pile]:
            continue
        vals = []
        for seed in seeds:
            before = data[(seed, LEVELS[-2], "1")].get(pile, {})
            after = data[(seed, LEVELS[-1], "1")].get(pile, {})
            names = [n for n, (m, z) in before.items() if z >= 1.0]
            if names:
                vals.append(sum(1 for n in names
                                if after.get(n, (0, 1.0))[1] < 1.0) / len(names))
        w(f"| {LABEL[pile]} | {rng_fmt(vals)} | {rng_fmt(deep[pile])} |")
    w("")

    w("## 5. The number this run exists for")
    w("")
    w(f"**The unplaced residue's own rescue rate over the last measured doubling is "
      f"{rng_fmt(unpl)}% (mate 1, {len(seeds)} seeds).** "
      f"`results/depth_rarefaction/` measured 44.0 (40.6–46.2) on a random block of the "
      f"library, and flagged that a random block is not the residue. This replaces that "
      f"figure for the pile the question was about; the library figure stands for the "
      f"library.")
    w("")
    if lib.exists():
        w(f"The library-wide run is [`../depth_rarefaction/TABLES.md`]"
          f"(../depth_rarefaction/TABLES.md), same instrument, same levels, same seeds.")
        w("")

    w("## Provenance")
    w("")
    w(f"- generator `{man['script']}`, content hash `{man['script_sha256'][:16]}`")
    w(f"- this report `experiments/report_depth_piles.py`, content hash "
      f"`{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}`")
    w(f"- git SHA `{man['git_sha']}`")
    w(f"- label table `{man['params']['classes_table']}`")
    w(f"- k = {man['params']['k']}, levels {man['params']['levels']}, "
      f"mate 2 trimmed {man['params']['trim_r2']} nt")
    w(f"- library sizes at each level: " + ", ".join(
        f"{f:g}x = {libpairs[f]:,}" for f in LEVELS))
    w(f"- total runtime {man['total_s']:.0f}s")
    w("")

    (d / "TABLES.md").write_text("\n".join(out) + "\n")
    print(f"wrote {d / 'TABLES.md'}")
    print(f"unplaced residue 0.5x->1x rescue: {rng_fmt(unpl)}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
