"""Turn the big-index run's rows into the three things it was launched for.

The instance searches and nothing else (`experiments/bigindex_search.py`), so everything
that needs a judgement happens here, where it can be re-run in seconds and a threshold can
be changed without another launch and another approval. Three outputs:

1. **IDENTITY DISTRIBUTIONS, residue against identified control.** Not a hit rate: the
   BLAST ground truth already established that nearly everything hits something, so a hit
   rate carries almost no information. What separates the piles, every time it has been
   measured, is continuous - identity AND the fraction of the read that aligns - so this
   reports the joint distribution in bins with Wilson intervals, each pile beside the other.

2. **REQUIREMENT 1.iii**: "evidence that no-match sequence really sits below ~80% identity."
   The honest version of that question needs both axes at once. A read whose best alignment
   covers a third of its length at 95% identity has NOT been shown to be a 95% match to
   anything, and a threshold on identity alone would count it as one. So the readout is the
   fraction of each pile with a strong alignment - high identity over most of the read - and
   the fraction with nothing above the floor, with the operating point stated on the line.

3. **THE RESIDUE NAME LISTS** the structure thread fixed the format of before launch, in
   `results/subtract_known/BIGINDEX_RUN_PLAN.md` section 4.1: cumulative (survivors of the
   existing tier-1 cascade PLUS this rung), one name per line as the FASTQ carries it, no
   mate suffix, deduplicated, with a JSON beside each file carrying the operating point a
   name list cannot itself carry.

PAIR LEVEL, AND BOTH MATES. A pair is removed by a rung if EITHER mate is identified, which
is why the block query sets carry both mates: a mate-1-only pass would under-identify and
every name list would be wrong in the direction that invents residue. The per-mate split is
reported too, because mate 2 is 93 nt against mate 1's 107 after the non-template trim and
this repo has already found that every mate-interchangeable field leans on that difference.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import platform
import statistics as st
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

#: Identity bins for the distribution. The edges are the ones the argument is about: 80% is
#: Imo's own number in requirement 1.iii, and the spike ladder puts recognition's wall
#: between 10 and 20% divergence, i.e. 80-90% identity, so these bins are where the two
#: independent lines of evidence have to meet.
IDENTITY_BINS = (0, 50, 60, 70, 80, 90, 95, 100)

#: Aligned-fraction bins. Separate axis, not folded into identity, because the contig
#: anatomy measured these as the two things that actually separate the piles (0.315 vs 0.567
#: aligned fraction at 35.0% vs 40.6% identity) and collapsing them loses the distinction.
ALN_FRAC_BINS = (0.0, 0.25, 0.5, 0.75, 1.01)


def wilson(k: int, n: int) -> tuple[float, float]:
    """95% interval. Quoted rather than the point estimate, per this repo's rule."""
    if n == 0:
        return (0.0, 0.0)
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def binof(value: float, edges) -> str:
    for lo, hi in zip(edges, edges[1:]):
        if lo <= value < hi:
            return f"{lo}-{hi}"
    return f">={edges[-1]}"


def pair_name(read: str) -> str:
    """`SRR23973193.4/1` -> `SRR23973193.4`. The name lists carry no mate suffix."""
    return read.rsplit("/", 1)[0] if "/" in read else read


def load(paths: list[str]) -> list[dict]:
    rows = []
    for p in paths:
        with open(p) as h:
            for line in h:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    # A half-written final line means the instance was killed mid-write,
                    # which is a fact about the run rather than odd data. Dropped here and
                    # counted in the report rather than silently averaged in.
                    rows.append({"_truncated": True})
    return rows


def sig_block(hits: list[dict], n: int, cutoff: float) -> dict:
    """What survives an E-value cutoff, and what shape it has.

    Reported as its own block rather than folded into the headline so that the two readings
    stay visible side by side: this report published the unfiltered one, and the difference
    between them is the finding, not a footnote.
    """
    sig = [r for r in hits if r.get("evalue") is not None and r["evalue"] < cutoff]
    part = [r for r in sig if r.get("aln_frac", 1.0) < 0.95]
    out = {"cutoff": cutoff, "reads": len(sig),
           "frac_of_queries": len(sig) / n if n else None,
           "frac_of_queries_ci95": wilson(len(sig), n),
           "partial_reads": len(part),
           "partial_frac_of_queries": len(part) / n if n else None}
    if sig:
        out["median_aln_frac"] = round(st.median([r["aln_frac"] for r in sig]), 4)
        out["median_pident"] = round(st.median([r["pident"] for r in sig]), 2)
    if part:
        out["partial_median_pident"] = round(st.median([r["pident"] for r in part]), 2)
        out["partial_median_aln_len"] = st.median([r["aln_len"] for r in part])
    if hits:
        out["median_evalue_all_hits"] = st.median([r["evalue"] for r in hits
                                                   if r.get("evalue") is not None])
    return out


def identified(row: dict, min_pident: float, min_aln: int,
               min_bitscore: float) -> bool:
    """Does this row clear the operating point?

    Both a length and an identity condition, never identity alone. A 20-base alignment at
    100% identity is noise at these read lengths, and an operating point that would count it
    is how a removal number gets inflated.
    """
    if not row.get("hit"):
        return False
    if row.get("bitscore", 0.0) < min_bitscore:
        return False
    return (row.get("pident", 0.0) >= min_pident
            and row.get("aln_len", 0) >= min_aln)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", nargs="+", required=True,
                    help="rows_*.jsonl from the run, one or more globs")
    ap.add_argument("--rung", required=True,
                    help="what this search adds to the cascade: core_nt, nr, core_nt_nr")
    ap.add_argument("--block", type=int, required=True,
                    help="which block (seed) these rows are from; names the outputs")
    ap.add_argument("--db-version", default="", help="copied into every operating point")
    ap.add_argument("--program", default="", help="blastn, blastx, diamond-blastx")
    ap.add_argument("--min-pident", type=float, default=90.0)
    ap.add_argument("--min-aln", type=int, default=40,
                    help="bases. 40 is the cascade's own nucleotide match floor, so a rung "
                         "added here is comparable with the rungs already in the curve")
    ap.add_argument("--min-bitscore", type=float, default=0.0)
    ap.add_argument("--strong-pident", type=float, default=80.0,
                    help="requirement 1.iii's own number")
    ap.add_argument("--strong-aln-frac", type=float, default=0.8)
    ap.add_argument("--significant", type=float, default=1e-5,
                    help="E-value below which a match is treated as real rather than as the "
                         "best chance alignment a short query finds in a large database")
    ap.add_argument("--residue-names", type=Path, default=None,
                    help="the cascade's own residue name list for this block, so the rung "
                         "is CUMULATIVE on it rather than measured against nothing. "
                         "Without it the output says so and lists only what these rows saw")
    ap.add_argument("--out", type=Path, default=ROOT / "results/bigindex")
    args = ap.parse_args()

    paths = sorted({p for g in args.rows for p in glob.glob(g)})
    if not paths:
        raise SystemExit(f"no files matched {args.rows}")
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    rows = load(paths)
    truncated = sum(1 for r in rows if r.get("_truncated"))
    rows = [r for r in rows if not r.get("_truncated")]

    # ---- distributions, per pile, per axis
    per_class: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        per_class[r.get("class") or "unlabelled"].append(r)

    dist: dict[str, dict] = {}
    for cls, rs in sorted(per_class.items()):
        n = len(rs)
        hits = [r for r in rs if r.get("hit")]
        ident = defaultdict(int)
        for r in hits:
            ident[binof(r.get("pident", 0.0), IDENTITY_BINS)] += 1
        frac = defaultdict(int)
        for r in hits:
            frac[binof(r.get("aln_frac", 0.0), ALN_FRAC_BINS)] += 1
        strong = [r for r in hits
                  if r.get("pident", 0.0) >= args.strong_pident
                  and r.get("aln_frac", 0.0) >= args.strong_aln_frac]
        dist[cls] = {
            "queries": n,
            "any_hit": len(hits) / n if n else None,
            "any_hit_ci95": wilson(len(hits), n),
            # 1.iii, stated as the two-axis quantity it has to be
            "strong": len(strong) / n if n else None,
            "strong_ci95": wilson(len(strong), n),
            "identity_bins": dict(sorted(ident.items())),
            "aln_frac_bins": dict(sorted(frac.items())),
            "median_pident": round(st.median([r["pident"] for r in hits]), 2) if hits else None,
            "median_aln_frac": round(st.median([r["aln_frac"] for r in hits]), 4)
            if hits else None,
            "identified_at_operating_point":
                sum(1 for r in rs if identified(r, args.min_pident, args.min_aln,
                                                args.min_bitscore)) / n if n else None,
            # SIGNIFICANCE, which the first version of this report did not read and which
            # inverts its conclusion. BLAST reports down to E = 10 by default, so "any hit"
            # against a 318.8 GB database counts the best CHANCE match a 107 nt query finds
            # - median E = 0.24 on the first pilot. These fields are therefore not an extra
            # cut of the same quantity, they are the quantity: a relative that is real
            # aligns end to end, and partial alignment is mostly what chance looks like.
            "significant": sig_block(hits, n, args.significant),
        }

    # ---- per mate, because the two mates are not interchangeable here
    per_mate: dict[str, dict] = {}
    for mate in ("1", "2"):
        rs = [r for r in rows if str(r.get("read", "")).endswith(f"/{mate}")
              and (r.get("class") == "residue")]
        if not rs:
            continue
        k = sum(1 for r in rs if identified(r, args.min_pident, args.min_aln,
                                            args.min_bitscore))
        per_mate[mate] = {"queries": len(rs), "identified": k / len(rs),
                          "ci95": wilson(k, len(rs))}

    # ---- the name lists: pairs of the residue that this rung does NOT identify
    residue_pairs: set[str] = set()
    removed_pairs: set[str] = set()
    for r in rows:
        if r.get("class") != "residue":
            continue
        name = pair_name(str(r.get("read", "")))
        residue_pairs.add(name)
        if identified(r, args.min_pident, args.min_aln, args.min_bitscore):
            removed_pairs.add(name)

    cascade_names: set[str] | None = None
    if args.residue_names and args.residue_names.exists():
        lines = args.residue_names.read_text().splitlines()
        if args.residue_names.suffix == ".tsv":
            # The structure thread's coverage table is the rung-zero reference that exists
            # in the repo, and its first column is the pair name in exactly this form. Read
            # it directly rather than asking somebody to maintain a second copy as a .txt:
            # a second copy of a reference list is how the codon-table bug happened.
            lines = [ln.split("\t")[0] for ln in lines[1:]]
        cascade_names = {ln.strip() for ln in lines if ln.strip()}
        # The rung-zero cross-check the plan asks for: if the pairs these rows cover and the
        # committed residue list disagree, the blocks drifted, and that is better found in a
        # file than in a result.
        drift = len(cascade_names ^ residue_pairs)
        survivors = sorted(cascade_names - removed_pairs)
    else:
        drift = None
        survivors = sorted(residue_pairs - removed_pairs)

    names_path = args.out / f"residue_names_{args.rung}_s{args.block}.txt"
    names_path.write_text("".join(f"{n}\n" for n in survivors))

    operating = {
        "file": names_path.name,
        "contains": f"pairs still unidentified after the tier-1 cascade PLUS {args.rung}; "
                    f"cumulative, not a delta",
        "rung": args.rung, "block": args.block,
        "program": args.program, "db_version": args.db_version,
        "min_pident": args.min_pident, "min_aln_bases": args.min_aln,
        "min_bitscore": args.min_bitscore,
        "pair_rule": "a pair is identified if EITHER mate clears the operating point",
        "count": len(survivors),
        "residue_pairs_in": len(cascade_names) if cascade_names is not None
        else len(residue_pairs),
        "removed_by_this_rung": len(removed_pairs),
        "name_drift_vs_committed_list": drift,
        "source_files": [Path(p).name for p in paths],
        "truncated_rows_dropped": truncated,
    }
    (args.out / f"residue_names_{args.rung}_s{args.block}.json").write_text(
        json.dumps(operating, indent=2))

    out = {"rung": args.rung, "block": args.block, "rows": len(rows),
           "truncated_rows_dropped": truncated,
           "distributions": dist, "residue_per_mate": per_mate,
           "operating_point": operating,
           "limits": [
               "identity is the identity of the ALIGNED PART only, so it is read beside "
               "aln_frac and never alone: a third of a read at 95% is not a 95% match",
               "the control pile is the cascade's already-identified reads, which is a "
               "contrast rather than a ground truth: a database this size has its own gaps",
               "one block is one seed; the repo rule is five, and the run adds blocks only "
               "while they fit its measured rate, so the seed count is whatever ran",
           ],
           "wall_seconds": time.perf_counter() - t0}
    (args.out / f"bigindex_{args.rung}_s{args.block}.json").write_text(
        json.dumps(out, indent=2))
    (args.out / f"manifest_{args.rung}_s{args.block}.json").write_text(json.dumps({
        "script": "experiments/bigindex_report.py", **provenance(__file__),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    print(f"{len(rows):,} rows from {len(paths)} file(s); rung {args.rung}, block {args.block}")
    for cls, d in dist.items():
        print(f"  {cls:12} n={d['queries']:6} any hit {d['any_hit']:.3f} "
              f"strong(>={args.strong_pident}% over >={args.strong_aln_frac:.0%}) "
              f"{d['strong']:.4f} "
              f"median pident {d['median_pident']} frac {d['median_aln_frac']}")
    for mate, d in per_mate.items():
        print(f"  residue mate {mate}: identified {d['identified']:.4f} "
              f"[{d['ci95'][0]:.4f}-{d['ci95'][1]:.4f}] on {d['queries']} reads")
    print(f"  {len(survivors):,} pairs survive {args.rung} -> {names_path}")
    if drift:
        print(f"  WARNING: {drift} names differ from the committed residue list for this "
              f"block; the blocks may have drifted")


if __name__ == "__main__":
    main()
