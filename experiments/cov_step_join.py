"""Does within-read coverage step mean "this molecule spans an abundant and a rare region"?

A cross-instrument test, requested by the thread measuring whether the residue is
structured. Its `cov_step` field - the largest step in exact 31-mer coverage along a read -
reaches 0.780 AUC separating residue from subtracted reads, and the chimaera reading of it
is dead (the contig-level anatomy fires on 0% of both unidentifiable and identified
contigs). The surviving explanation is coverage that is genuinely uneven along a molecule.

THE PREDICTION, stated so it can fail. If `cov_step` is uneven coverage along a molecule,
then within the residue, reads with a large step should be the ones whose two mates land on
assembled contigs of markedly different depth, and the rank correlation between `cov_step`
and |log2 of the ratio of the two mates' contig depths| should be clearly positive. If it
comes out flat, `cov_step` is measuring something else and the field list should say so.

WHY THIS IS WORTH RUNNING DESPITE A WEAK PRIOR. `mate_disagree`, the k-mer-level version of
the same idea, reaches only 0.762 AUC, so a strong relationship is not expected. The value
is that the two quantities come from different instruments: theirs from exact k-mer counting
over the whole library, this one from mapping to the assembly and reading SPAdes' per-contig
coverage out of the contig name. Agreement across instruments is evidence about the
mechanism; agreement within one instrument is evidence about the instrument.

CONTROLS, because a positive correlation here is easy to manufacture. `cov_spread` and
`mate_disagree` are correlated against the same depth ratio, so a correlation that is the
same for every field is a property of the mapping, not of `cov_step`. A shuffled join
(read names permuted) gives the null. And the pairs whose mates land on DIFFERENT contigs
are reported separately from those on the same contig, since the same-contig case can only
show depth difference through mapping error.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]
COV = re.compile(r"_cov_([0-9.]+)")


def read_fasta_mates(path: Path) -> dict[str, dict[str, str]]:
    """{read name: {"1": seq, "2": seq}} from a FASTA whose records are name/1 and name/2."""
    out: dict[str, dict[str, str]] = {}
    name = mate = None
    buf: list[str] = []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    out.setdefault(name, {})[mate] = "".join(buf)
                head = line[1:].strip()
                name, _, mate = head.partition("/")
                mate = mate or "1"
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out.setdefault(name, {})[mate] = "".join(buf)
    return out


def best_contig(idx, seq: str, min_mlen: int):
    best = None
    for h in idx.map(seq):
        if h.mlen < min_mlen:
            continue
        if best is None or h.mlen > best[1]:
            best = (h.ctg, h.mlen)
    return best


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    return num / (da * db) if da and db else 0.0


def spearman(a, b):
    return pearson(rank(a), rank(b))


def fisher_ci(rho: float, n: int):
    """95% interval on a rank correlation. Wide enough to matter at these n, and the point
    of reporting it is that "clearly positive" has to mean an interval clear of zero."""
    if n < 4 or abs(rho) >= 1:
        return [None, None]
    z = 0.5 * math.log((1 + rho) / (1 - rho))
    se = 1 / math.sqrt(n - 3)
    return [round(math.tanh(z - 1.96 * se), 4), round(math.tanh(z + 1.96 * se), 4)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--axes", type=Path,
                    default=ROOT / "results/residue_structure/coverage_axes_s1.tsv",
                    help="the sibling thread's per-read coverage axes, keyed by read name")
    ap.add_argument("--reads", type=Path,
                    default=Path("/tmp/claude-0/wk/recut/residue_s1.fa"))
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/covstep"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/cov_step_join")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    axes: dict[str, dict[str, float]] = {}
    with open(args.axes) as h:
        header = h.readline().rstrip("\n").split("\t")
        for line in h:
            f = line.rstrip("\n").split("\t")
            axes[f[0]] = {k: float(v) for k, v in zip(header[1:], f[1:])}

    mates = read_fasta_mates(args.reads)

    import mappy
    mmi = args.work / "contigs.mmi"
    if not mmi.exists():
        idx = mappy.Aligner(str(args.contigs), preset="sr", n_threads=args.threads,
                            fn_idx_out=str(mmi))
    else:
        idx = mappy.Aligner(str(mmi), preset="sr", n_threads=args.threads)
    if not idx:
        raise SystemExit("minimap2 index failed")

    rows = []
    for name, d in mates.items():
        a = axes.get(name)
        if a is None or "1" not in d or "2" not in d:
            continue
        h1 = best_contig(idx, d["1"], args.min_mlen)
        h2 = best_contig(idx, d["2"], args.min_mlen)
        if not h1 or not h2:
            continue
        c1, c2 = COV.search(h1[0]), COV.search(h2[0])
        if not c1 or not c2:
            continue
        v1, v2 = float(c1.group(1)), float(c2.group(1))
        if v1 <= 0 or v2 <= 0:
            continue
        rows.append({"name": name, "ctg1": h1[0], "ctg2": h2[0],
                     "same_contig": h1[0] == h2[0],
                     "depth_log2_ratio_abs": abs(math.log2(v1 / v2)),
                     **{k: a[k] for k in a}})

    fields = [k for k in ("cov_step", "cov_spread", "mate_disagree", "zero_frac", "log_cov")
              if rows and k in rows[0]]
    rng = random.Random(args.seed)

    def corr(sub, field):
        xs = [r[field] for r in sub]
        ys = [r["depth_log2_ratio_abs"] for r in sub]
        if len(sub) < 10:
            return None
        rho = spearman(xs, ys)
        perm = ys[:]
        rng.shuffle(perm)
        return {"n": len(sub), "spearman": round(rho, 4),
                "ci95": fisher_ci(rho, len(sub)),
                "shuffled": round(spearman(xs, perm), 4)}

    diff = [r for r in rows if not r["same_contig"]]
    same = [r for r in rows if r["same_contig"]]
    out = {
        "seed": args.seed,
        "axes_rows": len(axes),
        "pairs_both_mates_mapped": len(rows),
        "pairs_mates_on_different_contigs": len(diff),
        "depth_ratio_summary": {
            "all_median": round(st.median([r["depth_log2_ratio_abs"] for r in rows]), 4)
            if rows else None,
            "different_contig_median":
                round(st.median([r["depth_log2_ratio_abs"] for r in diff]), 4)
                if diff else None,
        },
        "correlations": {
            "all_mapped_pairs": {f: corr(rows, f) for f in fields},
            "mates_on_different_contigs": {f: corr(diff, f) for f in fields},
            "mates_on_same_contig": {f: corr(same, f) for f in fields},
        },
        "interpretation_rules": [
            "the prediction passes only if cov_step's interval is clear of zero AND larger "
            "than the same correlation for the control fields; a correlation shared by "
            "every field is a property of the mapping, not of cov_step",
            "the shuffled column is the null for each cell",
            "only pairs where BOTH mates map and both contigs carry a SPAdes cov_ field "
            "enter, which selects towards the assembled part of the residue - the rarest "
            "molecules are absent from this test by construction",
            "mates_on_same_contig is DEGENERATE BY CONSTRUCTION: both mates share one "
            "contig so the depth ratio is identically zero, and every correlation there "
            "is exactly 0.000 because a constant has no rank order. It is printed as the "
            "check that it came out that way, not as a measurement",
            "all_mapped_pairs is therefore A MIXTURE of that constant group and the "
            "informative one, so a correlation there mostly measures whether the mates "
            "landed on different contigs AT ALL, not how different the depths were. The "
            "test of the prediction is the mates_on_different_contigs scope alone",
        ],
        "wall_seconds": time.perf_counter() - t0,
    }
    (args.out / f"cov_step_join_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/cov_step_join.py", **provenance(__file__),
        "axes_source": str(args.axes),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    print(f"both mates mapped with depths: {len(rows)} of {len(axes)} residue pairs "
          f"({len(diff)} on different contigs)")
    for scope, d in out["correlations"].items():
        print(f"== {scope}")
        for f, c in d.items():
            if c:
                print(f"   {f:15} rho {c['spearman']:+.3f} ci {c['ci95']} "
                      f"shuffled {c['shuffled']:+.3f}  n={c['n']}")


if __name__ == "__main__":
    main()
