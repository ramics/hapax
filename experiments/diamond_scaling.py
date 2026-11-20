"""How does a DIAMOND search scale with query count and database size? The measurement
every cost estimate for a bigger index has to rest on.

WHY THIS EXISTS. The protein ladder's per-seed runs each put 8,000 read pairs against four
databases in about half an hour, and it is tempting to turn that into "N hours for the full
UniRef50" by multiplying. That extrapolation is wrong in a way worth measuring rather than
arguing about: at 8,000 queries a DIAMOND run is dominated by scanning the database, which
is a fixed cost paid once, so per-read throughput at a production query count (a whole
sample's residue is ~2 million pairs) is much better than at a ladder query count. A single
point cannot tell the two regimes apart - the same mistake as reading two index sizes as a
curve - so this script measures both axes:

1. QUERY SCALING at fixed database: the slope in seconds per added read is the marginal
   cost that a production estimate needs; the intercept is the database scan that a bigger
   index inflates.
2. DATABASE SCALING at fixed query count: how the intercept grows with index bytes, which
   is what buying a bigger database actually costs per pass.

Both in DIAMOND's default mode and in --sensitive, because the mode is a factor of several
and the cost of an index is not a property of the index alone.

Fitted on a least-squares line through the query ladder, reported with the residuals, and
the extrapolation to a full-sample residue stated as what it is: a linear extrapolation
past the largest point measured, which is exactly the thing this repo keeps catching.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import git_sha

ROOT = Path(__file__).resolve().parents[1]


def take(src: Path, n: int, dst: Path) -> int:
    """First n records. Order is the residue's own read order, which is flow-cell order and
    carries no relationship to sequence content, so a prefix is a fair sample here."""
    kept = 0
    with open(src) as h, open(dst, "w") as o:
        for line in h:
            if line.startswith(">"):
                if kept >= n:
                    break
                kept += 1
            o.write(line)
    return kept


def run(diamond: Path, db: Path, q: Path, threads: int, mode: str, work: Path) -> float:
    cmd = [str(diamond), "blastx", "-d", str(db), "-q", str(q), "-o",
           str(work / "hits.tsv"), "-p", str(threads), "-k", "1", "--quiet",
           "-f", "6", "qseqid", "sseqid", "pident", "bitscore"]
    if mode:
        cmd.append(mode)
    t = time.perf_counter()
    subprocess.run(cmd, check=True)
    return time.perf_counter() - t


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den if den else 0.0
    inter = my - slope * mx
    resid = [y - (inter + slope * x) for x, y in zip(xs, ys)]
    return slope, inter, resid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reads", type=Path,
                    default=Path("/tmp/claude-0/wk/recut/residue_s1.fa"))
    ap.add_argument("--query-ladder", type=int, nargs="+",
                    default=[2000, 8000, 32000])
    ap.add_argument("--dbs", type=Path, nargs="+",
                    default=[ROOT / "data/prot/uniref50_02pct.dmnd",
                             ROOT / "data/prot/uniref50_10pct.dmnd"])
    ap.add_argument("--modes", nargs="+", default=["", "--sensitive"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--target-pairs", type=int, default=2_000_000,
                    help="a whole sample's residue, for the extrapolation line")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/dscale"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/diamond_scaling")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    queries = {}
    for n in args.query_ladder:
        p = args.work / f"q{n}.fa"
        queries[n] = (p, take(args.reads, n, p))

    results = []
    for db in args.dbs:
        for mode in args.modes:
            for n in args.query_ladder:
                p, got = queries[n]
                secs = run(args.diamond_bin, db, p, args.threads, mode, args.work)
                results.append({"db": db.name, "db_bytes": db.stat().st_size,
                                "mode": mode or "default", "queries": got,
                                "seconds": secs,
                                "seconds_per_1k_reads": secs / got * 1000})
                print(f"{db.name:24} {mode or 'default':11} {got:6} reads "
                      f"{secs:8.1f}s  {secs / got * 1000:6.2f}s/1k", flush=True)

    fits = {}
    for db in args.dbs:
        for mode in args.modes:
            key = f"{db.name}|{mode or 'default'}"
            pts = [(r["queries"], r["seconds"]) for r in results
                   if r["db"] == db.name and r["mode"] == (mode or "default")]
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            slope, inter, resid = fit(xs, ys)
            fits[key] = {
                "marginal_seconds_per_read": slope,
                "fixed_seconds_db_scan": inter,
                "residuals": resid,
                "extrapolated_hours_for_target": (inter + slope * args.target_pairs) / 3600,
                "largest_measured_queries": max(xs),
                "extrapolation_factor": args.target_pairs / max(xs),
            }

    out = {"seed": args.seed, "threads": args.threads, "results": results, "fits": fits,
           "target_pairs": args.target_pairs,
           "caveats": [
               "the target line is a LINEAR EXTRAPOLATION far past the largest measured "
               "query count; it is a planning number, not a measurement",
               "a database subset is not the same shape as the full database: more "
               "distinct seeds per byte would make the full index slower per byte than "
               "these subsets",
               "wall time on a shared 4-core container; an 8 vCPU instance is not simply "
               "twice as fast",
           ],
           "wall_seconds": time.perf_counter() - t0}
    (args.out / f"scaling_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/diamond_scaling.py", "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) and v
                           and isinstance(v[0], Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    for k, f in fits.items():
        print(f"{k}: marginal {f['marginal_seconds_per_read'] * 1000:.2f}s/1k reads, "
              f"fixed scan {f['fixed_seconds_db_scan']:.0f}s, "
              f"{f['extrapolated_hours_for_target']:.1f}h for "
              f"{args.target_pairs:,} (x{f['extrapolation_factor']:.0f} extrapolation)")


if __name__ == "__main__":
    main()
