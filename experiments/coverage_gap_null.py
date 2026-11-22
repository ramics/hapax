#!/usr/bin/env python3
"""The coverage gap against its dinucleotide null, on the library with ground truth.

WHY THIS SCRIPT EXISTS, which is the only interesting thing about it. `results/
holdout_cascade/` reported a field - the gap between two definitions of how much of a
contig known proteins cover - that read +36.6 points on the withheld genome and 0.0 at
median, mean AND max across the 34 contigs of already-referenced genomes in the same
library. It was reported with its one missing control named: a union is the single
statistic that a pile of weak chance HSPs inflates without limit, because each one
contributes interval, and no shuffled arm had been run.

The gene-level thread then ran that arm on Berlin and the field FAILED, in the worst
available way: the null is not merely comparable, it is HIGHER. 14.47% of shuffled contigs
clearing the same significance bar showed a gap above 30 points against 5.81% of real
ones. The mechanism is exactly the one named when the control was requested.

So this measures the same thing on THIS library, where the positive control is a real
genome at a known distance rather than a population statistic, and tests the one
refinement proposed while that kill still stood: **a gap built only from long HSPs**. The
withheld genome's three genes align over 620, 493 and 135 aa, so a floor high enough to
exclude chance alignments also discards one of its real genes, and the sweep measures that
tension rather than assuming it away.

The test is still worth having now that the kill is gone, for two reasons. It fails on its
own terms rather than through the null - a floor moves the candidate from the 98.68th to the
99.27th percentile while discarding 96% of the real pile, which is shrinking the measurement
rather than sharpening a contrast. And with an empty null on this library, the sweep against
the REAL distribution is the only discriminating thing these data can offer: where the null
holds nothing, "the candidate is not chance" is as far as it goes, and the field's status is
UNTESTED rather than cleared.

Reported either way, and the thresholds are swept rather than chosen, because a floor
picked after seeing where one candidate sits is a threshold fitted to a positive set of
size one - which is the mistake the first version of this field already made once.

The shuffle is imported from `holdout_run.py` rather than rewritten: that one is the fixed
Altschul-Erikson walk that returns None on failure instead of silently handing back the
input, after a version that did the latter left 70.4% of a shuffled arm unshuffled.
Failures are counted and reported here for the same reason.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import random
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from holdout_run import dinuc_shuffle  # noqa: E402  the fixed walk, not a second one

ROOT = Path(__file__).resolve().parents[1]
SIGNIFICANCE_EVALUE = 1e-5
# Swept, not chosen. 0 is the field as originally reported.
AA_FLOORS = [0, 50, 100, 150, 200, 300, 400]


def log(*a) -> None:
    print(*a, file=sys.stderr, flush=True)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def diamond(binary: Path, query: Path, dbs: list[Path], out: Path, threads: int) -> Path:
    """One tabular DIAMOND pass per database, concatenated.

    Pooled across the viral and Swiss-Prot indexes because the cascade's protein tier is
    both of them, and because the thread this is answering pooled the same two: a gap
    computed over one index and compared with a gap computed over two is two instruments.
    """
    with open(out, "w") as w:
        for db in dbs:
            tmp = out.with_suffix(f".{Path(db).stem}.tsv")
            subprocess.run(
                [str(binary), "blastx", "--db", str(db), "--query", str(query),
                 "--out", str(tmp), "--threads", str(threads), "--outfmt", "6",
                 "qseqid", "sseqid", "pident", "length", "qstart", "qend",
                 "evalue", "bitscore", "qlen", "slen",
                 "--max-target-seqs", "5", "--evalue", "10", "--sensitive", "--quiet"],
                check=True)
            w.write(open(tmp).read())
    return out


def gaps(tsv: Path, aa_floor: int) -> dict[str, dict]:
    """Per contig: per-subject union, any-subject union, and the gap between them.

    `aa_floor` drops every HSP shorter than it BEFORE either union is formed, so the two
    definitions always see the same HSPs and the gap stays a property of how those HSPs
    are distributed across subjects rather than of which filter each side got.
    """
    per: dict[str, dict[str, list]] = collections.defaultdict(
        lambda: collections.defaultdict(list))
    allv: dict[str, list] = collections.defaultdict(list)
    qlen: dict[str, int] = {}
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 10 or float(f[6]) >= SIGNIFICANCE_EVALUE or int(f[3]) < aa_floor:
            continue
        q, s = f[0], f[1]
        a, b = sorted((int(f[4]), int(f[5])))
        per[q][s].append((a, b))
        allv[q].append((a, b))
        qlen[q] = int(f[8])

    def union(spans: list, L: int) -> float:
        spans = sorted(spans)
        m: list[list[int]] = []
        for a, b in spans:
            if m and a <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        return sum(b - a + 1 for a, b in m) / L if L else 0.0

    out = {}
    for q in allv:
        L = qlen[q]
        p = max(union(v, L) for v in per[q].values())
        a = union(allv[q], L)
        out[q] = {"per_subject": round(p, 4), "any_subject": round(a, 4),
                  "gap_points": round(100 * (a - p), 2), "subjects": len(per[q]),
                  "length": L}
    return out


def census(tsv: Path) -> dict:
    """How many contigs the arm returned ANY hit for, and how many a SIGNIFICANT one.

    Separated from the gap statistics because an arm with no significant hits and an arm
    that was never searched both produce an empty distribution, and this repo has now
    found three bugs in one day whose whole content was that confusion. The numbers below
    are what distinguish "the null is clean" from "the null did not run": the shuffled arm
    here returns 9,608 HSPs over 5,335 contigs and none of them significant, which is a
    measurement. An arm reporting 0 and 0 is a plumbing failure.
    """
    q_any: set[str] = set()
    q_sig: set[str] = set()
    best = None
    n = 0
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        n += 1
        q_any.add(f[0])
        e = float(f[6])
        if e < SIGNIFICANCE_EVALUE:
            q_sig.add(f[0])
        if best is None or e < best:
            best = e
    return {"hsps": n, "contigs_with_any_hit": len(q_any),
            "contigs_with_significant_hit": len(q_sig),
            "best_evalue": best}


def describe(vals: list[float]) -> dict:
    if not vals:
        return {"n": 0, "state": "no contig in this arm has a significant hit, so the "
                                 "arm has no gap distribution - this is a reading, not a "
                                 "missing measurement; see the census block"}
    v = sorted(vals)
    def q(p: float) -> float:
        return round(v[min(len(v) - 1, int(p * len(v)))], 2)
    return {"n": len(v), "median": round(st.median(v), 2), "p90": q(0.90),
            "p99": q(0.99), "max": round(v[-1], 2),
            "frac_above_10": round(sum(1 for x in v if x > 10) / len(v), 4),
            "frac_above_30": round(sum(1 for x in v if x > 30) / len(v), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=Path, required=True,
                    help="the identified, non-ribosomal pile the cascade removed")
    ap.add_argument("--candidate", required=True,
                    help="contig name of the withheld genome, the positive control")
    ap.add_argument("--leak-prefix", default=None,
                    help="optional contig name prefix for the near-reference arm")
    ap.add_argument("--dmnd", type=Path, nargs="+", required=True)
    ap.add_argument("--diamond-bin", type=Path, required=True)
    ap.add_argument("--min-contig", type=int, default=300)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--reuse-hits", action="store_true",
                    help="use the hit tables already in --work instead of "
                         "searching again. The searches are deterministic given "
                         "the same query file, index and seed, so this is a cache "
                         "and not a shortcut; the manifest records that it was used")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    t0 = time.time()

    seqs = {n: s for n, s in read_fasta(a.contigs).items() if len(s) >= a.min_contig}
    if a.candidate not in seqs:
        raise SystemExit(f"candidate {a.candidate} not in {a.contigs}")
    log(f"{len(seqs):,} contigs >= {a.min_contig} nt")

    real_fa = a.work / "real.fa"
    real_fa.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    shuf_fa = a.work / "shuffled.fa"
    failed = 0
    with open(shuf_fa, "w") as w:
        for n, s in seqs.items():
            sh = dinuc_shuffle(s, rng)
            if sh is None:
                failed += 1
                continue
            w.write(f">{n}\n{sh}\n")
    log(f"shuffle failed on {failed} contigs ({100 * failed / len(seqs):.2f}%), "
        f"dropped rather than returned unshuffled")

    hits = {}
    for arm, fa in (("real", real_fa), ("shuffled", shuf_fa)):
        out = a.work / f"hits_{arm}.tsv"
        if a.reuse_hits and out.exists() and out.stat().st_size > 0:
            log(f"  {arm}: reusing {out}")
            hits[arm] = out
            continue
        t = time.time()
        hits[arm] = diamond(a.diamond_bin, fa, a.dmnd, out, a.threads)
        log(f"  {arm}: {time.time() - t:.0f}s")

    # Read the census BEFORE the gaps, because it is what says whether an empty gap
    # distribution is a clean null or a search that never happened.
    cen = {arm: census(hits[arm]) for arm in hits}
    for arm, c in cen.items():
        log(f"  census {arm}: {c['hsps']:,} HSPs over "
            f"{c['contigs_with_any_hit']:,} contigs, "
            f"{c['contigs_with_significant_hit']:,} significant, "
            f"best E={c['best_evalue']:.3g}")

    result = {"census": cen, "arms": {}}
    for floor in AA_FLOORS:
        g = {arm: gaps(hits[arm], floor) for arm in hits}
        row = {}
        for arm in g:
            row[arm] = describe([v["gap_points"] for v in g[arm].values()])
        c = g["real"].get(a.candidate)
        if c:
            sh = sorted(v["gap_points"] for v in g["shuffled"].values())
            above = sum(1 for x in sh if x >= c["gap_points"])
            row["candidate"] = dict(c)
            rl = [v["gap_points"] for v in g["real"].values()]
            at_or_above = sum(1 for x in rl if x >= c["gap_points"])
            row["candidate"]["real_at_or_above_candidate"] = at_or_above
            row["candidate"]["real_percentile"] = (
                round(100 * sum(1 for x in rl if x <= c["gap_points"]) / len(rl), 2)
                if rl else None)
            # The only honest reading of one positive against a null: where in the null
            # it sits. A candidate at the null's p90 is not a detection.
            row["candidate"]["shuffled_at_or_above_candidate"] = above
            row["candidate"]["shuffled_fraction_at_or_above"] = (
                round(above / len(sh), 4) if sh else None)
        if a.leak_prefix:
            lk = [v["gap_points"] for n, v in g["real"].items()
                  if n.startswith(a.leak_prefix)]
            row["leak_arm"] = describe(lk)
        result["arms"][str(floor)] = row
        log(f"floor {floor:>3} aa: real median "
            f"{row['real'].get('median')} p90 {row['real'].get('p90')} | "
            f"shuffled median {row['shuffled'].get('median')} "
            f"p90 {row['shuffled'].get('p90')} | candidate "
            f"{row.get('candidate', {}).get('gap_points')} "
            f"(null frac at or above "
            f"{row.get('candidate', {}).get('shuffled_fraction_at_or_above')})")

    man = {
        "seed": a.seed,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "git_sha": subprocess.run(["git", "describe", "--always", "--dirty"],
                                  capture_output=True, text=True,
                                  cwd=ROOT).stdout.strip(),
        "contigs": str(a.contigs), "min_contig": a.min_contig,
        "candidate": a.candidate, "leak_prefix": a.leak_prefix,
        "dmnd": [str(d) for d in a.dmnd], "aa_floors": AA_FLOORS,
        "significance_evalue": SIGNIFICANCE_EVALUE,
        "contigs_scored": len(seqs), "shuffle_failures": failed,
        # Whether the hit tables were searched in this run or taken from --work. The
        # searches are deterministic given the same query file and index, but a reader
        # cannot check that from the output, so the output says which happened.
        "reused_hits": bool(a.reuse_hits),
        "seconds": round(time.time() - t0, 1),
        "reading": "with an empty shuffled null on this library the field cannot be "
                   "cleared here, only shown not to be chance; the floor sweep is "
                   "against the real distribution, and the floors are swept, not chosen",
        "note_on_the_withdrawn_kill":
            "the Berlin kill this script was written to test was WITHDRAWN: that null "
            "returned unshuffled input on dead-end walks, so every significant hit in it "
            "was real material compared against itself. This arm imports the fixed walk "
            "from holdout_run.py - 0 failures, 0 returned unchanged, 0 significant hits "
            "on 22,110 shuffled contigs. Second occurrence of that bug in this repo in "
            "two days; the shared module is experiments/shuffle_control.py.",
    }
    (a.out / f"coverage_gap_null_s{a.seed}.json").write_text(
        json.dumps(result, indent=1))
    (a.out / f"coverage_gap_null_manifest_s{a.seed}.json").write_text(
        json.dumps(man, indent=1))
    log(f"wrote {a.out}/coverage_gap_null_s{a.seed}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
