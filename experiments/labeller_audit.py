"""Audit of the per-window lineage labeller, and its strand-aware replacement.

Why this exists. Every mosaic count, lineage-switch count and join observable in
the two-lineage rabies study (experiments/codon_assembler/results/abundance_sweep)
comes from `lineage_windows()` in `two_lineage_rung.py`. That labeller reported a
mosaic in 13 of 42 single-lineage controls, which is impossible by construction.
Any flag trained on those labels inherits the error, so the labeller is audited
before anything is built on it.

The suspected cause: each 300 nt window is aligned to both lineages with edlib in
the forward orientation only. Assemblers emit contigs on either strand, and for a
reverse-strand window both edit distances are near-random, so the a/b call is a coin
flip and any contig of two or more windows can come out "mosaic".

What is measured, on the real genome pairs of the committed ladder (never mutated
genomes, which CLAUDE.md forbids for divergence rungs):
  * false-mosaic rate on pure single-lineage contigs, forward and reverse strand, at
    0, 0.5 and 1% added consensus error;
  * on synthetic one-join mosaics (lineage A then lineage B, joined at homologous
    coordinates), the fraction detected as mosaic and the join-position error.
Each is reported for the old labeller and the strand-aware one.

The strand-aware labeller scores each window against both lineages on both strands
and keeps the better strand per lineage. Calling per window rather than orienting the
whole contig once is deliberate: a chimaera can also be an inversion, and a per-contig
orientation would hide it.

Synthetic mosaics test sensitivity only in the easy case (one clean join, exact
sequence either side). They bound the labeller from above; they do not stand in for
real assembler mosaics.

Usage: python experiments/labeller_audit.py --genomes rabies_pairs.fa --seed 1 --n-seeds 5
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import random
import subprocess
import time
from pathlib import Path

import edlib

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "experiments" / "codon_assembler" / "experiments" / "two_lineage_rung.py"


def load_old():
    """Import the archived module for its PAIRS list (the genome pairs that ran)."""
    spec = importlib.util.spec_from_file_location("two_lineage_rung", OLD)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def lineage_windows_forward_only(contig: str, a: str, b: str, window: int = 300):
    """Verbatim copy of the labeller that produced the archived rows (two_lineage_rung.py
    up to commit d7b7783). Copied rather than imported because the archived file is being
    fixed in place, and this audit must keep testing what actually ran."""
    calls = []
    for i in range(0, len(contig) - window + 1, window):
        w = contig[i:i + window]
        da = edlib.align(w, a, mode="HW", task="distance")["editDistance"]
        db = edlib.align(w, b, mode="HW", task="distance")["editDistance"]
        if da < db:
            calls.append("a")
        elif db < da:
            calls.append("b")
    switches = sum(1 for x, y in zip(calls, calls[1:]) if x != y)
    return calls.count("a"), calls.count("b"), switches, calls


COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def lineage_windows_stranded(contig: str, a: str, b: str, window: int = 300):
    """Strand-aware per-window call, same return shape as the archived labeller."""
    calls = []
    for i in range(0, len(contig) - window + 1, window):
        w = contig[i:i + window]
        wr = rc(w)
        da = min(edlib.align(x, a, mode="HW", task="distance")["editDistance"] for x in (w, wr))
        db = min(edlib.align(x, b, mode="HW", task="distance")["editDistance"] for x in (w, wr))
        if da < db:
            calls.append("a")
        elif db < da:
            calls.append("b")
    switches = sum(1 for x, y in zip(calls, calls[1:]) if x != y)
    return calls.count("a"), calls.count("b"), switches, calls


def read_fasta(path: Path) -> dict[str, str]:
    seqs, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            seqs[key] = []
        elif key:
            seqs[key].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def add_errors(s: str, rate: float, rng: random.Random) -> str:
    if rate == 0:
        return s
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([c for c in "ACGT" if c != out[i]])
    return "".join(out)


def map_to_b(a: str, b: str, pos: int, flank: int = 150) -> int | None:
    """Homologous coordinate in b of position pos in a, via the left flank's end."""
    q = a[pos - flank:pos]
    r = edlib.align(q, b, mode="HW", task="locations")
    if not r["locations"] or r["editDistance"] > 0.2 * flank:
        return None
    return r["locations"][0][1] + 1


def is_mosaic(res) -> bool:
    wa, wb, sw, _ = res
    return bool(wa and wb and sw)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genomes", type=Path, required=True,
                    help="FASTA holding every accession in the archived PAIRS list")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--per-cell", type=int, default=100, help="contigs per (rung, condition, seed)")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "labeller_audit")
    args = ap.parse_args()

    old = load_old()
    genomes = read_fasta(args.genomes)
    labellers = {"archived": lineage_windows_forward_only, "stranded": lineage_windows_stranded}
    rows = []
    t0 = time.time()
    for seed in range(args.seed, args.seed + args.n_seeds):
        rng = random.Random(seed)
        for tag, acc_a, acc_b, div in old.PAIRS:
            sa, sb = genomes[acc_a], genomes[acc_b]
            # pure lineage-A contigs
            for strand in ("fwd", "rev"):
                for err in (0.0, 0.005, 0.01):
                    contigs = []
                    for _ in range(args.per_cell):
                        length = rng.choice([700, 1000, 1500, 3000])
                        i = rng.randrange(0, len(sa) - length)
                        c = add_errors(sa[i:i + length], err, rng)
                        contigs.append(rc(c) if strand == "rev" else c)
                    for name, fn in labellers.items():
                        fm = sum(is_mosaic(fn(c, sa, sb)) for c in contigs)
                        rows.append({"seed": seed, "pair": tag, "divergence": div, "kind": "pure",
                                     "strand": strand, "error": err, "labeller": name,
                                     "n": len(contigs), "called_mosaic": fm})
            # one-join mosaics, A then B at homologous coordinates
            for strand in ("fwd", "rev"):
                cases = []
                while len(cases) < args.per_cell:
                    half = rng.choice([450, 750, 1500])
                    i = rng.randrange(200, len(sa) - 2 * half - 200)
                    j = map_to_b(sa, sb, i + half)
                    if j is None or j + half > len(sb):
                        continue
                    c = sa[i:i + half] + sb[j:j + half]
                    cases.append((rc(c) if strand == "rev" else c, half))
                for name, fn in labellers.items():
                    hit, errs = 0, []
                    for c, half in cases:
                        res = fn(c, sa, sb)
                        if is_mosaic(res):
                            hit += 1
                            calls = res[3]
                            joins = [(k + 1) * 300 for k, (x, y) in enumerate(zip(calls, calls[1:])) if x != y]
                            # the archived convention: join at the first switching window boundary
                            errs.append(min(abs(p - half) for p in joins))
                    rows.append({"seed": seed, "pair": tag, "divergence": div, "kind": "mosaic",
                                 "strand": strand, "error": 0.0, "labeller": name,
                                 "n": len(cases), "called_mosaic": hit,
                                 "join_err_median": sorted(errs)[len(errs) // 2] if errs else None})
    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "rows.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
    except Exception:
        sha, dirty = "unknown", None
    (args.out / "manifest.json").write_text(json.dumps({
        "script": "experiments/labeller_audit.py", "git_sha": sha, "dirty": dirty,
        "seeds": list(range(args.seed, args.seed + args.n_seeds)), "per_cell": args.per_cell,
        "pairs": old.PAIRS, "window": old.WINDOW, "genomes": str(args.genomes),
        "edlib": getattr(edlib, "__version__", "unknown"), "python": platform.python_version(),
        "seconds": round(time.time() - t0, 1),
    }, indent=2))
    summarise(rows, args.out)


def summarise(rows, out: Path) -> None:
    """Pooled over seeds, with the per-seed range as spread."""
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        g[(r["kind"], r["strand"], r["error"], r["labeller"])].append(r["called_mosaic"] / r["n"])
    lines = ["| contigs | strand | added error | labeller | called mosaic, mean | range over rungs x seeds |",
             "|---|---|---|---|---|---|"]
    for (kind, strand, err, lab), v in sorted(g.items()):
        lines.append(f"| {kind} | {strand} | {err:.1%} | {lab} | {sum(v)/len(v):.1%} | {min(v):.1%}-{max(v):.1%} |")
    je = defaultdict(list)
    for r in rows:
        if r["kind"] == "mosaic" and r.get("join_err_median") is not None:
            je[(r["strand"], r["labeller"])].append(r["join_err_median"])
    lines += ["", "| mosaic strand | labeller | median join-position error (nt), median over cells |", "|---|---|---|"]
    for (strand, lab), v in sorted(je.items()):
        v = sorted(v)
        lines.append(f"| {strand} | {lab} | {v[len(v)//2]} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
