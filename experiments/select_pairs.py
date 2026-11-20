"""Choose several distinct real genome pairs per divergence rung.

Why this exists. The archived two-lineage ladder
(experiments/codon_assembler/experiments/two_lineage_rung.py) has one genome pair per
rung, so divergence and pair identity are collinear across the whole ladder: any
trend in collapse could be a property of seven particular pairs rather than of
divergence. Two or three distinct pairs per rung break that, so between-pair spread
within a rung can be compared with the trend across rungs.

Pairs come from the full pairwise matrix over the same 200 deposited rabies genomes the
archived matrix used (its accession list is committed in
experiments/codon_assembler/results/abundance_sweep/divergence_pairs.json), so the
selection is reproducible and still never mutates a genome to order: a dialled-in
divergence would be built by the same substitution process the assembler's k-mer
matching is being tested on, which makes the rung circular.

Selection rule, per rung: keep the archived pair (continuity with every committed
number), then add the closest pairs within the tolerance that share no genome with any
pair already chosen anywhere on the ladder, so no single genome can carry a rung or be
reused between rungs. A rung that cannot be filled is reported short, not padded.

Deterministic given the accession list, so it takes no seed. The full matrix is written
out so later selections are a lookup, not a recomputation.

Usage: python experiments/select_pairs.py --fasta data/pairs/rabies_200.fa
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path

import edlib

ROOT = Path(__file__).resolve().parents[1]
ARCHIVED = ROOT / "experiments/codon_assembler/results/abundance_sweep/divergence_pairs.json"
# (tag, a, b, divergence) exactly as the archived ladder ran them
ARCHIVED_PAIRS = [
    ("div02", "PZ227608.1", "PZ227621.1", 0.0198),
    ("div04", "PV491557.1", "PZ324370.1", 0.0381),
    ("div045", "MG458312.1", "PZ227553.1", 0.0450),
    ("div055", "PZ227534.1", "PZ227561.1", 0.0550),
    ("div065", "AB009663.2", "PZ227534.1", 0.0652),
    ("div07", "PZ227608.1", "AB044824.1", 0.0729),
    ("div10", "NC_001542.1", "AB009663.2", 0.0965),
]
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(accs: list[str], out: Path) -> None:
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for i in range(0, len(accs), 100):
            r = subprocess.run(["curl", "-sS", "-G", EFETCH, "--data-urlencode", "db=nuccore",
                                "--data-urlencode", "id=" + ",".join(accs[i:i + 100]),
                                "--data-urlencode", "rettype=fasta", "--data-urlencode", "retmode=text"],
                               capture_output=True, text=True, check=True)
            fh.write(r.stdout)
            time.sleep(0.5)


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.open():
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        out[name] = "".join(buf).upper()
    return out


SEQS: dict[str, str] = {}


def _init(seqs):
    global SEQS
    SEQS = seqs


def _dist(pair):
    a, b = pair
    d = edlib.align(SEQS[a], SEQS[b], task="distance")["editDistance"]
    # same normalisation as the archived find_divergence_pairs.py
    return a, b, d / max(len(SEQS[a]), len(SEQS[b]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fasta", type=Path, default=ROOT / "data/pairs/rabies_200.fa")
    ap.add_argument("--per-rung", type=int, default=3)
    ap.add_argument("--tol", type=float, default=0.003, help="max |divergence - rung target|")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, default=ROOT / "results/pair_ladder")
    args = ap.parse_args()

    accs = json.loads(ARCHIVED.read_text())["accessions"]
    fetch(accs, args.fasta)
    seqs = {k: v for k, v in read_fasta(args.fasta).items() if k in set(accs)}
    missing = sorted(set(accs) - set(seqs))
    print(f"{len(seqs)} of {len(accs)} genomes fetched; missing {missing}", flush=True)

    t0 = time.perf_counter()
    combos = list(itertools.combinations(sorted(seqs), 2))
    with Pool(args.threads, initializer=_init, initargs=(seqs,)) as pool:
        rows = pool.map(_dist, combos, chunksize=50)
    secs = time.perf_counter() - t0
    print(f"{len(rows)} pairs in {secs:.0f}s", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out / "matrix.tsv.gz", "wt") as fh:
        fh.write("a\tb\tdivergence\n")
        for a, b, d in rows:
            fh.write(f"{a}\t{b}\t{d:.5f}\n")
    lookup = {frozenset((a, b)): d for a, b, d in rows}

    used: set[str] = set()
    for _, a, b, _ in ARCHIVED_PAIRS:
        used |= {a, b}
    cands_by = {tag: sorted((abs(d - t), d, x, y) for x, y, d in rows if abs(d - t) <= args.tol)
                for tag, _, _, t in ARCHIVED_PAIRS}
    rungs = {}
    # Scarcest rungs pick first, so the dense ones do not use up the genomes a sparse
    # rung needs. A rung still short after the strict pass may reuse a genome from
    # ANOTHER rung (flagged), never one already in its own rung.
    order = sorted(ARCHIVED_PAIRS, key=lambda r: len(cands_by[r[0]]))
    for strict in (True, False):
        for tag, a, b, target in order:
            rung = rungs.setdefault(tag, [{"tag": tag, "a": a, "b": b,
                                           "divergence": round(lookup.get(frozenset((a, b)), target), 4),
                                           "archived": True, "reuses_genome": False}])
            mine = {g for p in rung for g in (p["a"], p["b"])}
            for _, d, x, y in cands_by[tag]:
                if len(rung) >= args.per_rung:
                    break
                if x in mine or y in mine or (strict and (x in used or y in used)):
                    continue
                rung.append({"tag": tag, "a": x, "b": y, "divergence": round(d, 4),
                             "archived": False, "reuses_genome": bool({x, y} & used)})
                used |= {x, y}
                mine |= {x, y}
    ladder = []
    for tag, a, b, target in ARCHIVED_PAIRS:
        rung = rungs[tag]
        ladder.append({"tag": tag, "target": target, "n_candidates_in_tol": len(cands_by[tag]),
                       "short": len(rung) < args.per_rung, "pairs": rung})
        print(tag, target, [(p["a"], p["b"], p["divergence"], p["reuses_genome"]) for p in rung],
              f"({len(cands_by[tag])} in tolerance)", flush=True)

    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        sha = "unknown"
    (args.out / "ladder.json").write_text(json.dumps(ladder, indent=2))
    (args.out / "manifest.json").write_text(json.dumps({
        "script": "experiments/select_pairs.py", "git_sha": sha, "genomes": len(seqs),
        "missing": missing, "pairs": len(rows), "per_rung": args.per_rung, "tol": args.tol,
        "seconds": round(secs, 1), "fasta": str(args.fasta),
    }, indent=2))


if __name__ == "__main__":
    main()
