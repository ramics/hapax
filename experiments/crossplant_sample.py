"""On-instance phase 1 of the cross-plant pass: cut each library to a fixed depth, draw queries.

Per library: the first `--depth` spots as split FASTQ (the subject, kept for phase 2), the
5' base composition of each mate (to check the trim from data), and per seed a sample of
query pairs drawn from the next `--pool` spots, which are never in the subject. Drawing
queries from outside the subject keeps the own-library count free of the read itself and of
size-biased selection, so growth against the panel compares like with like. Queries are adapter-cut, R2-trimmed by the library's kit rule and truncated to `--qlen` per
mate so every query offers the same number of 31-mers whatever the read length. A
mononucleotide shuffle of part of each seed's draw rides along as the null.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import subprocess
import sys
from pathlib import Path

ADAPTER = "AGATCGGAAGAGC"
COMP_POS = 30


def run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def dump(acc: str, depth: int, work: Path) -> tuple[Path, Path]:
    """First `depth` spots, split. Remote access first; stage the object if that fails."""
    base = ["fastq-dump", "--split-files", "-X", str(depth), "--outdir", str(work)]
    try:
        run(base + [acc])
    except subprocess.CalledProcessError:
        sra = work / f"{acc}.sra"
        for extra in (["--no-sign-request"], []):
            try:
                run(["aws", "s3", "cp", f"s3://sra-pub-run-odp/sra/{acc}/{acc}", str(sra),
                     *extra, "--only-show-errors"])
                break
            except subprocess.CalledProcessError:
                continue
        run(base + [str(sra)])
        sra.unlink(missing_ok=True)
    r1, r2 = work / f"{acc}_1.fastq", work / f"{acc}_2.fastq"
    if not (r1.exists() and r2.exists()):
        sys.exit(f"{acc}: no two mate files")
    return r1, r2


def records(path: Path):
    with open(path) as h:
        while True:
            head = h.readline()
            if not head:
                return
            seq = h.readline().rstrip()
            h.readline()
            h.readline()
            yield seq


def clean(seq: str, trim: int, qlen: int) -> str | None:
    seq = seq[trim:]
    cut = seq.find(ADAPTER)
    if cut >= 0:
        seq = seq[:cut]
    seq = seq[:qlen]
    if len(seq) < qlen or seq.count("N") > 0:
        return None
    return seq


def shuffled(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--acc", required=True)
    ap.add_argument("--trim-r2", type=int, required=True)
    ap.add_argument("--depth", type=int, default=2_000_000)
    ap.add_argument("--pool", type=int, default=200_000)
    ap.add_argument("--per-seed", type=int, default=200)
    ap.add_argument("--null-per-seed", type=int, default=20)
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--qlen", type=int, default=90)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--subjects", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    acc = args.acc
    args.work.mkdir(parents=True, exist_ok=True)
    args.subjects.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    done = args.out / f"sample_{acc}.json"
    if done.exists():
        print(f"= {acc} already sampled")
        return

    r1, r2 = dump(acc, args.depth + args.pool, args.work)
    comp = {m: [[0] * 5 for _ in range(COMP_POS)] for m in ("r1", "r2")}
    lens = {"r1": 0, "r2": 0}
    pairs = []
    for s1, s2 in zip(records(r1), records(r2)):
        for m, s in (("r1", s1), ("r2", s2)):
            lens[m] += len(s)
            for i, b in enumerate(s[:COMP_POS]):
                comp[m][i]["ACGTN".find(b) if b in "ACGT" else 4] += 1
        pairs.append((s1, s2))
    n = len(pairs)
    if n < args.depth + args.pool:
        sys.exit(f"{acc}: only {n} spots, under the {args.depth + args.pool} floor")
    subject, pairs = pairs[: args.depth], pairs[args.depth:]

    rows = []
    for seed in (int(s) for s in args.seeds.split(",")):
        rng = random.Random(f"{acc}|{seed}")
        order = list(range(len(pairs)))
        rng.shuffle(order)
        got = 0
        for i in order:
            q1 = clean(pairs[i][0], 0, args.qlen)
            q2 = clean(pairs[i][1], args.trim_r2, args.qlen)
            if q1 is None or q2 is None:
                continue
            rows.append((f"lib|{acc}|{args.depth + i}|s{seed}", q1, q2))
            if got < args.null_per_seed:
                rows.append((f"shuffled-lib|{acc}|{i}|s{seed}", shuffled(q1, rng),
                             shuffled(q2, rng)))
            got += 1
            if got == args.per_seed:
                break
        if got < args.per_seed:
            sys.exit(f"{acc}: only {got} clean pairs for seed {seed}")

    with open(args.out / f"queries_{acc}.tsv", "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    for src, mate in ((r1, 1), (r2, 2)):
        with open(src) as i, gzip.open(args.subjects / f"{acc}_{mate}.fastq.gz", "wt",
                                       compresslevel=1) as o:
            for _ in range(4 * args.depth):
                o.write(i.readline())
        src.unlink()
    del subject
    done.write_text(json.dumps({"acc": acc, "spots": n, "subject_spots": args.depth,
                                "query_pool_spots": len(pairs), "trim_r2": args.trim_r2,
                                "mean_len": {m: round(v / n, 1) for m, v in lens.items()},
                                "queries": len(rows), "composition_5p": comp}) + "\n")
    print(f"+ {acc} sampled {n} spots, {len(rows)} query rows")


if __name__ == "__main__":
    main()
