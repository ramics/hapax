"""Build the query sets the big-index run searches. Runs here, not on the instance.

WHY THE QUERIES ARE BUILT HERE. Deciding which reads are the residue is the cascade's job
and it has already been done, five times, by `experiments/subtract_known.py`: the residue
and a matched sample of SUBTRACTED (identified) pairs are written out as FASTA for exactly
this kind of reuse. Rebuilding those sets on the instance would mean shipping SILVA,
Kraken2's hash, the sample's own contigs and the labelling code to a machine that cannot be
debugged, to recompute something already on disk. So the instance gets FASTA and nothing
else, and the only thing that has to be right about it is the names.

THE CONTROL IS NOT OPTIONAL. Every pile in this project is reported beside the already-
identified pile, because that is what caught two bugs in one evening: a test firing on 100%
of both populations, and a shuffled control identical to the real reads. A big-index pass
that searched only the residue would produce a column of identities with nothing to compare
them to - and the protein ladder already had to re-label such a column "not searched".
Here the control costs the same per read as the residue and makes every number a contrast.

NAMES CARRY THE CLASS, as `<class>|<read>` with the read name exactly as the cascade wrote
it (`SRR23973193.<n>/<mate>`). That is what lets the instance-side script stay a search
wrapper with no reference data, and it is what lets the residue NAME LISTS the sibling
thread asked for be derived afterwards in the container, where the agreed format can still
be fixed without another launch.

BOTH MATES IN THE BLOCK FILES, mate 1 only in the samples. A pair survives a rung only if
NEITHER mate is identified, so a mate-1-only pass would under-identify and the name lists
would be wrong in the dangerous direction. It is also the rule this repo learned the hard
way: mate 2 is 93 nt against mate 1's 107 after the non-template trim, so anything
mate-interchangeable has to be reported per mate, which needs both mates present. The pilot
and the nr subsample are distributions rather than survival sets, so they take mate 1 only
and cost half as much.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]]) -> int:
    with open(path, "w") as h:
        for name, seq in records:
            h.write(f">{name}\n{seq}\n")
    return len(records)


def mate1(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return [(n, s) for n, s in records if n.endswith("/1")]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="draws the subsamples. The blocks themselves are not sampled")
    ap.add_argument("--blocks", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--cascade-work", type=Path,
                    default=Path("/tmp/claude-0/wk/recut"),
                    help="where subtract_known.py wrote residue_s*.fa and subtracted_s*.fa")
    ap.add_argument("--pilot", type=int, default=2000,
                    help="sequences in the timed pilot. Big enough that the rate is not "
                         "all database scan, small enough that the pilot is minutes")
    ap.add_argument("--scan-probe", type=int, default=10,
                    help="a query set small enough that its search time IS the database "
                         "scan. BLAST streams the whole subject database per invocation, "
                         "so this plus the pilot are two query counts against one database "
                         "and the difference between them is the marginal cost per read - "
                         "the split that a single rate cannot show and that this project "
                         "has already once got wrong by reading one")
    ap.add_argument("--nr-sample", type=int, default=3000,
                    help="sequences per pile for the blastx-against-nr pass. nr is only "
                         "readable by blastx, which is far slower per read than DIAMOND, "
                         "and the question it answers is a comparison of identity "
                         "distributions - which a few thousand per pile settles at the "
                         "percentage resolution this project reports anyway")
    ap.add_argument("--per-block-sample", type=int, default=0,
                    help="sequences per pile per block, 0 for all of them. The first "
                         "launch asked for all 429,004 and the pilot measured a sensitive "
                         "blastn against core_nt at roughly a second per read on eight "
                         "cores, which is days. A right-sized sample is not a compromise "
                         "here: every proportion in this project is already quoted with an "
                         "interval, and a few thousand reads per pile per block resolves a "
                         "distribution to the percentage level the report uses")
    ap.add_argument("--mate1-only", action="store_true",
                    help="halves the cost. Legitimate for a DISTRIBUTION and wrong for a "
                         "survival set: a pair survives a rung only if neither mate is "
                         "identified, so a name list built from mate 1 alone invents "
                         "residue. Use it for identity distributions, not for the lists")
    ap.add_argument("--eyeball-rows", type=Path, default=None,
                    help="rows from an earlier pass. With it, an eyeball.fa is written: the "
                         "reads whose alignments are to be rendered against the reference "
                         "they matched, stratified by E-value rather than sampled at "
                         "random. Proportionally the pile is ~90% chance matches, which is "
                         "the true composition and the wrong file - the populations worth "
                         "looking at are the small ones")
    ap.add_argument("--eyeball-n", type=int, default=1000)
    ap.add_argument("--eyeball-reads", type=Path, default=None,
                    help="the FASTA those rows were computed from. Required in practice: "
                         "this script re-samples its blocks every time it runs, so the "
                         "reads it would otherwise draw from share only a handful of names "
                         "with rows from an earlier set - which is how a 45-read selection "
                         "silently became 4")
    ap.add_argument("--out", type=Path, default=Path("/tmp/claude-0/wk/bigindex_queries"))
    ap.add_argument("--manifest-out", type=Path,
                    default=ROOT / "results/subtract_known")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    t0 = time.perf_counter()
    files: dict[str, dict] = {}

    pooled_residue: list[tuple[str, str]] = []
    pooled_control: list[tuple[str, str]] = []
    for b in args.blocks:
        res = args.cascade_work / f"residue_s{b}.fa"
        ctl = args.cascade_work / f"subtracted_s{b}.fa"
        for p in (res, ctl):
            if not p.exists():
                raise SystemExit(f"{p} is missing; the cascade's work dir is the input")
        residue = [(f"residue|{n}", s) for n, s in read_fasta(res)]
        control = [(f"identified|{n}", s) for n, s in read_fasta(ctl)]
        if args.mate1_only:
            residue, control = mate1(residue), mate1(control)
        if args.per_block_sample:
            # Sampled at random rather than taken as a prefix: the cascade writes its
            # residue in the library's own read order, and while that is flow-cell order
            # rather than anything biological, a prefix of it is one region of one lane.
            rb = random.Random(args.seed * 1000 + b)
            rb.shuffle(residue)
            rb.shuffle(control)
            residue = residue[:args.per_block_sample]
            control = control[:args.per_block_sample]
        out = args.out / f"block_s{b}.fa"
        # Interleaved rather than concatenated: a search killed part way through then has a
        # usable sample of BOTH piles instead of all of one and none of the other, and on a
        # machine with a dead man's switch that is not a hypothetical.
        merged: list[tuple[str, str]] = []
        for i in range(max(len(residue), len(control))):
            if i < len(residue):
                merged.append(residue[i])
            if i < len(control):
                merged.append(control[i])
        n = write_fasta(out, merged)
        files[out.name] = {"sequences": n, "residue": len(residue),
                           "identified": len(control), "bytes": out.stat().st_size,
                           "sha256_16": digest(out)}
        print(f"{out.name}: {n} sequences ({len(residue)} residue, "
              f"{len(control)} identified)", flush=True)
        pooled_residue += mate1(residue)
        pooled_control += mate1(control)

    rng.shuffle(pooled_residue)
    rng.shuffle(pooled_control)

    # The pilot is residue only and mate 1 only: it exists to measure a rate, and the rate
    # is what the full pass is projected from, so it has to be drawn from the same kind of
    # sequence the full pass is mostly made of.
    pilot = pooled_residue[:args.pilot]
    p = args.out / "pilot.fa"
    files["pilot.fa"] = {"sequences": write_fasta(p, pilot), "bytes": p.stat().st_size,
                         "sha256_16": digest(p)}

    sp = args.out / "scan_probe.fa"
    files["scan_probe.fa"] = {
        "sequences": write_fasta(sp, pooled_residue[:args.scan_probe]),
        "bytes": sp.stat().st_size, "sha256_16": digest(sp)}

    nr = pooled_residue[args.pilot:args.pilot + args.nr_sample] + \
        pooled_control[:args.nr_sample]
    rng.shuffle(nr)
    q = args.out / "nr_sample.fa"
    files["nr_sample.fa"] = {"sequences": write_fasta(q, nr), "bytes": q.stat().st_size,
                             "sha256_16": digest(q)}
    print(f"pilot.fa: {len(pilot)} residue mate-1; nr_sample.fa: {len(nr)} "
          f"({args.nr_sample} per pile)", flush=True)

    if args.eyeball_rows and args.eyeball_rows.exists():
        # The selection rule lives in partial_alignment_eyeball.py, which renders the
        # result, so it is imported rather than restated: two copies of a stratification
        # drift, and the stratum labels in the file have to mean what the quotas meant.
        from partial_alignment_eyeball import pick as eyeball_pick
        rows = [json.loads(l) for l in open(args.eyeball_rows)]
        if args.eyeball_reads:
            pooled = {n.split("|", 1)[-1]: sq
                      for n, sq in read_fasta(args.eyeball_reads)}
        else:
            pooled = {n.split("|", 1)[-1]: sq
                      for n, sq in pooled_residue + pooled_control}
        covered = sum(1 for r in rows if r.get("read") in pooled)
        if covered < 0.9 * len(rows):
            raise SystemExit(
                f"only {covered} of {len(rows)} rows name a read in the eyeball read set; "
                f"pass --eyeball-reads pointing at the FASTA those rows came from")
        names, avail = eyeball_pick(rows, pooled, args.eyeball_n, 30, 1e-5, 0.05, rng)
        eb = args.out / "eyeball.fa"
        files["eyeball.fa"] = {
            "sequences": write_fasta(eb, [(f"residue|{n}", pooled[n]) for n in names]),
            "bytes": eb.stat().st_size, "sha256_16": digest(eb),
            "available_per_stratum": avail, "from_rows": str(args.eyeball_rows)}
        print(f"eyeball.fa: {len(names)} reads, strata available {avail}", flush=True)

    manifest = {
        "script": "experiments/bigindex_queries.py", **provenance(__file__),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "files": files,
        "total_sequences": sum(f["sequences"] for f in files.values()),
        "total_bytes": sum(f["bytes"] for f in files.values()),
        "note": "the block files carry both mates and both piles; pilot and nr_sample are "
                "mate 1 only. Class is the part of each name before the pipe.",
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }
    args.manifest_out.mkdir(parents=True, exist_ok=True)
    (args.manifest_out / "bigindex_queries_manifest.json").write_text(
        json.dumps(manifest, indent=2))
    print(f"{manifest['total_sequences']:,} sequences, "
          f"{manifest['total_bytes'] / 1e6:.1f} MB in {args.out}")


if __name__ == "__main__":
    main()
