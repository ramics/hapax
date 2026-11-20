"""One search against one big database, on the instance. The whole of what EC2 does.

WHY THIS IS THE ONLY THING THAT RUNS THERE. The box exists for one reason: the databases
are too big for this container, and they sit in `us-east-1` where reading them is free.
Everything else this project does with a search result - assigning classes, Wilson
intervals, the residue name lists, the removal curve - already runs here in minutes and is
already tested here. Shipping any of it to the instance would buy nothing and risk the one
thing that cannot be retried cheaply: `UserData` is the whole run, with no SSH and no second
chance, so the payload is kept to the one step that has to be there.

So this script searches and reports. It makes no decisions and applies no thresholds: it
writes one row per query with the best alignment's identity, aligned length and bitscore,
and the identity DISTRIBUTION is read off those rows afterwards. That is deliberate - every
single-axis threshold this project has tried has died, and the deliverable it keeps arriving
at is identity to the closest relative as a continuous number.

WHAT A ROW CAN AND CANNOT SAY. `pident` is the identity of the aligned part only, so a 40%
alignment at 95% identity and a full-length one at 95% are the same number in that column
and wildly different claims. `aln_frac` is beside it for exactly that reason, and both are
needed to answer requirement 1.iii ("is the no-match sequence really below ~80% identity"):
a read whose best alignment covers a third of itself has not been shown to be at that
identity over its length.

TIMING IS A RESULT HERE, NOT INSTRUMENTATION. `results/diamond_scaling/` could only measure
two database sizes, and the linear and the superlinear law it fits disagree by a factor of
fifteen at the far end. The seconds this script records on a real 100 GB index are the third
point that separates them, so they are written to their own file and uploaded even when the
search itself fails.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import time
from pathlib import Path

#: Tab-separated fields asked of both programs, in this order. Identical for BLAST+ and
#: DIAMOND so one parser reads both: the formats differ in what they CAN report, not in
#: these columns, and a shared parser is one fewer thing to get wrong on a box that cannot
#: be debugged interactively.
FIELDS = ("qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore")


def command(program: str, db: Path, query: Path, out: Path, threads: int,
            max_targets: int, sensitive: bool, blast_bin: Path,
            diamond_bin: Path, blastn_task: str = "blastn") -> list[str]:
    """The search, as the exact argument list that will be recorded in the manifest.

    `blastn` defaults to the plain task rather than megablast, matching
    `experiments/remote_blast_check.py`: megablast is tuned for close matches, and the
    question here is whether a relative exists *however far*, so megablast would answer a
    different question and flatter the cascade.

    `blastn_task` exists because the plain task's cost turned out to be the binding
    constraint: measured on this run's own pilot, a sensitive blastn of 107 nt reads against
    core_nt is about a second per read on eight cores, which is 550+ hours for a whole
    sample's residue and therefore not a step any production tool can contain. So the task
    becomes a measured axis rather than a fixed choice: `dc-megablast` is built for
    cross-species nucleotide matching and is far cheaper, and running both on the SAME reads
    says what the cheap task fails to find. That difference is the cost-versus-completeness
    curve this project exists to publish, applied to itself.

    DIAMOND's `--sensitive` is not a free switch - measured in this repo it costs 10-14x
    default - so it is a flag, and which was used goes in the manifest beside the timing.
    """
    fmt = list(FIELDS)
    if program == "blastn":
        return [str(blast_bin / "blastn"), "-task", blastn_task, "-db", str(db),
                "-query", str(query), "-out", str(out), "-num_threads", str(threads),
                "-max_target_seqs", str(max_targets), "-evalue", "10",
                "-outfmt", " ".join(["6", *fmt])]
    if program == "blastx":
        return [str(blast_bin / "blastx"), "-db", str(db), "-query", str(query),
                "-out", str(out), "-num_threads", str(threads),
                "-max_target_seqs", str(max_targets), "-evalue", "10",
                "-outfmt", " ".join(["6", *fmt])]
    if program == "diamond-blastx":
        cmd = [str(diamond_bin), "blastx", "-d", str(db), "-q", str(query),
               "-o", str(out), "-p", str(threads), "-k", str(max_targets),
               "--quiet", "-f", "6", *fmt]
        if sensitive:
            cmd.append("--sensitive")
        return cmd
    raise SystemExit(f"unknown program {program!r}")


def count_queries(path: Path) -> int:
    n = 0
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                n += 1
    return n


def best_rows(hits: Path, queries: Path) -> list[dict]:
    """One row per QUERY, including the queries that matched nothing.

    The no-hit rows are the point. A hit table alone cannot distinguish "searched and found
    nothing" from "never searched", and this project has already published a table where
    that distinction was lost and a column of 0.00% had to be re-labelled "not searched".
    """
    best: dict[str, dict] = {}
    with open(hits) as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(FIELDS):
                continue
            q, s, pid, ln, qlen, ev, bits = f[:7]
            try:
                row = {"sseqid": s, "pident": float(pid), "aln_len": int(ln),
                       "qlen": int(qlen), "evalue": float(ev), "bitscore": float(bits)}
            except ValueError:
                continue
            row["aln_frac"] = row["aln_len"] / row["qlen"] if row["qlen"] else 0.0
            prev = best.get(q)
            if prev is None or row["bitscore"] > prev["bitscore"]:
                best[q] = row

    out = []
    with open(queries) as h:
        for line in h:
            if not line.startswith(">"):
                continue
            name = line[1:].strip().split()[0]
            # The class is carried in the query name as `<class>|<read>`, written by
            # experiments/bigindex_queries.py, so the instance needs no labelling code and
            # no reference data to say which pile a row belongs to.
            cls, _, read = name.partition("|")
            row = {"query": name, "class": cls if read else "", "read": read or name,
                   "hit": name in best}
            row.update(best.get(name, {}))
            out.append(row)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="the block this query set was drawn from; carried into the "
                         "manifest so a row can be traced to a block")
    ap.add_argument("--program", required=True,
                    choices=["blastn", "blastx", "diamond-blastx"])
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--db-version", default="",
                    help="the pinned NCBI prefix or UniProt release this db came from. "
                         "A database with no version is not citable and the whole point "
                         "of staging our own copy was to have one")
    ap.add_argument("--query", type=Path, required=True)
    ap.add_argument("--label", required=True, help="names the output files")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-targets", type=int, default=5)
    ap.add_argument("--sensitive", action="store_true")
    ap.add_argument("--blastn-task", default="blastn",
                    choices=["blastn", "dc-megablast", "megablast"],
                    help="which nucleotide task. Recorded in the timing file, because the "
                         "task is a bigger lever on cost than the database is")
    ap.add_argument("--blast-bin", type=Path, default=Path("/opt/hapax/blast/bin"))
    ap.add_argument("--diamond-bin", type=Path, default=Path("/opt/hapax/diamond"))
    ap.add_argument("--out", type=Path, default=Path("/opt/hapax/out"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    hits = args.out / f"hits_{args.label}.tsv"
    n_q = count_queries(args.query)
    cmd = command(args.program, args.db, args.query, hits, args.threads,
                  args.max_targets, args.sensitive, args.blast_bin, args.diamond_bin,
                  args.blastn_task)

    # The timing file is written BEFORE the search and rewritten after, so a run killed by
    # the dead man's switch still leaves behind what it was doing and with what arguments.
    # The reference project lost two fourteen-hour machines that left nothing at all.
    timing = {"label": args.label, "program": args.program, "db": str(args.db),
              "db_version": args.db_version, "queries": n_q, "threads": args.threads,
              "blastn_task": args.blastn_task if args.program == "blastn" else None,
              "sensitive": bool(args.sensitive), "command": cmd,
              "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "finished": None, "seconds": None, "exit_code": None,
              "platform": platform.platform(), "python": platform.python_version(),
              "cpu_count": os.cpu_count()}
    tpath = args.out / f"timing_{args.label}.json"
    tpath.write_text(json.dumps(timing, indent=2))

    t0 = time.perf_counter()
    proc = subprocess.run(cmd)
    secs = time.perf_counter() - t0
    timing.update({"finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   "seconds": secs, "exit_code": proc.returncode,
                   "seconds_per_1k_queries": secs / n_q * 1000 if n_q else None})
    tpath.write_text(json.dumps(timing, indent=2))
    if proc.returncode != 0:
        raise SystemExit(f"{args.program} exited {proc.returncode}")

    rows = best_rows(hits, args.query)
    rpath = args.out / f"rows_{args.label}.jsonl"
    with open(rpath, "w") as h:
        for r in rows:
            h.write(json.dumps(r) + "\n")

    hit = sum(1 for r in rows if r["hit"])
    print(f"{args.label}: {n_q} queries, {hit} with any hit, {secs:.1f}s "
          f"({secs / n_q * 1000:.2f}s/1k)" if n_q else f"{args.label}: no queries",
          flush=True)


if __name__ == "__main__":
    main()
