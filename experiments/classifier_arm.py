"""Classifier recall against contig length, with the classifier held fixed and
the assembler varied.

The sister thread's literature check found no published experiment that holds a
modern virus classifier fixed and varies the assembler on the same reads: Sutton
et al. 2019 varied 16 assemblers but scored by reference alignment, with no
classifier, and predates geNomad, VirSorter2 and ViraLM. So this arm is the
closure of that gap, and the shape of its output matters: recall is reported as
a curve against contig length, not as a count over 1500 bp, because the 1500 bp
figure in Wu et al. 2024 is that benchmark's harmonisation filter (the minimum
input length all nine of its tools accept) and not a property of any tool.
geNomad states no minimum at all; VIBRANT's four-ORF rule is the real floor.

The spike-derived contigs of every rung and arm are classified in one pass,
labelled by the rung they came from, so one geNomad run serves the whole sweep.
Background contigs from the same assemblies are classified alongside: they are
what a false positive would have to be drawn from, and in a wastewater sample
some of them are genuinely viral, which is why the output reports the rate
rather than calling it an error rate.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path

TAG = re.compile(r"^(s\d+)_n(\d+)_(.+?)\.(spike|bg)\.fa$")
BINS = [(0, 300), (300, 500), (500, 1000), (1000, 1500), (1500, 3000), (3000, 10 ** 9)]


def pool(contig_dir: Path, out: Path) -> tuple[Path, dict[str, dict]]:
    """One fasta of every kept contig, renamed so the rung travels with it.

    Renaming rather than keeping a side table: geNomad's output is keyed by
    contig name, and a name that carries its own provenance cannot be joined to
    the wrong rung later.
    """
    meta: dict[str, dict] = {}
    with out.open("w") as fh:
        for fa in sorted(contig_dir.glob("*.fa")):
            m = TAG.match(fa.name)
            if not m:
                continue
            seed, n, arm, kind = m.group(1), int(m.group(2)), m.group(3), m.group(4)
            i = 0
            for line in fa.open():
                if line.startswith(">"):
                    i += 1
                    name = f"{seed}|{n}|{arm}|{kind}|{i}"
                    meta[name] = {"seed": seed, "spike_pairs": n, "arm": arm,
                                  "kind": kind, "len": 0}
                    fh.write(f">{name}\n")
                    last = name
                else:
                    fh.write(line)
                    meta[last]["len"] += len(line.strip())
    return out, meta


def run_genomad(fasta: Path, db: Path, out: Path, threads: int, min_score: float) -> Path:
    res = out / f"{fasta.stem}_summary" / f"{fasta.stem}_virus_summary.tsv"
    if not res.exists():
        subprocess.run(["genomad", "end-to-end", "--cleanup", "--threads", str(threads),
                        "--min-score", str(min_score), str(fasta), str(out), str(db)],
                       check=True)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", type=Path, default=Path("results/abundance_sweep/contigs"))
    ap.add_argument("--db", type=Path, required=True, help="geNomad database directory")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-score", type=float, default=0.7, help="geNomad default")
    ap.add_argument("--work", type=Path, default=Path("/home/user/sweepdata/clf"))
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep"))
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    fasta, meta = pool(args.contigs, args.work / "pooled.fna")
    print(f"pooled {len(meta)} contigs", flush=True)

    summary = run_genomad(fasta, args.db, args.work / "genomad", args.threads, args.min_score)
    called: dict[str, float] = {}
    with summary.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            called[row["seq_name"]] = float(row.get("virus_score", 0) or 0)

    # Recall by length bin, spike contigs only; and the background rate, which
    # shares the bins so the two can be read against each other.
    by_bin: dict[tuple, dict] = defaultdict(lambda: {"n": 0, "called": 0})
    by_rung: dict[tuple, dict] = defaultdict(lambda: {"n": 0, "called": 0, "mass": 0, "called_mass": 0})
    for name, m in meta.items():
        lo_hi = next((b for b in BINS if b[0] <= m["len"] < b[1]), BINS[-1])
        hit = name in called
        key = (m["kind"], lo_hi)
        by_bin[key]["n"] += 1
        by_bin[key]["called"] += hit
        rk = (m["kind"], m["arm"], m["spike_pairs"])
        by_rung[rk]["n"] += 1
        by_rung[rk]["called"] += hit
        by_rung[rk]["mass"] += m["len"]
        by_rung[rk]["called_mass"] += m["len"] if hit else 0

    result = {
        "contigs": len(meta),
        "genomad_min_score": args.min_score,
        "wall": round(time.perf_counter() - t0, 1),
        "recall_by_length": [
            {"kind": k, "len_from": b[0], "len_to": b[1], "contigs": v["n"],
             "called_virus": v["called"],
             "rate": round(v["called"] / v["n"], 3) if v["n"] else None}
            for (k, b), v in sorted(by_bin.items())
        ],
        "by_rung": [
            {"kind": k, "arm": a, "spike_pairs": n, **v,
             "rate": round(v["called"] / v["n"], 3) if v["n"] else None}
            for (k, a, n), v in sorted(by_rung.items())
        ],
    }
    (args.out / "classifier.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result["recall_by_length"], indent=2))


if __name__ == "__main__":
    main()
