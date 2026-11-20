"""Split this sample's own assembly into what can be identified and what cannot.

WHY THIS EXISTS - a correction, not an addition. The first cut of the cascade
(`experiments/subtract_known.py`) subtracted the sample's own assembly wholesale, and that
was the most efficient filter measured by a wide margin: 35 MB of contigs removed 83.6% of
pairs alone. But Imo's three-tier architecture (2026-11-20) makes clear that this is the
wrong subtrahend:

  1. remove everything we can assemble AND identify - or identify at all, assembled or not;
  2. check everything we can assemble but CANNOT identify, for novelty;
  3. check the residual reads for genuinely new biology.

Subtracting the assembly wholesale swallows tier 2 entirely: a novel virus abundant enough
to assemble would be deleted as background before anyone looked at it. So the subtrahend
has to be defined by IDENTIFICATION, not by assembly. This script does that split, and the
cascade is re-run against the identified part only, with the unidentified contigs routed
to tier 2 rather than into the bin.

HOW A CONTIG COUNTS AS IDENTIFIED. Two independent routes, because they fail differently
and a contig only needs one of them:

  - nucleotide: Kraken2 assigns it a taxon at the stated confidence. Cheap, and the only
    route that can place non-coding sequence such as rRNA;
  - protein: DIAMOND blastx finds a hit at or above the stated bitscore. Reaches further
    down the divergence axis, which is the whole reason the project kept a protein tier.

Both are swept rather than fixed, and the partition is written at the operating point so
the downstream cascade is reproducible. A contig is "unidentifiable" only in the sense of
"not identified by these indexes at these thresholds", which is the same one-sided
statement the rest of this work makes about absence; the sweep is how that limit is
reported rather than hidden.

The tier-2 population is also described here rather than just counted: length and SPAdes
k-mer coverage from the contig headers, because tier 2's own success measure is
sensitivity on contigs where existing classifiers actually operate, and this project
measured that boundary at ~500 nt and ~2x coverage (PR #5, the abundance sweep).
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import dir_bytes, git_sha, kraken2

ROOT = Path(__file__).resolve().parents[1]

# SPAdes contig names carry length and k-mer coverage: NODE_1_length_12710_cov_727.441011
HEADER = re.compile(r"length_(\d+)_cov_([0-9.]+)")


def read_fasta(path: Path):
    name, buf = None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="unused: the input is a fixed assembly, but the repo requires it")
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--k2-db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--dmnd", type=Path, nargs="+",
                    default=[ROOT / "data/prot/sprot.dmnd",
                             ROOT / "data/prot/uniref50_10pct.dmnd"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--k2-confidence", type=float, nargs="+", default=[0.0, 0.05, 0.1, 0.2])
    ap.add_argument("--operating-confidence", type=float, default=0.1)
    ap.add_argument("--operating-bitscore", type=float, default=50.0)
    ap.add_argument("--min-len", type=int, default=200,
                    help="contigs shorter than this are not routed to tier 2: no "
                         "classifier in the 2024 benchmark operates there, so calling "
                         "them unidentifiable would overstate the tier-2 population")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/contigs"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/subtract_known")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()

    contigs = dict(read_fasta(args.contigs))
    lens = {k: len(v) for k, v in contigs.items()}
    cov = {}
    for k in contigs:
        m = HEADER.search(k)
        if m:
            cov[k] = float(m.group(2))

    # ---- nucleotide identification
    fa = args.work / "contigs.fa"
    with open(fa, "w") as h:
        for k, v in contigs.items():
            h.write(f">{k}\n{v}\n")
    k2_sweep, k2_hit = {}, set()
    timings = {}
    for conf in args.k2_confidence:
        out = args.work / f"k2_{conf}.out"
        t0 = time.perf_counter()
        cmd = [str(args.classify_bin), "-H", str(args.k2_db / "hash.k2d"),
               "-t", str(args.k2_db / "taxo.k2d"), "-o", str(args.k2_db / "opts.k2d"),
               "-M", "-p", str(args.threads), "-T", str(conf), "-O", str(out), str(fa)]
        subprocess.run(cmd, check=True, capture_output=True)
        timings[f"kraken2_conf{conf}"] = time.perf_counter() - t0
        hit = set()
        with open(out) as h:
            for line in h:
                _st, name, taxid, _l, _hits = line.rstrip("\n").split("\t")
                if taxid != "0":
                    hit.add(name)
        k2_sweep[f"{conf:.2f}"] = len(hit) / len(contigs)
        if conf == args.operating_confidence:
            k2_hit = hit
        out.unlink(missing_ok=True)

    # ---- protein identification
    prot_sweep, prot_hit = {}, {}
    for db in args.dmnd:
        if not db.exists():
            continue
        tsv = args.work / f"{db.stem}.tsv"
        t0 = time.perf_counter()
        subprocess.run([str(args.diamond_bin), "blastx", "--db", str(db), "--query", str(fa),
                        "--out", str(tsv), "--threads", str(args.threads), "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "bitscore", "--max-target-seqs", "1",
                        "--evalue", "10", "--sensitive", "--quiet"], check=True)
        timings[f"diamond_{db.stem}"] = time.perf_counter() - t0
        best: dict[str, float] = {}
        with open(tsv) as h:
            for line in h:
                q, _s, _pid, bits = line.rstrip("\n").split("\t")
                bits = float(bits)
                if bits > best.get(q, 0):
                    best[q] = bits
        prot_sweep[db.stem] = {str(b): sum(1 for v in best.values() if v >= b) / len(contigs)
                               for b in (30, 40, 50, 60, 80, 100)}
        prot_hit[db.stem] = {q for q, v in best.items() if v >= args.operating_bitscore}
        tsv.unlink(missing_ok=True)

    identified = set(k2_hit)
    for s in prot_hit.values():
        identified |= s
    unidentified = [k for k in contigs if k not in identified]
    tier2 = [k for k in unidentified if lens[k] >= args.min_len]

    id_fa = args.contigs.with_name(args.contigs.stem + "_identified.fa")
    un_fa = args.contigs.with_name(args.contigs.stem + "_unidentified.fa")
    with open(id_fa, "w") as h:
        for k in identified:
            h.write(f">{k}\n{contigs[k]}\n")
    with open(un_fa, "w") as h:
        for k in tier2:
            h.write(f">{k}\n{contigs[k]}\n")

    def desc(keys):
        ls = [lens[k] for k in keys]
        cs = [cov[k] for k in keys if k in cov]
        return {
            "contigs": len(keys),
            "bases": sum(ls),
            "len_median": st.median(ls) if ls else None,
            "len_p90": sorted(ls)[int(0.9 * len(ls))] if ls else None,
            "len_max": max(ls) if ls else None,
            "cov_median": st.median(cs) if cs else None,
            "frac_over_500nt": sum(1 for x in ls if x >= 500) / len(ls) if ls else None,
            "frac_over_2x": sum(1 for x in cs if x >= 2) / len(cs) if cs else None,
            "frac_over_500nt_and_2x": (sum(1 for k in keys
                                           if lens[k] >= 500 and cov.get(k, 0) >= 2)
                                       / len(keys)) if keys else None,
        }

    result = {
        "contigs_file": str(args.contigs),
        "contigs": len(contigs),
        "bases": sum(lens.values()),
        "operating_point": {"kraken2_confidence": args.operating_confidence,
                            "diamond_bitscore": args.operating_bitscore,
                            "kraken2_db": args.k2_db.name,
                            "protein_dbs": [d.stem for d in args.dmnd if d.exists()]},
        "kraken2_confidence_sweep": k2_sweep,
        "protein_bitscore_sweep": prot_sweep,
        "identified_frac": len(identified) / len(contigs),
        "identified_bases_frac": sum(lens[k] for k in identified) / sum(lens.values()),
        "route_nucleotide_only": len(k2_hit - set().union(*prot_hit.values())) / len(contigs)
        if prot_hit else None,
        "route_protein_only": (len(set().union(*prot_hit.values()) - k2_hit) / len(contigs)
                               if prot_hit else None),
        "all": desc(list(contigs)),
        "identified": desc(list(identified)),
        "unidentified": desc(unidentified),
        "tier2_population": desc(tier2),
        "tier2_min_len": args.min_len,
        "files": {"identified": str(id_fa), "unidentified_tier2": str(un_fa)},
        "index_bytes": {"kraken2": dir_bytes(args.k2_db / "hash.k2d"),
                        **{d.stem: dir_bytes(d) for d in args.dmnd if d.exists()}},
        "timings_seconds": timings,
    }
    (args.out / "contig_identification.json").write_text(json.dumps(result, indent=2))
    (args.out / "contig_identification_manifest.json").write_text(json.dumps({
        "script": "experiments/classify_contigs.py",
        "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }, indent=2))
    fa.unlink(missing_ok=True)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("files", "index_bytes", "timings_seconds")}, indent=2))


if __name__ == "__main__":
    main()
