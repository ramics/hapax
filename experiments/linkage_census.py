"""The linkage census: what does the residue sit NEXT TO?

THE QUESTION, in Imo's framing. Coverage says how often a read's sequence recurs. It cannot
say what that sequence is physically attached to. But a paired library answers adjacency for
free: find the other library pairs that cover the same sequence as a residue read, then look
at THEIR other mates. Each of those is a piece of sequence one insert away from something
the residue read also covers. Label them, and the residue splits into two populations that
mean quite different things:

  - **next to sequence we can name** - rRNA precursor spacers, UTRs, intergenic sequence,
    the untranslated flanks of contigs we identified. This is residue by annotation, not by
    biology. It is unidentified because the reference is a gene catalogue, not because the
    molecule is unknown.
  - **a neighbourhood as unknown as the read** - whole unknown molecules, which is the
    population everything else in this project is actually about.

Telling those apart changes what the residue IS, which is why it is worth a pass of its own.

THE CONTROL IS NOT OPTIONAL AND IS RUN IN THE SAME PASS. Already-identified pairs go through
the identical census. A residue read's neighbourhood looking unknown means nothing until it
is set beside what an identified read's neighbourhood looks like: in this project every
separator that turned out to be describing the database read the same on both piles, and
the tell was always the control moving with the test. There is a second control too, a
label permutation, which says what fraction of queries would be called "next to known" from
the label mix alone.

THREE BIASES, EACH NAMED AND NONE OF THEM REMOVABLE HERE.

 1. **A witness is found through shared sequence, so this is biased towards abundant
    residue.** A singleton molecule has no witnesses at all. Those queries are reported as
    `no_witnesses`, which is a different statement from "unknown neighbourhood" and is
    never folded into it.
 2. **About a third of this library's pairs overlap their own mate**, so for those the
    "other" mate is partly the same sequence as the matching one. Every row carries both
    mates' shared k-mer counts, and the census reports with and without the overlapping
    witnesses.
 3. **Labels come from the same indexes the cascade used**, so a neighbour labelled `none`
    is unidentified at the cascade's own operating point and not in any absolute sense. The
    operating point travels with the number.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import _mm_index, map_filter, provenance
from loud import Preflight

ROOT = Path(__file__).resolve().parents[1]

#: Neighbour classes, ordered from most to least identified. The order is the reading
#: order of the report and nothing more; it is not a score.
CLASSES = ["rrna", "identified_contig", "unidentified_contig", "none"]


def wilson(k: int, n: int, z: float = 1.96) -> dict | None:
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {"point": round(p, 5), "lo": round(max(0.0, c - h), 5),
            "hi": round(min(1.0, c + h), 5), "k": k, "n": n}


def read_fasta_pairs(path: Path, limit: int) -> list[tuple[str, str, str]]:
    """The cascade writes a pair as two records, `<name>/1` and `<name>/2`."""
    seqs: dict[str, dict[str, str]] = {}
    name = None
    buf: list[str] = []
    order: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    base, _, mate = name.rpartition("/")
                    seqs.setdefault(base, {})[mate] = "".join(buf)
                    if base not in order:
                        order.append(base)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        base, _, mate = name.rpartition("/")
        seqs.setdefault(base, {})[mate] = "".join(buf)
        if base not in order:
            order.append(base)
    out = []
    for base in order:
        d = seqs[base]
        if "1" in d and "2" in d:
            out.append((base, d["1"], d["2"]))
        if len(out) >= limit:
            break
    return out


def census(binary: Path, r1: Path, r2: Path, trim: int, queries: list, work: Path,
           tag: str, min_shared: int, max_witnesses: int, seed: int) -> Path:
    q = work / f"census_query_{tag}.tsv"
    with open(q, "w") as h:
        for name, s1, s2 in queries:
            h.write(f"{name}\t{s1}\t{s2}\n")
    rows = work / f"census_rows_{tag}.tsv"
    err = work / f"census_{tag}.err"
    with open(rows, "w") as out, open(err, "w") as e:
        proc = subprocess.run([str(binary), str(r1), str(r2), str(trim), str(q),
                               str(min_shared), str(max_witnesses), str(seed)],
                              stdout=out, stderr=e)
    if proc.returncode != 0:
        # The census writes its stderr to a file, so without this the log says only that
        # it failed.
        raise SystemExit(f"census {tag} exited {proc.returncode}:\n"
                         + "\n".join(err.read_text().splitlines()[-20:]))
    print(f"[{tag}] {err.read_text().strip().splitlines()[-1]}", flush=True)
    return rows


def label(neighbours: dict[str, str], refs: dict[str, object], min_mlen: int) -> dict:
    """Label each neighbour sequence at the cascade's own operating point.

    Order matters and is the cascade's: rRNA first with a NUCLEOTIDE index, because protein
    indexes fail on rRNA in the dangerous direction. A neighbour is labelled by the first
    reference it matches at `min_mlen` bases, so the classes are exclusive and the counts
    add up.
    """
    out = {k: "none" for k in neighbours}
    remaining = dict(neighbours)
    for cls in CLASSES[:-1]:
        idx = refs.get(cls)
        if idx is None or not remaining:
            continue
        hits = map_filter(idx, {k: (v, "") for k, v in remaining.items()}, min_mlen)
        for k in hits:
            out[k] = cls
        remaining = {k: v for k, v in remaining.items() if k not in hits}
    return out


def summarise(rows_path: Path, labels: dict[str, str], queries: list,
              drop_overlapping: bool) -> dict:
    """Per query, the mix of its neighbours' labels, then the piles those mixes fall into."""
    per: dict[str, Counter] = {name: Counter() for name, _s1, _s2 in queries}
    total_seen: dict[str, int] = {name: 0 for name, _s1, _s2 in queries}
    dropped = 0
    with open(rows_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        need = ["query", "witness", "shared_mate1", "shared_mate2", "witnesses_total",
                "other_mate", "other_seq"]
        if head != need:
            raise SystemExit(f"{rows_path}: unexpected header {head}")
        for line in fh:
            f = dict(zip(need, line.rstrip("\n").split("\t")))
            c1, c2 = int(f["shared_mate1"]), int(f["shared_mate2"])
            if drop_overlapping and min(c1, c2) > 0:
                # Both mates of the witness match the query, so its "other" mate is not a
                # neighbourhood - it is more of the same sequence.
                dropped += 1
                continue
            q = f["query"]
            if q not in per:
                continue
            total_seen[q] = int(f["witnesses_total"])
            per[q][labels.get(f["witness"] + "|" + f["other_mate"], "none")] += 1

    piles = Counter()
    rows = []
    for name, _s1, _s2 in queries:
        c = per[name]
        n = sum(c.values())
        if n == 0:
            pile = "no_witnesses"
        else:
            known = (c["rrna"] + c["identified_contig"]) / n
            asm = c["unidentified_contig"] / n
            unk = c["none"] / n
            pile = ("beside_known" if known >= 0.5 else
                    "beside_assembled_unidentified" if asm >= 0.5 else
                    "unknown_neighbourhood" if unk >= 0.5 else "mixed")
        piles[pile] += 1
        rows.append({"query": name, "pile": pile, "neighbours_sampled": n,
                     "witnesses_total": total_seen[name], **{k: c[k] for k in CLASSES}})
    n = len(queries)
    return {"piles": {k: wilson(v, n) for k, v in piles.items()},
            "pile_counts": dict(piles), "queries": n,
            "witnesses_dropped_as_overlapping": dropped,
            "neighbour_label_mix": {k: sum(p[k] for p in per.values()) for k in CLASSES},
            "rows": rows}


def permuted_null(summary: dict, rng: random.Random) -> dict:
    """What the piles would read if labels were dealt out at random.

    Not a bias statement: it says what the pile fractions would be from the label MIX
    alone, given each query's number of neighbours. A census whose piles match this is
    reporting the library's overall composition and nothing about adjacency.
    """
    mix = summary["neighbour_label_mix"]
    pool = [k for k, v in mix.items() for _ in range(v)]
    if not pool:
        return {}
    piles = Counter()
    for r in summary["rows"]:
        n = r["neighbours_sampled"]
        if n == 0:
            piles["no_witnesses"] += 1
            continue
        c = Counter(rng.choice(pool) for _ in range(n))
        known = (c["rrna"] + c["identified_contig"]) / n
        asm = c["unidentified_contig"] / n
        unk = c["none"] / n
        piles["beside_known" if known >= 0.5 else
              "beside_assembled_unidentified" if asm >= 0.5 else
              "unknown_neighbourhood" if unk >= 0.5 else "mixed"] += 1
    n = len(summary["rows"])
    return {k: wilson(v, n) for k, v in piles.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cascade-work", type=Path, default=Path("/tmp/claude-0/wk/recut"),
                    help="where subtract_known.py wrote residue_s*.fa and subtracted_s*.fa")
    ap.add_argument("--queries", type=int, default=5000,
                    help="pairs per pile. Both piles are always run: the control is what "
                         "makes the residue number readable")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--trim-r2", type=int, default=14)
    ap.add_argument("--binary", type=Path,
                    default=ROOT / "experiments/kmer_coverage/target/release/census")
    ap.add_argument("--min-shared", type=int, default=10,
                    help="31-mers a witness must share with the query. Ten 31-mers is "
                         "about 40 bases of overlap, the cascade's own match length")
    ap.add_argument("--max-witnesses", type=int, default=30,
                    help="witnesses kept per query, reservoir sampled over the library so "
                         "the cap is not a statement about one region of one lane")
    ap.add_argument("--min-mlen", type=int, default=40,
                    help="the cascade's operating point, so a neighbour labelled none is "
                         "unidentified at the same threshold the residue was defined by")
    # Accepted so old command lines still parse, and refused if given: nothing in this
    # script reads SILVA. It used to default to the SILVA paths, which put an index in
    # every manifest's parameters that the census never consulted.
    ap.add_argument("--silva", type=Path, nargs="*", default=[],
                    help="not read by this script; giving it is an error")
    ap.add_argument("--identified-contigs", type=Path,
                    default=ROOT / "data/geo/sample_contigs_identified.fa")
    ap.add_argument("--unidentified-contigs", type=Path,
                    default=ROOT / "data/geo/sample_contigs_unidentified.fa")
    ap.add_argument("--rrna-contigs", type=Path,
                    default=ROOT / "data/geo/sample_contigs_rrna.fa",
                    help="the sample's own ribosomal contigs: 1.4 MB of them removes more "
                         "of this library than all 272 MB of SILVA, so they are the rRNA "
                         "reference. There is no fallback: absent means refuse")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/census"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/linkage_census")
    ap.add_argument("--preflight", action="store_true",
                    help="check every input and tool, then exit: 0 if the run could "
                         "start, 2 if not")
    args = ap.parse_args()

    # A class with no reference used to print a WARNING and carry on, so every neighbour
    # of that class was labelled `none` - unidentified - which inflates exactly the
    # "neighbourhood as unknown as the read" number this census exists to measure.
    pf = Preflight("linkage_census.py")
    pf.exe(args.binary, "--binary (cargo build --release --bin census)")
    pf.file(args.r1, "--r1")
    pf.file(args.r2, "--r2")
    pf.file(args.rrna_contigs, "--rrna-contigs (class rrna)")
    pf.file(args.identified_contigs, "--identified-contigs (class identified_contig)")
    pf.file(args.unidentified_contigs, "--unidentified-contigs (class unidentified_contig)")
    pf.require(not args.silva, f"--silva {args.silva} given, but this script never reads "
                               f"SILVA; the rRNA class is labelled from --rrna-contigs")
    if not args.preflight:
        for tag in ("residue", "subtracted"):
            pf.file(args.cascade_work / f"{tag}_s{args.seed}.fa",
                    f"cascade work ({tag} pile)")
    pf.module("mappy", "labelling")
    pf.done(args.preflight)

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    piles = {}
    for tag, fa in (("residue", args.cascade_work / f"residue_s{args.seed}.fa"),
                    ("identified", args.cascade_work / f"subtracted_s{args.seed}.fa")):
        if not fa.exists():
            raise SystemExit(f"{fa} is absent; the cascade's work directory is the input "
                             f"and the control pile is not optional")
        piles[tag] = read_fasta_pairs(fa, args.queries)
        print(f"{tag}: {len(piles[tag])} query pairs from {fa.name}", flush=True)

    refs: dict[str, object] = {}
    refs["rrna"] = _mm_index(args.rrna_contigs, args.work / "rrna.mmi", args.threads)
    refs["identified_contig"] = _mm_index(args.identified_contigs,
                                          args.work / "identified.mmi", args.threads)
    refs["unidentified_contig"] = _mm_index(args.unidentified_contigs,
                                            args.work / "unidentified.mmi", args.threads)
    missing = [c for c in CLASSES[:-1] if c not in refs]
    if missing:
        raise SystemExit(f"no reference for classes {missing}: their neighbours would all "
                         f"be labelled none")

    rng = random.Random(args.seed)
    result = {}
    for tag, queries in piles.items():
        rows = census(args.binary, args.r1, args.r2, args.trim_r2, queries, args.work,
                      f"{tag}_s{args.seed}", args.min_shared, args.max_witnesses,
                      args.seed)
        # Label once per (witness, other mate), not once per row: a neighbour witnessing
        # several queries is the same sequence each time.
        seqs: dict[str, str] = {}
        with open(rows) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                seqs[f"{f[1]}|{f[5]}"] = f[6]
        print(f"[{tag}] labelling {len(seqs)} distinct neighbour mates", flush=True)
        labels = label(seqs, refs, args.min_mlen)
        for drop in (False, True):
            s = summarise(rows, labels, queries, drop_overlapping=drop)
            key = f"{tag}_without_overlapping_witnesses" if drop else tag
            s.pop("rows") if drop else None
            result[key] = s
        result[tag]["permuted_label_null"] = permuted_null(result[tag], rng)
        # The per-query rows are written out rather than kept in the summary: they are the
        # population the report's filters run over.
        with open(args.out / f"census_rows_{tag}_s{args.seed}.jsonl", "w") as h:
            for r in result[tag].pop("rows"):
                h.write(json.dumps(r) + "\n")

    out = {"script": "experiments/linkage_census.py", **provenance(__file__),
           "parameters": {k: (str(v) if isinstance(v, Path) else
                              [str(x) for x in v] if isinstance(v, list) else v)
                          for k, v in vars(args).items()},
           "classes": CLASSES,
           "references": {k: str(getattr(args, {"rrna": "rrna_contigs",
                                                "identified_contig": "identified_contigs",
                                                "unidentified_contig":
                                                    "unidentified_contigs"}[k]))
                          for k in refs},
           "result": result,
           "readout": "for each query pair, the labels of the OTHER mates of the library "
                      "pairs that cover the same sequence; then the pile each query's "
                      "neighbourhood falls into",
           "bias": "witnesses are found through shared sequence, so the census is blind to "
                   "singleton molecules, which land in no_witnesses and are never folded "
                   "into unknown_neighbourhood. About a third of pairs overlap their own "
                   "mate, so the census is reported with and without witnesses whose two "
                   "mates both match. Labels are at the cascade's own 40-base operating "
                   "point, so none means unidentified there and not in any absolute sense",
           "wall_seconds": round(time.perf_counter() - t0, 1),
           "python": platform.python_version(), "platform": platform.platform()}
    p = args.out / f"linkage_census_s{args.seed}.json"
    p.write_text(json.dumps(out, indent=2))

    for tag in ("residue", "identified"):
        r = result[tag]
        line = "  ".join(f"{k} {v['point'] * 100:.1f}%"
                         for k, v in sorted(r["piles"].items()))
        print(f"{tag:11} {line}", flush=True)
    print(f"written {p}", flush=True)


if __name__ == "__main__":
    main()
