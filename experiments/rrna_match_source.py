"""Which public protein entries do ribosomal reads match, and what are they?

Why this exists. Joining the sibling thread's SILVA calls onto the per-level
absence dumps (`results/absence_curve/`) produced the most surprising number of
this thread: P(absent | rRNA) is 0.998 against a viral-only index, exactly as it
should be for sequence with no protein product, and then *crosses below* the
non-ribosomal rate when UniRef50 enters the index - 0.473 against 0.720 at a
tenth of UniRef50, 0.143 against 0.252 with the sample's own proteins added.
Ribosomal RNA has no true peptide, so something in the public database matches
its translations, and the two explanations on the table both failed:

- spurious gene calls on the sample's own ribosomal contigs are 0.24% of the
  index and, measured directly, move P(absent | rRNA) only from 0.107 to 0.143
  when removed (`experiments/rrna_contigs.py`, §4.5);
- low complexity does nothing at all - restricting to k-mers with at least 5 of
  10 distinct residues leaves the figure at 0.144 against 0.143.

That leaves the database itself. UniRef50 includes TrEMBL, which is
automatically annotated, and a gene caller run over a genome assembly makes the
same mistake on ribosomal loci that Prodigal makes on our contigs - so the
artefact the sibling thread proposed may be real but living in the public
database at scale rather than in our 8 Maa of sample proteins.

This names the entries instead of inferring them. Every six-frame peptide 10-mer
of the ribosomal reads becomes a query, the UniRef50 subsample is streamed once,
and any entry containing a query k-mer is recorded with its description. The
non-ribosomal reads are carried through the same pass as a control, because the
description profile of "what wastewater reads match" is only interpretable
against "what ribosomal wastewater reads match". A database of real proteins
should look like real proteins; a tail of automatic annotation should look like
uncharacterised ones.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from insert_profile import read_pairs                     # noqa: E402
from peptide_absence import AA_CODE, hash_kmers, six_frames  # noqa: E402

UNINFORMATIVE = re.compile(
    r"uncharacteri[sz]ed|hypothetical|predicted protein|putative protein|"
    r"unnamed protein|ORF|open reading frame|DUF\d+", re.I)


def read_kmers(read: str, k: int) -> np.ndarray:
    hs = []
    for f in six_frames(read):
        for seg in f.split("*"):
            if len(seg) < k:
                continue
            codes = AA_CODE[np.frombuffer(seg.encode(), dtype=np.uint8)]
            h = hash_kmers(codes, k)
            if h.size:
                hs.append(h)
    return np.concatenate(hs) if hs else np.zeros(0, dtype=np.uint64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=20000)
    ap.add_argument("--block", type=int, default=400000)
    ap.add_argument("--r1", type=Path, default=Path("data/ww_berlin/SRR23973193_1.fastq.gz"))
    ap.add_argument("--r2", type=Path, default=Path("data/ww_berlin/SRR23973193_2.fastq.gz"))
    ap.add_argument("--trim-r2", type=int, default=14)
    ap.add_argument("--rrna-tsv", type=Path, required=True,
                    help="the sibling thread's rrna_s{seed}.tsv: read_id, mate, rrna")
    ap.add_argument("--db", type=Path,
                    default=Path("data/prot/uniref50_every10.fasta.gz"))
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--out", type=Path, default=Path("results/rrna_match_source"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    rrna = {}
    with args.rrna_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rrna[(r["read_id"], int(r["mate"]))] = int(r["rrna"])

    skip = (args.seed - 1) * args.block
    qs: list[np.ndarray] = []
    grp: list[np.ndarray] = []           # 1 = ribosomal read, 0 = not
    n = {0: 0, 1: 0}
    for name, a, b in read_pairs(args.r1, args.r2, args.pairs, args.trim_r2, skip):
        for mate, read in ((1, a), (2, b)):
            g = rrna.get((name, mate))
            if g is None:
                continue
            h = read_kmers(read, args.k)
            if h.size == 0:
                continue
            qs.append(h)
            grp.append(np.full(h.size, g, dtype=np.uint8))
            n[g] += 1
    q = np.concatenate(qs)
    g = np.concatenate(grp)
    order = np.argsort(q, kind="stable")
    q_sorted, g_sorted = q[order], g[order]
    uq, first = np.unique(q_sorted, return_index=True)
    # a unique k-mer is attributed to ribosomal reads if any ribosomal read has it
    bounds = np.append(first, uq.size and q_sorted.size)
    is_r = np.zeros(uq.size, dtype=bool)
    for i in range(uq.size):
        is_r[i] = bool(g_sorted[bounds[i]:bounds[i + 1]].any())
    print(f"queries {q.size} unique {uq.size} ribosomal-reads {n[1]} other {n[0]}",
          flush=True)

    desc_r: Counter = Counter()
    desc_o: Counter = Counter()
    hits_r = hits_o = entries = 0
    name = ""
    seq: list[str] = []

    def handle():
        nonlocal hits_r, hits_o
        if not seq or not name:
            return
        s = "".join(seq)
        codes = AA_CODE[np.frombuffer(s.encode(), dtype=np.uint8)]
        h = hash_kmers(codes, args.k)
        if h.size == 0:
            return
        pos = np.searchsorted(uq, h)
        np.clip(pos, 0, uq.size - 1, out=pos)
        m = uq[pos] == h
        if not m.any():
            return
        hit_r = bool(is_r[pos[m]].any())
        hit_o = bool((~is_r[pos[m]]).any())
        # description: UniRef headers are "ID description n=.. Tax=.."
        d = name.split(" ", 1)[1] if " " in name else ""
        d = d.split(" n=")[0].strip()
        if hit_r:
            desc_r[d] += 1; hits_r += 1
        if hit_o:
            desc_o[d] += 1; hits_o += 1

    with gzip.open(args.db, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                handle()
                entries += 1
                if entries % 500_000 == 0:
                    print(f"  {entries} entries, {hits_r} ribosomal-matching",
                          flush=True)
                name = line[1:].rstrip()
                seq = []
            else:
                seq.append(line.strip())
    handle()

    def profile(c: Counter) -> dict:
        tot = sum(c.values())
        unin = sum(v for d, v in c.items() if UNINFORMATIVE.search(d))
        return {"entries_matched": tot,
                "uninformative_frac": round(unin / tot, 4) if tot else None,
                "top": c.most_common(args.top)}

    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    res = {"seed": args.seed, "git_sha": sha, "db": str(args.db), "k": args.k,
           "db_entries": entries,
           "reads_ribosomal": n[1], "reads_other": n[0],
           "query_kmers_unique": int(uq.size),
           "ribosomal": profile(desc_r), "other": profile(desc_o),
           "seconds": round(time.perf_counter() - t0, 1)}
    (args.out / f"rrna_match_source_s{args.seed}.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items()
                      if k not in ("ribosomal", "other")}, indent=2))
    for lbl in ("ribosomal", "other"):
        print(f"\n{lbl}: {res[lbl]['entries_matched']} entries, "
              f"uninformative {res[lbl]['uninformative_frac']}")
        for d, c in res[lbl]["top"][:15]:
            print(f"  {c:7d}  {d[:80]}")


if __name__ == "__main__":
    main()
