"""Label every read of the whole library against SILVA and the sample's own assembly.

WHY THIS EXISTS. Every residue experiment so far looked at a residue read on its own: its
composition, its frames, its best hit. A 107 nt read carries very little, and the
continuum result (`results/residue_structure/`) says that per-read fields do not separate
the residue from anything. What a read does carry, and nothing has used for composition,
is its PLACE in the library: other reads cover the same sequence, and each of those has a
mate up to ~500 nt away. If those mates are ribosomal, the residue read sits on an rRNA
precursor (an external or internal transcribed spacer, which SILVA does not hold). If they
land on an identified contig, the residue read is the untranslated or intergenic flank of
something known. If nothing anywhere near it is known, it is a whole unknown molecule.

That question needs a label for every read in the library, not just the 100,000-pair
seed blocks, because a residue read's neighbours can be anywhere in the run. This script
makes that table once; `residue_linkage.py` reads it.

WHAT IS LABELLED, per mate, as continuous matching-base counts (never a boolean - the
operating point is applied by the reader, and the cascade's is >=40 matching bases):
  - `silva_mlen`: best minimap2 (sr) alignment to SILVA 138.2 SSU+LSU NR99.
  - `contig`, `contig_mlen`: best alignment to this sample's own SPAdes contigs
    (GSE228220, `WW_211013_20-22hmix_S5_`), all 82,947 of them, identified or not. Which
    contigs count as identified is a join made later against the cascade's own split, so
    this table does not depend on a classification it cannot reproduce here.

CALIBRATION. The table is checked against the cascade it must agree with: on the five seed
blocks, the fraction of pairs with either mate at silva_mlen >= 40 must reproduce the
rrna tier's 38.02% (37.95-38.13) in `results/subtract_recut/TABLES.md`, and every residue
pair must read below 40 on SILVA in both mates. A mismatch means this is a different
instrument and nothing downstream may be compared with the cascade.

Mate 2 is trimmed by 14 nt, the measured non-template prefix (`results/nontemplate_prefix/`),
the same operating point as the cascade's residue FASTAs.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import _silva_to_dna, provenance, _peak_rss  # noqa: E402

SAMPLE_PREFIX = "WW_211013_20-22hmix_S5_"


def extract_sample_contigs(geo: Path, out: Path) -> int:
    """This sample's contigs from the all-sample GEO file.

    The GEO file is truncated after sample S65 (data/MANIFEST.md), but this sample's
    records come before the break, so a read error after we have seen them and then left
    them is not a loss; one before is, and is raised.
    """
    n, inside = 0, False
    try:
        with gzip.open(geo, "rt") as g, open(out, "w") as h:
            for line in g:
                if line.startswith(">"):
                    was = inside
                    inside = line[1:].startswith(SAMPLE_PREFIX)
                    if was and not inside:
                        break  # the sample's records are contiguous; we have left them
                    n += inside
                if inside:
                    h.write(line)
    except EOFError:
        if n == 0 or inside:
            raise
    return n


def pairs(f1: Path, f2: Path, trim_r2: int, chunk: int):
    """Chunks of (pair_index, s1, s2). Pair index is the SRR spot number, which is how
    the cascade's read names (SRR23973193.N) join back to this table."""
    buf = []
    with gzip.open(f1, "rt") as h1, gzip.open(f2, "rt") as h2:
        while True:
            a = [h1.readline() for _ in range(4)]
            b = [h2.readline() for _ in range(4)]
            if not a[0]:
                break
            spot = int(a[0][1:].split()[0].split(".")[1])
            buf.append((spot, a[1].strip(), b[1].strip()[trim_r2:]))
            if len(buf) >= chunk:
                yield buf
                buf = []
    if buf:
        yield buf


def best(aligner, seq):
    """Best hit as (matching bases, target, target start, target length, strand).

    The position is kept because the composition questions downstream are positional: a
    read beside a ribosomal gene's 3' end is on a spacer, and a read beside a contig's
    end is on the flank the assembly did not reach.
    """
    top, ctg, st, tl, strand = 0, None, -1, 0, 0
    for h in aligner.map(seq):
        if h.mlen > top:
            top, ctg, st, tl, strand = h.mlen, h.ctg, h.r_st, h.ctg_len, h.strand
    return top, ctg, st, tl, strand


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--geo", type=Path, default=ROOT / "data/geo/GSE228220_Spades_RNA_all.fa.gz")
    ap.add_argument("--silva", type=Path, nargs="+",
                    default=[ROOT / "data/silva/SSU.fa.gz", ROOT / "data/silva/LSU.fa.gz"])
    ap.add_argument("--trim-r2", type=int, default=14)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--limit", type=int, default=0, help="pairs to label; 0 = all")
    ap.add_argument("--work", type=Path, default=ROOT / "data/work")
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded for the convention; the labelling is deterministic")
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    import mappy

    t0 = time.perf_counter()
    timings = {}
    contigs_fa = args.work / "sample_contigs_all.fa"
    if not contigs_fa.exists():
        n = extract_sample_contigs(args.geo, contigs_fa)
        print(f"extracted {n} contigs", flush=True)
    silva_fa = args.work / "silva.fa"
    if not silva_fa.exists():
        _silva_to_dna(args.silva, silva_fa)
    timings["prepare"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    silva = mappy.Aligner(str(silva_fa), preset="sr", fn_idx_out=str(args.work / "silva.mmi")) \
        if not (args.work / "silva.mmi").exists() else mappy.Aligner(str(args.work / "silva.mmi"), preset="sr")
    asm = mappy.Aligner(str(contigs_fa), preset="sr")
    names = list(asm.seq_names)
    cidx = {c: i for i, c in enumerate(names)}
    timings["index"] = time.perf_counter() - t0
    print(f"indexes ready {timings['index']:.0f}s, {len(names)} contigs", flush=True)

    def work(chunk):
        out = np.zeros((len(chunk), 15), dtype=np.int64)
        for j, (spot, s1, s2) in enumerate(chunk):
            row = [spot]
            for s in (s1, s2):
                a, _sc, ast, atl, astr = best(silva, s)
                c, cc, cst, _ctl, cstr = best(asm, s)
                row += [a, ast, atl * astr, cidx[cc] if cc else -1, c, cst, cstr]
            out[j] = row
        return out

    t0 = time.perf_counter()
    parts, done = [], 0
    with ThreadPoolExecutor(args.threads) as ex:
        it = pairs(args.r1, args.r2, args.trim_r2, args.chunk)
        pending = []
        for chunk in it:
            pending.append(ex.submit(work, chunk))
            if len(pending) >= args.threads * 2:
                r = pending.pop(0).result()
                parts.append(r)
                done += len(r)
                if done % 500000 < args.chunk:
                    print(f"{done} pairs {time.perf_counter() - t0:.0f}s", flush=True)
            if args.limit and done + len(pending) * args.chunk >= args.limit:
                break
        for p in pending:
            parts.append(p.result())
    tab = np.concatenate(parts)
    tab = tab[np.argsort(tab[:, 0])]
    timings["label"] = time.perf_counter() - t0

    # Per mate m in (1, 2): silva{m} matching bases, silva_st{m} start on the SILVA
    # record, silva_len{m} record length signed by strand; ctg{m} contig index into
    # contig_names, cm{m} matching bases, ctg_st{m} start on the contig, ctg_str{m} strand.
    cols = {}
    for m, base in ((1, 1), (2, 8)):
        cols[f"silva{m}"] = tab[:, base].astype(np.int16)
        cols[f"silva_st{m}"] = tab[:, base + 1].astype(np.int32)
        cols[f"silva_len{m}"] = tab[:, base + 2].astype(np.int32)
        cols[f"ctg{m}"] = tab[:, base + 3].astype(np.int32)
        cols[f"cm{m}"] = tab[:, base + 4].astype(np.int16)
        cols[f"ctg_st{m}"] = tab[:, base + 5].astype(np.int32)
        cols[f"ctg_str{m}"] = tab[:, base + 6].astype(np.int8)
    np.savez_compressed(args.work / "library_labels.npz", spot=tab[:, 0].astype(np.int32),
                        contig_names=np.array(names), **cols)
    man = {"script": "experiments/library_labels.py", **provenance(__file__),
           "parameters": {k: str(v) for k, v in vars(args).items()},
           "pairs_labelled": int(len(tab)), "timings": timings, "peak_rss_bytes": _peak_rss(),
           "mappy": mappy.__version__ if hasattr(mappy, "__version__") else "unknown"}
    (args.work / "library_labels_manifest.json").write_text(json.dumps(man, indent=2))
    print(json.dumps(man["timings"]), len(tab))


if __name__ == "__main__":
    main()
