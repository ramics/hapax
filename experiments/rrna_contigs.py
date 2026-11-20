"""Which contigs of the sample's own assembly are ribosomal RNA, and what that
does to the proteins predicted from them.

Why this exists. The index-coverage arm (`results/mate_width/REPORT.md` §4.3)
added proteins predicted by Prodigal from the sample's own assembly to the
peptide index, and the fully-absent background rate fell from 86% to 22%. The
rRNA arm in the sibling thread then found something that cannot happen if that
index holds no ribosomal peptides: reads identified as rRNA are *less* than half
as likely to be fully absent as other reads (0.125 against 0.286). Ribosomal RNA
has no true protein product, so a true index cannot match it. The explanation
that fits is a gene caller artefact - Prodigal is prokaryotic and runs here over
a mixed metatranscriptome, and ribosomal contigs are the highest-coverage,
longest contigs in the assembly, so spurious open reading frames called on them
enter the index as if they were proteins. If that is what happened, part of the
86-to-22 drop measures the gene caller rather than database incompleteness.

This labels the assembly rather than the reads, because a contig offers
thousands of k-mers where a 107 nt read offers 77, so the call is far less
threshold-sensitive. Exact canonical 31-mer containment against SILVA 138.1
SSURef_NR99 + LSURef_NR99, with one k-mer in eight retained by a hash residue:
membership stays exact for retained k-mers and only the number of chances falls,
which a contig can afford. SILVA sequences are stored in the sense orientation
of the gene, so the database is hashed on one strand and each contig on both.

Outputs the three numbers that decide the question: the share of predicted
proteins and of indexed residues sitting on ribosomal contigs, Prodigal's coding
density on those contigs against the assembly as a whole, and a split of the
protein FASTA so `absence_curve.py` can re-run the sample level with them
removed. The split is written rather than the verdict asserted: whether the
removal recovers the drop is the measurement, not this script's claim.
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from peptide_absence import COMP, hash_kmers, read_fasta  # noqa: E402

NT_CODE = np.zeros(256, dtype=np.uint64)
for _i, _b in enumerate("ACGT"):
    NT_CODE[ord(_b)] = _i + 1
    NT_CODE[ord(_b.lower())] = _i + 1
NT_CODE[ord("U")] = NT_CODE[ord("u")] = 4      # SILVA stores RNA


def nt_hashes(seq: str, k: int, keep_mod: int) -> np.ndarray:
    codes = NT_CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    h = hash_kmers(codes, k)
    return h[(h & np.uint64(keep_mod - 1)) == 0] if keep_mod > 1 else h


def silva_set(paths: list[Path], k: int, keep_mod: int, cache: Path) -> np.ndarray:
    if cache.exists():
        return np.load(cache)
    parts, buf, n = [], [], 0
    for p in paths:
        with gzip.open(p, "rt") as fh:
            seq: list[str] = []
            for line in fh:
                if line[0] == ">":
                    if seq:
                        buf.append("".join(seq)); n += 1; seq = []
                    if n % 50_000 == 0 and buf:
                        parts.append(np.unique(nt_hashes("N".join(buf), k, keep_mod)))
                        buf = []
                        print(f"  {n} sequences", flush=True)
                else:
                    seq.append(line.strip())
            if seq:
                buf.append("".join(seq)); n += 1
    if buf:
        parts.append(np.unique(nt_hashes("N".join(buf), k, keep_mod)))
    out = np.unique(np.concatenate(parts))
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache, out)
    print(f"SILVA: {n} sequences, {out.size} retained {k}-mers", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1, help="unused; the labelling is "
                    "deterministic, kept so every experiment takes --seed")
    ap.add_argument("--contigs", type=Path, default=Path("data/ww_berlin/spades_sample.fa"))
    ap.add_argument("--faa", type=Path, default=Path("data/prot/berlin_sample_proteins.faa"))
    ap.add_argument("--silva", type=Path, nargs="+", default=[
        Path("data/silva/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"),
        Path("data/silva/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz")])
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--keep-mod", type=int, default=8, help="retain 1 k-mer in N")
    ap.add_argument("--hit-frac", type=float, default=0.5,
                    help="contig called rRNA at this retained-k-mer hit fraction")
    ap.add_argument("--cds", type=Path,
                    default=Path("data/pepidx/cds_intervals_spades_sample.json"))
    ap.add_argument("--cache", type=Path, default=Path("data/silva/silva_31mers.npy"))
    ap.add_argument("--out", type=Path, default=Path("results/rrna_contigs"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    db = silva_set(args.silva, args.k, args.keep_mod, args.cache)

    cds = json.loads(args.cds.read_text()) if args.cds.exists() else {}
    rows = []
    for name, seq in read_fasta(args.contigs):
        ctg = name.split()[0]
        q = np.concatenate([nt_hashes(seq, args.k, args.keep_mod),
                            nt_hashes(seq.translate(COMP)[::-1], args.k, args.keep_mod)])
        if q.size == 0:
            continue
        pos = np.searchsorted(db, q)
        np.clip(pos, 0, db.size - 1, out=pos)
        frac = float((db[pos] == q).mean())
        iv = cds.get(ctg, [])
        rows.append({"contig": ctg, "len": len(seq), "kmers": int(q.size),
                     "hit_frac": round(frac, 4),
                     "coding_bases": int(sum(e - b + 1 for b, e in iv)),
                     "n_cds": len(iv)})

    rrna = {r["contig"] for r in rows if r["hit_frac"] >= args.hit_frac}

    # split the predicted proteins by whether their contig is ribosomal
    keep, drop, res_keep, res_drop = [], [], 0, 0
    cur: list[str] | None = None
    hdr = ""
    def flush():
        nonlocal res_keep, res_drop
        if cur is None:
            return
        body = "".join(cur)
        ctg = hdr[1:].rsplit("_", 1)[0]
        if ctg in rrna:
            drop.append(hdr + "\n" + body + "\n"); res_drop += len(body)
        else:
            keep.append(hdr + "\n" + body + "\n"); res_keep += len(body)
    for line in args.faa.read_text().splitlines():
        if line.startswith(">"):
            flush(); hdr = line; cur = []
        elif cur is not None:
            cur.append(line.strip())
    flush()

    out_keep = args.faa.with_name(args.faa.stem + "_nonrrna.faa")
    out_drop = args.faa.with_name(args.faa.stem + "_rrnacontigs.faa")
    out_keep.write_text("".join(keep))
    out_drop.write_text("".join(drop))

    tot_len = sum(r["len"] for r in rows)
    rr = [r for r in rows if r["contig"] in rrna]
    nr = [r for r in rows if r["contig"] not in rrna]

    def density(rs):
        L = sum(r["len"] for r in rs)
        return round(sum(r["coding_bases"] for r in rs) / L, 4) if L else None

    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    res = {
        "git_sha": sha, "k": args.k, "keep_mod": args.keep_mod,
        "hit_frac_threshold": args.hit_frac,
        "silva_kmers_retained": int(db.size),
        "contigs": len(rows), "contigs_rrna": len(rr),
        "assembly_bases": tot_len,
        "rrna_bases_frac": round(sum(r["len"] for r in rr) / tot_len, 4) if tot_len else None,
        "coding_density_rrna_contigs": density(rr),
        "coding_density_other_contigs": density(nr),
        "coding_density_all": density(rows),
        "proteins_total": len(keep) + len(drop),
        "proteins_on_rrna_contigs": len(drop),
        "proteins_on_rrna_contigs_frac": round(len(drop) / max(len(keep) + len(drop), 1), 4),
        "residues_total": res_keep + res_drop,
        "residues_on_rrna_contigs": res_drop,
        "residues_on_rrna_contigs_frac": round(res_drop / max(res_keep + res_drop, 1), 4),
        "hit_frac_histogram": {
            f"{lo:.1f}-{lo + 0.1:.1f}": int(sum(1 for r in rows if lo <= r["hit_frac"] < lo + 0.1))
            for lo in [i / 10 for i in range(10)]},
        "faa_nonrrna": str(out_keep), "faa_rrna_contigs": str(out_drop),
        "seconds": round(time.perf_counter() - t0, 1),
    }
    (args.out / "rrna_contigs.json").write_text(json.dumps(res, indent=2))
    with (args.out / "contig_rrna.tsv").open("w") as fh:
        fh.write("contig\tlen\tkmers\thit_frac\tcoding_bases\tn_cds\trrna\n")
        for r in rows:
            fh.write(f"{r['contig']}\t{r['len']}\t{r['kmers']}\t{r['hit_frac']}\t"
                     f"{r['coding_bases']}\t{r['n_cds']}\t{int(r['contig'] in rrna)}\n")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
