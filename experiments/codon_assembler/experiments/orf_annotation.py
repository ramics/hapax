"""Region annotation against GenBank CDS truth on real viral genomes.

hapax labels every consensus base noncoding, coding or dual-coding. Those
labels feed `regions.tsv` and the syn/nonsyn effect of each variant, so a
caller that over-calls dual-coding mislabels a large share of real genomes.
The wastewater probe found this on astrovirus. This experiment scores each
ORF caller base by base against GenBank CDS annotations:

* human astrovirus NC_001943: ORF1a/1b ribosomal frameshift plus ORF2,
  essentially single-coding (true dual-coding ~8 nt),
* HIV-1 HXB2 K03455: genuine overlaps (vif/pol, vpr/vif, tat/rev inside env),
* HBV NC_003977: the dual-coding extreme (S, C and X nested in or across P).

Truth frames come from the CDS features, with joins (including the
frameshift junctions) taken piece by piece. A base is truly dual when two
different frames cover it. The genomes are fetched with NCBI efetch into
`data/` and listed in data/MANIFEST.md. The script is deterministic, so
`--seed` is accepted only to match the other experiments' interface.

Usage: python experiments/orf_annotation.py --seed 0 --out results/orf_annotation
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
import urllib.request
from pathlib import Path

from common import HAPAX, ROOT, ensure_release_build, write_manifest

GENOMES = {"NC_001943": "human astrovirus 1", "K03455": "HIV-1 HXB2", "NC_003977": "hepatitis B virus (ayw)"}
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype={}&retmode=text"


def fetch(acc: str, data: Path) -> tuple[Path, Path]:
    data.mkdir(parents=True, exist_ok=True)
    fa, gb = data / f"{acc}.fa", data / f"{acc}.gb"
    for path, kind in ((fa, "fasta"), (gb, "gb")):
        if not path.exists():
            path.write_bytes(urllib.request.urlopen(EFETCH.format(acc, kind)).read())
            time.sleep(1.0)  # NCBI asks for <= 3 requests per second without a key
    return fa, gb


def cds_frames(gb: Path, length: int) -> list[set]:
    """Per base, the set of (strand, frame) codings from CDS features."""
    text = gb.read_text()
    locs = []
    for m in re.finditer(r"^     CDS             (.+?)\n(?=                     /)", text, re.M | re.S):
        locs.append(re.sub(r"\s+", "", m.group(1)))
    cover = [set() for _ in range(length)]
    for loc in locs:
        strand = "-" if loc.startswith("complement") else "+"
        parts = [(int(a), int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc)]
        if strand == "-":
            parts = parts[::-1]
        phase = 0  # position within the codon at the start of each part
        for a, b in parts:
            rng = range(a - 1, b) if strand == "+" else range(b - 1, a - 2, -1)
            for k, p in enumerate(rng):
                cod_pos = (phase + k) % 3
                start = p - cod_pos if strand == "+" else p + cod_pos
                cover[p % length].add((strand, start % 3 if strand == "+" else (length - 1 - start) % 3))
            phase = (phase + len(rng)) % 3
    return cover


def predicted(fa: Path, caller: str, min_orf: int, length: int) -> list[str]:
    out = subprocess.check_output([str(HAPAX), "orfs", "-i", str(fa), "--orf-caller", caller, "--min-orf", str(min_orf)], text=True)
    lab = ["noncoding"] * length
    for line in out.splitlines()[1:]:
        _, kind, s, e, v = line.split("\t")
        if kind == "region":
            for p in range(int(s) - 1, int(e)):
                lab[p] = v
    return lab


def score(truth: list[set], pred: list[str]) -> dict:
    t_dual = [len(x) >= 2 for x in truth]
    t_cod = [len(x) >= 1 for x in truth]
    p_dual = [x == "dual" for x in pred]
    p_cod = [x != "noncoding" for x in pred]
    tp = sum(a and b for a, b in zip(t_dual, p_dual))
    n = len(truth)
    return {
        "length": n,
        "true_dual_nt": sum(t_dual),
        "called_dual_nt": sum(p_dual),
        "dual_true_positive_nt": tp,
        "dual_precision": round(tp / sum(p_dual), 3) if sum(p_dual) else None,
        "dual_recall": round(tp / sum(t_dual), 3) if sum(t_dual) else None,
        "false_dual_fraction_of_genome": round((sum(p_dual) - tp) / n, 3),
        "coding_recall": round(sum(a and b for a, b in zip(t_cod, p_cod)) / sum(t_cod), 3),
        "noncoding_called_coding_nt": sum((not a) and b for a, b in zip(t_cod, p_cod)),
        "true_noncoding_nt": n - sum(t_cod),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-orf", type=int, default=300)
    ap.add_argument("--callers", nargs="+", default=["union", "atg", "scored"])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    ensure_release_build()
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    rows = []
    for acc, name in GENOMES.items():
        fa, gb = fetch(acc, ROOT / "data")
        seq = "".join(l.strip() for l in fa.read_text().splitlines()[1:])
        truth = cds_frames(gb, len(seq))
        for caller in a.callers:
            r = {"accession": acc, "virus": name, "caller": caller, **score(truth, predicted(fa, caller, a.min_orf, len(seq)))}
            rows.append(r)
    (a.out / "orf_annotation.json").write_text(json.dumps(rows, indent=2))
    cols = ["accession", "caller", "true_dual_nt", "called_dual_nt", "dual_precision", "dual_recall", "false_dual_fraction_of_genome", "coding_recall", "noncoding_called_coding_nt", "true_noncoding_nt"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(r[c]) for c in cols) + " |" for r in rows]
    (a.out / "orf_annotation.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    write_manifest(a.out, a.seed, {k: str(v) for k, v in vars(a).items()}, {"total": time.perf_counter() - t0}, {"accessions": GENOMES})


if __name__ == "__main__":
    main()
