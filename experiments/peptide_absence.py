"""Six-frame peptide k-mer absence: the index, a per-read score, and how far
apart two positions have to be before their novelty decorrelates.

Why this script exists in a thread about insert size. The wide-insert hypothesis
has two halves. The first - does spreading the mates apart raise the chance of
hitting a coding region - is settled by annotation arithmetic in
`experiments/mate_separation.py`. The second is the one that actually decides the
question: does spreading them apart make the two mates *independent observations
of novelty*? Every pair-level rule in results/readlevel_features/REPORT.md §3.1,
conjunction and disjunction alike, has its merit set by the mate correlation
rho, and insert length is the only experimental knob that moves rho. If novelty
is blocky at gene scale, a few hundred nucleotides of separation buys no
independence at all and the hypothesis is answered; if its correlation length is
short, width is worth paying for.

To measure that, novelty has to be a *continuous track along a genome* rather
than a verdict about a read. So this builds the one feature in the catalogue
with the right polarity (A12, six-frame peptide k-mer absence, the shortlist's
item 2) and uses it as that track:

- index: every amino-acid k-mer of RefSeq viral proteins, with one family
  removed by taxonomy so the target genome is genuinely novel to it. Stored as
  sorted 64-bit hashes; `numpy.searchsorted` answers membership. Hashing is
  vectorised over the concatenated proteome, because a Python loop over ~3e8
  k-mers is hours.
- score: for a nucleotide read, translate all six frames, and report the
  fraction of its peptide k-mers absent from the index. Absence rises with
  novelty, which is the polarity that database proximity (DIAMOND bitscore) gets
  backwards.
- track: slide a window along a genome and score it, giving novelty(position).
  The autocorrelation of that track against lag is the number the insert-size
  question reduces to, and it is read off at the lags the library can actually
  deliver (median 194 nt, p99 485 nt - measured in results/insert_profile/).

k = 10 is the specificity floor from the literature review (at k = 7 over 80% of
peptide k-mers are shared across phyla).

The honest caveats are in the report, but one belongs here: absence is produced
by sequencing error, low-complexity sequence and under-sequenced known organisms
as readily as by novelty, so a high score is not evidence of novelty on its own.
Nothing in this script is a detector; it is the track the insert-size question
needs, and the per-read scorer the rho measurement needs.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import subprocess
import time
from pathlib import Path

import numpy as np

CODON = {}
_B = "TCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODON[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_CODE = np.zeros(256, dtype=np.uint64)
for i, a in enumerate(AA):
    AA_CODE[ord(a)] = i + 1  # 0 is reserved for "not an ordinary residue"

COMP = str.maketrans("ACGTN", "TGCAN")
MASK = np.uint64(0xFFFFFFFFFFFFFFFF)
MULT = np.uint64(1099511628211)


def translate(nt: str) -> str:
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def six_frames(nt: str) -> list[str]:
    rc = nt.translate(COMP)[::-1]
    return [translate(nt[o:]) for o in (0, 1, 2)] + [translate(rc[o:]) for o in (0, 1, 2)]


def hash_kmers(codes: np.ndarray, k: int) -> np.ndarray:
    """FNV-style rolling hash of every length-k window, vectorised.

    `codes` holds residue codes 1..20 with 0 marking a break (protein boundary,
    stop codon, or an unusual residue). Windows containing a 0 are dropped,
    which is what keeps k-mers from spanning two proteins.
    """
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    h = np.zeros(n, dtype=np.uint64)
    ok = np.ones(n, dtype=bool)
    for off in range(k):
        c = codes[off:off + n]
        h = (h * MULT) & MASK
        h ^= c
        ok &= c != 0
    return h[ok]


def proteome_codes(faas: list[Path], exclude_accessions: set[str] | None = None,
                   exclude_names: tuple[str, ...] = ()) -> np.ndarray:
    """Concatenate the kept proteins into one code array, 0 between proteins.

    Several FASTA files are accepted so the index can span more than RefSeq
    viral. That matters for anything measured on real background: wastewater is
    mostly bacterial, so against a viral-only index almost every background read
    is fully absent and the score saturates (measured: 88% of Berlin reads score
    exactly 1.0 - results/mate_rho/, viral-only arm). Adding Swiss-Prot gives the
    score range on real background without pretending to be a complete index.
    """
    chunks: list[np.ndarray] = []
    kept = dropped = 0
    keep = True
    for faa in faas:
        opener = gzip.open if faa.suffix == ".gz" else open
        with opener(faa, "rt") as fh:
            for line in fh:
                if line[0] == ">":
                    acc = line[1:].split()[0]
                    low = line.lower()
                    keep = not ((exclude_accessions and acc in exclude_accessions)
                                or any(nm in low for nm in exclude_names))
                    kept += keep
                    dropped += not keep
                    if keep:
                        chunks.append(np.zeros(1, dtype=np.uint8))
                elif keep:
                    chunks.append(np.frombuffer(line.strip().upper().encode(), dtype=np.uint8))
    arr = np.concatenate(chunks)
    out = AA_CODE[arr]
    print(f"[index] proteins kept {kept}, dropped {dropped}, residues {out.size}")
    return out


def lineage_accessions(gpff: Path, family: str) -> set[str]:
    """Protein accessions whose taxonomy lineage contains `family`.

    Family-level holdout has to come from taxonomy, not from the species name in
    the FASTA header: a divergent relative of the target family is usually named
    after its host or its sample, so name matching leaves exactly the homologs
    that make the target look known.
    """
    want = family.lower()
    accs: set[str] = set()
    acc = None
    in_tax = False
    tax: list[str] = []
    with gzip.open(gpff, "rt") as fh:
        for line in fh:
            if line.startswith("VERSION"):
                acc = line.split()[1]
            elif line.startswith("  ORGANISM"):
                in_tax, tax = True, []
            elif in_tax:
                if line.startswith("REFERENCE") or line.startswith("COMMENT") or \
                        line.startswith("FEATURES"):
                    if acc and want in " ".join(tax).lower():
                        accs.add(acc)
                    in_tax = False
                else:
                    tax.append(line.strip())
    return accs


def build_index(faa: Path | list[Path], k: int, cache: Path, exclude_family: str | None,
                gpff: Path | None) -> np.ndarray:
    faas = [faa] if isinstance(faa, Path) else list(faa)
    src = "+".join(f.name.split(".")[0] for f in faas)
    tag = f"k{k}_{src}" + (f"_no_{exclude_family.lower()}" if exclude_family else "_all")
    path = cache / f"pepidx_{tag}.npy"
    if path.exists():
        return np.load(path, mmap_mode="r")
    cache.mkdir(parents=True, exist_ok=True)
    excl: set[str] = set()
    if exclude_family:
        if not gpff:
            raise SystemExit("--exclude-family needs --gpff for taxonomy")
        excl = lineage_accessions(gpff, exclude_family)
        print(f"[index] {len(excl)} proteins in {exclude_family} held out")
    codes = proteome_codes(faas, exclude_accessions=excl)
    h = hash_kmers(codes, k)
    del codes
    h = np.unique(h)
    np.save(path, h)
    print(f"[index] {h.size} distinct {k}-mers -> {path}")
    return np.load(path, mmap_mode="r")


class Scorer:
    def __init__(self, index: np.ndarray, k: int):
        self.index = np.asarray(index)
        self.k = k

    def _absent_fraction(self, peptides: list[str]) -> tuple[float, int, int]:
        hs: list[np.ndarray] = []
        for p in peptides:
            if len(p) < self.k:
                continue
            codes = AA_CODE[np.frombuffer(p.encode(), dtype=np.uint8)]
            hs.append(hash_kmers(codes, self.k))
        if not hs:
            return float("nan"), 0, 0
        q = np.concatenate(hs)
        if q.size == 0:
            return float("nan"), 0, 0
        pos = np.searchsorted(self.index, q)
        pos[pos >= self.index.size] = self.index.size - 1
        present = self.index[pos] == q
        return float(1.0 - present.mean()), int(q.size), int(present.sum())

    def score_nt(self, nt: str) -> tuple[float, int, int]:
        """Fraction of a read's six-frame peptide k-mers absent from the index.

        Pooled over all six frames, which is the naive version of the feature and
        is kept only for comparison: five of the six frames are nonsense, so
        their k-mers are absent whatever the read is, and they dominate the
        denominator. Measured on genomes whose own proteins are in the index, the
        pooled score sits at 0.70-0.75 rather than near zero (see
        results/peptide_absence/control_fullindex/). Use score_nt_best_frame.
        """
        return self._absent_fraction(
            [seg for f in six_frames(nt) for seg in f.split("*")])

    def score_nt_best_frame(self, nt: str) -> tuple[float, int, int]:
        """Lowest absent fraction across the six frames: the usable version.

        A read is "known" if *any* one of its frames translates into peptide
        k-mers the index has seen, so the statistic has to be a minimum over
        frames rather than an average. Returns (absent fraction of the best
        frame, k-mers in that frame, frames that could be scored).
        """
        best = (float("nan"), 0, 0)
        frames = 0
        for f in six_frames(nt):
            v, n, _ = self._absent_fraction([seg for seg in f.split("*")])
            if not n:
                continue
            frames += 1
            if n and (best[1] == 0 or v < best[0]):
                best = (v, n, 0)
        return best[0], best[1], frames

    def score_peptide_window(self, aa: str) -> float:
        return self._absent_fraction([seg for seg in aa.split("*")])[0]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    with path.open() as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip().upper())
    if name:
        out.append((name, "".join(buf)))
    return out


def novelty_track(seq: str, scorer: Scorer, win: int, step: int) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    for s in range(0, len(seq) - win + 1, step):
        v, n, _ = scorer.score_nt_best_frame(seq[s:s + win])
        if n:
            xs.append(s + win // 2)
            ys.append(v)
    return np.array(xs), np.array(ys)


def autocorr(x: np.ndarray, y: np.ndarray, lags_nt: list[int], step: int) -> dict:
    """Correlation of the novelty track with itself at a given separation.

    Reported as Pearson r at each lag plus the correlation of the *binarised*
    track (novel versus not, split at the track's own median), because the
    pair-level rules in the feature map are thresholded and it is the
    thresholded correlation that enters them as rho.
    """
    out = {}
    med = float(np.median(y))
    b = (y > med).astype(float)
    for lag in lags_nt:
        d = max(1, round(lag / step))
        if d >= y.size - 2:
            continue
        a1, a2 = y[:-d], y[d:]
        b1, b2 = b[:-d], b[d:]
        def r(u, v):
            if u.std() == 0 or v.std() == 0:
                return None
            return round(float(np.corrcoef(u, v)[0, 1]), 4)
        out[lag] = {"r_continuous": r(a1, a2), "r_thresholded": r(b1, b2), "n": int(a1.size)}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--faa", type=Path, nargs="+",
                    default=[Path("data/prot/viral.1.protein.faa.gz")])
    ap.add_argument("--gpff", type=Path, default=Path("data/prot/viral.1.protein.gpff.gz"))
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--cache", type=Path, default=Path("data/pepidx"))
    ap.add_argument("--genomes", type=Path, nargs="*", default=[],
                    help="FASTA genomes to build novelty tracks for")
    ap.add_argument("--exclude-family", default=None,
                    help="taxonomic family held out of the index (e.g. Rhabdoviridae)")
    ap.add_argument("--win", type=int, default=107, help="track window = read length")
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--lags", type=int, nargs="+",
                    default=[10, 20, 50, 100, 150, 200, 250, 300, 400, 500, 700,
                             1000, 1500, 2000, 3000])
    ap.add_argument("--out", type=Path, default=Path("results/peptide_absence"))
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    idx = build_index(args.faa, args.k, args.cache, args.exclude_family, args.gpff)
    scorer = Scorer(idx, args.k)
    res = {"seed": args.seed, "k": args.k, "index_kmers": int(np.asarray(idx).size),
           "exclude_family": args.exclude_family, "win": args.win, "step": args.step,
           "genomes": {}}

    for g in args.genomes:
        for name, seq in read_fasta(g):
            x, y = novelty_track(seq, scorer, args.win, args.step)
            if y.size == 0:
                continue
            res["genomes"][name.split()[0]] = {
                "file": str(g),
                "length": len(seq),
                "windows": int(y.size),
                "absent_fraction": {
                    "mean": round(float(y.mean()), 4),
                    "median": round(float(np.median(y)), 4),
                    "p05": round(float(np.quantile(y, 0.05)), 4),
                    "p95": round(float(np.quantile(y, 0.95)), 4),
                    "fraction_windows_fully_absent": round(float((y >= 0.999).mean()), 4),
                    "fraction_windows_fully_present": round(float((y <= 0.001).mean()), 4),
                },
                "autocorrelation": autocorr(x, y, args.lags, args.step),
                "track": [[int(a), round(float(b), 4)] for a, b in zip(x, y)],
            }
            print(name.split()[0], res["genomes"][name.split()[0]]["absent_fraction"])
            print("  autocorr", {l: v["r_thresholded"] for l, v in
                                 res["genomes"][name.split()[0]]["autocorrelation"].items()})

    (args.out / f"peptide_absence_s{args.seed}.json").write_text(json.dumps(res, indent=2))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/peptide_absence.py", "git_sha": sha, "seed": args.seed,
        "params": {k: str(v) for k, v in vars(args).items()},
        "wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2))


if __name__ == "__main__":
    main()
