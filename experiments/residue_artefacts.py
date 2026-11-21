"""Library-construction artefacts in the residue, measured against identified reads.

WHY THIS EXISTS. If almost everything real hits something, whatever matches nothing should
be enriched for artefact - a pile defined by matching nothing collects, by construction,
anything a library prep can make that no organism carries. The cascade's artefact tier only
looked for the DNA-seq kind (adapters, primer dimers, low entropy, homopolymers) and found
0.31% of the library, and `results/absence_residual/` found adapters at 0.0004 in the
residue against 0.0001 elsewhere and flat entropy, so "artefacts are not the story" was the
standing reading. But none of those looks for what an RNA-seq library actually makes:

- **fold-backs**: reverse transcriptase switching strand onto its own product makes a read
  that contains its own reverse complement (a hairpin). Real RNA has inverted repeats too
  (stems), so this is measured as a continuous longest-inverted-repeat length and read
  against identified reads and against shuffles, never thresholded alone.
- **same-strand mates**: in a normal pair mate 2 is the reverse complement of the far end
  of the fragment, so where mates overlap they overlap as reverse complements. Mates that
  share sequence in the SAME orientation cannot come from one linear fragment - that is a
  template switch or an inverted chimaera.
- **insert shorter than a read** (mate overlap starting before mate 1's 3' end), which in
  a trimmed library is the residue of adapter read-through.
- **poly(A)/poly(T)** runs, the oligo-dT / polyadenylation signature, which a database of
  genomes does not carry because a tail is not templated.
- **vector and PhiX** (UniVec_Core, NC_001422.1), the spike-in and cloning contaminants
  that Kraken2's genome indexes do not hold.

POLARITY, stated before the numbers: a field that is higher in the residue than in reads the
cascade identified points toward ARTEFACT only if it is also above the shuffled null. A
field that is higher in identified reads than in the residue is a property of known
biology (rRNA stems are inverted repeats), and is reported as that.

CONTROLS, in the same pass and through the same code:
  - `identified`: the cascade's subtracted pairs for the same seed (`subtracted_s{seed}.fa`),
    split further into rRNA and non-rRNA by the library label table, because rRNA's
    structure would otherwise dominate the fold-back comparison.
  - `shuffled`: every residue read mononucleotide-shuffled in place (keeps length and GC,
    destroys structure). Bias: a mononucleotide shuffle also destroys real dinucleotide
    bias, so its inverted-repeat length is slightly LOW; it is a null, not a baseline.

Mates are reported separately (CLAUDE.md: R2 is 93 nt against R1's 107 after the 14 nt
trim, and pair-symmetric fields average two opposite biases).
"""
from __future__ import annotations

import argparse
import json
import random
import statistics as st
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_pairs_fa(path: Path):
    out, name = {}, None
    with open(path) as h:
        for line in h:
            if line[0] == ">":
                name, mate = line[1:].strip().rsplit("/", 1)
            else:
                out.setdefault(name, ["", ""])[int(mate) - 1] = line.strip()
    return out


def longest_inverted_repeat(s: str, k: int = 10) -> int:
    """Longest exact stretch of `s` whose reverse complement also occurs in `s`, excluding
    a stretch matching itself (a palindrome is its own reverse complement and would score
    its full length trivially). Exact k-mer runs, so it is a lower bound on a stem with
    mismatches - the same bound for every pile, which is all the comparison needs."""
    if len(s) < k:
        return 0
    pos = {}
    for i in range(len(s) - k + 1):
        pos.setdefault(s[i:i + k], []).append(i)
    r = rc(s)
    n = len(s)
    best = 0
    # s[i:i+k] == r[j:j+k] is a direct match between s and its reverse complement, so
    # one inverted repeat is a run of consecutive (i, j), (i+1, j+1) on the diagonal i - j.
    runs = {}
    for j in range(n - k + 1):
        km = r[j:j + k]
        for i in pos.get(km, ()):
            src = n - k - j  # where this k-mer's reverse complement sits in s
            if abs(src - i) < k:  # the k-mer overlaps its own partner: palindrome core
                continue
            d = i - j
            prev = runs.get((d, j - 1))
            length = (prev + 1) if prev else 1
            runs[(d, j)] = length
            best = max(best, length)
    return best + k - 1 if best else 0


def mate_overlap(s1: str, s2: str, k: int = 15):
    """(normal, same_strand, insert). normal: mate 1 shares a k-mer with rc(mate 2), which
    is how an insert shorter than the two reads looks. same_strand: mate 1 shares a k-mer
    with mate 2 as written, which one linear fragment cannot produce. insert: for a normal
    overlap, the fragment's far end in mate 1 coordinates (mate 2's 5' end after its 14 nt
    trim), else None. The library is adapter-trimmed already, so read-through shows as
    short reads and short inserts, not as adapter sequence."""
    if len(s1) < k or len(s2) < k:
        return False, False, False
    k1 = {s1[i:i + k]: i for i in range(len(s1) - k + 1)}
    r2 = rc(s2)
    normal, short = False, None
    for j in range(len(r2) - k + 1):
        i = k1.get(r2[j:j + k])
        if i is not None:
            normal = True
            # r2 position j aligns to s1 position i, so r2's end (mate 2's 5' end, the
            # fragment's far end) sits at s1 coordinate i - j + len(r2).
            frag_end = i - j + len(r2)
            short = frag_end
            break
    if any(s2[j:j + k] in k1 for j in range(len(s2) - k + 1)):
        # a k-mer that is its own reverse complement can sit on both strands; ignore
        # those so a palindrome cannot fake a same-strand mate
        same = any(s2[j:j + k] in k1 and s2[j:j + k] != rc(s2[j:j + k])
                   for j in range(len(s2) - k + 1))
    else:
        same = False
    return normal, same, short


def poly_run(s: str) -> int:
    best = cur = 0
    prev = ""
    for c in s:
        cur = cur + 1 if c == prev and c in "AT" else (1 if c in "AT" else 0)
        prev = c
        best = max(best, cur)
    return best


def shuffle(s: str, rng: random.Random) -> str:
    x = list(s)
    rng.shuffle(x)
    return "".join(x)


def summarise(vals):
    a = np.asarray(vals, dtype=float)
    return {"n": int(a.size), "mean": float(a.mean()) if a.size else None,
            "p50": float(np.median(a)) if a.size else None,
            "p90": float(np.quantile(a, 0.9)) if a.size else None}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--recut", type=Path, default=ROOT / "data/recut")
    ap.add_argument("--labels", type=Path, default=ROOT / "data/work/library_labels.npz")
    ap.add_argument("--univec", type=Path, default=ROOT / "data/univec/UniVec_Core")
    ap.add_argument("--phix", type=Path, default=ROOT / "data/univec/phix.fa")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    ap.add_argument("--fold-thresholds", type=int, nargs="+", default=[14, 18, 24, 30])
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    import mappy
    t0 = time.perf_counter()

    res = read_pairs_fa(args.recut / f"residue_s{args.seed}.fa")
    known = read_pairs_fa(args.recut / f"subtracted_s{args.seed}.fa")

    # Split identified pairs into rRNA and not by the whole-library label table, at the
    # cascade's operating point (>=40 matching bases on SILVA, either mate).
    d = np.load(args.labels)
    spot = d["spot"]
    rr = np.maximum(d["silva1"], d["silva2"]) >= 40
    is_rrna = dict(zip(spot.tolist(), rr.tolist()))
    piles = {"residue": res,
             "identified_rrna": {k: v for k, v in known.items()
                                 if is_rrna.get(int(k.split(".")[1]), False)},
             "identified_other": {k: v for k, v in known.items()
                                  if not is_rrna.get(int(k.split(".")[1]), False)}}
    piles["residue_shuffled"] = {k: [shuffle(a, rng), shuffle(b, rng)] for k, (a, b) in res.items()}
    # Stated check: every residue pair must read below the operating point on SILVA,
    # or the label table is a different instrument from the cascade.
    bad = sum(1 for k in res if is_rrna.get(int(k.split(".")[1]), False))
    if bad:
        raise SystemExit(f"{bad} residue pairs read as rRNA in the label table: instrument mismatch")

    vec = mappy.Aligner(str(args.univec), preset="sr")
    phx = mappy.Aligner(str(args.phix), preset="sr")

    def best(al, s):
        return max((h.mlen for h in al.map(s)), default=0)

    out = {"seed": args.seed, "piles": {}}
    for pname, pile in piles.items():
        per_mate = {}
        for m in (0, 1):
            fold = [longest_inverted_repeat(v[m]) for v in pile.values()]
            poly = [poly_run(v[m]) for v in pile.values()]
            vecm = [best(vec, v[m]) for v in pile.values()]
            phxm = [best(phx, v[m]) for v in pile.values()]
            per_mate[f"mate{m + 1}"] = {
                "len": summarise([len(v[m]) for v in pile.values()]),
                "foldback_len": summarise(fold),
                "foldback_ge": {str(t): float(np.mean(np.asarray(fold) >= t)) for t in args.fold_thresholds},
                "polyAT_run": summarise(poly),
                "polyAT_ge12": float(np.mean(np.asarray(poly) >= 12)),
                "univec_ge30": float(np.mean(np.asarray(vecm) >= 30)),
                "phix_ge30": float(np.mean(np.asarray(phxm) >= 30)),
            }
        ov = [mate_overlap(a, b) for a, b in pile.values()]
        out["piles"][pname] = {
            "pairs": len(pile), **per_mate,
            "mates_overlap_normal": float(np.mean([o[0] for o in ov])),
            "mates_same_strand": float(np.mean([o[1] for o in ov])),
            "insert_from_overlap": summarise([o[2] for o in ov if o[2] is not None]),
            "insert_le_107": float(np.mean([o[2] is not None and o[2] <= 107 for o in ov])),
            "mate1_shorter_than_100": float(np.mean([len(a) < 100 for a, _ in pile.values()])),
        }
    out["seconds"] = time.perf_counter() - t0
    out["manifest"] = {"script": "experiments/residue_artefacts.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()},
                       "operating_points": {"rrna_split": "silva >=40 matching bases either mate",
                                            "vector_phix": ">=30 matching bases, minimap2 sr",
                                            "foldback": "exact 10-mer runs, palindrome cores excluded",
                                            "mate_overlap": "exact 15-mer, R2 trimmed 14"}}
    (args.out / f"artefacts_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    for p, v in out["piles"].items():
        print(p, v["pairs"], "fold>=18 m1/m2", round(v["mate1"]["foldback_ge"]["18"], 4),
              round(v["mate2"]["foldback_ge"]["18"], 4), "same-strand", round(v["mates_same_strand"], 4),
              "normal", round(v["mates_overlap_normal"], 3), "ins<=107", round(v["insert_le_107"], 4), "m1<100", round(v["mate1_shorter_than_100"], 4),
              "polyAT12 m1", round(v["mate1"]["polyAT_ge12"], 4), "univec m1", round(v["mate1"]["univec_ge30"], 4),
              "phix", round(v["mate1"]["phix_ge30"], 4))


if __name__ == "__main__":
    main()
