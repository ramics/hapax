"""How far does one spurious peptide derived from ribosomal RNA reach?

The mate-separation thread's cumulative index curve located where ribosomal reads stop
looking absent, and it is not our own assembly. Against a viral-only index P(absent | rRNA)
is 0.998, exactly what an index of real proteins must give for a molecule that encodes no
protein. It falls to 0.473 between the 1% and 10% UniRef50 subsamples, while non-ribosomal
background moves only 0.798 to 0.720 over the same step. Excluding the 388 proteins
predicted from our own ribosomal contigs moves it 0.107 to 0.143, so the mechanism measured
in `prodigal_on_rrna.py` is real but four orders of index size away from where it matters,
and a low-complexity filter changes nothing at all (0.144 against 0.143).

So the question is why a protein database contains anything a ribosomal read can match, and
why adding more of it hurts ribosomal reads several times as much as everything else. The
differential is the part a uniform sprinkling of badly annotated entries cannot explain.

**The hypothesis this measures: amplification by conservation.** UniRef50 includes TrEMBL,
which is automatically annotated, and an automated gene caller over a genome assembly makes
the same mistake on ribosomal loci that pyrodigal makes on our contigs. But ribosomal RNA is
the most conserved sequence there is, so a spurious entry derived from *one* organism's
ribosomal locus carries peptides that ribosomal reads from across the tree of life can
match, while a spurious entry from anywhere else reaches only its own neighbourhood.
Background sequence has no such multiplier. If the multiplier is large, a handful of bad
entries explains the whole differential and the number of them barely matters.

The measurement needs no protein database at all, which is why it belongs on this branch:
six-frame translate curated ribosomal RNA, and for each peptide 10-mer count how many
distinct SILVA sequences, genera and phyla it occurs in. That is exactly the number of
lineages one spurious entry's peptides would reach.

The control shuffles each sequence, preserving length and base composition while destroying
conservation and every shared motif. Reach there should collapse toward one. If real rRNA
does not greatly exceed it, amplification is wrong and the explanation lies elsewhere.

Reach is reported occurrence-weighted as well as per distinct k-mer, because what matters
is the reach of a peptide a read actually draws, not of an average catalogue entry.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
K_AA = 10

_B = np.full(256, 255, dtype=np.uint8)
for ch, code in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
    _B[ord(ch)] = code
    _B[ord(ch.lower())] = code

_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}
_AA = "ACDEFGHIKLMNPQRSTVWY"
# codon index = b0*16 + b1*4 + b2 with A,C,G,T = 0,1,2,3; 255 marks a stop
_TABLE = np.full(64, 255, dtype=np.uint8)
for codon, aa in _CODONS.items():
    i = "ACGT".index(codon[0]) * 16 + "ACGT".index(codon[1]) * 4 + "ACGT".index(codon[2])
    _TABLE[i] = 255 if aa == "*" else _AA.index(aa)


def peptide_codes(bases: np.ndarray) -> np.ndarray:
    """Every stop-free peptide K_AA-mer code over all six frames of a 2-bit base array."""
    out = []
    rc = (np.uint8(3) - bases[::-1])
    rc[bases[::-1] == 255] = 255
    for strand in (bases, rc):
        for off in range(3):
            n = (strand.size - off) // 3
            if n < K_AA:
                continue
            cod = strand[off : off + 3 * n].reshape(n, 3).astype(np.int32)
            bad = (cod == 255).any(axis=1)
            idx = np.where(bad, 0, cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2])
            aa = _TABLE[idx]
            aa[bad] = 255
            m = n - K_AA + 1
            if m <= 0:
                continue
            code = np.zeros(m, dtype=np.uint64)
            ok = np.ones(m, dtype=bool)
            valid = aa != 255
            a64 = np.where(valid, aa, 0).astype(np.uint64)
            for j in range(K_AA):
                code = code * np.uint64(25) + a64[j : j + m]
                ok &= valid[j : j + m]
            out.append(code[ok])
    return np.concatenate(out) if out else np.empty(0, dtype=np.uint64)


def reach(seqs: list[str], labels: list[int], n_groups: int) -> dict:
    """Distribution of how many distinct sequences and groups each peptide K_AA-mer spans."""
    codes, sids, gids = [], [], []
    for i, (s, g) in enumerate(zip(seqs, labels)):
        c = peptide_codes(_B[np.frombuffer(s.encode(), dtype=np.uint8)])
        if c.size:
            codes.append(c)
            sids.append(np.full(c.size, i, dtype=np.int32))
            gids.append(np.full(c.size, g, dtype=np.int32))
    if not codes:
        return {}
    code = np.concatenate(codes)
    sid = np.concatenate(sids)
    gid = np.concatenate(gids)

    def distinct_per_code(partner):
        pairs = np.unique(np.stack([code, partner.astype(np.uint64)], axis=1), axis=0)
        uc, cnt = np.unique(pairs[:, 0], return_counts=True)
        return uc, cnt

    uc, seq_reach = distinct_per_code(sid)
    _, grp_reach = distinct_per_code(gid)

    # occurrence weighting: the reach of a peptide a read actually draws
    order = np.searchsorted(uc, code)
    occ_seq = seq_reach[order]
    occ_grp = grp_reach[order]

    # The distribution matters more than the mean: a few 10-mers spanning the whole tree
    # while the rest span almost nothing is a different and more actionable claim than a
    # uniformly raised multiplier, and it predicts the shape of the matching entries' head.
    occ_per_code = np.bincount(order, minlength=uc.size).astype(np.int64)
    by_reach = np.argsort(-seq_reach)
    cum = np.cumsum(occ_per_code[by_reach])
    total_occ = int(code.size)
    head = {f"top_{n}_kmers_share_of_occurrences":
            (float(cum[min(n, cum.size) - 1] / total_occ) if cum.size else 0.0)
            for n in (1, 10, 100, 1000)}
    nseq = len(seqs)

    return {
        "sequences": len(seqs),
        "peptide_kmer_occurrences": int(code.size),
        "distinct_peptide_kmers": int(uc.size),
        "reach_distribution": {
            "occurrence_weighted_percentiles": {
                f"p{q}": float(np.percentile(occ_seq, q))
                for q in (10, 25, 50, 75, 90, 99, 99.9)
            },
            "per_distinct_kmer_percentiles": {
                f"p{q}": float(np.percentile(seq_reach, q))
                for q in (50, 90, 99, 99.9)
            },
            "distinct_kmers_reaching_at_least": {
                "1pc_of_sequences": float(np.mean(seq_reach >= 0.01 * nseq)),
                "10pc_of_sequences": float(np.mean(seq_reach >= 0.10 * nseq)),
                "50pc_of_sequences": float(np.mean(seq_reach >= 0.50 * nseq)),
            },
            "occurrences_from_kmers_reaching_at_least": {
                "1pc_of_sequences": float(
                    occ_per_code[seq_reach >= 0.01 * nseq].sum() / total_occ),
                "10pc_of_sequences": float(
                    occ_per_code[seq_reach >= 0.10 * nseq].sum() / total_occ),
                "50pc_of_sequences": float(
                    occ_per_code[seq_reach >= 0.50 * nseq].sum() / total_occ),
            },
            "head_share": head,
        },
        "per_distinct_kmer": {
            "mean_sequence_reach": float(seq_reach.mean()),
            "median_sequence_reach": float(np.median(seq_reach)),
            "fraction_reaching_one_sequence": float(np.mean(seq_reach == 1)),
            "p99_sequence_reach": float(np.percentile(seq_reach, 99)),
            "max_sequence_reach": int(seq_reach.max()),
        },
        "occurrence_weighted": {
            "mean_sequence_reach": float(occ_seq.mean()),
            "median_sequence_reach": float(np.median(occ_seq)),
            "fraction_reaching_one_sequence": float(np.mean(occ_seq == 1)),
            "mean_group_reach": float(occ_grp.mean()),
            "median_group_reach": float(np.median(occ_grp)),
            "groups_available": n_groups,
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sequences", type=int, default=8000)
    ap.add_argument("--min-length", type=int, default=1200)
    ap.add_argument("--group-depth", type=int, default=2,
                    help="taxonomy levels from the SILVA header that define a group; 2 is "
                         "roughly phylum, so group reach counts phyla a peptide spans")
    ap.add_argument("--silva", type=Path,
                    default=ROOT / "data" / "silva" / "SSURef_NR99.fasta.gz")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "rrna_content")
    args = ap.parse_args()

    t0 = time.perf_counter()
    rng = np.random.default_rng(args.seed)

    # one pass, reservoir sampled, keeping the taxonomy string
    res: list[tuple[str, str]] = []
    seen = 0
    name, cur = None, []
    with gzip.open(args.silva, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name and cur:
                    seen = _offer(res, args.sequences, seen, name,
                                  "".join(cur), args.min_length, rng)
                name, cur = line[1:].strip(), []
            else:
                cur.append(line.strip())
        if name and cur:
            seen = _offer(res, args.sequences, seen, name, "".join(cur),
                          args.min_length, rng)

    seqs = [s for _, s in res]
    taxa = []
    for h, _ in res:
        parts = h.split(" ", 1)
        lineage = parts[1] if len(parts) > 1 else ""
        taxa.append(";".join(lineage.split(";")[: args.group_depth]))
    uniq = sorted(set(taxa))
    gmap = {t: i for i, t in enumerate(uniq)}
    labels = [gmap[t] for t in taxa]
    print(f"{len(seqs):,} sequences, {len(uniq):,} distinct taxonomy groups at depth "
          f"{args.group_depth}")

    real = reach(seqs, labels, len(uniq))
    shuf = []
    for s in seqs:
        arr = np.frombuffer(s.encode(), dtype=np.uint8).copy()
        rng.shuffle(arr)
        shuf.append(arr.tobytes().decode())
    ctrl = reach(shuf, labels, len(uniq))

    mult = (real["occurrence_weighted"]["mean_sequence_reach"]
            / ctrl["occurrence_weighted"]["mean_sequence_reach"])
    out = {
        "seed": args.seed,
        "k_aa": K_AA,
        "reference": args.silva.name,
        "real_rrna": real,
        "composition_shuffled_control": ctrl,
        "conservation_multiplier": mult,
        "note": "reach is the number of distinct sequences or taxonomy groups a peptide "
                "10-mer occurs in, which is the number of lineages one spurious database "
                "entry derived from a ribosomal locus would match reads from",
    }
    (args.out / f"rrna_peptide_reach_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"rrna_peptide_reach_manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/rrna_peptide_reach.py",
        "git_sha": _git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
    }, indent=2))

    for label, d in (("real rRNA", real), ("shuffled", ctrl)):
        o = d["occurrence_weighted"]
        p = d["per_distinct_kmer"]
        print(f"{label:>10}: occurrence-weighted sequence reach mean "
              f"{o['mean_sequence_reach']:.2f} median {o['median_sequence_reach']:.0f}, "
              f"group reach mean {o['mean_group_reach']:.2f}, "
              f"singletons {o['fraction_reaching_one_sequence']:.4f} "
              f"| per distinct k-mer mean {p['mean_sequence_reach']:.2f} "
              f"max {p['max_sequence_reach']:,}")
    print(f"conservation multiplier: {mult:.1f}x")
    rd = real["reach_distribution"]
    print("  occurrence-weighted reach percentiles: " + ", ".join(
        f"{k} {v:.0f}" for k, v in rd["occurrence_weighted_percentiles"].items()))
    print("  share of occurrences from the top N distinct 10-mers by reach: " + ", ".join(
        f"{k.split('_')[1]} {v:.4f}" for k, v in rd["head_share"].items()))
    print("  distinct 10-mers reaching >=1/10/50% of sequences: " + ", ".join(
        f"{v:.5f}" for v in rd["distinct_kmers_reaching_at_least"].values())
        + "; occurrences they carry: " + ", ".join(
        f"{v:.4f}" for v in rd["occurrences_from_kmers_reaching_at_least"].values()))


def _offer(res, n, seen, name, seq, min_len, rng) -> int:
    seq = seq.upper().replace("U", "T")
    if len(seq) < min_len:
        return seen
    if len(res) < n:
        res.append((name, seq))
    else:
        j = int(rng.integers(0, seen + 1))
        if j < n:
            res[j] = (name, seq)
    return seen + 1


def _git_sha() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=ROOT) != 0
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
