"""When half a read matches well and the rest does not, what is the rest? (Imo, 1.ii)

THE OBSERVATION THIS ANSWERS. Against NCBI core_nt, residue reads do hit - 87% of them -
but over about a third of their length at ~89% identity, where a ribosomal read aligns end
to end at 100%. Partial alignment is therefore the NORMAL case for the residue, not an
anomaly, and "is it a chimaera?" is one of at least four explanations:

  A. CHIMAERA. The two parts come from different molecules, joined in the library or in
     the assembly. Tested by asking whether the unaligned part aligns somewhere else.
  B. DIVERGENT CONTINUATION. One molecule, whose conserved part is recognisable and whose
     variable part is not. Tested by asking whether the unaligned part still reaches the
     SAME protein at a threshold low enough to see it - the sharpest discriminator here,
     because a chimaera's second half has no reason to return to the first half's subject.
  C. ARTEFACT END. Adapter, poly-A/G, or low-complexity sequence in the unaligned part.
     Tested by entropy, adapter match and homopolymer run.
  D. SHORT TRUE HOMOLOGY. A genuinely short shared region - a conserved motif - with
     ordinary, unexplainable sequence either side.

WHY EVERY TEST IS RUN ON KNOWN READS TOO. Each test has a false-positive rate at 107 nt
and that rate is the whole answer to "which can you tell apart". This project has already
had a chimaera test fire on 26% of reads the cascade had identified, which makes it
useless as a per-read call; the same test is run here against the same control so the
number is in front of the reader rather than in a footnote. A category is only reported as
discriminating if it fires much more often on partial residue reads than on the known
control, and where it does not, this script says so.

The unit throughout is the single read, not the pair: the question is about one read's two
halves.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics as st
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import _mm_index, dir_bytes, provenance, shannon_entropy, spread

ROOT = Path(__file__).resolve().parents[1]

ADAPTERS = {
    "truseq_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "truseq_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "nextera": "CTGTCTCTTATACACATCT",
}


def load_reads(path: Path, limit: int, min_len: int = 0, seed: int = 1) -> dict[str, str]:
    """Any fasta: the per-mate files written by subtract_known.py, or a contig file.

    Sampled at random rather than taken from the top of the file, because an assembler
    writes contigs in length order and the first N would be the longest N - which is
    exactly the population whose partial alignments are easiest to explain.
    """
    import random
    seqs, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name and len("".join(buf)) >= min_len:
                seqs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name and len("".join(buf)) >= min_len:
        seqs.append((name, "".join(buf)))
    rng = random.Random(seed)
    rng.shuffle(seqs)
    return dict(seqs[:limit])


def best_nt(idx, seq: str, exclude: str | None = None):
    """Best nucleotide alignment: (query start, query end, matching bases, ref, ref pos).

    `exclude` drops alignments to the query's own entry. Without it the contig run is
    meaningless: the contigs ARE the reference, so every unaligned segment aligns
    perfectly to itself and the chimaera test fires on 100% of both populations - which
    is exactly what the first run of this script reported, and why it is worth leaving
    this note here rather than silently fixing it.
    """
    best = None
    for h in idx.map(seq):
        if exclude is not None and h.ctg == exclude:
            continue
        if best is None or h.mlen > best[2]:
            best = (h.q_st, h.q_en, h.mlen, h.ctg, h.r_st)
    return best


def diamond(query: dict[str, str], db: Path, binary: Path, work: Path, threads: int,
            evalue: float, mode: str = "--very-sensitive"):
    """blastx, keeping query start/end so the ALIGNED SPAN is available, not just a score."""
    fa = work / "q.fa"
    with open(fa, "w") as h:
        for k, v in query.items():
            h.write(f">{k}\n{v}\n")
    tsv = work / "q.tsv"
    subprocess.run([str(binary), "blastx", "--db", str(db), "--query", str(fa),
                    "--out", str(tsv), "--threads", str(threads), "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "bitscore", "qstart", "qend",
                    "--max-target-seqs", "5", "--evalue", str(evalue), mode, "--quiet"],
                   check=True)
    hits: dict[str, list] = {}
    with open(tsv) as h:
        for line in h:
            q, sid, pid, bits, qs, qe = line.rstrip("\n").split("\t")
            hits.setdefault(q, []).append(
                {"subject": sid, "pident": float(pid), "bitscore": float(bits),
                 "qstart": min(int(qs), int(qe)), "qend": max(int(qs), int(qe))})
    tsv.unlink(missing_ok=True)
    return hits


def artefact_end(seg: str) -> dict:
    ent = shannon_entropy(seg, 3) if len(seg) >= 12 else None
    longest_run = 1
    run = 1
    for i in range(1, len(seg)):
        run = run + 1 if seg[i] == seg[i - 1] else 1
        longest_run = max(longest_run, run)
    adapter = any(a[:12] in seg or seg in a for a in ADAPTERS.values())
    return {"entropy3": ent, "longest_homopolymer": longest_run, "adapter": adapter,
            "fires": bool(adapter or (ent is not None and ent < 0.6) or longest_run >= 12)}


def classify(seq: str, span: tuple[int, int], nt_idx, prot_hits: list, seg_hits: list,
             min_seg: int, exclude: str | None = None) -> dict:
    """Which of A-D fire for one read with a partial alignment."""
    qs, qe = span
    left, right = seq[:qs], seq[qe:]
    seg = left if len(left) >= len(right) else right
    out = {"aligned_fraction": (qe - qs) / len(seq), "unaligned_len": len(seg)}
    if len(seg) < min_seg:
        out["verdict"] = "unaligned part too short to test"
        return out
    # A: the unaligned part aligns somewhere else
    other = best_nt(nt_idx, seg, exclude) if nt_idx else None
    out["A_chimaera"] = bool(other and other[2] >= min_seg)
    out["A_other_ref"] = other[3] if other else None
    # B: the unaligned part reaches the same protein as the aligned part
    first_subjects = {h["subject"] for h in prot_hits}
    out["B_divergent_continuation"] = any(h["subject"] in first_subjects
                                          for h in seg_hits)
    out["B_best_seg_bitscore"] = max((h["bitscore"] for h in seg_hits), default=0.0)
    # C: artefact end
    art = artefact_end(seg)
    out["C_artefact_end"] = art["fires"]
    out["C_detail"] = art
    # D: nothing fires and the aligned part is short
    out["D_short_true_homology"] = not (out["A_chimaera"] or
                                        out["B_divergent_continuation"] or
                                        out["C_artefact_end"])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--residue", type=Path,
                    default=Path("/tmp/claude-0/wk/recut/residue_s1.fa"))
    ap.add_argument("--known", type=Path,
                    default=Path("/tmp/claude-0/wk/recut/subtracted_s1.fa"),
                    help="reads the cascade DID identify: the control every test is "
                         "scored against, because a test that fires here is not "
                         "measuring novelty")
    ap.add_argument("--max-reads", type=int, default=3000)
    ap.add_argument("--min-len", type=int, default=0,
                    help="length floor, used when the input is contigs rather than reads")
    ap.add_argument("--refs", type=Path, nargs="+",
                    default=[ROOT / "data/geo/sample_contigs.fa",
                             ROOT / "data/geo/all_contigs_1kb.fa"])
    ap.add_argument("--dmnd", type=Path,
                    default=ROOT / "data/prot/uniref50_10pct.dmnd")
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--evalue", type=float, default=10.0)
    ap.add_argument("--seg-evalue", type=float, default=1000.0,
                    help="the unaligned segment is searched at a deliberately permissive "
                         "threshold: the question is whether it reaches the SAME subject "
                         "at all, not whether that hit would be significant alone")
    ap.add_argument("--partial-lo", type=float, default=0.15)
    ap.add_argument("--partial-hi", type=float, default=0.75)
    ap.add_argument("--min-seg", type=int, default=30)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/anatomy"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/partial_alignment")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    nt_idx = None
    for ref in args.refs:
        if ref.exists():
            nt_idx = _mm_index(ref, args.work / f"{ref.stem}.mmi", args.threads)
            break

    result = {}
    for label, path in (("residue", args.residue), ("known", args.known)):
        reads = load_reads(path, args.max_reads, args.min_len, args.seed)
        prot = diamond(reads, args.dmnd, args.diamond_bin, args.work, args.threads,
                       args.evalue)
        # the span of the best protein hit defines "the part that matches"
        rows, segs = {}, {}
        for name, seq in reads.items():
            hits = sorted(prot.get(name, []), key=lambda h: -h["bitscore"])
            if not hits:
                rows[name] = {"status": "no protein hit at all"}
                continue
            top = hits[0]
            span = (top["qstart"] - 1, top["qend"])
            frac = (span[1] - span[0]) / len(seq)
            rows[name] = {"status": "hit", "aligned_fraction": frac,
                          "pident": top["pident"], "bitscore": top["bitscore"],
                          "subject": top["subject"], "span": span}
            if args.partial_lo <= frac <= args.partial_hi:
                left, right = seq[:span[0]], seq[span[1]:]
                seg = left if len(left) >= len(right) else right
                if len(seg) >= args.min_seg:
                    segs[name] = seg
        seg_hits = diamond(segs, args.dmnd, args.diamond_bin, args.work, args.threads,
                           args.seg_evalue) if segs else {}
        cats = Counter()
        detail = []
        for name, seg in segs.items():
            r = rows[name]
            c = classify(reads[name], r["span"], nt_idx,
                         [{"subject": r["subject"]}], seg_hits.get(name, []),
                         args.min_seg, exclude=name.split()[0])
            for k in ("A_chimaera", "B_divergent_continuation", "C_artefact_end",
                      "D_short_true_homology"):
                if c.get(k):
                    cats[k] += 1
            detail.append({"read": name, **{k: v for k, v in c.items()
                                            if k != "C_detail"}})
        hit_rows = [r for r in rows.values() if r.get("status") == "hit"]
        result[label] = {
            "reads": len(reads),
            "with_protein_hit": len(hit_rows),
            "aligned_fraction": spread([r["aligned_fraction"] for r in hit_rows])
            if hit_rows else None,
            "pident": spread([r["pident"] for r in hit_rows]) if hit_rows else None,
            "partial_reads": len(segs),
            "partial_frac_of_hits": len(segs) / len(hit_rows) if hit_rows else None,
            "categories": {k: v / len(segs) for k, v in cats.items()} if segs else {},
            "category_counts": dict(cats),
        }
        (args.out / f"{label}_detail_s{args.seed}.jsonl").write_text(
            "\n".join(json.dumps(d) for d in detail) + "\n")

    # the comparison that decides what is discriminating
    verdict = {}
    for k in ("A_chimaera", "B_divergent_continuation", "C_artefact_end",
              "D_short_true_homology"):
        r = result["residue"]["categories"].get(k, 0.0)
        n = result["known"]["categories"].get(k, 0.0)
        verdict[k] = {"residue": r, "known_control": n,
                      "ratio": (r / n) if n else None,
                      "discriminating": bool(n and r / n >= 2) or (r > 0.05 and n == 0)}

    out = {"seed": args.seed, "by_set": result, "verdict": verdict,
           "thresholds": {"partial_lo": args.partial_lo, "partial_hi": args.partial_hi,
                          "min_seg": args.min_seg, "seg_evalue": args.seg_evalue},
           "wall_seconds": time.perf_counter() - t0}
    (args.out / f"anatomy_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/partial_alignment_anatomy.py",
        **provenance(__file__),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "index_bytes": {str(args.dmnd): dir_bytes(args.dmnd)},
        "wall_seconds": time.perf_counter() - t0,
        "python": platform.python_version(), "platform": platform.platform(),
    }, indent=2))

    for label, r in result.items():
        print(f"{label}: {r['with_protein_hit']}/{r['reads']} with a protein hit, "
              f"aligned fraction median "
              f"{(r['aligned_fraction'] or {}).get('mean', 0):.2f}, partial "
              f"{r['partial_reads']} ({(r['partial_frac_of_hits'] or 0) * 100:.1f}% of hits)")
        for k, v in r["categories"].items():
            print(f"    {k:<28} {v * 100:5.1f}%")
    print("verdict:", json.dumps(verdict, indent=1))


if __name__ == "__main__":
    main()
