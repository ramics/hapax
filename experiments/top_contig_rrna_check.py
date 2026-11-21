"""Is the deepest unnameable sequence ribosomal that SILVA nearly caught, or not ribosomal at all?

WHY THIS EXISTS. `results/residue_contig_annotation/` finds the 200 highest-mass
unidentifiable contigs markedly LESS coding than length-matched unidentifiable ones (coding
density 0.604 against 0.818, half-coding 61% against 90%, Pfam 7.5% against 11.2%) while
leading every pile on Rfam at 3.0%, those hits being LSU and SSU rRNA. The reading offered
is that the pile is partly ribosomal sequence SILVA missed.

**That reading rests on 6 contigs of 200 against 2 in the control, and it should not be
believed on 6 contigs.** There is a direct test, and it costs one alignment: the contig
split called a contig ribosomal at >= 100 matching bases against SILVA SSU+LSU, so if these
contigs are ribosomal-but-missed they should pile up just under that line. If they sit at
zero they are not ribosomal at all, and the Rfam hits are either a profile reaching past
sequence identity on a handful of contigs or chance at the pile's own rate.

The two outcomes say different things about the tool. Just under the line means the operating
point is doing the damage and a second ribosomal instrument belongs in the cascade. At zero
means the deepest unnameable sequence is genuinely non-ribosomal non-coding material, and
"ranking by abundance finds rRNA" would be the wrong lesson to layer in.

CONTROLS, same code, same pass, from the same contig sets the annotation used:
  - `unidentified_lenmatched`: unidentifiable contigs drawn in the same 100 nt length bins,
    query names excluded, which is the comparison the annotation's own control makes;
  - `identified_norrna`: identified non-rRNA contigs, the known-positive for "this pile is
    ordinary non-ribosomal sequence";
  - `rrna_contigs`: the contigs the split DID call ribosomal, the known-positive for the
    quantity being measured. Without this arm a low number is unreadable, because nothing
    says what a true ribosomal contig scores on this instrument.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance  # noqa: E402


def read_fa(p: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name:
        out[name] = "".join(buf)
    return out


def best_mlen(al, seq: str) -> int:
    return max((h.mlen for h in al.map(seq)), default=0)


def summarise(v, thresholds=(40, 60, 100, 150)) -> dict:
    a = np.asarray(v, dtype=float)
    d = {"n": int(a.size), "p50": float(np.median(a)) if a.size else None,
         "p90": float(np.quantile(a, 0.9)) if a.size else None, "mean": float(a.mean()) if a.size else None,
         "frac_0": float((a == 0).mean()) if a.size else None}
    for t in thresholds:
        d[f"frac_ge_{t}"] = float((a >= t).mean()) if a.size else None
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--top", type=Path, default=ROOT / "results/residue_census/top_unidentifiable_contigs.tsv")
    ap.add_argument("--query-fa", type=Path, default=None,
                    help="score this FASTA instead of the top-mass contig list, e.g. the "
                         "residue-only assembly's new contigs. The arms are unchanged, so the "
                         "known-positive rRNA-contig reading stays comparable across queries.")
    ap.add_argument("--tag", default="top_mass_unidentifiable", help="name for the query pile")
    ap.add_argument("--preset", default="sr",
                    help="minimap2 preset. `sr` is the cascade's own operating point and "
                         "answers 'did the split's threshold nearly catch this'. It CANNOT "
                         "answer 'is this ribosomal at all': sr is tuned for short reads at "
                         "high identity, so divergent rRNA is invisible to it. Run `asm20` as "
                         "a second arm for that question - a covariance model scores structure "
                         "as well as sequence and reaches further than any identity aligner, so "
                         "a negative at sr alone must not be reported as 'not ribosomal'.")
    ap.add_argument("--own-all", type=Path, default=ROOT / "data/work/sample_contigs_all.fa")
    ap.add_argument("--sets", type=Path, default=ROOT / "data/geo/sets")
    ap.add_argument("--silva-mmi", type=Path, default=ROOT / "data/work/silva.mmi")
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_census")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    import mappy

    own = read_fa(args.own_all)
    if args.query_fa:
        query = read_fa(args.query_fa)
    else:
        names = [l.split("\t")[0] for l in open(args.top).read().strip().split("\n")[1:]]
        query = {n: own[n] for n in names if n in own}
    targets = [len(s) for s in query.values()]
    pools = {"unidentified_lenmatched": read_fa(args.sets / "sample_contigs_unidentified.fa"),
             "resasm_rederived": read_fa(ROOT / "data/work/resasm/resasm_rederived.fa"),
             "identified_norrna": read_fa(args.sets / "sample_contigs_identified_norrna.fa"),
             "rrna_contigs": read_fa(args.sets / "sample_contigs_rrna.fa")}

    al = (mappy.Aligner(str(args.silva_mmi), preset=args.preset) if args.preset == "sr"
          else mappy.Aligner(str(ROOT / "data/work/silva.fa"), preset=args.preset))
    if not al:
        raise SystemExit("could not open the SILVA index")

    piles = {args.tag: query}
    for tag, pool in pools.items():
        bins: dict[int, list[str]] = {}
        for k, s in pool.items():
            if k not in query:
                bins.setdefault(len(s) // 100, []).append(k)
        drawn, short = {}, 0
        for i, L in enumerate(targets):
            cand = bins.get(L // 100) or []
            if not cand:
                short += 1
                continue
            k = rng.choice(cand)
            drawn[f"{i}_{k}"] = pool[k]
        piles[tag] = drawn
        if short:
            piles.setdefault("_short", {})[tag] = short

    shortfall = piles.pop("_short", {})
    out = {"seed": args.seed, "preset": args.preset,
           "operating_point": "minimap2 sr against SILVA 138.2 SSU+LSU NR99; the contig split "
                              "called a contig ribosomal at >=100 matching bases",
           "query_contigs_excluded_from_control_pools": len(query),
           "length_bins_unfilled": shortfall,
           "piles": {k: summarise([best_mlen(al, s) for s in v.values()]) for k, v in piles.items()}}
    out["manifest"] = {"script": "experiments/top_contig_rrna_check.py", **provenance(__file__),
                       "parameters": {k: str(v) for k, v in vars(args).items()}}
    suffix = "" if args.preset == "sr" else f"_{args.preset}"
    (args.out / f"{args.tag}_rrna{suffix}_s{args.seed}.json").write_text(json.dumps(out, indent=1))
    for k, v in out["piles"].items():
        print(f"{k:28s} n={v['n']:4d} p50={v['p50']:6.0f} mean={v['mean']:8.1f} "
              f"zero={v['frac_0']*100:5.1f}% >=40 {v['frac_ge_40']*100:5.1f}% "
              f">=100 {v['frac_ge_100']*100:5.1f}% >=150 {v['frac_ge_150']*100:5.1f}%")
    if shortfall:
        print("unfilled length bins:", shortfall)


if __name__ == "__main__":
    main()
