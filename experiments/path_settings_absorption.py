#!/usr/bin/env python3
"""Which contigs stop being unidentifiable when the bubble-merging floor is lowered.

WHY THIS EXISTS. `results/path_settings/` measures, paired by seed, that `lax_bubble`
(`--merge-level 20,0.90`) takes 1.33 points of the library out of the residue against a
same-setting null of 0.01, and that essentially all of it comes out of tier 2 - reads that
mapped only to contigs nothing could identify now map to contigs something can. That is
Imo's objection made quantitative, and it is invisible both to a residue-only assembly and
to the ground-truth panel, where no genome loses breadth or purity in that arm.

But "reads moved between piles" does not say WHAT moved, and the direction matters: a lax
floor could be completing known organisms' contigs, which is a legitimate gain, or welding
unidentifiable sequence onto recognisable sequence so that the whole thing scores above the
protein tier's threshold, which is the failure this project is built around. This script
asks the contigs directly.

THE MEASUREMENT. Take the arm's identified contigs and the BASELINE's unidentifiable
contigs, align the first against the second, and report how much baseline-unidentifiable
sequence now sits inside something the cascade calls identified. Mass is weighted two ways:
by bases, and by bases times MEGAHIT's coverage header, because 10 kb of contig at 300x is a
different event from 10 kb at 2x and the read-level cascade is weighted the second way.
STATED BIAS: the coverage header understates mapped depth about fourfold
(`results/divergent_contigs/`), but it understates it for both arms, and this is a ratio
between arms rather than an absolute depth.

THE CONTROL IS THE WHOLE POINT, AND IT IS THE SAME CODE ON THE SAME PILES.
`baseline_repeat` is the baseline setting run again at a different thread count. Whatever
this measurement reads there is what it reads when nothing has changed, so an arm's number
means something only to the extent that it exceeds it.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from path_settings import merge, provenance, read_fasta  # noqa: E402


def cov_of(name: str) -> float:
    """MEGAHIT's multi= rewritten into the SPAdes-style name by path_settings_run.sh."""
    try:
        return float(name.split("_cov_")[1].split("_")[0])
    except (IndexError, ValueError):
        return 1.0


def len_of(name: str, fallback: int) -> int:
    try:
        return int(name.split("_length_")[1].split("_")[0])
    except (IndexError, ValueError):
        return fallback


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded; this stage is deterministic given the two contig sets")
    ap.add_argument("--baseline-unidentified", type=Path, required=True)
    ap.add_argument("--arms", nargs="+", required=True,
                    help="name=path_to_that_arm's_identified_contigs, repeated")
    ap.add_argument("--minimap2", type=Path, required=True)
    ap.add_argument("--min-identity", type=float, default=0.99,
                    help="the two assemblies are of the SAME reads, so sequence carried "
                         "across is near-identical; this is deliberately strict")
    ap.add_argument("--min-aln", type=int, default=200)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    t0 = time.perf_counter()

    base = read_fasta(args.baseline_unidentified)
    base_len = {n: len_of(n, len(s)) for n, s in base.items()}
    base_cov = {n: cov_of(n) for n in base}
    total_bases = sum(base_len.values())
    total_mass = sum(base_len[n] * base_cov[n] for n in base)

    rows = {}
    for spec in args.arms:
        name, _, path = spec.partition("=")
        cmd = [str(args.minimap2), "-t", str(args.threads), "-cx", "asm5",
               "--secondary=yes", "-N", "50", "-p", "0.01",
               str(args.baseline_unidentified), str(path)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"minimap2 failed for {name}: {proc.stderr[-2000:]}")
        spans: dict[str, list[tuple[int, int]]] = {}
        for line in proc.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            matches, blocklen = int(f[9]), int(f[10])
            if blocklen < args.min_aln or matches / max(1, blocklen) < args.min_identity:
                continue
            spans.setdefault(f[5], []).append((int(f[7]), int(f[8])))
        covered_bases = covered_mass = 0
        n_half = n_ninety = 0
        detail = []
        for n, sp in spans.items():
            cov_nt = merge(sp)
            frac = cov_nt / base_len[n]
            covered_bases += cov_nt
            covered_mass += cov_nt * base_cov[n]
            n_half += frac >= 0.5
            n_ninety += frac >= 0.9
            if frac >= 0.5:
                detail.append({"contig": n, "length": base_len[n], "cov": base_cov[n],
                               "covered_frac": frac})
        detail.sort(key=lambda r: -r["length"] * r["cov"])
        # detail holds only contigs at >=0.5 covered; top1 share is computed over it, so
        # it is an upper bound on how concentrated the mass is, not an estimate.
        rows[name] = {
            "baseline_unidentified_contigs": len(base),
            "baseline_unidentified_bases": total_bases,
            "contigs_touched": len(spans),
            "contigs_half_covered": n_half, "contigs_ninety_covered": n_ninety,
            "bases_absorbed": covered_bases,
            "absorbed_coverage_mass": covered_mass,
            "baseline_unidentified_coverage_mass": total_mass,
            # How much of the absorbed coverage-weighted mass is ONE contig. Measured
            # because it turned out to be nearly all of it in three arms - a 940 nt contig
            # at 6,972x and a 506 nt one at 15,975x - so the weighted percentage without
            # this column reads as a broad effect when it is one deep fragment.
            "top1_share_of_absorbed_mass": (
                max((r["length"] * r["covered_frac"] * r["cov"] for r in detail), default=0)
                / covered_mass
                if covered_mass else None),
            "bases_absorbed_frac": covered_bases / total_bases if total_bases else None,
            "coverage_weighted_mass_frac": covered_mass / total_mass if total_mass else None,
            "top_absorbed": detail[:25],
        }
        print(f"{name}: {n_ninety} of {len(base)} baseline-unidentifiable contigs are "
              f">=90% inside an identified contig; "
              f"{covered_bases / total_bases * 100:.2f}% of bases, "
              f"{covered_mass / total_mass * 100:.2f}% coverage-weighted", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"absorption_s{args.seed}.json").write_text(json.dumps(rows, indent=1))
    (args.out / f"absorption_manifest_s{args.seed}.json").write_text(json.dumps(
        {"script": "experiments/path_settings_absorption.py", **provenance(__file__),
         "argv": sys.argv[1:], "seconds": time.perf_counter() - t0,
         "minimap2": subprocess.run([str(args.minimap2), "--version"],
                                    capture_output=True, text=True).stdout.strip(),
         "parameters": {k: str(v) for k, v in vars(args).items()}}, indent=1))


if __name__ == "__main__":
    main()
