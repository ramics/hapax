"""Where does a relative stop being resolvable? A ladder in divergence.

In the mock community MEGAHIT separated every divergent relative from its
known neighbour except the closest pair, SARS-CoV-2 and RaTG13 at 4%
divergence, which came back only half covered. That is one point on a curve,
so this script walks the curve: the same two-genome mixture at a series of
divergences, everything else held fixed.

Rather than mutating a genome to hit a target divergence -- which would make
the answer a property of my mutation model -- the ladder is built from real
sequence pairs whose divergence is what it is, measured and reported:

* real pairs from the community set (SARS-CoV-2/RaTG13 0.040, HBV D/A 0.105,
  HIV-1 B/C 0.163, SARS-CoV-2/BM48-31 0.240, lambda/HK97 0.458),
* and, to fill the gap below 4% where no pair in that set sits, pairs of real
  within-host env haplotypes binned by their own pairwise divergence, which
  covers 0.0004 to 0.19 continuously.

Each rung is assembled on its own, as a two-member mixture at a fixed
coverage ratio, so the only thing changing along the ladder is divergence.
The number reported per rung is whether both members come back as separate
genomes, and at what coverage of each.
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
from pathlib import Path

import edlib


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def divergence(a, b):
    return edlib.align(a, b, task="distance")["editDistance"] / max(len(a), len(b))


REAL_PAIRS = [
    ("sars2_ratg13", "NC_045512.2", "MN996532.2"),
    ("hbv_d_a", "NC_003977.2", "X70185.1"),
    ("hiv1_b_c", "NC_001802.1", "AY772699.1"),
    ("sars2_bm4831", "NC_045512.2", "NC_014470.1"),
    ("lambda_hk97", "NC_001416.1", "NC_002167.1"),
]


def env_pairs(env_fa: Path, targets=(0.002, 0.005, 0.01, 0.02), tol=0.4):
    # targets are settable from the command line so the ladder can be filled in
    # where it matters: the interesting region is between 2.5% and 10%.
    """Real env haplotype pairs nearest each target divergence."""
    s = read_fasta(env_fa)
    s = {k: v for k, v in s.items() if set(v) <= set("ACGT")}
    ks = sorted(s)
    pairs = []
    for a, b in itertools.combinations(ks[:40], 2):
        d = divergence(s[a], s[b])
        if d > 0:
            pairs.append((d, a, b))
    pairs.sort()
    out = []
    for t in targets:
        best = min(pairs, key=lambda p: abs(p[0] - t))
        if abs(best[0] - t) <= t * tol:
            out.append((f"env_{t}", best[1], best[2], best[0], s[best[1]], s[best[2]]))
    return out


def run_rung(name, a_seq, b_seq, div, out: Path, art, megahit, cov_a, cov_b, read_len, seed, threads=4):
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "a.fa").write_text(f">member_a\n{a_seq}\n")
    (d / "b.fa").write_text(f">member_b\n{b_seq}\n")
    (d / "truth.fa").write_text(
        f">member_a depth={cov_a} host=member_a role=known coverage={cov_a} relative=member_b\n{a_seq}\n"
        f">member_b depth={cov_b} host=member_b role=novel coverage={cov_b} relative=member_a\n{b_seq}\n")
    r1, r2 = d / "r1.fq", d / "r2.fq"
    for p in (r1, r2):
        p.unlink(missing_ok=True)
    for i, (f, cov) in enumerate((("a.fa", cov_a), ("b.fa", cov_b))):
        pre = d / f"art_{i}_"
        subprocess.run([art, "-ss", "MSv3", "-i", str(d / f), "-p", "-l", str(read_len), "-f", str(cov),
                        "-m", "400", "-s", "50", "-rs", str(seed + i), "-na", "-o", str(pre)],
                       check=True, capture_output=True)
        for src, dst in ((f"{pre}1.fq", r1), (f"{pre}2.fq", r2)):
            with open(dst, "ab") as g, open(src, "rb") as h:
                g.write(h.read())
            Path(src).unlink()
    asm = d / "megahit"
    subprocess.run([megahit, "-1", str(r1), "-2", str(r2), "-o", str(asm),
                    "--presets", "meta-sensitive", "-t", str(threads)], check=True, capture_output=True)
    return dict(name=name, divergence=round(div, 5), contigs=asm / "final.contigs.fa", truth=d / "truth.fa")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refdir", type=Path, default=Path("data/community"))
    ap.add_argument("--env", type=Path, default=Path("data/env/ZM247F.fa"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--art", required=True)
    ap.add_argument("--megahit", required=True)
    ap.add_argument("--cov-a", type=float, default=100.0)
    ap.add_argument("--cov-b", type=float, default=25.0)
    ap.add_argument("--read-len", type=int, default=250)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--env-targets", type=float, nargs="+", default=[0.002, 0.005, 0.01, 0.02])
    ap.add_argument("--skip-real-pairs", action="store_true")
    a = ap.parse_args()

    rungs = []
    for name, acc_a, acc_b in ([] if a.skip_real_pairs else REAL_PAIRS):
        sa = max(read_fasta(a.refdir / f"{acc_a}.fa").values(), key=len)
        sb = max(read_fasta(a.refdir / f"{acc_b}.fa").values(), key=len)
        rungs.append((name, sa, sb, divergence(sa, sb)))
    for name, _, _, div, sa, sb in env_pairs(a.env, tuple(a.env_targets)):
        rungs.append((name, sa, sb, div))
    rungs.sort(key=lambda r: r[3])

    results = []
    for name, sa, sb, div in rungs:
        r = run_rung(name, sa, sb, div, a.out, a.art, a.megahit, a.cov_a, a.cov_b, a.read_len, a.seed)
        results.append({k: (str(v) if isinstance(v, Path) else v) for k, v in r.items()})
        print(f"{name}\tdiv={div:.5f}\t{r['contigs']}")
    (a.out / "ladder.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
