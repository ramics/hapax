"""A mock metagenome for the novelty question: can an assembler separate a
virus from its nearest relative when neither is in any reference?

The haplotype panels test resolution within one species. This builds the other
half: real bacterial and viral genomes mixed at realistic abundances, where
several viruses are present *as pairs* -- a "known" member and a relative at a
measured divergence. The question each pair asks is whether an assembler
returns two genomes or one, because collapsing a divergent relative into its
known neighbour is exactly how a novel virus gets missed, and joining them is
how a novel virus gets invented.

The novelty ladder, divergence measured by edit distance over the whole genome
(printed by `--report`):

| pair | divergence |
|---|---|
| SARS-CoV-2 / bat RaTG13 | 0.040 |
| HBV genotype D / genotype A | 0.105 |
| West Nile lineage 1 / lineage 2 | 0.205 |
| SARS-CoV-2 / bat CoV BM48-31 | 0.240 |
| phage lambda / phage HK97 | 0.458 |
| HPV16 / HPV11 | 0.506 |
| HHV-5 / HHV-7 | 0.546 |

Design choices:

* **Background dominates the reads.** Three bacterial genomes (13 Mb between
  them) at low coverage, as in a real stool or wastewater library where the
  virus is a few percent of the data. Without that, the viral contigs are
  unrealistically easy and the k-mer thresholds every assembler uses never
  come under pressure.
* **Abundance spans orders of magnitude.** A novel virus in real data is
  usually rare, so the second member of each pair sits below the first. The
  ladder of coverages is explicit in `COMMUNITY` rather than drawn at random,
  so that a failure can be read off against a coverage.
* **Phages are in there on purpose.** crAss-like and lambdoid phages carry
  the mosaic, recombined structure that produces false joins; a phage-free
  mock community flatters every assembler.
* **The within-host swarm stays in.** The ZM247F env panel goes in as its own
  component, so the sample contains both between-species novelty and
  within-species swarm at once, which is the real case.
* **Reads come from ART (Illumina) and Badread (nanopore)**, not from this
  repo's simulators: at metagenome scale the per-base Python loop is too slow,
  and these two are what the assemblers' own papers use, so the error model is
  not a thing this benchmark gets to choose.

Accessions and download dates are in data/MANIFEST.md. Nothing here designs
or modifies a sequence: it mixes published genomes at chosen abundances.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

# name, accession, coverage, role, nearest relative in the mix (or None)
COMMUNITY = [
    # background: bacterial chromosomes, low coverage, most of the read mass
    ("ecoli", "NC_000913.3", 3.0, "background", None),
    ("saureus", "NC_007795.1", 3.0, "background", None),
    ("pvulgatus", "NC_009614.1", 2.0, "background", None),
    # phages: mosaic genomes, moderate abundance
    ("phage_lambda", "NC_001416.1", 60.0, "known", "phage_hk97"),
    ("phage_hk97", "NC_002167.1", 20.0, "novel", "phage_lambda"),
    ("phage_t4", "NC_000866.4", 30.0, "known", None),
    ("phix174", "NC_001422.1", 150.0, "known", None),
    ("crassphage", "NC_024711.1", 40.0, "known", None),
    # known viruses with a divergent relative also present
    ("sars2", "NC_045512.2", 100.0, "known", "ratg13"),
    ("ratg13", "MN996532.2", 25.0, "novel", "sars2"),
    ("batcov_bm4831", "NC_014470.1", 10.0, "novel", "sars2"),
    ("hbv_d", "NC_003977.2", 120.0, "known", "hbv_a"),
    ("hbv_a", "X70185.1", 30.0, "novel", "hbv_d"),
    # West Nile lineage 1/2 stands where HIV-1 subtype B/C stood. Full HIV-1
    # genomes cannot be in the same truth set as the env swarm: a swarm
    # haplotype IS the env gene of an HIV-1 genome, so a contig running from
    # one into the other is a nesting in the truth set rather than a
    # misassembly, and the chimaera count becomes unreadable. The pair is
    # replaced by one at a comparable divergence with nothing nested in it.
    ("wnv1", "NC_009942.1", 80.0, "known", "wnv2"),
    ("wnv2", "NC_001563.2", 20.0, "novel", "wnv1"),
    ("hpv16", "NC_001526.4", 70.0, "known", "hpv11"),
    ("hpv11", "NC_001525.1", 15.0, "novel", "hpv16"),
    ("hhv5", "NC_006273.2", 15.0, "known", "hhv7"),
    ("hhv7", "NC_001716.2", 5.0, "novel", "hhv5"),
    # singletons with no relative in the mix: the easy cases, as controls
    ("ttv1", "NC_002076.2", 60.0, "known", None),
    ("norovirus", "NC_039477.1", 50.0, "known", None),
    ("rabies", "NC_001542.1", 40.0, "known", None),
    ("veev", "NC_075022.1", 30.0, "known", None),
]


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def one_seq(path):
    d = read_fasta(path)
    return max(d.values(), key=len)


def build(refdir: Path, out: Path, swarm_panel: Path | None, swarm_cov: float, art: str, badread: str,
          read_len: int, nanopore_reads: int, seed: int, threads: int = 4):
    out.mkdir(parents=True, exist_ok=True)
    comp_dir = out / "components"
    comp_dir.mkdir(exist_ok=True)
    truth = out / "truth_genomes.fa"
    rows = []
    with open(truth, "w") as tf:
        for name, acc, cov, role, rel in COMMUNITY:
            s = one_seq(refdir / f"{acc}.fa")
            (comp_dir / f"{name}.fa").write_text(f">{name}\n{s}\n")
            tf.write(f">{name} depth={cov} host={name} accession={acc} role={role} coverage={cov} relative={rel}\n{s}\n")
            rows.append(dict(name=name, accession=acc, coverage=cov, role=role, relative=rel, length=len(s)))
        if swarm_panel:
            for header, s in read_fasta(swarm_panel / "truth_haplotypes.fa").items():
                hap = header.split()[0]
                (comp_dir / f"swarm_{hap}.fa").write_text(f">swarm_{hap}\n{s}\n")
                tf.write(f">swarm_{hap} depth={swarm_cov} host=env_swarm accession={hap} role=swarm coverage={swarm_cov} relative=env_swarm\n{s}\n")
                rows.append(dict(name=f"swarm_{hap}", accession=hap, coverage=swarm_cov, role="swarm",
                                 relative="env_swarm", length=len(s)))

    # Illumina: ART per component, then concatenate
    r1, r2 = out / "illumina_1.fq", out / "illumina_2.fq"
    for p in (r1, r2):
        p.unlink(missing_ok=True)
    for i, row in enumerate(rows):
        pre = comp_dir / f"{row['name']}_art"
        cmd = [art, "-ss", "MSv3", "-i", str(comp_dir / f"{row['name']}.fa"), "-p", "-l", str(read_len),
               "-f", str(row["coverage"]), "-m", "400", "-s", "50", "-rs", str(seed + i), "-na",
               "-o", str(pre) + "_"]
        subprocess.run(cmd, check=True, capture_output=True)
        for src, dst in ((f"{pre}_1.fq", r1), (f"{pre}_2.fq", r2)):
            with open(dst, "ab") as f, open(src, "rb") as g:
                f.write(g.read())
            Path(src).unlink()

    # nanopore: Badread over the whole truth set, abundance from coverage x length
    total = sum(r["coverage"] * r["length"] for r in rows)
    nano = out / "nanopore.fq"
    with open(nano, "wb") as f:
        p = subprocess.run([badread, "simulate", "--reference", str(truth),
                            "--quantity", f"{nanopore_reads}x", "--seed", str(seed),
                            "--error_model", "nanopore2023", "--qscore_model", "nanopore2023",
                            "--length", "6000,3000"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True)
        f.write(p.stdout)

    meta = dict(components=rows, read_len=read_len, seed=seed, nanopore_quantity=f"{nanopore_reads}x",
                total_cov_bp=total, illumina_bp=sum(r["coverage"] * r["length"] for r in rows))
    (out / "community.json").write_text(json.dumps(meta, indent=2))
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refdir", type=Path, default=Path("data/community"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--swarm-panel", type=Path, default=None, help="a bench_panel output dir to fold in")
    ap.add_argument("--swarm-coverage", type=float, default=30.0)
    ap.add_argument("--art", default="art_illumina")
    ap.add_argument("--badread", default="badread")
    ap.add_argument("--read-len", type=int, default=250)
    ap.add_argument("--nanopore-reads", type=int, default=1, help="Badread quantity, as a coverage multiple")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    m = build(a.refdir, a.out, a.swarm_panel, a.swarm_coverage, a.art, a.badread, a.read_len,
              a.nanopore_reads, a.seed)
    print(f"{len(m['components'])} components, {m['illumina_bp']/1e6:.0f} Mbp Illumina")


if __name__ == "__main__":
    main()
