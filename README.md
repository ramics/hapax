<h1 align="center">hapax</h1>
<p align="center">Finding the viruses in wastewater that no database can name.</p>

<p align="center"><img src="docs/hapax-overview.svg" alt="Three-panel figure. Panel a: the seven stages of a run, with what each is given and what it writes — resolve, merge, assemble, identify, score the distance, report the residue, subtract other plants. Panel b: the three output classes, with what each reports and the qualifier that travels with it — divergent but related contigs, novel strains of known viruses, and the unplaceable residue. Panel c: the distance scale all three sit on, from far from published sequence to matched, and what every row carries — a distance, the closest relative with coverage and E-value, the operating point behind it, and its controls." width="100%"></p>

hapax subtracts the sequence a library already explains, assembles what survives, and reports
the rest as a distance: how far each contig sits from its closest published relative, which
relative that is, over how much of the contig, and at what significance. Identity comes per
gene, in nucleotide and protein, with the index it was measured against. Subtraction works on
read pairs, so every contig gets measured, named or not.

Three classes on one scale:

- **A** — contigs at real distance from anything published
- **B** — novel strains of named viruses, at the near end of that scale
- **C** — sequence nothing places, reported with the depth it was seen at

Recognition of a known virus holds to roughly 10% divergence and fades by 20%; detection needs
abundance above about one in a million. Verdicts are made per contig and per subject. A read is
too weak to carry one: at library scale a read recovers its own mate at chance
([`results/mate_friends/`](results/mate_friends/REPORT.md)).

## Setting up

```
python3 -m pip install -r requirements.txt
python3 -m pytest                              # needs no external tool
```

The binaries a run drives — BLAST, DIAMOND, Kraken2, foldseek — are pinned in
[`TOOLS.md`](TOOLS.md).

## Running it

```
python3 hapax.py run SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz --seed N --out RUNDIR \
    --library ACCESSION --trim-r2 K \
    --contigs IDENTIFIED_RRNA.fa IDENTIFIED_OTHER.fa --unidentified-contigs UNIDENTIFIED.fa \
    --identified-contigs IDENTIFIED_OTHER.fa --rrna-contigs IDENTIFIED_RRNA.fa \
    --assembly ASSEMBLY.fa --cross-sample-ref OTHER_SAMPLES.fa --spacers SPACERS.fa
python3 experiments/make_web.py results/web/pages/SAMPLE.json   # the run viewer
```

Per-sample values — trim, own contigs, assembly — are given on the command line and have no
defaults, because a default here is some other sample's file and a run would accept it in
silence. `--report-without cross_sample crispr ...` drops a field whose reference is missing
and records it as unsearched. Inputs, indexes and tools are checked before a run starts; a
missing one exits 2. `--plan` prints the stages and runs nothing.

`hapax.py` orchestrates only. Every number in `RUNDIR` comes from the `experiments/` script
that measured it, under that script's own seed and manifest, and the web view shows nothing
absent from a committed file.

## What a run writes

<p align="center"><img src="docs/hapax-outputs.svg" alt="Table of what a run writes, one row per unit. Library: one breakdown per run into identified, assembled but unnameable and unplaced, each fraction with its interval and operating point. Contig: one row per contig, whether or not an index named it, carrying distance to the closest published relative with identity and coverage per gene. Subject: one row per published relative the residue points at, with how many reads point at it and from how far. Read pair: one row per surviving pair, un-thresholded. Divergence ladder: one rung per divergence where a library carries a genome of known distance." width="100%"></p>

## Two earlier directions, archived with their results

- **A codon-aware assembler**, under [`experiments/codon_assembler/`](experiments/codon_assembler/).
  Negative: codon-space overlap scoring gained a little over nucleotide scoring at high
  simulated diversity, and lost to metaSPAdes and MEGAHIT on real multi-host swarms by merging
  sister swarms.
- **Strain separation during assembly**, which produced the benchmark in
  [`results/pair_ladder/`](results/pair_ladder/REPORT.md). Below about 3.8% divergence every
  assembler fuses lineages at any depth, while the reads retain enough information to separate
  them — so the limit lives in path selection. Live, and cited by the spec.

## Where things are

[`SPEC.md`](SPEC.md) is the spec, with the bounds on every claim in it.
[`results/README.md`](results/README.md) indexes the experiments, failures included.
[`CLAUDE.md`](CLAUDE.md) is how to work in the repo.
