<h1 align="center">hapax</h1>
<p align="center">Finding the viruses in wastewater that no database can name.</p>

<p align="center"><img src="docs/hapax-overview.svg" alt="Three-panel figure. Panel a: the seven stages of a run, with what each is given and what it writes — resolve, merge, assemble, identify, score the distance, report the residue, subtract other plants. Panel b: the three output classes, with what each reports and the qualifier that travels with it — divergent but related contigs, novel strains of known viruses, and the unplaceable residue. Panel c: the distance scale all three sit on, from far from published sequence to matched, and what every row carries — a distance, the closest relative with coverage and E-value, the operating point behind it, and its controls." width="100%"></p>

hapax removes the sequence in a library that is already known, assembles what survives, and
reports the rest as a distance. For every contig: how far it sits from its closest published
relative, which relative that is, over how much of the contig, and at what significance.
Identity is given per gene, nucleotide and protein, with the index it was measured against.
Removal applies to read pairs, so every contig is measured, named or not.

The output is three classes on one scale, ordered by distance from published sequence. A,
novel strains of known viruses, at the near end. B, contigs at real distance from anything
published, ranked on the gene that is conserved across a family and carrying the spread
between a contig's own genes, which is where the system's novelty work sits. C, sequence
nothing places, with the depth it was seen at.

[`SPEC.md`](SPEC.md) is the spec, including the bounds on every claim in it.
[`results/README.md`](results/README.md) indexes the experiments behind it, failures included.
[`CLAUDE.md`](CLAUDE.md) is how to work in the repo.

## Setting up

```
python3 -m pip install -r requirements.txt     # Python side, pinned where a number depends on it
python3 -m pytest                              # the suite; needs no external tool
```

The binaries a run drives (BLAST, DIAMOND, Kraken2, foldseek and the rest) are pinned in
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

One command builds the web view for any sample. What a page states about its sample that is
not in the run directory sits in that sample's descriptor under `results/web/pages/`: the
title, one sentence on the water, and every result file outside the run that a panel reads,
as a path or as `not_run` with a reason. Methods entries live in `results/web/methods.json`.
A row measured on one library is stored under that accession. On another sample's page it
reads "Not measured on this library".

Per-sample values are given on the command line: trim, own contigs, assembly. None of them
have defaults. A default here is another sample's file, and the run would use it without
complaint. A per-read field whose reference is unavailable is dropped by name
(`--report-without cross_sample crispr ...`) and recorded as not searched. Every input, index
and tool is checked before the run starts. A missing one stops it with exit status 2.

`hapax.py` only orchestrates. Every number in `RUNDIR` comes from the `experiments/` script
that measures it, under that script's own seed and manifest. The web view reads `RUNDIR` and
shows no number absent from a committed file. `--plan` prints the stages without running
them.

## What a run writes

<p align="center"><img src="docs/hapax-outputs.svg" alt="Table of what a run writes, one row per unit. Library: one breakdown per run into identified, assembled but unnameable and unplaced, each fraction with its interval and operating point. Contig: one row per contig, whether or not an index named it, carrying distance to the closest published relative with identity and coverage per gene. Subject: one row per published relative the residue points at, with how many reads point at it and from how far. Read pair: one row per surviving pair, un-thresholded. Divergence ladder: one rung per divergence where a library carries a genome of known distance." width="100%"></p>

## How far it reaches

Recognition of a known virus holds to roughly 10% divergence and fades by 20%. Detection
needs abundance above about one in a million. A read is the smallest unit the report works
at and the weakest: nearest-neighbour geometry over reads carries molecule-level
information, and the rate falls as the search space grows, so at library scale a read
recovers its own mate at chance ([`results/mate_friends/`](results/mate_friends/REPORT.md)).
A verdict is made at contig and subject level. These bounds are properties of the method,
measured across libraries (`SPEC.md` §7, §9).

## History

Two earlier directions are archived here with their results. What they ruled out determined
the shape of the current design.

- **The codon-aware assembler** was the first attempt, under
  [`experiments/codon_assembler/`](experiments/codon_assembler/) with its own spec and
  results. It is laid out as the repo root once was, so its scripts still run from there.
  The result is negative. Scoring overlaps in codon space gave a small gain over nucleotide
  scoring, and only at high simulated diversity. On real multi-host swarms it lost to
  metaSPAdes and MEGAHIT by merging sister swarms.
- **Strain separation during assembly** produced the benchmark in
  [`results/pair_ladder/`](results/pair_ladder/REPORT.md). Below about 3.8% divergence every
  assembler fuses lineages at any depth. The reads carry enough information to separate
  them, which places the limitation in path selection. That result is live and the spec
  cites it.
