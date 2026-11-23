<h1 align="center">hapax</h1>
<p align="center">Finding the viruses in wastewater that no database can name.</p>

<p align="center"><img src="docs/hapax-overview.svg" alt="Three-panel figure. Panel a: the seven stages of a run, with what each is given and what it writes — resolve, merge, assemble, identify, score the distance, report the residue, subtract other plants. Panel b: the three output classes, with what each reports and the qualifier that travels with it — divergent but related contigs, novel strains of known viruses, and the unplaceable residue. Panel c: the distance scale all three sit on, from far from published sequence to matched, and what every row carries — a distance, the closest relative with coverage and E-value, the operating point behind it, and its controls." width="100%"></p>

hapax removes the sequence in a library that is already known, assembles what survives, and reports the rest as a distance: for every contig, how far it sits from its closest published relative, which relative that is, over how much of the contig, and at what significance. Identity is given per gene, nucleotide and protein, with the index it was measured against. Removal applies to read pairs, so every contig is measured, named or not.

The output is three classes on one scale. **A**, contigs at real distance from anything published. **B**, novel strains of named viruses, at the near end of the same scale. **C**, sequence nothing places, with the depth it was seen at.

The signal the whole thing serves is **"new AND like X"**: a distance without a nearest-relative identity is half an answer.

## The rule the design is built around

**Searching only ever removes background. It is never used to score novelty.** An index can answer *is this already known*; it cannot answer *is this new*, because score ordering tracks database proximity. Measured, repeatedly. Any design that scores novelty with a search has reintroduced the failure this project is built around.

Two things follow, and they govern every number here. **Every number ships with its operating point**: thresholds dominate database size by an order of magnitude, and one Kraken2 database removes 42.8%, 21.7% or 3.5% of the same library depending only on its confidence setting. **Every field ships with its calibration**, a shuffled null beside a known-positive comparison; a field that reads the same on both is reported with that fact attached.

Which index you consult beats how fast you consult it, by orders of magnitude: 1.4 MB of a sample's own contigs removed more of the library than 272 MB of SILVA, and 15 MB of them removed 30 points where 16 GB of Kraken2 removed 0.4.

## How far it reaches

Recognition of a known virus holds to roughly 10% divergence and fades by 20%. Detection needs abundance above about one in a million. A verdict is made at contig and subject level. A read is the smallest unit the report works at and the weakest: the rate falls as the search space grows, so at library scale a read recovers its own mate at chance ([`results/mate_friends/`](results/mate_friends/REPORT.md)). These bounds are properties of the method, measured across libraries (`SPEC.md` §7, §9).

## What a run writes

<p align="center"><img src="docs/hapax-outputs.svg" alt="Table of what a run writes, one row per unit. Library: one breakdown per run into identified, assembled but unnameable and unplaced, each fraction with its interval and operating point. Contig: one row per contig, whether or not an index named it, carrying distance to the closest published relative with identity and coverage per gene. Subject: one row per published relative the residue points at, with how many reads point at it and from how far. Read pair: one row per surviving pair, un-thresholded. Divergence ladder: one rung per divergence where a library carries a genome of known distance." width="100%"></p>

## Running it

```
python3 -m pip install -r requirements.txt     # pinned where a number depends on it
python3 -m pytest                              # the suite; needs no external tool
```

The binaries a run drives (BLAST, DIAMOND, Kraken2, foldseek and the rest) are pinned in [`TOOLS.md`](TOOLS.md).

```
python3 hapax.py run SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz --seed N --out RUNDIR \
    --library ACCESSION --trim-r2 K \
    --contigs IDENTIFIED_RRNA.fa IDENTIFIED_OTHER.fa --unidentified-contigs UNIDENTIFIED.fa \
    --identified-contigs IDENTIFIED_OTHER.fa --rrna-contigs IDENTIFIED_RRNA.fa \
    --assembly ASSEMBLY.fa --cross-sample-ref OTHER_SAMPLES.fa --spacers SPACERS.fa
python3 experiments/make_web.py results/web/pages/SAMPLE.json   # the run viewer
```

Per-sample values are given on the command line and have no defaults, because a default here is another sample's file and the run would use it without complaint. Every input, index and tool is checked before the run starts; a missing one stops it with exit status 2. `--plan` prints the stages without running them.

`hapax.py` only orchestrates. Every number in `RUNDIR` comes from the `experiments/` script that measures it, under that script's own seed and manifest, and the web view shows no number absent from a committed file.

## History

Two earlier directions are archived here with their results, because what they ruled out determined the shape of the current design.

- **The codon-aware assembler**, under [`experiments/codon_assembler/`](experiments/codon_assembler/) with its own spec and results. The result is negative: scoring overlaps in codon space gave a small gain over nucleotide scoring, and only at high simulated diversity. On real multi-host swarms it lost to metaSPAdes and MEGAHIT by merging sister swarms.
- **Strain separation during assembly**, which produced the benchmark in [`results/pair_ladder/`](results/pair_ladder/REPORT.md). Below about 3.8% divergence every assembler fuses lineages at any depth. The reads carry enough information to separate them, which places the limitation in path selection. That result is live and the spec cites it.

[`SPEC.md`](SPEC.md) is the spec and the bounds on every claim in it, with §10 what is measured, §11 what is dead, and §12 what is open plus the experiment that would settle it. [`results/README.md`](results/README.md) indexes the experiments, failures included. [`CLAUDE.md`](CLAUDE.md) is how to work in the repo.
