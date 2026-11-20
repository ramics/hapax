# Real wastewater, second probe: a divergent astrovirus with the authors' own assembly as comparator

Run with `experiments/ww_berlin.sh`. Manifest beside this file. Companion to
`results/ww_cryptic/REPORT.md`, which recorded a dataset that could not discriminate.

## The dataset

Wyler et al. 2024, Environment International 190:108875. 17 months of total nucleic acid from
one Berlin treatment plant, published as GEO **GSE228220** / SRA **PRJNA948850**, 239 runs.
Unusually, the authors deposited their metaSPAdes and Trinity assemblies alongside the reads,
so a third assembly of the same sample is already on the table.

Probe sample: `SRR23973193` (RNA_WW_211013_20-22hmix), 8.7M read pairs, 2x107 bp, 825 MB.
Chosen as the astrovirus-richest sample whose published contigs are actually readable (see
the gotcha below). Target: the one near-complete astrovirus genome in the authors' SPAdes
output for this sample, 6,823 bp, **57% identical to HAstV-1** (3,900/6,801) — a real
divergent-virus case, not a reference-matching one.

## There is within-sample diversity here, and it is mostly synonymous

Counted from the reads with mpileup, independently of any assembler: of 5,000 positions on the
target genome at >=20x, **82 carry a minor allele between 20% and 80%** (1.6%). This is the
mixture that `SRR15202285` lacked.

Of the 536 such sites hapax calls across its astrovirus contigs, those falling in a single
reading frame split **254 synonymous to 60 non-synonymous — 81% synonymous**. That is
METHOD.md §1's premise ("most of the variation causing the tangle is synonymous") measured on
a real wastewater viral population rather than on the simulator.

## Assemblies against the authors' 6,823 bp contig

| | contigs touching it | covered | identity | misjoins |
|---|---|---|---|---|
| authors' SPAdes (whole sample) | 1 | 6,802 / 6,823 | 100% | 0 |
| MEGAHIT (whole sample) | 12 | 6,709 / 6,823 | 99.2% | 0 |
| hapax (5,171 recruited reads) | 4 | 6,663 / 6,823 | 97.0% | 0 |

hapax reaches 97.7% of the genome in 4 contigs, longest 4,147 bp, with no misjoin. On breadth
it is level with MEGAHIT and touches the target with far fewer pieces. **On identity it is
behind: 97.0% against MEGAHIT's 99.2%.** The gap is larger than the 1.6% polymorphism rate, so
it is not explained by hapax having picked a different haplotype. It needs chasing before any
accuracy claim is made — consensus error at low depth is the first suspect, since the recruited
set gives only ~43x.

**Dispatch:** 27.6% of candidate pairs fell through to the codon HMM, against 12.7% on the
SARS-CoV-2 wastewater data. Higher, but still a minority, so the two-tier design in
CLAUDE.md's order of work item 2 is holding on real input.

## The finding that matters most: `find_orfs` over-calls dual-coding

hapax labels **26% of its astrovirus contig bases (2,769 of 10,837) as dual-coding**, and 40%
on the narrower single-target run. The true layout of an astrovirus is essentially
single-coding: GenBank NC_001943 gives ORF1a/1b as `join(86..2845,2845..4380)` — two frames in
sequence via a ribosomal frameshift, not overlapping — and ORF2 at `4373..6733`, an 8 bp
overlap. Genuine dual-coding here is ~8 bases, not ~2,800.

Cause, at `src/consensus.rs:307`: `find_orfs` accepts every stop-free stretch of >=300 nt in
all six frames, with no start-codon requirement and no scoring, then `annotate_regions` marks
any base covered by two such stretches as dual. Chance stop-free runs on the reverse strand and
in alternate frames are common in real coding sequence, so they pile up.

This matters because METHOD.md §3 switches dual regions to nucleotide scoring, on the argument
that dual-coding regions have little synonymous variation to collapse. Applied to a spurious
call, that argument inverts: hapax turns off its own codon advantage over a quarter to a
third of a genome where the advantage should apply. A simulator emitting one clean ORF per
entity cannot surface this; only real coding sequence does.

## Caveats

- MEGAHIT and SPAdes assembled all 8.7M read pairs; hapax was given 5,171 reads recruited by
  mapping to astrovirus bait. Recruitment is not what limits hapax here — widening the bait
  from the SPAdes contig to that plus MEGAHIT's own astrovirus contigs plus four references
  moved the recruited set only from 5,171 to 5,186 reads, because the virus has ~5.2k reads in
  this run. But hapax's background contigs are therefore uncounted, and *finding* the virus in
  sewage is not being tested at all.
- One sample, one run, no replication. The identity figures come from minimap2 `asm10`
  alignment blocks and will move a little under other settings.
- The 6,823 bp target is itself an assembly, not a finished genome, so "100%" for SPAdes is
  circular by construction. It is a comparator, not truth.

## Dataset gotchas worth keeping

- `GSE228220_Spades_RNA_all.fa.gz` is ~1 GB and the GEO copy ends mid-stream: `gzip -t` fails,
  `zcat` yields 101 of 116 samples and errors at the end. Samples after WW_220616 are
  unreadable from it, which is why this probe uses 211013 and not the astrovirus-richest
  sample, 220629.
- The `S<N>` suffixes in the kaiju summary and in the SPAdes contig headers are *different* run
  indices. `S70` is WW_220629_22-24hmix in the kaiju table and WW_210527_22-24h in the SPAdes
  headers. Match samples by name.
- Trinity contigs are deposited as one combined file with no sample prefix, so they cannot be
  split per sample; only SPAdes supports a per-sample comparison.
