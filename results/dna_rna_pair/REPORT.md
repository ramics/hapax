# DNA and RNA from the same water, on the same morning

`experiments/dna_rna_pair.py`, one 200,000-pair block per library (seed 1; seeds 2-5
queued — **nothing here is a headline result yet** under this repo's five-seed rule).

Every figure this project has published comes from one library: SRR23973193, Berlin,
2021-10-13, RNA. PRJNA948850 turns out to hold 24 `GENOMIC` libraries beside its 215 RNA
ones, and 22 of those have an RNA library from the same date and the same collection window.
So this compares:

- `DNA_220616_20-22hmix_cer_BIOMcombi` (SRR23973375) against
  `RNA_WW_220616_20-22hmix` (SRR23973265) — **same water, same morning, same prep, one
  variable: which nucleic acid was extracted.**
- Both against `RNA_WW_211013_20-22hmix` (SRR23973193), our calibration library — **same
  nucleic acid, same plant, same naming series, eight months apart.**

Two clean axes rather than one muddled one, and the second doubles as the control that says
whether the instrument here matches the instrument that produced our published numbers.

## The table

| | DNA 2022-06-16 | RNA 2022-06-16 | RNA 2021-10-13 (ours) |
|---|---|---|---|
| rRNA, ≥40 aligned bases against SILVA SSU+LSU | **0.01%** | **67.12%** | 38.01% |
| identified, Kraken2 PlusPF-16 at confidence 0.1 | **0.87%** | 40.65% | 22.65% |
| identified, RefSeq viral proteins, bitscore ≥50 | 7.37% | 6.49% | 11.96% |
| identified, Swiss-Prot, bitscore ≥50 | **0.98%** | 5.86% | 6.51% |
| mean GC | **0.385** | 0.508 | 0.499 |
| duplicate pairs (first 50 bases of both mates) | 20.0% | 23.6% | 22.1% |
| modal read length, R1/R2 | 107/107 | 107/107 | 107/107 |

Wilson intervals are in `dna_rna_pair_s1.json`; at 200,000 pairs they are within ±0.3
points of every figure above, so the differences below are not sampling.

**The instrument agrees with our published numbers.** Kraken2 PlusPF-16 at confidence 0.1
removes 21.7% of the 2021 library in `results/subtract_known/`; measured here, in the same
pass as everything else, it reads **22.65%**. The rRNA figure needs more care: our published
34.5% came from exact 31-mer matching and this is alignment of ≥40 bases, deliberately more
permissive, which reads **38.01%** on the same library — consistent, and the published
figure was explicitly a lower bound. Without that column in the same run, a drift and an
instrument change could not have been told apart, and the first draft of this comparison
nearly reported one as the other.

## 1. The DNA library is a virome, not a metagenome

Three independent facts, and they only fit together one way.

- **It has essentially no ribosomal sequence**: 0.01% against 38-67% for the RNA libraries.
  A bacterial DNA metagenome cannot look like that — rDNA is a few kilobases of every
  genome, so a few tenths of a percent of reads at minimum. Something removed the cells.
- **Its most abundant classified taxon is `Darmviridae`**, a family in the `Crassvirales`:
  crAss-like phage, the standard signature of the human faecal DNA virome.
- **Viral proteins identify 7.5 times what cellular proteins do** (7.37% against 0.98%). In
  the RNA library of the same water the same ratio is 1.1 (6.49% against 5.86%). So it is
  not that this library is hard to identify in general; what is identifiable in it is viral.

The `cer` and `centri` in the library names are presumably ceramic-filter and
centrifugation concentration, which is consistent, but that is reading a file name and is
not evidence. What the sequence says is: virome.

## 2. Identifiability inverts depending on which instrument you bring

| instrument | size | DNA | RNA |
|---|---|---|---|
| Kraken2 PlusPF-16, nucleotide | 16 GB | 0.87% | 40.65% |
| RefSeq viral proteins, DIAMOND | 0.24 GB | 7.37% | 6.49% |

**For the DNA library, 0.24 GB of viral protein identifies 8.5 times what 16 GB of
nucleotide index does.** For the RNA library the same 16 GB index wins by a factor of six.
This is the project's standing result — protein reaches further than nucleotide — in its
sharpest form yet, and it says the right tier order is a property of the *library*, not of
the tool. A cascade tuned on an RNA library and pointed at a DNA virome would spend its
largest index to remove under one percent.

**99.1% of the DNA library is unidentified** at the operating point the rest of this project
uses. Our RNA residue is 22.9%. So the residue problem is roughly four times larger on the
DNA side, and it is a different residue: unidentifiable phage rather than unidentifiable
anything.

## 3. The largest subtrahend is not a stable property of a library

Same plant, same group, same naming series, eight months apart, both RNA:

- rRNA **38.01% → 67.12%**
- Kraken2 PlusPF-16 at 0.1 **22.65% → 40.65%**
- viral protein at bitscore ≥50 **11.96% → 6.49%**

The ribosomal fraction nearly doubled, the identifiable fraction nearly doubled with it (the
extra rRNA is itself classifiable), and the virus-identifiable fraction nearly halved as
ribosomal sequence crowded it out. **A removal figure is a statement about one sample on one
day, not about wastewater.** This is the first hard evidence for the per-sample calibration
run `SPEC.md` requires, and it means every cascade percentage in this repo should be read as
"on SRR23973193" rather than as a property of the method.

It also warns about the positive control: a holdout spike's survival depends on how much of
the library is ribosomal, which moved by 29 points between two samples of the same plant.

## 4. GC is a nucleic-acid difference, not a seasonal one

0.385 for DNA against 0.508 and 0.499 for the two RNA libraries. The two RNA libraries agree
across eight months while the same morning's DNA differs by 12 points, which is the
cleanest available statement that these two libraries sample different biology rather than
the same biology twice. Useful as a sanity axis: a pipeline reporting similar GC for a
DNA/RNA pair from one sample has mixed something up.

## 5. Do the two nucleic acids see the same organisms?

Of 720 taxa called across the pair, **81 are shared** (Jaccard 0.11), and those 81 account
for 62.2% of the DNA library's classified reads and 87.1% of the RNA library's. So the
classified material concentrates on a shared core while the tails disagree almost entirely.

Read this one with caution: the DNA library has only 1,736 classified reads in 200,000, so
its 268 taxa are a thin sample and the Jaccard is as much a statement about depth as about
biology. The read-weighted columns are the defensible half.

## 6. Does either library explain the other's unidentified reads?

This is the question a matched pair is actually for. If the DNA library held the molecules
our RNA residue cannot account for, a second extraction from the same water would be a
cheaper subtrahend than any public database — which is what this project's thesis predicts
for sample-local data, since 1.4 MB of the sample's own rRNA contigs beat all 272 MB of
SILVA.

`experiments/cross_library_witness.py`. 20,000 pairs of the 2022 RNA library, split by
whether Kraken2 identified them, scored against two whole libraries with
`experiments/kmer_coverage/`. The readout is the fraction of a pair's exact canonical
31-mers occurring at least once in the subject — deliberately the weakest criterion, because
the question is presence, not depth.

| subject library | library pairs | identified | unidentified | gap | unidentified sharing nothing |
|---|---|---|---|---|---|
| DNA virome, same water, same morning | 11,311,548 | 16.5% | **13.3%** | −3.1 | **60.7%** |
| our 2021 RNA library, same plant, 8 months earlier | 8,731,224 | 96.6% | **71.2%** | −25.4 | 14.1% |

**The DNA library does not hold our residue.** It is the larger of the two subjects and
accounts for almost none of the RNA library either way, and the gap runs the wrong
direction: unidentified reads are witnessed *less* than identified ones, which is what rarer
material does. **A second extraction from the same water is not a cheap subtrahend**, and
the DNA side is its own residue problem rather than a key to ours. That is a negative, and
it is the answer to whether sequencing both is worth it for this purpose.

**Nucleic acid separates these libraries far more than eight months does.** 71.2% of the
*unidentified* 2022 reads have most of their 31-mers present in a different sample of the
same plant from eight months earlier, and only 14.1% share nothing with it. This is the
recurrence axis measured directly rather than inferred: the residue is not one-off sequence,
it recurs across samples of one catchment at exact 31-mer identity. It is also a warning —
reads we describe as absent from everything are mostly not absent from the same plant in a
different season.

**Bias, stated rather than corrected:** the query reads are RNA, so the RNA subject shares
ribosomal sequence with them that a virome cannot. That inflates the second row, and 71.2%
is therefore an upper bound on the recurrence of non-ribosomal material. Correcting it needs
a third library this project does not have.

**Two bugs, both caught by the control.** The first version read the coverage table by
hard-coded column index, one off, and produced a fraction present of −0.71; an impossible
value is the cheapest bug report available and a plausible one would have shipped. The
header is now parsed into a name map, asserted against the tool's own header line, and every
fraction checked to lie in [0, 1]. Second, the tool's own library blocks cannot serve as the
control: it subtracts a library pair's own occurrences and does not do so for the
query-only path, so a library read is scored more harshly than a foreign one and the control
moves against the test. Both subjects here are scored through the identical query-only path.

## What is not measured here

- **One seed.** Five are queued. No figure above should be quoted until they land.
- ~~**The non-template prefix is unresolved.**~~ **RESOLVED, and the detector was wrong** —
  see [`../nontemplate_prefix/REPORT.md`](../nontemplate_prefix/REPORT.md). Mapping untrimmed
  reads to SILVA gives a 5' soft-clip mode of **14** in 5/5 seeds on both RNA libraries, with
  95% of anchored R2 reads clipping ≥10 and the same reads flush at their 3' end, against
  mate 1's mode of 0. The composition detector read **9 and 16**, and it was setting the trim
  in the run that produced the table above. `--r2-trim` now carries the measured 14, the
  detector is a diagnostic with a guard, and every figure here has been re-run on it.

  **But I first called 16 an over-trim, and that is withdrawn.** Re-running at 14 made the
  Kraken2 fraction go *down* — 0.4065 → 0.3899 on the 2022 library and 0.2265 → 0.2159 on the
  2021 one — from two more bases of what the mapping calls template. The prefix length is a
  distribution (mode 14 at 51% of anchored reads, 13 at 16%, 15 at 6.5%, 16 at 2.2%, with a
  tail), and Kraken2's confidence is a *fraction* of a read's minimizers, so two extra
  minimizers where that variability lives lower every read's confidence and push marginal
  ones under the threshold. 14 is the measured prefix; 16 is the better operating point for
  this classifier. The table above is at trim 14, and
  `../nontemplate_prefix/trim_sensitivity_*.json` is the curve that says what the choice
  costs — 1.7 points, the size of effects this repo reports as findings.

  Two caveats stay. The **DNA library's trim rests on the detector alone**, because SILVA
  anchored 0 to 11 of its reads per 100,000 and no reference here can replace it for a DNA
  virome. And the prefix is measured on the ribosomal fraction only, since SILVA is the
  anchor.
- **No cross-library subtraction yet.** The question this pair is really for — can the DNA
  library identify reads the RNA library's residue cannot explain, and vice versa — needs
  the DNA library's own assembly or k-mer index as a subtrahend. That is the next run, and
  it is the one that would say whether sequencing both is worth it.
- **Whether the 2022 RNA library has the same residue structure** as the 2021 one. Every
  read-level field in this repo was calibrated on the 2021 library; with 67% of this one
  ribosomal, the residue is a different mixture.
