# Where assembly stops and read-level evidence has to start

Reproduce with `experiments/abundance_sweep.py`, `experiments/classifier_arm.py` and
`experiments/readlevel_sweep.py`; parameters, git SHA and timings in `manifest.json`,
per-run rows in `rows.jsonl`, geNomad output in `classifier.json`, the read-level
arm in `readlevel_s1.json`.

## Why this run exists

Two directions for this project disagreed about the same unmeasured number. One
says viral classifiers need long contigs, so the assembler sets the sensitivity
ceiling and the work is downstream lift per assembler. The other says that at
surveillance abundance there is no contig at all, so only read-level evidence
exists. A Lander-Waterman calculation made before the run said the two are about
an order of magnitude apart in coverage, so one sweep can answer both.

The 1500 bp figure that framed the first argument does not survive checking. Wu
et al. 2024 (*Genome Biology* 25:97) used it because it is "the minimum length
of the input sequence that can be handled by all tools", dropping shorter
contigs with seqtk before any tool saw them: the compatibility floor of the
weakest tool in a panel of nine, not a property of any one of them. geNomad
states no minimum; ViraLM asks 500 bp and was trained on 300 bp-2 kb fragments.
So recall is reported here as a curve against contig length.

## Setup

- **Background**: 1,000,000 read pairs of real Berlin wastewater, SRR23973193
  (Wyler et al. 2024), mate 2 trimmed from base 15 to drop this library's ~14 bp
  non-template prefix. Held fixed across every rung and arm, so a difference
  along the curve cannot come from the background.
- **Spike**: rabies virus NC_001542 (11,932 nt), ART MSv3 paired reads, rungs as
  absolute read-pair counts that are nested subsamples of one pool per seed.
  Coverage, not relative abundance, is the x-axis: mean contig length depends on
  coverage and minimum overlap, and relative abundance here is an artefact of the
  background's size.
- **Held out** means Rhabdoviridae removed from the search database by NCBI
  taxonomy, verified: none of the spike's five proteins and no header naming
  Lyssavirus survive. Taxonomy rather than name matching, because name matching
  would drop a herpesvirus protein called "vesicular stomatitis virus matrix
  protein" and keep strain-level rhabdovirus names.
- **Nearest retained relative**: 35.3% identity over 2,158 aa (L, to an
  unclassified rhabdo-like metagenomic virus) and 27.4% (N). The other three
  proteins have no homolog at all, so about a third of the genome is invisible to
  translated search at any depth. This is a database property, and it is reported
  because it would otherwise look like a method property.
- 48 assemblies, mean ± sd over 2-3 seeds per rung, MEGAHIT 1.2.9 `--min-contig-len 200`.

## 1. The assembly curve

| Spike coverage | Breadth | Longest contig | Mass in contigs ≥1500 nt | Contigs |
|---|---|---|---|---|
| 0.17× | 0.0% | 0 | 0 | 0 |
| 0.50× | 0.0% | 0 | 0 | 0 |
| 1.68× | 4.2 ± 1.4% | 509 ± 165 | 0 | 1 |
| 2.51× | 15.8 ± 1.1% | 532 ± 48 | 0 | 4 ± 1 |
| 3.35× | 24.8 ± 5.2% | 684 ± 79 | 0 | 6 ± 1 |
| 5.03× | 64.5 ± 6.5% | 1,304 ± 555 | 848 ± 1,199 | 15 ± 4 |
| 7.54× | 91.1 ± 1.7% | 2,818 ± 856 | 5,689 ± 1,031 | 9 ± 3 |
| 10.06× | 97.6 ± 0.8% | 6,365 ± 2,462 | 10,359 ± 752 | 4 |
| 16.76× | 99.2 ± 0.9% | 11,846 ± 91 | 11,846 ± 91 | 1 |
| 50.28× | 99.5 ± 0.1% | 11,910 ± 12 | 11,910 ± 12 | 1 |

**Zero chimaera candidates in all 48 runs.**

The prediction made before the run was that 500 nt needs about 2.6× and 1,500 nt
about 4.8×. Observed: 500 nt is crossed between 1.7× and 2.5×, and 1,500 nt
between 5.03× and 7.54×, where the 5.03× rung straddles it with a standard
deviation larger than its mean. Real over-dispersed coverage is slightly worse
than the idealised calculation at the long end, which is the direction the
SecureBio spike-in (about 10× by Poisson, stopped at a zero-coverage gap) would
predict.

## 2. The minimum-overlap arm, which inverted

In a de Bruijn assembler the k-mer size *is* the minimum overlap two reads must
share, so it is the lever. It only works in one direction.

| Arm | 1.68× | 5.03× | 16.8× | 50.3× |
|---|---|---|---|---|
| default (k 21-141) | 509 nt | 1,304 nt | complete | complete |
| low k (21, 29, 39) | 444 nt | 1,296 nt | complete | complete |
| high k (71, 99, 127) | nothing | nothing | nothing | 2.6-4.8 kb in fragments |

Demanding long overlaps loses the genome outright up to 17×. Relaxing the
requirement gains nothing, because MEGAHIT's default already sweeps k from 21.
So on this assembler the setting is a way to lose a genome rather than a lever on
contig mass, and a reported gain from relaxing minimum overlap is probably a
statement about an assembler that does not sweep k.

## 3. geNomad, classifier fixed and assembler varied

One geNomad 1.9.9 pass at its default `--min-score 0.7` over all 62,365 kept
contigs. No published experiment holds a modern classifier fixed and varies the
assembler on the same reads: Sutton et al. 2019 varied 16 assemblers but scored
by reference alignment, with no classifier, and predates geNomad, VirSorter2 and
ViraLM.

Spike contigs, end to end:

| Coverage | Contigs | Called viral | Mass called of mass assembled |
|---|---|---|---|
| 1.68× | 2-4 | 0 | 0 of 1.0-1.4 kb |
| 5.03× | 30-33 | 8-9 (27%) | 3.6 of 16.2 kb |
| 16.8× | 1 | 1 | all |
| 50.3× | 1 | 1 | all |
| 50.3×, high k | 8 | 1 | 1.0 of 7.4 kb |

So assembly and classification fail in the same place rather than in different
ones: sensitivity is zero below about 2×, partial at 5×, complete by 17×. The
classifier floor is not the binding constraint — the contigs that reach it are
already long enough to be called.

**The background rate is the finding.** geNomad calls viral:

| Background contig length | Called viral |
|---|---|
| 300-500 nt | 16.9% (5,600 of 33,096) |
| 500-1,000 nt | 31.6% (6,697 of 21,173) |
| 1,000-1,500 nt | 48.2% (1,973 of 4,092) |
| 1,500-3,000 nt | 63.9% (1,467 of 2,297) |
| ≥3,000 nt | **98.3%** (1,594 of 1,622) |

This is a wastewater virome, so a high viral fraction is expected and most of
those calls are probably correct. That is exactly why it matters: in this sample
a viral call carries almost no information at the long end. A per-contig
false-positive budget cannot be built on "is it viral" where nearly everything
long is viral. The discriminating question is **new versus already known**, and
no tool surveyed answers that with a calibrated score.

## 4. The read-level arm: a clean negative result

DIAMOND blastx on single reads against the family-removed database — the arm a
real surveillance pipeline runs. 2,000,000 background reads, mate pairs scored
separately.

| Bitscore | Spike sensitivity | Background flagged per million |
|---|---|---|
| ≥40 | 8-10% | 93,479 |
| ≥50 | 2-3% | 58,934 |
| ≥80 | 0% | 1,489 |
| ≥150 | 0% | 0 |

**There is no operating point.** Where the spike is detectable the background is
four to five orders of magnitude above a one-per-million budget; where the
background goes quiet the spike has disappeared. Only 2,576 of 22,200 spike
reads (11.6%) get any hit at all, and `--very-sensitive` changes that to 2,613
(11.8%), so this is not a settings artefact.

The mechanism: a 107 nt read is about 35 aa, and at 27-35% identity to the
nearest retained relative those alignments score below ~50 bits, while the
background's own well-represented viruses hit their references hard. The score
ordering tracks database proximity, not novelty — the same lesson as the geNomad
background rate, reached from the other side.

## What this settles

1. **Below ~1× there is nothing to classify, and translated search cannot
   substitute.** Both halves of that are now measured rather than argued.
2. **Between 2× and 10× the assembler is the whole game**, and the classifier
   floor is not what binds; by 17× everything works.
3. **"Viral or not" is the wrong label for this problem.** Any claim about
   false positives per million in a virome background has to be about novelty,
   not viralness, or it is measuring a property of the sample.
4. **Minimum overlap is a way to lose genomes, not a lever**, at least on an
   assembler that already sweeps k.

## Caveats that would change the picture

- One spike genome, one background sample, 2-3 seeds per rung where the repo rule
  is five for headline numbers. The curve's shape is robust; the exact crossing
  points are not.
- One assembler. metaSPAdes and PenguiN were designed in but not run: PenguiN
  matters most, since its authors state it falls back to consensus below ~10×,
  which is exactly this regime.
- The read-level arm tests translated search only. Nothing here used composition,
  codon statistics or a protein language model, and nothing grouped reads, so the
  proposal those would support is untested — this sets its floor, not its ceiling.
- The spike is divergent (nearest retained relative 35% aa identity on L). A
  target with a 70%-identity relative in the database would make the read-level
  arm look far better and is the obvious next rung.

## 5. Two related lineages at fixed total coverage

Proposed by the sister thread as the rung that unifies this curve with the
stress test: if an assembler will join reads across a small difference, then two
lineages at 3.75× each are one genome at 7.5× to it, and 7.5× is where length
arrives. The stress test's "collapse" and this sweep's "length" would then be the
same event scored oppositely.

Four pairs of complete rabies genomes, divergence measured with edlib rather than
simulated, 3 seeds, same fixed background, same total budget of 450 read pairs
(≈7.5× of one genome) either split between two lineages or spent on one.
Lineage is called per 300 nt window, so a contig that switches is counted as the
mosaic it is. `two_lineage_rows.jsonl`.

| Divergence | Mode | Breadth A | Breadth B | Longest | Mass ≥1500 | Mosaic contigs |
|---|---|---|---|---|---|---|
| 2.0% | A alone, 3.8× | 41.9 ± 4.4% | — | 725 ± 189 | 0 | 0 |
| 2.0% | B alone, 3.8× | — | 41.9 ± 4.5% | 771 ± 152 | 0 | 0 |
| 2.0% | both, 3.8× each | 51.0 ± 11.5% | 41.1 ± 11.1% | 1,971 ± 239 | 3,605 ± 1,377 | 5 ± 2 |
| 3.8% | both | 52.9 ± 16.7% | 39.9 ± 11.6% | 2,143 ± 1,046 | 2,323 ± 2,708 | 3 |
| 7.3% | both | 35.3 ± 5.2% | 45.2 ± 10.4% | 1,053 ± 126 | 0 | 2 ± 1 |
| 9.7% | both | 43.2 ± 12.0% | 33.1 ± 3.6% | 1,652 ± 289 | 1,180 ± 1,042 | 3 ± 1 |

Single-lineage controls at half budget behave identically at every divergence
(41–47% breadth, longest 725–963, and **zero** mass over 1,500 nt in all 12
runs), so the comparison is clean.

**Pooling is real and it converts coverage into length.** At 2.0% the mixture's
longest contig is 1,971 nt against 725–771 for either lineage alone, and only the
mixture ever crosses 1,500 nt. Five of about twelve contigs carry windows from
both lineages, so the length is bought with mosaicism — an assembly that scores
better on every contiguity metric while being a chimaera of two real things.

**The strong form of the prediction is falsified.** Tolerant assembly was
predicted to return ~91% breadth of a mosaic, matching the one-lineage 7.54×
rung. It does not: per-lineage breadth in the mixture is 33–53%, barely above the
half-budget controls, and far from 91%. MEGAHIT pools enough to lengthen contigs
and fabricate mosaics, not enough to reconstruct a genome. So the two phenomena
are related but not identical: collapse is partial, and it costs lineage identity
without buying back the coverage it appeared to merge.

**The benefit is not monotone in divergence.** It is strongest at 2.0%, erratic
at 3.8% (one seed of three crossed 1,500 nt, hence the ± larger than the mean),
absent at 7.3%, and partially back at 9.7%. The 9.7% rung still produces three
mosaic contigs, which points at joins through conserved windows rather than
genome-wide pooling — the mechanism this project already measured when two HBV
hosts shared 709 exact 150-mers.

**Why this matters for what gets built.** Below about 4% an assembler silently
trades lineage identity for contig length, and nothing in its output says which
trade it made. That is a gain when a classifier needs length and a disaster when
the lineages are the answer. The same event cannot be prevented without losing
the length, which argues for labelling it: a per-contig flag saying "this contig
is supported by reads that disagree with each other in a structured way", which
is the one thing no assembler emits.
