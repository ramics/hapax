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

**Why the controls sit above §1's curve, and why it does not matter here.** The
controls run at 3.77×, where §1's curve interpolates to about 35% breadth, and
they give 40.6%. Two candidate explanations, separated by scoring the same
control contigs both ways (`--cross-score`, `two_lineage_controls_rescored.jsonl`):
§5 aligns to both lineages at `asm20` and 85% identity, §1 to one reference at
`asm10` and 90%, which is necessarily more permissive. Measured over 8 runs, the
looser scoring accounts for 1.5 ± 0.8 points — real but small. The rest is
between-genome variation: these are eight different deposited rabies genomes and
their individual control breadths run from 31.8% to 48.3%, a spread that
contains §1's interpolated 35% comfortably, while §1's curve is one genome
(NC_001542) at two seeds. So the two measurements agree within the variation
either can resolve. It matters only that the controls and the mixtures in this
section are scored identically to each other, which they are; the comparison
with §1 is a sanity check, not a term in the result.

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
absent at 7.3%, and partially back at 9.7%.

Read as noise, that is a wobble on a declining curve. Read as mechanism, it is
two processes with different ranges, and the 9.7% rung separates them. Pooling
that works genome-wide needs exact k-mer sharing along the genome, which falls
away by about 4%: this is the process that produces the length gain at 2.0% and
is gone by 7.3%. Joins through individual conserved windows need sharing only in
one place, so they happen at any divergence: this is what still makes three
mosaic contigs at 9.7% while contributing no breadth, and it is the mechanism
this project already measured directly when two HBV genomes from different hosts
shared 709 exact 150-mers. Only the second survives into the divergence range
where real surveillance samples sit, and only the first looks like an
improvement on a contiguity table.

The practical consequence is that a single threshold cannot be quoted. Pooling
is strong at 2.0%, already erratic at 3.8% where one seed of three crossed
1,500 nt and two did not, and absent at 7.3%, so "safe above 4%" is a claim
about a transition this ladder only brackets. The three added rungs in §5.1 are
there to say where inside the bracket it actually goes.

**Why this matters for what gets built.** Below about 4% an assembler silently
trades lineage identity for contig length, and nothing in its output says which
trade it made. That is a gain when a classifier needs length and a disaster when
the lineages are the answer. The same event cannot be prevented without losing
the length, which argues for labelling it: a per-contig flag saying "this contig
is supported by reads that disagree with each other in a structured way", which
is the one thing no assembler emits.

### 5.1 The 4–7% band, and why there is no threshold

The first ladder left 4–7% empty. That was an artefact of having downloaded 40 of
the ~2,840 deposited complete rabies genomes, not a property of the data: the full
pairwise matrix over 200 of them (`find_divergence_pairs.py`, 19,900 pairs) gives
real pairs at 4.50%, 5.50% and 6.52%. The band is thin — 570 of 19,900 pairs,
real lineage structure — but thin is not empty, so the rungs are real sequence
rather than a genome mutated to order. These are also the first rungs anywhere in
this project between 2.5% and 10%: `results/stress/ladder.md` has eight rungs
with nothing between 2.505% and 3.953%.

| Divergence | Mode | Breadth A | Breadth B | Longest | Mass ≥1500 | Mosaic |
|---|---|---|---|---|---|---|
| 4.50% | A alone | 38.2 ± 12.0% | — | 831 ± 326 | 0 | 0.0 |
| 4.50% | B alone | — | 40.5 ± 2.4% | 1,103 ± 175 | 0 | 1.0 |
| 4.50% | both | 39.3 ± 8.7% | 45.7 ± 4.2% | 1,622 ± 557 | 1,355 ± 2,346 | 2.7 |
| 5.50% | A alone | 38.8 ± 2.2% | — | 684 ± 21 | 0 | 0.7 |
| 5.50% | B alone | — | 47.3 ± 2.4% | 773 ± 176 | 0 | 0.0 |
| 5.50% | both | 44.5 ± 2.1% | 36.8 ± 4.9% | 1,636 ± 234 | 1,169 ± 1,019 | 3.3 |
| 6.52% | A alone | 40.0 ± 10.0% | — | 946 ± 132 | 0 | 0.7 |
| 6.52% | B alone | — | 38.5 ± 2.1% | 750 ± 164 | 0 | 0.0 |
| 6.52% | both | 35.1 ± 7.0% | 42.1 ± 4.4% | 1,055 ± 206 | 0 | 3.0 |

**It is a rate, not a threshold.** Across all 63 runs, the fraction of mixtures
in which any contig crosses 1,500 nt reads 3/3 at 2.0%, 2/3 at 3.8%, 1/3 at
4.5%, 2/3 at 5.5%, 0/3 at 6.5%, 0/3 at 7.3% and 2/3 at 9.7%. A single lineage at
half the budget crosses it in **0 of 42** control runs. So pooling is what
produces classifier-useful length, it does so in roughly half of runs over the
whole 3.8–9.7% range, and the rung-to-rung pattern is what sampling a rate of
that size with three seeds looks like rather than a transition to locate. The
earlier reading of §5 — a window that closes by about 4% — does not survive the
band being filled in.

The mechanism is the one the k-mer measurement implies (`kmer_survival.json`).
Shared 21-mer fraction sits within a few points of the uniform expectation
(1−d)^21 for every pair, so at 4.5–9.7% a 12 kb pair still shares 1,600–4,600
exact 21-mers: bridges are abundant, and what varies between seeds is whether one
falls where local coverage can use it. That is a property of the draw, which is
exactly what produces a rate. It also retires the two-mechanism paragraph above:
the longest exactly identical stretch is 85–111 nt for every pair at 4.5% and
above, against 107 nt reads, and even at 2.0% the stretches at actual join points
are 81–201 nt. Conserved windows long enough to join reads outright are real —
two HBV genomes from different hosts share 709 exact 150-mers — but they are not
what is operating here. One mechanism, isolated 21-mer bridges, at every rung.

**The binary readout says no dependence on divergence. The continuous one says
there is.** Pooling the seven rungs gives 10 of 21 mixtures crossing 1,500 nt, a
rate of 0.48 with a Wilson 95% interval of 0.28 to 0.68. An exact permutation
test on the trend statistic Σ kᵢdᵢ, enumerating all 352,716 arrangements of 10
crossings over 21 runs conditional on the margins, gives a two-sided P of 0.159.
On that evidence the rate looks constant across 2.0–9.7%, and the apparent dip at
6.5–7.3% is not evidence of anything.

But "did any contig cross 1,500 nt" throws away most of each run. Regressing log
longest-contig-length on divergence instead (`length_trend.py`) finds a clear
effect: over 46 mixture runs the slope is −0.0745 per percentage point of
divergence, 95% interval −0.119 to −0.030, permutation P = 0.0018. The longest
contig therefore shrinks by about 7% per point, roughly halving between 2% and
10%. Mass over 1,500 nt declines the same way (−0.74 per point, P = 0.011). The
42 single-lineage controls are the null anchor — divergence between two lineages
cannot affect a run containing one of them — and their slope is flat, +0.013
(−0.019 to +0.045), P = 0.42, so this is not a confound from different genomes
sitting at different rungs.

Two honest limits. The significant fit leans on the unbalanced rate study, which
contributes 18 runs at 4.5% and 13 at 7.3%; the balanced seven-rung ladder alone
gives the same direction and magnitude but does not reach significance (−0.056,
−0.117 to +0.004, P = 0.066). And the relationship is not clean: mean longest
contig at 9.7% (1,652 nt) exceeds 6.5% (1,055 nt), so a monotone slope is a
summary rather than a description.

So the earlier framing in this section was wrong in both directions. Collapse is
not a window that closes at 4%, and it is not a rate independent of divergence
either. Length declines steadily with divergence while the probability of
crossing any particular length threshold stays too noisy to order the rungs — a
threshold on divergence is still unquotable, but because the outcome is
continuous and the variance is large, not because divergence does not matter.

**A limit on the lineage caller, not on the mixtures.** Control runs contain a
single lineage, so a mosaic there is impossible by construction, yet 13 of 42
controls report one and the control mean is 0.33 (max 2). That is a measured
false-positive rate for the per-window lineage call itself — presumably where a
300 nt window genuinely matches both genomes or sits at a contig end — and it
means the mosaic *label* is noisy, not merely the count. Mixture runs average
3.19 with a minimum of 1, so the signal clears the floor, but two consequences
follow. A mosaic count here reads as "about three above a baseline of a third"
rather than as an exact count. And a collapse flag must be trained on the 1,500 nt
crossing event or on a label derived from alignment to the known truth, never on
this caller's mosaic count, with the control runs serving as the caller's measured
error rate. The crossing rate itself is unaffected: no control ever crossed.

**This is a dependency on the deliverable, not a footnote.** A per-contig collapse
flag cannot be calibrated more accurately than the labels it is trained against,
so a caller that invents a mosaic in 31% of runs where one is structurally
impossible is a ceiling on the flag, not a detail of this section. Replacing it —
base-resolution lineage assignment from the persisted alignments, scored against
single-lineage controls as the error rate — is unbudgeted work that comes before
any flag, and is listed here so it is not discovered later.
