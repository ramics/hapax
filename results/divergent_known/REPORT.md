# The distance a known thing was removed at

`experiments/divergent_known.py`, seed 1, one pass over the identified pile of
SRR23973193. Result file: [`divergent_known_s1.json`](divergent_known_s1.json). Run under
`hapax run` as the `divergent` stage, behind `--identified-contigs`.

This is the only experiment here that runs over what the cascade **removed** rather than
what it left. The cascade's subtrahend is one bit. A contig matching its reference at 99.9%
and one matching it at 35% record the same word, are both subtracted in tier 1, and neither
is reported again, because `strains` only looks at the residue. The identity is computed and
then discarded.

## Operating point

11,361 of 25,209 identified non-ribosomal contigs are at least 300 nt, which is the length
at which a 60 aa alignment can be most of the contig. DIAMOND blastx `--sensitive` against
RefSeq viral proteins and Swiss-Prot. A hit counts at **E < 1e-5 and at least 60 aa
aligned**. 10,172 contigs clear it. Gene level on the 600 longest divergent contigs and all
301 near-reference ones: Prodigal in metagenome mode, minimum 50 aa, DIAMOND blastp, Pfam at
the curators' gathering cutoffs. All 17,462,448 mates mapped back, mate 2 trimmed 16 nt.

## 1. The word "identified" spans eighty points of identity

| band | contigs | share | shuffled null |
|---|---|---|---|
| below 40%, no name survives | 3,355 | 32.98% | 106 |
| 40 to 70%, family level or looser | 5,820 | 57.22% | 176 |
| 70 to 90% | 563 | 5.53% | 15 |
| 90 to 97%, divergent strain | 387 | 3.80% | 7 |
| 97% and above, reference | 47 | 0.46% | 0 |

**434 of 10,172 identified contigs, 4.27%, sit at 90% or above.** 90.2% sit below 70%, where
the subject lends its name rather than establishing one. The cascade subtracts all five
bands with the same word.

This is not an error rate. Most of this pile is bacterial sequence whose own species is
absent from Swiss-Prot, and 40% identity to a characterised enzyme in a related organism is
the normal and correct reading of such a contig. The finding is that one bit cannot carry
the difference between that case and a 99.9% match, and that the pipeline currently makes
no record of which it was.

## 2. The null: existence separates, and nothing conditional on a hit can be tested

Each contig was dinucleotide-shuffled, preserving length and dinucleotide composition, and
put through the identical search in the same pass. The shuffle is
`experiments/shuffle_control.py`, imported.

**Of 11,361 shuffled contigs, 0 carry a significant alignment.** 89.53% of the real
contigs do. The whole shuffled library produces one alignment of any kind, 47 aa at
E=8.45e-4, which clears neither the E-value nor the length floor.

An earlier version of this section reported 2.68% for the null and a 33-fold ratio. That
number was an artefact: the shuffle it used returned its input unchanged when its walk
dead-ended, which happened for 324 contigs, and all 304 "shuffled" hits came from those
324. The correction is recorded in `CLAUDE.md` and in the commit that replaced the shuffle.

Two consequences:

- **Existence of a significant alignment is the calibrated readout, and the separation is
  total at this operating point.** E < 1e-5 over at least 60 aa against RefSeq viral
  proteins plus Swiss-Prot admits nothing from composition-matched noise. There is no
  ratio to quote, because the denominator is zero.
- **Every field conditional on a hit now has no null at all.** Identity, alignment length,
  subject coverage, bitscore and contig coverage were previously compared against a
  shuffled distribution; that distribution was real sequence. With a correct shuffle the
  comparison cannot be made, because the control pile is empty. Those fields are neither
  calibrated nor killed here. The real distributions are reported for what they are:
  identity p10 32.0, median 45.1, p90 69.6; aligned length 88 / 143 / 420 aa; contig
  coverage 0.407 / 0.815 / 0.987; subject coverage 0.151 / 0.342 / 0.905.

**Scope.** This calibrates whether a hit is real. It says nothing about whether identity
measures the distance between a real contig and its real relative. The holdout answers that
one directly and the answer is that identity works: a genome 30.28% diverged read 66.5%
protein identity over 620 aa where a referenced comparator from the same library read 99.9%,
ordered correctly and 30 points apart. The usable statement is the narrow one — a
significant alignment's identity does not tell you the query is real, and on short
alignments it is quantised and optimistic, so it ships with its coverage and its aligned
length every time.

Coverage and length floors used to be reported here as the lever that moves the null. With
the null at zero they move only the real pile, and that is a cost curve rather than a
calibration:

| floor | real contigs | shuffled |
|---|---|---|
| none | 89.53% | 0.00% |
| subject coverage 0.5, 100 aa | 26.32% | 0.00% |
| subject coverage 0.8, 150 aa | 9.92% | 0.00% |
| subject coverage 0.9, 200 aa | 6.45% | 0.00% |

A null of zero over 11,361 shuffles bounds the chance rate at roughly 1 in 11,000 at best;
it is not evidence that the rate is zero. The untested candidate remains low-complexity
sequence, which a dinucleotide shuffle preserves and which aligns to anything, and the
control for it is the aligner's own compositional masking, not run here.

## 3. One identity per contig is whichever gene won on bitscore

Over 600 divergent and 302 near-reference contigs, median per-contig spread between a
contig's own genes:

- divergent arm **12.8 points**, near-reference arm **2.45 points**.

`results/divergent_contigs/` measured 17.4 against 2.35 on the 89-taxon viral panel. This
pass reproduces the effect on a pile an order of magnitude larger and of different
composition.

**The arms are not matched as drawn, and the confound runs one way.** Divergent contigs carry
a median of 3 scorable genes over 3,820 nt; near-reference contigs carry 1 over 552 nt,
because a contig matching its reference at 95% is usually a short fragment of a well-sampled
organism. More genes and more length give a spread more room, which biases the gap upward.
Matched on the number of scorable genes, the gap holds at every band:

| scorable genes | divergent | control |
|---|---|---|
| 2 | 5.6 (n=167) | 1.3 (n=28) |
| 3 | 16.9 (n=347) | 2.65 (n=24) |
| 4 | 22.2 (n=27) | 3.4 (n=7) |

So a whole-contig identity is a property of the best gene, not of the contig, and the error
that introduces grows with distance. Rank on the conserved gene and show the spread.

## 4. What this changes in the pipeline

The flag worth having is already computed and thrown away. "Removed, and only barely
recognisable" is a distance with a coverage beside it, which is a different object from a
novelty score: it says how far the thing we removed was from the thing that removed it, and
says nothing about whether anything else in the world is closer. Ordering subjects by
identity would order them by database proximity, which is the failure this project is built
around, so this stage reports a distribution and per-contig rows and draws no ranking of
novelty.

`subtract_known.py` now records the same quantity at read level: the identity distribution
of what the protein tier removed, over alignments of at least 30 aa, with the shorter ones
counted separately because identity is quantised there.

## 4b. The coverage gap: the kill is withdrawn, and the background rate is the number

The holdout thread proposed a field and asked for the null it lacked. Two coverages of the
same contig from the same significant HSPs: **per-subject** is the union against whichever
single subject covers most of the contig, **any-subject** the union against all of them, and
the **gap** between them is the part that known proteins explain only collectively. On the
holdout library the withheld genome reads +36.6 points, the 95.85th percentile of 193
length-matched identified contigs whose gap is median 0.1, and 34 contigs from
already-referenced genomes read median, mean and maximum 0.0.

This section previously killed the field on a shuffled null that read worse than the real
pile. That null was not shuffled; with the shuffle corrected it contains no hits at all, so
**there is no shuffled gap distribution and the kill is withdrawn.** The field is untested
against noise, not cleared by it.

What the run does give is the field's rate on real identified background, which is the
comparison that matters for a flag meant to fire on a residue:

| gap over 10,172 identified contigs | |
|---|---|
| median | 0.38 points |
| p90 | 21.5 |
| p99 | 46.2 |
| above 10 points | 14.51% |
| above 30 points | 5.81% |

So 30 points is the 94.2nd percentile of already-identified contigs and the holdout's +36.6
sits between that and the 99th. A threshold there fires on a few hundred identified contigs
per library. That is consistent with the holdout thread's own length-matched percentile and
is the operating point the field has to be read at, rather than against the 0.0 that
referenced genomes read.

One direction of the field is solid without any null. A gap near zero means one subject
explains the contig, which is why 34 referenced contigs reading a maximum of 0.0 is a real
statement. A large gap is necessary for multi-gene architecture and not sufficient for it.
The version worth testing next is the gap conditioned on every contributing HSP exceeding a
length floor; the holdout candidate's two come from 620 aa and 493 aa, and a shuffled null
cannot reach those lengths at all.

## 5. The band is a continuum, and it is not short

Measured by the holdout thread on the same question with a different instrument: all 35,328
contigs against a RefSeq viral index with the holdout genome withheld, coverage taken as the
union of HSPs over the contig.

- **No gap.** Identity against coverage is a continuum and the holdout candidate sits inside
  the mass, not off its edge.
- Of 22,295 identified contigs, 7,635 clear E<1e-5 at protein level, median identity 40.0%
  (p5 27.1, p95 71.1). 75.5% of those sit below 50%.
- **The candidate's own cell is crowded.** At coverage ≥0.70 and identity 60 to 70%, where
  the holdout sits at 66.5% over 0.829 coverage, there are **304 contigs, 4.02% of the
  pile**, and its neighbours are more striking than it is: 4,156 nt at 68.0% over 95%
  coverage; 2,072 nt at 62.2% over 100%. A threshold drawn there fires about 300 times per
  library and the holdout is not top of its own list.
- **Identity fails in both directions.** The candidate is at the **93rd percentile of
  identity** among identified contigs: more similar to the database than 93% of what the
  cascade removed. No low-identity rule finds it.
- Coverage times identity does separate known from distant on known material: the
  already-referenced leak arm sits at coverage ≥0.7 and identity ≥90 in 31 of 74 significant
  contigs. It separates known from distant, not the holdout from its peers.

So the distance record is a measurement to report, not a filter to act on. The cascade
deletes 6,781 contigs whose nearest viral protein relative is under 60% identical, and
nothing in protein identity against a narrow index shortens that list. Whether a broader
index shortens it is now partly measured: against core_nt, 16.8% of the divergent contig
panel collapses to a strain-level identification and 55.7% stays distant
(`results/divergent_contigs/REPORT.md` §7). The protein equivalent against nr is still
running.

## 5b. This result has actually been reproduced

Two provenance changes landed on this script on the same day - the shuffle moved into
`experiments/shuffle_control.py`, and another thread made every named index required rather
than skipped when absent. They were merged before re-running rather than re-run once each,
so the result was recomputed from a materially different code state.

**Every number above came back identical**, to the digit: 11,361 contigs searched, 10,172
significant, 9,738 below 90%, 0 of 11,361 shuffles, the same identity histogram, the same
12.8 against 2.45 gene spread and the same gene-count-matched bands. Only the manifest
differs, which is what it is for.

That is a stronger claim than "reproducible from a committed script plus a seed", and almost
nothing in this repo can make it. It was not earned by intending to: the provenance rule in
`CLAUDE.md` forces a re-run whenever provenance changes, and the reproduction falls out of
that for free.

## 6. Blind spots

- Both reference sets are narrow, so each identity is an upper bound on identity and each
  divergence a lower bound on divergence.
- This stage measures the distance for one tier, the protein tier. The cascade's other tiers
  remove by Kraken2 and by nucleotide match and are not scored here.
- Protein identity understates nucleotide divergence: synonymous changes are invisible to
  it.
- Nothing here says whether a divergent contig is a new species. Demarcation criteria are
  per family and nothing in this stage fetches them.
- A contig is already a consensus over reads and hides within-sample strain mixtures below
  about 3.8% divergence.
- The gene-level arm is capped by contig length, so a short divergent contig can carry a
  whole-contig distance and no gene-level row.
- One library, one date. Read every fraction as specific to SRR23973193.

## Cost

951 s wall on 4 cores: DIAMOND blastx on the real contigs and on the null, blastp at gene
level, Pfam, and mapping 17.5 M mates. The null roughly doubles the search cost. It bought
one number here - zero - and it cost a day, because the first version of it was wrong in a
way that read as a passed control.
