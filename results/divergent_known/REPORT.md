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
aligned**. 10,170 contigs clear it. Gene level on the 600 longest divergent contigs and all
301 near-reference ones: Prodigal in metagenome mode, minimum 50 aa, DIAMOND blastp, Pfam at
the curators' gathering cutoffs. All 17,462,448 mates mapped back, mate 2 trimmed 16 nt.

## 1. The word "identified" spans eighty points of identity

| band | contigs | share | shuffled null |
|---|---|---|---|
| below 40%, no name survives | 3,350 | 32.94% | 106 |
| 40 to 70%, family level or looser | 5,821 | 57.24% | 175 |
| 70 to 90% | 564 | 5.55% | 16 |
| 90 to 97%, divergent strain | 388 | 3.82% | 7 |
| 97% and above, reference | 47 | 0.46% | 0 |

**435 of 10,170 identified contigs, 4.3%, sit at 90% or above.** 90.2% sit below 70%, where
the subject lends its name rather than establishing one. The cascade subtracts all five
bands with the same word.

This is not an error rate. Most of this pile is bacterial sequence whose own species is
absent from Swiss-Prot, and 40% identity to a characterised enzyme in a related organism is
the normal and correct reading of such a contig. The finding is that one bit cannot carry
the difference between that case and a 99.9% match, and that the pipeline currently makes
no record of which it was.

## 2. The null: existence separates, nothing about the hit does

Each contig was dinucleotide-shuffled, preserving length and dinucleotide composition, and
put through the identical search in the same pass.

- **Existence**: 89.52% of real contigs carry a significant alignment, against 2.68% of
  shuffles. A factor of 33.
- **Conditional on a hit existing, nothing separates them.** Median identity 45.1% real
  against 43.5% shuffled. Scored as the share of real hits above the shuffled p90, where a
  p90 puts 10% above it by construction: identity 11.98%, alignment length 7.10%, subject
  coverage 7.52%, bitscore 7.63%. The one field carrying anything is the fraction of the
  contig the alignment covers, 26.34%, with real median 0.814 against shuffled 0.496.

**Scope.** This calibrates whether a hit is real. It says nothing about whether identity
measures the distance between a real contig and its real relative. The holdout answers that
one directly and the answer is that identity works: a genome 30.28% diverged read 66.5%
protein identity over 620 aa where a referenced comparator from the same library read 99.9%,
ordered correctly and 30 points apart. The usable statement is the narrow one — a
significant alignment's identity does not tell you the query is real, and on short
alignments it is quantised and optimistic, so it ships with its coverage and its aligned
length every time.

Coverage and length floors are the only lever that moves the null, and they cost real hits:

| floor | real | shuffled |
|---|---|---|
| none | 89.52% | 2.68% |
| subject coverage 0.5, 100 aa | 26.53% | 1.14% |
| subject coverage 0.8, 150 aa | 9.96% | 0.49% |
| subject coverage 0.9, 200 aa | 6.47% | 0.23% |

2.68% is high for a null. The untested candidate is low-complexity sequence, which a
dinucleotide shuffle preserves and which aligns to anything; the control is the aligner's
own compositional masking, not run here. Until that is measured, 2.68% is an upper bound on
chance and 33 a lower bound on the ratio.

## 3. One identity per contig is whichever gene won on bitscore

Over 600 divergent and 301 near-reference contigs, median per-contig spread between a
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
index shortens it is what the running nr and core_nt pass decides.

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

911 s wall on 4 cores: 221 s of DIAMOND blastx on the real contigs, 213 s on the null, 162 s
of blastp, 108 s of Pfam, the rest mapping 17.5 M mates. The null doubles the search cost and
is the reason the result above can be stated at all.
