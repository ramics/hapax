# What the residue is made of

Berlin wastewater SRR23973193, five 100,000-pair blocks, mate 2 trimmed 14 nt. The residue
is the 22.90% of pairs the cascade (`results/subtract_recut/`) could not identify at its
operating points. This thread's question, Imo's words on 2026-11-21: *what on earth the
residue actually is*, rather than whether a planted holdout can be found in it.

Three measurements, each with its control in the same pass. **The residue is real
sequence, it is not a few hidden abundant genomes, and it is thin — mostly molecules seen
once.** None of that is good news for detection, and it is a cleaner problem than artefact
would have been.

## 1. It is not artefact

`experiments/residue_artefacts.py`, `artefacts_s{1..5}.json`. The cascade's artefact tier
only looked for the DNA-seq kind — adapters, primer dimers, low entropy, homopolymers — and
removed 0.31% of the library, and `results/absence_residual/` found adapter content flat
(0.0004 against 0.0001) with flat entropy. Neither looked for what an RNA-seq library
actually makes. Mean of 5 seeds, per pair, with min–max in the JSONs:

| field | residue | identified, rRNA | identified, other | residue shuffled |
|---|---|---|---|---|
| fold-back hairpin ≥18 nt, mate 1 | 0.18% | 0.06% | 0.02% | 0.00% |
| fold-back hairpin ≥18 nt, mate 2 | 0.17% | 0.07% | 0.05% | 0.00% |
| fold-back hairpin ≥24 nt, mate 1 | 0.10% | 0.01% | 0.01% | 0.00% |
| mates on the same strand | 0.55% | 0.14% | 0.12% | 0.00% |
| insert ≤107 nt (from mate overlap) | 7.10% | 5.80% | 3.32% | — |
| mate 1 shorter than 100 nt | 2.68% | 1.63% | 1.18% | 2.68% |
| poly(A/T) run ≥12, mate 1 | 0.02% | 0.00% | 0.01% | 0.02% |
| UniVec ≥30 bases, mate 2 | 0.00% | 0.19% | 0.13% | 0.00% |
| PhiX ≥30 bases | 0.00% | 0.00% | 0.00% | 0.00% |

**The template-switching signatures are enriched and still negligible.** Fold-backs run 3–9×
higher in the residue than in identified reads and 5–10× higher at the stricter ≥24 nt
threshold, and same-strand mates — which one linear fragment cannot produce — run about 4×
higher. Every one of them is under 1% of the residue. Shuffled reads score zero on all of
them, so the fields are measuring structure rather than composition.

**Polarity, measured on known material as the conventions require.** UniVec and poly(A/T)
run *higher on identified reads* than on the residue, so neither points at the residue; the
UniVec hits are vector-like sequence in real genomes, not contamination in the unknown pile.
PhiX is absent from this library at any pile. Short inserts are about twice as common in the
residue (7.1% against 3.3%), as is a short mate 1 — more degraded or short fragments, which
is a property of the material, not junk, and it is the one artefact-adjacent field with a
real difference.

**Verdict.** The standing worry — *if almost everything real hits something, whatever matches
nothing is mostly artefact* — does not hold for this library. Under 1% of the residue carries
any artefact signature tested. The chimaera question is unaffected and still sits where
`results/partial_alignment/` left it: after assembly, not on reads.

## 2. It does not hide abundant unknown genomes

`experiments/library_labels.py`, `library_label_classes.py`, `residue_assembly.py`.

Every one of the library's 8,731,224 pairs was labelled against SILVA 138.2 SSU+LSU and all
82,947 of this sample's own SPAdes contigs, per mate, as matching-base counts
(`data/work/library_labels.npz`, shared). Classes at the cascade's 40-base operating point,
using the pipeline thread's three-way contig split:

| class, per mate | mate 1 | mate 2 |
|---|---|---|
| rrna (SILVA ≥40 or an rRNA contig) | 38.76% | 38.06% |
| identified non-rRNA contig | 33.32% | 32.74% |
| unidentifiable contig | 8.67% | 8.26% |
| a contig outside the split (7,944 contigs) | 1.01% | 1.28% |
| nothing at ≥40 | 18.23% | 19.65% |

*Validation, before anything was read off it:* on the five seed blocks the table's rRNA
fraction is 37.96–38.13% against the cascade's rrna tier at 37.95–38.13%, and 0 of ~22.9k
residue pairs per seed reaches 40 bases on SILVA. It is the same instrument.
*Stated bias:* the table holds no read-level Kraken2 or DIAMOND, so "nothing" over-counts
against the cascade — it errs toward calling material unplaced.

**1,372,021 pairs (15.71% of the library) have neither mate matching anything at ≥40 bases.**
Assembled on their own with MEGAHIT (meta-sensitive, ≥200 nt), away from the rRNA and
abundant transcripts that crowd a whole-library graph, they give:

| | |
|---|---|
| contigs | 3,695 |
| total bases | 1,306,995 |
| median length | 331 nt |
| ≥500 nt / ≥1,000 nt | 288 / 7 |
| longest | 1,385 nt |

1.37 M pairs is ~0.29 Gbase of sequence returning 1.31 Mb of contigs. Mapping 200,000 of
those pairs back recruits 3.8% of them, and the ten best contigs take 4.2% of what recruits,
so there is no dominant piece. **Whatever the unplaced fraction is, it is many molecules each
too shallow to assemble, not a few abundant unknowns.**

## 3. It is thin: the unplaced pile sits at depth 1

`experiments/residue_depth.py` over `experiments/kmer_coverage` — exact canonical 31-mer
counts for every query read over the whole library, each pair's own occurrences subtracted,
so a count is *occurrences elsewhere in the run*. Five blocks, mean (min–max), mate 1:

| pile | pairs/seed | median depth | fully unwitnessed | depth <5 |
|---|---|---|---|---|
| identified, rRNA | 39,288 | 2,992 (2,876–3,162) | 0.68% (0.64–0.73) | 2.76% (2.67–2.88) |
| identified, other | 37,812 | 86.6 (85–90) | 4.50% (4.25–4.99) | 17.98% (17.62–18.69) |
| residue, on an unidentifiable contig | 8,386 | 6.0 (6–6) | 10.82% (9.90–11.71) | 44.80% (43.48–45.19) |
| residue, on a contig outside the split | 475 | 7,086 (5,047–10,620) | 6.88% (6.13–9.00) | 18.56% (16.73–21.46) |
| **residue, unplaced by anything** | **14,039** | **1.0 (1–1)** | **32.49% (30.92–35.61)** | **91.46% (90.89–92.29)** |

Mate 2 reads the same shape at its own scale (62 scorable 31-mers against 76), and both are
in the JSONs; a pooled figure is not given, per the conventions.

**The unplaced residue sits at a depth of one.** Over 90% of it occurs fewer than five times
in 8.7 M pairs and about a third occurs nowhere else at all, against 85–90 for reads the
cascade identified and ~3,000 for ribosomal ones. The two piles of the residue are an order
of magnitude apart from each other (6 against 1), which is why §4's split of the residue
matters and why "the residue" as one number hides the thing that governs it.

This explains §2 arithmetically rather than by inference: **you cannot assemble sequence you
hold one copy of.** It also places the residue against the walls already measured. A planted
genome's reads find each other by neighbourhood down to about 1 in 100,000 and not below
(`results/residue_structure/` §4); at 8.7 M pairs, 1 in 100,000 is ~87 pairs, and this
material is at 1. **The unplaced residue lies below the floor of every grouping method this
project has measured, and it is not a tuning problem — it is the sampling.**

*Polarity:* depth is reported as abundance and nothing else. The coverage-shape family was
measured to carry *familiarity* polarity (`results/cov_step_conservation/`), so a deep read
is not a more interesting read, and none of this is a novelty score.

*Incidental reproduction:* seed 1 was first computed from a single-block instrument run and
then recomputed from the five-block run; every figure above is identical to the published
digits, so the per-block numbers do not depend on which blocks share a pass.

## 4. The pile outside the split, reported rather than folded in

475 residue pairs per seed match only one of the 7,944 contigs that the cascade's three-way
split omits, and they sit at median depth ~7,000 — repetitive sequence, 30× deeper than
anything else in the residue. Folding them into either neighbour would move that neighbour's
median, which is why they have a row. They are 2% of the residue and worth someone's
attention: highly abundant own-assembly sequence that no index identified.

## 5. The spreads are five-block, not random

The five seeds are five contiguous regions of one flow cell, not a random sample
(`results/block_position_bias/`, companion thread; re-derived in
[`CONTIGS.md`](CONTIGS.md) §4 and `block_dispersion.json`). Location is unbiased - the blocks
sit within 0.16 points of the whole library on every class - but between-block dispersion runs
**1.44x and 1.61x** binomial on the rRNA and identified classes against **1.09x and 1.24x**
on the residue classes. So the residue ranges above are close to a simple random sample, and
the identified-pile ranges in §1 and §3 are narrow by about half again. No number changes and
nothing here turns on those widths.

## 6. What this says about where signal could live

- **Artefact rejection is not the lever.** Under 1% of the residue is artefact by any
  signature tested, so a better artefact filter buys at most that.
- **Depth is the governing variable, and it is not a parameter we control.** The unplaced
  residue is at depth 1. Anything that needs two observations of a molecule — assembly,
  neighbourhood grouping, mate-linkage kinship, recurrence within a run — has nothing to
  work with on that pile. The 8.4 points of the residue that sit on unidentifiable contigs
  at depth 6 are a different proposition and are where contig-level work can reach.
- **The honest target follows.** For the unplaced pile the deliverable can only be a
  per-molecule description with its own calibration — not a grouping, not a call — which is
  what SPEC.md already says. What this census adds is the reason: not that the axes are weak,
  but that most of the pile is a single observation.
- **Open, and handed on.** Coding density, Rfam RNA families and Pfam/RdRp over the 3,695
  residue-only contigs are the remaining composition questions; the contigs are the only
  place where a profile search has power, and that work sits with the companion thread.

## Provenance

| file | script | what |
|---|---|---|
| `artefacts_s{1..5}.json` | `experiments/residue_artefacts.py` | §1, 5 seeds |
| `assembly_manifest.json` | `experiments/residue_assembly.py` | §2 assembly |
| `depth_s{1..5}.json` | `experiments/residue_depth.py` | §3, 5 blocks |

Inputs not in the repo: the cascade's residue and subtracted FASTAs and contig split (from
the companion thread, via S3), `data/work/library_labels.npz` and
`library_label_classes.npz` (this thread, shared back), SILVA 138.2, UniVec_Core,
NC_001422.1. Each result carries its seed, git SHA, dirty flag and script content hash.
