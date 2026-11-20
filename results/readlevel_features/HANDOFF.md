# Handoff: what the two follow-on threads should not rediscover

The read-level work split on 2026-11-20 into two threads — peptide/protein fragment content,
and paired-read evidence. This is what has already been measured or checked, so neither
thread spends a morning re-deriving it. Everything here is backed by `REPORT.md` (measured,
5 seeds), `PRIOR_ART.md` (literature, with primary-source verification on the load-bearing
items), or `results/abundance_sweep/REPORT.md` (the earlier sweep).

## Shared ground, both threads

- **Below ~1x, assembly returns nothing**; 4.2% breadth at 1.68x, 64.5% at 5.03x. At 1e-6
  abundance in 1e6 pairs you expect **one pair**. Reaching the first assembled contig at that
  abundance needs ~1e8 pairs.
- **The budget: ≤1 false flag per million reads** is the target. Translated search's best
  non-zero point was 1,489/M at 0% sensitivity. That is the floor to clear.
- **"Is it viral" carries almost no information** in a wastewater virome (98.3% of background
  contigs ≥3 kb called viral). The question is new vs already known.
- **Do not claim the per-million unit as ours.** PACIFIC (*Sci Rep* 2021) Table 2 is headed
  `FP/1 M reads`. The surviving claim is sensitivity at a fixed per-million FP rate **on real
  background**. See `PRIOR_ART.md` §1.1.
- **Truth labels are free**: ART names spike reads after the accession, real background reads
  keep SRA names. Five seeds before quoting a number. Prefer continuous readouts to yes/no
  ones — a binary readout produced a null here that a continuous one overturned.

## For the protein/peptide-fragment thread

**Do not build the score from scratch.** Three published components:

- **Rarity index** (*NAR Genom Bioinform* 2024, PMC10993293) — a graded per-k-mer absence
  score over 45,785 genomes and 21,871 proteomes, with viral indices. This is the scoring
  primitive; take the weighting from it.
- **MetaProFi** (*Bioinformatics* 2023) — amino-acid k-mer Bloom index, six-frame queryable
  from nucleotide input, **zero false-negative rate**, which is the property an absence score
  requires.
- **k_aa = 10** is the specificity floor. At k_aa = 7, over 80% of peptide k-mers are shared
  across phyla, so a k=7 absence signal is noise.

**Already measured, do not re-measure** (`REPORT.md` §2):

- ORF-likeness has **no novelty polarity** and is weak anyway: at 107 nt, 33–72% of
  composition-matched shuffled reads are already stop-free in ≥1 of 6 frames.
- **47% of coding reads at 107 nt do not identify their own frame** (true frame unique in
  0.528 of reads; 0.943 at 250 nt). Any frame-dependent feature must marginalise over ~1.6
  candidate frames.
- Codon usage is ~35 codons for 61 categories; codon-pair bias ~34 observations for 3,721.

**Measure this first, it is unmeasured in the literature:** peptide-k-mer error
amplification. One substitution corrupts the whole codon and therefore every peptide k-mer
overlapping it — roughly k consecutive k-mers — and a frameshifting indel destroys the rest
of the frame. Error-driven absence is characterised for *nucleotide* k-mers and never for
peptide k-mers, where it is structurally worse. An absence score that has not accounted for
this is measuring its own error profile.

**Positions to know:** Kaiju already does six-frame translation with variable-length exact
amino-acid matches (not fixed 12-mers), and nobody scores its unclassified fraction as
novelty. sourmash supports protein/Dayhoff/hp encodings and `gather --output-unassigned`, but
reports unassigned fraction per **sample** and states it does not classify individual reads.
AGNOSTOS (*eLife* 2022) quantifies the unknown coding space per **gene**, after assembly.
geNomad's score calibration is the machinery a reviewer will ask why we did not extend.

**Read before claiming the slot is open:** Viromancer / Lyra-Virus (Feldman, Arze & Sabeti,
NeurIPS 2026 **workshop** oral, OpenReview `qQDFUh4sMJ`) scores individual reads under strict
family-level holdout across six viral realms. Abstract-only so far; pair aggregation and
calibration unknown.

## For the paired-read thread

**The model and its one free parameter** (`REPORT.md` §3.1). `P(both mates pass) = p² +
rho·p·(1−p)`, with separate correlations for true and false flags. At per-read sensitivity
0.3 and FP 1e-4, with `rho_true` = 0.9, the merit (FP reduction per unit sensitivity
surrendered) is 9,300 at `rho_fp` = 0, **9.3 at 0.1**, 3.7 at 0.25, 1.9 at 0.5, and 1.0 —
pointless — at 0.9. Four orders of magnitude on one unmeasured number.

**Known positions:**

- **The NAO already tried both rules** on real wastewater across thresholds 15–45 and chose
  the **disjunction** ("either mate ≥ 20"). They never computed a mate correlation. This is
  the most likely prior art to be raised, and it is a measurement of *their* `rho_fp`, not a
  refutation.
- **taxMaps** (*Genome Res* 2018) uses mate concordance as a precision proxy across seven
  real metagenomes, validated at Spearman ρ = 0.994 — never as a filter, never as a
  mate-vs-mate correlation. Its Fig 3B is the thing to calibrate against.
- **Classifiers all sum over mates, none conjoin**: Kraken 2 concatenates mates with an `N`
  and runs LCA once over the union; Kaiju gives only a disjunction rationale; Centrifuge sums
  over combined pair length; DIAMOND has no paired support at all.
- **Variant calling treats mates as fully dependent** and engineers the double-count away
  (Strelka's documented FP inflation, GATK capping overlap quality, fgbio `ClipBam`). ROVER
  (2014) is the one published conjunction, as an unanalysed design choice. **The
  partial-correlation middle is empty — that is the strongest claim either strand found.**
- **Merged mates / short inserts is dropped**, not deferred: wastewater RNA at 2×150 already
  mostly overlaps (58–76% merge by NAO's figures) and NAO's own analysis argues a positive gap
  beats any overlap. Do not revive it, and do not spend time measuring SRR23973193's inserts
  for that purpose.
- **One-sided anchor is demoted**: a 2010 primitive (NovelSeq), and ViFi (*NAR* 2018) already
  pairs the anchor with viral-family profile HMMs on the free mate. Its confound dominates —
  a dengue/Nextera study found most or all host–virus chimaeras artificial, amplicon chimaera
  rates exceed 70% for the rarest taxa (our target by construction), and they form
  reproducibly across replicates. If kept, the chimaera rate is the first experiment.

### Imo's wide-template hypothesis, assessed

> "A wide template with far apart mates would give max chance of hitting a coding region and
> max chance of novelty."

**The direction is right, and it has independent support we found rather than assumed.** The
NAO's read-length analysis argues on cost-per-informative-base grounds that a positive gap
beats any overlap — a 2×150 run on an 80 bp fragment yields only 80 bp of unique sequence.
And Kaiju's own discussion gives exactly this rationale for pairing: *"the chance of one mate
overlapping with a protein-coding gene is higher than for single-end sequencing."* A wide
template samples ~2L unique bases instead of ~L, so it doubles the sequence in which
something informative can land. This is also the same evidence that killed the merged-mate
idea, pointing the other way, so the two positions are consistent, not contradictory.

**Three qualifications, each cheap to check and worth checking before building on it:**

1. **Coding density caps the coding-hit gain.** HXB2 is ≥81% coding by contiguous CDS
   (7,884 of 9,719 nt, a lower bound since spliced CDS are excluded), and most compact RNA
   viruses are higher. A single 107 nt mate already lands in coding ~4 times in 5, so going
   from one mate to two far-apart mates raises P(at least one in coding) from ~0.81 to ~0.96
   at best. Real but modest. **The gain from a wide template is mostly the doubling of unique
   bases, not the coding-hit probability** — worth separating, because the two have different
   ceilings and only the first keeps improving with insert length.
2. **Wastewater may not offer the choice.** That 58–76% of viral read pairs already overlap at
   2×150 means the fragments in that substrate are mostly under ~300 nt. Insert width there is
   set by how degraded the RNA is, not by library design, so the hypothesis is largely
   untestable in wastewater and genuinely available in intact clinical material. This is the
   same constraint that makes nanopore unhelpful for this regime.
3. **Insert width moves BOTH correlations, in opposite directions — and this is the part that
   changes the planned experiment.** Wider templates make the mates more independent, which
   is good for `rho_fp` (lower, so the conjunction suppresses false positives far better:
   merit 1.9 → 9.3 as `rho_fp` falls from 0.5 to 0.1) but bad for `rho_true` (also lower, so
   the conjunction costs more sensitivity: the pair-sensitivity factor falls from 0.93x at
   `rho_true` = 0.9, to 0.65x at 0.5, to 0.30x at full independence). **So there is an optimum
   insert width, and `rho_fp` should be measured as a function of insert size rather than as a
   single number.** That is a strictly better version of the experiment than the one
   `REPORT.md` §3.1 proposes, and it comes directly from Imo's hypothesis.

A wide template also raises the chance that the two mates fall in *different genes*, which
matters for a one-sided anchor (more chance one mate is conserved and placeable while the
other is novel) and is the mechanism by which `rho_true` falls. Both effects are the same
underlying fact and should be measured together.

## Status of PR #8

Two commits: the measured feature map with its script, and the prior-art review with the
three corrections it forced. No CI in this repo. Reviewable as it stands; the open question
for a reviewer is whether the shortlist ordering in `REPORT.md` §5 is right, given that the
`rho_fp` measurement and the peptide-absence build are now owned by separate threads.
