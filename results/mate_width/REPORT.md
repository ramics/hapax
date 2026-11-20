# Does a wide library template help you see novelty?

The hypothesis, in Imo's words (2026-11-20): *"quite a wide template with far apart
mates would give max chance of hitting a coding region and max chance of novelty."*

This tests it instead of assuming it, and splits it into the four things that have to
be true for it to be worth doing:

1. the material and the protocol can deliver wide inserts at all (§1);
2. mate separation raises the chance of hitting **coding** sequence (§2);
3. mate separation raises the chance of hitting **novel** sequence, and decorrelates
   the two mates' novelty (§3);
4. what that buys survives the mate correlation of **false** flags (§4).

Headline: **(2) is false, (3) is true and measured, (1) is capped by the library kit
rather than by the sewage, and (4) finds the both-mates rule worthless at public-database
index coverage but worth 1.8–2.6× once the score has range, with width helping both
rules.** The idea survives as a library-design
recommendation with a measured effect size, not as a way to re-use public data. And at
Imo's operating point of 1 false flag per million, no pair rule gets there at any insert
length (§4.1): the pair buys sensitivity, not orders of magnitude of specificity.

Chasing (4) opened a second half of the report that was not in the plan, and it is the
half with the more portable findings. §4.4–§4.7: what the still-unmatched reads actually
are; the absence curve against real UniRef50 indexes up to 1.73 G peptide 10-mers, which
**has not flattened**; and the discovery that a third of this library is ribosomal RNA
which a large public protein index scores as *less* novel than everything else, because
about 63 misannotated public entries match ribosomal reads across the whole tree.

Reproduce with (all five seeds each):

```
python3 experiments/insert_profile.py  --seed S --r1 … --r2 … --ref … --pairs 200000
python3 experiments/mate_separation.py --seed S --contigs data/ww_berlin/spades_sample.fa
python3 experiments/peptide_absence.py --seed S --k 10 [--exclude-family F] --genomes …
python3 experiments/novelty_blocks.py  --seed S --accession A --exclude-family F
python3 experiments/mate_rho.py        --seed S --r1 … --r2 … --ref … --faa … --pairs 40000
python3 experiments/subsample_uniref.py --every 10 100        # streams, never stores
python3 experiments/absence_residual.py --seed S
python3 experiments/rrna_contigs.py                           # labels the assembly
python3 experiments/absence_curve.py   --seed S
python3 experiments/rrna_match_source.py --seed 1 --rrna-tsv …/rrna_s1.tsv
```

Every number is marked **measured** (a script here, 5 seeds unless stated), **prior
run** (committed elsewhere in this repo), **literature** (cited), or **model**
(closed-form arithmetic). Data: SRR23973193 (Wyler et al. 2024, Berlin wastewater,
PRJNA948850/GSE228220) and the authors' own SPAdes contigs for that sample; RefSeq
viral proteins, Swiss-Prot and streamed UniRef50 subsamples for the protein index;
SILVA 138.1 SSU+LSU NR99 for ribosomal labelling; GenBank reference genomes.

---

## 1. What inserts the library actually delivers, and what it costs to ask for width

This had never been measured here. The ~250 nt figure in earlier write-ups was an ART
simulation parameter in the abundance sweep, not a property of the library (flagged as
an open caveat in `results/readlevel_features/REPORT.md` §6).

Two independent estimators, 200,000 pairs per seed, mate 2 trimmed of its 14 nt
non-template prefix first (**prior run**, `results/ww_berlin_consensus/REPORT.md`):

| | median | p75 | p95 | p99 | max |
|---|---|---|---|---|---|
| mapping to the authors' contigs ≥2 kb | **207.8 ± 1.3** | 263 ± 2.4 | 387 ± 5.2 | 515 ± 8.6 | 1,017–1,158 |
| mate overlap (reference-free, sees <200 only) | 147 | 164 | 177 | 179 | 180 |

(**measured**, 5 seeds, mean ± sd. The overlap estimator is truncated by construction
at read1+read2 = 200 nt, so it is a check on the left half, and it agrees there.)

**38.4% ± 0.5% of all pairs have a detectable mate overlap** (**measured**), and 55.6%
of the pairs whose insert could be measured are under 200 nt. That is the same regime
NAO reports for other wastewater runs (58–76% of viral pairs merging at 2×150,
**literature**), and it is why the merged-mate mechanism was dropped as having no
headroom.

Short contigs cannot hold long inserts, so they would truncate the tail. Restricting
the reference from all 82,947 contigs to the 1,416 above 2 kb moves the median from 194
to 208 and p95 from 363 to 392 (**measured**): the shape is real, not an artefact of
contig length.

**The cost of asking for width, on this library** (fraction of measurable molecules
retained by a minimum-insert requirement, **measured**, 5 seeds):

| minimum insert | 150 | 200 | 250 | 300 | 350 | 400 | 500 | 600 |
|---|---|---|---|---|---|---|---|---|
| molecules retained | 0.825 | 0.541 | 0.302 | 0.155 | 0.084 | 0.045 | 0.0123 | 0.0038 |

At one read in a million, molecules sampled is the binding resource (**prior run**,
§1 of the feature map), so *selecting* wide pairs out of an existing run costs one to
two orders of magnitude of sensitivity to buy the effect measured in §3. That route is
dead. The live question is whether a library can be *built* wide, which is §1.1.

### 1.1 Can the wet lab deliver it? The kit, not the sewage, is the cap

A literature pass (sources in `PRIOR_ART_INSERTS.md`) changes the framing I started
with. I assumed wastewater RNA was too fragmented for wide inserts. Below ~1 kb that is
not the binding constraint:

- Wyler et al. used the **SMARTer Stranded Total RNA-Seq v3 Pico** kit; its manual
  states libraries span 200–1,000 bp with a maximum at 300–400 bp and that prep adds
  150 bp, i.e. an insert mode of 150–250 nt (**literature**). Our measured median of
  208 is the kit's own output distribution.
- Illumina's TruSeq Stranded Total RNA manual tabulates measured inserts against
  fragmentation time: **with fragmentation skipped entirely, median insert is still
  200 bp (range 130–350)** on intact high-RIN RNA, because random-hexamer priming is
  internal (**literature**). Every total-RNA kit in wastewater use is random-primed, so
  400 nt+ inserts are not reachable by tuning fragmentation.
- Contradicting the "too fragmented" story: encapsidated viral RNA in sewage has T90 of
  8–12 days, 1,200 nt tiled amplicons work on wastewater (if less evenly than 400 nt
  ones), and in this very dataset the authors argue a plant virus behaved as intact RNA
  while SARS-CoV-2 was already fragmented (**literature**). Material is not uniformly
  short; the protected fraction is a minority of copies (<30%, **literature**) but it
  exists.
- Ceiling: Illumina documents tunable median inserts to 550 bp on the DNA side, cores
  decline much past 700 bp, clustering is biased against long fragments, and ≥1 kb is
  poorly supported (**literature**). Nobody has published a long-insert or mate-pair
  library for an RNA virome, and nobody has published an insert distribution for a
  wastewater virome at all.

So the honest statement is: **a 400–550 nt insert library of wastewater RNA is a
bespoke prep (full-length RT, DNA-side size selection, ideally on a nuclease-protected
fraction) that nobody has built; 1 kb is not available on Illumina.** The claim I made
in this thread's first reply - that the sample type has no wide templates to give - was
too strong, and is withdrawn in that form.

---

## 2. Coding hits: the gain is real but negligible (hypothesis half 1 fails)

Slide a template across real sequence, put a 107 nt mate at its 5' end and a 93 nt mate
at its 3' end, and ask how often at least one mate has ≥60 nt inside a CDS. Annotation
is GenBank CDS for references and Prodigal (metagenomic mode) for 300 of the Berlin
sample's own contigs ≥3 kb, which is the composition a wastewater method actually meets.

Berlin contigs, coding density **0.9227** pooled (**measured**, 5 seeds):

| insert | 107 | 150 | 200 | 250 | 300 | 400 | 500 | 1,000 | 3,000 |
|---|---|---|---|---|---|---|---|---|---|
| P(≥1 mate coding) | 0.9454 | 0.9647 | 0.9793 | 0.9863 | 0.9896 | 0.9926 | 0.9927 | 0.9956 | 0.9948 |
| ± sd | 0.0027 | 0.0025 | 0.0024 | 0.0025 | 0.0021 | 0.0021 | 0.0019 | 0.0020 | 0.0044 |
| P(both coding) | 0.9344 | — | 0.9143 | — | 0.9062 | — | 0.9028 | 0.8904 | — |

Reference genomes agree, with the ceiling set by each one's coding density: astrovirus
0.976 → P(≥1) reaches 1.000 by insert 150; SARS-CoV-2 0.979 → 1.000 by 400; phiX174
0.960 → 1.000 by 250; *E. coli* 0.863 → 0.965 at 3,000; HIV-1 0.833 → 0.952 at 500
(**measured**).

Two things to take from this. First, the whole available gain is bounded by the
non-coding fraction: going from the real median insert (208) to 500 nt buys **1.3
percentage points** on the Berlin contigs, and going to infinity buys under 2. Second,
**P(both mates coding) falls** as the mates separate, 0.934 → 0.903, so any rule that
needs both mates coding is hurt by width. Half 1 of the hypothesis is true in direction
and too small to pay for.

---

## 3. Novelty hits: the gain is real and large (hypothesis half 2 holds)

Novelty is not codingness, and it behaves completely differently. "Novel" here means
what a detection method would have to see: with the target's whole family removed from
the protein database, which positions still have a retained homolog and which have
none. Built per residue with phmmer (HMMER3 via pyhmmer, E ≤ 1e-3) of each annotated
protein against RefSeq viral proteins minus that family, mapped back to genome
coordinates through the CDS location (joins and complement strands handled), then the
same template sweep as §2.

Four genomes, four family holdouts, 5 seeds each (**measured**):

| target (family removed) | coding | novel share of coding | novel blocks: n / median / max (nt) |
|---|---|---|---|
| rabies NC_001542 (*Rhabdoviridae*) | 0.909 | 0.276 | 8 / 381 / 891 |
| HIV-1 K03455 (*Retroviridae*) | 0.832 | 0.839 | 7 / 369 / 2,568 |
| astrovirus NC_001943 (*Astroviridae*) | 0.975 | 0.635 | 3 / 1,197 / 2,864 |
| SARS-CoV-2 NC_045512 (*Coronaviridae*) | 0.978 | 0.775 | 14 / 580 / 9,848 |

Novelty is blocky but the blocks are **not** simply whole genes: rabies' phosphoprotein
and M2 keep no homolog at all, its L protein is 92.6% covered and its glycoprotein
58.4%, so there is real within-gene patchiness, and the median novel run is 369–1,197
nt depending on the genome.

P(at least one mate in novel sequence), and the correlation of the two mates' novelty:

| insert | rabies | HIV-1 | astrovirus | SARS-CoV-2 |
|---|---|---|---|---|
| 107 | 0.244 / ρ 0.977 | 0.704 / ρ 0.966 | 0.625 / ρ 0.986 | 0.762 / ρ 0.993 |
| 200 | 0.290 / ρ 0.725 | 0.760 / ρ 0.737 | 0.657 / ρ 0.846 | 0.781 / ρ 0.904 |
| 300 | 0.319 / ρ 0.578 | 0.809 / ρ 0.536 | 0.671 / ρ 0.752 | 0.800 / ρ 0.798 |
| 500 | 0.367 / ρ 0.348 | 0.875 / ρ 0.252 | 0.697 / ρ 0.610 | 0.818 / ρ 0.677 |
| 700 | 0.391 / ρ 0.247 | 0.906 / ρ 0.132 | 0.720 / ρ 0.468 | 0.824 / ρ 0.642 |
| 1,000 | 0.400 / ρ 0.263 | 0.940 / ρ −0.044 | 0.753 / ρ 0.244 | 0.848 / ρ 0.525 |
| 2,000 | 0.433 / ρ 0.001 | 0.942 / ρ −0.138 | 0.843 / ρ −0.170 | 0.875 / ρ 0.332 |

(**measured**, 5 seeds; seed sd on P(≥1) is 0.001–0.005 and on ρ is ≤0.010 throughout,
so the trends are far outside the spread.)

Both halves of Imo's instinct are visible here:

- **More novelty per pair.** Going from the overlap regime to a 500 nt insert raises the
  per-pair chance of catching novel sequence by 50% relative in rabies (0.244 → 0.367),
  24% in HIV, 12% in astrovirus, 7% in SARS-CoV-2. **The gain is biggest exactly where
  novelty is a patchy minority of the genome**, which is the hard case; where most of
  the genome is already novel there is little headroom.
- **Two more-independent looks.** The mates' novelty correlation falls from ~0.98 at the
  overlap regime to 0.25–0.68 at 500 nt. This is the first number in the project for the
  correlation that every pair-level rule's merit depends on, measured rather than
  assumed, and insert length is the only knob that moves it.

Note the contrast with §2: codingness decorrelates to ~0 by insert 300 (its φ on the
Berlin contigs goes 0.89 → 0.22 → 0.03), while novelty is still correlated at 0.35–0.68
at 500 nt. Novelty blocks are longer than untranslated gaps, so novelty is the quantity
for which separation is the scarce resource.

---

## 4. The false-flag side, and a prediction of mine that was wrong

A pair rule's merit needs both correlations: ρ_true (§3) and **ρ_fp**, the correlation
of *false* flags between mates. Every read in SRR23973193 is background, so a flag on
any of them is a false flag and the correlation measured between real mates is ρ_fp
directly. Score: six-frame peptide 10-mer absence, the only single-read feature with
novelty polarity (`results/readlevel_features/REPORT.md` A12), taking the **minimum
over the six frames** - the pooled-over-frames version is dominated by the five
nonsense frames and is not usable (measured: a genome whose own proteins are in the
index still scores 0.70–0.75 pooled, versus 0.083 best-frame;
`results/peptide_absence/`).

**Scorer validation** (**measured**): rabies against an index containing its own
family scores mean absent fraction 0.083 with 89.4% of windows fully present; against
the index with *Rhabdoviridae* removed it scores 0.9945 with 94.6% of windows fully
absent. The feature has the right polarity and responds to the holdout.

**But exact peptide 10-mer absence saturates for a novel family.** With the family
removed, 94.6% of rabies windows have *every* peptide 10-mer absent (k=8: 32.0%; k=6:
0.0%, mean 0.251) (**measured**). At the divergence that matters - the nearest retained
relative is 27–35% identical (**prior run**) - exact 10-mers are all absent, so the
score cannot grade novelty at k=10 and cannot distinguish "novel" from "known but
divergent". That is a real constraint on shortlist item 2, and it says any usable
version of A12 needs either a smaller k, graded/inexact matching, or the published
rarity-index weighting rather than a plain absence count.

ρ_fp on real pairs with mate identity kept, 40,000 pairs per seed. Two indexes: RefSeq
viral alone (87.6M peptide 10-mers) and RefSeq viral + Swiss-Prot (221.1M), the second
run because of the saturation noted below. Pearson r on the continuous score
(**measured**, 5 seeds for the combined index, mean ± sd):

| separation | n | r, viral+Swiss-Prot | r, viral only (seed 1) |
|---|---|---|---|
| gap 0 (two halves of one mate) | 39,879 | 0.783 ± 0.004 | 0.793 |
| insert 107–150 | 2,531 | 0.964 ± 0.007 | 0.975 |
| insert 150–200 | 4,043 | 0.908 ± 0.009 | 0.904 |
| insert 200–250 | 3,381 | 0.898 ± 0.007 | 0.910 |
| insert 250–300 | 1,962 | 0.885 ± 0.014 | 0.874 |
| insert 300–400 | 1,416 | 0.872 ± 0.014 | 0.883 |
| insert 400–600 | 478 | 0.848 ± 0.018 | 0.855 |
| **gap ∞ (shuffled pairing)** | 39,923 | **0.001 ± 0.005** | −0.006 |

**This contradicts my own prediction.** The feature map estimated that absence-based
scores would have materially *lower* ρ_fp than homology scores, because error- and
low-complexity-driven absence is mate-specific. Measured, ρ_fp is 0.85–0.96 and falls
only slowly with separation. Withdraw that estimate.

**The saturation objection, and its answer.** Against viral+Swiss-Prot, **86.3% of
background reads still have every peptide 10-mer absent** (88.3% viral-only,
**measured**), because no index available here covers environmental protein space. A
correlation over a saturated score could just be the correlation of "neither mate
matches anything", and the thresholded ρ is degenerate (every percentile lands on a
threshold of 1.0). So the same statistics were repeated on the pairs where the score is
below 1.0 and has range (**measured**, 5 seeds):

| subset | n | r | same subset, mates re-paired |
|---|---|---|---|
| at least one mate non-saturated | 7,453 | **0.770 ± 0.009** | −0.004 ± 0.012 |
| both mates non-saturated | 3,447 | **0.773 ± 0.006** | −0.009 ± 0.013 |

The correlation survives at 0.77 against a 0.00 floor from its own shuffled control, so
the high ρ_fp is a property of the flags, not of the ceiling. And the mechanism is the
one the feature map named: **these false flags are organism-level, not read-level.** Both
mates come off the same molecule of the same under-sequenced organism, so both are
equally unknown to the index. Accidents would have been independent; organisms are not.

Raw per-pair scores are kept in `pair_scores_s*.npz` beside each result so further
questions need no re-run.

### 4.1 What this does to the pair rules at 1 false positive per million

Imo fixed the operating point at **1 false flag per million reads** (2026-11-20).
Feeding the measured ρ_true (§3) and ρ_fp (above) into the two rules gives
(`experiments/pair_rule_merit.py`, computation over measured inputs, per-read FP 1e-5
and per-read sensitivity 0.5 on a read that is in novel sequence):

| genome | insert | ρ_true | ρ_fp | either-mate: sens / FP per M | both-mates: sens / FP per M | both-mates merit |
|---|---|---|---|---|---|---|
| rabies | 107 | 0.977 | 0.974 | 0.123 / 10.3 | 0.121 / 9.7 | 1.02 |
| rabies | 300 | 0.578 | 0.884 | 0.193 / 11.2 | 0.126 / 8.8 | 0.89 |
| rabies | 500 | 0.348 | 0.848 | **0.243** / 11.5 | 0.124 / 8.5 | 0.80 |
| HIV-1 | 500 | 0.252 | 0.848 | 0.601 / 11.5 | 0.274 / 8.5 | 0.74 |
| astrovirus | 500 | 0.610 | 0.848 | 0.416 / 11.5 | 0.280 / 8.5 | 0.95 |
| SARS-CoV-2 | 500 | 0.680 | 0.848 | 0.475 / 11.5 | 0.344 / 8.5 | 0.99 |

Three conclusions, in order of how much they matter:

1. **The both-mates rule is dead for organism-level absence, which is what dominates at
   every threshold measured here.** Merit (false-flag reduction per unit of sensitivity
   given up) is 0.74–1.02, i.e. at or below the break-even of 1.0: at the measured ρ_fp
   it suppresses false flags by only 1.05–1.2× while cutting sensitivity by more. This
   is the §3.1 model evaluated at its measured parameter, and it agrees with NAO having
   chosen the disjunction empirically. §4.2 gives the measured pair rate and reconciles
   this with the error-floor result of PR #11, which is about a different term and is
   not in conflict.
2. **Either-mate plus a wide insert is the combination that pays.** Pair sensitivity
   rises by a factor of 1.2–1.7 from insert 107 to 500 while pair false flags stay flat
   (10.3 → 11.5 per million), because the gain comes from ρ_true falling, not from the
   threshold moving.
3. **Neither rule reaches 1 per million, and no insert length changes that.** Under an
   either-mate rule the pair false-flag rate is roughly twice the per-read rate, so the
   pair level hits 1/M only if the per-read rate is already ~5e-7. **The mate pair buys
   sensitivity at a fixed operating point; it does not buy orders of magnitude of
   specificity.** Reaching 1/M has to come from aggregation across molecules -
   cross-sample recurrence (C3) or grouping (C2) - with the pair as a multiplier on
   that rather than a substitute for it.

### 4.2 Reconciliation with the error floor (PR #11)

The peptide-absence thread measured the *instrument-error* component of a false flag and
found, correctly, that it has **zero** mate correlation by construction: an error in mate
1 and an error in mate 2 are separate events, so that component squares exactly under the
conjunction. At SRR23973193's measured error rate it takes 26 of 26 peptide 10-mers absent
to hold 1 FP per million on a single read, and 16 of 26 on a pair.

Both results are right and they are about different terms. A background read's absence has
two causes: **instrument error** (independent between mates, squares) and **organism-level
absence** - an under-sequenced known strain, a low-complexity or fast-evolving protein, a
non-coding fragment (shared between mates, does not square). The error term sets a floor;
the organism term sets the binding constraint, and it is the one measured here.

How binding, measured directly on real Berlin pairs rather than modelled (**measured**,
5 seeds, viral + Swiss-Prot index):

| flag threshold | per-read rate | **both mates** | if mates were independent | either mate |
|---|---|---|---|---|
| every peptide 10-mer absent | 0.863 | **0.813** (813,300/M) | 0.745 | 0.914 |
| ≥95% absent | 0.885 | 0.844 | 0.784 | 0.923 |
| ≥90% absent | 0.901 | 0.867 | 0.812 | 0.930 |
| ≥80% absent | 0.920 | 0.892 | 0.847 | 0.939 |

So at the *strictest* threshold this score can express, the measured pair-level false-flag
rate on real background is **8.1 × 10^5 per million**, 1.09× the independence prediction.
The error-floor budget of "16 of 26 absent on a pair" is a more lenient threshold than
"26 of 26", and real background crosses it in 89% of pairs. Error independence is real and
cannot rescue this: the conjunction squares a term that is not what produces the false
flags.

The honest qualification on my side: this measurement is against RefSeq viral + Swiss-Prot
(221M peptide 10-mers), which is a small fraction of real protein space, so **the organism
term measured here is inflated by index incompleteness** and a complete index would lower
it - by how much is unmeasured. The structural point survives that: any index deployable on
an environmental sample is incomplete for that sample, and incompleteness is shared between
mates off the same molecule, which is exactly the correlation that defeats the conjunction.
What would settle the size of the term is the same score against a complete known-protein
index (UniRef or nr), which neither thread has run.

Agreed wording with PR #11: the conjunction squares the error term and does nothing for
the organism term; at measured thresholds on real background the organism term dominates,
so the conjunction does not deliver a usable false-positive reduction **at the index
coverage available from public databases**. §4.3 shows that last clause is load-bearing.

### 4.3 How much of this is the index? An index-coverage arm, and a revision

The saturation caveat above deserved a measurement rather than a caveat, so here is one.
The question: is 86% of background being fully absent a statement about the *reads* or
about the *database*? The sharpest available test is to index the proteins predicted from
**the sample's own assembly** - 79,667 proteins, 8.0 Maa, from Prodigal over the authors'
72,877 contigs. That is the one "database" guaranteed to contain this sample's abundant
organisms. It is circular as a detector and says nothing about what a real pipeline could
know; as a decomposition of the absence rate it is decisive.

Adding it grows the index by 3% in k-mers and changes everything (**measured**, 5 seeds,
40,000 pairs per seed):

| index | 10-mers | background fully absent | mean absence | r (all pairs) | r (insert 107–150) | r (insert 400–600) |
|---|---|---|---|---|---|---|
| RefSeq viral | 87.6M | 0.883 | 0.941 | 0.878 | 0.975 | 0.855 |
| + Swiss-Prot | 221.1M | 0.863 | 0.938 | 0.860 | 0.964 | 0.848 |
| **+ the sample's own assembly** | **228.0M** | **0.222** | **0.486** | **0.650** | **0.855** | **0.284** |

Three things follow, and one of them revises §4.1.

**1. The 86% absence is database incompleteness, not novelty.** 8 Maa drawn from this
sample removes three quarters of it, where 200 Maa of Swiss-Prot removed two points. So
"absence in wastewater means wastewater" (the protein thread's phrasing) is sharper than
that: it means *not in public protein databases*, and the organisms are sitting in the
sample's own contigs.

**2. The mate correlation becomes more visible, not less.** In ratio terms the excess over
independence grows as the score gains range: both mates fully absent is 3.26× the
independent prediction with the sample-assembly index (0.161 vs 0.049), against 1.09× with
the public-only index (0.813 vs 0.745). The correlation was always there; saturation was
hiding it by pushing every rate towards 1.

**3. REVISION to §4.1: the both-mates rule is not dead in general - it is dead at the
index coverage public databases give you.** With range in the score, ρ_fp falls steeply
with separation (0.855 at insert 107–150 → 0.284 at 400–600), and the conjunction's merit
goes above break-even (`pair_rule_merit_sample_index.json`, per-read FP 1e-5):

| insert | ρ_true (rabies) | ρ_fp | both-mates merit | either-mate pair sensitivity |
|---|---|---|---|---|
| 107 | 0.977 | 0.870 | 1.14 | 0.123 |
| 300 | 0.578 | 0.351 | **2.25** | 0.193 |
| 500 | 0.348 | 0.380 | **1.77** | 0.243 |

So at better index coverage **a wide insert helps both rules**: the either-mate rule through
ρ_true and the both-mates rule through ρ_fp. That strengthens the wide-insert case rather
than weakening it, and it was my own flagged caveat that turned out to carry the verdict.
What does *not* change: the absolute rate. At the strictest threshold the pair false-flag
rate is still 160,774 per million (231,230 at ≥80% absent), five orders of magnitude from
Imo's target, so no pair rule rescues peptide-10-mer absence as a detector.

**A retraction about the environment, which was my error and not the environment's.**
An earlier version of this section reported that the clean version of this arm - the same
measurement against a comprehensive public index - could not be run here, because the
egress proxy would not carry UniRef50's 9.2 GB to disk. The file cannot be *stored*, and
that part is true: one attempt arrived with the gzip stream broken, a second was refused,
and a 512 MB-chunked ranged download failed permanently on chunk 3 of 18. But it
*streams*, and a subsample never needs the file. `experiments/subsample_uniref.py` streams
it and keeps nested every-10th and every-100th subsamples, and the completeness check is
decisive rather than hopeful: UniProt's REST API reports 41,243,524 clusters for this
release (`x-total-results` for `query=identity:0.5`), which is exactly the record count the
stream reached - including on the first attempt, whose output had been discarded as
prefix-biased on the strength of a `BadGzipFile` raised *after* the last record. The lesson
is the general one for this project: a failure at the end of a stream is not evidence about
its middle, and the cheap check (ask the source how many records there should be) was
available the whole time.

§4.4 below is that measurement. The remaining limit is real but much smaller: the far end
of the curve, a full UniRef50 or nr index, is ~160 GB of sorted hashes and still does not
fit, so the largest point here is a tenth of UniRef50 rather than all of it.

---

---

### 4.4 What the residue is, once the index stops being the explanation

§4.3 left 23.0% ± 0.4% of background reads with no known peptide 10-mer, and that residue
is what an absence-based detector would actually flag, so it is worth knowing what it is.
Two mechanisms predict it and they point opposite ways:

1. **rare and unsequenced** - organisms neither in public databases nor abundant enough to
   assemble here. If this dominates, the feature separates "abundant or known" from "rare
   and unknown", the detection target sits on the same side as the false positives, and no
   index fixes it.
2. **abundant but non-coding** - a read off ribosomal RNA, a UTR or a structural RNA can
   never match a protein index however complete, and in a total-RNA library rRNA is the
   most abundant thing present, not the rarest.

Both are testable from data already on disk, which is why this cost an afternoon rather
than a run: every read is mapped to the sample's own assembly, and SPAdes writes per-contig
coverage into the contig name, so abundance comes free with the mapping
(`experiments/absence_residual.py`, 5 seeds, 39.8k reads per seed,
`results/absence_residual/`).

| | fully absent | the rest |
|---|---|---|
| maps to the sample's own assembly | 0.391 (0.007) | 0.928 (0.002) |
| inside a predicted CDS, of those that map | 0.180 (0.012) | 0.725 (0.005) |
| contig coverage, median | 115 (7) | 140 (3) |
| on a contig below 2x | 0.055 | 0.060 |
| mate also maps | 0.442 (0.009) | 0.912 (0.003) |
| unmapped but mate mapped | 0.072 (0.002) | 0.028 (0.001) |
| adapter present | 0.0004 | 0.0001 |
| dinucleotide entropy, mean | 3.804 | 3.839 |

The residue splits two ways rather than one.

**39% maps, and where it maps it is not rare.** Median coverage 115 and only 5.5% below 2x,
indistinguishable from the rest of the library, but it sits outside gene calls at 0.18
against 0.73. That is mechanism 2's signature.

**61% does not map to the sample's own assembly at all, and the mate agrees.** Only 7.2% of
residual reads are unmapped with a mapped mate, so these are whole molecules missing from
the assembly rather than unalignable ends of abundant ones. That is mechanism 1's
signature, and it is the larger half.

**Artefacts are not the story.** Adapter content is enriched fourfold and is 0.0004 in
absolute terms; entropy and homopolymer distributions are flat between the groups. Worth
recording, because a residue defined by matching nothing is exactly where this project
expected artefacts to collect.

**My advance prediction was wrong** and is recorded as such: I predicted a mapped fraction
"well above half", on the implicit assumption that the assembly captures the library. It
captures 93% of ordinary reads and 39% of these.

### 4.5 A challenge to §4.3, and what it cost to check

The sibling rRNA thread then measured something that §4.3 cannot survive unexamined: 34.5%
of this library is ribosomal RNA by exact 31-mer containment against SILVA 138.1, and reads
called rRNA are *less* than half as likely to be fully absent as other reads (0.125 against
0.286). Ribosomal RNA has no protein product, so a clean index cannot match it. The
proposed explanation was that Prodigal - a prokaryotic caller, run here over a mixed
metatranscriptome - calls spurious genes on the highest-coverage contigs in the assembly,
which are ribosomal, and that those spurious peptides are a large part of what made the
sample-self index look so effective.

Checked directly, by labelling the assembly rather than the reads
(`experiments/rrna_contigs.py`, SILVA SSU+LSU NR99, exact 31-mer containment; a contig
offers thousands of k-mers where a 107 nt read offers 77, so the call barely depends on the
threshold):

| | |
|---|---|
| contigs called ribosomal (hit fraction ≥ 0.5) | 5,303 of 82,856 |
| share of assembly bases | 1.72% |
| Prodigal coding density on them | 0.281, against 0.831 elsewhere |
| predicted proteins on them | 388 of 79,667 |
| **indexed residues on them** | **19,137 of 8,042,848 = 0.24%** |

**What 0.24% does and does not settle.** It is a share of the *index*, so it answers
"could spurious ribosomal calls have moved 86% to 22%" - clearly not. It does not by itself
settle P(absent | rRNA) = 0.125, because a read needs exactly one matching 10-mer, so the
quantity that matters there is the chance an rRNA read falls within reach of a spurious
call. That is the coding density on ribosomal contigs, 0.281, two orders of magnitude
larger than the index share, and it agrees with the sibling thread's read-level figure of
about 0.38 for mapped rRNA reads inside a call. The gap's direction is the expected one:
partly-ribosomal contigs carry real genes beside ribosomal stretches, so read-level sits
above contig-level - which also means 1.72% of bases is a floor twice over, once for
divergent lineages and once for mixed contigs.

So the mechanism is real and small *as a share of the index*. Spurious gene calls on ribosomal contigs do enter the
index, their coding density is a third of the rest of the assembly's, and they are a
quarter of one percent of the residues that moved the absence rate from 86% to 22%. The
conclusion does not depend on the threshold: loosening it to 0.3, 0.2, 0.1 and an absurd
0.05 gives 0.62%, 0.77%, 0.96% and 1.10%. The sibling thread's own generalised version
agrees from the other side - pyrodigal over SILVA sequences calls genes at density 0.164,
and its composition-shuffled control calls *more*, so the calling is generic short-ORF
behaviour rather than something ribosomal structure provokes.

The challenge therefore stands as a correctly identified contaminant of the index and a
rejected explanation of the drop. What it leaves live is a different route to the same
observation: translated rRNA may be low-complexity and match the low-complexity stretches
of real proteins, which no cleaning of gene calls would remove. §4.6 measures that, because
every query k-mer in the curve carries its distinct-residue count.

Two limits of this labelling, stated because they bound it in opposite directions: realised
k-mer retention was 1 in 65 rather than the intended 1 in 8 (the low bits of an FNV-style
hash over a 4-letter alphabet are not uniform), so each contig had fewer chances than
designed and `hit_frac` must not be read as a coverage estimate; and exact 31-mers miss
divergent environmental rRNA on the contig side exactly as on the read side, so 1.72% of
bases is a floor.

**The nucleotide-side comparison.** The subtraction thread measures the same subtraction in
nucleotide space and finds the sample's own 35 MB assembly removes 83.6% of pairs on its
own, and 40.5% on top of 0.7 GB of Kraken2 viral plus 16 GB of PlusPF
(`results/subtract_known/`, seed 1 of 5). Both sides therefore land in the same place from
different alphabets: the sample's own contigs are worth more than any public database we
can reach, which is a statement about what has been deposited rather than about what is
novel.


### 4.6 The curve, against real public indexes

`experiments/absence_curve.py`, 5 seeds, 20,000 pairs per seed, cumulative levels so each
is a superset of the last. No index is built: the reads' peptide k-mers are hashed once and
each proteome is streamed past them, so memory is set by the reads rather than the database
and nothing of UniRef50 is stored.

| level | cumulative db 10-mers | fully absent | mean absence | both mates absent | r (mates) | excess over independence |
|---|---|---|---|---|---|---|
| RefSeq viral | 0.169 G | 0.8834 (0.0009) | 0.9386 | 0.8478 | 0.884 | 1.087 |
| + Swiss-Prot | 0.373 G | 0.8625 (0.0011) | 0.9341 | 0.8119 | 0.867 | 1.091 |
| + UniRef50, 1% | 0.496 G | 0.8019 (0.0021) | 0.9211 | 0.7160 | 0.815 | 1.113 |
| + UniRef50, 10% | 1.731 G | 0.6346 (0.0022) | 0.8716 | 0.4992 | 0.678 | 1.240 |
| + the sample's own proteins, ribosomal contigs excluded | 1.738 G | 0.2148 (0.0025) | 0.4887 | 0.1317 | 0.654 | 2.854 |
| + the 388 ribosomal-contig proteins | 1.738 G | 0.2016 (0.0030) | 0.4663 | 0.1238 | 0.639 | 3.045 |

**This revises the framing of §4.3, though not its arithmetic.** I wrote that 8 Maa from the
sample "did roughly forty times what 200 Maa of Swiss-Prot did". Per residue that is still
true and now understated. But the implied reading - that public databases are a dead end
for this subtraction - is wrong: a tenth of UniRef50 moves the fully-absent rate 22.8 points
(0.8625 → 0.6346), where Swiss-Prot moved 2.1. The earlier curve looked flat because both
its points were small and viral-weighted, which is exactly the mistake the process notes in
this project warn about: two points are not a curve.

**The curve has not flattened where we must stop.** Each tenfold step in index size is still
buying ~20 points at the largest size reachable here, so the far end - full UniRef50, or nr
- is not an extrapolation anyone should make from these three public points. That is the
headline deliverable the project's thesis asks for (background removed per GB of index, and
where it flattens), and the honest version of it is: *it has not flattened yet*.

**The sample's own assembly remains the largest single step**, and it is 8 Maa against
1,270. Both facts stand together, and together they say the same thing: absence in
wastewater is a statement about what has been deposited, not about what is novel.

**The mate correlation tracks index coverage, not the reads.** r falls 0.884 → 0.639 as
coverage grows, while the ratio excess over independence rises 1.087 → 3.045. Those move in
opposite directions because they are different quantities: saturation was suppressing the
correlation's *visibility* while inflating r. The pair-rule merit in §4.1 should therefore
be read off the best-covered index available, not the cheapest.

### 4.7 The ribosomal inversion, and where it actually comes from

Joining the sibling thread's SILVA calls onto the per-level dumps (40,000 rows per seed,
100% matched on read_id + mate) gives the sharpest diagnostic in this thread:

| level | P(absent \| rRNA) | P(absent \| not rRNA) |
|---|---|---|
| RefSeq viral | 0.998 (0.000) | 0.823 (0.002) |
| + Swiss-Prot | 0.960 (0.002) | 0.811 (0.002) |
| + UniRef50, 1% | 0.810 (0.003) | 0.798 (0.002) |
| + UniRef50, 10% | 0.473 (0.003) | 0.720 (0.002) |
| + sample, ribosomal contigs excluded | 0.143 (0.003) | 0.252 (0.003) |
| + the 388 ribosomal-contig proteins | 0.107 (0.004) | 0.251 (0.003) |

At a viral-only index ribosomal reads behave exactly as sequence with no protein product
should: 0.998 fully absent, the highest figure anywhere in this report. The inversion -
ribosomal reads being *less* absent than other reads - is created between the 1% and 10%
UniRef50 subsamples, not by our assembly. And it is not explained by either candidate:

- **spurious calls on our own ribosomal contigs**: removing all 388 of them moves the figure
  from 0.107 to 0.143. The sibling thread predicted 0.4 in advance and I predicted 0.25-0.35;
  both are wrong, and wrong the same way, so the error was not the index-share reasoning I
  had already conceded.
- **low complexity**: restricting to k-mers with at least 5 of 10 distinct residues changes
  nothing at any level - 0.144 against 0.143 for ribosomal reads, 0.261 against 0.252 for
  the rest.

What is left is the public database itself, and the mechanism the sibling thread proposed is
the right one in the wrong place: UniRef50 includes TrEMBL, which is automatically annotated,
and a gene caller run over a genome assembly makes the same mistake on ribosomal loci that
Prodigal makes on our contigs - at a scale of millions of entries rather than 388.
`experiments/rrna_match_source.py` tests that by naming the entries instead of inferring
them: every entry containing a ribosomal read's peptide 10-mer is recorded with its
description, with the non-ribosomal reads carried through the same pass as a control.

**The entries, named** (`results/rrna_match_source/`, seed 1, the whole 4,124,353-entry
UniRef50 10% subsample streamed once, non-ribosomal reads as a control in the same pass).
The sibling thread sharpened the prediction before the run: a uniform sprinkling of bad
automatic annotations cannot produce the *differential* in the table above, but conservation
can - rRNA is the most conserved sequence there is, so one spurious entry from a single
organism's ribosomal locus matches ribosomal reads across the whole tree, while a spurious
entry elsewhere matches only its own neighbourhood. That predicts a heavily concentrated
match distribution. It is:

| | ribosomal reads (13,580) | other reads (26,380) |
|---|---|---|
| entries matched | 4,209 | 36,218 |
| read-matches held by the top 1 entry | 2.97% | 0.11% |
| by the top 10 | 14.3% | 0.44% |
| **by the top 100** | **52.2%** | **1.62%** |
| by the top 1,000 | 79.5% | 10.9% |
| entries each matching ≥1% of the group's reads | 63 | 1 |
| most reads matched by one entry | 1,330 | 283 |
| description uninformative (uncharacterised/hypothetical/ORF/DUF) | 66.0% | 46.7% |

**Amplification confirmed, and the head is tiny.** Sixty-three entries out of four million
reach one percent or more of the ribosomal reads each, and a hundred entries carry over half
of all ribosomal read-matches, against 1.6% for the control. Both my advance prediction
(top-100 above half for ribosomal, under a tenth for the control) and the sibling thread's
conservation mechanism are confirmed on the same table. The annotation signal is there too -
66.0% of ribosomal-matching entries are uncharacterised against 46.7% of control entries -
but it is the weaker readout, exactly as predicted, because a spurious ORF can inherit a
plausible name: the single largest offender is "Cell wall-associated hydrolase" at 1,330
ribosomal reads.

Two caveats. The control is not matched for abundance, so some of its concentration is the
library rather than the database; the ratio between the two heads is the signal, not the
ribosomal head alone. And this arm is seed 1 only - the five-seed statement is the
per-level join table above, and nothing here depends on the one seed beyond the shape of a
distribution whose head is three orders clear of its control.

**Why this matters beyond a curiosity.** A third of this library is ribosomal, and against a
large public protein index a third of the library stops looking novel for a reason that has
nothing to do with biology. Any absence-based score calibrated on a big protein database is
partly measuring annotation artefacts, and the fix is not a better threshold: it is to
subtract rRNA explicitly, with a nucleotide index, before a protein index is consulted at
all. The cheaper mitigation, if a protein index must carry it, is a blocklist: sixty-three
entries explain most of it.

---

## 5. What follows

**The recommendation.** If a library is being designed, target the widest insert the
chemistry and the flowcell allow (400–550 nt), and do not expect it from a random-primed
total-RNA kit. The measured payoff is a 7–50% relative rise in the chance a pair
contains novel sequence, and a fall in the mates' novelty correlation from ~0.98 to
0.25–0.68, which is what makes two mates two observations rather than one. Do **not**
try to get this by filtering wide pairs out of existing runs: that costs 1–2 orders of
magnitude of molecules (§1) to buy a factor of 1.1–1.5.

**What would kill the recommendation.** Three things, none of them tested here:

1. The disjunction's measured gain (§4.1, 1.2–1.7× pair sensitivity from insert 107 to
   500) failing to survive a real per-read score with a real operating point. §4.1
   assumes per-read sensitivity 0.5 on a novel read, which nothing here measures.
2. The protected (less fragmented) fraction of sewage viral RNA turning out not to
   support 500 nt inserts in practice; §1.1 found no one has measured it.
3. Clustering bias eating the gain: Illumina's own documentation says long fragments
   cluster worse, so a library built at median 500 may sequence at median 350.

**Not established here.** No detection method is demonstrated; §3 measures the novelty
structure a method would exploit, not a method. The novelty track depends on the
database and on phmmer's sensitivity at E ≤ 1e-3 - a more sensitive search would shrink
the novel fraction and could lengthen or shorten blocks. Four genomes is four, not a
survey, and all four are human-associated RNA viruses with annotated CDS. The insert
distribution is one run from one plant with one kit.

**What the later arms add.** §4.6 and §4.7 were run after the above and change two of its
framings:

4. **Public protein databases are not a dead end, and the curve has not flattened.** A tenth
   of UniRef50 takes the fully-absent background rate from 0.862 to 0.635, where Swiss-Prot
   moved it 0.021. Anything stated here about what "unknown" means in this library is a
   statement at the index size reached, not a limit of the approach - and the biggest step
   is still the sample's own 8 Maa, which is a fact about deposition rather than novelty.
5. **A third of this library cannot be scored by a protein index at all, and scoring it
   anyway measures annotation quality.** rRNA reads are 0.998 absent against real proteins
   and 0.473 against a tenth of UniRef50, because ~63 misannotated public entries reach
   ribosomal reads across the whole tree. Any pipeline built on the features in this report
   should subtract rRNA with a nucleotide index first. This was not in the plan for this
   thread and is the most portable thing in it.

**The residue, for whoever builds the detector.** Of the 23% that still matches nothing with
the best index reachable here, 61% is not in the sample's own assembly and its mate agrees -
material too rare or too divergent to assemble - and 39% maps to abundant contigs outside
any gene call. The first is where the target lives; it is also where the irreducible false
positives live, which is §4.4's uncomfortable point and unchanged by anything since.

