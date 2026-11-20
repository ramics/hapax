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
rather than by the sewage, and (4) is bad news for the conjunction rule.** The idea
survives as a library-design recommendation with a measured effect size, not as a way
to re-use public data.

Reproduce with (all five seeds each):

```
python3 experiments/insert_profile.py  --seed S --r1 … --r2 … --ref … --pairs 200000
python3 experiments/mate_separation.py --seed S --contigs data/ww_berlin/spades_sample.fa
python3 experiments/peptide_absence.py --seed S --k 10 [--exclude-family F] --genomes …
python3 experiments/novelty_blocks.py  --seed S --accession A --exclude-family F
python3 experiments/mate_rho.py        --seed S --r1 … --r2 … --ref … --faa … --pairs 40000
```

Every number is marked **measured** (a script here, 5 seeds unless stated), **prior
run** (committed elsewhere in this repo), **literature** (cited), or **model**
(closed-form arithmetic). Data: SRR23973193 (Wyler et al. 2024, Berlin wastewater,
PRJNA948850/GSE228220) and the authors' own SPAdes contigs for that sample; RefSeq
viral proteins and Swiss-Prot for the protein index; GenBank reference genomes.

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

ρ_fp on real pairs, viral-only index, 40,000 pairs (**measured**, seed 1):

| separation | n | r (continuous score) |
|---|---|---|
| gap 0 (two halves of one mate) | 39,885 | 0.793 |
| insert 107–150 | 2,558 | 0.975 |
| insert 150–200 | 3,996 | 0.904 |
| insert 200–250 | 3,283 | 0.910 |
| insert 250–300 | 1,939 | 0.874 |
| insert 300–400 | 1,456 | 0.883 |
| insert 400–600 | 479 | 0.855 |
| gap ∞ (shuffled pairing) | 39,928 | −0.006 |

**This contradicts my own prediction.** The feature map estimated that absence-based
scores would have materially *lower* ρ_fp than homology scores, because error- and
low-complexity-driven absence is mate-specific. Measured, ρ_fp is 0.85–0.98 and barely
moves with separation, which puts the conjunction's merit at ≈1.0 - i.e. pointless - on
the model's own table. Withdraw that estimate.

One limitation keeps this from being the final word, and it is the reason for the
second arm: against a viral-only index **88% of background reads score exactly 1.0**
(**measured**), because wastewater is mostly bacterial and nothing matches whatever the
read is. A saturated score can look correlated for the wrong reason, and the thresholded
ρ is degenerate (every percentile lands on the same threshold of 1.0). The arm under
`results/mate_rho/viral_sprot/` repeats this with Swiss-Prot added to the index so the
score has range on real background; §4.1 reports it, and until it is in, the
"conjunction is dead" reading is provisional.

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

1. ρ_fp staying near 0.9 after §4.1's index fix *and* the disjunction rule gaining
   nothing from decorrelation either - at that point width improves only the single-read
   hit rate, not the pair logic, and a 1.1–1.5× factor does not justify a bespoke prep.
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
