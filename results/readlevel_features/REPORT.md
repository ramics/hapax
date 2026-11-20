# Read-level features for novel viral families at sub-1x coverage

What can be extracted from Illumina reads when coverage is too low to assemble, what
each feature actually measures, and what would kill it. Breadth first, then depth on
the paired-end axis.

Reproduce the measured sections with:

```
python3 experiments/readlevel_features.py --n-reads 4000 --n-seeds 5
```

Outputs `results/readlevel_features/readlevel_features.json` and `manifest.json`
(seed, git SHA, parameters, timings). Genomes come from
`/mnt/project-files/stress/community/`; CDS coordinates and sequence both from
`/mnt/project-files/ref/K03455.gb`.

Every number below is marked **measured** (this script, 5 seeds), **prior run**
(committed elsewhere in this repo, cited), **model** (closed-form arithmetic in this
script), or **estimate** (reasoned, not computed). Nothing here fits a classifier or
claims a detection method works; this bounds what the features can do, which is the
cheap half of the question and the half that kills candidates before they are built.

---

## 1. The regime, and what actually binds

Two constraints define the problem, and most feature proposals fail on one of them
before any model is trained.

**Constraint 1: the organism contributes almost no reads.** For a 12 kb genome, 200 nt
of observed sequence per pair (107 nt mate 1 plus a 93 nt trimmed mate 2, matching
SRR23973193):

| Run size (pairs) | Abundance | Pairs from organism | Coverage | P(≥2 pairs) |
|---|---|---|---|---|
| 1e6 | 1e-6 | 1.0 | 0.017x | 0.26 |
| 1e6 | 1e-5 | 10 | 0.168x | >0.999 |
| 1e6 | 1e-4 | 100 | 1.68x | >0.999 |
| 1e8 | 1e-7 | 10 | 0.168x | >0.999 |
| 1e8 | 1e-6 | 100 | 1.68x | >0.999 |

(**model**.) The assembly arm measured nothing recoverable below 1x, 4.2% breadth at
1.68x, 64.5% at 5.03x (**prior run**, `results/abundance_sweep/REPORT.md`). Reading the
two together: at one-in-a-million abundance you need ~1e8 pairs before the organism even
reaches the coverage where assembly produces its first 509 nt contig. In a 1e6-pair run
at that abundance you expect **one pair**, and a 26% chance of two.

That is the whole reason read-level work exists here, and it also sets a hard rule:
**any feature requiring two reads from the same organism is unavailable at the bottom of
the range.** Grouping features are not wrong, they are simply out of scope below ~1e-5
in a single modest run, and they come back as soon as the run is large.

**Constraint 2: the false-positive budget is brutal, and must be in reads per million.**
With ~20 true reads available, here is what a per-read false-positive rate costs
(**model**):

| Run size | FP rate | False flags | False:true (20 true reads) |
|---|---|---|---|
| 1e6 pairs | 0.1/M | 0.2 | 0.01 : 1 |
| 1e6 pairs | 1/M | 2 | 0.1 : 1 |
| 1e6 pairs | 10/M | 20 | 1 : 1 |
| 1e6 pairs | 1,489/M | 2,978 | 149 : 1 |
| 1e8 pairs | 1/M | 200 | 10 : 1 |
| 1e8 pairs | 1,489/M | 297,800 | 14,900 : 1 |

A usable method needs roughly **≤1 false flag per million reads**. The best non-zero
operating point translated search reached was 1,489/M at 0% sensitivity (**prior run**).

The gap must be stated carefully, and narrower than this project has been stating it.
PACIFIC (*Sci Rep* 2021) already reports false positives in this unit: its Table 2 columns
are headed `FP/1 M reads (Balanced)` and `FP/1 M reads (Unbalanced)`, giving 81.6 and 66.6
per million for Coronaviridae, 2 and 1 for SARS-CoV-2. So "nobody expresses false positives
per million reads" is false. What survives checking its full text is that sensitivity
(Table 1) and false positives (Table 2) come from different experiments at one fixed
threshold, with no curve from which sensitivity could be read off at a chosen per-million
rate, and that both FP tables are ART-simulated, with no real-background FP rate computed.
**Safe claim: no published read-level detector reports sensitivity at a fixed false-positive
rate per million background reads, measured on real background.** See `PRIOR_ART.md` §1.1.

### The polarity problem

This governs the whole catalogue. Almost every published read-level feature measures
*is this coding / is this viral*. In a wastewater virome that question is nearly
answered before it is asked: geNomad called 98.3% of background contigs ≥3 kb viral
(**prior run**). The discriminating question is **new versus already known**, and only
a minority of features below have that polarity at all. A feature with the wrong
polarity cannot be rescued by a better classifier on top of it.

---

## 2. Feature catalogue

Verdict column: **dead** (ruled out, reason given), **wrong polarity** (measures
viralness/codingness, not novelty), **open** (worth testing), **needs scale** (works,
but only above ~1e-5 abundance or ~1e8 pairs).

### A. Single-read, sequence-intrinsic

| # | Feature | Measures | Verdict |
|---|---|---|---|
| A1 | Translated homology search (DIAMOND blastx, BATH, MAIRA) | Distance to database | **dead** as a novelty score |
| A2 | Nucleotide k-mer composition (3–6-mers) | Organism identity | **wrong polarity** |
| A3 | ORF-likeness / stop-codon density | Codingness | **wrong polarity**, and weak at 107 nt |
| A4 | Frame identifiability | Prerequisite for A5, A9 | **measured constraint**, not a feature |
| A5 | Codon usage bias | Organism/host identity | **wrong polarity**, data-starved |
| A6 | Low-complexity / Shannon entropy | Artefact | **open** (as a filter, not a detector) |
| A7 | GC fraction, GC skew | Organism identity | **wrong polarity** |
| A8 | Dinucleotide bias (CpG, UpA) | Host adaptation | **open**, weak alone |
| A9 | Codon-pair bias | Host adaptation | **dead** at this length (data-starved) |
| A10 | RNA secondary structure propensity | Structured RNA | **dead** at this length |
| A11 | Nucleotide k-mer **absence** from a known-sequence index | Novelty | **open** — correct polarity |
| A12 | Six-frame peptide k-mer **absence** from a protein index | Novelty, divergence-tolerant | **open** — the strongest single-read candidate |
| A13 | Protein language model embedding (ESM-2) on ~35 aa | Protein-likeness, family | **open**, expensive, needs A4 |
| A14 | Quality / error profile | Artefact | **open** (as a filter) |

**A1 — translated homology. Dead, measured.** 8–10% sensitivity at 93,478 background
flags/M; 0% once background falls to 1,489/M; only 11.6% of target reads get any hit
(**prior run**). The cause is structural, not a tuning failure: a 107 nt read is ~35 aa,
the nearest retained relative was 27–35% identical, and the score therefore ranks reads
by **database proximity**, which is the opposite of novelty. Known-but-divergent
background scores above genuinely novel sequence.
*Killed by:* already killed. Retained as the floor any candidate must clear.

**A3 — ORF-likeness. Wrong polarity, and weaker than usually assumed.** Fraction of
reads with no stop codon in at least one of six frames, 5 seeds, 4,000 reads per cell
(**measured**):

| Genome | GC | Real, 107 nt | Shuffled, 107 nt | Gap | Real, 250 nt | Shuffled, 250 nt |
|---|---|---|---|---|---|---|
| *S. aureus* | 0.33 | 0.840 | 0.332 | +0.508 | 0.642 | 0.009 |
| SARS-CoV-2 | 0.38 | 0.974 | 0.437 | +0.536 | 0.906 | 0.025 |
| HIV-1 | 0.42 | 0.967 | 0.586 | +0.381 | 0.803 | 0.066 |
| rabies | 0.45 | 0.975 | 0.588 | +0.388 | 0.864 | 0.066 |
| phiX174 | 0.45 | 0.939 | 0.605 | +0.334 | 0.721 | 0.063 |
| HBV | 0.48 | 0.986 | 0.656 | +0.330 | 0.921 | 0.099 |
| phage λ | 0.50 | 0.922 | 0.704 | +0.218 | 0.694 | 0.103 |
| *E. coli* | 0.51 | 0.933 | 0.721 | +0.212 | 0.733 | 0.121 |

Shuffled is a mononucleotide shuffle of the same genome, so it is the composition-matched
null. Three things follow. First, at 107 nt even *random* sequence is stop-free in some
frame 33–72% of the time, so the feature carries well under one bit. Second, the gap
shrinks as GC rises, which is the wrong direction: high-GC bacterial background is
exactly where it fails worst. Third and fatally, **a wastewater background is real
coding sequence too**, so the relevant contrast is not real-versus-shuffled but
novel-coding-versus-known-coding, where this feature is identical by construction.
The closed-form model agrees with the measurement (at GC 0.50, per-frame stop-free
0.186 at 107 nt, ~0.71 across six frames).
*Killed by:* already dead as a novelty feature. Would only matter if someone proposed it
as one.

**A4 — frame identifiability. A constraint, not a feature.** Reads drawn inside
annotated HXB2 CDS, so the true frame is known. The true frame is stop-free in
1.0000 of reads at every length (the invariant that must hold, and the check that caught
two bugs in this script). The useful number is whether it is the *only* stop-free frame,
i.e. whether the read identifies its own frame (**measured**, 5 seeds):

| Read length | True frame unique | Mean stop-free frames |
|---|---|---|
| 75 nt | 0.244 | 2.17 |
| 107 nt | 0.528 | 1.61 |
| 150 nt | 0.767 | 1.26 |
| 250 nt | 0.943 | 1.06 |

At the Berlin read length, **47% of coding reads do not identify their frame**. Every
codon-level feature (A5, A9, A13, and the RAMICS-style codon HMM) must therefore either
marginalise over ~1.6 candidate frames or accept a wrong frame half the time. This is
also the sharpest argument for the merged-mate protocol change in §3.

**A5, A9 — codon usage and codon-pair bias. Data-starved at this length.** A 107 nt read
gives ~35 codons to estimate a 61-category distribution, over ~1.6 candidate frames
(A4). Codon-pair bias needs ~3,721 categories from ~34 observations. Both also measure
host adaptation, i.e. organism identity, not novelty.
*Killed by:* the arithmetic above. A novel virus in a known host has known-looking codon
usage, which is the point.

**A11, A12 — k-mer and peptide absence. The only single-read features with the right
polarity.** Instead of asking how close a read is to the database, ask how much of it is
*missing* from the database: the count of nucleotide k-mers (A11) or six-frame peptide
k-mers (A12) absent from an index of all known sequence. Novel sequence should be
enriched in absent k-mers. A12 is the better bet because peptide k-mers survive
synonymous divergence, which is exactly the regime where a 27–35% identical relative
defeats A1.

This is the interesting corner, and it is also where the honest difficulties are:
absence is produced by sequencing error, low-complexity sequence, and under-sequenced
*known* organisms just as readily as by novelty. The calibration problem is the
contribution, not the feature.
*Killed by:* absent-k-mer counts failing to separate target from background at a fixed
per-million budget, once error and low-complexity reads are excluded by A6/A14. Concrete
kill: if the background's absent-k-mer distribution has a tail reaching the target's
median, the feature cannot reach 1/M whatever the threshold.

**A13 — protein language models.** ESM-2 embeddings of ~35 aa fragments, as Rider does
for ORFs. Needs A4's frame (47% ambiguous at 107 nt), and is the most expensive option
per read by orders of magnitude. Worth testing only after A12, and mainly as the
comparison that says whether a learned representation beats k-mer absence.
*Killed by:* failing to beat A12 at equal per-million FP rate, which would make the
compute indefensible.

**A6, A14 — complexity and quality. Filters, not detectors.** Their role is to stop
artefacts entering the numerator. They have no novelty polarity and should never be
reported as detection features, but omitting them is how an absence-based score (A11/A12)
gets fooled: a quality-collapsed tail is maximally novel to any database.

### B. Pair-level — the under-explored axis

| # | Feature | Measures | Verdict |
|---|---|---|---|
| B1 | Mate conjunction of any per-read score | FP suppression | **open** — hinges on one unmeasured number |
| B2 | One-sided anchor (one mate known, one unplaceable) | Novelty adjacent to known | **demoted** — occupied since 2010, chimaera-dominated |
| B3 | Pooled composition across both mates | Doubles observation count | **open**, weak |
| B4 | Mate overlap / merging (short inserts) | Effective read length | **dead** — no headroom (§3.3) |
| B5 | Insert-size discordance | Structural variation | **dead** here (needs a reference) |
| B6 | Cross-mate frame continuity | — | **dead** (frames are unlinked without the insert) |
| B7 | Mate-pair graph connectivity among unplaceable reads | Grouping | **needs scale** |

See §3 for the depth treatment.

### C. Run-level aggregation

| # | Feature | Measures | Verdict |
|---|---|---|---|
| C1 | Overlap grouping / micro-assembly | Multi-read support | **needs scale** (~1 pair expected at 1e-6) |
| C2 | Statistical grouping, LSA eigengenome partitioning | Multi-read support without overlap | **open**, prior art — required baseline |
| C3 | Cross-sample recurrence (sites, dates) | Reality of a weak signal | **open** — best evidence-per-read in the catalogue |
| C4 | Replicate concordance | Artefact rejection | **open**, cheap |
| C5 | Abundance-ratio consistency across samples | Coherence of a candidate group | **needs scale** |
| C6 | PCR-duplicate family structure | Library complexity | **dead** as novelty evidence |

**C3 deserves emphasis.** A single read flagged once is uninterpretable at 1/M. The same
unplaceable sequence recurring in independent samples converts a weak per-read score
into strong evidence, because the background's false flags are drawn afresh each time
while a real organism persists. This multiplies across samples in a way no single-read
feature can, and it is the natural home for a calibrated score. It needs longitudinal
data, which Berlin (one run) does not provide — so this is the one direction that
requires a dataset the project does not yet have.

### D. Features that look available but are not

Listed because leaving them out would read as an oversight. Strand bias, coverage
evenness, and k-mer spectrum shape are all population-level statistics that need many
reads from the organism; at ~1 read they are undefined. Mapping-based features
(insert-size outliers, soft-clip profiles, split reads) need a reference to map to, which
by construction a novel family lacks — they are characterisation features, available
after a candidate exists, not detection features. Methylation and modified-base signals
do not exist in Illumina data.

---

## 3. Depth: the paired-end axis

Paired-end is the one axis where something is both available at sub-1x and not used for
novelty in the surveyed literature. Three mechanisms, in descending order of how much I
believe them.

### 3.1 The mate conjunction, and the number it hinges on

Requiring *both* mates to pass a threshold is the obvious way to buy false-positive
suppression. The flattering version — "the FP rate squares" — assumes the mates are
independent, and they are not: both come from the same molecule of the same organism.
Modelling each mate as passing with marginal probability `p` and correlation `rho`,
`P(both) = p² + rho·p·(1-p)`, and keeping the correlations for true and false flags
separate, at per-read sensitivity 0.3 and FP 1e-4 with `rho_true` = 0.9 (**model**):

| `rho_fp` | Pair sensitivity | Pair FP | FP reduction | Merit |
|---|---|---|---|---|
| 0.00 | 0.279 (0.93x) | 1.0e-8 | 10,000x | 9,300 |
| 0.10 | 0.279 (0.93x) | 1.0e-5 | 10x | 9.3 |
| 0.25 | 0.279 (0.93x) | 2.5e-5 | 4x | 3.7 |
| 0.50 | 0.279 (0.93x) | 5.0e-5 | 2x | 1.9 |
| 0.75 | 0.279 (0.93x) | 7.5e-5 | 1.3x | 1.2 |
| 0.90 | 0.279 (0.93x) | 9.0e-5 | 1.1x | 1.0 |

Merit is FP reduction per unit sensitivity surrendered; 1.0 means the conjunction is
pointless. The sensitivity cost is small and barely moves, because a novel organism's
two mates are both novel. **Everything rides on `rho_fp`, and the merit swings four
orders of magnitude between `rho_fp` = 0 and 0.1.**

That reduces the paired-end question to one falsifiable statement:

> The mate conjunction is worth having exactly when background false flags are
> **read-level accidents** (low complexity, adapter remnant, chimeric molecule, quality
> collapse) rather than **organism-level properties** (a genuinely divergent astrovirus
> in the sample). Accidents are roughly independent between mates; organisms are not.

My estimate is that `rho_fp` is high for homology-based scores — the Berlin background's
flags are overwhelmingly real divergent viruses, an organism-level property — and
materially lower for absence-based scores (A11/A12), where low-complexity and error-driven
absence is mate-specific. If that is right, the conjunction is nearly worthless for A1
and valuable for A12, which is a testable prediction and not a convenient one.

**This has been tried, informally, and the other way won.** The Nucleic Acid Observatory's
open lab notebook implements both rules on real wastewater across thresholds 15-45 and
**chooses the disjunction** - pass a pair if *either* mate scores at least 20 - with its
whole false-positive reduction trajectory computed under that rule. It plots the two mates'
scores against each other and never computes a correlation. Unpublished, and the most
likely thing to be raised as prior art against this section.

That sharpens the experiment rather than killing it. If the conjunction's value swings four
orders of magnitude on `rho_fp`, a group finding the disjunction better on their data is a
data point about their `rho_fp`, not a refutation. The claim that survives is narrower:
**the mate-score correlation has never been measured**, taxMaps shows mate concordance
tracks precision at rho = 0.994 without explaining why, and this model predicts both why
the disjunction won there and when a conjunction would win instead.

**The measurement is cheap and has not been done.** One translated-search run with mate
identity retained gives `rho_fp` directly; the existing read-level arm scored reads as an
unpaired pool and discarded it. Same for an absence-based score once A12 exists. This is
the single highest value-per-hour experiment in this document. See `PRIOR_ART.md` §1.3.

**Measure it against insert size, not as one number.** Wider templates make the mates more
independent, which lowers `rho_fp` (good: merit rises from 1.9 to 9.3 as `rho_fp` falls from
0.5 to 0.1) and lowers `rho_true` as well (bad: the pair-sensitivity factor falls from 0.93x
at `rho_true` 0.9 to 0.65x at 0.5). Both move together with insert width, so there is an
optimum, and a single `rho_fp` would hide it. See `HANDOFF.md` for the derivation and for
the coding-density ceiling that bounds the related wide-template argument.

### 3.2 The one-sided anchor

A pair where one mate places confidently against known sequence and the other is
unplaceable is a different kind of evidence from either mate alone: it localises novel
sequence *next to* something identifiable. It has the right polarity, it needs no
assembly, and it works at exactly one pair. It is how an integrated element, a
reassortant segment, or a divergent gene beside a conserved one would present.

**Demoted after the prior-art review.** The mechanism is not new: "one-end-anchored" (OEA)
pairs are a named primitive since NovelSeq (*Bioinformatics* 2010) and standard in insertion
discovery, and ViFi (*NAR* 2018) already does the divergence-tolerant version - keep a pair
if one end maps to human and the other maps to a virus **or scores well against viral-family
profile HMMs**. Applying it to free-living novel viruses in a metagenome is not published,
but it is a re-application rather than a new mechanism.

The confound is also worse than I assumed. A dengue/mosquito Nextera study found most, if
not all, detected host-virus chimaeras were artificial, with unrelated control RNA showing
about the same rate, via template switching during PCR. No shotgun-specific chimaera rate is
published; amplicon figures reach **>70% for less abundant species**, and our target is by
construction the rarest thing in the library. Worse, chimaeras form reproducibly across
independent amplifications, so replicate agreement will not filter them.
*Killed by:* the above, largely. If kept, the chimaera rate is the first experiment and not a
caveat to add later. See `PRIOR_ART.md` §3.

### 3.3 Merged mates: the protocol lever

This is the cheapest real win available and it is a library-design choice, not an
algorithm. §2 A4 measured frame identifiability at 0.528 for 107 nt and 0.767 for 150 nt;
§2 A3 measured the random-sequence stop-free rate falling from 0.72 to 0.12 (*E. coli*
shuffle) between 107 and 250 nt. Effective read length is the strongest lever on every
codon-level feature in the catalogue, and short-insert libraries supply it for free: if
the insert is shorter than twice the read length, the mates overlap and merge into one
contiguous read of up to ~200 nt.

SRR23973193 has 107 nt reads (**prior run** — the sweep set `--read-len 107` to match it).
Its insert distribution is **not** established: the sweep simulated spike fragments at
`-m 250 -s 40`, which was a parameter choice, not a measurement of the library, and I have
not checked it. If inserts really centre near 250 nt the mates do not overlap and this
signal is absent by construction, but that needs confirming before the paragraph below is
load-bearing; the mate-overlap rate is directly measurable with `bbmerge` on the real
reads, which is minutes of work. A library targeting ~180 nt inserts at 2x107 would merge most pairs and
roughly double effective read length, moving frame identifiability from ~0.53 toward
~0.85 and cutting random-sequence stop-free rates several-fold.

**Killed by the prior-art review — drop this mechanism.** Two findings, either sufficient.
The NAO's own read-length analysis argues on cost-per-informative-base grounds that a
positive gap beats any overlap and that more overlap is worse, since a 2x150 run on an 80 bp
fragment yields only 80 bp of unique sequence and short inserts read through into adapters.
And wastewater RNA at 2x150 **already mostly overlaps**: NAO reports viral read pairs not
collapsed at 42% (Yang 2020) and 24% (Spurbeck 2023), so roughly 58-76% merge already and
there is no design headroom to capture. The paragraph above assumed room that does not
exist. See `PRIOR_ART.md` §1.2.

This also supersedes the insert-size caveat I flagged: the mechanism fails whatever
SRR23973193's inserts turn out to be, so measuring them is no longer on the path.

### 3.4 What paired-end does not give

Frames cannot be linked across mates without knowing the insert length to the base, so
there is no "codon continuity across the pair" feature (B6). Insert-size discordance
needs a reference the novel family does not have (B5). And the pair does not double the
information in the way a 214 nt read would: two 107 nt chunks with an unknown gap are two
observations, not one longer one, for every frame-dependent feature.

---

## 4. How each candidate should be scored

One protocol, so results are comparable across features and against the measured floor:

1. **Truth labels are free.** ART names simulated spike reads after the reference
   accession; real background reads keep their SRA names (**prior run**, the existing
   read-level arm's trick). No manual curation needed.
2. **Report the full curve**, sensitivity against flagged background reads per million,
   not a single operating point. The shape is what decides usability, and the existing
   arm's curve is only interesting because it shows sensitivity hitting zero before the
   FP rate becomes acceptable.
3. **State the detectable fraction separately.** With the target's family removed from
   the database, part of its genome has no retained homolog at all; mixing that database
   property into a method's sensitivity makes a database limitation look like a method
   limitation.
4. **Run a spike-free arm of the same background.** The background genuinely contains
   divergent viruses, so "fraction of background flagged" is not a false-positive rate
   until the spike-free rate is subtracted.
5. **Five seeds before any number is quoted** (repo convention; the existing read-level
   arm has 2–3 and is under-powered, which is noted in its own report).
6. **Prefer continuous readouts.** Yesterday a binary readout (did any contig reach
   1,500 nt) produced a null that a continuous one overturned. Report scores, not flags,
   and threshold afterwards.

---

## 5. Shortlist

In order, cheapest discriminating experiment first:

1. **Measure `rho_fp` from the existing data** (§3.1). One DIAMOND run with mate identity
   retained, on reads already identified. Settles whether paired-end conjunction is worth
   anything at all, for a few hours of compute, and the answer is informative either way.
2. **Build A12, six-frame peptide k-mer absence**, and curve it against the 1,489/M floor.
   This is the only single-read feature with the right polarity and divergence tolerance.
   If it cannot beat translated search at a fixed per-million rate, the single-read
   direction is finished and that is a publishable negative. **Build it from existing
   parts**: the rarity index (*NAR Genom Bioinform* 2024) already publishes a graded
   per-k-mer absence score over 45,785 genomes and 21,871 proteomes, MetaProFi
   (*Bioinformatics* 2023) is an amino-acid k-mer Bloom index with a zero false-negative
   rate, and k_aa = 10 is the specificity floor (k = 7 is shared across phyla for >80% of
   k-mers). Quantify peptide-k-mer error amplification first - one substitution corrupts
   every peptide k-mer overlapping its codon, roughly k consecutive, and that is unmeasured
   in the literature.
3. **Read the Viromancer / Lyra-Virus OpenReview paper** (NeurIPS 2026 workshop, Sabeti lab)
   before claiming per-read novelty scoring is unoccupied. It scores individual reads under
   family-level holdout across six viral realms; whether it has pair aggregation or any
   calibration is unknown from the abstract.
4. **LSA as the grouping baseline** (C2) before anything is built on grouping, per the
   prior-art check.

Demoted: the one-sided anchor (§3.2), which needs a chimaera-rate measurement first.
Dropped: merged mates (§3.3).

Not recommended: anything codon-level at 107 nt (A5, A9, and the codon HMM) until
effective read length is addressed; A13 before A12.

## 6. What this document does not establish

- No feature here has been shown to work. Sections 2 and 3 bound ceilings and kill
  candidates; the only measured detection result in the project remains the
  translated-search negative.
- The ORF-likeness and frame measurements use whole genomes from the stress community as
  a stand-in for a wastewater background's composition. Real wastewater is a mixture with
  different GC structure, and the direction of the GC effect (§2 A3) says the mixture
  matters.
- `rho_fp` is a model parameter here, not a measurement. Every statement about the value
  of the mate conjunction is conditional on it, including my own estimate that it differs
  between homology and absence scores.
- The merged-mate mechanism (§3.3) is withdrawn: the prior-art review found both that
  wastewater RNA at 2x150 already mostly overlaps and that the one group to analyse insert
  length argued the opposite. SRR23973193's insert distribution remains unverified but no
  longer matters, since the mechanism fails either way.
- `PRIOR_ART.md` records what the review cost us: the per-million claim narrowed (PACIFIC
  already publishes the unit), the mate conjunction shown to have been tried informally and
  rejected in favour of the disjunction, and one-sided anchoring found to be a 2010
  primitive whose confound dominated real signal in the one controlled study.
- Cross-sample recurrence (C3) is the strongest aggregation idea in the catalogue and
  cannot be tested on the data the project holds.
