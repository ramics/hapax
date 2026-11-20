# The noise floor under peptide k-mer absence

**Script:** `experiments/peptide_absence_noise.py` · **Seeds:** 5 · **Date:** 2026-11-20
**Data:** 6 genomes, 95 usable CDS (see `data/MANIFEST.md`)

`results/readlevel_features/PRIOR_ART.md` left one quantity unmeasured and described it as
"structurally worse than the nucleotide case": how far a single sequencing error propagates
through peptide k-mer space. It is the noise floor under the whole peptide-absence
proposal, because absence is the signal we intend to read as novelty and sequencing error
manufactures absence for free. This measures it, before any proteome index is downloaded,
because a failure here would kill the direction at zero data cost.

**It does not fail. The concern was wrong-signed, and a different constraint bites instead.**

---

## 1. Headline

1. **Peptide absence is more error-robust than nucleotide absence at matched k**, by roughly
   the synonymous fraction of the genetic code. The PRIOR_ART.md worry is retired.
2. **Indels are survivable at Illumina rates and fatal at nanopore rates**, which is a
   second and independent reason this feature belongs on short reads.
3. **A single 107 nt read cannot carry a one-in-a-million call on peptide absence alone.**
   At Q30 the one-per-million budget requires every k-mer in the read to be absent — zero
   headroom. Aggregation is forced, not optional.
4. At the thousand-per-million budget that translated search already set, there is real
   headroom: 12 of 26 k-mers.
5. **Berlin is measured, not assumed, onto the 0.001 row** (§8): SRR23973193's own quality
   strings give a mean error probability of 0.00089 on R1 and 0.00117 on R2, at a modal
   read length of 107. So point 3 applies to the real data rather than to a hypothetical.
6. **The mate pair relaxes the error term and does NOT rescue the target** (§9). Errors in
   two mates are independent, so that component squares: on error alone the budget moves
   from an impossible T ≥ 26 of 26 on a single read to T ≥ 16 of 26 on a pair. **But the
   error term is not what sets the floor.** PR #10 measured the total mate correlation of
   false flags on real Berlin pairs at ρ_fp = 0.77–0.98, with a 0.00 shuffled control, and
   found the both-mates rule at or below break-even. Measured against real Berlin
   background with a viral + Swiss-Prot index, **§9's own T ≥ 16 threshold is crossed by
   925,036 pairs per million against a target of 1** (§10). The 16-of-26 figure is an
   error-only bound and is not an operating point.
7. **So peptide 10-mer absence does not reach one per million on wastewater, and is not
   close** (§10). 86% of background reads have no known peptide 10-mer at all: absence
   there means "wastewater", not "novel virus". This is the ORF-likeness polarity failure
   one level up — the background has the property we were treating as the signal.

## 2. Synonymous shielding

A substitution that changes the nucleotide but not the amino acid costs peptide space
nothing. Measured on real CDS rather than taken from the codon table, because codon usage
is not uniform:

| | synonymous fraction |
|---|---|
| all substitutions | **22.7%** (22.6%–22.8% across seeds) |
| first codon position | 3.2% |
| second codon position | 0.2% |
| third codon position | 64.8% |

This is the mechanism behind §3. It is also the reason the intuition in PRIOR_ART.md went
the wrong way: that intuition counted the *reach* of a change (one changed amino acid
corrupts up to k_aa k-mers) and forgot the *rate* (about a fifth of substitutions never
produce a changed amino acid at all). Reach and k-mer count scale together, so they cancel,
and the rate term is what is left.

## 3. Corruption per read, peptide against nucleotide

Reads of 107 nt drawn inside a CDS, true frame known by construction. Peptide k_aa = 10
against nucleotide k_nt = 30, so the two cover the same span of sequence. Corruption is the
fraction of a read's k-mers that differ from the error-free truth, compared positionally.

| substitution rate | peptide corrupted | nucleotide corrupted | difference | reads with none corrupted |
|---|---|---|---|---|
| 0.0005 | 1.16% (1.12–1.21) | 1.49% (1.45–1.53) | -0.33 pp | 96.0% |
| 0.001 | 2.31% (2.27–2.37) | 2.99% (2.91–3.05) | -0.68 pp | 92.2% |
| 0.002 | 4.48% (4.34–4.56) | 5.78% (5.66–5.84) | -1.30 pp | 85.3% |
| 0.005 | 10.90% (10.76–11.07) | 13.94% (13.79–14.16) | -3.04 pp | 66.7% |
| 0.01 | 20.75% (20.62–20.89) | 26.04% (25.96–26.15) | -5.29 pp | 44.5% |

Peptide corruption is lower at every error rate tested, and the gap widens with error rate.
The ramped error profile, with errors concentrated toward the 3' end as real Illumina data
has them, is reported beside the flat one in the JSON; it moves the peptide number by
little, because clustered errors corrupt fewer *distinct* k-mers than dispersed ones. The
flat model is therefore the pessimistic one and the one quoted here.

## 4. Indels

A substitution is local. An indel frameshifts everything downstream and destroys every
peptide k-mer after it, while leaving nucleotide k-mers only locally disturbed. Substitution
rate held at 0.001:

| indel rate | peptide corrupted | nucleotide corrupted |
|---|---|---|
| 0.0 | 2.35% (2.29–2.46) | 3.00% (2.93–3.09) |
| 1e-05 | 2.32% (2.24–2.38) | 2.99% (2.90–3.07) |
| 0.0001 | 2.95% (2.88–3.02) | 3.66% (3.57–3.75) |
| 0.001 | 8.70% (8.61–8.86) | 9.40% (9.32–9.55) |
| 0.01 | 46.96% (45.87–47.61) | 47.70% (46.68–48.33) |
| 0.05 | 91.03% (90.89–91.22) | 91.50% (91.38–91.68) |

At Illumina indel rates the term is negligible. At nanopore rates it dominates and the
feature stops working. This is independent of the cost argument recorded in
`hapax-readlevel-detection-pivot` and points the same way.

## 5. Headroom: the number that shapes the design

A read from a known organism is falsely flagged as novel if at least T of its peptide
k-mers come back absent. Under the idealisation below, the absent k-mers are exactly the
corrupted ones, so the error-induced false-flag rate at threshold T is P(corrupted ≥ T).
Headroom is how many k-mers remain above T for a genuinely novel read to spend.

107 nt at frame offset 1 or 2 carries 26 peptide 10-mers, and the worst case is quoted.

| substitution rate | ≤93,478/M | ≤1,489/M | ≤100/M | ≤1/M |
|---|---|---|---|---|
| 0.0005 | T ≥ 1, headroom 25 | T ≥ 11, headroom 15 | T ≥ 19, headroom 7 | T ≥ 24, headroom 2 |
| 0.001 | T ≥ 1, headroom 25 | T ≥ 14, headroom 12 | T ≥ 21, headroom 5 | T ≥ 26, headroom 0 |
| 0.002 | T ≥ 9, headroom 17 | T ≥ 19, headroom 7 | T ≥ 22, headroom 4 | not reachable \* |
| 0.005 | T ≥ 11, headroom 15 | T ≥ 22, headroom 4 | T ≥ 26, headroom 0 | not reachable \* |
| 0.01 | T ≥ 15, headroom 11 | T ≥ 26, headroom 0 | not reachable \* | not reachable \* |

Worst threshold and smallest headroom across 5 seeds, 10,000,000 samples per seed (resolution 0.10 per million). An entry marked \* hit a zero count in at least one seed, so it is a bound rather than an estimate.

**Reading this.** At the 1,489 per million that translated search set as its cleanest
operating point (`results/abundance_sweep/REPORT.md`), Q30 leaves 12 of 26 k-mers of
headroom — a real working margin. At 100 per million the margin is 5. At one per million
the threshold reaches the full k-mer count: the entire read must be absent, and a single
synonymous-free error anywhere in it produces a false flag. At error 0.002 the one-per-
million budget is not reachable at all.

So the design consequence is not subtle. **At sub-1x coverage with a one-in-a-million prior,
a single read's peptide absence cannot clear the budget on its own.** Whatever the index
says about how absent a novel read really is, the error floor alone forecloses it. Pairing
and grouping are therefore load-bearing rather than refinements, which is consistent with
the pooling arithmetic in `results/readlevel_features/REPORT.md` and makes the sibling
thread's mate-correlation measurement a dependency of this one rather than a parallel track.

## 6. A correction made during this run

The first pass computed §5 from the 20,000-read direct simulation. With 20,000 reads the
smallest non-zero rate observable is 50 per million, so every threshold reported as meeting
the 100- and 1-per-million budgets was a zero count rather than a measurement. It gave
headroom 5 at one per million where the honest answer is 0.

The published numbers come instead from a sampler that works through the structure of the
problem: given which amino-acid positions changed, the corrupted k-mer count is determined,
so only the substitution count (Poisson) and the measured non-synonymous probability are
needed. That runs 10 million samples per condition, resolving 0.1 per million, and it
reproduces the direct simulation's mean corrupted fraction to within 0.003 at every error
rate — the validation is in the JSON under `validation_against_direct_simulation`. Each
threshold in §5 carries the false-positive rate actually achieved and a `resolved` flag,
which is false where the count hit zero and the number is a bound rather than an estimate.

This is the fourth time on this project that a claim has turned on a zero count or an
unbalanced arm. The check that caught it was asking what rate the sample size could
resolve, before reading the table.

## 7. What this does not establish

**The false-positive side only, and optimistically.** The idealisation is that a known
organism's error-free peptide k-mers are all present in the database, so absence comes only
from error. Real absence also comes from undersequenced known strains, low-complexity
peptides, reads from non-coding sequence, and the five wrong frames. Every rate here is
therefore a floor and the real feature can only be noisier, which makes a pass here
permissive rather than encouraging.

**Nothing about sensitivity.** How absent a genuinely novel read is depends on the distance
to its nearest retained relative and needs the index. That is step two.

**The true frame is known by construction.** 47% of coding reads at 107 nt cannot identify
their own frame (`results/readlevel_features/REPORT.md`), so a frame-blind scorer does
worse than every number here. Quantifying that gap needs the index too.

**The error rate is swept, not measured.** SRR23973193's actual quality profile is not on
disk in this container, so which row of §5 applies to the Berlin data is unknown. Measuring
it is cheap once the run is fetched and should happen before step two is costed.

**The genome set is narrow.** 95 CDS from 6 genomes, of which lambda contributes 72, so the
codon-usage breadth behind §2 is thin. The measured 22.7% sits close to the uniform-codon
expectation, so it is unlikely to move much, but it has not been shown across a wide range
of codon usage and should not be quoted as if it had.

## 8. Which row Berlin is actually on

§7 recorded that the error rate was swept rather than measured, and that which row applies
to SRR23973193 was unknown. `experiments/berlin_quality_profile.py` closes that, reading
the run's own quality strings from ENA over HTTP (no SRA toolkit needed), 400,000 reads per
mate:

| | reads | modal length | mean error probability | first position → last |
|---|---|---|---|---|
| R1 | 400,000 | 107 | **0.00089** | 0.00175 → 0.00092 |
| R2 | 400,000 | 107 | **0.00117** | 0.00120 → 0.00136 |

**Berlin sits on the 0.001 row**, so §5's conclusion applies to the real data: at a
one-per-million budget a single read must have every one of its 26 peptide 10-mers absent.

Two things worth noting in the shape. The profile is close to flat rather than ramped —
R1's error probability actually *falls* across the read — which is what NovaSeq's binned
quality scores produce, and it means the flat model used in §3 is the right one here rather
than merely the pessimistic one. And these are Phred scores, the instrument's estimate of
its own error rate, which is usually optimistic; the true rate is at or above this, and the
headroom table gets worse as the rate rises.

## 9. What a mate pair buys

A sequencing error in one mate and a sequencing error in the other are separate events, so
for the error component of a false flag the mate correlation is zero by construction. This
is the one component of the absence score whose correlation does not have to be measured.
Requiring both mates to clear the threshold therefore squares the error-induced rate.

It also happens to be the only aggregation guaranteed to exist at the target abundance. At
one in a million a rare organism contributes about one fragment, so there is nothing to
group read-to-read (`results/readlevel_features/REPORT.md`); the two mates of a single
fragment are always both there.

Threshold needed to hold error-induced false flags at one per million, 26 k-mers per read:

| substitution rate | single read | mate pair | per-read rate the pair needs |
|---|---|---|---|
| 0.0005 | T ≥ 24, headroom 2 | **T ≥ 11, headroom 15** | 645/M |
| 0.001 | T ≥ 26, headroom 0 | **T ≥ 16, headroom 10** | 882/M |
| 0.002 | not reachable | **T ≥ 20, headroom 6** | 721/M |
| 0.005 | not reachable | **T ≥ 23, headroom 3** | 797/M |
| 0.01 | not reachable | **T ≥ 26, headroom 0** | 818/M |

At Berlin's measured error rate the budget moves from **impossible on a single read** to
T ≥ 16 of 26 on a pair, leaving 10 k-mers of headroom, with each mate individually allowed
a rate near 900 per million.

**Read that as an error-only bound, not an operating point.** See the reconciliation
immediately below: the organism term, which this calculation excludes, has since been
measured and it dominates.

**What this does not cover, and it is the whole of step two.** Squaring applies to absence
caused by error. Absence caused by the organism — an undersequenced known strain, a
low-complexity protein, a fragment of non-coding sequence — is a property of the molecule
and appears in both mates together, so it does not square at all. That correlated term is
unmeasured, it is what will actually set the floor, and it can only be seen on real
background with a real index. This result therefore says the error floor is not the binding
constraint once mates are used; it does not say the budget is met.

It also sharpens the conjecture in the sibling thread's `rho_fp` work from the other side:
that analysis predicted the mate conjunction would be near-worthless for homology scores
and valuable for absence scores, because error-driven absence is mate-specific. The error
component of rho_fp for peptide absence is not merely low, it is zero. What remains to
measure there is the organism-driven component, which is the same quantity step two needs.


## 10. Reconciliation with PR #10, which measured the total correlation

§9 was written before `results/mate_width/REPORT.md` (PR #10, the mate-separation thread)
landed its ρ_fp measurement. That measurement supersedes the practical conclusion here, and
this section is the agreed joint reading rather than a competing one.

**What they measured.** Six-frame peptide 10-mer absence, best-frame statistic, scored on
real SRR23973193 pairs with mate identity kept. Every read in that run is background, so
every flag is a false flag and the correlation between real mates is ρ_fp directly. They
get **0.85–0.98 across insert bins**, against **0.00** from a shuffled-pairing control. Aware
that 88.3% of background reads saturate at a score of 1.0 against a viral-only index, they
repeated it on the pairs where the score has genuine range: **0.770 ± 0.009** with at least
one mate non-saturated, 0.773 ± 0.006 with both, each against a shuffled control at
−0.004 ± 0.012. So the correlation is a property of the flags, not of the ceiling. Their
pair-rule arithmetic puts the both-mates merit at 0.80–1.02, at or below break-even, and
their headline is that at one false flag per million no pair rule gets there at any insert.

**Why both results are right.** They decompose the same quantity. A false flag's absence
has an instrument component and an organism component. The instrument component is
independent between mates — that is §9, and it is exact. The organism component is a
property of the molecule: both mates come off the same fragment of the same
under-sequenced organism, so both are equally unknown to the index, and it does not square
at all. Their measurement is of the **total**, which is what merit depends on, and it says
the organism component dominates. §9 therefore bounds a term that turns out to be the
smaller one.

**The correction to this report.** The 16-of-26 pair threshold is an error-only bound. It
is not an achievable operating point and should not be quoted as one. The honest statement
is: sequencing error stops being the binding constraint once mates are used, and the
binding constraint is organism-level absence, which the pair rule does not help with.

**The Swiss-Prot arm landed, and it settles it.** PR #10 has since run the same measurement
against a viral + Swiss-Prot index of 221M peptide 10-mers, 2.5x the viral-only one. I
re-derived the numbers from their per-pair scores with
`experiments/absence_threshold_check.py` (5 seeds, 40,000 real Berlin pairs per seed, every
read background, so every flag false):

| | measured |
|---|---|
| background reads fully absent, single read | **86.35%** (88.3% against viral-only) |
| background pairs with both mates fully absent | 81.33% |
| mate correlation, unconditional, no subsetting | 0.860 |

And applying §9's own threshold to that background:

| threshold | single read | both mates | if independent | conjunction gain |
|---|---|---|---|---|
| T ≥ 16 of 26 (§9's pair threshold) | 938,046/M | **925,036/M** | 879,930/M | 1.05x |
| T ≥ 21 of 26 | 919,490/M | 891,205/M | 845,462/M | 1.05x |
| T ≥ 26 of 26 (every k-mer absent) | 862,891/M | 813,300/M | 744,581/M | 1.09x |

**The §9 threshold is exceeded by real background in 925,036 pairs per million against a
target of 1.** That is not erosion at the margin, it is six orders of magnitude, and the
conjunction buys a factor of 1.05 rather than the squaring the error term promised. The
honest verdict is that **peptide 10-mer absence as a thresholded detector on wastewater
background does not reach one per million, and is not close.**

**The mechanism, and why it is familiar.** 86% of Berlin background reads have no known
peptide 10-mer at all. Absence does not mean "novel virus" there, it means "wastewater":
the background is itself largely unknown to protein databases. That is the same polarity
failure that killed ORF-likeness in `results/readlevel_features/REPORT.md` — the background
has the property we were treating as the signal — arriving one level up. It also predicts,
from our side, the transfer problem recorded for Lyra-Virus in
`results/readlevel_features/VIROMANCER_REVIEW.md`: a clinical respiratory background is
human plus known pathogens, where "unknown" is informative; an environmental virome is not.

**What remains open is how much of that 86% is index incompleteness.** Expanding the index
2.5x moved absence by 2 points and the tightest-bin correlation from 0.975 to 0.964, so the
returns look flat, but two points do not make a curve. A complete index is not feasible in
this sandbox — UniRef50 is ~2e10 peptide 10-mers, roughly 160 GB of sorted hashes against
54 GB free, and a usable Bloom filter needs more RAM than the box has. PR #10 is therefore
building an index-size scaling curve with deterministic UniRef50 subsamples, and the
extrapolation decides whether any achievable index changes this. My own expectation, stated
in advance so it can be scored: closing six orders of magnitude needs essentially every
background read to match, and nothing in the first two points suggests that, so I expect
the curve to confirm the kill rather than overturn it.

**Previously open, now closed.** The index used was viral-only,
and no index available in this sandbox covers environmental protein space — which is why
88.3% of background saturates. Two things follow. The 0.770 subset is conditioned on the
score being non-saturated, and conditioning on an outcome can bias a correlation in either
direction, so it is a different estimand from the unconditional one rather than a cleaner
version of it. And the saturation itself is a property of index coverage, not of the
feature: against a comprehensive index, background reads from known organisms should match,
and both the absence rate and its mate correlation could move substantially. That arm has now landed and is reported above.

**One more of their results lands on this feature directly.** They measured six-frame
peptide 10-mer absence under a family holdout and found it **saturates**: 94.6% of windows
fully absent at k = 10. Novel reads being maximally absent is good for detection and bad
for grading — a plain absence count cannot rank degrees of novelty at the divergence that
matters. Any design built on this feature has to be a detector with a calibrated threshold,
not a novelty score, or has to find a statistic with range where the count has none.
