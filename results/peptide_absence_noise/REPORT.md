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
   headroom: 12 of 26 k-mers. That is the regime worth building for.

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
