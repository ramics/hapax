# Cross-sample recurrence of the residue

**What was measured.** For every read pair in a depth-stratified query set drawn from
SRR23973193, the fraction of each mate's canonical 31-mers present in each of the **115 other
paired RNA libraries** of PRJNA948850 — the same Berlin plant, **84 distinct waters over 40
dates spanning 502 days**. Presence is exact k-mer witnessing, which costs one scan per
subject library whatever the query count, so all five seeds travel in one query file
(198,284 pairs, 18.45 M distinct 31-mers) and the seed rides in each read's name.

Everything here is **on SRR23973193 as the query**. The drift caveat applies to the subjects
too: this plant's rRNA fraction went 38.0% → 67.1% over eight months, so the panel is not a
stationary background.

**Why this exists.** The residue census put the unplaceable pile at **median within-run depth
1** (`results/residue_census/REPORT.md`), which rules out by arithmetic every method needing
a second observation. The 115 other libraries are the only place a second observation can
come from, and they convert one observation into an occupancy pattern. It replaces the old
71.2% cross-library figure (`results/dna_rna_pair/cross_witness_s1.json`), which pooled rRNA
and had no depth axis, so it was an upper bound rather than a measurement.

Scripts: `experiments/panel_inventory.py`, `recurrence_queries.py`, `recurrence_pass.py`
(+ `ec2_recurrence.py`, `recurrence_userdata.sh`), `recurrence_report.py`,
`recurrence_self_control.py`, `recurrence_grouping.py`. Run: `i-02e76104598299792`,
c6i.2xlarge on demand, 3.5 h, **$1.78**; pilot over 8 libraries `i-0a1cd44806a8b4c31`, $1.07.

## 1. Operating points, because recurrence has no natural one

| name | presence of a mate in a subject |
|---|---|
| `any` | at least one of its 31-mers present — the most sensitive reading the instrument has |
| `most` | at least half of its 31-mers present — close to "the molecule is there" |

Three denominators, reported throughout because they differ: **115 runs**, **84 waters**
(date + 2 h composite window; preparation variants of one water are not independent),
**40 dates**. No threshold is applied to any read: occupancy is reported as a distribution.

## 2. Controls, read before the test

| control | what it bounds | reading |
|---|---|---|
| own library (`self_control.json`) | the sensitive end | every real pile **1.0000** (min 1.000, fully witnessed) |
| composition-preserving shuffle | the specific end | every shuffled pile **0.0000** at both points, in 0 of 115 |
| rRNA | the familiarity reference | depth ≥10: on **all 40 dates** for 97.8% of reads |
| subject sequencing depth | occupancy as a panel property | Spearman ρ = **−0.02** against subject read count |
| identified reads at matched own depth | the panel's detection limit | depth 0: on a median **37 of 40 dates** |

The last row is the one that matters. Residue reads are rare *here*, and rarity here is the
confounder, so the comparison is against reads of the **same within-run depth** that a
database could name. Per-subject witnessing spans 0.157–0.493 (median 0.379) of the query
set, so no subject is inert and none is doing all the work.

**Stratified, not matched.** Recurrence and within-run depth are one variable (unplaced
median depth 1, identified 87), so the query set takes a fixed number of pairs per
(pile, depth bin) cell and recurrence is a curve against depth with the identified and rRNA
curves beside it. A single matched number is a weighted average away; the curve additionally
answers whether the piles differ **at the same depth**, which is the question.

## 3. Occupancy over the 40 dates, mate 1, `any`

Fraction of reads present on **no** other date, with the five-seed spread, and the
occupancy quantiles.

| pile | own depth | n | in none | spread | p10 | p50 | p90 | in all 40 |
|---|---|---|---|---|---|---|---|---|
| unplaced | 0 | 9,873 | **21.5%** | 21.0–22.2 | 0 | 6 | 30 | 1.6% |
| unplaced | 1–2 | 10,000 | **19.5%** | 18.6–20.4 | 0 | 7 | 32 | 2.6% |
| unplaced | 3–9 | 10,000 | 14.4% | 13.5–14.9 | 0 | 14 | 39 | 8.0% |
| unplaced | 10–99 | 2,036 | 1.4% | 0.5–2.9 | 15 | 40 | 40 | 57.7% |
| unident. contig | 0 | 4,884 | 1.6% | 0.9–2.0 | 7 | 29 | 40 | 10.7% |
| unident. contig | 1–2 | 9,615 | 1.8% | 1.4–2.4 | 8 | 30 | 40 | 12.9% |
| identified | 0 | 5,577 | **0.9%** | 0.5–1.1 | 17 | 37 | 40 | 26.3% |
| identified | 1–2 | 10,000 | 0.5% | 0.2–0.7 | 19 | 38 | 40 | 30.2% |
| identified | 10–99 | 10,000 | 0.0% | 0.0–0.1 | 22 | 40 | 40 | 56.0% |
| rRNA | 0 | 1,385 | 0.6% | 0.0–1.1 | 33 | 40 | 40 | 80.5% |
| rRNA | 10–99 | 10,000 | 0.0% | 0.0–0.0 | 40 | 40 | 40 | 97.8% |
| shuffled-unplaced | any | 9,909 | 100.0% | 100–100 | 0 | 0 | 0 | 0.0% |

**At matched within-run depth the unplaced residue is absent from the whole plant 24× as
often as identified reads** (21.5% vs 0.9% at depth 0). Both piles are equally rare in our
library, so rarity here cannot explain it. The gap narrows with depth and is nearly gone by
depth 10–99 (1.4% vs 0.0%), so **depth 1 is two populations**: thin material the plant does
hold, and molecules it does not hold at all.

Mate 2 reads the same direction and slightly stronger (23.7% vs 1.4% at depth 0), as its 93
nt against mate 1's 107 predicts — 62 scorable 31-mers against 76, so fewer chances to be
witnessed. Pooling the mates would average a real asymmetry; both are tabulated in
`recurrence_full_s1.json`.

At `most`, every pile moves and the ordering survives: unplaced depth 0 **43.2%** in none
against identified **17.6%**, a 2.5× ratio rather than 24×. The sensitive point separates the
piles; the strict point separates rare material from common material in both piles at once.
**Neither is the right one, which is why both ship.**

Over 84 waters instead of 40 dates the `in none` column is unchanged by construction
(absence is absence) while `in all` collapses — identified depth 0 goes 26.3% → 5.9%, rRNA
80.5% → 69.5% — so an occupancy *fraction* must name its denominator or it is unreadable.

## 4. Error or missing molecule

A read sharing nothing with 115 libraries could be common material carrying errors. Errors
strike the two mates independently, so error predicts **exactly one** mate unwitnessed;
a molecule the panel does not hold takes **both**. Same instrument as the mate-disagreement
error floor in `results/residue_structure/`.

| pile, own depth 0 | both unwitnessed | exactly one | ratio |
|---|---|---|---|
| rRNA | 0.00% | 1.59% | **0.00** |
| identified | 0.43% | 1.47% | 0.29 |
| unident. contig | 1.13% | 1.80% | 0.62 |
| unplaced | **17.91%** | 8.82% | **2.03** |
| shuffled | ~100% | ~0.1% | ~1,000 |

The two ends are both present and both behave: on rRNA, where a fully unwitnessed read has
to be an artefact, the field reads **0.00**; on shuffled sequence, which no library holds,
it saturates. The residue sits an order of magnitude from the error end. **The residue has
more missing molecules, not more errors** — the same conclusion read-level mate disagreement
reached by a different route (2.33% known vs 2.46% residue errors, 3.0% vs 7.0% absence).

## 5. Transient or resident

At **matched date-occupancy**, so the comparison is between reads the panel holds equally
often, the residue's shape in time is indistinguishable from known material.

| pile | dates occupied | n | span p10 | p50 | p90 | median months | within 14 days |
|---|---|---|---|---|---|---|---|
| unplaced | 2–3 | 2,839 | 29 | **191** | 384 | 2 | 3.8% |
| identified | 2–3 | 243 | 28 | **223** | 396 | 2 | 4.1% |
| unident. contig | 2–3 | 592 | 21 | 184 | 346 | 2 | 7.3% |
| unplaced | 4–8 | 4,235 | 218 | **359** | 460 | 5 | 0.0% |
| identified | 4–8 | 807 | 217 | **339** | 461 | 5 | 0.0% |
| unplaced | 9–20 | 6,408 | 372 | 448 | 502 | 9 | 0.0% |
| identified | 9–20 | 3,391 | 339 | 460 | 502 | 10 | 0.0% |

A burst — one shedding event, a few days — would show a small span inside one month. **There
is no burst signature anywhere in the residue.** What recurs at all recurs like resident
material: months apart, across seasons.

So the two legible populations in the residue are *resident at low level* and *recurring
across seasons*, and both look like known material once occupancy is matched. **What is
distinct is the 21.5% that is nowhere else in 502 days of the same sewer.**

**The limit this cannot pass.** Absence from 115 libraries means below the panel's limit in
all of them, not non-existent, and a molecule below one copy per library is invisible by
arithmetic. The limit is characterised rather than assumed — identified reads of the same
within-run rarity are found on a median 37 of 40 dates — but a single shedding event and
something permanently sub-single-copy produce the same reading here, and this instrument
cannot separate them.

## 6. Grouping by where, not by what

Reads at depth 1 cannot be grouped by sequence: no overlap, no neighbourhood, no second copy
(`results/mate_friends/`, `results/residue_structure/` §4). Occupancy over 84 waters is a
different space, and it needs no second copy in this run. The label-free ground truth is
free here: a pair's two mates came off one molecule, so rank every mate 2 by occupancy
similarity to a given mate 1 and ask whether the own mate is in the top 10 — the same
instrument `results/mate_friends/` used on sequence neighbourhoods, so the numbers compare.

Chance is `k/n`. The null is a **curveball swap** preserving both marginals (each read's
sample count and each sample's read count), because an abundant read is present nearly
everywhere and a deeply sequenced subject witnesses nearly everything, and a null that
destroyed either would make those look like structure.

**The kill criterion was stated first**: vectors with no entropy — in every sample or in
one — cannot group anything. 79–81% of unplaced vectors are informative (neither empty nor
full) over 84 dimensions, so the question is live.

### 6a. The control that decides it, read first

The median insert is 194 nt against 107 + 93 = 200 nt of read, so **many pairs overlap**:
their mates are partly the same bases, witnessed by the same k-mers, with identical
occupancy vectors for a trivial reason. An own-mate rate over those pairs measures overlap,
not co-occurrence — the same shape as the chimaera test whose query sat inside its own
reference. So the headline rate is the one over pairs whose mates share **no** 31-mer, with
candidates restricted to the same subset so chance moves with it. Disjoint mates are 62–64%
of unplaced pairs, 72–75% of unidentified-contig pairs, 76–77% of identified pairs. The
residue appears to overlap more, and that is itself a selection effect rather than a fact
about its inserts: the control runs on pairs with **both** mates witnessed somewhere, and an
overlapping pair shares the k-mers that do the witnessing, so it clears that bar more often
in the pile where witnessing is scarce. It is the direction that inflates the uncontrolled
number most, which is why the controlled one is the headline.

### 6b. Own mate found by occupancy alone

4,000 reads per pile per seed, k = 10, 84 water dimensions, five seeds.

| pile | own mate in top 10 | chance | curveball null | × chance | × null |
|---|---|---|---|---|---|
| unplaced, disjoint mates | **0.302–0.341** | 0.0051 | 0.012–0.022 | ~61× | ~19× |
| unident. contig, disjoint | 0.220–0.249 | 0.0034 | 0.010–0.014 | ~68× | ~19× |
| identified, disjoint | 0.192–0.205 | 0.0033 | 0.008–0.017 | ~60× | ~16× |
| unplaced, all pairs | 0.451–0.480 | 0.0032 | 0.016–0.019 | ~145× | ~27× |

**Occupancy over 84 waters identifies a read's own molecule, with no sequence comparison at
all.** A third of residue pairs whose mates share not one 31-mer still find their own partner
in a top-10 list of 4,000 candidates. The multiplier over the marginal-preserving null is
~19×, so it is not an artefact of abundant reads being everywhere or deep subjects witnessing
everything.

**Reproduced from clean.** The overlap control landed after the first grouping run, so the
run was repeated from a committed tree (`c4beb51`, `tree_dirty: false`, same
`script_sha256`) and every rate above came back byte-identical. That is the re-run the
provenance rule forces, and here it is also an actual reproduction.

**Two things this does not say.** A multiplier is not transportable across candidate-pool
size — the absolute rate falls roughly as 1/n while the multiplier rises, measured in
`results/mate_friends/` — so these numbers belong to n = 4,000 and nothing else. And the
rate for the residue is *higher* than for identified reads, which is a selection lean, not a
discovery: all-empty and all-full vectors are excluded as uninformative, and what remains of
the residue is the sparse middle, which is the most distinctive kind of vector there is.

### 6c. Structure, and the absence of clusters

Mean Jaccard to the ten nearest neighbours, informative vectors only, against the curveball
null on the same reads:

| pile | kNN Jaccard | null | excess |
|---|---|---|---|
| unplaced | 0.587–0.601 | 0.543–0.564 | **~0.04** |
| unident. contig | 0.783–0.791 | 0.742–0.750 | ~0.04 |
| identified | 0.840–0.849 | 0.806–0.815 | ~0.035 |

The excess over the null is real, small, and **the same size in the residue as in named
material**. So occupancy space carries structure beyond the marginals, and the residue is no
more clustered in it than the known material is — which is what the joint-view work already
found on ten other axes (PRs 17–19: one continuum, no clusters). **Report it as a density,
not a clustering**; the own-mate control is what licenses using occupancy to link reads
pairwise, and nothing here licenses cutting the pile into groups.

## 7. Polarity, and what this is not

Co-occurrence similarity is a **familiarity** field, like every other field in this repo:
reads that recur together recur, and recurrence is the opposite of novel. Nothing here
scores novelty, and the searching that produced the piles only ever removed background.

What recurrence buys is the other thing: **it is the one axis that works past the depth-1
wall**, because it needs no second copy in this library. A read seen once here can still be
seen on 40 dates — and 21.5% of the unplaced pile is seen on none of them, which is a
quantity where there was previously an adjective.
