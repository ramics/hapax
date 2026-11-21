# Profile reach is a negative; the pair as one molecule is a positive

Two questions about what can be said about a molecule seen exactly once, measured on the
same divergence ladder in the same pass, 5 seeds. Both were cheap and local, and both
were chosen because SPEC.md already leans on them without having measured them.

The tables are [`TABLES.md`](TABLES.md), generated from `cells.json` and
`paired_merged.json` by `experiments/report_profile_reach.py`.

## Why these two, and for which pile

The ~11 points of this library that nothing can place sit at median 31-mer depth 1
(`results/residue_census/` §3): 91% of that pile occurs fewer than five times in 8.7 M
pairs and about a third occurs nowhere else. Everything that needs a second observation —
assembly, neighbourhood grouping, mate-linkage kinship, within-run recurrence — has
nothing to work with there, and that is sampling rather than tuning. So for this pile the
only available deliverable is a *description of one molecule*, which is what SPEC.md §4.1
says. The two instruments that description leans on hardest are the ones measured here.

## The design, and the control that decides each claim

The holdout is D00239, human rhinovirus A1B — the only clean holdout of
`results/spike_holdout/`, absent from this sample's own assembly, unlike the other two
genomes spiked there. It is walked away from itself by uniform substitution at
0/5/10/20/30%, exactly as that ladder does, and sequenced through this library's own
measured parameters: mates of 107 and 93 nt after the 14 nt non-template prefix trim
(`results/nontemplate_prefix/`), inserts from a lognormal fitted to the measured median
194 and p95 363 (`results/insert_profile/`), and substitution errors at the
reference-free rates 0.00068 along the body and 0.0019 at the final base
(`results/quality_profile/`, which supersedes the Q-derived rates in
`results/berlin_quality/` as upper bounds conservative by about 2.6x).

Three piles are searched in the same pass, and the order they are read in is the point:

1. **composition-matched shuffles** of the same reads — what each instrument reads on
   sequence carrying no information;
2. **real library reads** — known material, which the cascade's indexes do identify, and
   which gives every field its direction;
3. the **holdout** at each rung.

A field that reads the same on (2) as on (3) describes the database rather than the
query. That has now happened to six single-read axes in this project, and in all three
failures of 2026-11-20 the tell was the control moving with the test.

The shuffle is searched **at every length the comparison uses** — 106 nt, 92 nt and the
149 nt merged length — rather than only at the shortest. That exists because of §2's
claim: a longer query can reach a fixed bitscore floor by chance more easily than a short
one, so a merged-against-mate gain measured at a fixed floor could be arithmetic rather
than recognition. The control answers it cleanly. Shuffled reads score **0.00% at
DIAMOND bitscore >= 40 and >= 50, and 0.00% at Pfam E <= 1e-2, at all three lengths**.
Chance is not what is being measured.

One by-product worth keeping, because it decides how results from different runs can be
compared. Pfam's gathering thresholds gave identical figures across two passes of this
experiment, while the per-domain E-value cells moved by up to 0.42 points when the target
set grew from 711,101 to 775,559 peptides — an E-value depends on the size of the set it
was computed against and `GA` does not. **So `GA` is the operating point that is stable
between runs and an E-value cut is not**, which matters for a tool whose §6 requirement is
per-sample recalibration.

## 1. Profile reach: a negative, and a correction to SPEC.md

SPEC.md §4.2 lists "profile reach" as one of the five axes of the signal-versus-noise
measurement, and §7 states it as a blind spot: *"Pfam-A profiles are the only axis still
alive past 20% divergence, where DIAMOND and every nucleotide index read zero."* That
sentence was never measured. Both halves of it are false.

**Pfam-A is weaker than translated sequence search at every rung of the ladder.** Mate 1,
5 seeds, DIAMOND at bitscore >= 50 against Pfam at E <= 1e-2: 91.5% against 39.9% at 0%
divergence, 81.8% against 27.6% at 5%, 61.3% against 18.6% at 10%, 10.0% against 5.8% at
20%, 0.1% against 0.5% at 30%. Only at the last rung do they meet, and there both are
within noise of nothing. This is not a threshold artefact: the profile side is reported as
a curve over per-domain E-value with the curated gathering thresholds as one point on it,
and the ordering does not change anywhere on that curve (`TABLES.md` §1).

**DIAMOND does not read zero past 20% divergence.** At a bitscore floor of 40 it
recognises 39.9% of reads at 20% divergence and 5.7% at 30%, where Pfam reads 5.8% and
0.5%. The claim that it reads zero came from the removal table in
`results/spike_holdout/`, which is a *cascade* figure at the cascade's operating point —
chosen for subtraction, where a false removal is expensive — and not a measurement of how
far a search can see. Those are different questions and the spec collapsed them.

**The gap carries nothing.** The reportable quantity was never the profile hit rate — a
domain hit means "carries a recognisable domain", and this background is itself mostly
coding (`results/readlevel_features/`), so a hit rate is a familiarity signal. The
quantity is the gap: no significant translated hit, yes significant profile hit. It peaks
at **3.8% (2.2–4.8) on the holdout at 20% divergence against 2.6% (2.2–3.1) on reads the
cascade already identifies** — overlapping seed ranges, a difference of about one point,
and the control moving with the test for the sixth time in this project (`TABLES.md` §2).

**What survives.** Profile search is not useless — it is a legitimate `motif` field in the
per-read report, and the shuffle null confirms it is measuring structure rather than
composition. What does not survive is the claim that it reaches further than anything
else, and SPEC.md §4.2 and §7 both have to be corrected, because the axis was listed
there as the one thing alive in the band this project cares about.

## 2. The pair as one molecule: a positive, and a wall that moves

A pair is one molecule of ~194 nt median insert, and the mates of a short-insert pair
overlap, so those pairs merge into a single read. Measured on this library with the
merge criterion used here, 35.1% of real pairs merge; the published self-overlap figure
is 26-27% (`results/quality_profile/`), measured with a different criterion, and the two
are not compared as if they were the same instrument.

The comparison is **paired**: merged read against mate 1 *of the same molecule*,
restricted to pairs that merge. Comparing merged reads to all mate-1 reads would compare
insert lengths and call it a result.

The effect is large and it grows with divergence (`TABLES.md` §4). Recognition at DIAMOND
bitscore >= 50, merged against mate 1 of the same molecule:

| divergence | merged (~140 nt) | mate 1 (~105 nt) |
|---|---|---|
| 0% | 91.1 (89.8–92.0) | 90.9 (89.6–92.0) |
| 5% | 87.2 (85.5–88.9) | 80.1 (79.2–81.5) |
| 10% | 75.5 (72.0–78.9) | 59.6 (55.7–64.3) |
| 20% | **38.4 (33.2–43.5)** | **10.5 (8.9–11.6)** |
| 30% | **7.1 (5.1–9.2)** | **0.0 (0.0–0.2)** |
| real library reads | 16.0 (14.1–17.7) | 8.9 (7.2–10.6) |
| shuffled null, same lengths | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) |

At low divergence merging buys nothing, because a 107 nt mate is already enough. In the
20-30% band — the band where every index in this project reads zero or nearly zero — it is
the difference between nothing and something: 10.5% to 38.4% at 20%, and 0.0% to 7.1% at
30%. On real library reads it is a 1.8x gain, and the profile side moves the same way
(8.6% against 5.5%).

**So part of the published 10-20% recognition wall is a read-length wall, not an index
wall.** That matters more than the size of the effect, because a read-length wall is one
we can move without new sequencing, new indexes or new money: the information is already
in the library we have, in about a third of its pairs.

Three things this does not say. It is one genome, a picornavirus with a single long
polyprotein ORF, which is favourable to both instruments. The merged reads here carry no
error-correction credit for their overlap — we have a measured error rate but no quality
model, and awarding merging an accuracy benefit we cannot falsify would be worse than
leaving it out — so the measured gain is a lower bound against a real merger. And the
simulated merge rate is about 43% against the library's 35%, so the simulation's merged
pile samples slightly shorter inserts than the real one; direction, that makes the
measured gain conservative.

## 3. What this changes

- **SPEC.md §4.2 and §7 need the profile-reach correction.** The axis stays as a reported
  field and loses its claim to reach.
- **Mate merging belongs in the pipeline before anything clever does.** It is free, it is
  measured, and it moves the one wall this project has described as a property of the
  indexes.
- **Nothing here helps the depth-1 pile find relatives it does not have.** Both results
  are about how far an instrument can see, not about abundance. The pile is still at depth
  1, and 85% of it still has no significant relative in core_nt. What merging buys is a
  better *description* of each single molecule, which is exactly the deliverable SPEC.md
  §4.1 specifies for this pile.

## Provenance

| file | script | what |
|---|---|---|
| `cells.json` | `experiments/profile_reach.py` | every cell of the design, 5 seeds |
| `paired_merged.json` | same | the paired merged-against-mate comparison |
| `TABLES.md` | `experiments/report_profile_reach.py` | generated from the two above |
| `manifest.json` | `experiments/profile_reach.py` | seed, git SHA, script content hash, operating points, timings |

Inputs: D00239 (ENA fasta) as the holdout; SRR23973193 (ENA fastq) as the known-material
control; RefSeq `viral.1.protein.faa.gz` plus UniProt Swiss-Prot as the translated-search
index; Pfam-A as the profile index. Versions and sizes are in `manifest.json`.
