# What assembly path selection costs: measured on a whole library, against 25 real genomes

**Run:** `experiments/path_settings_run.sh`, which calls `experiments/path_settings.py`,
`experiments/holdout_indexes.py`, `experiments/classify_contigs.py` and
`experiments/subtract_known.py`. Every table in this file is printed by
`experiments/report_path_settings.py` from the committed outputs; none is typed. Tables live
in [`TABLES.md`](TABLES.md).

## The question, and the objection that shaped it

`results/pair_ladder/` left a standing claim: below about 3.8% divergence every assembler
tested fuses lineages, PenguiN's raw output holds clean multi-kb pieces of distinct lineages
among mosaics, so *the reads carry the information and path selection is what is missing*.
That had only ever been tested on simulated two-lineage mixtures.

The obvious follow-up - re-assemble the residue at different path settings and see whether
it comes out better - is the wrong experiment, and Imo said why before it was run:

> splitting the residue from the rest feels like the wrong move when wiggling path settings,
> because the residue might be absorbed into the rest under laxer settings, no?

Residue membership is decided by the identification step, which ran on a previous assembly,
so a residue-only arm bakes in that assembly's path choices. Worse, the failure mode worth
hunting is invisible to it: under laxer settings residue reads can be pulled into an
abundant known contig instead of forming contigs of their own, and a residue-only assembly
cannot see that happen. So **every arm here assembles the whole library, and everything about
the residue is an output of the arm rather than an input to it.**

## What was run

Ten whole-library assemblies of SRR23421544 (4,290,429 single-end reads, the library the
post-cutoff holdout comes from), one per setting, then identification against indexes rolled
back to 2022-05-10, then the read-level cascade with each arm's own contigs as the
own-assembly tier. The two levers this project has already measured as dead - `--min-count`
and the k-list - are held at MEGAHIT's defaults in every arm, so the ladder varies path
selection and nothing else. `meta_sensitive` is carried separately because the preset also
changes `--min-count`, and no path conclusion rests on it.

Cost: zero. Everything ran in one container.

## The ground truth, which nobody chose

The depositors assembled 25 virus genomes out of this same water and deposited them on
2023-07-31, and `results/holdout_candidates/depth.tsv` records each one's depth in this run,
measured by mapping back: **0.4x to 700x, at 0.08% to 30.4% divergence from public
reference**. 53 further genomes in the same committed panel have **exactly zero** reads
mapped in this run, which makes a false-join control free and drawn from the same pool.

**The scorer reproduces the committed run before it is believed.** Against the two attributed
contigs of `results/holdout_cascade/`, it returns 3,990 nt covering **96.21%** of OQ722337.1
at purity 1.0, and 6,383 nt covering **99.41%** of OQ722330.1 - the published figures, from
different code.

## Two of my own measurement defects, found and fixed before any of this was believed

Both are the shape CLAUDE.md warns about - a flag that silently decides the result - and both
first produced a *more interesting* answer than the truth.

1. **`--secondary=no` credited only one member of a near-identical pair.** The panel contains
   at least six pairs above 99.4% identity, several co-present here at high depth
   (OQ722321/OQ722330 99.98% over 6,063 nt, at 499x and 700x). With secondary alignments
   suppressed, each contig was credited to its single best genome, so one member of each pair
   read 0.00 breadth. Read as data that says *a 499x genome failed to assemble*, and it also
   made the all-against-all pass over the panel report that no two panel genomes were related
   at all. With secondaries on, recovery at >=90% breadth rises from 9/25 to 20/25 at
   baseline, and the apparent run-to-run instability (one genome moving 32 points between two
   runs of the *same* setting) disappears entirely.
2. **Purity against one accession misread a partial deposit as absorption.** OQ722323 is a
   5,048 nt partial deposit 99.84% identical to the 9,467 nt OQ722327. A contig carrying the
   full molecule explains 0.53 of itself with the short record and 0.92 with the group's, so
   per-accession purity reported a correct recovery as absorption into foreign sequence.
   Purity is now measured against the near-identical group, with the accession-level figure
   kept beside it; where a genome has no near-identical relative the two are the same number.

## The answer

**Nothing in the ladder beats the default, and the levers that move anything move it down.**
At >=90% breadth of the genome: baseline 20/25, and its own repeat 20/25. `lax_bubble`,
`strict_bubble`, `keep_thin` and `prune_hard` are all 20/25 and move **no genome by as much
as 10 points in either direction**. `no_bubble` is 18/25, `keep_local` 16/25, `no_rescue`
15/25.

**The aggregate hides a trade, and the trade is the interesting part.** Counting genomes
whose breadth moves at least 10 points against baseline: `no_bubble` +1/-2, `keep_local`
+1/-4, `no_rescue` 0/-6, `meta_sensitive` +1/0, everything else 0/0. So path settings do
move specific molecules - they just move more out than in. The one real gain is worth naming:
**`no_bubble` is the only arm that sees OQ722315 at all**, a 20.7x genome 27.2% divergent
from public reference, going from 0.00 to 0.33 of its length, while costing 604 kb of
assembly, 89 of the 322 contigs over 3 kb, and two other genomes.

**The null is zero and that is what makes the rest readable.** `baseline_repeat` is the same
setting at a different thread count: 2 contigs different out of 35,328, 6,838 bases, 15 of
the 3,968 contigs over 1 kb, and not one genome moved by 10 points. So the four arms that
read 0/0 are genuinely flat, not noise-limited.

**Absorption was looked for explicitly and did not happen.** In every arm, including
`lax_bubble`, which exists to cause it: **zero** genomes recovered at >=90% breadth sit at
below 0.9 purity, **zero** contigs span two unrelated panel molecules, and **zero** of the 53
zero-read control genomes attract a contig. Whatever the bubble-merging similarity floor is
doing, at 0.90 it is not pulling these molecules into longer foreign contigs.

**What the ladder costs and buys at the thin end is the one large effect.** `no_rescue`
(`--no-mercy --no-local`) drops the assembly from 22.04 Mb to 14.13 Mb, loses 40% of the
contigs over 1 kb and six genomes, and runs in 3.9 minutes instead of 11.8. So MEGAHIT's two
low-coverage rescue mechanisms account for **36% of assembled bases** in this library. That is
a measurement of how much of this assembly exists only because of path rescue - and the
levers that would push further in that direction (`keep_thin`, `keep_local`) buy nothing:
`keep_local` adds 1,017 contigs and 172 kb of sequence and **loses four genomes**.

## What this cannot say

- **One library, one sample.** SRR23421544 is a treated-wastewater RNA library; Berlin's
  drift caveat applies here as everywhere. **The Berlin arm was not run.** The design called
  for the same ladder over the discriminating dataset, so that residue size on real
  background came out as an output of each setting; it was about eight hours of assembly and
  Imo chose to skip it once the holdout arm showed the levers flat. So **how much the residue
  of SRR23973193 moves with path settings is unmeasured**, and nothing here should be read as
  answering it.
- **It cannot answer whether both members of a near-identical pair come out separately.** At
  99.9% identity a 200 nt alignment at >=95% identity cannot say which member a contig came
  from. That question needs variant-level attribution, not alignment breadth, and so the
  `pair_ladder` fusion result is neither confirmed nor contradicted here.
- **The panel is the recoverable end of the problem.** These 25 genomes are what the
  depositors could already assemble. Nothing in this panel sits at the median depth 1 where
  `results/residue_census/` puts most of the residue, and no assembly setting can create
  coverage that is not there.
- **Identification here is not the production cascade.** The two 8 GB Kraken2 tiers are not
  staged: they are identical across arms and `results/subtract_recut/` measures their
  marginal contribution at 0.27 and 0.00 points. So a residue fraction in this directory is
  comparable across these arms and to nothing else.
- **The held-out index is not byte-identical to PR #29's**, for the same reason recorded
  there: upstream RefSeq viral keeps dropping records, so the rebuild reads 19,202 kept
  nucleotide records against 19,211 recorded and was built with `--allow-mismatch`. Direction:
  fewer pre-cutoff relatives make a candidate look more novel, the dangerous way.

## The read-level readout, which is where Imo's objection lands

Identification and the read-level cascade were run per arm, so residue size is an output of
each setting. Paired by seed - five seeds are five disjoint blocks of one pool, so the same
seed is the same reads - the baseline's residue is 35.42% of pairs and its own repeat moves
it by **-0.00 points (range -0.01 to +0.01)**. That is the instrument's noise floor, and it
is tight enough to read small effects against.

**`lax_bubble` takes 1.33 points of the library out of the residue** (range -1.30 to -1.36),
and **1.24 of those points come out of tier 2** - reads that previously matched only contigs
nothing could identify. Nothing in the ground-truth panel moves in that arm: no genome loses
breadth, no genome loses purity, no contig spans two molecules. This is exactly the shape the
objection predicted: on a residue-only design the setting would have looked like an
improvement, and the panel would not have contradicted it.

**But the contig-level test says the mechanism is not novel sequence being swallowed.**
`experiments/path_settings_absorption.py` aligns each arm's identified contigs against the
*baseline's* unidentifiable pile and asks how much of it now sits inside something the
cascade calls identified. Coverage-weighted: null `baseline_repeat` **0.18%**, `lax_bubble`
**0.26%**, `strict_bubble` **0.25%**, `no_bubble` **0.67%**. So the arm with the largest
contig-level absorption is the one that keeps paths *apart*, and `lax_bubble`'s figure is
1.4x a null that exists only because MEGAHIT's output shifts slightly with thread count.
0.26% of that pile cannot account for 1.24 points of reads.

The likelier mechanism, which this instrument cannot separate from absorption, is read
recruitment: a merged contig is a better magnet for a swarm than the variant-specific contigs
it replaces, and the cascade's own-assembly tier subtracts a read on 40 matching bases. Under
that reading the reads leaving the residue are variants of material the sample already holds,
and the residue shrinks without anything new being found. **Either way the number to report
is the pair, not the residue alone**: 1.33 points left the residue and nothing in the truth
panel improved.
