# Is a contiguous block of the FASTQ a random sample of the library?

`experiments/block_position_bias.py`, 5 seeds (the block statistics are exact over the
whole library and do not vary; the seed varies the random-sample control). **Every number
is in the per-seed JSONs** beside this file.

## Why it was asked

Every percentage this project reports comes from a contiguous block. `subtract_known.py`
skips `(seed - 1) * 1,500,000` pairs and reads the next `--pairs` in file order, so "five
seeds" means five regions of one flow cell, not five random draws of the library. The
per-read report then samples at random **within** that block, and the web view samples at
random within the report: the randomness is nested inside a selection that is not random,
and nothing here had checked the outer layer.

Two things could be wrong, and only the second turns out to be real.

## The instrument

The residue thread labelled all 8,731,224 pairs against SILVA 138.2 and all 82,947 of the
sample's own contigs, in file order. So "does content vary with position in the file" can
be answered exactly, over every pair in the library, with no alignment and no sampling of
its own. Pairs are scored by mate 1 alone, because the mates are not interchangeable.

The control is uniform random samples of the same size from the same library. It must read
1.0 on the dispersion statistic or the instrument is wrong rather than the blocks; **it
reads 0.92 to 1.10 across every class and every seed.**

## Location: the blocks are not in an unrepresentative place

| class | whole library | across all 87 blocks | block 1, z against binomial | the 5 cascade blocks, minus the library |
|---|---|---|---|---|
| `none` (residue upper bound) | 18.230% | 17.875 - 18.596% | +0.3 | **+0.024 points** |
| `rrna` | 38.760% | 38.249 - 39.275% | -1.1 | **-0.141 points** |
| `identified` | 33.323% | 32.743 - 33.860% | +0.8 | **+0.101 points** |
| `unident` (tier 2) | 8.675% | 8.445 - 8.954% | +0.5 | **+0.021 points** |
| `other_ctg` | 1.013% | 0.894 - 1.136% | -0.7 | **-0.005 points** |

Block 1 was the prime suspect — the start of the file is the first tile, where edge
effects and a flow-cell corner's optical duplicates live. It is not an outlier on any
class: every z is inside ±1.1, which is ordinary sampling noise.

**The five blocks the cascade actually uses sit within 0.14 points of the whole library on
every class, and within 0.03 points on both residue classes.** So no headline percentage
in this repo is meaningfully displaced by the block design. That is the question answered,
and it is answered in the direction that costs nothing.

## Dispersion: the intervals are too narrow, by a measured factor

The blocks do differ from each other by more than sampling noise, which is the half that
matters, because a Wilson interval assumes a simple random sample.

| class | overdispersion | lag-1 autocorrelation between adjacent blocks |
|---|---|---|
| `none` | **1.15x** | +0.04 |
| `unident` | **1.13x** | +0.02 |
| `other_ctg` | 1.28x | -0.02 |
| `rrna` | **1.43x** | +0.24 |
| `identified` | **1.57x** | +0.28 |

Overdispersion is the observed between-block standard deviation divided by the binomial
standard deviation a Wilson interval assumes. **1.57x on `identified` means the published
interval on that class is 57% too narrow**; for `identified` at n = 100,000 the half-width
goes from 0.29 points to 0.46.

The two readings agree with each other and say where the structure is. The abundant
classes carry it: adjacent blocks resemble one another (+0.24, +0.28), which is what
spatial structure along a flow cell looks like. **The two residue classes do not** — their
autocorrelation is within noise of zero and their overdispersion is 1.13-1.15x, so a block
is very nearly a random sample for exactly the classes this project is about.

## Reproduced independently, and one trap worth more than the result

The residue thread re-derived this from the same label table with its own script
(`experiments/block_dispersion.py`, `results/residue_census/block_dispersion.json`). Over
the same 87 windows it reads 1.44 and 1.61 on `rrna` and `identified` against the 1.43 and
1.57 here, with the autocorrelation splitting the same way, and five-block means within
0.16 points of the library. One figure differs: `unident` reads 1.24 there against 1.13
here, so treat that class's factor as 1.1-1.25 rather than a point value. Everything the
conclusion rests on is reproduced.

**The trap.** Computed from the five cascade blocks alone rather than from all 87 windows,
the same statistic reads **0.65 for `rrna` and 0.44 for `identified`** - "no inflation, if
anything less" - where the full series says 1.44 and 1.61. Five numbers cannot measure a
variance: the sampling error of a five-sample standard deviation swamps the effect and here
it reverses its sign.

That bears directly on this repo's convention. A headline is reported as a min-max across
five seeds, and that five-sample spread is **an honest observed range and a useless
dispersion estimate**. Any correction factor must come from the window series over the
whole library, which is where `hapax.py`'s `DESIGN_EFFECT` takes it from. A tool that
estimated its own design effect from its own five seeds would have confidently reported
that block sampling made its intervals too *wide*.

## What to do with it

1. **Do not switch the cascade to random sampling.** The location bias it would remove is
   at most 0.14 points, and reading a random sample out of a gzipped FASTQ costs a full
   pass over 8.7 million pairs per run. The measurement does not justify the cost.
2. **Widen the intervals by the design effect, per class**, or state it. 1.15x on the
   residue, 1.43-1.57x on rRNA and identified. This is a correction that can be applied
   to numbers already published without re-running anything, which is why it was worth
   measuring this way. Take the factor from the full window series and never from a run's
   own seeds, for the reason above.
3. **A five-seed spread across blocks is not five independent draws.** With adjacent-block
   autocorrelation at +0.28 on `identified`, five blocks spaced 1.5 M pairs apart are
   closer to independent than adjacent ones would be, but the spread they report is still
   a spread over one flow cell's structure, not over the sampling distribution.

## What it cannot say

It is blind to any bias that does not move these five classes: a block could hold a
different mix of viral families at identical ribosomal content and this would not see it.
The labels are a coarser instrument than the cascade — no Kraken2, no DIAMOND — so `none`
is an upper bound on the residue rather than the residue, and the design effect measured
on labels is carried across to cascade numbers as the best available estimate rather than
as the same measurement. And it describes SRR23973193. A second library would need its
own labels.
