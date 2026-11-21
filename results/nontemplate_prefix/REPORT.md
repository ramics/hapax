# The 14 nt trim, tested against a reference

`experiments/nontemplate_prefix_map.py`, 5 seeds of 100,000 pairs each, disjoint blocks of
the file. One commit (`6490dfb`), clean tree, one script hash in all five manifests.

## Why this was worth running

Every experiment in this repo trims 14 bases from the 5' end of mate 2. That trim is
load-bearing: it is why R2 is 93 nt against R1's 107, and that length difference is the
single cause of every pair asymmetry `CLAUDE.md` warns about — asymmetric absence leaning to
R2 at 64:36 while a fixed-width error signature leans to R1 at 58:42, from one mechanism
pointing two ways.

**A correction first.** The trim had been checked against a reference once before, and my
commit message for this script said it had not. The earlier check is in
`experiments/codon_assembler/results/ww_berlin_consensus/REPORT.md`: 2,455 of 2,572 mate-2
reads recruited to one astrovirus bait were soft-clipped at the 5' end, mode 14, at full base
quality. That is where the number came from, and it was right.

What was genuinely open is narrower, and `results/dna_rna_pair/REPORT.md` records it as
UNRESOLVED: that check measured only that mate-2 reads clip. It did not have mate 1 as a
control, did not have the other end of the same read as a control, was 2,572 reads recruited
to a single virus, and said nothing about any other library. Meanwhile the composition
detector in `dna_rna_pair.py` fired on all four mate-and-library combinations and read
**9 for R1 and 16 for R2** — numbers nothing had tested.

**Composition cannot settle it and mapping can.** Illumina's 5' priming bias and a genuine
non-template prefix both make the first bases unlike the body of the read, which is all
composition sees. They differ in one respect it cannot see: priming bias changes *which*
molecules get sampled, not whether their first bases came from the template. Primed bases
still align to the genome; non-template bases cannot. Map the untrimmed read and it becomes
a clipping question.

## The table

Soft-clip at each end of the best anchored SILVA alignment, untrimmed reads, ≥60 matching
bases required before an end is read. Range over 5 seeds.

| | R1 5' | R1 3' | R2 5' | R2 3' |
|---|---|---|---|---|
| **RNA 2021** (SRR23973193) mode | **0** (5/5) | 0 | **14** (5/5) | 0 |
| mode share | 0.569–0.574 | — | 0.507–0.518 | — |
| clip ≥ 10 bases | 0.087–0.093 | — | **0.946–0.948** | — |
| clip exactly 0 | 0.569–0.574 | 0.904–0.910 | ~0.02 | 0.945–0.948 |
| mean 5' minus mean 3' | **1.89–2.07** | | **13.53–13.69** | |
| **RNA 2022** (SRR23973265) mode | **0** (5/5) | 0 | **14** (5/5) | 0 |
| clip ≥ 10 bases | 0.085–0.088 | — | **0.952–0.955** | — |
| mean 5' minus mean 3' | **1.20–1.28** | | **13.54–13.63** | |

Anchored reads per seed: 34,986–35,456 (R1) and 33,851–34,330 (R2) for 2021; 61,821–62,181
and 60,005–60,249 for 2022.

## 1. The trim is real, and the control is the other end of the same read

**Mate 2's 5' end clips at 14 bases in half of all anchored reads and at ten or more in 95%
of them. The same reads align flush at their 3' end — 94.6% clip zero there.** A divergent
reference entry would clip both ends; a quality problem would clip the 3' end harder, not
the 5'. One end of one mate, at one length, is what a non-template prefix looks like and
nothing else here does.

Mate 1 is the second control and it separates cleanly: mode 0, 57% of reads clipping nothing
at all, and a 5'-minus-3' excess of about 2 bases against mate 2's 13.6.

## 1a. A correction to this report's own first conclusion

I first read §2 below as meaning `dna_rna_pair.py` had been **over-trimming mate 2 by 2
bases**, corrected it to 14 and re-ran. **The Kraken2 removal fraction went down** — 1.66
points on the 2022 RNA library (0.4065 → 0.3899) and 1.06 on the 2021 one (0.2265 →
0.2159) — from being given two *more* bases of what the mapping calls template. That is the
wrong direction for a correction, so "over-trim" is withdrawn.

Two bases of real template cannot make a read harder to classify. What can is that **the
prefix length is a distribution, not a constant**: mode 14 at 51% of anchored reads, with 13
at 16%, 12 at 6.5%, 15 at 6.5%, 16 at 2.2% and a tail. Kraken2's confidence is the *fraction*
of a read's minimizers falling in the called clade, so two extra minimizers at the end of the
read where the prefix variability lives lower every read's confidence a little and push
marginal ones below the threshold. A 5,000-pair check reproduces it directly: 0.3844 at trim
14 against 0.4048 at trim 16, mate 1 untouched and the same block at both rungs.

**So the mode answers one question and not the other.** Trimming at 14 is the best estimate
of the prefix length, which is what §1 measures. Trimming above it is the better operating
point for a fraction-thresholded classifier. The trim is therefore an operating point in this
repo's sense, and `experiments/trim_sensitivity.py` measures its curve
(`trim_sensitivity_*.json`) rather than leaving it as a decision somebody made once — 1.7
points is the size of effects reported here as findings.

## 2. The composition detector was wrong in both directions, and this is why it is a flag

`dna_rna_pair.py` read **9 bases on R1 and 16 on R2**. Mapping reads **0 and 14**, in both
libraries, in all five seeds. So the detector is not reading the prefix length: there is no
fixed-length prefix on R1 at all, and R2's is 14.

But as §1a records, that does **not** make its 16 a worse trim — on a fraction-thresholded
classifier it is a better one. The detector was right by accident and for a reason it cannot
state, which is the definition of a flag. Treat `nontemplate_prefix` as saying that something
is compositionally odd at a read's start, which is true and useful, and never as a length.

## 3. Mate 1 is not perfectly clean, and the difference in shape is the point

R1 carries a real 5' excess — about 2 bases on average, with 8.7–9.3% of reads clipping ten
or more. But it has **no mode**: the distribution runs 0, 3, 5, 6 with no spike. A fixed-length
non-template prefix produces a spike; 5' priming and end-quality effects produce a tail. R1
has the tail and R2 has the spike on top of a tail. This is exactly the distinction
composition could not make, and it is why trimming R1 by a fixed amount would be wrong.

## 4. It is a property of the prep, not of the run

The 2022 library, eight months later at the same plant, gives the same mode at the same
length with the same 95% of reads clipping ≥10. **Every library in this series should be
assumed to carry it**, which matters because the 2022 libraries are new to this project and
`experiments/dna_rna_pair.py` trimmed R2 by 16 on a detector reading rather than by 14 on a
measurement. It now takes `--r2-trim`, defaulting to the measured 14, with a guard that stops
the run if the detector disagrees by more than 4 — but see §1a: 14 is the measured prefix and
16 classifies better, so the default is the honest estimate rather than the optimum, and the
curve beside it says what the choice costs.

## 5. The coordinate assumption was checked, not assumed

This reads mappy's `q_st` as a coordinate on the query as supplied, so that the 5' clip is
`q_st` on both strands. If that were wrong, minus-strand hits would carry the mirror image
and the two ends would average into looking alike — the failure mode `CLAUDE.md`'s "read the
control first" rule exists for. Every 5' distribution is therefore split by strand.
**R2's mode is 14 on both strands** (plus n≈32,000, minus n≈1,600) and R1's is 0 on both. No
mirror, so the coordinates are what the code assumes.

The strand counts are themselves consistent with a stranded library: R1 maps almost entirely
minus and R2 almost entirely plus, as mates of one fragment must.

## 6. Two things this run also happens to calibrate

- **The anchored fraction is an independent rRNA measurement.** 35.0–35.5% of 2021 pairs and
  61.8–62.2% of 2022 pairs anchor ≥60 bases to SILVA. The same libraries read 38.01% and
  67.12% at the ≥40-base threshold in `results/dna_rna_pair/`, and 34.5% by exact 31-mers in
  `results/subtract_known/`. Three instruments at three thresholds, ordered as thresholds
  predict. The ~29-point drift between the two libraries is reproduced here by a third
  method.
- **The mapq floor does not matter.** A sensitivity run at `--min-mapq 30`
  (`mapq30/nontemplate_prefix_s1.json`) keeps only reads hitting a unique region of a highly
  redundant reference, which costs 92% of the n: 2,940 anchored R1 reads instead of 35,000.
  The answer does not move — R2 mode 14, 5'-minus-3' 13.64 against R1's 1.63. So the clipping
  is a property of the reads' ends and not of which rRNA entry was hit.

## What this cannot say

- **The DNA library has no answer here and the row says so.** SILVA anchored 0 to 11 reads
  of 100,000 across the five seeds, because that library is 0.01% ribosomal. Its modes swing
  between 0, 12 and 42 on single-digit n, which is the expected behaviour of a measurement
  with no data behind it, and is reported with its n so it cannot be read as one. No
  nucleotide reference in this container identifies enough of a DNA virome to replace SILVA.
- **Ribosomal reads only.** The anchor is SILVA, so this measures the prefix on the
  ribosomal fraction. There is no reason for a library-prep artefact to care which molecule
  it is on, but this run does not prove it does not.
- **Mate-level, not pair-level.** Each mate is mapped independently and no pair rule is used.
- **The bias runs one way.** Adapter read-through clips 3' ends and inflates the within-read
  control, which can only make a real 5' excess harder to see. Divergence from the SILVA
  entry clips both ends and cancels in the difference.

## The trim is an operating point, and 16 beats the measured 14

`experiments/trim_sensitivity.py`, 5 seeds per library, 100,000 pairs per seed, Kraken2
PlusPF-16 at confidence 0.1, mate 1 untouched at every rung. Mean across seeds
(min–max), as a percentage of read pairs identified.

| mate 2 trimmed by | RNA 2021-10-13 (ours) | RNA 2022-06-16 |
|---|---|---|
| 0 (untrimmed) | 19.93 (19.89–20.02) | 35.86 (35.74–35.95) |
| 10 nt | 21.83 (21.79–21.87) | 39.29 (39.02–39.45) |
| 12 nt | 21.67 (21.63–21.69) | 39.00 (38.71–39.16) |
| 14 nt | 21.67 (21.63–21.73) | 38.86 (38.58–39.04) |
| 16 nt | 22.74 (22.71–22.78) | 40.59 (40.39–40.73) |
| 18 nt | 22.33 (22.30–22.37) | 39.84 (39.58–39.97) |
| 20 nt | 21.91 (21.85–21.97) | 39.12 (38.85–39.24) |
| 24 nt | 21.44 (21.39–21.47) | 38.23 (38.03–38.31) |

**The best rung is 16 in 5 of 5 seeds in both libraries**, and the ladder spans 2.7–2.9
points on the 2021 library and 4.65–4.92 on the 2022 one. Moving from the measured prefix
mode of 14 to 16 buys 1.03–1.11 points here and 1.60–1.84 points there.

This is the measurement behind §1a's withdrawal. The prefix is a distribution with mode 14,
so 14 is the right answer to "how long is the prefix" and the wrong answer to "what should
Kraken2 be given": confidence is a fraction of a read's minimizers, so the residual
non-template bases left by a modal trim lower every read's confidence, and cutting two
template bases costs less than leaving two non-template ones. Over-trimming past 16 does
turn around, as it must.

**So the trim is two operating points, not one.** Anything that reads SEQUENCE — coverage,
composition, folding, the per-read fields — takes 14, because those two bases are real
template. Anything that reads a CLASSIFIER's confidence takes 16. `hapax run` currently
passes 14 to everything, which costs about one point of Kraken2 identification on this
library and is recorded in every run manifest rather than left implicit.

**The ladder is not monotonic on either side of the peak**, and that is in the data rather
than in the noise: trim 10 reads 21.83% against 12 and 14 at 21.67%, a gap of 0.16 points
with a seed-to-seed range under 0.06. Minimizers are selected in windows, so two bases can
move whether a window's minimizer is template or not, and the effect need not be smooth in
the trim. It means the peak must be found by sweeping rather than reasoned to, which is
what this experiment is for, and it means a trim chosen on one library should be re-swept
on the next.

