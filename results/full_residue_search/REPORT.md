# The whole unplaced residue, searched

`experiments/full_residue_search.py`, seed 1, SRR23973193, local, no spend. 900 s wall on
4 cores. Manifest and counts: `manifest_s1.json`; the ranking built on top of it:
`ranking_summary.json`. The 41,829-row subject table and the per-read output stay in the
run directory, not here.

## What was searched

Every unplaced pair in the library: **1,372,021 pairs, 2,744,042 mates**, against
`viral.dmnd` and `uniprot_sprot.dmnd` at `--sensitive`, reporting at `--evalue 100` and
`--max-target-seqs 1`, with significance applied afterwards at **E < 1e-5**. A control of
**200,000 already-placed pairs** went through the same search in the same pass.

This replaces a 6,000-pair sample. The leaderboard it feeds had leaders carrying two and
three reads, and medians of a handful.

The pile comes from `/mnt/project-files/shared/resasm/library_label_classes.npz`, the
whole-library nucleotide label table, not from re-running the cascade: a warm 100k-pair
cascade did not finish in 35 minutes here, so 8.7 M pairs is tens of hours. **Direction of
the resulting bias: the table is nucleotide-only, so its unplaced pile is a little more
generous than the cascade's residue and the named-relative rate below is, if anything, a
small overestimate. One direction, not two.** Its unplaced count reproduces
`results/residue_census/` §2 to the pair.

## The numbers, each with its control

| | leftovers | already-placed control |
|---|---|---|
| mates searched | 2,744,042 | 400,000 |
| any hit reported at E<100 | 1,273,084 (46.39%) | 183,897 (45.97%) |
| **significant at E<1e-5** | **214,502 (7.82%)** | **39,992 (10.00%)** |
| distinct subjects | 41,829 | — |

**92.18% of the leftovers point at nothing we can name**, and the control says what that
figure is read against: material the cascade already identified clears the same threshold
only 10.00% of the time, because most of it was identified by a nucleotide index that
protein search does not reproduce. The gap between the two piles is 2.2 points, not the
order of magnitude the words "identified" and "leftover" suggest.

**The "any hit" row is the chance band, measured here by accident and cleanly: 46.39%
against 45.97%.** Two piles that differ by 2.2 points at a significance threshold are
indistinguishable at the reporting threshold, which is what `results/weak_hits/` found
with a shuffled null (92.3% of shuffled reads hit at this cutoff). The row is reported
so that nobody reads a database-size curve of it as reach. It is not evidence of
anything.

## What the subjects are

Of the 214,502 named leftover mates: **166,448 (6.07% of all leftovers) point at a
putative phage**, 16,530 (0.60%) at a virus with no phage signal, 13,914 (0.51%) at
something living whose kind could not be placed, 10,628 (0.39%) at something putatively
bacterial or archaeal, 6,982 (0.25%) at something putatively eukaryotic.

**The host axis carries nothing here and that is the finding.** Host confidence is a
ramp over identity and aligned length. Over all 41,829 subjects its **maximum is 0.40
and its median 0.17**: not one subject in the whole residue is close enough to a known
organism for that organism's host to transfer. The host groups are labels on the matched database record, not claims about
the water, and the page says so above the colours.

## Two classification bugs this scale exposed

Both were invisible at 243 subjects and obvious at 41,829, which is the argument for
running the whole pile even when a sample would do.

1. **`Bacteriophage X` was grouped `unknown_host`.** The keyword `phage` is anchored at a
   left word boundary - the boundary that stops it firing inside `phagocyte` - so it
   cannot match inside `Bacteriophage`, which is RefSeq's commonest phage naming. Fixed by
   listing the `-phage` compounds explicitly rather than relaxing the boundary, in both
   the host rules and the kind rules, which had disagreed about the same word.
2. **`Rodent astrovirus` was grouped `human_associated`.** `Astrovir` is in the human list
   and the rules were tried in file order. Rules now carry a `priority`: a record that
   names its host animal beats a virus-family prefix, because the family says what the
   virus is and the host word says who it was found in. It only fires when the name
   carries the animal - `Norovirus GV` is a mouse virus and nothing here can tell.

## What this does not say

Nothing here scores novelty. An index can say whether something is already known; it
cannot say whether something is new, and the ordering this feeds is a triage order for a
person's attention built from four measured quantities. The 92.18% is a statement about
what the references can answer, not about what is in the water.
