# Where hapax's 97.0% on the Berlin astrovirus came from

Run with `experiments/ww_berlin_consensus.sh` after `experiments/ww_berlin.sh`. Per-base
attribution in `diag_*/` (from `experiments/consensus_diag.py`), table in `table.md`. Follows
`results/ww_berlin/REPORT.md`, which found hapax 97.0% identical to the authors' 6,823 bp SPAdes
astrovirus contig against MEGAHIT's 99.2%, and suspected the consensus step.

## Short answer

The per-column vote was not the problem. When every disagreeing base is checked against the reads,
two other causes account for nearly all of the gap:

1. **The library carries ~14 bp of non-template sequence on mate 2.** 2,455 of 2,572 recruited
   mate-2 reads are soft-clipped by minimap2 at the 5' end, mode 14 bp, at full base quality, so
   quality trimming would not remove it. hapax aligned every read end to end, so that sequence went
   into overlaps and into the pileup as mismatches and indels. Effects: junk indels of 9 to 241 bp
   in the consensus, fragmentation, and most of the high HMM fall-through the earlier report
   noted (31% of candidate pairs; 7.8% once the prefix is removed).
2. **One contig holds several astrovirus lineages.** The sample has at least five lineages that the
   authors' SPAdes assembled separately, 10 to 14% apart (e.g. NODE_1338 is 86.5% identical to the
   target over 2.1 kb, NODE_8170 is 89.9%). hapax's default overlap floor (`min_ident` 0.70, tier-1
   reject only above 30% mismatch) joins them into one component, and a per-column majority over
   several lineages follows whichever is locally commoner. On trimmed reads, **107 of 151 remaining
   substitutions are bases carried by another lineage's reads, while every read of the target's own
   lineage at that site carries the target's base**; 21 more sit where no target-lineage read maps
   to the hapax contig at all.

Shown cleanly: given only the target lineage's reads with the prefix trimmed, hapax's consensus is
**99.77%** identical to the authors' contig (6,300 bp covered, 3 contigs, longest 3,829 bp).

## Numbers

hapax on the 5,171 recruited reads, scored against the authors' contig (minimap2 `asm10`, identity =
matching bases over alignment block length, so indels count). "Before" is the binary at the start of
this work; every other row includes the consensus fixes below.

| run | contigs | longest block | covered | identity | HMM fraction |
|---|---|---|---|---|---|
| before, raw reads (reproduces the earlier 97.0%) | 4 | 3,901 | 6,722 | 96.97% | 0.309 |
| raw reads | 4 | 3,897 | 6,720 | 97.53% | 0.309 |
| raw reads, `--mode nucleotide` | 3 | 3,882 | 6,702 | 98.01% | 0.309 |
| mate-2 prefix trimmed | **2** | **6,079** | **6,809** | 97.14% | **0.078** |
| trimmed, `--mode nucleotide` | 2 | 6,079 | 6,809 | 97.82% | 0.078 |
| trimmed, `--min-ident 0.90` | 3 | 6,082 | 6,809 | 97.65% | 0.078 |
| trimmed, `--min-ident 0.95` | 2 | 5,985 | 6,300 | **98.88%** | 0.078 |
| trimmed, `--min-ident 0.95 --mode nucleotide` | 2 | 5,985 | 6,300 | 98.88% | 0.078 |
| target lineage only, raw | 5 | 2,023 | 6,208 | 95.06% | 0.276 |
| target lineage only, trimmed | 3 | 3,829 | 6,300 | **99.77%** | 0.037 |

For comparison (from `results/ww_berlin`): MEGAHIT on the whole raw sample 99.2% over 12 contigs,
6,709 covered; authors' SPAdes 100% by construction. MEGAHIT on the trimmed sample: see the end.

Reading the table:

- Trimming is what buys contiguity: one 6,079 bp block covering 99.8% of the genome, against
  4 contigs before. Identity does not rise with it, because longer contigs absorb more lineages.
- The identity floor is what buys identity: 97.1% → 98.9% at 0.95. Remaining substitutions at 0.95
  are still mostly lineage (51 of 71); lineages at ~10% can share stretches above 95%.
- **Codon mode costs identity on this sample at the default floor** (97.14% codon vs 97.82%
  nucleotide, trimmed). That is synonymous collapse doing what it is designed to do, joining reads
  that differ synonymously, applied across lineages that are different genomes. At a 0.95 floor the
  two modes are identical, because almost nothing reaches the codon HMM. Recorded per CLAUDE.md:
  on this real sample the codon advantage is either harmful or absent.
- Target-lineage reads alone, untrimmed, are *worse* (95.1%) than the mixture, because the prefix
  junk dominates when depth is lower: a 241 bp deletion and 9 to 23 bp indels. Trimmed, 99.77%.

## Per-base attribution (`diag_*/manifest.json`)

| run | substitutions | from another lineage | no target-lineage reads | target lineage agrees with hapax | other | indel events (bases) |
|---|---|---|---|---|---|---|
| raw | 134 | 81 | 31 | 14 | 8 | 8 (39) |
| trimmed | 151 | 107 | 21 | 13 | 10 | 1 (45) |
| trimmed, 0.95 | 71 | 51 | 9 | 7 | 4 | 0 |
| target lineage only, raw | 9 | n/a | | | | 8 (307) |

"Target lineage agrees with hapax" means the comparator's own base is the minority among its own
lineage's reads: a polymorphic site or an error in the comparator, not in hapax.

## What changed in the code

- `consensus::align_to_ref` may now leave either read end unaligned for a fixed cost (`CLIP` = 10,
  minimap2's short-read end bonus). A tail is clipped only if aligning it scores below -10, which a
  genuine SNP near a read end (-3) never does. The read is kept whole; clipped bases are not voted.
  Unit test `non_template_read_ends_are_clipped_not_voted`.
- Reads that fail to realign in one polishing round keep their placement and are retried against
  the next draft. Before, a read missed in round 0 (the noisiest draft) was gone for good.
- `orf::call_orfs` returns no ORFs for a sequence shorter than a codon or than `min_len`. A component
  whose reads all failed to realign left an empty draft and panicked the run; any tighter overlap
  setting (e.g. `--min-ident 0.90`) hit it on this sample.
- `--rounds` exposes the number of polishing rounds.

Effect on the simulated diversity ladder (5 seeds × 4 diversities × codon/nucleotide × tau 0.9/1.0,
80 runs, against `results/ladder_v2`): contigs, contigs per entity, chimaeric contigs and missed
entities identical in all 80; mean recovered fraction unchanged in 68 and higher in 12
(mean 0.9792 → 0.9802). The simulator has no non-template ends, so little was expected.

## Tried and dropped

A read-weighted vote (each read weighted by `exp(-divergence/0.02)` against the current draft, so
the vote follows one lineage through a read's length). On raw reads it made identity worse
(97.5% → 94.9%): a read that disagrees with a draft error is exactly the read that gets
down-weighted, so draft errors lock themselves in. On trimmed reads it moved identity by 0.2 points
at most. Not kept. Splitting lineages needs to happen on linked variant sites, not on whole-read
divergence against a mosaic draft.

## What this leaves open

- **The real decision is the join criterion, not the consensus.** CLAUDE.md says reads are joined
  only when they share nucleotide identity, but the default floor (0.70) and the simulator's
  notion of one population (up to 0.20 subs/site within, sisters at 0.30) let 10-14% lineages
  join. On real wastewater those are different genomes, assembled separately by SPAdes. Either
  the floor rises (costs coverage: 6,809 → 6,300 here) or the consensus has to split lineages
  within a component on linked sites. The ladder should gain a rung with distinct lineages at
  10-15% inside a mixture, since it cannot currently see this failure.
- hapax has no read-end preprocessing. The trimmed runs here strip a fixed 14 bp from mate 2,
  chosen from the measured soft-clip mode; a general fix would detect non-template ends from the
  overlaps themselves.
- One sample, one comparator that is itself an assembly, lineage labels from best-bait mapping
  (which favours the authors' contig for reads near it). The "target lineage only" 99.77% is
  therefore an upper bound on what a perfect lineage split would give, not a fair head-to-head.
