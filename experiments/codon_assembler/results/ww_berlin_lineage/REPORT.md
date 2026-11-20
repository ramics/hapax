# Splitting lineages on linked sites (`--split-lineage`)

Follows `results/ww_berlin_consensus/REPORT.md`. That report found that most of the remaining
identity gap on the Berlin astrovirus (SRR23973193) came from one contig holding several lineages
10 to 14% apart, which the per-column vote then blends. Imo chose to fix this by keeping the
permissive overlap join and splitting a component into lineages using variants linked within reads.
This report covers that split, `src/lineage.rs`, which is **off by default**.

Berlin: `experiments/ww_berlin_lineage.sh`, tables `target.tsv`, `lineages.tsv`,
`multi_lineage.tsv`. Ladder: `experiments/diversity_ladder.py`, 5 seeds × 4 diversities, tau 0.9,
results in `results/ladder_lineage/` (`comparison.tsv`, one directory per setting with its manifest).

## Short answer

- **On Berlin it works.** The target lineage comes out as one 6,075 bp contig at **99.63%**
  against the authors' contig (hapax default 97.81%, MEGAHIT on the same trimmed reads 98.92%).
  Nine SPAdes lineages are recovered as separate contigs, against three with the split off.
- **On the simulated ladder it over-splits.** A single within-host swarm comes out as 2.0 to 4.4
  contigs per entity, against 1.0 to 1.7 with the split off. It does remove every chimaera (0 in
  all 160 runs, against up to 0.07 per run with it off) and every missed entity.
- **No floor serves both.** `--lineage-split-div` (merge split lineages closer than this) reduces
  the over-splitting only from 0.09 up, and at 0.09 Berlin's lineages merge again (target 98.66%,
  NODE_5252 95.19%). The swarms in the ladder have 5 to 17% divergence between haplotype chains, and
  Berlin's lineages are 7 to 14% apart. The ranges overlap, so a divergence threshold cannot tell
  "one swarm" from "two genomes". That is the same definitional conflict the consensus report
  raised, now measured.

So it ships off by default. It is worth turning on for a wastewater sample where separate genomes
are expected. It is not worth turning on for a within-host swarm.

## Berlin (SRR23973193, recruited astrovirus reads, mate-2 prefix trimmed)

Against the target (`target.tsv`, asm10, identity over alignment blocks, indels count):

| run | blocks | identity |
|---|---|---|
| hapax default | 2 | 97.81% |
| `--split-lineage 1` | 2 | **99.63%** |
| `--split-lineage 1 --lineage-split-div 0.06` | 2 | 99.63% |
| `--split-lineage 1 --lineage-split-div 0.09` | 3 | 98.66% |
| `--split-lineage 1 --mode nucleotide` | 2 | 99.63% |
| MEGAHIT, trimmed | 2 | 98.92% |

The default row is 97.81% rather than the 97.14% in the consensus report because the consensus
now fills the holes that clipped read ends leave (below).

Per lineage (`lineages.tsv`; each contig's primary alignment against all SPAdes astrovirus contigs):

| lineage | split: blocks, identity | MEGAHIT: blocks, identity |
|---|---|---|
| NODE_37 (target) | 2, 99.63% | 4, 97.44% |
| NODE_1118 | 4, 99.72% | 1, 99.92% |
| NODE_1338 | 1, 97.78% (one 46 bp insertion; 1 substitution) | 1, 99.94% |
| NODE_1471 | 2, 99.78% | 1, 99.84% |
| NODE_1845 | 1, 96.13% | 1, 98.56% |
| NODE_5252 | 1, 99.87% | 1, 99.57% |
| NODE_8170 | 1, 100% | 1, 99.64% |
| NODE_7719, 13532, 60396 | 100% each | 100% (7719 only) |

hapax with the split wins clearly on the target lineage. MEGAHIT wins on contiguity for the minor
lineages and on NODE_1845. One hapax contig (`multi_lineage.tsv`) runs through NODE_7719 →
NODE_13532 → NODE_60396 → NODE_1845, each piece 100% or close. That is either a genome SPAdes broke
into four or a chimaera; the comparator is itself an assembly, so this sample cannot tell which.
It is reported as a candidate chimaera.

## Ladder (32-tip swarms, one population per entity; `results/ladder_lineage/comparison.tsv`)

Contigs per found entity, mean over 5 seeds (codon mode; nucleotide mode is within 0.1 except at
0.20, where it reads 4.35/4.10/3.55):

| diversity | split off | split, floor 0 | floor 0.06 | floor 0.09 |
|---|---|---|---|---|
| 0.05 | 1.00 | 2.00 | 1.50 | 1.05 |
| 0.10 | 1.05 | 3.25 | 2.95 | 1.65 |
| 0.15 | 1.25 | 3.10 | 3.00 | 2.55 |
| 0.20 | 1.63 | 3.70 | 3.45 | 3.00 |

Chimaeric contigs per run go from 0 to 0.07 (split off) to 0 in every split run. Missed entities
go from 0 to 0.2 to 0. Recovered fraction is 0.995 or higher in every split cell except 0.05 at
floors 0.06 and 0.09 (0.975, against 0.999 off).

The split is doing what it is built to do: haplotype chains in these swarms really are 5 to 17%
apart over hundreds of bases, and the reads link each chain's sites. The simulator calls them one
entity, and Berlin's SPAdes calls the same structure several genomes.

## How the split works (`src/lineage.rs`)

1. **Informative sites:** columns where the minor allele is carried by at least 3 reads and 10% of
   depth.
2. **Read groups:** reads that overlap and agree at every shared informative site (at least 2) are
   joined, best agreement first. A join must also pass a profile test against the group so far
   (conflicts at most 4% of the overlap, or 30% of the shared sites). A site counts only where both
   sides have two or more reads, so one read's error cannot veto a join.
3. **Colouring:** groups, largest first, join the lineage they fit best on shared sites. A site
   conflicts only when the group's allele is carried by under a fifth of the lineage's reads, so
   sister haplotypes inside one lineage are not counted as conflicts. A group that shares no site
   with any lineage joins the biggest lineage whose span it lies in. Otherwise it starts a piece.
4. **Stitching:** pieces are joined end to end across a stretch with no informative site only when
   the join is unique both ways, counting only the nearest continuations. Where two lineages meet
   and no read says which continues which, the pieces stay apart.
5. **Reads without sites** go to the nearest lineage, then the one with fewer conflicts, then the
   one the read lies deeper inside. They are placed nearest first, so each lineage grows from its
   own side. Stitching then runs again, because the grown spans can make a junction unique.
6. Lineages below 30 reads or 2% of the component are folded back. Each set is rebuilt into its
   own consensus. `--lineage-depth` repeats the split on each piece (no change on Berlin).

Two earlier versions were dropped. The first stitched by chains, which fragmented one lineage
across stretches with no sites. The second coloured by span overlap, which chained minor lineages
across gaps and blended NODE_1338 with another lineage (95.1%).

The consensus change that ships with it: when clipped read tails or heads agree on sequence that
the draft lacks (at least 4 reads and 25% of local depth), that sequence is inserted or bridged.
The clip fix in PR #4 had left a 45 bp hole in the target contig, and this fills it. With the split off, the simulated
ladder (tau 0.9) matches `results/ladder_v2` on contigs, contigs per entity and chimaeras in every cell.
The recovered fraction is equal or slightly higher (by up to 0.004).

## What this leaves open

- **The real fix is a definition, not a threshold.** "One entity" in the simulator spans the same
  divergence as "separate genomes" in wastewater. Either the ladder gains a rung with distinct
  lineages mixed together (so the split can be scored where it should fire), or the split uses a
  signal other than divergence. Candidates are coverage steps at lineage boundaries, or whether
  the chains co-vary along the whole genome.
- The NODE_1338 46 bp insertion and NODE_1845's 15 bp of indels have not been traced to reads.
- One real sample, comparators that are themselves assemblies, labels from best-bait mapping.
