# Near-relative separation ladder: where incumbent assemblers stop keeping two lineages apart

Scripts: `experiments/select_pairs.py` (pair choice), `experiments/pair_ladder.py` (runs),
`experiments/report_pair_ladder.py` (tables, `summary.md`). Rows: `rows.jsonl`, 2,205 runs,
0 tool errors. Manifests: `manifest_rows_part1.json` (first 1,667 runs) and `manifest_rows.json`
(the resumed remainder after a worker restart; same code, same seeds). Contigs and per-contig
labels for every run: `/mnt/project-files/strain_ladder/contigs/`.

## Design

- **Pairs.** Seven divergence rungs from 2.0% to 9.7%, three real rabies genome pairs per rung,
  genome-disjoint within the ladder except one flagged reuse at 5.5% (`ladder.json`). The
  archived ladder had one pair per rung, so divergence was collinear with pair identity.
- **Reads.** ART HS25, 2x150, insert 400±40, spike only (no background), equal abundance,
  5x, 15x or 40x per lineage. Seeds 1-5. Every tool sees byte-identical reads.
- **Tools at defaults.** MEGAHIT 1.2.9, metaSPAdes 4.3.0 (`--only-assembler`), PenguiN
  `nuclassemble` (raw) and `guided_nuclassemble`, Haploflow. MEGAHIT also assembles each lineage
  alone as a control (630 runs).
- **Labels.** The strand-aware 300 nt window caller audited in `results/labeller_audit`. A
  contig is *pure* if every called window is one lineage. Headline readout: **separated breadth
  of the weaker lineage**, the fraction of that genome covered by pure contigs of it. Mosaic mass
  and lineage switches per 10 kb come from the same runs.

## The labeller's null on real assembler output

**0 of 630** single-lineage MEGAHIT assemblies contain a contig called mosaic (mean breadth
0.975). The archived caller reported a mixture in 13 of 42. So mosaic counts here can be read
directly.

## Headline, 15x per lineage

Mean over 3 pairs × 5 seeds. Brackets give the range of the three per-pair means.

| divergence | MEGAHIT sep. breadth | MEGAHIT mosaic mass | metaSPAdes sep. | PenguiN guided sep. | PenguiN raw sep. (longest pure piece, median) |
|---|---|---|---|---|---|
| 2.0% | 0.00 [0.00-0.00] | 0.98 | 0.00 | 0.00 | 0.98 (2.9 kb) |
| 3.8% | 0.08 [0.08-0.09] | 0.78 | 0.00 | 0.00 | 0.99 (4.2 kb) |
| 4.4-4.5% | 0.52 [0.23-0.95] | 0.34 | 0.14 | 0.13 | 0.99 (7.4 kb) |
| 5.2-5.6% | 0.97 [0.93-0.99] | 0.01 | 0.04 | 0.33 | 1.00 |
| 6.5-6.8% | 0.99 [0.99-1.00] | 0.00 | 0.04 | 0.77 | 0.99 |
| 7.3% | 0.88 [0.80-1.00] | 0.12 | 0.00 | 0.70 | 0.99 |
| 9.7% | 0.99 [0.99-0.99] | 0.00 | 0.66 | 0.99 | 0.99 |

Going to 40x per lineage improves separation for no tool: every 40x cell lies within the
spread of its 15x counterpart (`summary.md`). At 5x nobody separates fully: MEGAHIT's separated breadth rises steadily from 0.00 at
2.0% to 0.86 at 9.7%, limited by breadth as much as by collapse.

## What it says

1. **MEGAHIT collapses two lineages at 3.8% and below, at any depth tested.** The output is
   one consensus that switches lineage about every 0.8-1.2 kb (8.6-12.9 switches per 10 kb).
   From 5.2% up, with 15x or more, it returns both lineages separately. **4.5% is the
   transition, and the pair matters as much as the divergence there:** per-pair means run from
   0.23 to 0.95 at 15x on pairs within 0.1 points of each other. That spread is exactly what
   the one-pair-per-rung ladder could not show.
2. **metaSPAdes collapses up to 7.3%** and only partly separates at 9.7%. That fits its
   design, which removes strain bubbles on purpose, and the stress test's earlier finding.
3. **PenguiN, the closest published claim, does not separate either in its usable mode.**
   Guided mode, with its redundancy reduction, returns a median of 3 contigs holding about two
   genomes' worth of sequence, and collapses at 3.8% and below. Raw mode is different and is the interesting result. It emits
   a median of 1,016 contigs per run, about 113 times the mass of the two genomes, and 70-83% of that mass is
   mosaic at 2-3.8%. But among them are pure pieces of *both* lineages, a median 2.9 kb long
   at 2.0% and 4.2 kb at 3.8%. Their union covers 98-99% of each genome. Nothing in its output
   says which contigs those are.
4. **So below about 4.5% the reads carry enough to separate the lineages over kilobases**, as
   PenguiN's raw paths show, **but no tool here selects those paths.** The open problem is
   choosing a consistent path, not getting the information. That is the target for this
   direction: 2-4.5% divergence at 15x or more per lineage, where every incumbent returns a
   mosaic and the right pieces demonstrably exist.

## Real-background parity (MEGAHIT, 2.0% and 4.5%, 15x)

Does the spike-only shortcut change the answer? The same 30 MEGAHIT mixture runs (6 pairs × 5
seeds, byte-identical spike reads) were repeated with 1,000,000 real read pairs from the Berlin
wastewater run SRR23973193 added, the first million in file order, with the ~14 nt non-viral
mate-2 prefix trimmed (PR #4). Before labelling, contigs that do not align to either lineage over
≥80% of their length are dropped, so background contigs cannot pick up an arbitrary lineage call.
Rows: `rows_background.jsonl`, manifest `manifest_rows_background.json`. The manifest's git SHA is
the commit before the `--background` code was committed; the code is in this PR. 0 tool errors.

| rung | separated breadth, spike only | with Berlin background | mosaic mass, spike / bg | switches per 10 kb, spike / bg |
|---|---|---|---|---|
| 2.0% | 0.00 [0.00-0.00] | 0.00 [0.00-0.00] | 0.98 / 0.98 | 12.9 / 12.5 |
| 4.5% | 0.52 [0.23-0.95] | 0.44 [0.23-0.79] | 0.34 / 0.38 | 2.6 / 2.6 |

**The background does not change the conclusion.** At 2.0% both conditions collapse completely.
At 4.5%, 9 of 15 runs give separated breadth within 0.01 of each other with and without the
background, 4 lose some and 2 gain some (mean change −0.08, SD 0.24). Nearly all of the drop is
one run: pair MG458312.1/PZ227553.1, seed 1, went from 0.96 to 0.05. That is not a collapse. Both
lineages are still assembled end to end, but each long contig has picked up one 300 nt window of
the other lineage at its end, and the strict all-windows-pure rule then counts the whole 8 kb as
mosaic. Switch density, which does not have that cliff, is identical (2.6 per 10 kb in both
conditions). So the background moves where MEGAHIT ends its contigs, and that matters only for
the strict readout. It does not move where the lineages collapse.

## Not established, and what would change it

- **Haploflow is not judged here.** At defaults it is mostly mosaic at every rung, including
  9.7%, and returns nothing at 5x. That is more likely a configuration problem than a verdict,
  and its `--strict` and `--error-rate` options were never swept. Do not quote it.
- **Spike only, checked on two rungs.** The real-background parity check below covers 2.0% and
  4.5% for MEGAHIT at 15x only. Other tools and rungs are untested with a background.
- **Equal abundance only.** A minor lineage at 1:4 or 1:10 is the case the stress test found
  hardest, and it is not covered.
- **One virus family**, 21 pairs, simulated reads with a single error profile. HS25 was chosen
  because under MSv3 PenguiN could not build one 1 kb contig from a single lineage at 60x, which
  would have made the benchmark unfair to it.
- **Labels are window-resolution.** A contig with one swapped 300 nt window counts as mosaic
  in the strict readout. Switch density is the measure that tells rare local swaps (metaSPAdes
  at 7.3%: 1.4 per 10 kb) apart from collapse (MEGAHIT at 2%: 12.9). Join positions are known
  only to ±150 nt.
- **PenguiN raw's separated breadth is a union over hundreds of contigs.** It is an upper bound
  on what a path-selection step could recover, not output a user could use.
