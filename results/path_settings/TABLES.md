## 1. What each setting is, and which way it was expected to move

| arm | flags added to the baseline | stated direction |
|---|---|---|
| `baseline` | `(none)` | the committed setting; everything else is read against it |
| `baseline_repeat` | `(none)` | same setting, different thread count: the same-instrument difference, which bounds what an arm has to beat to mean anything |
| `keep_local` | `--low-local-ratio 0.05 --disconnect-ratio 0.02` | keep low-coverage local branches that the defaults remove |
| `keep_thin` | `--prune-level 0 --prune-depth 1` | keep low-depth unitigs. Expect more low-depth contigs and more junk |
| `lax_bubble` | `--merge-level 20,0.90` | merge bubbles down to 90% similarity. THE ABSORPTION ARM: expect longer contigs, more mosaics, and a residue that shrinks wrongly |
| `meta_sensitive` | `--presets meta-sensitive` | the preset the residue assembly used, never compared to the default until now. NOT a path arm: it also sets --min-count 1, a lever this project already measured as dead, so it is reported separately and no path conclusion rests on it |
| `no_bubble` | `--bubble-level 0` | no bubble merging at all. The extreme of keeping paths apart |
| `no_rescue` | `--no-mercy --no-local` | remove both low-coverage rescue mechanisms. The least-path-rescue end |
| `prune_hard` | `--prune-level 3` | prune harder. Expect a larger residue and little change deep down |
| `strict_bubble` | `--merge-level 20,0.98` | merge only near-identical bubbles. Expect more, shorter contigs and less lineage fusion; breadth may fall |

## 2. The assemblies

| arm | contigs | bases | N50 | max | >=1 kb | >=3 kb | minutes | threads |
|---|---|---|---|---|---|---|---|---|
| `baseline` | 35,328 | 22,037,571 | 649 | 20,695 | 3,968 | 322 | 11.8 | 4 |
| `baseline_repeat` | 35,330 | 22,030,733 | 649 | 20,695 | 3,953 | 319 | 21.1 | 2 |
| `keep_local` | 36,345 | 22,209,442 | 634 | 21,798 | 3,878 | 259 | 11.9 | 4 |
| `keep_thin` | 35,510 | 22,019,623 | 645 | 21,798 | 3,902 | 300 | 10.3 | 4 |
| `lax_bubble` | 35,204 | 22,012,860 | 651 | 20,695 | 3,954 | 332 | 12.1 | 4 |
| `meta_sensitive` | 53,906 | 29,809,319 | 547 | 15,436 | 4,095 | 304 | 36.8 | 4 |
| `no_bubble` | 35,715 | 21,433,865 | 616 | 20,695 | 3,577 | 233 | 12.0 | 4 |
| `no_rescue` | 23,684 | 14,132,183 | 614 | 14,854 | 2,363 | 189 | 3.9 | 4 |
| `prune_hard` | 26,975 | 18,283,273 | 737 | 20,695 | 3,772 | 298 | 10.8 | 4 |
| `strict_bubble` | 35,567 | 21,954,518 | 639 | 20,695 | 3,892 | 296 | 11.7 | 4 |

## 3. Ground truth: 25 genomes the depositors assembled from this same water

| arm | >=90% recovered | >=50% | mean breadth | median purity | p10 purity | contig joins | false joins (control) |
|---|---|---|---|---|---|---|---|
| `baseline` | 20/25 | 21/25 | 0.810 | 1.000 | 0.983 | 0 | 0/53 |
| `baseline_repeat` | 20/25 | 21/25 | 0.808 | 1.000 | 0.943 | 0 | 0/53 |
| `keep_local` | 16/25 | 21/25 | 0.785 | 1.000 | 0.981 | 0 | 0/53 |
| `keep_thin` | 20/25 | 21/25 | 0.811 | 1.000 | 0.995 | 0 | 0/53 |
| `lax_bubble` | 20/25 | 21/25 | 0.809 | 1.000 | 0.983 | 0 | 0/53 |
| `meta_sensitive` | 20/25 | 21/25 | 0.819 | 1.000 | 0.989 | 0 | 0/53 |
| `no_bubble` | 18/25 | 21/25 | 0.804 | 1.000 | 0.985 | 0 | 0/53 |
| `no_rescue` | 15/25 | 20/25 | 0.747 | 1.000 | 1.000 | 0 | 0/53 |
| `prune_hard` | 20/25 | 21/25 | 0.807 | 1.000 | 0.994 | 0 | 0/53 |
| `strict_bubble` | 20/25 | 20/25 | 0.807 | 1.000 | 0.986 | 0 | 0/53 |

### 3a. The trade, which the aggregate hides

Genomes whose breadth moves by at least 10 points against the `baseline` arm. The `baseline_repeat` row is the null: it is the same setting, so whatever it shows is what a difference has to beat.

| arm | gained >=10 pts | lost >=10 pts | net | recovered but embedded (breadth >=0.9, purity <0.9) |
|---|---|---|---|---|
| `baseline` | 0 | 0 | +0 | 0 |
| `baseline_repeat` | 0 | 0 | +0 | 0 |
| `keep_local` | 1 | 4 | -3 | 0 |
| `keep_thin` | 0 | 0 | +0 | 0 |
| `lax_bubble` | 0 | 0 | +0 | 0 |
| `meta_sensitive` | 1 | 0 | +1 | 0 |
| `no_bubble` | 1 | 2 | -1 | 0 |
| `no_rescue` | 0 | 6 | -6 | 0 |
| `prune_hard` | 0 | 0 | +0 | 0 |
| `strict_bubble` | 0 | 0 | +0 | 0 |

### 3b. Per genome, baseline against every other arm

| genome | depth | divergence % | baseline: breadth/purity | baseline_repeat: breadth/purity | keep_local: breadth/purity | keep_thin: breadth/purity | lax_bubble: breadth/purity | meta_sensitive: breadth/purity | no_bubble: breadth/purity | no_rescue: breadth/purity | prune_hard: breadth/purity | strict_bubble: breadth/purity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| OQ722330.1 | 700.7x | 0.08 | 0.99/1.00 | 0.99/1.00 | 1.00/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/0.99 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 |
| OQ722321.1 | 499.0x | 0.07 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/0.99 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 |
| OQ722337.1 | 452.4x | 30.28 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 | 0.96/1.00 |
| OQ722331.1 | 354.7x | 9.84 | 0.99/0.99 | 0.99/0.99 | 0.99/1.00 | 0.99/0.99 | 0.99/0.99 | 0.99/1.00 | 0.91/0.98 | 0.94/1.00 | 0.99/0.99 | 0.99/0.99 |
| OQ722319.1 | 282.8x | 9.94 | 0.99/0.99 | 0.99/0.99 | 0.99/1.00 | 0.99/0.99 | 0.99/0.99 | 0.99/1.00 | 0.91/0.98 | 0.94/1.00 | 0.99/0.99 | 0.99/0.99 |
| OQ722336.1 | 181.8x | 7.03 | 0.97/1.00 | 0.97/0.99 | 0.81/0.92 | 1.00/1.00 | 0.97/1.00 | 1.00/1.00 | 0.94/0.99 | 0.92/1.00 | 0.97/1.00 | 0.97/1.00 |
| OQ722325.1 | 142.3x | 15.66 | 0.98/0.98 | 0.98/0.92 | 0.63/1.00 | 0.98/0.99 | 0.98/0.98 | 0.98/0.92 | 0.98/0.99 | 0.98/0.99 | 0.98/0.99 | 0.98/0.99 |
| OQ722335.1 | 63.4x | 17.09 | 1.00/0.99 | 1.00/0.99 | 1.00/0.98 | 1.00/0.99 | 1.00/0.99 | 1.00/0.99 | 1.00/0.99 | 1.00/1.00 | 1.00/0.99 | 1.00/0.99 |
| OQ722324.1 | 58.0x | 17.09 | 0.99/0.99 | 0.99/0.99 | 0.99/0.98 | 0.99/0.99 | 0.99/0.99 | 0.99/0.99 | 0.99/0.99 | 0.99/1.00 | 0.99/0.99 | 0.99/0.99 |
| OQ722334.1 | 44.7x | 1.9 | 1.00/1.00 | 0.98/1.00 | 0.96/0.99 | 1.00/1.00 | 1.00/1.00 | 0.99/1.00 | 0.88/1.00 | 0.98/1.00 | 1.00/1.00 | 0.99/0.99 |
| OQ722328.1 | 43.9x | 4.12 | 0.98/1.00 | 0.97/1.00 | 0.76/1.00 | 0.98/1.00 | 0.98/1.00 | 0.98/1.00 | 0.95/1.00 | 0.91/1.00 | 0.98/1.00 | 0.98/1.00 |
| OQ722317.1 | 42.2x | 4.25 | 0.98/1.00 | 0.97/1.00 | 0.75/1.00 | 0.97/1.00 | 0.98/1.00 | 0.98/1.00 | 0.92/1.00 | 0.90/1.00 | 0.97/1.00 | 0.98/1.00 |
| OQ722332.1 | 33.4x | 1.1 | 0.54/1.00 | 0.52/1.00 | 0.84/0.95 | 0.54/1.00 | 0.51/1.00 | 0.63/1.00 | 0.57/1.00 | 0.14/1.00 | 0.54/1.00 | 0.45/0.95 |
| OQ722338.1 | 32.9x | 27.19 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 0.87/1.00 | 0.89/1.00 | 1.00/1.00 | 1.00/1.00 |
| OQ722327.1 | 27.9x | 24.48 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 |
| OQ722318.1 | 21.0x | 2.69 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.86/1.00 | 0.99/1.00 | 0.99/1.00 |
| OQ722315.1 | 20.7x | 27.18 | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.33/0.99 | 0.00/- | 0.00/- | 0.00/- |
| OQ722323.1 | 16.6x | 19.42 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 |
| OQ722333.1 | 12.4x | 0.48 | 0.95/0.94 | 0.95/0.94 | 0.97/0.98 | 0.95/1.00 | 0.95/0.94 | 0.92/1.00 | 0.95/1.00 | 0.73/1.00 | 0.91/1.00 | 0.96/1.00 |
| OQ722320.1 | 10.6x | 0.47 | 0.97/0.94 | 0.97/0.94 | 1.00/0.98 | 0.97/1.00 | 0.97/0.94 | 0.93/1.00 | 0.96/1.00 | 0.75/1.00 | 0.93/1.00 | 0.98/1.00 |
| OQ722316.1 | 8.9x | 30.05 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.99/1.00 | 0.91/1.00 | 0.99/1.00 | 0.99/1.00 |
| OQ722329.1 | 3.7x | 2.72 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 0.87/1.00 | 1.00/1.00 | 1.00/1.00 |
| OQ722326.1 | 3.4x | 16.05 | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- |
| OQ722314.1 | 0.7x | 26.08 | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.04/1.00 | 0.00/- | 0.00/- | 0.00/- | 0.00/- |
| OQ722322.1 | 0.4x | 30.38 | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.00/- | 0.11/1.00 | 0.00/- | 0.00/- | 0.00/- | 0.00/- |

## 3c. What each arm leaves unidentifiable, which is tier 2

The absorption question the truth panel cannot answer. A setting that fuses novel sequence into abundant known contigs shows up here as identified bases rising while tier-2 mass falls, with no change in the panel. Identification is against indexes rolled back to 2022-05-10, the day this water was collected, at the committed operating point (Kraken2 confidence 0.1, protein bitscore 50).

| arm | identified contigs % | identified bases % | tier-2 contigs | tier-2 bases | tier-2 bases >=500 nt and >=2x |
|---|---|---|---|---|---|
| `baseline` | 63.11 | 70.32 | 13,034 | 6,541,620 | 31.0% / 96.1% |
| `baseline_repeat` | 63.08 | 70.29 | 13,043 | 6,544,532 | 30.9% / 96.1% |
| `keep_local` | 63.71 | 70.24 | 13,188 | 6,608,745 | 31.0% / 96.1% |
| `keep_thin` | 63.18 | 70.22 | 13,075 | 6,557,602 | 31.0% / 96.1% |
| `lax_bubble` | 63.06 | 70.33 | 13,006 | 6,531,285 | 31.0% / 96.2% |
| `meta_sensitive` | 59.46 | 65.90 | 21,853 | 10,166,306 | 24.8% / 73.5% |
| `no_bubble` | 63.34 | 69.71 | 13,093 | 6,491,269 | 30.0% / 94.9% |
| `no_rescue` | 65.31 | 71.39 | 8,216 | 4,043,909 | 28.7% / 99.7% |
| `prune_hard` | 64.81 | 72.15 | 9,492 | 5,092,576 | 36.8% / 100.0% |
| `strict_bubble` | 63.15 | 70.15 | 13,106 | 6,553,184 | 30.7% / 95.8% |

## 4. Residue, as an output of each arm

NOTE: these fractions are not comparable to the 77.10% of `results/subtract_known/`. The two 8 GB Kraken2 tiers are not staged here (fixed across arms, marginal 0.27 and 0.00 points there), and this is a different library. They are comparable across the arms below and nowhere else.

| arm | removed % | residue % (mean) | residue min-max | tier 2 % (assembled, unidentifiable) | tier 3 % (maps nothing) | seeds |
|---|---|---|---|---|---|---|
| `baseline` | 64.58 | 35.42 | 34.89-35.83 | 15.38 | 20.04 | 5 |
| `baseline_repeat` | 64.59 | 35.41 | 34.89-35.81 | 15.68 | 19.73 | 5 |
| `keep_local` | 64.31 | 35.69 | 35.23-36.08 | 13.89 | 21.80 | 5 |
| `keep_thin` | 65.27 | 34.73 | 34.21-35.14 | 13.11 | 21.62 | 5 |
| `lax_bubble` | 65.91 | 34.09 | 33.60-34.49 | 14.14 | 19.95 | 5 |
| `meta_sensitive` | 66.16 | 33.84 | 33.30-34.19 | 14.58 | 19.26 | 5 |
| `no_bubble` | 62.97 | 37.03 | 36.47-37.42 | 16.06 | 20.97 | 5 |
| `no_rescue` | 57.37 | 42.63 | 41.91-43.01 | 14.35 | 28.28 | 5 |
| `prune_hard` | 63.80 | 36.20 | 35.68-36.59 | 13.59 | 22.61 | 5 |
| `strict_bubble` | 64.13 | 35.87 | 35.34-36.25 | 15.41 | 20.46 | 5 |

### 4b. Residue against baseline, PAIRED BY SEED

Five seeds are five disjoint blocks of one read pool, so the same seed is the same reads and the paired difference is the comparison with the subsampling noise taken out. The `baseline_repeat` row is the null.

| arm | residue change (pts, mean) | per-seed range | tier 2 change | tier 3 change | seeds paired |
|---|---|---|---|---|---|
| `baseline` | +0.00 | +0.00 to +0.00 | +0.00 | +0.00 | 5 |
| `baseline_repeat` | -0.00 | -0.01 to +0.01 | +0.31 | -0.31 | 5 |
| `keep_local` | +0.27 | +0.20 to +0.33 | -1.49 | +1.76 | 5 |
| `keep_thin` | -0.69 | -0.71 to -0.66 | -2.26 | +1.58 | 5 |
| `lax_bubble` | -1.33 | -1.36 to -1.30 | -1.24 | -0.09 | 5 |
| `meta_sensitive` | -1.57 | -1.63 to -1.52 | -0.80 | -0.78 | 5 |
| `no_bubble` | +1.61 | +1.58 to +1.66 | +0.68 | +0.93 | 5 |
| `no_rescue` | +7.22 | +7.01 to +7.37 | -1.03 | +8.24 | 5 |
| `prune_hard` | +0.78 | +0.76 to +0.80 | -1.79 | +2.57 | 5 |
| `strict_bubble` | +0.45 | +0.42 to +0.47 | +0.03 | +0.42 | 5 |

## 4c. Absorption at the contig level

How much of the BASELINE's unidentifiable pile now sits inside a contig the cascade calls identified, per arm. `baseline` is 0 by construction (the two piles are disjoint) and `baseline_repeat` is the null: the same setting at a different thread count, so it is what this reads when nothing changed.

| arm | contigs >=90% inside an identified contig | of bases | coverage-weighted | of which the single largest contig |
|---|---|---|---|---|
| `baseline` | 0 of 13,034 | 0.00% | 0.00% | - |
| `baseline_repeat` | 13 of 13,034 | 0.12% | 0.18% | 54% |
| `lax_bubble` | 38 of 13,034 | 0.40% | 0.26% | 29% |
| `strict_bubble` | 44 of 13,034 | 0.35% | 0.25% | 10% |
| `no_bubble` | 71 of 13,034 | 0.67% | 0.67% | 14% |
| `keep_thin` | 24 of 13,034 | 0.26% | 7.61% | 95% |
| `keep_local` | 76 of 13,034 | 0.69% | 8.25% | 79% |
| `prune_hard` | 84 of 13,034 | 0.83% | 0.47% | 20% |
| `no_rescue` | 28 of 13,034 | 0.40% | 0.31% | 22% |
| `meta_sensitive` | 720 of 13,034 | 5.42% | 10.63% | 73% |

## 5. Provenance

- script `experiments/path_settings.py` at `d4a01adccdea`, git `6668bf4`
- MEGAHIT v1.2.9, minimap2 2.28-r1209
- reads sha256 `a199c914af6e`
