Seeds: 5; 200,000 pairs per library per seed; script `e5051a199e16`, commit `ae0bfd1a925e`. 3 of 5 seeds ran from a dirty tree, so their commit names HEAD rather than the code; the script hash above is identical across all 5 and is what pins them.

Mate 2 trimmed by [14] bases — the measured prefix mode. It moves the Kraken2 row by about 1.7 points between 14 and 16; see `../nontemplate_prefix/`.

| | dna (DNA_220616_20-22hmix_cer_BIOMcombi) | rna (RNA_WW_220616_20-22hmix) | rna_2021 (RNA_WW_211013_20-22hmix (calibration)) |
|---|---|---|---|
| rRNA, ≥40 aligned bases against SILVA | 0.01 (0.00-0.01) | 67.19 (67.05-67.38) | 38.11 (38.02-38.28) |
| identified, Kraken2 PlusPF-16 at confidence 0.1 | 0.88 (0.86-0.92) | 38.93 (38.68-39.19) | 21.71 (21.54-21.86) |
| identified, RefSeq viral proteins, bitscore ≥50 | 7.32 (7.02-7.47) | 6.41 (6.09-6.56) | 11.91 (11.43-12.33) |
| any hit, RefSeq viral proteins | 66.70 (66.55-66.84) | 54.02 (53.44-54.55) | 65.58 (65.03-66.06) |
| identified, Swiss-Prot, bitscore ≥50 | 0.89 (0.81-0.98) | 5.96 (5.64-6.26) | 6.49 (6.12-6.84) |
| any hit, Swiss-Prot | 58.34 (57.37-58.86) | 63.92 (63.59-64.53) | 68.60 (68.12-69.06) |
| mean GC | 0.385 (0.384-0.385) | 0.508 (0.508-0.508) | 0.498 (0.498-0.499) |
| duplicate pairs (first 50 bases of both mates) | 21.81 (19.96-24.80) | 25.02 (23.57-25.99) | 24.21 (22.09-25.74) |
| modal length R1/R2 | 107/107 | 107/107 | 107/107 |

Mean across seeds (min-max). Per-seed Wilson intervals are in the JSONs; at this depth they are inside ±0.3 points of every figure, so the spread above is seed-to-seed variation and not sampling.

## Viral against cellular protein identification

| library | viral | Swiss-Prot | ratio |
|---|---|---|---|
| dna | 7.32 (7.02-7.47) | 0.89 (0.81-0.98) | 8.277 (7.520-9.166) |
| rna | 6.41 (6.09-6.56) | 5.96 (5.64-6.26) | 1.075 (1.048-1.105) |
| rna_2021 | 11.91 (11.43-12.33) | 6.49 (6.12-6.84) | 1.836 (1.730-1.918) |

A library whose identifiable material is viral reads high here; the RNA libraries sit near 1.

## Do the two nucleic acids see the same organisms?

| | mean (min-max) |
|---|---|
| taxa called in both | 84.400 (76.000-95.000) |
| Jaccard over called taxa | 0.121 (0.109-0.137) |
| DNA classified reads on shared taxa | 64.86 (61.87-69.49) |
| RNA classified reads on shared taxa | 85.67 (79.67-87.59) |

The DNA library has few classified reads, so its taxon count is as much a statement about depth as about biology; the read-weighted rows are the defensible half.

**Confound, carried from the run:** the RNA member is 2022-06-16 and the project's calibration library is 2021-10-13, so DNA-vs-RNA here is clean but either against 2021 mixes nucleic acid with eight months

