# Tables: residue structure

Generated from `structure_s1.json`; see `REPORT.md` for what each one means.

## Separability — residue against subtracted, 5-fold CV, mean (min–max) over 5 seeds

| space | AUC |
|---|---|
| all eleven axes | 0.9053 (0.9042–0.9077) |
| best single axis (`log_cov`) | 0.8921 (0.8911–0.8934) |
| five coverage axes | 0.9043 (0.9030–0.9068) |
| six composition and quality axes | 0.6864 (0.6783–0.6914) |

## Per-axis AUC, seed 1

| axis | AUC | | axis | AUC |
|---|---|---|---|---|
| log_cov | 0.892 | | entropy3 | 0.683 |
| cov_step | 0.780 | | UpA_oe | 0.618 |
| cov_spread | 0.767 | | gc | 0.525 |
| mate_disagree | 0.762 | | CpG_oe | 0.521 |
| zero_frac | 0.617 | | homopolymer | 0.519 |
| | | | mean_error_prob | 0.507 |

## Gap statistic against a covariance-matched reference, seed 1

| k | 1 | 2 | 3 | 4 | 5 | 6 | 10 | 20 |
|---|---|---|---|---|---|---|---|---|
| gap | −0.0015 | −0.0017 | 0.041 | 0.070 | 0.107 | 0.133 | ~0.22 | ~0.33 |

Tibshirani rule: k = 1 in all five seeds. s_k ≈ 0.0025 throughout.

## Axis means, block 1

| axis | residue (22,839) | subtracted (76,906) | shuffled null (19,858) | planted 1e-5 |
|---|---|---|---|---|
| log_cov | 1.640 | 6.170 | 0.000 | 0.964 |
| cov_spread | 1.132 | 2.861 | 0.000 | 1.150 |
| zero_frac | 0.276 | 0.050 | 1.000 | 0.186 |
| mate_disagree | 0.506 | 1.418 | 0.001 | 0.704 |
| cov_step | 1.136 | 3.428 | 0.000 | 0.994 |
| gc | 0.504 | 0.497 | 0.504 | 0.416 |
| entropy3 | 0.939 | 0.952 | 0.941 | 0.948 |
| CpG_oe | 0.953 | 0.941 | 0.955 | 0.716 |
| UpA_oe | 0.775 | 0.897 | 0.777 | 0.946 |

## 10-NN purity of a planted genome, 5 simulations per abundance

| abundance | pairs | all axes but quality | coverage only | composition only | random residue subset | chance |
|---|---|---|---|---|---|---|
| 1e-3 | 8,731 | 0.9959 | 0.9960 | 0.7903 | 0.3960 | 0.277 |
| 1e-4 | 873 | 0.8940 | 0.9179 | 0.2920 | 0.0578 | 0.037 |
| 1e-5 | 87 | 0.1041 | 0.0607 | 0.0395 | 0.0248 | 0.004 |
| 1e-6 | 9 | 0.0267 | 0.0266 | 0.0089 | 0.0178 | 0.0004 |

The control column is the same-sized random residue subset in the matching space
(`all_but_quality`); per-space control values are in the JSON.

## What cov_step is: block 1 residue pairs with cov_step >= 2 (n = 4,290 of 22,975)

| readout | value | same readout, all residue pairs |
|---|---|---|
| no unwitnessed 31-mer at all | 45.0% | — |
| step survives dropping unwitnessed k-mers (>=1) | 72.0% | — |
| longest zero-run 25–40 k-mers (one-error signature) | 20.4% | 7.1% |
| mean longest zero-run | 16.3 | 21.8 |

Rank correlations over all block-1 residue pairs: cov_step with the same step computed on
witnessed k-mers only 0.804, with longest zero-run −0.128, with zero_frac −0.029, with
log_cov 0.592.
