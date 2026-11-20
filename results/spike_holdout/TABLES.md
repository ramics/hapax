Seeds: 5; 2,000 simulated pairs per genome per rung; genomes: NC_003630.1, NC_001943.1, NC_038311.1.

## Spike survival through tier 1, by divergence

| divergence | NC_003630.1 | NC_001943.1 | NC_038311.1 | all |
|---|---|---|---|---|
| 0% | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) |
| 5% | 0.5 (0.0-1.6) | 0.3 (0.0-0.9) | 1.0 (0.1-2.0) | 0.6 (0.0-2.0) |
| 10% | 13.9 (7.7-17.2) | 11.5 (6.4-19.2) | 12.5 (8.8-15.8) | 12.7 (6.4-19.2) |
| 20% | 85.0 (79.6-88.6) | 81.8 (74.8-85.5) | 79.8 (73.8-89.8) | 82.2 (73.8-89.8) |
| 30% | 99.7 (98.8-100.0) | 99.1 (98.2-100.0) | 99.6 (98.4-100.0) | 99.5 (98.2-100.0) |

Mean % of simulated pairs surviving every index the cascade subtracts (min-max across seeds).

## Which tier removes the spike

| divergence | rrna | selfasm_identified | selfasm_unidentified | cross_sample | k2_viral | k2_pluspf16 | prot_viral | prot_sprot |
|---|---|---|---|---|---|---|---|---|
| 0% | 0.0 | 65.5 | 0.0 | 66.5 | 100.0 | 64.0 | 97.5 | 97.5 |
| 5% | 0.0 | 51.6 | 0.0 | 54.8 | 86.6 | 5.9 | 94.6 | 94.7 |
| 10% | 0.0 | 21.8 | 0.0 | 24.3 | 28.1 | 0.4 | 82.9 | 83.0 |
| 20% | 0.0 | 0.2 | 0.0 | 0.2 | 0.0 | 0.0 | 17.4 | 17.6 |
| 30% | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.5 | 0.5 |

## At the real abundance

- A virus at 1 in a million contributes about 8.7 read pairs in this 8,731,224-pair run.
- Per-read survival above is measured at 2,000 pairs per rung, because a simulation at the real rate cannot resolve the quantity; the sentence above is arithmetic on top of it.
- Library recurrence is measured separately, with exactly the number of spike pairs 1-in-a-million implies, in results/spike_abundance/.
