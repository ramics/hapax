Seeds: 5; 2,000 simulated pairs per genome per rung.

**Read the per-genome tier table, not the pooled one.** Of the spiked genomes, `NC_003630.1` pepper mild mottle virus is resident, 285x; `NC_001943.1` human astrovirus 1 is resident, 11.7x as a 6,823 nt contig; `NC_038311.1` human rhinovirus 1 is absent from this sample's assembly. A genome already in the sample is removed by the sample's own assembly for a reason that has nothing to do with detecting a novel virus, so only the holdout column (NC_038311.1) answers the question the spike is for.

## Spike survival through tier 1, by divergence

| divergence | NC_003630.1 | NC_001943.1 | NC_038311.1 | all |
|---|---|---|---|---|
| 0% | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) | 0.0 (0.0-0.0) |
| 5% | 0.5 (0.0-1.6) | 0.3 (0.0-0.9) | 1.0 (0.1-2.0) | 0.6 (0.0-2.0) |
| 10% | 13.9 (7.7-17.2) | 11.5 (6.4-19.2) | 12.5 (8.8-15.8) | 12.7 (6.4-19.2) |
| 20% | 85.0 (79.6-88.6) | 81.8 (74.8-85.5) | 79.8 (73.8-89.8) | 82.2 (73.8-89.8) |
| 30% | 99.7 (98.8-100.0) | 99.1 (98.2-100.0) | 99.6 (98.4-100.0) | 99.5 (98.2-100.0) |

Mean % of simulated pairs surviving every index the cascade subtracts (min-max across seeds).

## Which tier removes the spike, per genome

The column that answers the spike's question is the holdout's. `selfasm_identified` and `cross_sample` read 0.0 on it at every rung, because a virus absent from the sample cannot be in the sample's assembly.

| divergence | genome | rrna | selfasm_identified | selfasm_unidentified | cross_sample | k2_viral | k2_pluspf16 | prot_viral | prot_sprot |
|---|---|---|---|---|---|---|---|---|---|
| 0% | NC_003630.1 | 0.0 | 100.0 | 0.0 | 100.0 | 100.0 | 64.3 | 99.5 | 99.5 |
| 0% | NC_001943.1 | 0.0 | 96.4 | 0.0 | 99.4 | 100.0 | 59.6 | 99.5 | 99.4 |
| 0% | NC_038311.1 **(holdout)** | 0.0 | 0.0 | 0.0 | 0.0 | 100.0 | 68.0 | 93.6 | 93.6 |
| 5% | NC_003630.1 | 0.0 | 93.8 | 0.0 | 94.8 | 87.5 | 6.6 | 96.1 | 96.1 |
| 5% | NC_001943.1 | 0.0 | 61.1 | 0.0 | 69.7 | 85.8 | 5.2 | 96.4 | 96.6 |
| 5% | NC_038311.1 **(holdout)** | 0.0 | 0.0 | 0.0 | 0.0 | 86.4 | 5.9 | 91.3 | 91.3 |
| 10% | NC_003630.1 | 0.0 | 44.4 | 0.0 | 46.6 | 28.0 | 0.6 | 81.1 | 81.1 |
| 10% | NC_001943.1 | 0.0 | 20.9 | 0.0 | 26.4 | 28.8 | 0.5 | 82.8 | 83.0 |
| 10% | NC_038311.1 **(holdout)** | 0.0 | 0.0 | 0.0 | 0.0 | 27.3 | 0.1 | 84.8 | 84.9 |
| 20% | NC_003630.1 | 0.0 | 0.4 | 0.0 | 0.5 | 0.0 | 0.0 | 14.5 | 14.4 |
| 20% | NC_001943.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 17.8 | 18.1 |
| 20% | NC_038311.1 **(holdout)** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 20.0 | 20.2 |
| 30% | NC_003630.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.3 | 0.3 |
| 30% | NC_001943.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.7 | 0.8 |
| 30% | NC_038311.1 **(holdout)** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.4 | 0.4 |

## Which tier removes the spike, pooled over genomes (**inflated — kept for comparison only**)

Two of the three genomes are resident, so this table credits the sample's own assembly with removals that a novel virus would not have given it. It is the table this repo quoted before the split above existed.


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
