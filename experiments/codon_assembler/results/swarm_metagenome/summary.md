# Multi-swarm metagenome: 5 seed(s)

Eight real within-host swarms (see experiments/fetch_swarm_panels.py) on a
bacterial and phage background. Contiguity and chimaeras from the same runs.

| tool | swarms found /8 | contigs per found swarm | completeness | cross-swarm chimaeras | all chimaeras | haplotypes resolved | total contigs | variant recall | variant precision | wall (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hapax-codon | 4.80 ± 0.45 | 1.25 ± 0.09 | 0.57 ± 0.05 | 3.60 ± 0.89 | 11.80 ± 1.64 | 0.00 ± 0.00 | 7406.00 ± 40.88 | 0.47 ± 0.05 | 0.55 ± 0.06 | 283.95 ± 105.03 | 2.29 ± 0.05 |
| hapax-codon-t1groups | 5.00 ± 0.00 | 1.00 ± 0.00 | 0.62 ± 0.00 | 3.80 ± 0.45 | 5.80 ± 0.84 | 1.00 ± 0.00 | 7555.60 ± 27.13 | 0.48 ± 0.03 | 0.52 ± 0.03 | 284.16 ± 109.31 | 2.31 ± 0.13 |
| hapax-nucleotide | 4.80 ± 0.45 | 1.25 ± 0.09 | 0.57 ± 0.05 | 3.60 ± 0.89 | 7.00 ± 1.73 | 0.00 ± 0.00 | 7477.60 ± 37.38 | 0.47 ± 0.05 | 0.55 ± 0.06 | 280.96 ± 107.33 | 2.33 ± 0.13 |
| haploflow | 8.00 ± 0.00 | 3.95 ± 0.78 | 0.87 ± 0.03 | 12.40 ± 3.13 | 12.60 ± 2.88 | 1.00 ± 0.00 | 106.60 ± 9.76 | 0.65 ± 0.06 | 0.75 ± 0.04 | 126.18 ± 4.14 | 3.46 ± 0.00 |
| megahit | 7.60 ± 0.89 | 3.28 ± 0.43 | 0.73 ± 0.04 | 1.80 ± 0.84 | 1.80 ± 0.84 | 1.00 ± 0.00 | 4170.80 ± 41.41 | 0.53 ± 0.08 | 0.86 ± 0.02 | 300.53 ± 120.53 | 0.39 ± 0.00 |
| metaspades | 7.00 ± 0.00 | 1.63 ± 0.22 | 0.81 ± 0.04 | 2.00 ± 0.71 | 2.40 ± 0.55 | 1.40 ± 0.55 | 3369.00 ± 50.02 | 0.56 ± 0.04 | 0.82 ± 0.02 | 373.80 ± 13.25 | 1.31 ± 0.00 |
