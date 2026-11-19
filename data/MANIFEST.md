# Data manifest

`data/` is gitignored except for this file. The development sandbox cannot reach
NCBI, so nothing below has been downloaded yet. Record accession, source URL and
download date when a file is fetched.

| Purpose | Accession / source | Downloaded | Notes |
|---|---|---|---|
| HIV-1 reference (coordinates for *env*, *gag/pol*) | HXB2, GenBank K03455 | not yet | thesis EvolveAGene seed was HXB2 *env* |
| HIV-1 *env* population data with ground truth | Imo's curated alignments (thesis ch. 2/4) | not yet | ask Imo for location |
| HBV reference (dual-coding stress test) | NC_003977 | not yet | ~2/3 of the genome is dual-coding |
| Public wastewater metagenome | to choose (SRA) | not yet | pick after HIV results |

## Within-host env panels for the incumbent stress test

Downloaded 2026-11-19 from NCBI nuccore via eutils, by querying
`HIV-1[orgn] AND env[gene] AND 2400:2700[SLEN] AND single genome amplification`
and grouping by host isolate. These are published single-genome-amplification
(SGA) env sequence sets from the Zambia transmission cohort, so each file is a
sample of the haplotypes that co-existed in one person.

| file | host | sequences | mean pairwise divergence | max | length range (nt) |
|---|---|---|---|---|---|
| `env/ZM231F.fa` | ZM231F | 48 | 0.0028 | 0.0157 | 2513-2545 |
| `env/ZM247F.fa` | ZM247F | 65 | 0.0300 | 0.1434 | 2412-2604 |
| `env/ZM206F.fa` | ZM206F | 35 | 0.0030 | 0.0109 | 2538-2560 |
| `env/ZM214M.fa` | ZM214M | 35 | 0.0082 | 0.0164 | 2544-2568 |
| `env/hxb2_env.fa` | HXB2 reference | 1 | - | - | 2571 (K03455:6225-8795) |

Two ZM247F sequences carry IUPAC ambiguity codes and are dropped by
`bench_panel.build_panel` rather than resolved.
