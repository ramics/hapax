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

## Mock metagenome components

Downloaded 2026-11-19 from NCBI nuccore. Background: `NC_000913.3` (E. coli
K-12 MG1655), `NC_007795.1` (S. aureus NCTC 8325), `NC_009614.1` (Phocaeicola
vulgatus ATCC 8482). Phages: `NC_001416.1` (lambda), `NC_002167.1` (HK97),
`NC_000866.4` (T4), `NC_001422.1` (phiX174), `NC_024711.1` (Carjivirus
communis / crAssphage). Viruses: `NC_045512.2` (SARS-CoV-2), `MN996532.2`
(RaTG13), `NC_014470.1` (bat CoV BM48-31), `NC_003977.2` (HBV genotype D),
`X70185.1` (HBV genotype A), `NC_001802.1` (HIV-1 subtype B), `AY772699.1`
(HIV-1 subtype C, 04ZASK146B1), `NC_001526.4` (HPV16), `NC_001525.1` (HPV11),
`NC_006273.2` (HHV-5 Merlin), `NC_001716.2` (HHV-7), `NC_002076.2` (TTV1),
`NC_039477.1` (norovirus GII), `NC_001542.1` (rabies), `NC_075022.1` (VEEV).

Measured whole-genome divergence of the known/relative pairs (edit distance
over the longer genome): SARS-CoV-2/RaTG13 0.040, HBV D/A 0.105, HIV-1 B/C
0.163, SARS-CoV-2/BM48-31 0.240, lambda/HK97 0.458, HPV16/HPV11 0.506,
HHV-5/HHV-7 0.546.
