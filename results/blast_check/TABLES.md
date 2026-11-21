# BLAST ground truth against the public databases

5 seed(s): 1, 2, 3, 4, 5. 60 sequences per class per seed, labelled by the same references and thresholds as the cascade.

**Batches lost to malformed replies and correctly dropped**, per seed: seed 1: 1. Those queries are dropped rather than scored as misses - scoring them as misses is how an earlier version of this table reported 50% any-hit for a class that sits at 100%.

## blastn against nt

| class | what it is | any hit | strong hit | median identity | median query coverage |
|---|---|---|---|---|---|
| `rrna` | removed by SILVA | 100.0 (100.0-100.0) | **89.3 (86.7-93.3)** | 99.8 | 1.0 |
| `public_nucleotide` | removed by Kraken2, not by SILVA | 100.0 (100.0-100.0) | **94.7 (93.3-96.7)** | 100.0 | 1.0 |
| `selfasm_identified` | removed only by this sample's own identified contigs | 89.0 (85.0-93.3) | **26.7 (16.7-33.3)** | 93.0 | 0.6 |
| `residue_assembled` | survived, maps to an unidentifiable contig of this sample (tier 2) | 82.3 (78.3-85.0) | **8.7 (5.0-11.7)** | 90.3 | 0.4 |
| `residue_missing` | survived, maps to no contig at all (tier 3) | 89.0 (86.7-91.7) | **5.7 (3.3-6.7)** | 88.2 | 0.4 |
| `tier2_contigs` | unidentifiable contigs >=500 nt, >=2x coverage | 98.8 (95.6-100.0) | **5.2 (3.3-10.0)** | 87.1 | 0.1 |

## How big each class is in the library

The sampled 60 per class says nothing about abundance, so the proportions come from the labelled block itself.

| class | share of pairs |
|---|---|
| `rrna` | 38.0 (37.6-38.3) |
| `public_nucleotide` | 2.7 (2.5-2.9) |
| `selfasm_identified` | 34.3 (33.9-34.7) |
| `residue_assembled` | 8.8 (8.6-9.1) |
| `residue_missing` | 16.2 (16.1-16.4) |

## What the table says

- **Any hit is not a test.** Read the strong-hit column instead; the two disagree by tens of points on exactly the classes the architecture cares about.
- The two residue classes separate from the removed ones continuously, in identity and coverage together, which is the same shape the contig anatomy and the coverage-shape work both arrived at independently.
- 60 per class per seed gives percentage-level resolution and no better. Every figure here is a proportion with an interval in the per-seed files; the range across seeds is in brackets above.
