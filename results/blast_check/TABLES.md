# BLAST ground truth against the public databases

3 seed(s): 3, 4, 5. 60 sequences per class per seed, labelled by the same references and thresholds as the cascade.

**2 seed(s) excluded as pre-fix** (blast_check_s1.json (seed 1), blast_check_s2.json (seed 2)): written before failed batches were excluded from the denominator, so their misses and their lost queries are the same number and cannot be separated now. Queued to re-run rather than averaged in.

**Batches lost to malformed replies and correctly dropped**, per seed: seed 3: 2, seed 4: 1. Those queries are dropped rather than scored as misses - scoring them as misses is how an earlier version of this table reported 50% any-hit for a class that sits at 100%.

## blastn against nt

| class | what it is | any hit | strong hit | median identity | median query coverage |
|---|---|---|---|---|---|
| `rrna` | removed by SILVA | 100.0 (100.0-100.0) | **90.0 (86.7-93.3)** | 100.0 | 1.0 |
| `public_nucleotide` | removed by Kraken2, not by SILVA | 100.0 (100.0-100.0) | **94.4 (93.3-95.0)** | 100.0 | 1.0 |
| `selfasm_identified` | removed only by this sample's own identified contigs | 89.1 (85.0-95.6) | **25.2 (17.8-33.3)** | 92.4 | 0.6 |
| `residue_assembled` | survived, maps to an unidentifiable contig of this sample (tier 2) | 81.1 (78.3-85.0) | **9.3 (6.7-11.1)** | 89.4 | 0.4 |
| `residue_missing` | survived, maps to no contig at all (tier 3) | 90.0 (88.3-91.7) | **5.6 (3.3-6.7)** | 87.5 | 0.4 |
| `tier2_contigs` | unidentifiable contigs >=500 nt, >=2x coverage | 100.0 (100.0-100.0) | **3.9 (3.3-5.0)** | 87.1 | 0.1 |

## How big each class is in the library

The sampled 60 per class says nothing about abundance, so the proportions come from the labelled block itself.

| class | share of pairs |
|---|---|
| `rrna` | 38.1 (38.0-38.3) |
| `public_nucleotide` | 2.8 (2.7-2.9) |
| `selfasm_identified` | 34.1 (33.9-34.3) |
| `residue_assembled` | 8.7 (8.6-8.7) |
| `residue_missing` | 16.3 (16.2-16.4) |

## What the table says

- **Any hit is not a test.** Read the strong-hit column instead; the two disagree by tens of points on exactly the classes the architecture cares about.
- The two residue classes separate from the removed ones continuously, in identity and coverage together, which is the same shape the contig anatomy and the coverage-shape work both arrived at independently.
- 60 per class per seed gives percentage-level resolution and no better. Every figure here is a proportion with an interval in the per-seed files; the range across seeds is in brackets above.
