# What is in this water, asked the way a public-health reader asks it

`experiments/public_health_profile.py` and `experiments/known_strain_divergence.py`,
seed 1, SRR23973193, local, no spend. 377 s and 61 s. Committed outputs:
`profile_s1.json`, `strain_divergence_s1.json`; panel definition
`experiments/pathogen_panel.json`.

## Why this exists

Everything else in this repo asks what the indexes CANNOT name, because that is the
research question. A public-health reader asks the opposite question first, and the run
directory could not answer it: **the cascade records whether a read was removed and by
which tier, never what it was removed AS.** "What fraction of this library is bacterial"
was not a number this pipeline could produce.

## The whole library, classified, at three thresholds

Kraken2 2.17.2, PlusPF-08, paired, mate 2 trimmed 14, over all **8,731,224 read pairs**
in 20 s per pass.

| | confidence 0 | **0.1** | 0.2 |
|---|---|---|---|
| classified | 37.73% | **6.88%** | 0.50% |
| Bacteria | 31.87% | **5.21%** | 0.45% |
| Viruses | 4.01% | **0.54%** | 0.01% |
| Eukaryota | 0.66% | **0.02%** | 0.001% |
| Archaea | 0.03% | **0.002%** | - |

**A 75-fold range across one parameter**, which is the standing finding of this repo
restated on the whole library instead of a block. Of the part that is classified at all
at 0.1, 75.7% is bacterial, 7.9% viral.

The most abundant virus is **Tobamovirus** (pepper mild mottle and relatives, 34,542
pairs at 0.1, 0.40% of the library), which is diet-derived and the classic faecal
indicator - a sanity check that the instrument is reading what it should.

**Most bacterial reads never reach a species.** 321,748 pairs at 0.1 are placed no deeper
than the domain. The page labels those rows with the rank they reached rather than
printing "Bacteria" as if it were an organism, which the first version did.

## The screening panel, and the null that sets its threshold

44 taxa, chosen for what a public-health reader asks of sewage and **not** by looking at
this sample. Every entry is reported found or not found; a not-found row carries the
rule-of-three 95% upper bound, which is a bound on what this sequencing could have seen.

**17 found, 27 not found** at confidence 0.1. Strongest rows: Mamastrovirus (27 pairs,
1 in 323,000), Norovirus (16, 1 in 546,000), Leptospira (14), Cryptosporidium (12),
*Pseudomonas aeruginosa* (8). Not found includes SARS-CoV-2, influenza A and B, measles,
mumps, rubella, rotavirus A, hepatitis A and E, *M. tuberculosis* and mpox.

**The threshold is not a preference, it is the calibration.** The same panel run over
this library with every read shuffled lights up **24 of 44 rows at confidence 0** and
**0 of 44 at 0.1** (43,151 pairs classified on the null at 0, 27 at 0.1). So a panel read
at confidence 0 is more than half noise, and 0.1 is where the null goes quiet.

**What the null does not control, and nothing here does:** misassignment between close
relatives. A gut *Bacteroides* read landing on a *Shigella* k-mer is a real k-mer match
and a shuffle says nothing about it. Several found rows sit at one or two pairs and
should be read as "worth a targeted assay", never as a detection.

## How far the known viruses have drifted

On the sample's own SPAdes contigs, because at read level this measurement is
impossible: 107 nt is ~35 amino acids, where strain variation is mostly synonymous and
identity sits in the chance band. 377 contigs are viral or under a watchlist taxon, 368
have a protein hit at E<1e-5, **187 taxa have an alignment of at least 60 aa and 89 of
them sit below 90% identity** to their closest reference protein.

Examples: Steitzviridae 47.6% over 499 aa; *Hohrdovirus edaphohabitans* 70.2% over 721
aa; Cripavirus drosophilae 81.2% over 1,850 aa; **Human astrovirus 92.2% over 762 aa**
with a 6,823 nt contig; Mamastrovirus (on the panel) 98.2% over 318 aa.

**A finding worth carrying: Kraken2 confidence is length-dependent, and not in the
direction a reader assumes.** At confidence 0.1 the viral calls on these contigs collapse
from **1,021 to 11**, and the survivors are the SHORT contigs - a 9 kb divergent viral
contig has thousands of minimizers and matches few of them. The cascade's read-level
operating point does not transfer to contigs. This script therefore assigns at
confidence 0 and requires a protein alignment before a row exists: Kraken2 proposes, the
alignment disposes.

## Two numbers that look like a contradiction

The classifier names 6.88% of the library; the cascade calls 77.03% of it identifiable.
Different questions. **Naming** a read means placing it on a taxon. **Recognising** one
means matching it to sequence we have seen: 37.99% of this library matches ribosomal RNA
and a further 30.30% marginal matches contigs assembled from this very sample, and
neither carries a name. Most of what we can remove as known, we remove by matching
sequence nobody has named. The two instruments agree where they overlap: the same Kraken2
index removes 6.81% inside the cascade on its 100k block against 6.88% here on the whole
library.

## Blind spots

- A metatranscriptome. Every fraction is of the RNA that was sequenced, not of the
  organisms in the water.
- PlusPF-08 and not PlusPF-16, because the 15 GB hash does not fit this container
  alongside the taxonomy. The cascade measures PlusPF-16 at 21.67% alone against
  PlusPF-08's 6.81%, so the classified fractions here are low by roughly that factor.
- Protein identity understates nucleotide divergence; a nucleotide comparison needs
  reference genomes this sandbox cannot fetch.
- Assembly fuses lineages below ~3.8% divergence, so the strain table measures drift from
  the DATABASE and is blind to strain mixtures WITHIN the sample.
- One sample on one day. `results/dna_rna_pair/` measures this plant's rRNA fraction
  moving 38.0% to 67.1% in eight months.
