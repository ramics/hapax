# How many measurements are in the per-read report?

30,000 residue reads over 5 seeds, 13 numeric fields.

**Participation ratio 10.5 of 13 fields.** That is the effective number of independent axes: the report is measuring about 10 things, not 13.

## Groups of fields that move together (|rho| >= 0.5)

- fold_mfe, gc
- entropy3  (independent)
- homopolymer  (independent)
- dup_kmer16  (independent)
- longest_orf_aa, stop_free_frames
- CpG_oe  (independent)
- UpA_oe  (independent)
- fold_z  (independent)
- library_abundance  (independent)
- cross_sample_samples  (independent)
- closest_bitscore  (independent)

## Does the field separate residue from the reads the cascade identified?

| field | residue mean | known mean | standardised difference |
|---|---|---|---|
| library_abundance | 8.89 | 89.07 | -1.12 |
| UpA_oe | 0.7783 | 0.8911 | -0.40 |
| stop_free_frames | 0.809 | 0.4602 | +0.39 |
| fold_mfe | -26.98 | -24.59 | -0.28 |
| longest_orf_aa | 58.77 | 56.45 | +0.24 |
| entropy3 | 0.8739 | 0.8838 | -0.23 |
| fold_z | -0.4934 | -0.3994 | -0.11 |
| gc | 0.5034 | 0.4964 | +0.09 |
| CpG_oe | 0.9567 | 0.9408 | +0.05 |
| dup_kmer16 | 0.002516 | 0.002046 | +0.01 |
| homopolymer | 4.494 | 4.509 | -0.01 |
| closest_bitscore | 31.6 | *not searched* | *cannot be computed* |

**`closest_bitscore`** read identically zero across the control, which means the measurement was never made there - the protein search ran on the residue only. Reported as *not searched* rather than as a difference: scored naively, `closest_bitscore` showed the largest separation in this table (+9.5 standardised) and it was an artefact of a population nobody searched.

Independence and usefulness are different properties: a field can be its own axis and still say nothing about whether a read is interesting, and both columns are needed to read the report honestly.

Two results from the gate experiment (PR #17) belong beside this table. First, the composition and quality fields were tested JOINTLY on residue against subtracted reads and reach AUC 0.686, where the full space reaches 0.905 and the coverage axes alone reach 0.904: the axes PR #8 killed one at a time also fail together, so a high participation ratio among them is dimensionality without discrimination. Second, the residue is a continuum - no cluster count beats k=1 against a covariance-matched reference over five seeds - so nothing here should be read as carving groups.

A measured caveat on `library_abundance`: it comes from a count-min sketch that only ever overstates and that includes the read's own occurrences, so it cannot resolve coverage SHAPE. PR #17's exact 31-mer counter (`experiments/kmer_coverage/`, whole library, own occurrences subtracted) does, at 0.76-0.78 AUC per shape field. Those fields belong in THIS matrix when they are adopted rather than being assumed independent of abundance.
