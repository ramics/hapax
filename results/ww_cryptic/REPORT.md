# Real wastewater, first probe: cryptic SARS-CoV-2 lineages (SRR15202285)

Run with `experiments/ww_cryptic.sh SRR15202285`. Manifest beside this file.

## The dataset and why it was chosen

`SRR15202285`, wastewater WGS, New York, 28 June 2021 (BioProject `PRJNA746354`; one of
the public runs Gregory et al. 2022 re-screened, their Table 2). 332k read pairs, 2x100 bp,
2x10 MB. Picked because the authors of that paper and of Smyth et al. 2022 state the
failure in hapax's own terms: "it is not possible to reconstruct individual genomes using
standard methods"; "we cannot link mutations unless they are found on the same amplicon";
"mutations appearing on different reads cannot be linked together".

## The result: this dataset does not discriminate

| | contigs | aligning to SARS-CoV-2 | ref bases covered | candidate misjoins |
|---|---|---|---|---|
| MEGAHIT (raw pair) | 161 | 25 | 12,097 / 29,903 | 1 (k119_85) |
| hapax (mapped reads) | 147 | 18 | 10,612 / 29,903 | 1 (contig_1) |
| hapax (spike reads only) | 3 | 3 | 2,057 | 0 |

A misjoin here is an alignment jump >1 kb on the reference inside one contig.

**The ceiling is the reads, not the assembler.** 61.6% of reads map to SARS-CoV-2 at mean
depth 1349x, but only 11,994 of 29,903 positions reach 10x. Coverage is island-shaped
(amplicon-like dropout). MEGAHIT's 12,097 covered bases *is* the 11,994 the coverage
allows. There is no contiguity left on the table for a better overlap criterion to win.

**Both assemblers recover the published cryptic constellation, identically.** In the RBD
island (ref 22805-23205) hapax's `contig_3` and MEGAHIT's `k119_87` give byte-identical
difference strings against Wuhan-Hu-1: 19 amino-acid changes plus a 3 bp deletion at
codon 482, in linkage on one contig —

```
K444S V445K G446S Y449R L452Q N460K K462R S477N T478E E484P
F486I F490P G496S Q498Y P499S N501T Y505H V511I  (+ L425 synonymous, del nt23008-23010)
```

Gregory Table 2 describes this run as 100% cryptic with ~20 RBD changes including Y449R,
Q498Y, N501T, Y505H. All four are present. (Table 2 was read via a model summary of the
PMC page, not verbatim; the row should be checked against the PDF before it is cited.)

So the linkage the authors could not establish *is* established here — by plain MEGAHIT
as well as by hapax. The reason is visible in the same table: the cryptic lineage is 100%
of the viral reads in this run. There is no mixture to deconvolve, so the RBD assembles
as a clean single haplotype. The authors' difficulty was spanning the gaps *between*
islands and attributing mutations across them, which is a library-prep problem.

**Both make the same misjoin.** MEGAHIT `k119_85` (2,400 bp) and hapax `contig_1`
(2,192 bp) each weld ref ~400-870 to ref ~15.3-16.3 kb. Two independent assemblers
converging on the same false join points at chimaeric molecules in the library rather
than at either graph; Gregory's reviewers raised PCR chimaerism from ~50 cycles. This is
a caution for the whole dataset family: apparent long-range linkage in these runs needs
evidence that it is not an amplification artefact.

## What this means for the dataset shortlist

Keep the cryptic-lineage runs as a *correctness* fixture: a real-data case with a
published constellation that an assembler must reproduce without inventing linkage.
They are not a contiguity benchmark, and they will not show a codon-space advantage,
because the hard part is coverage breadth and the viral population in any one run is
near-clonal.

A discriminating real dataset needs within-sample diversity at the same locus and
coverage that is broad as well as deep. The runs with mixtures (Gregory's `SRR15202284`,
98% cryptic; `SRR15128983`, 99% cryptic, PANGO A.29) are the next ones to try, and the
Berlin deep total-nucleic-acid dataset (Wyler et al. 2024) is the better-shaped candidate
because it tracks entero/noro/astrovirus "down to strain or individual nucleotide
variation levels" in untargeted libraries.

## Caveats

- hapax was given reads pre-filtered by mapping to the reference; MEGAHIT was given the
  raw pair. The comparison is therefore about assembly of the virus, not about finding it.
  hapax's background contigs are not in this table because it never saw the background.
- One run, one seed-free deterministic pass. Nothing here is averaged over anything.
- hapax ran 663k reads in 64 s on 8 threads; 12.7% of candidate pairs fell through to the
  codon HMM.
