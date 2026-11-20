# Data manifest

`data/` is gitignored except for this file. The development sandbox cannot reach
NCBI at times; experiments/orf_annotation.py fetches what it needs. Record accession, source URL and
download date when a file is fetched.

| Purpose | Accession / source | Downloaded | Notes |
|---|---|---|---|
| HIV-1 reference (coordinates for *env*, *gag/pol*) | HXB2, GenBank K03455 | 2026-11-19 (efetch fasta + gb) | thesis EvolveAGene seed was HXB2 *env* |
| HIV-1 *env* population data with ground truth | Imo's curated alignments (thesis ch. 2/4) | not yet | ask Imo for location |
| HBV reference (dual-coding stress test) | NC_003977 | 2026-11-19 (efetch fasta + gb) | about half the genome is dual-coding (1,591 of 3,215 nt by CDS) |
| Public wastewater metagenome | to choose (SRA) | not yet | pick after HIV results |
| Human astrovirus 1 (ORF annotation truth) | NC_001943 | 2026-11-19 (efetch fasta + gb) | used by experiments/orf_annotation.py |

## Multi-swarm metagenome (experiments/fetch_swarm_panels.py)

| Purpose | Accession / query | Downloaded | Notes |
|---|---|---|---|
| swarm panel hiv_env_ZM247F | `Human immunodeficiency virus 1[orgn] AND ZM247F AND env AND 2400:3300[slen]` (EU166779.1..JX239326.1) | 2026-11-19 | 99 records, 72 unique ACGT-only haplotypes |
| swarm panel hiv_env_ZM246F | `Human immunodeficiency virus 1[orgn] AND ZM246F AND env AND 2400:3300[slen]` (EU166681.1..EU166758.1) | 2026-11-19 | 78 records, 61 unique ACGT-only haplotypes |
| swarm panel hiv_gagpol_2139 | `KC156179:KC156209[accn]` (KC156179.1..KC156209.1) | 2026-11-19 | 31 records, 19 unique ACGT-only haplotypes |
| swarm panel hcv_Pt3P / hcv_Pt3S | `MT416885:MT416975[accn]` split on isolate Pt3.Primary / Pt3.Super | 2026-11-19 | 52 / 39 records, 49 / 39 unique haplotypes |
| swarm panel hcv_Pt7 | `MT417081:MT417136[accn]` (MT417081.1..MT417136.1) | 2026-11-19 | 56 records, 54 unique ACGT-only haplotypes |
| swarm panel hbv_p03 | `KY881820:KY881843[accn]` (KY881820.1..KY881843.1) | 2026-11-19 | 24 records, 24 unique ACGT-only haplotypes |
| swarm panel hbv_p10 | `KY881980:KY882003[accn]` (KY881980.1..KY882003.1) | 2026-11-19 | 24 records, 24 unique ACGT-only haplotypes |
| background NC_000913.3 | Escherichia coli str. K-12 substr. MG1655, complete genome | 2026-11-19 | 4641652 nt |
| background NC_003228.3 | Bacteroides fragilis NCTC 9343, complete sequence | 2026-11-19 | 5205140 nt |
| background NC_004461.1 | Staphylococcus epidermidis ATCC 12228, complete sequence | 2026-11-19 | 2499279 nt |
| background NC_001416.1 | Enterobacteria phage lambda, complete genome | 2026-11-19 | 48502 nt |
| background NC_000866.4 | Enterobacteria phage T4, complete genome | 2026-11-19 | 168903 nt |
| background NC_024711.1 | Carjivirus communis, complete genome | 2026-11-19 | 97065 nt |
| background NC_001422.1 | Escherichia phage phiX174, complete genome | 2026-11-19 | 5386 nt |
| background K03455.1 | Human immunodeficiency virus type 1 (HXB2), complete genome; | 2026-11-19 | 9719 nt |


## Labeller audit (experiments/labeller_audit.py)

`data/labeller_audit/rabies_pairs.fa`: the 11 rabies genomes of the archived two-lineage
ladder (PAIRS in experiments/codon_assembler/experiments/two_lineage_rung.py), fetched from
NCBI efetch (nuccore, fasta) on 2026-11-20: PZ227608.1, PZ227621.1, PV491557.1, PZ324370.1,
MG458312.1, PZ227553.1, PZ227534.1, PZ227561.1, AB009663.2, AB044824.1, NC_001542.1.

## Mate-width thread (experiments/insert_profile.py, mate_separation.py,
## peptide_absence.py, novelty_blocks.py, mate_rho.py)

| Purpose | Accession / source | Downloaded | Notes |
|---|---|---|---|
| Berlin wastewater reads | SRR23973193 (PRJNA948850) via ENA fastq | 2026-11-20 | 8,731,224 pairs; 2x107 before trimming mate 2's 14 nt non-template prefix |
| the authors' own assembly | GSE228220_Spades_RNA_all.fa.gz (GEO) | 2026-11-20 | sample WW_211013_20-22hmix_S5_L001: 82,947 contigs, 29.8 Mb; subsets >=1 kb (3,530) and >=2 kb (1,416) used as insert-size references |
| protein index | RefSeq release viral.1.protein.faa.gz + .gpff.gz | 2026-11-20 | 722,823 proteins / 175.8 Maa; the gpff supplies taxonomy for family holdout |
| protein index (background range) | UniProt Swiss-Prot uniprot_sprot.fasta.gz (current release) | 2026-11-20 | added because 88% of Berlin reads score 1.0 against a viral-only index |
| annotated genomes | NC_001542.1, NC_001943.1, NC_045512.2, K03455.1, NC_001422.1, NC_000913.3 (efetch gbwithparts) | 2026-11-20 | CDS annotation for the coding and novelty tracks |

## Pair-ladder background (experiments/pair_ladder.py --background)

`data/berlin/SRR23973193_{1,2}.fastq.gz`: Berlin wastewater RNA virome (PRJNA948850), from ENA
(ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/093/SRR23973193/) on 2026-11-20; 404,863,127 and
419,810,100 bytes. `bg_1.fq`/`bg_2.fq` are the first 1,000,000 pairs (`seqkit head`), mate 2
trimmed with `seqkit subseq -r 15:-1`.

## Peptide absence noise floor (experiments/peptide_absence_noise.py)

Cached under `data/peptide_noise/`, fetched from NCBI efetch (nuccore, gbwithparts) on
2026-11-20. Six small genomes with simple CDS annotation, chosen for spread in GC and
codon usage rather than depth, since the quantity measured is a property of the genetic
code and the error model rather than of any one organism. 95 usable CDS in total, of
which lambda contributes 72 — the codon-usage breadth is accordingly narrow and that
caveat travels with the synonymous-fraction number.

| Accession | Organism | Usable CDS |
|---|---|---|
| K03455.1 | HIV-1 HXB2 | 6 |
| NC_001542.1 | Rabies lyssavirus | 5 |
| NC_001943.1 | Human astrovirus 1 | 1 |
| NC_003977.2 | Hepatitis B virus | 4 |
| NC_001416.1 | Enterobacteria phage lambda | 72 |
| NC_001422.1 | Escherichia phage phiX174 | 7 |

"Usable" means a simple (non-joined, non-fuzzy) CDS location, in frame, unambiguous, with
no internal stop — the conditions under which the annotation can be trusted to define a
true reading frame.

## Berlin quality profile (experiments/berlin_quality_profile.py)

`data/berlin/sample_R{1,2}.fastq`: the first 400,000 read pairs of SRR23973193, streamed
from ENA over HTTP (`ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/093/SRR23973193/`) on 2026-11-20.
No SRA toolkit is installed in this sandbox and none is needed for a prefix. Run metadata
from the ENA portal API the same day: 8,731,224 pairs, 1,852,591,930 bases, Illumina
NovaSeq 6000, modal read length 107 nt.

## Mate-width thread, later arms (absence_curve.py, rrna_contigs.py, rrna_match_source.py,
## subsample_uniref.py, absence_residual.py)

| Purpose | Accession / source | Downloaded | Notes |
|---|---|---|---|
| public protein index, large | UniRef50 (uniref50.fasta.gz, current release) | 2026-11-20 | **streamed, never stored**: nested every-10th and every-100th subsamples kept (4,124,353 and 412,436 records; 1.27 and 0.127 Gaa). Completeness checked against UniProt REST `x-total-results` for `query=identity:0.5` = 41,243,524 clusters, which the stream reaches exactly |
| ribosomal reference | SILVA 138.1 SSURef_NR99 + LSURef_NR99 (arb-silva.de Exports) | 2026-11-20 | 197 MB + 67 MB gzipped, 605,794 sequences; used to label assembly contigs as ribosomal by exact 31-mer containment |
| sample-self protein index | Prodigal (pyrodigal, meta) over the GSE228220 contigs for this sample | 2026-11-20 | 79,667 proteins / 8.04 Maa, split into `berlin_sample_proteins_{nonrrna,rrnacontigs}.faa` by the SILVA contig label |

