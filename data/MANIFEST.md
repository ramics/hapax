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

## Subtraction cascade (experiments/subtract_known.py, residue_report.py, subsample_uniref.py)

All fetched 2026-11-20 into the gitignored `data/`. Sizes are what the index-size axis of
`results/subtract_known/REPORT.md` prices, so they are recorded as bytes on disk, not as
record counts alone.

| Purpose | Source | Downloaded | Notes |
|---|---|---|---|
| Berlin wastewater reads | SRR23973193 (PRJNA948850) via ENA fastq | 2026-11-20 | 8,731,224 pairs; ENA records library_strategy=RNA-Seq, library_selection=cDNA, source TRANSCRIPTOMIC, NovaSeq 6000. No rRNA depletion is recorded, which the rRNA tier's 38% bears out |
| rRNA reference | SILVA 138.2 SSURef_NR99 + LSURef_NR99 (ftp.arb-silva.de) | 2026-11-20 | SSU 510,495 records / 747.4 Mnt; U->T converted, headers cut to the accession |
| Kraken2 viral | genome-idx k2_viral_20261022 | 2026-11-20 | 0.65 GB hash + taxo |
| Kraken2 Standard-8 | genome-idx k2_standard_08_GB_20261022 | 2026-11-20 | 8.0 GB |
| Kraken2 PlusPF-8 | genome-idx k2_pluspf_08_GB_20261022 | 2026-11-20 | 8.0 GB |
| Kraken2 PlusPF-16 | genome-idx k2_pluspf_16_GB_20261022 | 2026-11-20 | 16.0 GB |
| this sample's own assembly | GSE228220_Spades_RNA_all.fa.gz, records prefixed `WW_211013_20-22hmix_S5_` | 2026-11-20 | 82,947 contigs / 35.2 MB, matching the count recorded for this sample in PR #10 |
| all 116 samples' assembly, >=1 kb | same file, every sample | 2026-11-20 | 270,626 contigs / 577 MB. **The GEO file is truncated**: it decompresses to 14,454,331 records and then ends mid-stream (`unexpected end of file`) at sample `WW_220616_20-22hmix_S65`, although the served Content-Length matches the bytes received. This sample's own contigs are complete; the all-sample tier is therefore a lower bound on what a complete file would subtract |
| RefSeq viral proteins | refseq/release/viral/viral.1.protein.faa.gz | 2026-11-20 | DIAMOND db 0.24 GB |
| Swiss-Prot | uniprot_sprot.fasta.gz, current release | 2026-11-20 | DIAMOND db 0.29 GB |
| UniRef50 nested subsets | streamed from ftp.uniprot.org without storing the file (experiments/subsample_uniref.py, seed 1) | 2026-11-20 | 41,243,524 records streamed; kept 2% (823,978 records / 254.6 Maa, DIAMOND db 0.36 GB) and 10% (4,121,360 / 1,270.9 Maa, db 1.82 GB). A 25% subset (10,307,898 / 3,180.3 Maa) was written and then deleted for disk; a DIAMOND db of the full 12.7 Gaa would be ~30 GB against ~15 GB free |
| Pfam-A | ftp.ebi.ac.uk Pfam current_release Pfam-A.hmm.gz | 2026-11-20 | 3.65 GB uncompressed, used for the motif field |

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

## Ribosomal content of the library (experiments/rrna_content.py, rrna_control.py,
## prodigal_on_rrna.py, rrna_peptide_reach.py)

The SILVA download is the same one the section above records; it serves both threads, which
label different things with it — contigs there, reads here.

`data/silva/SSURef_NR99.fasta.gz` and `data/silva/LSURef_NR99.fasta.gz`: SILVA release
138.1 small- and large-subunit ribosomal RNA references, NR99 (clustered at 99% identity),
downloaded from `ftp.arb-silva.de/release_138_1/Exports/` on 2026-11-20. 197,133,995 and
67,498,420 bytes; 510,508 and 95,286 sequences; 747,391,099 and 279,545,413 bases. Release
138.2 is not present at that path and returns 404. Used by
`experiments/rrna_content.py` to identify ribosomal reads positively rather than by the
absence of a gene call. The derived index (`data/silva/index_k31.npy`, 15,821,248 canonical
31-mer codes, 1-in-8 hash subsample) is a cache and is rebuilt by the script if deleted;
the retained count is far below an eighth of the 1.027 Gbp of positions because ribosomal
RNA is conserved enough to be highly redundant at k = 31.

`data/berlin/prefix_R{1,2}.fastq`: the first 1,620,000 read pairs of SRR23973193 from the
same ENA path, 2026-11-20, extending the 400,000-pair prefix above so the five read blocks
used by `experiments/absence_residual.py` on branch `claude/project-thread-clcuow` are all
available. The R1 stream closed early twice and was retried until the full line count
landed; a truncated prefix is detectable because every block is addressed by line offset.

## Residue structure (experiments/residue_structure.py, experiments/kmer_coverage/)

| Purpose | Source | Downloaded | Notes |
|---|---|---|---|
| Berlin wastewater reads, complete run | SRR23973193 via ENA fastq (`ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/093/`) | 2026-11-20 | both mates complete, 8,731,224 pairs each, line counts verified against the manifest of PR #16's cascade. The whole run is needed because this experiment counts every query k-mer across the entire library rather than a prefix |
| planted-control genome | Human rhinovirus A1B, D00239, ENA fasta | 2026-11-20 | 7,152 nt. The clean holdout of the spike ladder: absent from this sample's own assembly, unlike the other two genomes used there |

Intermediates written to the gitignored `data/residue_structure/`: `coverage.tsv` and
`composition.tsv` (568,500 pairs, 24 MB each), `spike.tsv`, `null.tsv`, and the residue
read names taken from `results/subtract_recut/residue_rows_s{1..5}.jsonl`. The points the
figure plots are committed to `results/residue_structure/joint_points_s1.json` so the figure
has provenance without the intermediates.

## The matched DNA/RNA pair (experiments/dna_rna_pair.py)

Fetched 2026-11-21 into the gitignored `data/berlin_pair/`, from ENA, each file verified
against the MD5 ENA publishes for it because ENA's FTP front end has served this project
truncated gzip before.

| Purpose | Accession | Downloaded | Notes |
|---|---|---|---|
| matched DNA library | SRR23973375 (`DNA_220616_20-22hmix_cer_BIOMcombi`, GSM7117187) | 2026-11-21 | PRJNA948850, `library_source=GENOMIC`, `library_strategy=OTHER`; 11,311,548 pairs / 2.37 Gbp; NovaSeq 6000. Measured to be a VIROME rather than a metagenome (0.01% rRNA, top taxon *Darmviridae*) - see `results/dna_rna_pair/REPORT.md` |
| matched RNA library | SRR23973265 (`RNA_WW_220616_20-22hmix`, GSM7117161) | 2026-11-21 | Same plant, same date, same 20-22h composite, same extraction and prep as the DNA member; 11,627,981 pairs / 2.42 Gbp. Same `RNA_WW_` naming series as our calibration sample, eight months later |

The project holds 239 runs in all, 215 `TRANSCRIPTOMIC` and 24 `GENOMIC`, and 22 of the 24
DNA libraries have an RNA library of the same date and collection window. The DNA set spans
2022-05-25 to 2022-07-19, so there is no DNA library for our 2021-10-13 calibration sample.
