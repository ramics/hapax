# Prior art: subtracting the known, and what is left

Survey done 2026-11-20 for the thread "subtract the known, then characterise the residue".
The point of the survey is to say which parts of Imo's thesis - *find the most efficient
way to extract the known, then score the unknown* - are already solved engineering, so the
contribution is not an index structure that three groups have already built.

Read this before proposing to build anything here.

## 1. Subtraction as a pipeline step is twenty years old and settled

Host and background subtraction followed by homology search on the remainder is the
standard viral-discovery workflow, not a new idea.

- **VIP** (Li et al., *Sci Rep* 2016) filters background reads, then classifies the
  remainder by nucleotide and then remote amino-acid homology, then assembles. The
  ordering - nucleotide first, protein second, because protein reaches further - is the
  same ordering this thread's cascade uses, and it was published a decade ago.
- **IMSA** (Dimitrov et al. 2013) and the host-subtraction validation study of Daly et al.
  (*PLoS ONE* 2015, e0129059) are the reference points for "subtract host, assemble the
  rest".
- A systematic review of novel-virus metagenomics (Kinsella et al. 2022, PMC9784588)
  found **little agreement between groups on the optimal pipeline**, and that more than
  half of successful viral-discovery papers used *de novo* assembly first with no
  subtraction step at all. So there is no canonical ordering to beat, which is part of why
  a measured curve is worth having.
- **SURPI+** (Chiu lab, *Nat Commun* 2024) tests novel-virus detection against a database
  with all known human viral sequences removed - the family-holdout design this project
  used in PR #8 and PR #10.

Known cost of subtraction, which the write-up must quote rather than discover: the SURPI
evaluation measured that pre-processing removed **0.43% of true viral reads**. Subtraction
is not free, and "removed as known" is never exactly "known".

## 2. Cheap index structures are solved, by several groups, and we should use them

Every one of these is a better engineered answer to "store a lot of sequence in few bytes"
than anything this project would write.

| Structure | What it is | Reported cost |
|---|---|---|
| **Kraken2** compact hash (Wood et al. 2019) | minimizer -> LCA, with `--max-db-size` subsampling to a byte budget | the prebuilt capped builds used here: 0.6 GB (viral), 8 GB, 16 GB; the uncapped 2020 Standard was 47 GB |
| **COBS** (Bingmann et al. 2019) | compact bit-sliced signature index, successor to **BIGSI** | ~3.0 GiB on the benchmark where HowDe-SBT was ~1.9 GiB; the only one of its cohort that scaled to the full 100,000-microbe set |
| **MetaProFi** (Srikakulam et al., *Bioinformatics* 2023, btad101) | chunked Bloom filter over **amino-acid** k-mers, queryable with nucleotide input via internal six-frame translation; exact or threshold-T approximate search | UniProtKB + Tara Oceans + 2,585 RNA-seq datasets indexed in <1,290 min using <60 GiB RAM |
| **sourmash / branchwater** | FracMinHash sketches, not an index of every k-mer | a branchwater database of ~7.9 GB |
| **Logan** (Chikhi et al. 2024) | assembly of the entire SRA: 4.59 Pbp of unitigs, 0.90 Pbp of contigs; `logan-search.org` is a k=31 index over 23.4M accessions | one DIAMOND pass over all Logan contigs: ~11 h on 60k vCPUs, ~$10k |
| **Serratus** (Edgar et al., *Nature* 2022) | 5.7M SRA runs, 10.2 petabases, searched for RdRp; 881,167 palmprints, 260,808 unique | the palmprint is explicitly offered as the cheap index |

Two things follow. First, **MetaProFi is the right tool for a peptide-k-mer absence
index** and the sibling thread should not write one. Second, the honest framing of our
disk limit is not "we cannot index everything" - Logan and Serratus can, for money - but
"what does each increment buy", which nobody in this table reports as a curve against a
real background.

## 3. The closest prior art: WVDB, the wastewater virome database

Kantor, Shakya et al., *medRxiv* 10.64898/2026.05.19.26353600 v2 (18 June 2026), "A
genome-resolved view of the wastewater RNA virome". This is the paper that occupies the
most ground and it appeared five months ago.

- 321 US wastewater metagenomes, 11 sites, 6 cities, median 1.6 billion reads per sample.
  MEGAHIT per sample, geNomad 1.11.0 + checkV 1.0.3, vclust at 95% ANI, two independent
  complete assemblies required per vOTU: **21,015 vOTUs**, 79% ssRNA bacteriophage.
- **The number we have to beat or match:** Kraken2 2.17.1 against NCBI core-nt
  (2025-10-15) left **on average 56% of reads unclassified (range 23.1-83.6%)**, and of
  those, **10.9% (range 0.0-40.3%)** mapped to WVDB by bowtie2. Separately, 18.3% of
  clean non-rRNA reads mapped to WVDB.
- They state plainly that **"read-based methods cannot currently separate novel from
  previously seen unclassified sequences"**, and their own novelty argument is at genome
  level: just 407 of 21,015 vOTUs had a species-level BLASTN hit to core-nt, 3,845 (18%)
  had a genus-level hit to any database, 12.6% had no VOG protein hit.
- Their collector's curve plateaus in *new genomes* at around 100 samples.

What they do **not** do, checked in the full text: no curve of unclassified-read removal
against database size; no per-read characterisation of the residual unclassified reads (no
GC, coding, repeat, adapter, chimaera or motif fields); no read-level procedure for
removing WVDB-mapped reads from the unclassified pool and no bowtie2 identity threshold
given. Residue handling is left as future work ("co-assembly of reads that don't map to
WVDB").

One internal discrepancy to be careful about when citing them: the abstract says WVDB
"captured around one fifth of the reads that could not be classified by Kraken2" while
Results 2.2 gives a mean of 10.9%. The 18.3% figure is of *all* clean non-rRNA reads, not
of the unclassified ones, which is the most likely source of the abstract's number. Quote
the body, not the abstract.

## 4. What unclassified fractions other people report

Useful as sanity bounds on our own numbers, with the caveat that each is a different
sample type, database and threshold.

- 56% mean unclassified against core-nt on US wastewater RNA (WVDB, above).
- ~72% unassigned (range 34-97%) against Kraken2 Standard 2024 on a wastewater **DNA**
  dataset.
- 41% no hits, ~55% bacteria, ~2% viruses in one wastewater set (METAGENE-1 preprint).
- Kraken2's classified fraction on a mock community fell from ~51% at confidence 0.05 to
  ~5% at the strictest threshold tested (*Sci Rep* 2025), and a SILVA-based Kraken2 build
  classified under 2% of reads even at the most lenient threshold. **This is why the
  confidence sweep in our cascade is not optional**: a single classified-fraction number
  is a threshold choice wearing a measurement's clothes.
- NAO's own wastewater notebooks: ribosomal fraction **0.67%** in one DNA run but
  **8-21%** in unenriched RNA runs, human 0.4-6.5%, viral 0.6-8.6%.

## 5. The chimaera and artefact constraint

Not a footnote in a residue defined by matching nothing.

- A dengue/Nextera study found most or all host-virus chimaeras were artificial (PCR
  template switching).
- Amplicon chimaera rates exceed **70% for the rarest taxa** (Haas et al. 2011), which is
  our target by construction, and they form reproducibly across replicates, so replicate
  agreement cannot filter them.
- `ViFi` (NAR 2018) already does the divergence-tolerant one-sided-anchor read, so that
  primitive is occupied (recorded in PR #8's prior art).

## 6. Where our claim actually lies, after the survey

Dropped as not ours:

- subtraction-then-residue as a workflow (VIP, IMSA, SURPI+, and every virome pipeline);
- cheap indexes (Kraken2, COBS/BIGSI, MetaProFi, sourmash, Logan, Serratus);
- a wastewater-specific catalogue of unclassified sequence (WVDB, five months ago);
- the per-million false-positive unit (PACIFIC 2021, corrected in PR #8).

What survives, and is what this thread should deliver:

1. **The removal-per-gigabyte curve on a real wastewater library**, with every tier's
   threshold swept rather than fixed, including the point that nobody else reports: the
   sample's *own* assembly as a subtraction index, priced in bytes against the public
   catalogues.
2. **A per-read evidence report on the residue with every field calibrated** against a
   composition-matched null and against reads the cascade did subtract - which is exactly
   the thing WVDB names as missing and leaves as future work.
3. **An honest extrapolation to a comprehensive index**, stated as an extrapolation,
   because the comprehensive index does not fit on the disk and saying so is a result.

Nothing in this list is a novelty classifier, and that is deliberate: this project has
measured four separate single-read features failing as novelty scores for the same reason
(PR #8), and WVDB independently states the read-level separation does not currently exist.
