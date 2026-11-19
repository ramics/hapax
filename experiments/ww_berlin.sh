#!/usr/bin/env bash
# Second real-data probe: a divergent astrovirus in Berlin wastewater, with the
# authors' own assembly as the comparator.
#
# Why this dataset. Wyler et al. 2024 (Environment International 190:108875) deep-sequenced
# total nucleic acid from a Berlin treatment plant for 17 months and published, alongside the
# reads, their metaSPAdes and Trinity assemblies (GEO GSE228220, SRA PRJNA948850). That is the
# thing the cryptic-lineage runs could not offer: a real sample where someone else's assembly
# of the same reads is on the table, and where within-sample diversity actually exists.
#
# The target. In sample RNA_WW_211013_20-22hmix their SPAdes run produced one 6,823 bp
# astrovirus genome, 57% identical to HAstV-1 (NC_001943) -- divergent enough to be a novelty
# case, complete enough to score against.
#
# Honesty note on the comparison. MEGAHIT and SPAdes assembled the whole sample. hapax is given
# reads recruited by mapping to astrovirus bait, because it has no front end for separating a
# virus from 8.7M read pairs of sewage. The recruitment is NOT the reason hapax finds less:
# widening the bait from the SPAdes contig to the SPAdes contig plus MEGAHIT's own astro contigs
# plus four references moved the recruited set from 5,171 to 5,186 reads. The astrovirus simply
# has ~5.2k reads in this run. What the asymmetry does mean is that hapax's background contigs
# are not counted here, and finding the virus is not being tested.
#
# No seed: nothing here is stochastic.
set -euo pipefail

RUN=${RUN:-SRR23973193}            # RNA_WW_211013_20-22hmix
SAMPLE=${SAMPLE:-WW_211013_20-22hmix_S5_L001}
OUT=data/ww_berlin
THREADS=${THREADS:-4}
HAPAX=${HAPAX:-target/release/hapax}

mkdir -p "$OUT"; cd "$OUT"

# 1. reads, and the authors' assemblies
for n in 1 2; do
  until gzip -t "${RUN}_$n.fastq.gz" 2>/dev/null; do
    curl -sfL -C - -O "https://$(curl -s -G https://www.ebi.ac.uk/ena/portal/api/filereport \
      --data-urlencode "accession=$RUN" --data-urlencode 'result=read_run' \
      --data-urlencode 'fields=fastq_ftp' --data-urlencode 'format=tsv' |
      tail -1 | cut -f2 | tr ';' '\n' | sed -n "${n}p")"
  done
done
S=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228220/suppl
[ -s GSE228220_Spades_RNA_all.fa.gz ] || curl -sfL -O "$S/GSE228220_Spades_RNA_all.fa.gz"
[ -s GSE228220_RNA_kaiju_summary_annotated.tsv.gz ] || curl -sfL -O "$S/GSE228220_RNA_kaiju_summary_annotated.tsv.gz"
# NB: that SPAdes file is ~1 GB and its GEO copy ends mid-stream -- zcat recovers 101 of the
# 116 samples and errors at the end. Samples after WW_220616 are unreadable from it, which is
# why this probe uses 211013 rather than the astrovirus-richest sample (220629). The S-numbers
# in the kaiju table and in the SPAdes headers are different run indices: match by sample name.

# 2. the authors' contigs for this sample, and the astrovirus ones among them
zcat GSE228220_Spades_RNA_all.fa.gz 2>/dev/null |
  awk -v s="^>${SAMPLE}_NODE_" '/^>/{keep=($0 ~ s)} keep' > "spades_${SAMPLE}.fa"
for acc in NC_001943.1 NC_013060.1 NC_013443.1 JQ403108.1; do
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text"
done > astro_refs.fa
minimap2 -x asm20 -t "$THREADS" astro_refs.fa "spades_${SAMPLE}.fa" 2>/dev/null |
  cut -f1 | sort -u > astro_contig_ids.txt
seqkit grep -f astro_contig_ids.txt "spades_${SAMPLE}.fa" > bait_astro.fa
# the single near-complete genome is the scoring target
seqkit sort -l -r bait_astro.fa | seqkit head -n 1 > node37.fa

# 3. the incumbent, on everything
[ -d "megahit_$RUN" ] || megahit -1 "${RUN}_1.fastq.gz" -2 "${RUN}_2.fastq.gz" \
  -o "megahit_$RUN" -t "$THREADS" --min-contig-len 300 -m 0.7

# 4. hapax, on the recruited astrovirus reads
cat bait_astro.fa astro_refs.fa > bait_all.fa
minimap2 -ax sr -t "$THREADS" bait_all.fa "${RUN}_1.fastq.gz" "${RUN}_2.fastq.gz" 2>/dev/null |
  samtools sort -@2 -o astro.bam -
samtools index astro.bam
samtools view -b -F 0x904 astro.bam | samtools fasta - > astro_reads.fa
"../../$HAPAX" assemble -i astro_reads.fa -o hapax_astro -t "$THREADS" --min-contig 300

# 5. is there anything to phase? Counted from the reads, independently of any assembler.
echo "=== within-sample diversity at the target locus"
samtools mpileup -f node37.fa -d 1000 -Q 20 astro.bam 2>/dev/null | awk '
{ bases=toupper($5); gsub(/\^./,"",bases); gsub(/\$/,"",bases);
  ref=gsub(/[.,]/,"",bases); mm=0;
  for(i=1;i<=length(bases);i++){c=substr(bases,i,1); if(c~/[ACGT]/) mm++}
  tot=ref+mm; if(tot>=20){sites++; f=mm/tot; if(f>=0.2 && f<=0.8) het++} }
END{printf "  sites at >=20x: %d, minor allele 20-80%%: %d (%.2f%%)\n", sites, het, 100*het/sites}'

# 6. score each assembly against the authors' contig
echo "=== assemblies vs the authors' 6,823 bp astrovirus contig"
for f in "hapax_astro/contigs.fa" "megahit_$RUN/final.contigs.fa" "node37.fa"; do
  minimap2 -c -x asm10 -t "$THREADS" node37.fa "$f" 2>/dev/null > "$f.paf"
  printf "  %-40s contigs=%s covered=%s/6823 identity=%s misjoins=%s\n" "$f" \
    "$(cut -f1 "$f.paf" | sort -u | wc -l)" \
    "$(awk '{print $8"\t"$9}' "$f.paf" | sort -n | awk 'BEGIN{e=0;t=0}{if($1>e){t+=$2-$1;e=$2} else if($2>e){t+=$2-e;e=$2}} END{print t}')" \
    "$(awk '{m+=$10; b+=$11} END{if(b)printf "%.1f%%", 100*m/b}' "$f.paf")" \
    "$(awk '{print $1"\t"$8"\t"$9}' "$f.paf" | sort -k1,1 -k2,2n |
       awk '{if($1==p && $2-pe>1000) print $1; p=$1; pe=$3}' | sort -u | wc -l)"
done

# 7. what kind of variation is it, and where does hapax think it is coding?
echo "=== hapax variant effects at sites with depth>=20 and 20-80% minor allele"
awk -F'\t' 'NR>1 && $3>=20 && $6>=0.2 && $6<=0.8 {print $13}' hapax_astro/variants.tsv |
  sed 's/[0-9]//g' | sort | uniq -c | sort -rn
echo "=== region calls (see REPORT.md: 'dual' is over-called here)"
awk -F'\t' 'NR>1{l=$3-$2+1; t[$4]+=l; tot+=l} END{for(k in t) printf "  %s: %d (%.0f%%)\n", k, t[k], 100*t[k]/tot}' \
  hapax_astro/regions.tsv
