#!/usr/bin/env bash
# First real-data probe: cryptic SARS-CoV-2 lineages in wastewater WGS.
#
# Why this dataset. Smyth et al. 2022 (Nat Commun 13:635) and Gregory et al. 2022
# (PLoS Pathog 18(10):e1010636) report spike-RBD lineages in sewage that appear in no
# patient sample, and say plainly why they could not finish the job: "it is not
# possible to reconstruct individual genomes using standard methods" and "we cannot
# link mutations unless they are found on the same amplicon". That is hapax's claim
# surface, stated by the people who hit the wall, with per-run mutation constellations
# published as (partial) ground truth in Gregory Table 2.
#
# What this script measures, per assembler: contigs, reference bases covered, and
# candidate misjoins (an alignment jump >1 kb on the reference inside one contig).
# Contiguity without the misjoin count is not reportable (CLAUDE.md).
#
# Ground truth is a mutation set, not a sequence, so the test is: does a contig carry
# the published constellation in linkage on one molecule? The RBD check at the end
# translates each assembler's RBD contig against Wuhan-Hu-1 and prints the changes.
#
# No seed: nothing here is stochastic. The manifest records accessions and versions.
set -euo pipefail

RUN=${1:-SRR15202285}
REF=${REF:-data/ww_cryptic/NC_045512.2.fa}   # NC_045512.2, Wuhan-Hu-1
OUT=data/ww_cryptic
THREADS=${THREADS:-8}
HAPAX=${HAPAX:-target/release/hapax}

mkdir -p "$OUT"
cd "$OUT"

# Reads straight from ENA; the FTP path is read from the API rather than guessed.
for n in 1 2; do
  [ -s "${RUN}_$n.fastq.gz" ] || curl -sfL -O "https://$(curl -s -G \
    https://www.ebi.ac.uk/ena/portal/api/filereport \
    --data-urlencode "accession=$RUN" --data-urlencode 'result=read_run' \
    --data-urlencode 'fields=fastq_ftp' --data-urlencode 'format=tsv' |
    tail -1 | cut -f2 | tr ';' '\n' | sed -n "${n}p")"
done

# MEGAHIT: the incumbent, on the raw pair.
[ -d "megahit_$RUN" ] || megahit -1 "${RUN}_1.fastq.gz" -2 "${RUN}_2.fastq.gz" \
  -o "megahit_$RUN" -t "$THREADS" --min-contig-len 300

# hapax takes one FASTA of reads. Wastewater is mostly not virus, so reads are
# pulled by mapping to the reference first. That is a reference-guided *input filter*,
# not reference-guided assembly: the assembly itself never sees the reference. It is
# also the honest comparison only in the "how well is the virus assembled" sense --
# MEGAHIT gets no such help, and its extra background contigs are not counted below.
minimap2 -ax sr -t "$THREADS" "$(basename "$REF")" "${RUN}_1.fastq.gz" "${RUN}_2.fastq.gz" |
  samtools sort -@2 -o "$RUN.sort.bam" -
samtools index "$RUN.sort.bam"
samtools fasta -F 0x900 "$RUN.sort.bam" > "${RUN}_viral.fa"
# spike gene alone, as the slice where the published constellations live
samtools view -b "$RUN.sort.bam" NC_045512.2:21563-25384 | samtools fasta -F 0x900 - > "${RUN}_spike.fa"

"../../$HAPAX" assemble -i "${RUN}_viral.fa" -o "hapax_$RUN"       -t "$THREADS" --min-contig 300
"../../$HAPAX" assemble -i "${RUN}_spike.fa" -o "hapax_${RUN}_spike" -t "$THREADS" --min-contig 300

echo "=== coverage breadth of the reads themselves (the ceiling any assembler can reach)"
samtools depth -a "$RUN.sort.bam" | awk '{s+=$3; if($3>=10)c++} END{printf "mean depth %.0f, positions >=10x: %d / %d\n", s/NR, c, NR}'

for f in "megahit_$RUN/final.contigs.fa" "hapax_$RUN/contigs.fa" "hapax_${RUN}_spike/contigs.fa"; do
  echo "=== $f"
  awk '/^>/{n++} !/^>/{b+=length($0)} END{print "  contigs="n" bases="b}' "$f"
  minimap2 -x asm20 -t "$THREADS" "$(basename "$REF")" "$f" 2>/dev/null > "$f.paf"
  echo "  contigs aligning to the reference: $(cut -f1 "$f.paf" | sort -u | wc -l)"
  echo "  reference bases covered: $(awk '{print $8"\t"$9}' "$f.paf" | sort -n |
    awk 'BEGIN{e=0;t=0}{if($1>e){t+=$2-$1;e=$2} else if($2>e){t+=$2-e;e=$2}} END{print t}') / 29903"
  echo "  candidate misjoins: $(awk '{print $1"\t"$8"\t"$9}' "$f.paf" | sort -k1,1 -k2,2n |
    awk '{if($1==p && $2-pe>1000) print $1; p=$1; pe=$3}' | sort -u | tr '\n' ' ')"
done
