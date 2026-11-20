#!/usr/bin/env bash
# Does --split-lineage pull the Berlin astrovirus lineages apart (results/ww_berlin_lineage)?
#
# Run after experiments/ww_berlin_consensus.sh, which leaves astro_reads_trim.fa (mate-2
# non-template prefix removed) in data/ww_berlin. The comparators are the authors' own SPAdes
# contigs (bait_astro.fa): NODE_37 is the 6,823 bp target, the others are lineages 10-14% from it.
#
# Two scores per assembly:
#  - target: every contig against NODE_37 alone (minimap2 asm10), as in ww_berlin_consensus.
#  - per lineage: each contig's primary alignment against all SPAdes astrovirus contigs (asm20,
#    no secondaries), so a contig is scored against the lineage it is closest to. A contig that
#    blends two lineages scores low against both; one whose pieces map to different contigs
#    end to end is listed as a multi-lineage contig (a candidate chimaera, or a genome SPAdes
#    itself broke: the comparator is an assembly, not truth).
# MEGAHIT is the whole-sample assembly on the same trimmed reads, from ww_berlin_consensus.
#
# No seed: nothing here is stochastic.
set -euo pipefail

OUT=data/ww_berlin
THREADS=${THREADS:-4}
HAPAX=${HAPAX:-target/release/hapax}
RES=results/ww_berlin_lineage
mkdir -p "$RES"
ROOT=$(pwd)
cd "$OUT"

lineage_score() { # dir name
  minimap2 -c -x asm20 --secondary=no -t "$THREADS" bait_astro.fa "$1/contigs.fa" 2>/dev/null \
    | awk '$12>=0' > "$1/vs_lineages.paf"
  minimap2 -c -x asm10 -t "$THREADS" node37.fa "$1/contigs.fa" 2>/dev/null > "$1/vs_target.paf"
  awk -v n="$2" -v OFS='\t' '{t=$6; sub(/^WW_211013_20-22hmix_S5_L001_/,"",t); sub(/_length.*/,"",t);
      m[t]+=$10; b[t]+=$11; k[t]++; ind=$0; gsub(/.*cg:Z:/,"",ind); while (match(ind,/[0-9]+[ID]/)) {
      x=substr(ind,RSTART,RLENGTH); if (x+0>=10) big[t]+=x+0; ind=substr(ind,RSTART+RLENGTH)}
      L[t]+=$9-$8}
    END{for (t in m) print n, t, k[t], L[t], sprintf("%.2f%%", 100*m[t]/b[t]), big[t]+0}' "$1/vs_lineages.paf" \
    | sort -k2,2 >> "$ROOT/$RES/lineages.tsv"
  awk -v n="$2" -v OFS='\t' '{m+=$10; b+=$11; k++} END{print n, k, sprintf("%.2f%%", 100*m/b)}' \
    "$1/vs_target.paf" >> "$ROOT/$RES/target.tsv"
  awk -v n="$2" -v OFS='\t' '{t=$6; sub(/^WW_211013_20-22hmix_S5_L001_/,"",t); sub(/_length.*/,"",t);
      if (!(($1,t) in seen)) {seen[$1,t]=1; c[$1]++; s[$1]=s[$1]" "t}}
    END{for (q in c) if (c[q]>1) print n, q, s[q]}' "$1/vs_lineages.paf" >> "$ROOT/$RES/multi_lineage.tsv"
}

printf "run\tlineage\tblocks\tcovered\tidentity\tindel_bases_ge10\n" > "$ROOT/$RES/lineages.tsv"
printf "run\tblocks\tidentity\n" > "$ROOT/$RES/target.tsv"
printf "run\tcontig\tlineages\n" > "$ROOT/$RES/multi_lineage.tsv"
declare -A T
run() { # name args...
  local name=$1; shift
  local t0=$(date +%s.%N)
  "$ROOT/$HAPAX" assemble -i astro_reads_trim.fa -o "hxl_$name" -t "$THREADS" --min-contig 300 "$@" 2>/dev/null
  T[$name]=$(echo "$(date +%s.%N) - $t0" | bc)
  lineage_score "hxl_$name" "$name"
}
run default
run split           --split-lineage 1
run split_div0.06   --split-lineage 1 --lineage-split-div 0.06
run split_div0.09   --split-lineage 1 --lineage-split-div 0.09
run split_nt        --split-lineage 1 --mode nucleotide
mkdir -p megahit_trim_lin && cp megahit_trim/final.contigs.fa megahit_trim_lin/contigs.fa
lineage_score megahit_trim_lin megahit_trimmed

cat > "$ROOT/$RES/manifest.json" <<JSON
{"script": "experiments/ww_berlin_lineage.sh", "git_sha": "$(git -C "$ROOT" rev-parse HEAD)",
 "input": "data/ww_berlin/astro_reads_trim.fa (SRR23973193 recruited astrovirus reads, mate 2 trimmed 14 bp)",
 "comparators": "authors' SPAdes astrovirus contigs (bait_astro.fa), target NODE_37",
 "seed": null, "threads": $THREADS,
 "wall_s": {$(for k in "${!T[@]}"; do printf '"%s": %s, ' "$k" "${T[$k]}"; done | sed 's/, $//')}}
JSON
cat "$ROOT/$RES/target.tsv" "$ROOT/$RES/lineages.tsv" "$ROOT/$RES/multi_lineage.tsv"
