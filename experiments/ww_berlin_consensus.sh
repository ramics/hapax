#!/usr/bin/env bash
# Where hapax's 97.0% identity on the Berlin astrovirus comes from (results/ww_berlin_consensus).
#
# Run after experiments/ww_berlin.sh, which leaves the reads, the recruited astrovirus reads,
# the authors' 6,823 bp contig (node37.fa) and the read-to-bait BAM in data/ww_berlin.
#
# Three questions, each answered from the reads rather than from an identity number:
#  1. Is the per-column vote wrong? -> consensus_diag.py asks the reads at every disagreeing base.
#  2. Is it the library? -> 96% of mate-2 reads start with ~14 bp of non-template sequence at full
#     base quality (minimap2 soft-clips it). Measured below, then removed, and hapax re-run.
#  3. Is it lineage mixing? -> each read is labelled with its best bait contig (the authors' SPAdes
#     contigs: NODE_37 is the target lineage, NODE_1338/1845/1118/5252/8170 and references are
#     others), and every substitution is attributed to the lineage whose reads carry hapax's base.
#
# No seed: nothing here is stochastic.
set -euo pipefail

OUT=data/ww_berlin
THREADS=${THREADS:-4}
HAPAX=${HAPAX:-target/release/hapax}
TRIM=${TRIM:-14}
RES=results/ww_berlin_consensus
mkdir -p "$RES"
ROOT=$(pwd)
cd "$OUT"

echo "=== mate-2 5' soft-clip length (minimap2, reads vs baits): count by length"
samtools view -F 0x904 -f 128 astro.bam | python3 -c '
import sys, re, collections
c = collections.Counter()
for l in sys.stdin:
    f = l.split("\t"); rev = int(f[1]) & 16
    m = re.search(r"(\d+)S$", f[5]) if rev else re.match(r"^(\d+)S", f[5])
    c[int(m.group(1)) if m else 0] += 1
n = sum(c.values())
print("  reads", n, "clipped>=10:", sum(v for k, v in c.items() if k >= 10), "mode:", c.most_common(1)[0])
' | tee "$ROOT/$RES/r2_clip.txt"

# read -> lineage (best bait contig), for attribution
samtools view -F 0x904 astro.bam | awk -v OFS='\t' '{
  m = (int($2/64)%2) ? "/1" : "/2"; b = $3; sub(/^WW_211013_20-22hmix_S5_L001_/, "", b); sub(/_length.*/, "", b);
  print $1 m, b }' > read2bait.tsv

# recruited reads with the mate-2 prefix removed
awk -v t="$TRIM" '/^>/{r2 = ($0 ~ /\/2$/); print; next} {print (r2 ? substr($0, t + 1) : $0)}' \
  astro_reads.fa > astro_reads_trim.fa

score() {
  minimap2 -c -x asm10 -t "$THREADS" node37.fa "$1/contigs.fa" 2>/dev/null > "$1/vs_target.paf"
  printf "| %s | %s | %s | %s/6823 | %s | %s |\n" "$2" \
    "$(cut -f1 "$1/vs_target.paf" | sort -u | wc -l)" \
    "$(awk '{print $11}' "$1/vs_target.paf" | sort -n | tail -1)" \
    "$(awk '{print $8"\t"$9}' "$1/vs_target.paf" | sort -n | awk 'BEGIN{e=0;t=0}{if($1>e){t+=$2-$1;e=$2} else if($2>e){t+=$2-e;e=$2}} END{print t}')" \
    "$(awk '{m+=$10; b+=$11} END{printf "%.2f%%", 100*m/b}' "$1/vs_target.paf")" \
    "$(grep hmm_fraction "$1/stats.json" | tr -dc '0-9.')"
}

echo "=== hapax configurations vs the authors' contig"
echo "| run | contigs touching | longest block | covered | identity | HMM fraction |" | tee "$ROOT/$RES/table.md"
echo "|---|---|---|---|---|---|" | tee -a "$ROOT/$RES/table.md"
run() { # name input args...
  local name=$1 inp=$2; shift 2
  "$ROOT/$HAPAX" assemble -i "$inp" -o "hx_$name" -t "$THREADS" --min-contig 300 "$@" 2>/dev/null
  score "hx_$name" "$name" | tee -a "$ROOT/$RES/table.md"
}
run raw            astro_reads.fa
run raw_mode_nt    astro_reads.fa      --mode nucleotide
run raw_ident95    astro_reads.fa      --min-ident 0.95
run trim           astro_reads_trim.fa
run trim_mode_nt   astro_reads_trim.fa --mode nucleotide
run trim_ident90   astro_reads_trim.fa --min-ident 0.90
run trim_ident95   astro_reads_trim.fa --min-ident 0.95
run trim_ident95_nt astro_reads_trim.fa --min-ident 0.95 --mode nucleotide

echo "=== per-base attribution"
for name in raw trim trim_ident95; do
  inp=astro_reads.fa; [ "$name" = raw ] || inp=astro_reads_trim.fa
  python3 "$ROOT/experiments/consensus_diag.py" --contigs "hx_$name/contigs.fa" --target node37.fa \
    --reads "$inp" --labels read2bait.tsv --target-lineage NODE_37 --threads "$THREADS" \
    --out "$ROOT/$RES/diag_$name" > /dev/null
  rm -f "$ROOT/$RES/diag_$name"/*.bam*
done

echo "=== the target lineage alone (reads whose best bait is NODE_37)"
awk '$2=="NODE_37"{print $1}' read2bait.tsv > n37.ids
seqkit grep -f n37.ids astro_reads.fa 2>/dev/null > n37_reads.fa
seqkit grep -f n37.ids astro_reads_trim.fa 2>/dev/null > n37_reads_trim.fa
run n37_raw  n37_reads.fa
run n37_trim n37_reads_trim.fa
python3 "$ROOT/experiments/consensus_diag.py" --contigs hx_n37_raw/contigs.fa --target node37.fa \
  --reads n37_reads.fa --threads "$THREADS" --out "$ROOT/$RES/diag_n37_raw" > /dev/null
rm -f "$ROOT/$RES/diag_n37_raw"/*.bam*
