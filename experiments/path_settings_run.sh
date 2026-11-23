#!/bin/bash
# The exact invocations behind results/path_settings/. Same role as
# experiments/holdout_cascade_run.sh: the scripts it calls each write their own manifest
# with their own content hash, and this records the sequence and the arguments, which the
# manifests alone do not make readable.
#
# WHAT RUNS HERE AND WHY IN THIS ORDER
#
#   1. experiments/path_settings.py --stages assemble score
#      Ten whole-library assemblies of SRR23421544, one per path setting, then recovery and
#      purity of all 78 committed panel genomes against each one. This stage needs no index
#      at all: the ground truth is 25 genomes the depositors assembled out of this same
#      water and deposited, and the 53 panel genomes with zero reads mapped in this run are
#      the false-join control. So the decisive arm costs nothing but CPU.
#
#   2. the held-out indexes, rebuilt from the committed accession list
#      Identification has to be made against indexes rolled back to 2022-05-10, the day
#      this water was collected, or the experiment is decided in the subtrahend: today's
#      RefSeq contains the withheld genome itself. experiments/holdout_indexes.py applies
#      results/holdout_run/dropped_accessions.json and refuses to build unless six counts
#      reproduce. They do not, by nine nucleotide records and four proteins, because
#      upstream RefSeq viral keeps dropping records (19,641 served today against 19,650 in
#      PR #29 and 19,650 -> 19,202 kept after filtering). --allow-mismatch is therefore
#      passed and recorded, exactly as results/holdout_cascade/ did. DIRECTION: fewer
#      pre-cutoff relatives can only make a candidate look MORE novel, which is the
#      dangerous direction, so every per-arm identification is also run against today's
#      unfiltered index as the drift control.
#
#   3. experiments/classify_contigs.py, once per arm
#      Splits that arm's own assembly into what can be identified and what cannot. This is
#      the step that defines the residue, and running it per arm is the whole point of
#      Imo's objection: residue membership is an OUTPUT of the arm, never an input.
#
#   4. experiments/subtract_known.py, once per arm, five seeds
#      The read-level cascade with that arm's identified contigs as the own-assembly tier.
#      TWO TIERS ARE DELIBERATELY ABSENT and this is the one judgement call in the file:
#      the 8 GB Kraken2 Standard and PlusPF databases are not staged here. They are
#      identical across arms, so they cannot move a difference between arms, and
#      results/subtract_recut/ measures their marginal contribution at 0.27 and 0.00 points
#      against the own-contig tier's 30.33. So a residue fraction from this script is NOT
#      comparable to the 77.10% in results/subtract_known/ - it is comparable only across
#      the arms here, which is the only comparison this experiment makes.
set -euo pipefail

S="${S:?set S to the staging directory}"
OUT="${OUT:-results/path_settings}"
THREADS="${THREADS:-4}"
READS="${READS:-$S/reads/SRR23421544.fastq.gz}"
ARMS="${ARMS:-baseline,baseline_repeat,lax_bubble,no_bubble,strict_bubble,keep_thin,keep_local,prune_hard,no_rescue,meta_sensitive}"
# Five seeds are five disjoint 200k-read blocks of ONE pool, so the spread between them is
# subsampling noise and nothing else. BLOCK IS A PROPERTY OF THE LIBRARY: 1,000,000 puts
# the five blocks at 0, 1 M, 2 M, 3 M and 4 M of 4,290,429 reads. subtract_known.py's
# default of 1,500,000 is Berlin's stride and would make seeds 4 and 5 read nothing, which
# is the defect results/holdout_cascade/ found the hard way.
BLOCK="${BLOCK:-1000000}"
PAIRS="${PAIRS:-200000}"
SEEDS="${SEEDS:-1 2 3 4 5}"
CASCADE_ARMS="${CASCADE_ARMS:-$ARMS}"
K2BIN="${K2BIN:-$S/k2pkg/root/share/kraken2-2.17.2/libexec/classify}"

mkdir -p "$OUT" "$S/wk"

# ---- 1. the ladder and the ground-truth scoring --------------------------------------
python3 -I experiments/path_settings.py --seed 1 \
  --reads "$READS" \
  --megahit "$S/bin/megahit" --minimap2 "$S/bin/minimap2" \
  --threads "$THREADS" --memory 0.72 --arms "$ARMS" \
  --work "$S/wk" --out "$OUT"

# ---- 2. the held-out indexes ----------------------------------------------------------
if [ ! -s "$S/holdout/db_holdout_prot.pin" ] && [ ! -s "$S/holdout/db_holdout_prot.00.pin" ]; then
  python3 -I experiments/holdout_indexes.py --seed 1 \
    --nt-fa "$S/refseq/viral.1.1.genomic.fna.gz" \
    --prot-fa "$S/refseq/viral.1.protein.faa.gz" \
    --gbff "$S/refseq/viral.1.genomic.gbff.gz" \
    --blastbin "$S/bin" --out "$S/holdout" --allow-mismatch
fi
# DIAMOND databases from the same two FASTAs the BLAST databases were built from, so the
# protein tier and the protein drift control read the same sequence. No committed script
# built these before now; results/holdout_cascade/ took them as pre-existing paths.
for pair in "holdout_prot:$S/holdout/holdout_prot.fa" "full_prot:$S/refseq/viral.1.protein.faa.gz"; do
  name="${pair%%:*}"; fa="${pair#*:}"
  [ -s "$S/prot/$name.dmnd" ] || "$S/bin/diamond" makedb --in "$fa" \
    --db "$S/prot/$name.dmnd" --threads "$THREADS" >/dev/null
done
[ -s "$S/prot/sprot.dmnd" ] || "$S/bin/diamond" makedb --in "$S/db/uniprot_sprot.fasta.gz" \
  --db "$S/prot/sprot.dmnd" --threads "$THREADS" >/dev/null

# ---- 2b. the SILVA index classify_contigs.py reads ------------------------------------
# classify_contigs.py takes --silva-mmi and, if the file is absent, quietly labels nothing
# ribosomal - the exact failure shape results/fail_loudly_audit/ is about. The index is
# normally a by-product of subtract_known.py, which cannot run before identification here,
# so it is built FROM subtract_known's own two functions rather than from a local copy of
# them: a second SILVA cleaner or a different minimap2 preset would be a second instrument.
if [ ! -s "$S/wk/silva.mmi" ]; then
  python3 -I - "$S" "$THREADS" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, "experiments")
import subtract_known as sk
S, threads = Path(sys.argv[1]), int(sys.argv[2])
fa = S / "wk/silva_dna.fa"
if not fa.exists():
    sk._silva_to_dna([S / "db/SSURef_NR99_1381.fasta.gz",
                      S / "db/LSURef_NR99_1381.fasta.gz"], fa)
sk._mm_index(fa, S / "wk/silva.mmi", threads)
print("silva.mmi built from subtract_known._silva_to_dna + _mm_index")
PY
fi

# ---- 3. identification, one pass per arm ---------------------------------------------
# The header rewrite is the same pure relabelling holdout_cascade_run.sh uses: MEGAHIT
# writes `flag= multi= len=` where classify_contigs.py and subtract_known.py read length
# and coverage out of a SPAdes-style name. Rewriting the name keeps both committed scripts
# unmodified, where a second header parser would be a second instrument.
for arm in ${ARMS//,/ }; do
  src="$S/wk/asm_$arm/final.contigs.fa"
  dst="$S/wk/asm_$arm/contigs_renamed.fa"
  test -s "$src" || { echo "arm $arm has no assembly at $src" >&2; exit 2; }
  if [ ! -s "$dst" ]; then
    python3 -I - "$src" "$dst" "$arm" <<'PY'
import re, sys
src, dst, arm = sys.argv[1:4]
with open(src) as fh, open(dst, "w") as w:
    i = 0
    for line in fh:
        if line.startswith(">"):
            i += 1
            L = re.search(r"len=(\d+)", line)
            C = re.search(r"multi=([0-9.]+)", line)
            orig = line[1:].split()[0]
            w.write(f">{arm}_NODE_{i}_length_{L.group(1)}_cov_{C.group(1)}_megahit_{orig}\n")
        else:
            w.write(line)
print(f"{arm}: {i} contigs relabelled", file=sys.stderr)
PY
  fi
  if [ ! -s "$OUT/identification_$arm.json" ]; then
    python3 -I experiments/classify_contigs.py --seed 1 \
      --contigs "$dst" \
      --k2-db "$S/db/k2_viral_20210517.d" \
      --classify-bin "$K2BIN" \
      --dmnd "$S/prot/holdout_prot.dmnd" "$S/prot/sprot.dmnd" \
      --diamond-bin "$S/bin/diamond" \
      --silva-mmi "$S/wk/silva.mmi" \
      --threads "$THREADS" --work "$S/wk/contigs_$arm" --out "$S/wk/idout_$arm"
    cp "$S/wk/idout_$arm/contig_identification.json" "$OUT/identification_$arm.json"
    cp "$S/wk/idout_$arm/contig_identification_manifest.json" \
       "$OUT/identification_${arm}_manifest.json"
  fi
done

# ---- 4. the read-level cascade, one pass per arm per seed -----------------------------
for arm in ${CASCADE_ARMS//,/ }; do
  base="$S/wk/asm_$arm/contigs_renamed"
  for sd in $SEEDS; do
    [ -s "$OUT/cascade_${arm}_s$sd.json" ] && continue
    python3 -I experiments/subtract_known.py --seed "$sd" \
      --r1 "$READS" --r2 "$READS" --trim-r2 0 \
      --pairs "$PAIRS" --block "$BLOCK" \
      --k2-dbs "$S/db/k2_viral_20210517.d" --classify-bin "$K2BIN" \
      --no-k2-memory-map \
      --silva "$S/db/SSURef_NR99_1381.fasta.gz" "$S/db/LSURef_NR99_1381.fasta.gz" \
      --selfasm "${base}_rrna.fa" "${base}_identified_norrna.fa" \
      --observe "${base}_unidentified.fa" \
      --dmnd "$S/prot/holdout_prot.dmnd" "$S/prot/sprot.dmnd" \
      --diamond-bin "$S/bin/diamond" \
      --threads "$THREADS" --work "$S/wk/cascade_$arm" --out "$S/wk/cascout_${arm}_s$sd"
    cp "$S/wk/cascout_${arm}_s$sd/subtract_s$sd.json" "$OUT/cascade_${arm}_s$sd.json" 2>/dev/null \
      || cp "$S/wk/cascout_${arm}_s$sd"/*.json "$OUT/" 2>/dev/null || true
  done
done

echo "done: $OUT"
