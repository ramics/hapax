#!/bin/bash
# The exact invocations behind results/holdout_cascade/. This file is the provenance of
# the run: the three committed scripts it calls each write their own manifest with their
# own content hash, and this records the arguments they were given, which those manifests
# alone do not make readable as a sequence.
#
# WHY A SHELL SCRIPT AND NOT A STAGE IN hapax.py. Two of the three steps ARE hapax.py (it
# is called once per seed below). The other two - splitting this library's own assembly
# into what can be identified and what cannot, and scoring the attributed contigs - sit
# either side of it, because `hapax run` deliberately takes an assembly as an INPUT and
# refuses to pretend it assembled one itself. So the ordering lives here.
#
# STAGE 0, the one worth reading. The sample's own identified contigs are the most
# efficient subtrahend this project has measured, and they are also the trap: if the
# withheld genome's own contig counts as "identified" it is subtracted as background in
# tier 1 and never reaches the part of the report anyone reads. That is exactly what this
# run is testing, so the split must be made by IDENTIFICATION against the date-held-out
# indexes and never by assembly. classify_contigs.py is the committed script that does it.
#
# SINGLE-END. This library has one mate. r2 is r1 again with --trim-r2 0, which is exact
# for every cascade tier because each is a min, max or union over the two mates; the
# per-read report is told to drop mate 2 rather than duplicate it, and the two
# mate-linkage stages refuse to run. See hapax.py --single-end.
set -euo pipefail

S="${S:?set S to the staging directory}"
OUT="${OUT:-results/holdout_cascade}"
THREADS="${THREADS:-4}"
PAIRS="${PAIRS:-200000}"
SEEDS="${SEEDS:-1 2 3 4 5}"
# hapax run's default is 20,000 residue pairs through the per-read report. MEASURED HERE:
# the `fold` field alone costs 290 ms per read (ViennaRNA, 20 dinucleotide shuffles per
# read for its own null), so 20,000 pairs is ~97 minutes per seed and ~8 hours over five -
# not runnable in this container, and the cost is in one field rather than in the search.
# 2,000 is the same per-arm sample size results/holdout_run/ used. The per-read half is
# not where this experiment's question lives: the withheld genome reaches a reader as a
# contig, and tier 3 here is read-level background. Recorded rather than silently lowered.
REPORT_PAIRS="${REPORT_PAIRS:-2000}"
K2BIN="$S/k2pkg/root/share/kraken2-2.17.2/libexec/classify"

mkdir -p "$OUT" "$S/wk"

# ---- 0a. MEGAHIT header -> the SPAdes-like form the committed scripts parse ------------
# classify_contigs.py and residue_report.py read length and coverage out of the contig
# NAME (SPAdes writes NODE_1_length_12710_cov_727.4). MEGAHIT writes the same two numbers
# as `flag=.. multi=.. len=..` instead. Renaming is a pure relabelling of the same contig
# and keeps both committed scripts unmodified; the alternative is a second header parser,
# which is a second instrument.
if [ ! -s "$S/asm/contigs_renamed.fa" ]; then
  python3 -I - "$S/asm/final.contigs.fa" "$S/asm/contigs_renamed.fa" <<'PY'
import re, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh, open(dst, "w") as w:
    i = 0
    for line in fh:
        if line.startswith(">"):
            i += 1
            L = re.search(r"len=(\d+)", line)
            C = re.search(r"multi=([0-9.]+)", line)
            orig = line[1:].split()[0]
            w.write(f">SRR23421544_NODE_{i}_length_{L.group(1)}_cov_{C.group(1)}"
                    f"_megahit_{orig}\n")
        else:
            w.write(line)
print(f"{i} contigs relabelled", file=sys.stderr)
PY
fi
CONTIGS="$S/asm/contigs_renamed.fa"

# ---- 0b. identified vs unidentifiable, against the HELD-OUT indexes -------------------
# The Kraken2 database here is the pre-cutoff viral one (2021-05-17), and the protein
# index is the date-filtered RefSeq viral release. Using today's indexes would decide the
# experiment's result in the subtrahend.
if [ ! -s "$OUT/contig_identification.json" ]; then
  python3 -I experiments/classify_contigs.py --seed 1 \
    --contigs "$CONTIGS" \
    --k2-db "$S/db/k2_viral_20210517.d" \
    --classify-bin "$K2BIN" \
    --dmnd "$S/prot/viral_holdout.dmnd" "$S/prot/sprot.dmnd" \
    --diamond-bin "$S/bin/diamond" --classify-bin "$K2BIN" \
    --silva-mmi "$S/wk/silva.mmi" \
    --threads "$THREADS" --work "$S/wk/contigs" --out "$OUT"
fi
# classify_contigs.py writes its splits beside the contig file it was given.
ID_RRNA="$S/asm/contigs_renamed_rrna.fa"
ID_NORRNA="$S/asm/contigs_renamed_identified_norrna.fa"
ID_ALL="$S/asm/contigs_renamed_identified.fa"
UNID="$S/asm/contigs_renamed_unidentified.fa"

# The per-read report now refuses to start when a field's reference is absent, rather than
# dropping the field and writing 0.0 for it. The committed run_s1-3 were made without the
# cross-sample assemblies, Rfam, the CRISPR spacers, Pfam and the chimaera reference
# (results/fail_loudly_audit/), so those five are left out BY NAME here and recorded as
# "not searched: --without".
# ---- 1. the committed cascade, one run directory per seed -----------------------------
# Five seeds are five disjoint 200k-read blocks of ONE read pool, so the spread between
# them is subsampling noise and nothing else. It is not library variation and it is not
# biological variation, and it is narrower than either by an unknown factor.
for sd in $SEEDS; do
  if [ -s "$OUT/run_s$sd/manifest.json" ]; then echo "seed $sd already done"; continue; fi
  python3 hapax.py run "$S/reads/SRR23421544.fastq.gz" "$S/reads/SRR23421544.fastq.gz" \
    --seed "$sd" --out "$OUT/run_s$sd" --work "$S/wk/hapax_shared" \
    --single-end --trim-r2 0 --pairs "$PAIRS" --threads "$THREADS" \
    --report-pairs "$REPORT_PAIRS" \
    --report-without cross_sample rfam crispr pfam chimaera \
    --k2-no-memory-map \
    --k2-dbs "$S/db/k2_viral_20210517.d" "$S/db/k2_standard08.d" "$S/db/k2_pluspf08.d" \
    --silva "$S/db/SSURef_NR99_1381.fasta.gz" "$S/db/LSURef_NR99_1381.fasta.gz" \
    --dmnd "$S/prot/viral_holdout.dmnd" "$S/prot/sprot.dmnd" \
    --diamond-bin "$S/bin/diamond" --classify-bin "$K2BIN" \
    --contigs "$ID_RRNA" "$ID_NORRNA" \
    --unidentified-contigs "$UNID" \
    --stages cascade contigs reads intervals layers strains \
    || echo "seed $sd exited non-zero; its manifest records which stage"
done

# ---- 2. what the output says about the withheld genome --------------------------------
# Resumable, because this step is deterministic given the assembly and the indexes and
# costs over an hour: `blastn -task blastn` of 22 Mb of contigs is the slow part, and
# -task blastn rather than megablast is not negotiable here - megablast would miss exactly
# the divergent matches the whole experiment is about.
if [ ! -s "$OUT/contig_rows_s1.tsv" ]; then
python3 -I experiments/holdout_cascade.py --seed 1 \
  --contigs "$CONTIGS" \
  --holdout-nt-db "$S/holdout/db_holdout_nt" \
  --holdout-prot-db "$S/prot/viral_holdout.dmnd" \
  --full-nt-db "$S/full/db_full_nt" \
  --full-prot-db "$S/prot/viral_full.dmnd" \
  --sprot-db "$S/prot/sprot.dmnd" \
  --blastbin "$S/bin" --diamond-bin "$S/bin/diamond" \
  --classify "$OUT/contig_identification.json" \
  --identified-fa "$ID_ALL" \
  --unidentifiable-fa "$UNID" \
  --run-dir "$OUT/run_s1" \
  --threads "$THREADS" --work "$S/wk/score" --out "$OUT"
fi

echo "done: $OUT"
