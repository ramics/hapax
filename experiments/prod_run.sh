#!/bin/bash
# Production end to end over one staged library: assemble, identify, then `hapax run`.
#
# WHAT THIS IS FOR. SPEC 13.1 says the command line is the contract and `hapax run` is the
# whole tool, and that is true of the part that measures. It is not the whole of a sample:
# `hapax run` deliberately takes an assembly as an INPUT and refuses to pretend it
# assembled one itself, so two steps sit either side of it - assembling the library, and
# splitting that assembly into what can be identified and what cannot. `holdout_cascade_run.sh`
# is the same ordering written for one library, with that library's facts in it. This file
# is that ordering with every sample-specific value taken from the environment and nothing
# defaulted, so a library the repo has never seen runs without a flag being edited.
#
# NO VALUE THAT IS A PROPERTY OF THE LIBRARY IS DEFAULTED, and the ones that are properties
# of the FILE are not asked for either: they are measured from it. The Berlin defaults
# silently cut 14 real bases off another kit's mate 2 and left two of five seeds of a
# smaller library reading nothing, and in both cases the output directory still looked
# right, so there is no fallback here for the same reason `library_facts.py` raises. But a
# refusal is not the goal: the pair count, the read lengths, the layout and the seed stride
# all come off one pass over the reads (`measure_library.py`), so a library nobody has
# touched needs no commit before it runs. What remains for a person is what reading the
# file cannot answer, and on the first new library that was exactly one thing: TRIM_R2, the
# non-template prefix, which needs the reads mapped onto something already known.
#
# THE ONE THING THIS SCRIPT CHECKS ITSELF is the block arithmetic, against the library's
# counted reads rather than a remembered count, because that failure is invisible: each
# seed fails loudly on its own and the directory still looks like a five-seed run. It is
# checked even when the stride was measured, because a stated stride overrides it.
set -xeuo pipefail

ACC="${ACC:?set ACC to the run accession, e.g. SRR23702504}"
S="${S:?set S to the staging directory}"
OUT="${OUT:?set OUT to the output directory}"
R1="${R1:?set R1 to mate 1 of the library}"
R2="${R2:?set R2 to mate 2 of the library}"
# No default: it is measured per library (experiments/nontemplate_prefix_map.py) and an
# inherited value is the bug this project has already paid for twice.
TRIM_R2="${TRIM_R2:?set TRIM_R2 to the measured non-template prefix of mate 2}"
# MEASURED, not required. The stride is arithmetic on the library's own pair count, so
# asking a caller for it turned "a sample nobody has touched" into "a sample nobody has
# touched plus a commit". Stating it still wins and is recorded as stated.
BLOCK="${BLOCK:-}"
PAIRS="${PAIRS:-200000}"
SEEDS="${SEEDS:-1 2 3 4 5}"
THREADS="${THREADS:-8}"
# The interpreter, because `python3` is not it everywhere: Amazon Linux
# 2023's python3 is 3.9 and every script here needs 3.11+, so an
# unqualified `python3` on an instance runs the wrong interpreter and dies
# inside a stage rather than at the top. Overridable, and `python3` locally,
# where python3 IS the interpreter.
PY="${PY:-python3}"
MEMORY="${MEMORY:-0.85}"
# `hapax run` defaults to 20,000 residue pairs through the per-read report. The `fold`
# field alone costs 290 ms/read (ViennaRNA, 20 dinucleotide shuffles per read for its own
# null), so 20,000 is ~97 min per seed and ~8 h over five. 2,000 is what
# results/holdout_cascade/ and results/holdout_run/ used, recorded rather than quietly
# lowered: it is the per-read arm's sample size, and the residue table's intervals widen
# with it.
REPORT_PAIRS="${REPORT_PAIRS:-2000}"
# `panel` and `calibration` are VALIDATION stages and are deliberately absent: both score
# the run against sequence whose distance was fixed before it - the depositors' own
# genomes from that water, and a spiked divergence ladder. Production water has no
# depositor, so a production run has nothing to score against and must not pretend the
# question applies (SPEC 13.4). They stay runnable by name on a fixture library.
STAGES="${STAGES:-cascade contigs reads linkage census intervals layers strains shortlist crossplant}"
MEGAHIT="${MEGAHIT:?set MEGAHIT to the megahit binary}"
MINIMAP2="${MINIMAP2:?set MINIMAP2 to the minimap2 binary}"
DIAMOND="${DIAMOND:?set DIAMOND to the diamond binary}"
K2BIN="${K2BIN:?set K2BIN to the Kraken2 classify binary}"
# The index set. Space-separated lists, so the operating point is visible in the log.
K2_DBS="${K2_DBS:?set K2_DBS to the Kraken2 databases}"
K2_DB_CONTIGS="${K2_DB_CONTIGS:?set K2_DB_CONTIGS to the one database the contig split uses}"
SILVA="${SILVA:?set SILVA to the SILVA SSU and LSU fasta files}"
DMND="${DMND:?set DMND to the DIAMOND protein databases}"
# Per-read fields whose reference this run does not stage. They are dropped BY NAME, which
# is recorded as "not searched", and never as 0.0 (results/fail_loudly_audit/). `chimaera`
# is deliberately NOT in this list: it reads against this sample's own assembly, which this
# run has, and CLAUDE.md requires a chimaera number beside every contiguity claim - beside
# its own false-positive rate, which for split-mapping is 39-46% and makes it a description
# of how the assembler joins genes rather than a filter.
REPORT_WITHOUT="${REPORT_WITHOUT:-cross_sample rfam crispr pfam}"

mkdir -p "$OUT" "$S/wk"
cd "$(dirname "$0")/.."            # the repo root, whatever the caller's cwd

# ---- 0. the blocks must fit the library, checked before anything runs -----------------
LAST_SEED=$(for s in $SEEDS; do echo "$s"; done | sort -n | tail -1)
NSEEDS=$(for s in $SEEDS; do echo "$s"; done | wc -l)
# One pass over both mates gives the pair count, the modal read length per mate, the
# layout and a stride that fits the seeds disjointly, written beside the run with
# provenance. The facts that are properties of the FILE are read off the file; the two
# that are not - the non-template prefix and what the water was - are named as such by
# the measurement itself and are the only things a person supplies.
if [ ! -s "$OUT/measured_library.json" ]; then
  "$PY" -I experiments/measure_library.py --r1 "$R1" --r2 "$R2" \
    --pairs-per-seed "$PAIRS" --seeds "$NSEEDS" --accession "$ACC" \
    --out "$OUT/measured_library.json" > /dev/null
fi
READS=$("$PY" -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["pairs"])' \
  "$OUT/measured_library.json")
if [ -z "$BLOCK" ]; then
  BLOCK=$("$PY" -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["block"])' \
    "$OUT/measured_library.json")
  echo "=== block stride $BLOCK measured from $READS pairs in this library ==="
else
  echo "=== block stride $BLOCK stated; the measurement is recorded beside it ==="
fi
echo "=== read length $("$PY" -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["read_length"])' "$OUT/measured_library.json") measured ==="
NEED=$(( (LAST_SEED - 1) * BLOCK + PAIRS ))
echo "library $READS pairs; $LAST_SEED seeds x $PAIRS pairs at stride $BLOCK needs $NEED"
if [ "$NEED" -gt "$READS" ]; then
  echo "REFUSING: seed $LAST_SEED would start at $(( (LAST_SEED - 1) * BLOCK )) in a" \
       "library of $READS pairs. Lower BLOCK to at most" \
       "$(( (READS - PAIRS) / (LAST_SEED - 1) )) or lower PAIRS." >&2
  exit 2
fi

# ---- 1. mate 2's non-template prefix comes off before the assembler sees it -----------
# Read-level stages trim as they stream; MEGAHIT reads files, so it needs a trimmed copy.
if [ ! -s "$S/asm_in/${ACC}_2.fq" ]; then
  "$PY" -I experiments/trim_mate2.py --r1 "$R1" --r2 "$R2" \
    --library "$ACC" --trim-r2 "$TRIM_R2" \
    --out-dir "$S/asm_in" --prefix "$ACC" \
    --manifest "$OUT/trim_manifest.json"
fi

# ---- 2. assemble the whole library, at the committed baseline setting -----------------
# Through path_settings.py rather than a bare megahit line, because that script is the one
# place in this repo that records the full argv, the MEGAHIT version and the timing of an
# assembly, and `results/path_settings/` measured the whole ladder of alternatives against
# this arm. `--arms baseline` is MEGAHIT's defaults plus `--min-contig-len 200`.
if [ ! -s "$S/wk/asm_baseline/final.contigs.fa" ]; then
  "$PY" -I experiments/path_settings.py --seed 1 \
    --stages assemble --arms baseline \
    --reads "$S/asm_in/${ACC}_1.fq" --r2 "$S/asm_in/${ACC}_2.fq" \
    --megahit "$MEGAHIT" --minimap2 "$MINIMAP2" \
    --threads "$THREADS" --memory "$MEMORY" \
    --no-truth-panel \
    --work "$S/wk" --out "$OUT/assembly"
fi
ASM="$S/wk/asm_baseline/final.contigs.fa"
test -s "$ASM" || { echo "the assembly is missing at $ASM" >&2; exit 2; }
# The trimmed copies are an intermediate of the assembly and nothing downstream reads them.
rm -f "$S/asm_in/${ACC}_1.fq" "$S/asm_in/${ACC}_2.fq"

# ---- 2a. MEGAHIT header -> the SPAdes-like form the committed scripts parse -----------
# NAMED sample_contigs.fa DELIBERATELY. crossplant_stage.py's `queries` arm
# requires `sample_contigs_unidentified.fa` and `sample_contigs_identified.fa`
# in --contig-dir BY NAME, because that is what the Berlin pass had: the
# depositors' assembly was staged as data/geo/sample_contigs.fa. Naming the
# relabelled assembly anything else makes a stage that would otherwise run
# refuse for a reason that has nothing to do with this sample. Found by
# running hapax's own preflight in the container before launching.
CONTIGS="$S/wk/asm_baseline/sample_contigs.fa"
if [ ! -s "$CONTIGS" ]; then
  "$PY" -I experiments/megahit_rename.py --contigs "$ASM" --out "$CONTIGS" --tag "$ACC" \
    --manifest "$OUT/rename_manifest.json"
fi

# ---- 2b. the SILVA index classify_contigs.py reads -----------------------------------
# Built from subtract_known.py's own two functions rather than from a copy of them: a
# second SILVA cleaner or a different minimap2 preset would be a second instrument reading
# the same quantity. The cascade builds the same index in its own work directory.
if [ ! -s "$S/wk/silva.mmi" ]; then
  "$PY" -I - "$S" "$THREADS" $SILVA <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, "experiments")
import subtract_known as sk
S, threads, silva = Path(sys.argv[1]), int(sys.argv[2]), [Path(p) for p in sys.argv[3:]]
fa = S / "wk/silva_dna.fa"
if not fa.exists():
    sk._silva_to_dna(silva, fa)
sk._mm_index(fa, S / "wk/silva.mmi", threads)
print("silva.mmi built from subtract_known._silva_to_dna + _mm_index")
PY
fi

# ---- 3. identified vs unidentifiable, by IDENTIFICATION and never by assembly ---------
# The sample's own identified contigs are the most efficient subtrahend this project has
# measured (+30.3 points) and also the trap: a divergent genome's own contig that counts as
# "identified" is subtracted as background in tier 1 and never reaches the part of the
# report anyone reads (results/holdout_cascade/). classify_contigs.py makes the split, and
# `removal_evidence.tsv` is what keeps the removal from being silent about how far away the
# thing it removed was (SPEC 4.5).
if [ ! -s "$OUT/contig_identification.json" ]; then
  "$PY" -I experiments/classify_contigs.py --seed 1 \
    --contigs "$CONTIGS" \
    --k2-db "$K2_DB_CONTIGS" --classify-bin "$K2BIN" \
    --dmnd $DMND --diamond-bin "$DIAMOND" \
    --silva-mmi "$S/wk/silva.mmi" \
    --threads "$THREADS" --work "$S/wk/contigs" --out "$OUT"
fi
ID_RRNA="$S/wk/asm_baseline/sample_contigs_rrna.fa"
ID_NORRNA="$S/wk/asm_baseline/sample_contigs_identified_norrna.fa"
UNID="$S/wk/asm_baseline/sample_contigs_unidentified.fa"
EVIDENCE="$OUT/removal_evidence.tsv"
for p in "$ID_RRNA" "$ID_NORRNA" "$UNID"; do
  test -f "$p" || { echo "classify_contigs.py wrote no $p" >&2; exit 2; }
done

# ---- 4. the cascade and the report, one run directory per seed ------------------------
# Five seeds are five disjoint blocks of ONE read pool, so the spread between them is
# subsampling noise: not library variation and not biological variation, and narrower than
# either by an unknown factor.
for sd in $SEEDS; do
  if [ -s "$OUT/run_s$sd/manifest.json" ]; then echo "seed $sd already done"; continue; fi
  "$PY" hapax.py run "$R1" "$R2" \
    --seed "$sd" --out "$OUT/run_s$sd" --work "$S/wk/hapax_shared" \
    --library "$ACC" --trim-r2 "$TRIM_R2" --pairs "$PAIRS" --block "$BLOCK" \
    --threads "$THREADS" --report-pairs "$REPORT_PAIRS" \
    --report-without $REPORT_WITHOUT \
    --k2-dbs $K2_DBS \
    --silva $SILVA \
    --dmnd $DMND \
    --diamond-bin "$DIAMOND" --classify-bin "$K2BIN" --minimap2 "$MINIMAP2" \
    --contigs "$ID_RRNA" "$ID_NORRNA" \
    --unidentified-contigs "$UNID" \
    --identified-contigs "$ID_NORRNA" \
    --rrna-contigs "$ID_RRNA" \
    --assembly "$CONTIGS" \
    --removal-evidence "$EVIDENCE" \
    --stages $STAGES \
    || echo "seed $sd exited non-zero; its manifest records which stage"
  # Uploaded as soon as it exists: a run that is cut short by its own timer should leave
  # every finished seed behind it, not nothing.
  if [ -n "${AFTER_SEED:-}" ]; then eval "$AFTER_SEED"; fi
done

echo "done: $OUT"
