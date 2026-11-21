#!/bin/bash
# The whole cross-sample recurrence run. No SSH, no SSM, no inbound access: this script is
# everything that will ever happen on the instance, so a mistake in it costs another launch.
# Tokens are upper-case names between double at-signs; every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, from Imo's launch guide and unchanged from the big-index
# run that proved it: the dead man's switch first, before anything that can fail; a
# heartbeat to S3 because console output is unreliable; one exit helper on every path that
# cancels the pending timer first, since `shutdown -h now` is a no-op while one is
# scheduled; python3.12 by name because this AMI's python3 is 3.9; and `test` on every
# input, because `aws s3 cp` over a missing key can leave a run scanning nothing.
#
# WHAT IS DIFFERENT FROM THE BIG-INDEX RUN. No database is staged: the subjects are SRA
# objects fetched one at a time and deleted after each scan, so the volume is sized by the
# largest single library rather than by a corpus. The searcher is this repo's own Rust tool,
# built here from the payload rather than shipped as a binary, because the container's glibc
# is newer than this AMI's and a binary built there dies on the instance with a loader error
# that looks nothing like its cause.
set -xeuo pipefail
mkdir -p /opt/hapax/out
exec > >(tee /dev/console /opt/hapax/log-userdata.txt) 2>&1
echo "=== hapax recurrence run starting $(date -Is) commit @@COMMIT@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/recurrence-$IID"
OUT="$S3/results/recurrence/$IID"

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --include "*.json" --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --region "$REGION" || true
  shutdown -c || true
  shutdown -h now || systemctl --no-block poweroff || halt -p || true
  exit "${2:-1}"
}

# A hard stop powers the machine off from outside this script, so `die` never runs for it.
# This loop is what makes a run that never finishes diagnosable, and it is armed before
# anything slow for the same reason the timer is.
(
  while true; do
    sleep @@HEARTBEAT_SECONDS@@
    df -h / > /opt/hapax/log-disk.txt 2>&1 || true
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --include "*.json" --region "$REGION" >/dev/null 2>&1 || true
    aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --region "$REGION" >/dev/null 2>&1 || true
  done
) &

dnf install -y tar gzip xz cargo >/dev/null || die "dnf could not install the build tools"

PYBIN=""
for cand in python3.12 python3.11; do
  dnf install -y "$cand" "$cand-pip" >/dev/null 2>&1 || true
  if command -v "$cand" >/dev/null 2>&1; then PYBIN="$cand"; break; fi
done
[ -n "$PYBIN" ] || die "no python 3.11+ available; nothing to run"
"$PYBIN" -m pip install --quiet numpy || die "numpy did not install"
echo "=== interpreter: $PYBIN ($($PYBIN -V)) ==="

# ---------------------------------------------------------------- the SRA toolkit, pinned
# Pinned because `fasterq-dump` is the one step whose output format this run depends on, and
# "latest" is not a version anyone can cite six months from now.
cd /opt/hapax
curl -fsSL --retry 5 -o sratoolkit.tar.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/@@SRA_VERSION@@/sratoolkit.@@SRA_VERSION@@-centos_linux64.tar.gz" \
  || die "the SRA toolkit @@SRA_VERSION@@ did not download"
tar xzf sratoolkit.tar.gz
export PATH="/opt/hapax/sratoolkit.@@SRA_VERSION@@-centos_linux64/bin:$PATH"
# Without a settings file the toolkit prompts, and a prompt on a machine nobody can log
# into is a hang that the dead-man timer eventually bills for in full.
mkdir -p /root/.ncbi
printf '/LIBS/GUID = "@@VDB_GUID@@"\n/libs/cloud/report_instance_identity = "false"\n' \
  > /root/.ncbi/user-settings.mkfg
fasterq-dump --version || die "fasterq-dump will not run"

# ---------------------------------------------------------------- code, queries, panel
aws s3 cp "$S3/code/@@CODE_KEY_NAME@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
test -f experiments/recurrence_pass.py || die "the payload has no recurrence_pass.py"
test -f experiments/kmer_coverage/Cargo.toml || die "the payload has no kmer_coverage crate"

( cd experiments/kmer_coverage && cargo build --release --offline 2>/dev/null \
    || cargo build --release ) || die "kmer_coverage did not build on the instance"
TOOL=/opt/hapax/experiments/kmer_coverage/target/release/kmer_coverage
test -x "$TOOL" || die "kmer_coverage built nothing executable"

aws s3 cp "$S3/scratch/recurrence-@@QUERY_SET@@/queries.tsv" /opt/hapax/queries.tsv \
  --region "$REGION" || die "the query set did not arrive"
aws s3 cp "$S3/scratch/recurrence-@@QUERY_SET@@/panel.tsv" /opt/hapax/panel.tsv \
  --region "$REGION" || die "the panel did not arrive"
test -s /opt/hapax/queries.tsv || die "the query set is empty"
cp /opt/hapax/queries.tsv /opt/hapax/out/queries_used.tsv

# ---------------------------------------------------------------- which subjects, in which order
# SPREAD ACROSS THE PANEL, NOT THE FIRST N. The panel is sorted by date, so a pilot that
# took the first eight would measure eight libraries from one fortnight of 2021 - and
# `results/dna_rna_pair/` already showed this plant's composition moving 29 points across
# eight months, so the pilot would answer a seasonal question by accident. An even stride
# over the dates gives a pilot whose recurrence figure is at least about the same thing the
# full pass will be about.
"$PYBIN" - <<'PY' > /opt/hapax/subjects.txt || die "the subject list could not be built"
import csv, sys
rows = [r for r in csv.DictReader(open("/opt/hapax/panel.tsv"), delimiter="\t")
        if r["capture_enriched"] == "0" and r["paired"] == "1"
        and r["run"] != "@@QUERY_RUN@@"]
rows.sort(key=lambda r: (r["date"], r["window"], r["prep"]))
n = @@SUBJECT_COUNT@@
if n and n < len(rows):
    stride = len(rows) / n
    rows = [rows[int(i * stride)] for i in range(n)]
for r in rows:
    print(r["run"])
PY
test -s /opt/hapax/subjects.txt || die "no subjects were selected"
wc -l < /opt/hapax/subjects.txt
cp /opt/hapax/subjects.txt /opt/hapax/out/subjects.txt

# ---------------------------------------------------------------- the pass, in shards
# WHY SHARDS AND NOT THREADS. The scan is one single-threaded sweep of a library's reads,
# and `fasterq-dump` is the only threaded step, so one subject at a time leaves most of
# eight vCPU idle. Shards are whole subjects handed to separate processes, which parallelises
# both steps without touching the tool. Each shard gets its own work directory, because they
# delete their own FASTQ as they go and a shared directory would have one shard's cleanup
# remove another's input.
#
# ROUND ROBIN, NOT CONTIGUOUS BLOCKS. The subject list is ordered by date, so contiguous
# blocks would put one season in each shard and a shard that died would remove a season
# rather than a sample of the panel. Interleaving makes a partial run a thinner panel
# instead of a biased one.
mkdir -p /opt/hapax/work
"$PYBIN" - <<PY || die "the shards could not be built"
from pathlib import Path
subs = Path("/opt/hapax/subjects.txt").read_text().split()
n = max(1, @@SHARDS@@)
for i in range(n):
    Path(f"/opt/hapax/shard{i}.txt").write_text("\n".join(subs[i::n]) + "\n")
    print(f"shard {i}: {len(subs[i::n])} subjects")
PY

set +e
PIDS=""
for i in $(seq 0 $(( @@SHARDS@@ - 1 ))); do
  mkdir -p "/opt/hapax/work/shard$i"
  "$PYBIN" experiments/recurrence_pass.py \
    --queries /opt/hapax/queries.tsv \
    --subjects $(cat "/opt/hapax/shard$i.txt") \
    --trim-r2 @@TRIM_R2@@ --threads @@SHARD_THREADS@@ \
    --tool "$TOOL" --work "/opt/hapax/work/shard$i" --out /opt/hapax/out \
    --timings-name "pass_timings_shard$i.json" \
    --budget-hours @@BUDGET_HOURS@@ > "/opt/hapax/log-pass$i.txt" 2>&1 &
  PIDS="$PIDS $!"
done
RC=0
for p in $PIDS; do wait "$p" || RC=$?; done
set -e
for i in $(seq 0 $(( @@SHARDS@@ - 1 ))); do
  tail -15 "/opt/hapax/log-pass$i.txt" || true
  cp "/opt/hapax/log-pass$i.txt" /opt/hapax/out/ || true
done
# A shard that failed is reported and does not discard the shards that worked: the output is
# one npz per subject, so a short panel is a smaller denominator and not a broken table.
[ "$RC" -eq 0 ] || echo "=== at least one shard exited $RC; see the per-shard logs ==="
ls /opt/hapax/out/witness_*.npz | wc -l > /opt/hapax/out/subjects_scored.txt
[ -s /opt/hapax/out/subjects_scored.txt ] || die "no subject was scored"

die "the recurrence pass finished" 0
