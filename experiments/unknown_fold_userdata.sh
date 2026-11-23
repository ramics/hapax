#!/bin/bash
# The whole structure-on-the-unnameable-pile run. No SSH, no SSM, no inbound access, so this
# script is the only thing that will ever happen on the instance and a second attempt costs
# another launch. It lives in the repo rather than inside the launcher's Python because it is
# the artefact worth reviewing before launch; `--print-run-script` renders it. Tokens are
# upper-case names between double at-signs, so every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, from Imo's launch guide and unchanged from
# `reach_userdata.sh`: the dead man's switch first, before anything that can fail; a
# heartbeat to S3, because console output is unreliable; one exit helper on every path,
# cancelling the pending timer first because `shutdown -h now` is a no-op while one is
# scheduled; python3.12 by name, since this AMI's python3 is 3.9; and a decision on the
# FILES after every download, because `aws s3 cp --recursive` over an empty prefix exits 0
# and `foldseek databases` exits 0 on a 17-byte rate-limit page.
#
# WHAT THIS RUN BUYS THAT THE CONTAINER CANNOT, which is different from every other run in
# this repo. It is not the databases: the only staged thing is a 2.4 GB model and a public
# CDN download. It is the cycles. ProstT5 encodes about 50 residues a second on the
# container's four cores and this run has of the order of a million residues of ORFs to
# encode, so the container is most of a day and eight vCPU is a few hours. **The money buys
# the cycles, not the databases** - the first run in this repo for which that is true.
#
# NOTHING HERE IS GRADED. The reach run dropped an arm whose database failed to stage,
# because a six-arm table with an honest gap beat a relaunch. This run is one arm, so a
# staging failure has nothing to fall back to and every stage is fatal. An empty subject
# library or a missing null would read as "structure places nothing", which is the result
# this run is trying to measure honestly.
set -xeuo pipefail
mkdir -p /opt/hapax/out /opt/hapax/bin /opt/hapax/work
exec > >(tee /dev/console /opt/hapax/log-run.txt) 2>&1
echo "=== hapax unknown_fold run starting $(date -Is) commit @@COMMIT@@ ==="

# Re-armed here rather than inherited: the bootstrap's timer is against the bootstrap, and
# this script `exec`s into its place, so arming again is what keeps the bill bounded if the
# bootstrap's own timer is ever removed.
shutdown -c || true
shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/unknown_fold-$IID"
OUT="$S3/results/unknown_fold/$IID"
THREADS=@@THREADS@@
SHARDS=@@SHARDS@@

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  # Logs and the exit reason go up on every path, including the failures. An object in the
  # store is the only record; GetConsoleOutput has come back empty for the first ten minutes
  # after exactly the instance somebody needed to read.
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --region "$REGION" || true
  shutdown -c || true
  shutdown -h now || systemctl --no-block poweroff || halt -p || true
  exit "${2:-1}"
}

# A hard stop powers the machine off from outside this script, so neither the uploads in
# `die` nor anything else here runs. This loop is what makes a run that never finishes
# diagnosable at all, and it is armed before anything slow for the same reason the timer is.
(
  while true; do
    sleep @@HEARTBEAT_SECONDS@@
    { df -h /; echo; ls /opt/hapax/work/subjects/cif 2>/dev/null | wc -l; echo "cif files";
      ls -la /opt/hapax/work/ 2>/dev/null; } > /opt/hapax/log-progress.txt 2>&1 || true
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --region "$REGION" >/dev/null 2>&1 || true
  done
) &

dnf install -y tar gzip xz findutils >/dev/null || die "dnf could not install the basics"

PYBIN=""
for cand in python3.12 python3.11; do
  dnf install -y "$cand" "$cand-pip" >/dev/null 2>&1 || true
  if command -v "$cand" >/dev/null 2>&1; then PYBIN="$cand"; break; fi
done
[ -n "$PYBIN" ] || die "no python 3.11+ available; nothing to run"
echo "=== interpreter: $PYBIN ($($PYBIN -V)) ==="

"$PYBIN" -m pip install --quiet --upgrade pip >/dev/null 2>&1 || true
"$PYBIN" -m pip install --quiet pyrodigal==3.7.1 \
  || die "pyrodigal did not install; the ORFs cannot be called"

# ---------------------------------------------------------------- foldseek
# PINNED, and the pin is the whole arm: release 9-427df8a reads ProstT5 weights only as
# `model.safetensors` or `model.gguf` and SEGFAULTS on the f16 gguf this account holds,
# which on an instance surfaces as a bare exit 1 with no message and cost the reach run four
# launches. 10-941cd33 reads `prostt5-f16.gguf`, which is what that file is.
cd /opt/hapax
curl -fsSL --retry 5 -o foldseek.tar.gz \
  "https://github.com/steineggerlab/foldseek/releases/download/@@FOLDSEEK_VERSION@@/foldseek-linux-avx2.tar.gz" \
  || die "foldseek @@FOLDSEEK_VERSION@@ did not download"
tar xzf foldseek.tar.gz
cp /opt/hapax/foldseek/bin/foldseek /opt/hapax/bin/foldseek || die "no foldseek binary in the tarball"
FSVER=$(/opt/hapax/bin/foldseek version) || die "foldseek will not run"
echo "=== foldseek $FSVER ===" | tee /opt/hapax/log-foldseek.txt
[ "$FSVER" = "@@FOLDSEEK_COMMIT@@" ] \
  || die "foldseek reports $FSVER, not the pinned @@FOLDSEEK_COMMIT@@; an unpinned build is how the segfault came back"

# ---------------------------------------------------------------- the model
mkdir -p /opt/hapax/model
aws s3 cp "s3://@@BUCKET@@/@@PROSTT5_KEY@@" /opt/hapax/model/prostt5-f16.gguf \
  --region "$REGION" || die "the ProstT5 weights did not stage"
# Decided on the FILE, not the exit code, and on a size rather than existence: a truncated
# 2.4 GB copy is the failure that would read as the encoder producing nothing.
SZ=$(stat -c%s /opt/hapax/model/prostt5-f16.gguf)
[ "$SZ" -gt 2000000000 ] || die "the staged model is $SZ bytes, not the expected 2.4 GB"

# ---------------------------------------------------------------- code and queries
aws s3 cp "$S3/code/@@CODE_KEY@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
for f in experiments/unknown_fold.py experiments/shuffle_control.py experiments/provenance.py; do
  test -f "$f" || die "the payload is short: no $f"
done

mkdir -p /opt/hapax/piles
aws s3 cp "$S3/scratch/unknown-fold-piles-@@PILE_SET@@/" /opt/hapax/piles/ \
  --recursive --region "$REGION" || die "the contig piles did not copy"
for f in sample_contigs_unidentified.fa sample_contigs_identified_norrna.fa; do
  test -s "/opt/hapax/piles/$f" || die "no $f; the run would start on nothing"
  grep -c '^>' "/opt/hapax/piles/$f" >> /opt/hapax/log-piles.txt
done

# Imports before the long read. On the instance a short payload costs hours and a launch;
# locally it costs a second, and this is the check that catches the module added to a
# script's imports and not to the launcher's payload list.
"$PYBIN" - <<'PRECHECK' || die "the shipped payload cannot import"
import sys
sys.path.insert(0, "experiments")
import importlib.util
spec = importlib.util.spec_from_file_location("uf", "experiments/unknown_fold.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
import pyrodigal                                                     # noqa: F401
print("=== payload imports cleanly, bands:", m.BANDS, "===")
PRECHECK

# ---------------------------------------------------------------- the run
# One process, which does the library, the ORFs, the encoding and the search, and writes one
# manifest. `--shards` is passed through to the encoder: foldseek's own threading goes
# sublinear well before eight cores, so the encoding is split across independent processes
# instead, each with its own copy of the model in memory.
cd /opt/hapax
set +e
"$PYBIN" experiments/unknown_fold.py \
  --seeds "@@SEEDS@@" \
  --contigs-per-seed @@CONTIGS_PER_SEED@@ \
  --unidentified /opt/hapax/piles/sample_contigs_unidentified.fa \
  --identified /opt/hapax/piles/sample_contigs_identified_norrna.fa \
  --foldseek /opt/hapax/bin/foldseek \
  --model /opt/hapax/model \
  --work /opt/hapax/work \
  --out /opt/hapax/out \
  --threads "$THREADS" \
  --shards "$SHARDS" \
  > /opt/hapax/log-unknown-fold.txt 2>&1
RC=$?
set -e
tail -40 /opt/hapax/log-unknown-fold.txt || true
[ "$RC" -eq 0 ] || die "unknown_fold.py exited $RC; see log-unknown-fold.txt" "$RC"

test -s /opt/hapax/out/manifest.json || die "the run wrote no manifest, so nothing is traceable"
test -s /opt/hapax/out/hits.tsv.gz || die "the run wrote no hits file"

aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --region "$REGION" \
  || die "the results did not upload, which is the only copy"
die "finished cleanly" 0
