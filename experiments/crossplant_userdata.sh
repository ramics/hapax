#!/bin/bash
# The whole cross-plant run (crossplant_sample.py, then crossplant_pass.py), adapted from
# recurrence_userdata.sh: same timer, heartbeat, exit helper and toolkit pin.
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
echo "=== hapax crossplant run starting $(date -Is) commit @@COMMIT@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/crossplant-$IID"
OUT="$S3/results/crossplant/$IID"

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
test -f experiments/crossplant_sample.py || die "the payload has no crossplant_sample.py"
test -f experiments/crossplant_pass.py || die "the payload has no crossplant_pass.py"
test -f experiments/crossplant_reduce.py || die "the payload has no crossplant_reduce.py"
test -f experiments/kmer_coverage/Cargo.toml || die "the payload has no kmer_coverage crate"

( cd experiments/kmer_coverage && cargo build --release --offline 2>/dev/null \
    || cargo build --release ) || die "kmer_coverage did not build on the instance"
TOOL=/opt/hapax/experiments/kmer_coverage/target/release/kmer_coverage
test -x "$TOOL" || die "kmer_coverage built nothing executable"

aws s3 cp "$S3/scratch/recurrence-s1-5/queries.tsv" /opt/hapax/berlin_queries.tsv \
  --region "$REGION" || die "the Berlin query set did not arrive"
aws s3 cp "$S3/scratch/crossplant-@@QUERY_SET@@/panel.tsv" /opt/hapax/panel.tsv \
  --region "$REGION" || die "the panel did not arrive"
test -s /opt/hapax/berlin_queries.tsv || die "the Berlin query set is empty"
test -s /opt/hapax/panel.tsv || die "the panel is empty"

# ---------------------------------------------------------------- subjects: ACC:TRIM_R2
# Berlin's SMARTer Pico mate 2 keeps its 14-nt trim; every other kit gets none. Libraries
# under the fixed depth are left out rather than scanned at a lower depth. A pilot takes an
# even stride over (project, date) so it touches every project.
"$PYBIN" - <<'PY' > /opt/hapax/subjects.txt || die "the subject list could not be built"
import csv
rows = [r for r in csv.DictReader(open("/opt/hapax/panel.tsv"), delimiter="\t")
        if r["paired"] == "1" and int(r["read_count"]) >= @@DEPTH@@ + 200000]
rows.sort(key=lambda r: (r["project"], r["date"], r["run"]))
n = @@SUBJECT_COUNT@@
if n and n < len(rows):
    stride = len(rows) / n
    rows = [rows[int(i * stride)] for i in range(n)]
for r in rows:
    print(f'{r["run"]}:{14 if r["project"] == "PRJNA948850" else 0}')
# Our own sample at the same depth, so its piles have an own-library abundance and it can be
# held out like any other library. It is not a panel member and never enters a union.
print("@@QUERY_RUN@@:14")
PY
test -s /opt/hapax/subjects.txt || die "no subjects were selected"
cp /opt/hapax/subjects.txt /opt/hapax/out/subjects.txt

# ---------------------------------------------------------------- phase 1: subsets + queries
mkdir -p /opt/hapax/work /opt/hapax/subsets /opt/hapax/out/sample
set +e
cat /opt/hapax/subjects.txt | xargs -P 6 -I{} bash -c '
  spec="{}"; acc="${spec%%:*}"; trim="${spec##*:}"
  '"$PYBIN"' experiments/crossplant_sample.py --acc "$acc" --trim-r2 "$trim" \
    --depth @@DEPTH@@ --work "/opt/hapax/work/$acc" --subjects /opt/hapax/subsets \
    --out /opt/hapax/out/sample >> /opt/hapax/log-sample.txt 2>&1 \
    || echo "SAMPLE FAILED $acc" >> /opt/hapax/log-sample.txt
  rm -rf "/opt/hapax/work/$acc"'
set -e
ls /opt/hapax/out/sample/sample_*.json | wc -l > /opt/hapax/out/subjects_sampled.txt
[ "$(cat /opt/hapax/out/subjects_sampled.txt)" -gt 0 ] || die "no subject was sampled"

# ---------------------------------------------------------------- one query file
# Berlin piles truncated to the same 90 nt per mate as every library query; shorter pairs
# are dropped and counted.
"$PYBIN" - <<'PY' || die "the query file could not be built"
import glob, json
L = 90
kept = dropped = 0
with open("/opt/hapax/queries.tsv", "w") as out:
    for line in open("/opt/hapax/berlin_queries.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f[1]) < L or len(f[2]) < L:
            dropped += 1
            continue
        out.write(f"berlin|{f[0]}\t{f[1][:L]}\t{f[2][:L]}\n")
        kept += 1
    libs = 0
    for p in sorted(glob.glob("/opt/hapax/out/sample/queries_*.tsv")):
        out.write(open(p).read())
        libs += 1
json.dump({"berlin_kept": kept, "berlin_dropped_short": dropped, "libraries": libs},
          open("/opt/hapax/out/query_build.json", "w"))
PY
wc -l < /opt/hapax/queries.tsv

# ---------------------------------------------------------------- phase 2: scan, in shards
ls /opt/hapax/out/sample/sample_*.json | sed 's#.*/sample_##; s#\.json##' > /opt/hapax/sampled.txt
"$PYBIN" - <<'PY' || die "the shards could not be built"
from pathlib import Path
ok = set(Path("/opt/hapax/sampled.txt").read_text().split())
subs = [s for s in Path("/opt/hapax/subjects.txt").read_text().split() if s.split(":")[0] in ok]
n = max(1, @@SHARDS@@)
for i in range(n):
    Path(f"/opt/hapax/shard{i}.txt").write_text("\n".join(subs[i::n]) + "\n")
    print(f"shard {i}: {len(subs[i::n])} subjects")
PY
set +e
PIDS=""
for i in $(seq 0 $(( @@SHARDS@@ - 1 ))); do
  "$PYBIN" experiments/crossplant_pass.py --queries /opt/hapax/queries.tsv \
    --subjects $(cat "/opt/hapax/shard$i.txt") --subject-dir /opt/hapax/subsets \
    --tool "$TOOL" --work "/opt/hapax/work/shard$i" --out /opt/hapax/out/witness \
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
cp /opt/hapax/log-sample.txt /opt/hapax/out/ || true
gzip -c /opt/hapax/queries.tsv > /opt/hapax/out/queries_used.tsv.gz || true
[ "$RC" -eq 0 ] || echo "=== at least one shard exited $RC ==="
ls /opt/hapax/out/witness/witness_*.npz | wc -l > /opt/hapax/out/subjects_scored.txt
"$PYBIN" experiments/crossplant_reduce.py --witness /opt/hapax/out/witness \
  --queries /opt/hapax/queries.tsv --out /opt/hapax/out/packed > /opt/hapax/log-reduce.txt 2>&1 \
  || echo "=== reduce failed; the per-subject npz files are still uploaded ==="
die "the crossplant pass finished" 0
