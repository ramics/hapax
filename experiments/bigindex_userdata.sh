#!/bin/bash
# The whole run. There is no SSH, no SSM and no inbound access to this instance, so this
# script is the only thing that will ever happen on it and there is no second attempt that
# does not cost another launch and another approval.
#
# It lives in the repo rather than inside the launcher's Python because it is the artefact
# worth reviewing before launch: `experiments/ec2_bigindex.py --print-user-data` renders it
# with its tokens filled in, and that rendering is what gets read line by line. A token is
# an upper-case name between double at-signs rather than Python's {braces}, so that every
# brace below is bash's own and the renderer can refuse a token it failed to replace.
#
# ORDER IS THE SAFETY PROPERTY, from Imo's launch guide:
#   1. the dead man's switch first, before anything that can fail;
#   2. a heartbeat to S3, because console output is unreliable and was empty or truncated
#      every time the reference project needed it;
#   3. one exit helper on every path, which cancels the pending timer before shutting down
#      because `shutdown -h now` is a no-op while a timer is scheduled;
#   4. python3.12 by name, since this AMI's python3 is 3.9;
#   5. `test -f` every input, because `aws s3 cp --recursive` over an empty prefix exits 0
#      and the workers would then start on nothing.
set -xeuo pipefail
mkdir -p /opt/hapax/out
exec > >(tee /dev/console /opt/hapax/log-userdata.txt) 2>&1
echo "=== hapax big-index run starting $(date -Is) commit @@COMMIT@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/bigindex-$IID"

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  # Logs and the exit reason go up on every path, including the failures. An object in the
  # store is the only record; GetConsoleOutput has come back empty for the first ten
  # minutes after exactly the instance somebody needed to read.
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --include "*.json" --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$S3/results/bigindex/$IID/" --recursive \
    --exclude "*" --include "rows_*.jsonl" --include "timing_*.json" \
    --region "$REGION" || true
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
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --include "*.json" --region "$REGION" >/dev/null 2>&1 || true
    aws s3 cp /opt/hapax/out/ "$S3/results/bigindex/$IID/" --recursive --exclude "*" \
      --include "rows_*.jsonl" --include "timing_*.json" --region "$REGION" \
      >/dev/null 2>&1 || true
    df -h / > /opt/hapax/log-disk.txt 2>&1 || true
  done
) &

stage_on() { case " @@STAGES@@ " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

dnf install -y tar gzip xz >/dev/null || die "dnf could not install the archive tools"

PYBIN=""
for cand in python3.12 python3.11; do
  dnf install -y "$cand" >/dev/null 2>&1 || true
  if command -v "$cand" >/dev/null 2>&1; then PYBIN="$cand"; break; fi
done
[ -n "$PYBIN" ] || die "no python 3.11+ available; nothing to run"
echo "=== interpreter: $PYBIN ($($PYBIN -V)) ==="

# ---------------------------------------------------------------- the searchers
# Both pinned. BLAST+ because the staged NCBI databases are BLAST volumes and nothing else
# reads them; DIAMOND at 2.1.10 specifically because that is the version every timing in
# `results/diamond_scaling/` was measured with, and the whole point of the third scaling
# point is that it be comparable to the first two.
cd /opt/hapax
curl -fsSL --retry 5 -o blast.tar.gz \
  "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/@@BLAST_VERSION@@/ncbi-blast-@@BLAST_VERSION@@+-x64-linux.tar.gz" \
  || die "BLAST+ @@BLAST_VERSION@@ did not download"
tar xzf blast.tar.gz && mv "ncbi-blast-@@BLAST_VERSION@@+" blast
curl -fsSL --retry 5 -o diamond.tar.gz \
  "https://github.com/bbuchfink/diamond/releases/download/v@@DIAMOND_VERSION@@/diamond-linux64.tar.gz" \
  || die "DIAMOND @@DIAMOND_VERSION@@ did not download"
tar xzf diamond.tar.gz
test -x /opt/hapax/diamond || die "the DIAMOND binary is not where it was expected"
/opt/hapax/blast/bin/blastn -version || die "BLAST+ will not run"
/opt/hapax/diamond version || die "DIAMOND will not run"

# ---------------------------------------------------------------- code and queries
aws s3 cp "$S3/code/@@CODE_KEY_NAME@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
test -f experiments/bigindex_search.py || die "the payload is short: no bigindex_search.py"

mkdir -p /opt/hapax/queries
aws s3 cp "$S3/scratch/bigindex-queries-@@QUERY_SET@@/" /opt/hapax/queries/ \
  --recursive --region "$REGION" || die "the query sets did not copy"
test -f /opt/hapax/queries/pilot.fa || die "no pilot.fa in the query sets"
ls -la /opt/hapax/queries/

# Imports before the long read. On the instance a short payload costs hours and a launch;
# here it costs a second.
"$PYBIN" - <<'PRECHECK' || die "the shipped payload cannot import"
import importlib.util
spec = importlib.util.spec_from_file_location("w", "experiments/bigindex_search.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
print("=== payload imports cleanly, fields:", m.FIELDS, "===")
PRECHECK

search() {  # label program db db_version query extra...
  local label="$1" program="$2" db="$3" ver="$4" query="$5"; shift 5
  "$PYBIN" experiments/bigindex_search.py --program "$program" --db "$db" \
    --db-version "$ver" --query "$query" --label "$label" --threads @@THREADS@@ \
    "$@" > "/opt/hapax/log-$label.txt" 2>&1 \
    || die "search $label failed; see log-$label.txt"
  aws s3 cp /opt/hapax/out/ "$S3/results/bigindex/$IID/" --recursive --exclude "*" \
    --include "rows_$label.jsonl" --include "timing_$label.json" --region "$REGION" || true
}

# Project a pass from a measured rate, and say so in the log either way. This is the
# guide's "check the first chunk against a measured baseline" with the only baseline that is
# honest here: no prior blastn-against-core_nt figure exists in this project, so the pilot
# measures the rate and this decides what still fits the time left.
#
# It does not kill the run when the whole set does not fit. A pass that cannot finish
# produces no table and spends the money anyway - but so does a run that refuses everything
# because block five would have overrun. So the caller uses this per block and stops adding
# blocks, which leaves a smaller table that is complete rather than a large one that is
# truncated, and the projection is in the log so the scope reduction is visible.
fits() {  # timing_label queries budget_hours
  "$PYBIN" - "$1" "$2" "$3" <<'FITS'
import json, sys
label, queries, budget = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
t = json.load(open(f"/opt/hapax/out/timing_{label}.json"))
rate = t["seconds_per_1k_queries"]
hours = rate * queries / 1000 / 3600
verdict = "fits" if hours <= budget else "DOES NOT FIT"
print(f"projection from {label}: {rate:.2f}s/1k x {queries:.0f} queries = {hours:.2f}h "
      f"against {budget:.2f}h -> {verdict}")
sys.exit(0 if hours <= budget else 3)
FITS
}

# ---------------------------------------------------------------- stage core_nt
if stage_on core_nt; then
  mkdir -p /opt/hapax/db
  # Unsigned, because the bucket is public and this account has no grant on it: a SIGNED
  # request is the one that fails here, with a bare AccessDenied that reads like a bucket
  # policy. That is the same trap as the regional endpoint, from the other side.
  aws s3 cp "s3://@@NCBI_BUCKET@@/@@NCBI_PREFIX@@/" /opt/hapax/db/ --recursive \
    --no-sign-request --exclude "*" --include "core_nt.*" \
    --include "core_nt-nucl-metadata.json" --region "$REGION" \
    || die "core_nt did not stage"
  test -f /opt/hapax/db/core_nt.nal || die "core_nt staged without its alias file"
  du -sh /opt/hapax/db > /opt/hapax/log-db.txt
  df -h / >> /opt/hapax/log-db.txt
fi

# ---------------------------------------------------------------- pilot
if stage_on pilot; then
  search pilot_core_nt blastn /opt/hapax/db/core_nt "@@NCBI_PREFIX@@" \
    /opt/hapax/queries/pilot.fa
  cat /opt/hapax/out/timing_pilot_core_nt.json
  # Reported, not enforced: the whole set not fitting is a scope answer, and the per-block
  # gate below is what acts on it. A pilot that cannot even project is a different thing,
  # and `search` has already died by then.
  fits pilot_core_nt @@FULL_QUERIES@@ @@FULL_BUDGET_HOURS@@ \
    || echo "=== the full @@FULL_QUERIES@@-query pass does not fit; blocks will be added while they fit ==="
fi

# ---------------------------------------------------------------- full core_nt pass
# One block at a time, each gated on the time left, so the run ends with whole blocks
# rather than a truncated one. A block is a self-contained measurement - its own residue
# and its own identified control - so three blocks is a three-seed result and this repo
# reports spread across what it actually ran.
if stage_on core_nt_full; then
  BUDGET=@@FULL_BUDGET_HOURS@@
  START=$SECONDS
  for q in /opt/hapax/queries/block_*.fa; do
    label="core_nt_$(basename "$q" .fa)"
    n=$(grep -c '^>' "$q" || echo 0)
    left=$("$PYBIN" -c "import sys;print(max(0.0,float(sys.argv[1])-float(sys.argv[2])/3600))" \
      "$BUDGET" "$((SECONDS - START))")
    if ! fits pilot_core_nt "$n" "$left"; then
      echo "=== stopping before $label: $n queries do not fit the ${left}h left ==="
      echo "stopped before $label with ${left}h left of a ${BUDGET}h budget" \
        > /opt/hapax/log-scope.txt
      break
    fi
    search "$label" blastn /opt/hapax/db/core_nt "@@NCBI_PREFIX@@" "$q"
  done
fi

# ---------------------------------------------------------------- nr, on a subsample
# nr is a BLAST protein database, so the only thing that reads it is blastx, which is far
# slower per read than DIAMOND. The question nr has to answer - does a 100 GB UniRef100
# index see what a 928 GB nr sees - is a comparison of identity DISTRIBUTIONS, and a few
# thousand reads per pile answers that at percentage resolution. So this pass is
# deliberately a subsample, which is also what makes "thousands per class instead of 60"
# true rather than aspirational.
if stage_on nr; then
  aws s3 cp "s3://@@NCBI_BUCKET@@/@@NCBI_PREFIX@@/" /opt/hapax/db/ --recursive \
    --no-sign-request --exclude "*" --include "nr.*" \
    --include "nr-prot-metadata.json" --region "$REGION" || die "nr did not stage"
  test -f /opt/hapax/db/nr.pal || die "nr staged without its alias file"
  du -sh /opt/hapax/db >> /opt/hapax/log-db.txt
  test -f /opt/hapax/queries/nr_sample.fa || die "no nr_sample.fa in the query sets"
  search nr_sample blastx /opt/hapax/db/nr "@@NCBI_PREFIX@@" \
    /opt/hapax/queries/nr_sample.fa
fi

echo "=== hapax big-index run finished $(date -Is) ==="
die "run complete" 0
