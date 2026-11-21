#!/bin/bash
# The whole run: no SSH, no SSM, no inbound access, so this script is the only thing that
# will ever happen on the instance and a second attempt costs another launch and another
# approval. It lives in the repo rather than in the launcher's Python because it is the
# artefact worth reviewing before launch; `--print-user-data` renders it. Tokens are
# upper-case names between double at-signs, so every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, from Imo's launch guide: the dead man's switch first,
# before anything that can fail; a heartbeat to S3, because console output is unreliable;
# one exit helper on every path, cancelling the pending timer first because `shutdown -h
# now` is a no-op while one is scheduled; python3.12 by name, since this AMI's python3 is
# 3.9; and `test -f` on every input, because `aws s3 cp --recursive` over an empty prefix
# exits 0 and the searches would start on nothing.
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
  echo "=== search $label starting $(date -Is) ==="
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
# TWO query counts against the same database, because one rate cannot be read. BLAST
# streams the whole subject database per invocation, so a tiny query set measures the SCAN
# and nothing else; the difference between it and the 500-read pilot is the marginal cost
# per read. That split is the number this run exists to produce - and it is also what sizes
# the chunks, since a chunk pays the scan again and a reclaim costs a chunk.
if stage_on pilot; then
  search scan_probe blastn /opt/hapax/db/core_nt "@@NCBI_PREFIX@@" \
    /opt/hapax/queries/scan_probe.fa
  search pilot_core_nt blastn /opt/hapax/db/core_nt "@@NCBI_PREFIX@@" \
    /opt/hapax/queries/pilot.fa
  "$PYBIN" - @@CHUNK_TARGET_SECONDS@@ > /opt/hapax/chunk.txt <<'SIZE'
import json, sys
target = float(sys.argv[1])
s = json.load(open("/opt/hapax/out/timing_scan_probe.json"))
p = json.load(open("/opt/hapax/out/timing_pilot_core_nt.json"))
scan = s["seconds"]
marginal = max((p["seconds"] - scan) / p["queries"], 1e-6)
chunk = int(max(100, min(8000, (target - scan) / marginal)))
print(chunk)
print(f"=== scan {scan:.0f}s on {s['queries']} queries; marginal {marginal:.3f}s/read "
      f"from {p['queries']} queries; chunk {chunk} to fill {target:.0f}s ===",
      file=sys.stderr)
SIZE
  CHUNK=$(head -1 /opt/hapax/chunk.txt)
  echo "=== chunk size $CHUNK ==="
  cp /opt/hapax/out/timing_scan_probe.json /opt/hapax/out/timing_pilot_core_nt.json /opt/hapax/ || true
fi
CHUNK="${CHUNK:-@@CHUNK@@}"

# ------------------------------------------------- the passes, in small chunks
# CHUNKS, NOT BLOCKS, AND SPOT IS THE REASON. Two instances were reclaimed in four hours on
# 2026-11-21 (`Server.SpotInstanceTermination`), the second of them 2.4 hours into a
# 6,000-read block, which left nothing: `search` uploads a block's rows when the block
# finishes, so a reclaim mid-block loses every read in it. The fix is the unit of work. A
# few hundred reads is ten to twenty minutes, so a reclaim now costs one chunk, and the
# chunk's rows are in the store before the next one starts.
#
# The chunks are cut from the block files, which interleave residue and identified control
# (bigindex_queries.py), so a run that stops early is a smaller sample of BOTH piles rather
# than all of one - the property that makes a partial run still a measurement.
split_blocks() {  # chunk_size
  mkdir -p /opt/hapax/queries/chunks
  "$PYBIN" - "$1" <<'SPLITB'
import sys
from pathlib import Path
size = int(sys.argv[1])
for src in sorted(Path("/opt/hapax/queries").glob("block_*.fa")):
    recs = [r for r in src.read_text().split(">") if r.strip()]
    stem = src.stem
    for i in range(0, len(recs), size):
        out = Path(f"/opt/hapax/queries/chunks/{stem}_c{i // size:03d}.fa")
        out.write_text("".join(">" + r for r in recs[i:i + size]))
    print(f"=== {stem}: {len(recs)} queries in {(len(recs) + size - 1) // size} chunk(s) ===")
SPLITB
}

# label_prefix budget_hours extra_args...
chunked_pass() {
  local tag="$1" budget="$2"; shift 2
  local start=$SECONDS first_label=""
  for q in /opt/hapax/queries/chunks/*.fa; do
    local label="${tag}_$(basename "$q" .fa)" n
    n=$(grep -c '^>' "$q" || echo 0)
    if [ -n "$first_label" ]; then
      local left
      left=$("$PYBIN" -c "import sys;print(max(0.0,float(sys.argv[1])-float(sys.argv[2])/3600))" \
        "$budget" "$((SECONDS - start))")
      if ! fits "$first_label" "$n" "$left"; then
        echo "=== stopping $tag before $label: ${left}h left of ${budget}h ==="
        echo "$tag stopped before $label with ${left}h left" >> /opt/hapax/log-scope.txt
        break
      fi
    fi
    search "$label" blastn /opt/hapax/db/core_nt "@@NCBI_PREFIX@@" "$q" "$@"
    [ -n "$first_label" ] || first_label="$label"
  done
}

if stage_on core_nt_full || stage_on core_nt_fast; then
  split_blocks "$CHUNK"
fi

# The sensitive task first: if only one of the two finishes, it is the one worth having.
if stage_on core_nt_full; then
  chunked_pass core_nt @@FULL_BUDGET_HOURS@@
fi

# `dc-megablast` on the SAME reads. It is built for cross-species nucleotide matching and is
# far cheaper, so running both says what the cheap one FAILS to find - this project's own
# cost-versus-completeness curve applied to itself.
if stage_on core_nt_fast; then
  chunked_pass fast @@FAST_BUDGET_HOURS@@ --blastn-task dc-megablast
fi

# ---------------------------------------------------------------- nr, on a subsample
# nr is a BLAST protein database, so only blastx reads it, far slower per read than
# DIAMOND. The question - does a 100 GB UniRef100 index see what a 928 GB nr sees - is a
# comparison of identity DISTRIBUTIONS, which a few thousand reads per pile settles at
# percentage resolution. So this pass is deliberately a subsample.
if stage_on nr; then
  aws s3 cp "s3://@@NCBI_BUCKET@@/@@NCBI_PREFIX@@/" /opt/hapax/db/ --recursive \
    --no-sign-request --exclude "*" --include "nr.*" \
    --include "nr-prot-metadata.json" --region "$REGION" || die "nr did not stage"
  test -f /opt/hapax/db/nr.pal || die "nr staged without its alias file"
  du -sh /opt/hapax/db >> /opt/hapax/log-db.txt
  test -f /opt/hapax/queries/nr_sample.fa || die "no nr_sample.fa in the query sets"
  # In chunks, gated the same way the blocks are. There is no measured blastx-against-nr
  # rate anywhere in this project, and the plausible range spans two orders of magnitude
  # per read, so a single call either finishes or leaves nothing. The first chunk measures
  # the rate; the rest are added while they fit. Chunks are interleaved residue and
  # control (bigindex_queries.py shuffles the two piles together), so a short run is a
  # smaller sample of both rather than all of one.
  mkdir -p /opt/hapax/queries/nrchunks
  "$PYBIN" - <<'SPLIT'
from pathlib import Path
src = Path("/opt/hapax/queries/nr_sample.fa").read_text().split(">")
recs = [r for r in src if r.strip()]
size = @@NR_CHUNK@@
for i in range(0, len(recs), size):
    out = Path(f"/opt/hapax/queries/nrchunks/nr_{i // size:03d}.fa")
    out.write_text("".join(">" + r for r in recs[i:i + size]))
print(f"=== {len(recs)} nr queries in {(len(recs) + size - 1) // size} chunk(s) ===")
SPLIT
  NRBUDGET=@@NR_BUDGET_HOURS@@
  NRSTART=$SECONDS
  first=1
  for q in /opt/hapax/queries/nrchunks/nr_*.fa; do
    label="nr_$(basename "$q" .fa)"
    n=$(grep -c '^>' "$q" || echo 0)
    if [ "$first" -eq 0 ]; then
      left=$("$PYBIN" -c "import sys;print(max(0.0,float(sys.argv[1])-float(sys.argv[2])/3600))" \
        "$NRBUDGET" "$((SECONDS - NRSTART))")
      if ! fits nr_nr_000 "$n" "$left"; then
        echo "=== stopping nr before $label: ${left}h left ==="
        echo "nr stopped before $label with ${left}h left" >> /opt/hapax/log-scope.txt
        break
      fi
    fi
    search "$label" blastx /opt/hapax/db/nr "@@NCBI_PREFIX@@" "$q"
    first=0
  done
fi

echo "=== hapax big-index run finished $(date -Is) ==="
die "run complete" 0
