#!/bin/bash
# Production end to end over a library this repo has never seen: stage everything, measure
# the one per-library fact nothing can inherit, then run experiments/prod_run.sh.
#
# No SSH, no SSM, no inbound access: this script is everything that will ever happen on the
# instance, so a mistake in it costs another launch. Tokens are upper-case names between
# double at-signs; every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, unchanged from crossplant_userdata.sh and the big-index run
# that proved it: the dead man's switch first, before anything that can fail; a heartbeat to
# S3 because console output is unreliable; one exit helper on every path that cancels the
# pending timer first, since `shutdown -h now` is a no-op while one is scheduled; python by
# name because this AMI's python3 is 3.9; and `test` on every input, because `aws s3 cp`
# over a missing key can leave a run reading nothing.
#
# WHAT IS DIFFERENT FROM EVERY EARLIER RUN HERE. This one stages the cascade's whole index
# set (SILVA, four Kraken2 databases, two protein indexes) and three binaries that no
# committed script has ever downloaded - Kraken2, MEGAHIT and minimap2 were hand-placed in
# a staging directory on somebody's machine, which is why `$S/k2pkg/root/share/...` appears
# as a bare path in nine scripts with nothing in the repo that builds it. Those three are
# the gap this file closes, and each is pinned.
#
# THE EXPENSIVE ORDER DECISION. The cheap, fail-fast things come first: the binaries, which
# are a few MB and prove they run, and the measurement that decides the trim, which needs
# only SILVA and 500,000 pairs. The assembly - hours, and the one step that cannot be
# resumed from nothing - comes after all of them. A launch that is going to die for a
# missing tool should die at minute ten, not at hour four.
set -xeuo pipefail
mkdir -p /opt/hapax/out
exec > >(tee /dev/console /opt/hapax/log-userdata.txt) 2>&1
echo "=== hapax prod run starting $(date -Is) commit @@COMMIT@@ for @@ACC@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/prod-@@ACC@@-$IID"
# The S3 side is OUTS3 throughout, because prod_run.sh takes OUT as its LOCAL output
# directory and one name for both would upload into a directory or write into a URL.
OUTS3="$S3/results/prod/@@ACC@@/$IID"
ACC="@@ACC@@"
THREADS="@@THREADS@@"

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --include "*.json" --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$OUTS3/" --recursive --region "$REGION" || true
  shutdown -c || true
  shutdown -h now || systemctl --no-block poweroff || halt -p || true
  exit "${2:-1}"
}

# A hard stop powers the machine off from outside this script, so `die` never runs for it.
# This loop is what makes a run that never finishes diagnosable, and it is armed before
# anything slow for the same reason the timer is. It is also what makes a timed-out run
# worth paying for: every finished seed is already in the store.
(
  while true; do
    sleep @@HEARTBEAT_SECONDS@@
    { df -h /; free -g; uptime; ls -la /opt/hapax/out/ 2>/dev/null; } \
      > /opt/hapax/log-machine.txt 2>&1 || true
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --include "*.json" --region "$REGION" >/dev/null 2>&1 || true
    aws s3 cp /opt/hapax/out/ "$OUTS3/" --recursive --region "$REGION" >/dev/null 2>&1 || true
  done
) &

# ---------------------------------------------------------------- system packages
# cargo, because two cascade stages are this repo's own Rust: `linkage` and `census` call
# experiments/kmer_coverage/target/release/{kmer_coverage,census}, and hapax preflights
# every selected stage, so a missing binary refuses the WHOLE run at the top rather than
# skipping those two. Built here from the payload rather than shipped as a binary, because
# the container's glibc is newer than this AMI's and a binary built there dies on the
# instance with a loader error that looks nothing like its cause.
dnf install -y tar gzip bzip2 xz unzip zstd pigz gcc gcc-c++ make zlib-devel cargo \
  >/dev/null || die "dnf could not install the build tools"

PYBIN=""
for cand in python3.12 python3.11; do
  dnf install -y "$cand" "$cand-pip" "$cand-devel" >/dev/null 2>&1 || true
  if command -v "$cand" >/dev/null 2>&1; then PYBIN="$cand"; break; fi
done
[ -n "$PYBIN" ] || die "no python 3.11+ available; nothing to run"
echo "=== interpreter: $PYBIN ($($PYBIN -V)) ==="

# ---------------------------------------------------------------- code, then its deps
cd /opt/hapax
aws s3 cp "$S3/code/@@CODE_KEY_NAME@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
for f in hapax.py requirements.txt experiments/prod_run.sh experiments/trim_mate2.py \
         experiments/trim_from_prefix.py experiments/megahit_rename.py \
         experiments/nontemplate_prefix_map.py experiments/classify_contigs.py \
         experiments/path_settings.py experiments/subtract_known.py \
         experiments/libraries.json experiments/stage_references.json; do
  test -f "/opt/hapax/$f" || die "the payload has no $f"
done
chmod +x experiments/prod_run.sh
# IN A VENV, and the first launch of this template is why: `python3.12 -m pip install -r
# requirements.txt` returned 0 on this AMI and left `import mappy` failing under
# `python3 -I`, which is how every script in this repo runs. The install went somewhere -I
# does not look, and --quiet hid it. A venv makes the interpreter and its packages one
# thing, and $PY is passed down to prod_run.sh so that no stage can reach the AMI's own
# python3, which is 3.9 and cannot run this code at all.
#
# One install, as TOOLS.md says a fresh container needs: deliberately the committed
# requirements file and not a hand-picked subset, because a stage that dies for a missing
# import three hours in costs more than the wheels do.
#
# zstandard is not a repo dependency and is not in requirements.txt: it is needed only to
# unpack the two conda packages below, on the instance, and declaring it in the repo's
# requirements would claim that something here imports it.
"$PYBIN" -m venv /opt/hapax/venv || die "the venv could not be created"
PY=/opt/hapax/venv/bin/python
"$PY" -m pip install --upgrade pip setuptools wheel > /opt/hapax/log-pip.txt 2>&1 \
  || die "pip itself could not be upgraded in the venv; see log-pip.txt"
"$PY" -m pip install zstandard >> /opt/hapax/log-pip.txt 2>&1 \
  || die "zstandard did not install; the conda packages cannot be unpacked"
"$PY" -m pip install -r requirements.txt >> /opt/hapax/log-pip.txt 2>&1 \
  || die "pip install -r requirements.txt failed; see log-pip.txt"
# One module at a time, so the die message names the module rather than the list.
for m in mappy numpy pyrodigal pyhmmer RNA sklearn; do
  "$PY" -I -c "import $m" || die "the venv cannot import $m; see log-pip.txt"
done
"$PY" -m pip freeze > /opt/hapax/out/pip_freeze.txt 2>&1 || true
echo "=== venv: $("$PY" -V) at $PY ==="

# ---------------------------------------------------------------- this repo's own Rust
( cd experiments/kmer_coverage && cargo build --release --offline 2>/dev/null \
    || cargo build --release ) || die "kmer_coverage did not build on the instance"
for b in kmer_coverage census; do
  test -x "/opt/hapax/experiments/kmer_coverage/target/release/$b" \
    || die "the crate built no $b; the linkage and census stages would refuse the run"
done

# ---------------------------------------------------------------- binaries, pinned
# Each one is proved to RUN here, not merely to download: a conda-built binary that cannot
# find a library on this AMI fails in a way that looks like a missing database later on.
mkdir -p /opt/hapax/bin

curl -fsSL --retry 5 -o diamond.tar.gz \
  "https://github.com/bbuchfink/diamond/releases/download/v@@DIAMOND_VERSION@@/diamond-linux64.tar.gz" \
  || die "DIAMOND @@DIAMOND_VERSION@@ did not download"
tar xzf diamond.tar.gz -C /opt/hapax/bin diamond || die "the DIAMOND tarball has no diamond"
/opt/hapax/bin/diamond version || die "DIAMOND will not run"

# Kraken2 from the bioconda package, because the cascade calls the `classify` binary
# directly rather than the wrapper, and that binary lives inside the package at
# share/kraken2-<version>/libexec/. Extracted under k2pkg/root so the path matches the one
# nine committed scripts already name.
curl -fsSL --retry 5 -o k2.conda "@@KRAKEN2_URL@@" || die "Kraken2 did not download"
mkdir -p /opt/hapax/k2pkg/root
"$PY" -I - /opt/hapax/k2.conda /opt/hapax/k2pkg/root <<'PY' || die "Kraken2 did not unpack"
import io, sys, tarfile, zipfile
src, dst = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(src) as z:
    name = [n for n in z.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst")][0]
    data = z.read(name)
try:
    import zstandard
except ImportError:
    sys.exit("zstandard is not importable; the package cannot be unpacked")
raw = zstandard.ZstdDecompressor().decompress(data, max_output_size=512 << 20)
with tarfile.open(fileobj=io.BufferedReader(io.BytesIO(raw)), mode="r|") as t:
    t.extractall(dst)
print("kraken2 unpacked")
PY
K2BIN="/opt/hapax/k2pkg/root/share/kraken2-@@KRAKEN2_VERSION@@/libexec/classify"
test -x "$K2BIN" || die "Kraken2's classify is not at $K2BIN"
"$K2BIN" 2>&1 | head -2 | grep -q "Usage" || die "Kraken2's classify will not run on this AMI"

for spec in "megahit @@MEGAHIT_URL@@" "minimap2 @@MINIMAP2_URL@@"; do
  name="${spec%% *}"; url="${spec#* }"
  curl -fsSL --retry 5 -o "$name.pkg" "$url" || die "$name did not download"
  mkdir -p "/opt/hapax/$name"
  case "$url" in
    *.conda) "$PY" -I - "/opt/hapax/$name.pkg" "/opt/hapax/$name" <<'PY' || die "unpack failed"
import io, sys, tarfile, zipfile, zstandard
src, dst = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(src) as z:
    name = [n for n in z.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst")][0]
    data = z.read(name)
raw = zstandard.ZstdDecompressor().decompress(data, max_output_size=512 << 20)
with tarfile.open(fileobj=io.BufferedReader(io.BytesIO(raw)), mode="r|") as t:
    t.extractall(dst)
PY
      ;;
    *.tar.bz2) tar xjf "/opt/hapax/$name.pkg" -C "/opt/hapax/$name" || die "$name unpack failed" ;;
    *) die "no unpacker for $url" ;;
  esac
  test -x "/opt/hapax/$name/bin/$name" || die "$name is not at /opt/hapax/$name/bin/$name"
done
MEGAHIT=/opt/hapax/megahit/bin/megahit
MINIMAP2=/opt/hapax/minimap2/bin/minimap2
"$MEGAHIT" --version || die "MEGAHIT will not run"
"$MINIMAP2" --version || die "minimap2 will not run"
{ "$MEGAHIT" --version; "$MINIMAP2" --version; /opt/hapax/bin/diamond version; } \
  > /opt/hapax/out/tool_versions.txt 2>&1

# ---------------------------------------------------------------- the SRA toolkit, pinned
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

# ---------------------------------------------------------------- the reads
# Staged where experiments/libraries.json says this library's reads live, so every script
# resolves them through the committed entry rather than through a path passed by hand.
READS=/opt/hapax/data/yeongdok
mkdir -p "$READS" /opt/hapax/sra
if [ ! -s "$READS/${ACC}_2.fastq.gz" ]; then
  prefetch --max-size 100G -O /opt/hapax/sra "$ACC" || die "prefetch $ACC failed"
  fasterq-dump --split-files --threads "$THREADS" --outdir "$READS" \
    --temp /opt/hapax/sra "/opt/hapax/sra/$ACC" || die "fasterq-dump $ACC failed"
  test -s "$READS/${ACC}_2.fastq" || die "$ACC did not come out as two mates; the entry in \
experiments/libraries.json says paired and the stages below would read one mate twice"
  pigz -p "$THREADS" "$READS/${ACC}_1.fastq" "$READS/${ACC}_2.fastq" || die "pigz failed"
  rm -rf /opt/hapax/sra
fi
ls -l "$READS" > /opt/hapax/log-reads.txt

# ---------------------------------------------------------------- SILVA, then the trim
# 138.1 Exports, the same two files and the same byte counts data/MANIFEST.md records for
# every committed ribosomal number: a different release would be a different operating
# point for the largest tier in the cascade.
mkdir -p /opt/hapax/data/silva
for u in "@@SILVA_SSU_URL@@" "@@SILVA_LSU_URL@@"; do
  f="/opt/hapax/data/silva/$(basename "$u")"
  [ -s "$f" ] || curl -fsSL --retry 5 -o "$f" "$u" || die "SILVA did not download: $u"
done
SSU="/opt/hapax/data/silva/$(basename "@@SILVA_SSU_URL@@")"
LSU="/opt/hapax/data/silva/$(basename "@@SILVA_LSU_URL@@")"
printf '%s %s\n' "$(stat -c%s "$SSU")" "$SSU" >> /opt/hapax/log-silva.txt
printf '%s %s\n' "$(stat -c%s "$LSU")" "$LSU" >> /opt/hapax/log-silva.txt
[ "$(stat -c%s "$SSU")" = "197133995" ] \
  || die "the SILVA SSU file is $(stat -c%s "$SSU") bytes, not the 197,133,995 every \
committed ribosomal number was measured against. A changed release is a changed operating \
point and this run would not be comparable to any of them"
[ "$(stat -c%s "$LSU")" = "67498420" ] || die "the SILVA LSU file is the wrong size"

# The DNA-converted copy the prefix measurement reads, built from subtract_known's own two
# functions: a second SILVA cleaner would be a second instrument.
cd /opt/hapax
if [ ! -s /opt/hapax/data/silva/silva_dna.fa ]; then
  "$PY" -I - "$SSU" "$LSU" <<'PY' || die "silva_dna.fa could not be built"
import sys
from pathlib import Path
sys.path.insert(0, "experiments")
import subtract_known as sk
sk._silva_to_dna([Path(p) for p in sys.argv[1:]], Path("data/silva/silva_dna.fa"))
print("silva_dna.fa written")
PY
fi

# ---- the one number nothing can inherit ----------------------------------------------
# SPEC 4.1: the mate-2 non-template prefix is per library and the system refuses to guess.
# Five seeds, because a mode that is not the same in every seed is not a length to cut, and
# experiments/trim_from_prefix.py refuses in that case rather than choosing. A refusal here
# stops the run BEFORE the expensive half, which is the whole point of measuring first.
mkdir -p /opt/hapax/out/nontemplate_prefix
for sd in @@PREFIX_SEEDS@@; do
  [ -s "/opt/hapax/out/nontemplate_prefix/nontemplate_prefix_s$sd.json" ] && continue
  "$PY" -I experiments/nontemplate_prefix_map.py --seed "$sd" \
    --pairs @@PREFIX_PAIRS@@ --libraries "$ACC" \
    --silva /opt/hapax/data/silva/silva_dna.fa --threads "$THREADS" \
    --work /opt/hapax/wk/ntprefix --out /opt/hapax/out/nontemplate_prefix \
    || die "the prefix measurement failed on seed $sd"
done
TRIM=$("$PY" -I experiments/trim_from_prefix.py --library "$ACC" \
  --prefix-json /opt/hapax/out/nontemplate_prefix/nontemplate_prefix_s*.json \
  --out /opt/hapax/out/nontemplate_prefix/trim_from_prefix.json) \
  || die "no single mate-2 trim can be named for $ACC; see \
out/nontemplate_prefix/trim_from_prefix.json. The run stops here rather than inheriting \
another library's prefix or picking a mode by eye"
echo "=== measured trim_r2 for $ACC: $TRIM ==="
echo "$TRIM" > /opt/hapax/out/trim_r2.txt

# ---------------------------------------------------------------- the index set
mkdir -p /opt/hapax/db /opt/hapax/prot
# Kraken2 databases, streamed out of the public genome-idx bucket so the volume never holds
# both the archive and the database. Names are the ones data/MANIFEST.md records, which is
# what makes this run's tier table readable beside the committed one.
for db in @@K2_DBS@@; do
  d="/opt/hapax/db/$db.d"
  if [ ! -s "$d/hash.k2d" ]; then
    mkdir -p "$d"
    curl -fsSL --retry 5 "https://genome-idx.s3.amazonaws.com/kraken/$db.tar.gz" \
      | tar xz -C "$d" || die "the Kraken2 database $db did not stage"
    test -s "$d/hash.k2d" || die "$db unpacked without a hash.k2d"
  fi
done
du -sh /opt/hapax/db/* > /opt/hapax/log-db.txt 2>&1 || true

# The two protein indexes of the cascade's last two tiers, in that order: RefSeq viral
# proteins then Swiss-Prot. `divergent` requires exactly these two and in this order.
if [ ! -s /opt/hapax/prot/viral.dmnd ]; then
  curl -fsSL --retry 5 -o /opt/hapax/prot/viral.1.protein.faa.gz \
    "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz" \
    || die "RefSeq viral proteins did not download"
  /opt/hapax/bin/diamond makedb --in /opt/hapax/prot/viral.1.protein.faa.gz \
    --db /opt/hapax/prot/viral.dmnd --threads "$THREADS" >/dev/null \
    || die "diamond makedb failed on the viral proteins"
fi
if [ ! -s /opt/hapax/prot/sprot.dmnd ]; then
  curl -fsSL --retry 5 -o /opt/hapax/prot/uniprot_sprot.fasta.gz \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
    || die "Swiss-Prot did not download"
  /opt/hapax/bin/diamond makedb --in /opt/hapax/prot/uniprot_sprot.fasta.gz \
    --db /opt/hapax/prot/sprot.dmnd --threads "$THREADS" >/dev/null \
    || die "diamond makedb failed on Swiss-Prot"
fi

# ---------------------------------------------------------------- the run itself
# Every sample-specific value is passed, none defaulted: the accession, the measured trim,
# the block sized to this library's read count, and the index set.
S=/opt/hapax/stage
mkdir -p "$S"
ACC="$ACC" S="$S" OUT="/opt/hapax/out/prod_$ACC" \
  R1="$READS/${ACC}_1.fastq.gz" R2="$READS/${ACC}_2.fastq.gz" \
  TRIM_R2="$TRIM" BLOCK="@@BLOCK@@" PAIRS="@@PAIRS@@" SEEDS="@@SEEDS@@" PY="$PY" \
  REPORT_PAIRS="@@REPORT_PAIRS@@" THREADS="$THREADS" MEMORY="@@MEGAHIT_MEMORY@@" \
  STAGES="@@STAGES@@" \
  MEGAHIT="$MEGAHIT" MINIMAP2="$MINIMAP2" DIAMOND=/opt/hapax/bin/diamond K2BIN="$K2BIN" \
  K2_DBS="$(for db in @@K2_DBS@@; do printf '/opt/hapax/db/%s.d ' "$db"; done)" \
  K2_DB_CONTIGS="/opt/hapax/db/@@K2_DB_CONTIGS@@.d" \
  SILVA="$SSU $LSU" \
  DMND="/opt/hapax/prot/viral.dmnd /opt/hapax/prot/sprot.dmnd" \
  AFTER_SEED='aws s3 cp /opt/hapax/out/ '"$OUTS3"'/ --recursive --region '"$REGION"' || true' \
  bash experiments/prod_run.sh || die "prod_run.sh exited non-zero; its manifests say where" 1

# The assembly itself, so a later pass over this sample never re-assembles: it is the
# expensive half and nothing about it depends on the stages that follow.
cp -f "$S/wk/asm_baseline/sample_contigs.fa" /opt/hapax/out/ 2>/dev/null || true
gzip -f /opt/hapax/out/sample_contigs.fa 2>/dev/null || true

die "finished" 0
