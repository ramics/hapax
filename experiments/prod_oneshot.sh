#!/bin/bash
# The documented install and the documented command, on a library this repo has never
# seen. Nothing else.
#
# WHY THIS FILE EXISTS BESIDE prod_userdata.sh. That template is 400 lines, of which about
# 200 stage indexes and binaries by hand and then call a 253-line shell driver. Everything
# those 450 lines do is now a command in the package: `conda env create -f
# environment.yml && pip install .` is the install (PLAN.md 4), `hapax fetch-db --all` is
# the index set, and `hapax run --config run.toml --report` is the pass. PLAN.md is judged
# on the report a user gets from one install and one command on an input nobody here has
# looked at, and the only way to find out whether that is true is to run those commands
# and nothing else on a machine with none of it installed.
#
# So this template deliberately does NOT reimplement any of it. Where the old template has
# a pinned curl, this has the command that the README tells a reader to type. If a command
# cannot do what the hand-rolled version did, that is the finding, and it belongs in the
# run's report rather than in a workaround here.
#
# WHAT STILL HAPPENS BY HAND, and why each one is not a gap in the package:
#
#   - The reads. A library is fetched from the archive with the archive's own toolkit;
#     `hapax` takes FASTQ paths and has no business holding SRA credentials or a prefetch
#     cache.
#   - The mate-2 trim. It is a property of the sequencing run and the package refuses to
#     assume one (SPEC 4.1). The instrument that measures it needs ribosomal anchors, so
#     it runs here, after `fetch-db` has produced SILVA, and its answer is written into
#     the run configuration WITH ITS BASIS. A refusal stops the run before the expensive
#     half, which is the point of measuring first.
#   - The panel pointers. `experiments/libraries.json` has to say whether a truth panel
#     and a cross-plant panel apply to this library, because an unknown library and a
#     library with no panel are different facts and the second must not be read off the
#     first. That is the one committed file an unseen library still needs a line in, and
#     the run's report says so.
#
# No SSH, no SSM, no inbound access: this script is everything that will ever happen on
# the instance, so a mistake in it costs another launch. Tokens are upper-case names
# between double at-signs; every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, unchanged from prod_userdata.sh and the big-index run that
# proved it: the dead man's switch first, before anything that can fail; a heartbeat to S3
# because console output is unreliable; one exit helper on every path that cancels the
# pending timer first, since `shutdown -h now` is a no-op while one is scheduled; and
# `test` on every input, because `aws s3 cp` over a missing key can leave a run reading
# nothing.
set -xeuo pipefail
mkdir -p /opt/hapax/out
exec > >(tee /dev/console /opt/hapax/log-userdata.txt) 2>&1
echo "=== hapax one-command run starting $(date -Is) commit @@COMMIT@@ for @@ACC@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/oneshot-@@ACC@@-$IID"
OUTS3="$S3/results/oneshot/@@ACC@@/$IID"
ACC="@@ACC@@"
THREADS="@@THREADS@@"
# experiments/provenance.py raises rather than writing "unknown", correctly, and there is
# no git and no checkout here. The launcher refuses to ship a dirty tree, so this commit
# describes the code in the payload exactly.
export HAPAX_GIT_SHA="@@COMMIT@@"
# The package finds the stage scripts from the checkout it is installed from; this
# instance installs from an unpacked tarball, so it is told where rather than left to
# guess (src/hapax/cli.py refuses instead of guessing).
export HAPAX_REPO=/opt/hapax

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  # --exclude the environments: `--include "*.json"` over /opt/hapax walks every installed
  # package, and an earlier run wrote some 1,500 botocore paginator files into a bucket
  # with no delete.
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --include "run.toml" --include "*.json" \
    --exclude "mamba/*" --exclude "*/site-packages/*" --exclude "db/*" \
    --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$OUTS3/" --recursive --region "$REGION" || true
  shutdown -c || true
  shutdown -h now || systemctl --no-block poweroff || halt -p || true
  exit "${2:-1}"
}

# A hard stop powers the machine off from outside this script, so `die` never runs for it.
# This loop is what makes a run that never finishes diagnosable, and it is what makes a
# timed-out run worth paying for: every finished seed is already in the store.
(
  while true; do
    sleep @@HEARTBEAT_SECONDS@@
    { df -h /; free -g; uptime; ls -laR /opt/hapax/out/ 2>/dev/null | head -200; } \
      > /opt/hapax/log-machine.txt 2>&1 || true
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --include "run.toml" --include "*.json" --exclude "mamba/*" \
      --exclude "*/site-packages/*" --exclude "db/*" --region "$REGION" \
      >/dev/null 2>&1 || true
    aws s3 cp /opt/hapax/out/ "$OUTS3/" --recursive --region "$REGION" \
      >/dev/null 2>&1 || true
  done
) &

# ---------------------------------------------------------------- system packages
# Only what conda cannot provide: an archiver, and a compiler for this repo's own Rust.
# Everything the pipeline drives comes from environment.yml.
dnf install -y tar gzip bzip2 xz unzip zstd pigz gcc gcc-c++ make zlib-devel cargo \
  >/dev/null || die "dnf could not install the build tools"

# ---------------------------------------------------------------- the code
cd /opt/hapax
aws s3 cp "$S3/code/@@CODE_KEY_NAME@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
for f in pyproject.toml environment.yml requirements.txt hapax.py \
         experiments/libraries.json experiments/stage_references.json \
         experiments/nontemplate_prefix_map.py experiments/trim_from_prefix.py \
         experiments/subtract_known.py experiments/make_web.py; do
  test -f "/opt/hapax/$f" || die "the payload has no $f"
done

# ---------------------------------------------------------------- the documented install
# `conda env create -f environment.yml && pip install .` is PLAN.md 4's whole install.
# micromamba rather than a full conda, because it is one static binary and the solve is
# the same; the environment file is the committed one, unedited.
curl -fsSL --retry 5 "https://micro.mamba.pm/api/micromamba/linux-64/@@MICROMAMBA_VERSION@@" \
  | tar -xj -C /opt/hapax bin/micromamba \
  || die "micromamba @@MICROMAMBA_VERSION@@ did not download"
export MAMBA_ROOT_PREFIX=/opt/hapax/mamba
MM=/opt/hapax/bin/micromamba
"$MM" create -y -f /opt/hapax/environment.yml > /opt/hapax/log-conda.txt 2>&1 \
  || die "conda env create -f environment.yml failed; see log-conda.txt"
ENV=/opt/hapax/mamba/envs/hapax
export PATH="$ENV/bin:$PATH"
PY="$ENV/bin/python"
"$PY" -V || die "the environment has no python"
# Every binary the pipeline drives, proved to RUN and recorded before anything expensive
# starts. A number produced by an unrecorded version is a number without provenance, and
# an install that resolved to something other than TOOLS.md's pins should be visible in
# the log of the run that used it rather than inferred afterwards.
{ "$PY" -V; megahit --version; diamond --version; kraken2 --version; minimap2 --version;
  node --version; } > /opt/hapax/log-versions.txt 2>&1 \
  || die "the environment is missing a tool the pipeline drives; see log-versions.txt"
cat /opt/hapax/log-versions.txt
# `classify` is Kraken2's own helper and is on no distribution's PATH; hapax.tools finds
# it beside the wrapper. Proved here so a missing one costs ten minutes, not four hours.
"$PY" - <<'PYEOF' > /opt/hapax/log-classify.txt 2>&1 || die "hapax cannot find Kraken2's classify helper; see log-classify.txt"
from pathlib import Path
from hapax import tools
from hapax.config import Config
print(tools.resolve(Config(path=Path("preflight")), "classify"))
PYEOF
cat /opt/hapax/log-classify.txt

"$PY" -m pip install /opt/hapax > /opt/hapax/log-pip.txt 2>&1 \
  || die "pip install . failed; see log-pip.txt"
command -v hapax || die "pip install . left no hapax on PATH"
hapax --version || die "the installed hapax will not start"
"$PY" -m pip freeze > /opt/hapax/out/pip_freeze.txt 2>&1 || true

# ---------------------------------------------------------------- this repo's own Rust
# The linkage and census stages call experiments/kmer_coverage/target/release/*, and hapax
# preflights every selected stage, so a missing binary refuses the WHOLE run at the top
# rather than skipping those two. Built here rather than shipped, because the container's
# glibc is newer than this AMI's and a binary built there dies on the instance with a
# loader error that looks nothing like its cause.
( cd /opt/hapax/experiments/kmer_coverage && cargo build --release --offline 2>/dev/null \
    || cargo build --release ) || die "kmer_coverage did not build on the instance"
for b in kmer_coverage census; do
  test -x "/opt/hapax/experiments/kmer_coverage/target/release/$b" \
    || die "the crate built no $b; the linkage and census stages would refuse the run"
done

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
READS=/opt/hapax/data/reads
mkdir -p "$READS" /opt/hapax/sra
if [ ! -s "$READS/${ACC}_2.fastq.gz" ]; then
  prefetch --max-size 100G -O /opt/hapax/sra "$ACC" || die "prefetch $ACC failed"
  fasterq-dump --split-files --threads "$THREADS" --outdir "$READS" \
    --temp /opt/hapax/sra "/opt/hapax/sra/$ACC" || die "fasterq-dump $ACC failed"
  test -s "$READS/${ACC}_2.fastq" || die "$ACC did not come out as two mates; the archive \
record says paired and the stages below would read one mate twice"
  pigz -p "$THREADS" "$READS/${ACC}_1.fastq" "$READS/${ACC}_2.fastq" || die "pigz failed"
  rm -rf /opt/hapax/sra
fi
ls -l "$READS" > /opt/hapax/log-reads.txt
R1="$READS/${ACC}_1.fastq.gz"
R2="$READS/${ACC}_2.fastq.gz"

# ---------------------------------------------------------------- the documented fetch
# `hapax fetch-db --all` is the whole index set: SILVA converted to DNA, the four Kraken2
# databases in tier order, the two protein indexes, and Pfam-A. It records each download's
# URL, byte count and digest, and refuses a SILVA release whose size has moved, because a
# changed release is a changed operating point for the cascade's largest tier.
#
# The printed `[databases]` table is captured rather than retyped: a second copy of these
# paths here could disagree with what was actually fetched.
hapax fetch-db --dry-run --all > /opt/hapax/log-fetchdb-plan.txt 2>&1 || true
cat /opt/hapax/log-fetchdb-plan.txt
hapax fetch-db --all --dest /opt/hapax/db --threads "$THREADS" \
  > /opt/hapax/db_block.toml 2> /opt/hapax/log-fetchdb.txt \
  || die "hapax fetch-db --all failed; see log-fetchdb.txt"
cp -f /opt/hapax/db/fetch_manifest.json /opt/hapax/out/fetch_manifest.json || true
grep -q "^\[databases\]" /opt/hapax/db_block.toml \
  || die "hapax fetch-db printed no [databases] table, so there is nothing to run against"
du -sh /opt/hapax/db/* > /opt/hapax/log-db.txt 2>&1 || true

# ---- the one number nothing can inherit ----------------------------------------------
# SPEC 4.1: the mate-2 non-template prefix is per library and the system refuses to guess.
# Five seeds, because a mode that is not the same in every seed is not a length to cut, and
# experiments/trim_from_prefix.py refuses in that case rather than choosing.
mkdir -p /opt/hapax/out/nontemplate_prefix
for sd in @@PREFIX_SEEDS@@; do
  [ -s "/opt/hapax/out/nontemplate_prefix/nontemplate_prefix_s$sd.json" ] && continue
  "$PY" -I /opt/hapax/experiments/nontemplate_prefix_map.py --seed "$sd" \
    --pairs @@PREFIX_PAIRS@@ --r1 "$R1" --r2 "$R2" \
    --silva /opt/hapax/db/silva/silva_dna.fa --threads "$THREADS" \
    --work /opt/hapax/wk/ntprefix --out /opt/hapax/out/nontemplate_prefix \
    || die "the prefix measurement failed on seed $sd"
done
set +e
MEASURED=$("$PY" -I /opt/hapax/experiments/trim_from_prefix.py --library "$ACC" \
  --prefix-json /opt/hapax/out/nontemplate_prefix/nontemplate_prefix_s*.json \
  --out /opt/hapax/out/nontemplate_prefix/trim_from_prefix.json)
MEASURED_RC=$?
set -e
TRIM="@@TRIM_R2_STATED@@"
TRIM_BASIS=""
if [ -n "$TRIM" ]; then
  # A stated trim always wins and is never checked, because stating it is the override.
  # The measurement still ships beside it; what changes is whether its answer is USED.
  TRIM_BASIS="stated on the launcher; the instance's own measurement exited $MEASURED_RC \
with answer '${MEASURED:-none}' and was not used"
else
  [ "$MEASURED_RC" = "0" ] || die "no single mate-2 trim can be named for $ACC; see \
out/nontemplate_prefix/trim_from_prefix.json. The run stops here rather than inheriting \
another library's prefix or picking a mode by eye. State it on the launcher's --trim-r2 \
with a basis if this library cannot be measured this way"
  TRIM="$MEASURED"
  TRIM_BASIS="measured on this instance by experiments/nontemplate_prefix_map.py over \
@@PREFIX_SEEDS@@ seeds of @@PREFIX_PAIRS@@ pairs against SILVA, mode agreed by every seed \
(experiments/trim_from_prefix.py refuses otherwise)"
fi
echo "=== mate-2 trim for $ACC: $TRIM ($TRIM_BASIS) ==="

# ---------------------------------------------------------------- the run configuration
# One file, assembled from what was measured here and what fetch-db actually downloaded.
# The block stride is deliberately NOT stated: `hapax run` measures it off the read count,
# and a stride sized to another library leaves the later seeds reading nothing while the
# run directory still looks like a five-seed result.
cat > /opt/hapax/run.toml <<TOMLEOF
[inputs]
r1 = "$R1"
r2 = "$R2"
accession = "$ACC"

[inputs.sample_facts]
@@SAMPLE_FACTS@@

[library]
trim_r2 = $TRIM
trim_r2_basis = "$TRIM_BASIS"

[run]
seeds = @@SEEDS_COUNT@@
pairs = @@PAIRS@@
report_pairs = @@REPORT_PAIRS@@
threads = $THREADS
stages = [@@STAGES_TOML@@]
report_without = [@@REPORT_WITHOUT_TOML@@]

TOMLEOF
cat /opt/hapax/db_block.toml >> /opt/hapax/run.toml
cp -f /opt/hapax/run.toml /opt/hapax/out/run.toml
cat /opt/hapax/run.toml

# ---------------------------------------------------------------- the documented command
# `hapax run --config run.toml --out DIR --report`, once per seed, into one directory tree.
# --work is shared so the assembly is paid for once: it is the expensive half and nothing
# about it depends on the seed. --report means the pass only reports success when a page
# was built in which every declared section came from a file this run wrote.
for sd in @@SEEDS@@; do
  OUT="/opt/hapax/out/run_s$sd"
  [ -s "$OUT/report.html" ] && continue
  hapax run --config /opt/hapax/run.toml --out "$OUT" --seed "$sd" \
    --work /opt/hapax/wk --report \
    > "/opt/hapax/log-run-s$sd.txt" 2>&1 \
    || die "hapax run failed on seed $sd; its manifest and log-run-s$sd.txt say where"
  test -s "$OUT/report.html" || die "seed $sd exited 0 and wrote no report, which is the \
one outcome that must not be possible"
  aws s3 cp /opt/hapax/out/ "$OUTS3/" --recursive --region "$REGION" || true
done

# The assembly itself, so a later pass over this sample never re-assembles.
cp -f /opt/hapax/wk/asm_baseline/sample_contigs.fa /opt/hapax/out/ 2>/dev/null || true
gzip -f /opt/hapax/out/sample_contigs.fa 2>/dev/null || true

die "finished" 0
