#!/bin/bash
# The whole reach run: no SSH, no SSM, no inbound access, so this script is the only thing
# that will ever happen on the instance and a second attempt costs another launch. It lives
# in the repo rather than in the launcher's Python because it is the artefact worth reviewing
# before launch; `--print-user-data` renders it. Tokens are upper-case names between double
# at-signs, so every brace below is bash's own.
#
# ORDER IS THE SAFETY PROPERTY, from Imo's launch guide and unchanged from
# `bigindex_userdata.sh`: the dead man's switch first, before anything that can fail; a
# heartbeat to S3, because console output is unreliable; one exit helper on every path,
# cancelling the pending timer first because `shutdown -h now` is a no-op while one is
# scheduled; python3.12 by name, since this AMI's python3 is 3.9; and `test -f` on every
# input, because `aws s3 cp --recursive` over an empty prefix exits 0 and the arms would
# start on nothing.
#
# WHAT IS DIFFERENT FROM THE BIG-INDEX RUN. That run was one program against one enormous
# database, so its risk was time. This one is seven programs against six small databases, so
# its risk is setup: any one tool or download failing takes an arm with it. So the staging is
# graded. A failure in the shared protein set is fatal, because every arm is scored against
# it and a table missing its baseline is not a table. A failure in Pfam, in the structure
# databases or in the viral genomes drops exactly that arm, records why, and lets the rest
# run - a six-arm table with one honest gap beats a relaunch.
set -xeuo pipefail
mkdir -p /opt/hapax/out /opt/hapax/bin /opt/hapax/db
exec > >(tee /dev/console /opt/hapax/log-userdata.txt) 2>&1
echo "=== hapax reach run starting $(date -Is) commit @@COMMIT@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

S3="s3://@@BUCKET@@/@@PREFIX@@"
REGION="@@REGION@@"
IID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo unknown)
LOGS="$S3/logs/reach-$IID"
OUT="$S3/results/reach/$IID"
THREADS=@@THREADS@@

die() {
  echo "=== shutting down: $1 ==="
  echo "$1" > /opt/hapax/why.txt 2>/dev/null || true
  sync || true
  # Logs and the exit reason go up on every path, including the failures. An object in the
  # store is the only record; GetConsoleOutput has come back empty for the first ten minutes
  # after exactly the instance somebody needed to read.
  aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
    --include "why.txt" --region "$REGION" || true
  aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --exclude "*" --include "rows_*.jsonl" \
    --include "timing_*.json" --include "*.json" --include "orfs.faa" \
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
    df -h / > /opt/hapax/log-disk.txt 2>&1 || true
    aws s3 cp /opt/hapax/ "$LOGS/" --recursive --exclude "*" --include "log-*.txt" \
      --region "$REGION" >/dev/null 2>&1 || true
    aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --exclude "*" --include "rows_*.jsonl" \
      --include "timing_*.json" --include "*.json" --region "$REGION" >/dev/null 2>&1 || true
  done
) &

dnf install -y tar gzip xz gcc make findutils >/dev/null || die "dnf could not install the build tools"

PYBIN=""
for cand in python3.12 python3.11; do
  dnf install -y "$cand" "$cand-pip" >/dev/null 2>&1 || true
  if command -v "$cand" >/dev/null 2>&1; then PYBIN="$cand"; break; fi
done
[ -n "$PYBIN" ] || die "no python 3.11+ available; nothing to run"
echo "=== interpreter: $PYBIN ($($PYBIN -V)) ==="

# pyrodigal calls genes and hands back the CDS nucleotide beside the protein, which arm 5
# needs and a GFF plus a contig would make an off-by-one away from wrong. parasail is the
# vectorised Smith-Waterman arm 2 is. Both are wheels, so no compiler is involved and a
# build failure cannot take the run with it.
"$PYBIN" -m pip install --quiet --upgrade pip >/dev/null 2>&1 || true
"$PYBIN" -m pip install --quiet pyrodigal==3.7.1 parasail==1.3.4 \
  || die "pyrodigal/parasail did not install; arms 2 and the ORF calling cannot run"

# ---------------------------------------------------------------- the searchers
cd /opt/hapax
# BLAST+ because the staged NCBI databases are BLAST volumes and nothing else reads them;
# DIAMOND at 2.1.10 because that is the version every timing in `results/diamond_scaling/`
# was measured with, and a comparison against the production baseline has to be against the
# production binary.
curl -fsSL --retry 5 -o blast.tar.gz \
  "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/@@BLAST_VERSION@@/ncbi-blast-@@BLAST_VERSION@@+-x64-linux.tar.gz" \
  || die "BLAST+ @@BLAST_VERSION@@ did not download"
tar xzf blast.tar.gz && mv "ncbi-blast-@@BLAST_VERSION@@+" blast
curl -fsSL --retry 5 -o diamond.tar.gz \
  "https://github.com/bbuchfink/diamond/releases/download/v@@DIAMOND_VERSION@@/diamond-linux64.tar.gz" \
  || die "DIAMOND @@DIAMOND_VERSION@@ did not download"
tar xzf diamond.tar.gz
cp /opt/hapax/diamond /opt/hapax/bin/diamond
for prog in blastp tblastn blastdbcmd makeblastdb; do
  cp "/opt/hapax/blast/bin/$prog" "/opt/hapax/bin/$prog" || die "BLAST+ is missing $prog"
done
/opt/hapax/bin/blastp -version || die "BLAST+ will not run"
/opt/hapax/bin/diamond version || die "DIAMOND will not run"

# HMMER from source: there is no dnf package on this AMI and the arm that needs it is the
# one the literature says should win, so it is fatal rather than graded.
curl -fsSL --retry 5 -o hmmer.tar.gz \
  "http://eddylab.org/software/hmmer/hmmer-@@HMMER_VERSION@@.tar.gz" \
  || die "HMMER @@HMMER_VERSION@@ did not download"
tar xzf hmmer.tar.gz && cd "hmmer-@@HMMER_VERSION@@"
./configure --prefix=/opt/hapax/hmmer >/opt/hapax/log-hmmer-configure.txt 2>&1 \
  || die "HMMER configure failed; see log-hmmer-configure.txt"
make -j"$THREADS" >/opt/hapax/log-hmmer-make.txt 2>&1 || die "HMMER build failed"
make install >>/opt/hapax/log-hmmer-make.txt 2>&1 || die "HMMER install failed"
cd /opt/hapax
for prog in hmmsearch jackhmmer hmmpress; do
  cp "/opt/hapax/hmmer/bin/$prog" "/opt/hapax/bin/$prog" || die "HMMER is missing $prog"
done
/opt/hapax/bin/hmmsearch -h >/dev/null || die "hmmsearch will not run"

# NO FOLDSEEK BINARY IS STAGED. The structure arm moved to the container (see the note
# further down where its databases used to be built), so downloading it here would be an
# unused dependency whose failure could still drop a run. `FOLDSEEK_SOURCE` is kept in the
# manifest, reading `not staged`, rather than removed: a field that disappears looks like a
# field nobody recorded.
FOLDSEEK_SOURCE="not staged: the structure arm runs in-container"
echo "=== foldseek: $FOLDSEEK_SOURCE ===" | tee /opt/hapax/log-foldseek.txt

# ---------------------------------------------------------------- code and panel
aws s3 cp "$S3/code/@@CODE_KEY@@" code.tar.gz --region "$REGION" \
  || die "the code payload did not arrive"
tar xzf code.tar.gz
test -f experiments/reach_arms.py || die "the payload is short: no reach_arms.py"

mkdir -p /opt/hapax/queries
aws s3 cp "$S3/scratch/reach-panel-@@PANEL_SET@@/" /opt/hapax/queries/ \
  --recursive --region "$REGION" || die "the panel did not copy"
test -f /opt/hapax/queries/panel.fa || die "no panel.fa; the arms would start on nothing"
# One null pile per seed. Five, because the real arms are fixed sequence and the shuffle is
# the only place a seed means anything in this comparison - so a null that reads zero reads
# zero five times or it is one draw.
NULLS=$(ls /opt/hapax/queries/null_s*.fa 2>/dev/null | tr '\n' ' ')
[ -n "$NULLS" ] || die "no null_s*.fa; there would be no shuffled control"
grep -c '^>' /opt/hapax/queries/panel.fa $NULLS

# Imports before the long read. On the instance a short payload costs hours and a launch;
# locally it costs a second, and this is the check that catches the module that was added to
# a tool's imports and not to the launcher's payload list.
"$PYBIN" - <<'PRECHECK' || die "the shipped payload cannot import"
import importlib.util
spec = importlib.util.spec_from_file_location("arms", "experiments/reach_arms.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
import pyrodigal, parasail                                           # noqa: F401
print("=== payload imports cleanly, arms:", m.ARMS, "===")
PRECHECK

# ---------------------------------------------------------------- the shared subject set
# FATAL, because every arm is scored against it. Swiss-Prot as BLAST volumes from the public
# NCBI bucket, unsigned: the bucket is public and this account has no grant on it, so a
# SIGNED request is the one that fails, with a bare AccessDenied that reads like a bucket
# policy. RefSeq viral proteins on top, because that is the pair of sets
# `results/divergent_known/` searched and this run has to be comparable to it.
aws s3 cp "s3://@@NCBI_BUCKET@@/@@NCBI_PREFIX@@/" /opt/hapax/db/ --recursive \
  --no-sign-request --exclude "*" --include "swissprot.*" --region "$REGION" \
  || die "swissprot did not stage"
/opt/hapax/bin/blastdbcmd -db /opt/hapax/db/swissprot -entry all -out /opt/hapax/db/prot.faa \
  || die "swissprot staged but will not dump; the shared subject set cannot be built"

VIRAL_PROT=0
for n in 1 2 3 4 5 6 7 8; do
  if curl -fsSL --retry 3 -o "/opt/hapax/db/viral.$n.protein.faa.gz" \
      "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.$n.protein.faa.gz"; then
    gzip -dc "/opt/hapax/db/viral.$n.protein.faa.gz" >> /opt/hapax/db/prot.faa
    VIRAL_PROT=$((VIRAL_PROT + 1))
  fi
done
[ "$VIRAL_PROT" -gt 0 ] || die "no RefSeq viral protein release files downloaded"
grep -c '^>' /opt/hapax/db/prot.faa > /opt/hapax/log-subjects.txt
/opt/hapax/bin/makeblastdb -in /opt/hapax/db/prot.faa -dbtype prot \
  -out /opt/hapax/db/prot -title "swissprot+refseq_viral" \
  || die "the shared subject set will not format for BLAST"
/opt/hapax/bin/diamond makedb --in /opt/hapax/db/prot.faa --db /opt/hapax/db/prot.dmnd \
  --threads "$THREADS" --quiet || die "the shared subject set will not format for DIAMOND"

# ---------------------------------------------------------------- graded databases
ARMS="diamond blastp exhaustive"

# Pfam-A, for the profile-library arm. Graded: if the release moves, the arm drops.
if curl -fsSL --retry 3 -o /opt/hapax/db/Pfam-A.hmm.gz \
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"; then
  gzip -df /opt/hapax/db/Pfam-A.hmm.gz
  curl -fsSL --retry 3 -o /opt/hapax/db/Pfam-A.clans.tsv.gz \
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz" \
    && gzip -df /opt/hapax/db/Pfam-A.clans.tsv.gz || true
  curl -fsSL --retry 3 -o /opt/hapax/db/Pfam.version.gz \
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" \
    && gzip -df /opt/hapax/db/Pfam.version.gz || true
  if /opt/hapax/bin/hmmpress /opt/hapax/db/Pfam-A.hmm; then
    ARMS="$ARMS hmmsearch"
  fi
fi
# jackhmmer needs only the shared set, which is fatal above, so it is always in.
ARMS="$ARMS jackhmmer"

# hh-suite, for the profile-to-profile arm. Graded, and ordered after jackhmmer because the
# query profile is built from the alignment jackhmmer writes - deliberately, rather than by
# an hhblits pass against UniRef30, which would be a 50 GB download AND a different subject
# universe behind the query side from the one every other arm here searches.
if curl -fsSL --retry 3 -o hhsuite.tar.gz \
    "https://github.com/soedinglab/hh-suite/releases/download/v@@HHSUITE_VERSION@@/hhsuite-@@HHSUITE_VERSION@@-AVX2-Linux.tar.gz"; then
  mkdir -p /opt/hapax/hhsuite && tar xzf hhsuite.tar.gz -C /opt/hapax/hhsuite
  for prog in hhsearch hhmake; do
    find /opt/hapax/hhsuite -name "$prog" -type f -perm -u+x -exec cp {} /opt/hapax/bin/ \; || true
  done
  mkdir -p /opt/hapax/db/hh_pfam
  if [ -x /opt/hapax/bin/hhsearch ] && [ -x /opt/hapax/bin/hhmake ] \
      && curl -fsSL --retry 3 -o pfam_hh.tar.gz "@@HH_PFAM_URL@@"; then
    tar xzf pfam_hh.tar.gz -C /opt/hapax/db/hh_pfam --strip-components=0
    # The archive's internal prefix is not fixed between releases, so the database name is
    # discovered from what landed rather than assumed. An assumed name is a wasted launch.
    HHBASE=$(find /opt/hapax/db/hh_pfam -name "*_hhm.ffindex" | head -1 | sed 's/_hhm.ffindex//')
    if [ -n "$HHBASE" ]; then
      ln -sf "$HHBASE"_hhm.ffindex /opt/hapax/db/hh_pfam/pfam_hhm.ffindex || true
      for suffix in hhm.ffdata hhm.ffindex a3m.ffdata a3m.ffindex cs219.ffdata cs219.ffindex; do
        [ -e "${HHBASE}_${suffix}" ] && ln -sf "${HHBASE}_${suffix}" \
          "/opt/hapax/db/hh_pfam/pfam_${suffix}" || true
      done
      ARMS="$ARMS hhsearch"
    fi
  fi
fi

# RefSeq viral genomes, for the codon arm: tblastn against these gives the subject's
# nucleotide coordinates, which is what makes a codon alignment reconstructable without a
# network round trip per subject.
# TWO NAMING PATTERNS, AND THE FIRST RUN FOUND OUT WHY. RefSeq's viral release names its
# protein files `viral.N.protein.faa.gz` but its genomic files `viral.N.1.genomic.fna.gz` -
# the extra index is part of the genomic name and not of the protein one. Assuming the
# protein pattern held gave eight 404s on run i-0672a2ca4a782aa21 and dropped the codon arm
# from that run's table. Both patterns are tried here rather than the corrected one alone,
# because a release that renames again should degrade to the other spelling rather than to
# a silent zero.
VIRAL_NT=0
for n in 1 2 3 4 5 6 7 8; do
  for name in "viral.$n.1.genomic.fna.gz" "viral.$n.genomic.fna.gz"; do
    if curl -fsSL --retry 3 -o "/opt/hapax/db/$name" \
        "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/$name"; then
      gzip -dc "/opt/hapax/db/$name" >> /opt/hapax/db/viral_nt.fna
      VIRAL_NT=$((VIRAL_NT + 1))
      break
    fi
  done
done
[ "$VIRAL_NT" -gt 0 ] || echo "=== no viral genomic release files; the codon arm will be \
absent from this table rather than wrong ===" | tee -a /opt/hapax/log-arms.txt
if [ "$VIRAL_NT" -gt 0 ] && /opt/hapax/bin/makeblastdb -in /opt/hapax/db/viral_nt.fna \
    -dbtype nucl -out /opt/hapax/db/viral_nt -title "refseq_viral_genomic"; then
  ARMS="$ARMS codon"
fi

# THE STRUCTURE ARM DOES NOT RUN HERE ANY MORE, and that is a deliberate move rather than a
# step skipped. It is now `experiments/reach_structure.py`, which runs in the container for
# nothing against every solved virus structure in the PDB. Three things pushed it off the
# instance, in increasing order of importance.
#
# One: `foldseek databases` exits 0 on a failed download. Its downloads come from a single
# Cloudflare worker, that worker answered HTTP 429 with a 17-byte body from three networks
# across several hours on 2026-11-22, and foldseek untarred the error body, printed
# `gzip: stdin: not in gzip format`, moved nothing, and exited 0. Two runs added structure
# arms on that exit code and then failed on `Input .../afdb_sprot does not exist`.
#
# Two: the only bulk-servable alternative, EBI's 28.6 GB `swissprot_pdb_v6.tar`, came down
# at 6.7 MB/s and then stopped advancing entirely at 33 GB of a 200 GB disk - no error, no
# socket timeout on the curl - and took a whole run's dead-man timer with it.
#
# Three, and this is the one that matters: AlphaFold/Swiss-Prot was chosen because it is
# Swiss-Prot folded, which would have made the arm these same subjects in another alphabet.
# But the shared database here is Swiss-Prot PLUS the RefSeq viral proteins, and
# `report_reach.py` §5b counts 47 of the 51 ORFs finding their nearest relative in the
# RefSeq viral half, which is folded nowhere. So that target could not have answered the
# question it was chosen for, and paying to stage it again would buy a weak null.
#
# The ProstT5 weights are therefore not staged here either. Nothing in this script needs
# them, and `reach_structure.py` pulls them from the same S3 key into the container.
echo "=== structure arm runs in-container via experiments/reach_structure.py, not here \
===" | tee -a /opt/hapax/log-arms.txt

du -sh /opt/hapax/db > /opt/hapax/log-db.txt
df -h / >> /opt/hapax/log-db.txt
echo "=== arms that will run: $ARMS ===" | tee /opt/hapax/log-arms.txt

# ---------------------------------------------------------------- the run
# `reach_arms.py` wraps each arm itself, so one tool dying inside it does not stop the
# others and its traceback lands in that arm's timing file. This only fails the run if every
# arm failed, which means the setup is wrong rather than one tool.
# shellcheck disable=SC2086
"$PYBIN" experiments/reach_arms.py --seed @@SEED@@ \
  --panel /opt/hapax/queries/panel.fa --null $NULLS \
  --out /opt/hapax/out --threads "$THREADS" \
  --bin-dir /opt/hapax/bin --db-dir /opt/hapax/db \
  --arms $ARMS > /opt/hapax/log-arms-run.txt 2>&1 \
  || die "every arm failed; see log-arms-run.txt"

# Provenance the report cannot reconstruct from the rows: which binary, which release, which
# subject counts. Written here because this script is the only thing that knows.
cat > /opt/hapax/out/stage_manifest.json <<MANIFEST
{
 "instance_id": "$IID",
 "instance_type": "@@INSTANCE_TYPE@@",
 "commit": "@@COMMIT@@",
 "seed": @@SEED@@,
 "threads": $THREADS,
 "arms_staged": "$ARMS",
 "blast_version": "@@BLAST_VERSION@@",
 "diamond_version": "@@DIAMOND_VERSION@@",
 "hmmer_version": "@@HMMER_VERSION@@",
 "hhsuite_version": "@@HHSUITE_VERSION@@",
 "null_piles": "$NULLS",
 "foldseek": "$FOLDSEEK_SOURCE",
 "ncbi_blast_snapshot": "@@NCBI_PREFIX@@",
 "refseq_viral_protein_files": $VIRAL_PROT,
 "refseq_viral_genomic_files": $VIRAL_NT,
 "shared_subject_sequences": $(cat /opt/hapax/log-subjects.txt 2>/dev/null || echo 0),
 "pfam_version": "$(head -3 /opt/hapax/db/Pfam.version 2>/dev/null | tr '\n' ' ' || echo unknown)"
}
MANIFEST

aws s3 cp /opt/hapax/out/ "$OUT/" --recursive --region "$REGION" \
  || echo "=== the final upload failed; the heartbeat copies are what there is ==="
aws s3 cp /opt/hapax/queries/panel.json "$OUT/panel.json" --region "$REGION" || true
die "the run finished" 0
