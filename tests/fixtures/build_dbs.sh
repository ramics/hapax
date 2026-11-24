#!/usr/bin/env bash
# Build the fixture's toy indexes from its committed source sequence.
#
# WHY THE INDEXES ARE BUILT AND NOT COMMITTED. A .dmnd file, a Kraken2 hash and a
# minimap2 .mmi are build artefacts of one tool version; the versions are pinned in
# TOOLS.md and a committed binary would silently stop matching them. The source FASTAs
# are committed, this script turns them into indexes, and it is the same script CI runs
# and a developer runs, so there is one recipe rather than one per machine. In step 7 of
# PLAN.md this becomes `hapax fetch-db` for the real databases.
#
# WHY THE TAXONOMY IS WRITTEN BY HAND. kraken2-build --build needs nodes.dmp and
# names.dmp, and --download-taxonomy fetches about 400 MB from NCBI, which is the wrong
# dependency for a test that has to run on every pull request. Four nodes are enough for
# a two-sequence library, and writing them here keeps the fixture offline.
#
# Usage: tests/fixtures/build_dbs.sh [OUTDIR]     (default: a dbs/ beside this script)
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${1:-$here/dbs}"
mkdir -p "$out"

for tool in diamond kraken2-build minimap2; do
  command -v "$tool" >/dev/null || { echo "$tool is not on PATH; see TOOLS.md" >&2; exit 1; }
done

echo "== protein index (DIAMOND) =="
diamond makedb --in "$here/db/prot_known.faa" --db "$out/prot_known" --quiet

# A SECOND, BROADER PROTEIN TIER. The distance stage takes an ordered pair of protein
# indexes, viral then Swiss-Prot, and refuses any other number: it reports which of the
# two held the nearest relative, and one index cannot answer that. The second tier here is
# the open reading frames of the fixture's bacterial background, found on the committed
# slice at build time rather than fetched, which keeps the fixture offline.
#
# WHAT IT MUST NOT CONTAIN is the divergent arm's own source. The arm is a mutated PMMoV
# and PMMoV's real proteins are tier one; if they were in both tiers the stage would read
# the same identity whichever tier won, which is the trap of a query sitting inside the
# reference it is searched against - it reads 100% on both piles and the control moves
# with the test.
python3 - "$here/genomes/background_ecoli.fa" "$out/prot_broad.faa" <<'PY'
import sys

CODONS = {}
BASES = "TCAG"
AMINO = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for i, b1 in enumerate(BASES):
    for j, b2 in enumerate(BASES):
        for k, b3 in enumerate(BASES):
            CODONS[b1 + b2 + b3] = AMINO[i * 16 + j * 4 + k]

src, dst = sys.argv[1], sys.argv[2]
seq = "".join(l.strip() for l in open(src) if not l.startswith(">")).upper()
rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))

def translate(s):
    return "".join(CODONS.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

# Six frames, each cut at stops, keeping anything 60 aa or longer. 60 because a DIAMOND
# subject shorter than that gives alignments too short to clear an E-value floor on a
# 150 nt read, and a subject nothing can significantly hit is a subject that is not in
# the database as far as any measurement is concerned.
peptides = []
for strand, s in (("plus", seq), ("minus", rc)):
    for frame in range(3):
        for piece in translate(s[frame:]).split("*"):
            if len(piece) >= 60:
                peptides.append((strand, frame, piece))
if not peptides:
    raise SystemExit("no open reading frame of 60 aa in the background slice")
with open(dst, "w") as fh:
    for n, (strand, frame, piece) in enumerate(peptides, 1):
        fh.write(f">fixture_background_orf_{n} {strand} frame {frame + 1}\n{piece}\n")
print(f"{len(peptides)} background ORFs")
PY
diamond makedb --in "$out/prot_broad.faa" --db "$out/prot_broad" --quiet

echo "== nucleotide index (minimap2) =="
minimap2 -d "$out/nt_known.mmi" "$here/db/nt_known.fa" 2>/dev/null
minimap2 -d "$out/silva_toy.mmi" "$here/db/silva_toy.fa" 2>/dev/null

echo "== classifier library (Kraken2) =="
k2="$out/k2_toy"
rm -rf "$k2"
mkdir -p "$k2/taxonomy"

# A minimal NCBI-format taxonomy: root, the two domains the fixture's named material
# falls in, and one species under each. The taxids are the real NCBI ones (562
# Escherichia coli, 329852 Emesvirus zinderi / phage MS2) so that a name the classifier
# returns is the name the organism actually has, which is what the report prints.
cat > "$k2/taxonomy/nodes.dmp" <<'DMP'
1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	1	|	0	|	0	|	0	|
2	|	1	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|
10239	|	1	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|
562	|	2	|	species	|	EC	|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|
329852	|	10239	|	species	|	MS	|	0	|	1	|	11	|	1	|	0	|	1	|	0	|	0	|
DMP
cat > "$k2/taxonomy/names.dmp" <<'DMP'
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
10239	|	Viruses	|		|	scientific name	|
562	|	Escherichia coli	|		|	scientific name	|
329852	|	Emesvirus zinderi	|		|	scientific name	|
DMP

# Kraken2 takes the taxid from the sequence header, so the library is the same sequence
# as the nucleotide index with a taxid written onto each record.
#
# AND IT HAS TO BE ADDED WITH --add-to-library, NOT COPIED INTO library/. Dropping the
# FASTA into that directory looks like it should work, and the build then aborts with "no
# preliminary seqid/taxid mapping files found": the seqid-to-taxid mapping is derived from
# the headers by --add-to-library, which writes prelim_map.txt, and --build reads that
# rather than the sequence. --no-masking because low-complexity masking calls dustmasker,
# and a toy index should not depend on a tool nothing else in the fixture needs.
python3 - "$here/db/nt_known.fa" "$out/k2_library.fna" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
taxid = {"NC_000913.3": 562, "NC_001417.2": 329852}
with open(src) as fh, open(dst, "w") as out:
    for line in fh:
        if line.startswith(">"):
            acc = line[1:].split()[0].split(":")[0]
            if acc not in taxid:
                raise SystemExit(f"no taxid for {acc}; add it to this script")
            out.write(f">{line[1:].split()[0]}|kraken:taxid|{taxid[acc]}\n")
        else:
            out.write(line)
PY

# k=31 because 31 is the k this repo's every k-mer operating point is stated at, and the
# minimizer length follows from it. --minimizer-spaces must be at most a quarter of the
# minimizer length, which kraken2-build enforces: its default of 7 is legal under the
# default l=31 and not under l=25, so leaving it out fails the build rather than being
# clamped. Stated here rather than inherited.
kraken2-build --add-to-library "$out/k2_library.fna" --db "$k2" --no-masking
kraken2-build --build --db "$k2" --threads "${THREADS:-2}" \
  --kmer-len 31 --minimizer-len 25 --minimizer-spaces 6

echo "built into $out:"
ls -la "$out"
