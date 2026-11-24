# The end-to-end fixture

A 9,987-pair library carrying all three of the report's output classes, with toy indexes
built from its own sequence. It exists so that `hapax run` can be a test instead of
something somebody runs by hand on an instance and reads afterwards.

**It is not a benchmark, and no number from it belongs in `results/`.** Three genomes, ten
thousand pairs, and indexes built from the fixture's own sequence: recall against them
says nothing about recall against RefSeq. It answers one question, which is whether every
section of the report is filled from a file the run wrote.

## What is in it

| arm | sequence | depth | what it is for |
|---|---|---|---|
| known virus | phage MS2, `NC_001417.2`, 3,569 nt, unmodified | 300x | in the toy nucleotide index and the toy classifier library, so it must be subtracted as tier 1 |
| divergent virus | pepper mild mottle virus, `NC_003630.1`, mutated here | 150x | the **real** PMMoV is in the toy protein index and nowhere else, so this arm is divergent-but-related |
| bacterial background | `NC_000913.3:4164000-4180000`, the E. coli rrnB window | 60x | the classifier and SILVA arms; rRNA is about a third of the discriminating library and an rRNA stage with no rRNA has nothing to do |
| unplaceable | 4,000 nt, composition-matched to the fixture, in no index | 3x | the residue report and the depth-1 pile |

MS2's four ORFs are maturation, coat, lysis and replicase, which are the roles the
gene-level asymmetry work is stated in (SPEC 2.2), so the fixture exercises per-gene
ranking and not only a whole-contig identity.

## The divergent arm's distance is measured, not stated

`make_fixture.py` substitutes nucleotides at 15%, rejecting any substitution that would
make a stop codon inside an annotated CDS or touch a start codon, and then **measures**
the protein identity of each ORF against the real PMMoV and writes it into
`manifest.json`. The seed-1 build gives 68.4%, 68.9%, 68.9% and 71.3% per ORF, which sits
beside the holdout positive control's 66.5% — the one distance in this repo nobody chose.

15% was reached by measuring: 22% put the arm at 56%, below the holdout. The mapping from
a nucleotide rate to a protein identity runs through the genetic code and the rejection
rule, so it is not worth predicting.

`tests/test_end_to_end.py` asserts the **band** rather than the number. An arm at 99%
would be subtracted as tier 1 and an arm at 30% would be unplaceable, and either way the
section under test would be measuring something else.

## Error model

0.00068 per base along the body and 0.0019 at the final base: measured reference-free from
self-overlapping pairs of SRR23973193 (`results/berlin_quality/`, CLAUDE.md). A fixture
with no errors lets an assembler look better than it is, and an invented rate is a number
without provenance.

## Rebuilding

```
python tests/fixtures/make_fixture.py --seed 1                 # fetches from NCBI
python tests/fixtures/make_fixture.py --seed 1 --refs DIR      # from local copies
tests/fixtures/build_dbs.sh OUTDIR                             # the toy indexes
```

The committed files are the output of `--seed 1`, and a rebuild is byte-identical: the
gzip headers are written with `mtime=0` so a checksum does not change with the clock.
`manifest.json` carries the seed, the accessions, the download date, this script's content
hash and the measured divergences, and `test_the_fixture_is_present_and_pinned` fails if
the tree and the manifest disagree.

The indexes are **built, not committed**. A `.dmnd`, a Kraken2 hash and a `.mmi` are build
artefacts of one tool version, the versions are pinned in `TOOLS.md` and requested by
`environment.yml`, and a committed binary would quietly stop matching them. The Kraken2
taxonomy is four nodes written by `build_dbs.sh`, because `--download-taxonomy` fetches
about 400 MB from NCBI and that is the wrong dependency for a test that runs on every pull
request.

There are **two protein indexes**, in order: `prot_known.dmnd` from the committed viral
proteins, and `prot_broad.dmnd` from the open reading frames of the bacterial background
slice, found at build time so that nothing is fetched. The distance stage takes the pair
and reports which of the two held the nearest relative, and refuses a list of any other
length. The second tier deliberately excludes the divergent arm's source: the real PMMoV
proteins are in tier one only, because a query that sits inside the reference it is
searched against reads the same identity on both piles.

## What ground truth is allowed to do here

Every read's name carries the arm it came from. That is usable only by a test: it tells
the end-to-end test which output class a missing section lost, rather than only that the
section is missing. It never reaches a report. A known distance validates this tool and is
not output of it (SPEC 13.4), and the report has no ground-truth section on any sample.
