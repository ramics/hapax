# A doubling of depth buys recurrence, not assembly

Tables: [`TABLES.md`](TABLES.md), generated from the committed per-level TSVs by
`experiments/report_depth_rarefaction.py`. Five seeds, each with its own query block and
its own nested subsample.

## The question, and why it was worth asking locally

`results/residue_census/` §3 put the unplaced residue at median 31-mer depth 1 and
concluded that this is sampling rather than tuning: *you cannot assemble sequence you
hold one copy of*. That makes depth the governing variable of the whole residue problem,
and depth is the one variable money can change. "Should we sequence this water deeper"
is therefore the most expensive open question in the project — and it does not have to be
answered by spending. A fixed set of real molecules can be scored against nested
subsamples of the library we already have, and the slope read off directly.

## What was measured

A fixed query block of 10,000 real pairs per seed, scored against nested libraries at
1/8, 1/4, 1/2 and 1 of the run with exact whole-library 31-mer counts
(`experiments/kmer_coverage`). Nested, so the ladder is a rarefaction of one library
rather than four unrelated samples. The query pairs are excluded from every library at
every level, so a depth is always "occurrences elsewhere" and no level is contaminated by
a read witnessing itself.

"Fully unwitnessed" is the strict definition — every scorable 31-mer of the read having
no witness anywhere else. The tool's zero-run field exists precisely to separate a rare
molecule, which is unwitnessed end to end, from a substitution error, which destroys the
up to 31 k-mers spanning it inside otherwise covered sequence. The strict definition
excludes the error case by construction, which is why it is the one reported.

## The two numbers, and they point different ways

**Seen once becomes seen twice, fast.** Of mate 1's fully unwitnessed reads at half the
library, 44.0% (40.6–46.2 over five seeds) have a witness at the full library. The rate
rises along the ladder — 24.4%, then 34.0%, then 44.0% — because each doubling adds twice
the absolute number of reads of the one before it while witnessing depends on absolute
depth. So the last measured transition is a floor for the next doubling rather than a
ceiling. This is measurement and not a fit: the same named reads are scored against both
libraries and the join is on the name.

**Assembly-grade depth barely moves.** Of those same rescued reads, 0.3% (0.1–0.5) reach
depth 5 or more. Across the whole ladder the fraction of the library below depth 5 falls
from 38.5% to 24.4% over an eightfold increase in effort — about 4.7 points per doubling,
against 5.6 points per doubling for the fully-unwitnessed fraction. A doubling therefore
turns a molecule seen once into a molecule seen twice, and almost never into a molecule
deep enough to assemble.

## What that means for the depth-1 pile

The instinct behind "sequence it deeper" is assembly: get the pile to a depth where
contigs form, then characterise the contigs. **That instinct is wrong at the price of one
doubling.** Nothing here reaches the ~2x breadth coverage where assembly and
classification start working at all (`results/spike_abundance/`), and the ~10x where the
spike assembles properly is three and a half doublings away.

What a doubling does buy is the thing the project has been unable to get any other way:
**a second observation.** Every grouping and recurrence method measured here needs one —
assembly, neighbourhood kinship, mate-linkage kinship, within-run recurrence — and on
this pile they all have nothing to work with because the molecule is there once. Turning
44% of "once" into "twice" does not create assemblies, but it does create the minimum
input for *recurrence*, which is the axis the companion thread is measuring across the
other 115 samples, and it moves reads out of the one category nothing can distinguish
from error.

So the honest recommendation is narrow and it is not "sequence deeper to assemble the
residue": a doubling of this library would roughly halve the fraction of it that is
unwitnessed, and that is worth having for confirmation and recurrence rather than for
contigs. Whether it is worth the money is Imo's call, and this is the curve to make it on.

## Stated biases, with direction

The full list with directions is in `TABLES.md`. The three that matter:

1. **This models re-sequencing the same library deeper, not sequencing new water.** The
   same plant's rRNA fraction went 38.0% to 67.1% in eight months
   (`results/dna_rna_pair/`), so new water samples a different population. Direction:
   optimistic about new water, roughly right about deeper sequencing of this library.
2. **Duplicates read as depth.** A read rescued by a PCR duplicate of its own molecule has
   gained confirmation that it is not an error, but no new sequence. Direction: inflates
   the rescue rate as a measure of assembly gain, honest as a measure of "seen twice".
3. **A random query block is not the residue.** These are real library reads of every
   kind, so the ladder describes the library. The unplaced residue sits at the bottom of
   this distribution by definition, and its own rescue rate will be lower than the library
   average. The ladder bounds the shape; it does not give that pile's value. Measuring the
   residue's own transition needs the cascade's residue read list, which lives with the
   pipeline thread.

## Provenance

| file | script | what |
|---|---|---|
| `depth_s{1..5}_f{0.125,0.25,0.5,1.0}.tsv` | `experiments/depth_rarefaction.py` | per-read depths at each level |
| `curve.json` | same | the curve and its linear-fit residual |
| `manifest.json` | same | seeds, git SHA, script content hash, levels, timings |
| `TABLES.md` | `experiments/report_depth_rarefaction.py` | generated from the TSVs above |

The analysis lives in a separate script from the generator on purpose, so it can be
extended without invalidating the manifest's hash of the code that produced the data.

The run shared the machine with the `profile_reach` searches, so the timings in
`manifest.json` are operational rather than benchmarks.
