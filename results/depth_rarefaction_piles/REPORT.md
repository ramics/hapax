# The depth ladder on the residue itself: 41%, and still not assembly

Tables: [`TABLES.md`](TABLES.md), generated from the committed per-level TSVs by
`experiments/report_depth_piles.py`. Five seeds, 4,000 query pairs per pile per seed.

## Why this run exists

[`../depth_rarefaction/`](../depth_rarefaction/REPORT.md) measured what a doubling of
sequencing depth would buy — 44.0% of fully unwitnessed reads gaining a witness — and
flagged its own limit in the same breath: the query block was a random sample of the
library, not the residue, so the ladder bounded the shape rather than giving the value for
the pile the question was about. This closes that gap. Same instrument, same levels, same
seeds; the query block is drawn per pair-level **pile** instead.

## The piles, validated before anything was read off them

The piles come from the residue thread's shared per-mate label table
(`library_label_classes.npz`: one class per mate against SILVA 138.2 and this sample's own
82,947 contigs at the cascade's 40-base operating point). A pair needs one class, and the
rule here is by priority — ribosomal, then identified, then a contig outside the
three-way split, then an unidentifiable contig, then nothing. That rule is this script's
choice, not the table's, and it is stated in `experiments/depth_rarefaction.py`.

The five piles partition the library exactly, and `unplaced` — both mates matching nothing
— comes out at **1,372,021 pairs, 15.71% of the library**, which is
[`../residue_census/`](../residue_census/REPORT.md) §2's figure to the pair. The script
refuses to run if it does not.

Query blocks are drawn at random across the whole library rather than from a contiguous
block, because the label table shows the file order is spatially structured and a
contiguous block is therefore not a sample of the library.

## The answer

**The unplaced residue's own rescue rate over the last measured doubling is 41.0%
(39.5–42.8), against 44.0% (40.6–46.2) for a random block of the library.** The library
figure was a good proxy after all — three points optimistic, in the direction predicted
when it was published.

The ladder underneath it is the sharper picture. The unplaced pile is 75.6% fully
unwitnessed at an eighth of the library and still 26.7% at the full library, against 0.1%
for ribosomal reads and 0.3% for identified ones; its median depth is 0 at every level
below the full library and 1 there. Nothing else in the library is remotely in this
condition, which is the census's finding reproduced on an independent draw.

**And the qualifier survives the change of pile: a witness is not a depth.** Of the
unplaced reads rescued over that doubling, **0.3% (0.2–0.5) reach median depth 5 or more**,
the depth where assembly and classification begin to work at all
(`../spike_abundance/`). Every other pile reads the same way — the highest is 2.2%, on the
repetitive contig-outside-the-split pile. So a doubling of this library moves the unplaced
residue from a median depth of 1 toward 2, and not toward anything assemblable.

**Two piles, two different propositions.** The residue on unidentifiable contigs — tier 2,
at median depth 6 — rescues at 65.6% (59.7–68.2) and is only 2.5% unwitnessed at full
depth. It is already past the wall the unplaced pile sits behind, which is why the two
halves of the residue were split between threads in the first place.

## What it does not say

The rates for the ribosomal and identified piles rest on 6 and 41 unwitnessed reads per
seed respectively, because reads from abundant material are almost never unwitnessed. They
are in the table with their n and should not be read as rates.

The stated bias of the label table carries over: it holds no read-level Kraken2 or
DIAMOND, so `unplaced` is more generous than the cascade's residue and errs toward calling
material unplaced. Direction: this pile is slightly closer to the library average than the
cascade's true residue, so **41% is, if anything, a small overestimate** of what a doubling
would rescue from the cascade's own residue.

The two biases named in the library-wide report also still hold. Subsampling models
re-sequencing *this* library deeper rather than sequencing new water, and drift makes new
water a different sample. And PCR duplicates read as depth, so a rescue is honest as "seen
twice" and inflated as "assembly gain" — which matters less here than it would have, since
the assembly gain is 0.3% either way.

## Provenance

| file | script | what |
|---|---|---|
| `depth_s{1..5}_f{0.125,0.25,0.5,1.0}.tsv.gz` | `experiments/depth_rarefaction.py --classes` | per-read depths per pile at each level |
| `manifest.json` | same | seeds, git SHA, script content hash, levels, label table, timings |
| `TABLES.md` | `experiments/report_depth_piles.py` | generated from the TSVs above |

Input not in the repo: `library_label_classes.npz` from the residue thread, via
`/mnt/project-files/shared/resasm/`, whose own manifest records the git SHA and script
hash that produced it.
