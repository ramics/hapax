# The weak-hit headroom is the chance band — but merging converts, and noise never does

Tables: [`TABLES.md`](TABLES.md), generated from `bands.json` and `conversion.json` by
`experiments/report_weak_hits.py`. Five seeds, 100,912 queries.

## The claim under test

The protein-index ladder established that bigger databases buy hits rather than
identifications: the residue's "any hit at all" climbed towards 100% across four
cumulative indexes while identifications stayed near 1% and were still tripling at the
largest. The standing reading of that was **headroom** — almost every residue read already
has a protein hit, so the prize is turning those weak hits into calls.

That reading needs two things the ladder did not have in the same pass: a shuffled null,
and the same molecule at two lengths. This run has both, over the query set of
[`../profile_reach/`](../profile_reach/REPORT.md) searched at a permissive reporting cutoff.

## 1. The headroom, as stated, is not there

**The weak-hit band on real library reads is indistinguishable from the band on
composition-matched shuffles.** Median bitscore 27.7 against 27.7, median identity 47.2%
against 47.7%, and 92.3% of *shuffled* reads produce a "hit" at the reporting cutoff at
all (`TABLES.md` §1).

So "any hit at all" is not a measure of anything about the read. It is the chance band,
and a curve of it against database size is a curve of how many chances the database
offers. This project has already withdrawn one claim for exactly this — a nucleotide BLAST
identity quoted without its E-value, where the chance band read identically on residue and
already-identified reads. The same lesson now holds for **translated search against a
protein index**, measured with the null in the same pass rather than inferred.

It goes further than the null. At 30% divergence the holdout's own weak hits (28.1 bits,
49.6% identity) are also indistinguishable from the shuffle's. A genuinely related read
that is far enough away produces a weak hit that looks exactly like noise, so **the weak
hit's own numbers cannot tell you which it is.** Something outside that hit has to decide.

## 2. Merging is something outside that hit, and it discriminates perfectly

Among pairs whose mate 1 produced only a weak hit, the fraction whose merged read produces
a significant one (`TABLES.md` §2):

| pile | n weak on mate 1 | converted |
|---|---|---|
| holdout, 10% divergence | 87 | 35.6 (27.5–45.0) |
| holdout, 20% divergence | 246 | 34.6 (29.9–40.5) |
| holdout, 30% divergence | 413 | 18.2 (16.2–21.5) |
| real library reads | 330 | 6.6 (5.4–8.0) |
| **composition-matched shuffle (null)** | **387** | **0.0 (0.0–0.0)** |

Already-significant hits are kept at 98.7–100% throughout, so nothing is traded away.

**Zero of 387 shuffled weak hits per seed convert, in 5 of 5 seeds.** That is what makes
this a gain rather than arithmetic, and the mechanism is in the last column of the table.
The length-arithmetic null — a merged read's bitscore over its mate's scaled by the length
ratio — is **at or below 1.0 everywhere, and it falls as the material gets less related**:
0.987 on the identical holdout, 0.952 at 20% divergence, 0.873 at 30%, 0.794 on real
library reads, 0.745 on shuffles. A merged read gains *less* score than its extra length
would predict, and least of all on noise.

That is the whole story. An E-value penalises a longer query, so a chance hit cannot be
promoted by lengthening the query — the penalty rises and the chance score does not. A
true alignment extends into the added sequence and clears the penalty. **Merging is
therefore not a longer ruler; it is a test that only real alignments pass.**

## 3. What to lean on, and what this costs to say honestly

- **Do not quote "any hit" against a protein index, at any database size, without the
  null.** 92% of shuffled reads have one.
- **The remaining read-level gain in identification is narrow and real.** On real library
  reads the conversion is 6.6% of the weak subset, which is about one point of the
  overlapping pairs — on top of the +5.3 points merging already buys at the call line
  (22.1% against 16.8% at E <= 1e-3, both on overlapping pairs).
- **It is free.** No new index, no new sequencing, no threshold to tune: the gain is
  measured at DIAMOND's own default significance cutoff.

Two numbers that must not be compared, because the denominators differ: 20.0% of *all*
background mate-1 reads are significant, against 16.8% of the *overlapping* subset.
Overlapping pairs are short inserts, which this library's residue is enriched for, and
they identify slightly less well. Every merged-against-mate figure here is within the
overlapping subset.

## 4. A correction to results/profile_reach/

That run searched with DIAMOND's **default** `--evalue 0.001`, so no weak hit was ever
reported and its `diamond_any_hit` field in `cells.json` counted only significant hits.
Nothing published from it depended on a weak-hit population — the reach ladder and the
merged-against-mate comparison are all at or above significance, and its bitscore floors
of 30 and 40 sit inside that cutoff rather than below it — but the field is misnamed for
what it holds, and the weak-hit question could not be asked of that run at all. This is
the same instrument at a second operating point, on a byte-identical query set.

## Provenance

| file | script | what |
|---|---|---|
| `bands.json` | `experiments/weak_hit_conversion.py` | §1, per seed, per pile, per mate |
| `conversion.json` | same | §2, paired within a molecule, with the length null |
| `TABLES.md` | `experiments/report_weak_hits.py` | generated from the two above |
| `manifest.json` | `experiments/weak_hit_conversion.py` | seeds, git SHA, script hash, operating points |

The search intermediates live in `results/profile_reach/work/` and are gitignored; the
manifest records the query set and both cutoffs.
