# Mate linkage: does a read find its own mate?

`experiments/mate_friends.py`, 5 seeds, 20,000 pairs each (39,9xx mate points), one script
hash across all five. **Every number below is in
[`TABLES.md`](TABLES.md)**, generated from the committed per-seed JSONs.

## Why this exists

`results/residue_structure/` §4 asked whether one genome's reads find each other in the
joint space, by planting a genome at a known abundance and measuring the fraction of a
planted read's 10 nearest neighbours that came from the same source. It works to about 1 in
100,000 of the library and dies there. That measurement needed a planted arm to know the
answer, it treated a PAIR as one point, and so it never used the one piece of ground truth
the real library already carries: **the two mates of a pair came from the same molecule.**

So here each MATE is a point, and the question is whether a read's own partner is among its
10 nearest neighbours. No planting, no simulation, real sequence, every abundance present in
the library — and it is the one lever the kin wall was never given.

## What it reads

| population | recovery at k=10 | × chance |
|---|---|---|
| permuted-partner null | 0.022% | 0.9× |
| **mates sharing no 31-mer** | **0.580%** | **23.2×** |
| mates that overlap | 6.576% | 262× |
| all mates | 2.483% | 99× |

Chance is 0.0250%. **Read the null first**: permuting which mate is whose partner lands on
chance, so the geometry is not leaking partner identity through anything but the axes.

31.7% of pairs share at least one canonical 31-mer between their mates, and those pairs
contain the same bases twice, so their coverage profiles agree by construction. The disjoint
row is therefore the honest figure, and the test is conservative in the safe direction: a
shared 31-mer can also come from repeat content, so the split over-assigns to the
overlapping half.

**Mate linkage carries real information about molecule identity: 23× chance with no shared
sequence at all, across 5 seeds with a 21–26× range.**

## It is a bound, not an estimate

Two mates share an insert, so they are far more alike than two reads a kilobase apart in the
same genome. This is an **upper bound** on friendship resolution between reads of one
genome, and the planted-genome purity in `results/residue_structure/` is the corresponding
lower bound. Quoting it as the grouping rate would be the error; stating it as a bound is
the point.

## What happens as the space fills, which is the part that decides anything

| pairs in the space | disjoint recovery | × chance |
|---|---|---|
| 2,496 | 3.044% | 15.2× |
| 4,992 | 1.661% | 16.6× |
| 9,984 | 0.969% | 19.3× |
| 19,968 | 0.579% | 23.1× |

The level falls roughly as 1/n while the multiplier **rises**, about 1.5× for every eightfold
increase in the number of points. Both halves matter and they point different ways:

- The rising multiplier says the axes genuinely locate the molecule rather than running a
  fixed-size lottery. A pure lottery would hold enrichment flat.
- The falling level says that at library scale a mate's partner is still not findable. One
  read in a 35-million-pair library does not find its own mate by neighbourhood, let alone
  find a stranger from the same genome.

**The far end is not extrapolated.** What the multiplier does at 35 million points is
unmeasured, and the measurement stops at 20,000 pairs because that is what was run.

## The two mates are not interchangeable, and here the difference is entirely length

| | mate 1 → mate 2 | mate 2 → mate 1 |
|---|---|---|
| as sequenced (107 nt and 93 nt) | 2.311% | 2.653% |
| mate 1 truncated to 93 nt | 4.093% | 4.087% |

Truncating mate 1 to mate 2's length collapses the asymmetry to the third decimal — 4.093
against 4.087, close but not identical, so this is an agreement rather than a duplicated
column. The direction difference in the raw row is length, not biology.

The second thing in that table is larger and was not the question: **matching the lengths
nearly doubles recovery**, 2.48% to 4.09% overall. Unequal mate lengths are themselves noise
in this space, because the coverage axes are computed over different numbers of scorable
31-mers. Anything downstream that compares two reads' coverage profiles should match their
lengths first.

## Where it works best is the opposite of where the kin wall is

| decile | log coverage | recovery | × chance |
|---|---|---|---|
| 1 | 0.69–1.10 | 2.73% | 109× |
| 2 | 1.10–1.61 | 4.08% | 163× |
| 3 | 1.61–2.30 | 3.99% | 159× |
| 5 | 3.66–4.99 | 2.13% | 85× |
| 10 | 9.58–12.76 | 1.04% | 42× |

Recovery is highest in the second and third coverage deciles and lowest at the top. That is
crowding, not biology: an abundant read sits in a dense region where thousands of reads are
equally close, and a rare one sits in a sparse corner where its partner is one of few
candidates. It is a property of the local density of this space, and at library scale the
sparse corner fills too, so it is not a claim that linkage beats the abundance wall.

The degenerate case runs the other way and is handled: mates with no witnesses anywhere
share an exact coordinate with thousands of others, so **ties are broken by a random
permutation before every query**. sklearn returns neighbours in index order, and a pair's two
mates are adjacent in any natural ordering, so index-order ties would hand every degenerate
mate its own partner and manufacture the entire result. That is the failure that inflated a
nearest-neighbour control in this project once already.

## Coverage and composition contribute about equally here

| axes | recovery | × chance |
|---|---|---|
| the four coverage axes alone | 1.183% | 47× |
| the six composition axes alone | 0.974% | 39× |

This is a contrast with the planted measurement, where kinship was mostly coverage (0.918
purity on coverage alone against 0.292 on composition alone at 1 in 10,000). Two mates of
one molecule share composition directly; two reads from distant parts of one genome share
abundance but need not share local composition. The gap between the two measurements is
itself the statement that this is an upper bound.

## What this changes

1. The grouping step in `hapax run` is calibrated on this, and the run directory carries
   `linkage.json` so no grouping can be read off the web view without its operating point.
2. Unequal mate lengths should be matched before any pairwise comparison of coverage
   profiles. That is a cheap change with a measured effect of the same size as the signal.
3. Mate linkage does not move the 1-in-100,000 wall. It is a real multiplier on the
   geometry and not a locator, and the honest version of "reads find their friends" keeps
   both facts attached to it.

## Operating point

31-mers, exact counts over the whole library by `experiments/kmer_coverage` (query-only path,
so a pair's own occurrences are not subtracted); mate 2 trimmed by 14 nt, the measured
non-template prefix; 10 axes, each standardised across the points in the space being
queried; Euclidean neighbourhoods with k = 1, 5, 10 and 50; ties permuted at random before
every query; blocks are the same seed-determined disjoint blocks the rest of the repo uses.
