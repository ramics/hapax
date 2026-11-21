# The linkage census: what the residue sits next to

`experiments/linkage_census.py` driving `experiments/kmer_coverage/src/bin/census.rs`,
5 seeds, 3,000 pairs per pile per seed. **Every number is in
[`TABLES.md`](TABLES.md)**, generated from the committed per-seed JSONs.

## The question, and why coverage could not answer it

Coverage says how often a read's sequence recurs in the library. It cannot say what that
sequence is attached to. A paired library answers adjacency for free: find the other pairs
covering the same sequence as a residue read, then look at THEIR other mates. Each of those
is a piece of sequence one insert away from something the residue read also covers. Label
them and the residue splits into two populations that mean quite different things:

- **unidentified because the reference is a gene catalogue** — rRNA precursor spacers, UTRs,
  intergenic sequence, the untranslated flanks of things we did identify;
- **a neighbourhood as unknown as the read** — whole unknown molecules, which is the
  population the rest of this project is about.

## What it reads

| the query sits... | residue | already identified | permuted-label null |
|---|---|---|---|
| beside sequence we can name | **3.27%** | 94.67% | 2.87% |
| beside contigs that assembled but nothing identified | **30.50%** | 0.94% | 56.97% |
| a neighbourhood as unknown as the read | **35.26%** | 3.27% | 7.79% |
| mixed | 0.77% | 0.30% | 2.18% |
| no witnesses at all | 30.19% | 0.82% | - |

Read the controls first.

- **The identified pile separates sharply**: 94.7% of already-subtracted reads sit beside
  sequence we can name, against 3.3% of residue reads. The instrument works.
- **The permuted-label null changes which rows survive.** `unknown_neighbourhood` reads 4.5
  times its null, so unknown neighbours genuinely cluster on the same queries — which is
  what "whole molecule" means rather than "scattered unknown k-mers". But `beside_known` at
  3.27% is barely above its null of 2.87%, so **the honest claim is that at most about 3%
  of the residue is unidentified for annotation reasons, and the null already accounts for
  most of that**. The row that reads below its null,
  `beside_assembled_unidentified` at 30.5% against 57.0%, is the same clustering seen from
  the other side: a query's neighbours agree with each other more than a random draw from
  the mix would.

## The most robust number here is the mix, not the piles

The piles depend on a half-of-neighbours rule. The neighbour label mix does not, and it is
where the two populations are furthest apart:

| neighbour class | residue | already identified |
|---|---|---|
| rRNA | 1.31% | 52.12% |
| identified contig | 6.77% | 41.74% |
| unidentified contig | 64.08% | 1.23% |
| none | 27.84% | 4.91% |

A residue read's neighbourhood is 92% sequence that is itself unidentified. An identified
read's neighbourhood is 94% sequence we can name. **The residue is not an annotation
artefact: it lives in a different part of the library.**

## What the census is blind to, said plainly

**30.2% of residue pairs have no witnesses at all** — nothing else in 8.7 million pairs shares
ten 31-mers with them. Those are reported as `no_witnesses` and are never folded into
`unknown_neighbourhood`, because a singleton molecule and a molecule with an unknown
neighbourhood are different claims and this method cannot tell them apart. The census is
found through shared sequence, so it is structurally blind to exactly the rarest material.
The identified pile's figure for the same row is 0.82%, which is the size of that blind
spot: it is almost entirely a residue phenomenon.

## The overlap control moves the numbers and not the direction

About a third of this library's pairs overlap their own mate, so for those the "other" mate
is partly the same sequence as the matching one. Dropping every witness whose two mates both
match the query:

| the query sits... | residue | already identified |
|---|---|---|
| beside sequence we can name | 4.17% | 93.61% |
| beside contigs that assembled but nothing identified | 26.55% | 0.97% |
| a neighbourhood as unknown as the read | 16.21% | 2.45% |
| no witnesses at all | 52.27% | 2.70% |

The unknown-neighbourhood pile halves and `no_witnesses` rises to 52.3%, because dropping
overlapping witnesses throws away most of a residue read's evidence. The separation from
the control survives intact, and `beside_known` actually rises. So the direction is robust
and the magnitude of the unknown pile is not: between 16% and 35% of the residue has a
measurably unknown neighbourhood, depending on whether overlapping witnesses are allowed to
count, and the rest of the residue is either beside assembled-unidentified sequence or
invisible to the method.

## What this changes

1. **The residue is mostly not an annotation artefact.** At most about 3.3% of it sits beside
   sequence we can name, and the null explains most of that. The old worry that the residue
   is largely untranslated flanks of known genes is not supported.
2. **A third of the residue is beside contigs that assembled without being identified**,
   which makes tier 2 the right place to look next for those reads: they have assembled
   neighbours, so they are not unassemblable.
3. **A third has no witnesses at all**, which is the population no neighbourhood method can
   reach, and it is six points larger than the whole-library singleton rate would suggest.

## Operating point

31-mers; witnesses share at least ten of them (about 40 bases, the cascade's own match
length); at most 30 witnesses per query, reservoir sampled over the whole library so the cap
is not a statement about one region of one lane; neighbours labelled by minimap2 at 40
aligned bases against the sample's own ribosomal contigs, then its identified contigs, then
its unidentified contigs, in that order, so the classes are exclusive and `none` means
unidentified at the same threshold that defined the residue; mate 2 trimmed by 14 nt.
