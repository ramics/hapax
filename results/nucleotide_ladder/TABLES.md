Seeds: 5; 100,000 pairs per seed; database k2_pluspf16 (16.1 GB).

Axis is REFERENCE CONTENT, not bytes: the library was classified once against the largest database available and nested random subsets of the reference taxa were then applied to the calls. Both biases run the same way and are stated in the script's docstring; the stronger one is that deleting calls understates what a real smaller database would remove, which flatters the left-hand side and makes this a conservative test of whether the curve has flattened.

## Kraken2 confidence 0.0

| reference content | taxa | removed % of library | marginal | removed % of non-rRNA | marginal |
|---|---|---|---|---|---|
| 0.3% | 6 | 0.05 (0.01-0.07) | 0.05 (0.01-0.07) | 0.01 (0.00-0.02) | 0.01 (0.00-0.02) |
| 1% | 20 | 0.24 (0.10-0.54) | 0.20 (0.04-0.53) | 0.03 (0.01-0.05) | 0.02 (0.00-0.04) |
| 3% | 57 | 1.00 (0.32-3.09) | 0.75 (0.11-2.74) | 0.08 (0.05-0.11) | 0.05 (0.03-0.06) |
| 10% | 202 | 3.68 (2.06-5.23) | 2.68 (1.52-4.58) | 1.31 (0.38-4.84) | 1.23 (0.27-4.77) |
| 30% | 642 | 9.42 (6.14-12.34) | 5.75 (3.86-7.22) | 2.66 (1.02-5.97) | 1.35 (0.63-2.60) |
| 50% | 1086 | 16.88 (10.53-22.41) | 7.46 (3.26-16.27) | 4.65 (1.86-7.42) | 1.99 (0.84-5.85) |
| 100% | 2155 | 42.82 (42.62-43.00) | 25.94 (20.54-32.08) | 11.61 (11.31-11.83) | 6.96 (4.41-9.45) |

Read concentration over reference taxa: top taxon 14.9%, top 10 51.9%, top 100 81.0% of classified reads. The tail is long and the head is heavy, so a content fraction is a weak proxy for what a subset buys: whether one dominant taxon falls inside the subset moves a rung more than hundreds of rare ones, which is why the spread across seeds is wide and is reported rather than averaged away.

## Kraken2 confidence 0.1

| reference content | taxa | removed % of library | marginal | removed % of non-rRNA | marginal |
|---|---|---|---|---|---|
| 0.3% | 1 | 0.00 (0.00-0.02) | 0.00 (0.00-0.02) | 0.00 (0.00-0.00) | 0.00 (0.00-0.00) |
| 1% | 2 | 0.05 (0.00-0.16) | 0.04 (0.00-0.16) | 0.00 (0.00-0.00) | 0.00 (0.00-0.00) |
| 3% | 9 | 0.15 (0.03-0.61) | 0.11 (0.01-0.45) | 0.01 (0.00-0.02) | 0.01 (0.00-0.01) |
| 10% | 31 | 1.49 (0.19-3.20) | 1.34 (0.15-2.59) | 0.59 (0.04-2.58) | 0.58 (0.03-2.58) |
| 30% | 99 | 2.60 (0.73-5.19) | 1.11 (0.45-1.99) | 0.89 (0.13-2.83) | 0.29 (0.09-0.57) |
| 50% | 166 | 6.37 (1.58-16.00) | 3.77 (0.39-15.28) | 1.72 (0.25-3.55) | 0.83 (0.12-3.27) |
| 100% | 338 | 21.67 (21.63-21.73) | 15.30 (5.72-20.10) | 4.30 (4.14-4.49) | 2.59 (0.95-3.90) |

Read concentration over reference taxa: top taxon 49.0%, top 10 86.4%, top 100 98.1% of classified reads. The tail is long and the head is heavy, so a content fraction is a weak proxy for what a subset buys: whether one dominant taxon falls inside the subset moves a rung more than hundreds of rare ones, which is why the spread across seeds is wide and is reported rather than averaged away.

## Does it flatten?

- At confidence 0.0, the last doubling of reference content (50% to 100%) still buys 25.9 points of the library, and no plateau is visible. Read that next to the concentration above, though: most of that step is one dominant taxon being present rather than absent, so this ladder is a weak instrument for the SHAPE of the curve and a strong one only for the claim that it has not visibly saturated.
- At confidence 0.1, the last doubling of reference content (50% to 100%) still buys 15.3 points of the library, and no plateau is visible. Read that next to the concentration above, though: most of that step is one dominant taxon being present rather than absent, so this ladder is a weak instrument for the SHAPE of the curve and a strong one only for the claim that it has not visibly saturated.

The stronger nucleotide evidence that the curve has not flattened is not this ladder but the remote BLAST in the main report: against NCBI core_nt, two orders of magnitude more sequence than anything local, 87% of residue reads find a partial, distant match. Sequence our indexes miss is findable in a database 70 times larger.

This matches the protein side (PR #10 sections 4.6-4.7), where each tenfold step in index size was still buying about 20 points at the largest index reachable. So 'removal per gigabyte and where it flattens' has, so far, no flattening point on either axis, and this report does not invent one.
