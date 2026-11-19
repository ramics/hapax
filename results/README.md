# Results

Everything here is reproducible from a committed script plus a seed. Each
directory has a `manifest.json` with the git SHA, parameters and timings.
These runs used commit `4e55e54`.

## What the simulation is (and is not)

`python/swarm.py`. There are four entities per world: two **sister populations**
whose roots are 0.30 subs/site apart, plus two unrelated families. Each entity is
a 32-tip coalescent swarm scaled to the target mean pairwise diversity. Selection
acts at the codon level: ω = 0.23, 10% of codons positive and 10% purifying, as in
the thesis EvolveAGene runs. Each genome is about 4.4 kb with three genes, a
300-nt dual-coding overlap, and UTRs. Reads are single-end 150 nt at 50× per
entity with modern Illumina error rates.

The sister populations are what make chimaeras possible. Unrelated families would
flatter any assembler.

Not yet covered:
- real data
- paired-end reads
- the external assemblers (MEGAHIT, metaSPAdes, PLASS, PenguiN). They are not
  installable in the development sandbox. The ladder runs them when they are on
  PATH (`--external`).

## 1. Diversity ladder: does codon scoring buy anything?

`experiments/diversity_ladder.py --seed 1 --n-seeds 5`. The figure is
`ladder/ladder.png`. Paired by seed in `ladder/paired.md`.

Codon mode and nucleotide mode run the **identical pipeline**. The only
difference is whether the tier-2 HMM has codon states. That makes this the clean
test of the thesis. Results at τ = 0.9 (mean ± sd over 5 seeds):

| diversity | config | contigs per found entity | entities missed | chimaera rate | recovered | HMM fraction | wall (s) |
|---|---|---|---|---|---|---|---|
| 0.05 | codon | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.71 ± 0.19 |
| 0.05 | codon-spaced | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.82 ± 0.45 |
| 0.05 | nucleotide | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.15 ± 0.18 |
| 0.05 | nucleotide-spaced | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.12 ± 0.09 |
| 0.1 | codon | 1.05 ± 0.11 | 0.20 ± 0.45 | 0.07 ± 0.15 | 0.94 ± 0.13 | 0.38 ± 0.04 | 2.21 ± 0.53 |
| 0.1 | codon-spaced | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.38 ± 0.04 | 2.11 ± 0.20 |
| 0.1 | nucleotide | 1.05 ± 0.11 | 0.20 ± 0.45 | 0.07 ± 0.15 | 0.94 ± 0.13 | 0.38 ± 0.04 | 1.45 ± 0.64 |
| 0.1 | nucleotide-spaced | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.38 ± 0.04 | 1.32 ± 0.46 |
| 0.15 | codon | 1.25 ± 0.31 | 0.00 ± 0.00 | 0.04 ± 0.09 | 1.00 ± 0.01 | 0.39 ± 0.03 | 2.32 ± 0.75 |
| 0.15 | codon-spaced | 1.10 ± 0.14 | 0.00 ± 0.00 | 0.09 ± 0.12 | 1.00 ± 0.01 | 0.39 ± 0.03 | 2.18 ± 0.41 |
| 0.15 | nucleotide | 1.35 ± 0.28 | 0.00 ± 0.00 | 0.04 ± 0.09 | 1.00 ± 0.00 | 0.39 ± 0.03 | 1.12 ± 0.04 |
| 0.15 | nucleotide-spaced | 1.20 ± 0.21 | 0.00 ± 0.00 | 0.08 ± 0.11 | 1.00 ± 0.01 | 0.39 ± 0.03 | 1.28 ± 0.45 |
| 0.2 | codon | 1.63 ± 0.41 | 0.20 ± 0.45 | 0.04 ± 0.09 | 0.94 ± 0.13 | 0.45 ± 0.06 | 2.49 ± 1.04 |
| 0.2 | codon-spaced | 1.65 ± 0.28 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.45 ± 0.06 | 2.09 ± 0.16 |
| 0.2 | nucleotide | 1.70 ± 0.37 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.45 ± 0.06 | 1.36 ± 0.57 |
| 0.2 | nucleotide-spaced | 1.85 ± 0.28 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.45 ± 0.06 | 1.13 ± 0.06 |

The wall times are noisy because other jobs shared the CPU. Use
`scaling/` for timing.

**Findings**

1. **The codon gain is real but small, and it only appears at high diversity.**
   Up to 0.10 subs/site, both modes give essentially one contig per entity, so
   there is nothing to gain. At 0.15 to 0.20, codon mode gives 0.05 to 0.35
   fewer contigs per entity than nucleotide mode on the same worlds. With spaced
   seeds it is never worse on any seed (`paired.md`), and at 0.20 / τ = 1.0
   it is better on 4 of 5 seeds. That is consistent, but it is a long way from
   "one population, one contig".
2. **With contiguous seeds at 0.20, codon mode added one chimaera on one seed out
   of five.** This is the METHOD.md §5 risk showing up. The numbers are too
   small to call it a trend, but the effect points the wrong way.
3. **Seeding matters more than scoring.** The codon-spaced seed, which ignores
   wobble positions, removed every entity miss and every 0.10 chimaera that
   contiguous seeds produced, and it did so in *both* modes. Most of the
   improvement available so far comes from the index rather than the HMM. That
   matches CLAUDE.md's point that indexing is the bigger lever.
4. **The remaining fragmentation is not an overlap-scoring problem.** At 0.20,
   the extra contigs are mostly deep clades of one swarm. In the coalescent tree
   they sit as far apart as the sister populations do. They connect to the rest
   of their entity by a handful of overlaps, while the sister populations connect
   by one or two. Telling those apart is a population-structure question. A
   per-pair score cannot answer it.
5. **Kill criterion "too slow" is borderline.** The HMM fraction is 0.24 at 0.05
   and 0.45 at 0.20 diversity: not "most overlaps", but close to half. Codon
   mode costs about 1.8× nucleotide mode end to end at 50×.

## 2. Scaling

`experiments/scaling.py --seed 1 --n-seeds 5`. Same world, rising coverage,
diversity 0.10, 4 threads, mean of 5 seeds.

| coverage | mode | reads | HMM fraction | candidates (s) | overlap (s) | consensus (s) | wall (s) |
|---|---|---|---|---|---|---|---|
| 25 | codon | 2942 | 0.502 | 0.02 | 0.65 | 0.42 | 1.12 ± 0.04 |
| 25 | nucleotide | 2942 | 0.502 | 0.02 | 0.11 | 0.42 | 0.58 ± 0.01 |
| 50 | codon | 5888 | 0.361 | 0.07 | 1.04 | 0.84 | 2.00 ± 0.08 |
| 50 | nucleotide | 5888 | 0.361 | 0.07 | 0.19 | 0.84 | 1.15 ± 0.02 |
| 100 | codon | 11778 | 0.187 | 0.24 | 1.16 | 1.63 | 3.12 ± 0.19 |
| 100 | nucleotide | 11778 | 0.187 | 0.23 | 0.20 | 1.61 | 2.13 ± 0.08 |
| 200 | codon | 23554 | 0.097 | 0.84 | 1.23 | 3.31 | 5.55 ± 0.18 |
| 200 | nucleotide | 23554 | 0.097 | 0.84 | 0.23 | 3.36 | 4.60 ± 0.04 |
| 400 | codon | 47119 | 0.048 | 3.40 | 1.35 | 6.63 | 11.74 ± 0.12 |
| 400 | nucleotide | 47119 | 0.048 | 3.35 | 0.26 | 6.70 | 10.65 ± 0.24 |

**Findings**

- Total wall clock is close to linear in read count. Consensus and overlap
  scale linearly. Overlap time stays flat because the partner cap per read keeps
  the near-identical neighbours.
- **Candidate generation is not linear.** Reads grow 16×, but candidate time
  grows about 170×, because every minimizer bucket grows with depth. This is the
  "fast quadratic assembler" CLAUDE.md warns about, and it is fixed next.
- The HMM fraction *falls* with depth. The partner cap fills up with
  near-identical reads that tier 1 settles. So the speed kill-criterion number
  depends on coverage and should always be quoted with it.

---

## 3. Dispatch fixes, v2 (commit `5f8a7e2`)

`ladder_v2/` and `scaling_v2/` hold the results after two changes to the
dispatch.

- **Bucket cap.** For each read, each seed bucket now returns only the 48
  entries per strand whose seed positions are nearest the read's own. For
  fixed-length reads those are the partners with the largest overlaps, and they
  are the nearest partners in every bucket the pair shares, so their anchors
  still add up.
- **Lazy tier 2.** An ambiguous pair goes to the HMM only if its two reads are
  in different tier-1 components (clusters joined by near-identical,
  indel-free overlaps). A pair already joined through those overlaps cannot
  change the grouping, so it keeps its cheap ungapped estimate as a "tier 3"
  edge.

A first version dropped those skipped pairs entirely. On 1 of 5 seeds that let a
weak cross-population edge into a read's top k, which merged the sister
populations. Keeping the pairs as tier-3 edges fixed it, and the fix is in this
run.

The v2 worlds differ slightly from v1 because the simulator's indel code
changed. The `-v0` configs rerun the old dispatch on the same v2 worlds, so
each comparison against v0 is paired. Results at τ = 0.9, mean ± sd over 5 seeds:

| diversity | config | contigs per found entity | entities missed | chimaera rate | recovered | HMM fraction | wall (s) |
|---|---|---|---|---|---|---|---|
| 0.05 | codon | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 1.05 ± 0.05 |
| 0.05 | codon-spaced | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 1.01 ± 0.03 |
| 0.05 | codon-spaced-t1groups | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 1.02 ± 0.02 |
| 0.05 | codon-spaced-v0 | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.69 ± 0.17 |
| 0.05 | nucleotide | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 0.98 ± 0.02 |
| 0.05 | nucleotide-spaced | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 0.98 ± 0.01 |
| 0.05 | nucleotide-spaced-t1groups | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.02 ± 0.01 | 0.99 ± 0.03 |
| 0.05 | nucleotide-spaced-v0 | 1.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.24 ± 0.06 | 1.09 ± 0.04 |
| 0.1 | codon | 1.05 ± 0.11 | 0.20 ± 0.45 | 0.07 ± 0.15 | 0.95 ± 0.11 | 0.05 ± 0.01 | 1.08 ± 0.07 |
| 0.1 | codon-spaced | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.05 ± 0.01 | 1.09 ± 0.05 |
| 0.1 | codon-spaced-t1groups | 1.15 ± 0.22 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.05 ± 0.01 | 1.12 ± 0.07 |
| 0.1 | codon-spaced-v0 | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.38 ± 0.04 | 2.03 ± 0.11 |
| 0.1 | nucleotide | 1.05 ± 0.11 | 0.20 ± 0.45 | 0.07 ± 0.15 | 0.95 ± 0.11 | 0.05 ± 0.01 | 0.97 ± 0.04 |
| 0.1 | nucleotide-spaced | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.05 ± 0.01 | 0.98 ± 0.02 |
| 0.1 | nucleotide-spaced-t1groups | 1.20 ± 0.33 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.05 ± 0.01 | 0.99 ± 0.01 |
| 0.1 | nucleotide-spaced-v0 | 1.05 ± 0.11 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.38 ± 0.04 | 1.15 ± 0.03 |
| 0.15 | codon | 1.25 ± 0.31 | 0.00 ± 0.00 | 0.04 ± 0.09 | 1.00 ± 0.01 | 0.09 ± 0.03 | 1.20 ± 0.08 |
| 0.15 | codon-spaced | 1.10 ± 0.14 | 0.00 ± 0.00 | 0.09 ± 0.12 | 1.00 ± 0.00 | 0.09 ± 0.03 | 1.19 ± 0.11 |
| 0.15 | codon-spaced-t1groups | 1.40 ± 0.22 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.09 ± 0.03 | 1.21 ± 0.11 |
| 0.15 | codon-spaced-v0 | 1.10 ± 0.14 | 0.00 ± 0.00 | 0.09 ± 0.12 | 1.00 ± 0.01 | 0.39 ± 0.03 | 1.99 ± 0.07 |
| 0.15 | nucleotide | 1.35 ± 0.28 | 0.00 ± 0.00 | 0.04 ± 0.09 | 1.00 ± 0.00 | 0.09 ± 0.03 | 0.98 ± 0.06 |
| 0.15 | nucleotide-spaced | 1.20 ± 0.21 | 0.00 ± 0.00 | 0.08 ± 0.11 | 1.00 ± 0.00 | 0.09 ± 0.03 | 0.97 ± 0.06 |
| 0.15 | nucleotide-spaced-t1groups | 1.55 ± 0.21 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.09 ± 0.03 | 0.99 ± 0.03 |
| 0.15 | nucleotide-spaced-v0 | 1.20 ± 0.21 | 0.00 ± 0.00 | 0.08 ± 0.11 | 1.00 ± 0.01 | 0.39 ± 0.03 | 1.10 ± 0.06 |
| 0.2 | codon | 1.63 ± 0.41 | 0.20 ± 0.45 | 0.04 ± 0.09 | 0.94 ± 0.12 | 0.15 ± 0.02 | 1.31 ± 0.10 |
| 0.2 | codon-spaced | 1.65 ± 0.28 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.15 ± 0.02 | 1.34 ± 0.09 |
| 0.2 | codon-spaced-t1groups | 2.10 ± 0.45 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.15 ± 0.02 | 1.34 ± 0.09 |
| 0.2 | codon-spaced-v0 | 1.65 ± 0.28 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.45 ± 0.06 | 2.13 ± 0.20 |
| 0.2 | nucleotide | 1.70 ± 0.37 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.15 ± 0.02 | 1.01 ± 0.01 |
| 0.2 | nucleotide-spaced | 1.80 ± 0.27 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.15 ± 0.02 | 1.00 ± 0.01 |
| 0.2 | nucleotide-spaced-t1groups | 2.25 ± 0.47 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.15 ± 0.02 | 1.01 ± 0.02 |
| 0.2 | nucleotide-spaced-v0 | 1.85 ± 0.28 | 0.00 ± 0.00 | 0.00 ± 0.00 | 1.00 ± 0.00 | 0.45 ± 0.06 | 1.14 ± 0.05 |

**Findings**

1. **Lazy tier 2 costs nothing in accuracy and removes most of the HMM work.**
   `codon-spaced` matches `codon-spaced-v0` exactly on contigs, misses,
   chimaeras and recovery at every diversity. The HMM fraction falls from
   0.24–0.45 to 0.02–0.15, and codon mode's overhead over nucleotide mode falls
   from about 1.9× to about 1.3× end to end. The §5 "too slow" kill criterion is
   now comfortably met. At 0.20 diversity and 50× coverage, 15% of candidates
   reach the HMM.
2. **The codon gain is unchanged (`ladder_v2/paired.md`).** With spaced
   seeds, codon mode gives 0.10 to 0.35 fewer contigs per entity at 0.15 to 0.20
   diversity, is never worse on any seed, and adds no chimaeras. The gain is
   small but consistent.
3. **Using tier-1 components as the first-stage clusters (`-t1groups`) gives
   zero chimaeras at every diversity, but more contigs** (2.1 vs 1.65 at 0.20).
   Best-k clustering lets a few weak edges through at 0.15, where one seed
   produced a chimeric contig in every best-k config, codon and nucleotide
   alike. This is exactly the contiguity-versus-chimaera trade CLAUDE.md asks
   to be reported. Best-k stays the default for now, and the choice is open.
4. **Candidate generation is better but still not linear.** Per read it costs
   8 µs at 25× and 43 µs at 400×. In v1 it rose from 7 to 72 µs. The rest of
   the pipeline is linear. The next step is profiling (likely index-size cache
   effects and the global candidate sort), not another heuristic.

| coverage | mode | reads | HMM fraction | candidates (s) | overlap (s) | consensus (s) | wall (s) | candidates µs/read |
|---|---|---|---|---|---|---|---|---|
| 25 | codon | 2942 | 0.130 | 0.02 | 0.18 | 0.41 | 0.63 ± 0.04 | 8 |
| 25 | nucleotide | 2942 | 0.130 | 0.02 | 0.03 | 0.42 | 0.50 ± 0.03 | 8 |
| 50 | codon | 5888 | 0.066 | 0.08 | 0.19 | 0.84 | 1.15 ± 0.05 | 13 |
| 50 | nucleotide | 5888 | 0.066 | 0.07 | 0.04 | 0.83 | 0.99 ± 0.02 | 12 |
| 100 | codon | 11778 | 0.024 | 0.25 | 0.15 | 1.64 | 2.13 ± 0.16 | 21 |
| 100 | nucleotide | 11778 | 0.024 | 0.25 | 0.04 | 1.64 | 2.01 ± 0.08 | 21 |
| 200 | codon | 23554 | 0.013 | 0.78 | 0.19 | 3.37 | 4.51 ± 0.05 | 33 |
| 200 | nucleotide | 23554 | 0.013 | 0.77 | 0.05 | 3.36 | 4.35 ± 0.04 | 33 |
| 400 | codon | 47119 | 0.012 | 2.05 | 0.39 | 6.72 | 9.48 ± 0.26 | 43 |
| 400 | nucleotide | 47119 | 0.012 | 2.08 | 0.10 | 6.67 | 9.17 ± 0.17 | 44 |

## Where this leaves the kill criteria (METHOD.md §5)

| criterion | status |
|---|---|
| No gain over nucleotide assemblers | **Weak gain, internal baseline only.** Codon scoring beats the same pipeline in nucleotide mode by 0.1 to 0.35 contigs per entity, and only at ≥ 0.15 diversity. It has not yet been compared against MEGAHIT or the haplotype assemblers (SAVAGE, Haploflow, VStrains), which the stress-test thread is covering. |
| Collapsing erases differences in overlapping genes | **Not yet tested properly.** The simulator has a dual-coding region, and the HMM falls back to nucleotide scoring when two offsets both win. Variants keep their nucleotide calls in the output. The real test is HIV *gag/pol* and HBV. |
| Too slow | **Passed after v2.** At most 15% of candidates reach the HMM, and codon mode costs about 1.3× nucleotide mode. |
| Chimaeras increase | **No increase from codon scoring in v2.** Chimaeras come from the grouping step and appear equally in both modes. Tier-1 grouping removes them at a contiguity cost. |
