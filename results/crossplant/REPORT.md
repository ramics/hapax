# Cross-plant backgrounds: how different they are, and what a union of them leaves

**Two questions from Imo (2026-11-22):** how different are the wastewater backgrounds from
each other, and does subtracting a union of them all from held-out samples leave much
behind. Plus a third point: if something in the subtracted background stops being a
background and shows appreciable coverage above the panel, report the outsize growth instead
of subtracting it.

**Answers, one line each.** Plant identity is real but small: a read recurs in another water
of its own plant at 55.8% and in a plant on another continent at 37.6%, a paired gap of
+14.3 points that holds in all 32 plants. A union of 512 other libraries removes 87.5% of a
held-out library's reads and is **still falling at 512**, so it has not saturated — but what
it removes is abundance, not identity: every pile loses its deep sequence to the union and
keeps its depth-1 sequence, and *identified* reads at depth 1 survive at 70.5% against the
unplaced residue's 79.0%. Outsize growth is measurable and well separated from its control,
but what it finds is sequence most of the panel also carries blooming locally, not something
leaving the background.

## 1. What was measured

531 subject libraries, each cut to exactly **2,000,000 read pairs** so a difference between
plants is not a difference in library size (`controls.json`,
`subject_pairs_min_max: [2000000, 2000000]`). 530 come from the panel in `panel.tsv` — five
projects, 32 plants, four continents — and the 531st is SRR23973193 itself, held out and
never part of any union.

Witnessing is exact canonical 31-mer presence, the same instrument as
`results/recurrence/`. Two operating points throughout: **`any`** (at least one of a query's
31-mers present in the subject) and **`most`** (half or more). Nothing is thresholded into a
verdict; the rates are the result.

**The query axis is sampled, the subject axis is not.** Queries come from 95 libraries, three
per plant at a stride over that plant's dates, plus the Berlin piles: 145,040 rows in all
(`query_build.json`). All 530 panel libraries remain subjects, because the matrix and the
union are what need the full list. The query axis had to be sampled because the scan's cost
is not linear in the query count — see §6.

Every query is truncated to **90 nt per mate**, so a 151 nt library's reads offer no more
31-mers than a 101 nt library's. Reads are reported **per mate** throughout, as the repo
requires.

## 2. Controls, read first

| control | mate 1 | mate 2 | what it bounds |
|---|---|---|---|
| shuffled query reads in other libraries (`most`) | **0.0000** | **0.0016** | the floor |
| own library, spots outside its own subject (`most`), median of 95 libraries | **0.7100** | 0.6800 | the ceiling |
| own library, same, worst library | 0.0900 | 0.1000 | |
| own library, same, best library | 0.9800 | 0.9700 | |

The floor is clean: mononucleotide shuffles of real reads are witnessed essentially never, so
nothing below is composition. Mate 2's 0.0016 against mate 1's 0.0000 is small but real and
in the expected direction — mate 2 is the one carrying non-template sequence (§5).

The ceiling matters for reading §3. Queries are drawn from spots **outside** the 2M-pair
subject of their own library, so the own-library rate is an honest measure of how often this
instrument re-finds a read at this depth: typically 0.71, but from 0.09 to 0.98 across
libraries. An earlier version of this table quoted only the worst library, which read 0.09
beside a within-plant rate of 0.56 and looked like a cross-library rate above its own
ceiling. It was the minimum of 95 libraries, not a ceiling.

## 3. How different are the backgrounds?

`plant_summary.tsv`, median over five seeds, paired per query plant.

| stratum | `most1` | `most2` | plants |
|---|---|---|---|
| own library, other spots (ceiling) | 0.7100 | 0.6800 | 95 libs |
| **within plant, other water** | **0.5577** | 0.5497 | 32 |
| other plant, same project | 0.4986 | 0.4802 | 30 |
| **other plant, other project** | **0.3761** | 0.3620 | 32 |
| shuffled (floor) | 0.0000 | 0.0016 | — |

A monotone ladder between a measured floor and a measured ceiling. The ordering is what you
would expect and the *size* is the finding: most of what is in one plant's water is in
everyone's.

**Paired, not pooled.** The within-plant rate itself runs from 0.1530 (Seocheon, Korea) to
0.9593 (plant A, Ohio) — a six-fold spread — so a pooled within-minus-between difference
would be comparing plants rather than comparing within against between. Each query plant is
therefore its own comparison: the paired within-minus-other-project gap is **+0.1432**
(`most1`, range +0.0500 to +0.3791) and **+0.1485** (`most2`), and it is **positive in 32 of
32 plants**. The two mates agree to within 0.005, which is the check that this is not a
mate-length artefact.

**Direction of the bias, and the control that bounds it.** Subject read length is not matched
across projects (101 nt SoCal, 106 Berlin, 118 Ohio, 147 India, 151 Korea and
PRJNA1543313), and a longer subject read offers more 31-mers per pair at fixed depth, so it
should witness slightly more. It does, in that direction: Pearson r = +0.071 between a
plant's median read length and the rate it receives, and plants under 120 nt receive 0.4147
against 0.4518 for plants at 120 nt or more — a 3.7-point effect against the 14.3-point
effect being measured. At a fixed 151 nt the received rate still spans 0.3234 (Korea) to
0.4927 (SoCal SJ), four times the length effect, so length is a nuisance here and not the
story.

## 4. Does a union of them all leave much behind?

### A held-out library against a growing union

`saturation.tsv`, `heldout = ALL`, `most1`: the fraction of a held-out library's own reads
left unwitnessed by a union of k other libraries.

| pool | k=1 | k=4 | k=16 | k=64 | k=256 | k=512 |
|---|---|---|---|---|---|---|
| own plant | 0.3650 | 0.2550 | 0.1825 | 0.1050 | — | — |
| other plants | 0.5950 | 0.4000 | 0.3150 | 0.2350 | 0.1650 | 0.1450 |
| everything | 0.5850 | 0.4050 | 0.3000 | 0.2450 | 0.1550 | **0.1250** |

**It does not saturate.** The last three doublings still pay: 0.1850 at k=128, 0.1550 at
k=256, 0.1250 at k=512. Adding plants keeps removing material, slowly, with no knee.

Libraries from the held-out library's **own plant** are worth about 2.6 times more each: at
k=128 an own-plant union leaves 0.0800 where an other-plant union leaves 0.2100.

Per project at k=512 (`pool = all`, `most1`): Korea **0.3225**, Berlin 0.2750,
PRJNA1543313 0.1550, Ohio 0.0625, SoCal 0.0550, India 0.0550. Read this with the pool's
composition in mind rather than as plant distinctiveness: SoCal contributes 273 of the 530
panel libraries, so a held-out SoCal library has many near-siblings in any union, while
Korea contributes 32 and Berlin 6. The own-plant row above is the version of this that is
not confounded by pool composition.

### Our own sample's piles against the same union

`saturation.tsv`, `heldout = berlin:<pile>|<depth bin>`, `pool = all`, `most1`. The depth bin
is the read's occurrence depth in its own library.

| pile | depth 100+ | depth 10–99 | depth 3–9 | depth 1–2 | depth 0 |
|---|---|---|---|---|---|
| rRNA | **0.0000** | 0.0050 | 0.0650 | 0.1700 | 0.2850 |
| identified | 0.1700 | 0.2050 | 0.3750 | 0.5550 | **0.7050** |
| unidentified contig (tier 2) | 0.1000 | 0.1750 | 0.4300 | 0.6650 | 0.7100 |
| unplaced residue (tier 3) | **0.0000** | 0.2150 | 0.6150 | 0.7950 | **0.7900** |
| every shuffled pile | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

**Depth decides this, and the pile label barely matters.** Read down a column rather than
across a row. At depth 100+ the union takes everything, including all of the abundant part
of our own unplaced residue. At depth 1 it takes almost nothing from any pile — and the
critical comparison is **identified 0.7050 against unplaced 0.7900**. The union is not
failing on the residue because the residue is novel; it fails on anything present once,
including sequence we can already name, and the whole known-versus-unknown difference at
matched depth is 8.5 points.

The shuffled rows reading 1.0000 at every union size is the floor doing its job: 530
libraries of real wastewater remove **nothing at all** from composition-matched nonsense.

## 5. Outsize growth, reported instead of subtracted

A read's own-library 31-mer count against the panel's 99th percentile for that same read,
both mates, taking the smaller of the two ratios. Because queries come from spots outside
their library's subject, the own-library count does not include the read itself.

The control is a **swap**: a random panel library stands in for the own library and is scored
against the rest of the panel identically. If the own library were just another panel draw,
that is what growth would read.

| | ≥2× | ≥5× | ≥10× | ≥100× | p99 |
|---|---|---|---|---|---|
| real | **0.1143** | 0.0674 | 0.0470 | 0.0122 | **142.15** |
| swap control | 0.0013 | 0.0002 | 0.00002 | 0.0000 | 0.96 |

48,716 of 95,000 library queries qualify as background. The test fires 88 times harder than
exchangeability predicts at ≥2× and the control stays flat at p99 0.96, so this is a real
readout and not a restatement of library-to-library variation. (An earlier version of it was
not: queries were drawn from the very spots that formed the subject, so every read witnessed
itself, and it read 54% against a 6% control. The fix was to draw queries from a separate
pool of spots.)

**What it finds, stratified by how much of the panel carries the read at all:**

| panel occupancy | background reads | capped | ≥2× | ≥100× |
|---|---|---|---|---|
| 10–50% | 2,125 | 0.005 | 0.1144 | 0.0009 |
| ≥50% | 46,591 | 0.048 | 0.1143 | **0.0127** |

The two strata are indistinguishable at ≥2× and differ 14-fold at ≥100×: **the extreme
growth lives in sequence that most of the panel also carries.** That is a bloom of shared
material, which is the "already seen" direction this repo has twice been caught by, not
something leaving the background. The under-10% stratum is empty **by construction** — a
read must be witnessed in at least 10% of the panel to count as background at all — and that
gate is right rather than a gap: a read cannot outgrow a background it was never part of.
The honest middle case is the 10–50% band, and there the ≥100× tail is 0.09%.

Counts are capped at 65535 by `crossplant_reduce.py`, so `growth_top.tsv` carries
`own_count_capped` and the ranking above the cap is not meaningful — the top 500 rows are all
capped and their ratios are lower bounds. The cap touches 0.5% of the mid-occupancy stratum
and 4.8% of the high one, well under the 1.22% of reads reading ≥100×, so the pooled tail is
not a cap artefact.

## 6. What the pass cost, and the trim

**The mate-2 trim, decided from data in the same pass.** Berlin's SMARTer Pico protocol
leaves ~14 nt of non-template sequence on mate 2 and every experiment here trims it; no other
project uses that kit, so all 530 were run untrimmed and the 5′ base composition was measured
alongside (`composition_5p.tsv`, median max-base fraction over positions 0–14 minus positions
20–29):

| project | mate-2 excess | trim applied |
|---|---|---|
| SRR23973193 | +0.294 | 14 |
| PRJNA948850 (Berlin) | +0.267 | 14 |
| PRJNA924011 (Ohio) | **+0.134** | 0 |
| PRJNA1543313 | +0.088 | 0 |
| PRJNA729801 (SoCal) | +0.072 | 0 |
| PRJNA941081 (Korea) | +0.044 | 0 |
| PRJNA842541 (India) | +0.013 | 0 |

Mate 1 carries about +0.07 everywhere, so +0.07 is the baseline, not zero. Berlin's +0.267 is
far above it and the trim is justified. **Ohio's +0.134 is about double the baseline and was
not trimmed**, so Ohio's rates carry an untrimmed 5′ bias on mate 2; it is a caveat on that
project's row, not on the matrix as a whole.

**The scan has a cliff in the query count, not a rate.** 201,526 queries scanned a 2M-pair
library in 15–44 s; 776,826 queries, same libraries, took 1,400–1,800 s — 40× the time for
3.9× the queries, measured both ways on the same subjects. The first full pass was launched
at the larger size, needed ~74 hours, and was killed at two. The launcher now refuses a
query file over `--query-ceiling` rows. A separate lesson: the pilot's 15–44 s came from eight
libraries whose subject files were still in page cache, so it measured processor time with
the reads free and was not a rate to plan from; at 531 libraries nothing is cached and the
real figure at 145,040 queries is 265 s per subject on eight shards.

Three EC2 runs, about $5 in all. Eight shards on eight vCPU beat three by 1.8× in aggregate
throughput even though each shard individually slowed, and each shard's subject order is
shuffled under a fixed seed so that a run stopped by its budget yields a smaller panel spread
over every project and date rather than the first two projects and none of the rest.

## 7. Provenance

| file | what |
|---|---|
| `panel.tsv`, `panel_manifest.json` | the 536-row panel; `experiments/crossplant_panel.py` |
| `plant_matrix.tsv` | every plant × plant cell, four operating points |
| `plant_summary.tsv` | §3, the three strata and the paired difference |
| `saturation.tsv` | §4, both the held-out libraries and our own piles |
| `growth_top.tsv` | §5, the 500 largest ratios with their occupancy and cap flag |
| `composition_5p.tsv` | §6, the 5′ composition of all 531 libraries |
| `controls.json` | §2 and §5, floor, ceiling and swap control |
| `report_manifest.json` | seeds, git SHA, `script_sha256` |

Scripts: `experiments/crossplant_panel.py` builds the panel,
`crossplant_sample.py` cuts each library to a fixed depth and draws queries,
`crossplant_pass.py` scans, `crossplant_reduce.py` packs, `crossplant_report.py` reports.
`experiments/crossplant_contig_fixture.py` plants a run with known answers and the report
returns them — a within-plant 0.596 against a planted 0.60, an off-diagonal 0.150 against
0.15, zero on the shuffles, growth firing at every read while its swap control stays flat.
It was built before the real data landed and it caught two bugs in this report: a field
compared as a string that is an integer in memory, which silently found zero plants for the
paired comparison, and the own-library ceiling quoting only its minimum.

**Caveats.** Every number is at `most` or `any` on exact 31-mers: this instrument sees
identity, not homology, so a divergent relative of a panel sequence counts as absent.
Percentages are of query reads, not of a library. The five seeds resample the query draw, so
the min-max is an honest range over that draw and not a dispersion estimate. Drift applies
here as everywhere: the panel spans 2020 to 2025 and SRR23973193 is one morning in Berlin.

**The tier 2 contig arm runs over this same panel at the same fixed depth** and is reported
in §8 when it lands; contigs are the stronger probe, because a contig offers many independent
windows of one molecule where a 90 nt read offers close to a coin flip.
