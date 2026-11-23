# hapax — production system spec

**The goal, in Imo's words (2026-11-23): flag a novel virus that might be pointed at
humans.** This document is what the measurements of the last weeks say a system aimed at
that goal can and cannot be. It is an account of what was measured, not a design sketch:
every number in it is reproducible from a committed script plus a seed, and anything that
has not been measured is in §12 marked as open rather than asserted anywhere else.

The system emits **three output classes on one distance continuum**, and for each class it
states what it reports and what it refuses to claim:

| class | what it is | §|
|---|---|---|
| **A — divergent but related contigs** | assembled, identified, at real distance from anything published | §5 |
| **B — novel strains of known viruses** | the near end of the same continuum; a first-class case, not an edge case | §6 |
| **C — the unplaceable residue** | reads and contigs nothing names, at median depth 1 | §7 |

A and B are one table with a column, not two products, because *"the line between 'new
strain of old thing' and 'new thing' doesn't exist in evolution"* (Imo, 2026-11-22). They
are separated below because they carry different measured constraints and different
refusals, not because the system scores them differently.

**Two things the system never emits.** There is no novelty score and no novel-or-not call:
score ordering tracks database proximity, measured repeatedly, so a search can say *is this
already known* and can never say *is this new*. And there is **no threat score**. Over
41,829 subjects reached by the residue search, host confidence maxes at 0.40 with a median
of 0.17 and nothing is human-associated (`results/full_residue_search/`), so a host label
describes **a database entry, not the water**, and the system prints it that way.

Provenance for the headline numbers: `results/subtract_known/REPORT.md`,
`results/subtract_recut/`, `results/holdout_cascade/REPORT.md`, `results/holdout_run/`,
`results/divergent_contigs/`, `results/reach/`, `results/path_settings/`,
`results/crossplant/`, `results/residue_census/`, `results/recurrence/`. The one-line index
to all of them, including everything killed, is [`results/README.md`](results/README.md).

---

## 1. What a positive output looks like

One row, and it is a sentence with two halves:

> **This contig is 66.5% identical at protein level over 620 aa to *X*, covering 46% of
> itself, at 452× even mapped depth, assembled from this water on this date — and its three
> genes span 24.3 identity points with the polymerase at 66.5% and the maturation protein at
> 53.5%.**

That is the withheld genome of `results/holdout_cascade/`, a real virus deposited after the
library that contains it, and it is the system's positive control. The distance half says
*how new*; the identity half says *like what*.

**Both halves are required, and that is Imo's correction of 2026-11-23: "what it IS is as
important as whether it's new — 'new AND like X' is the combined signal that says
trouble."** Novelty alone is a curiosity. A distance beside a named nearest relative and a
gene role is the sentence somebody can act on. So every row of every output class carries its
nearest-relative identity beside its distance, and where one half is missing the system says
which half is missing rather than reporting the other half alone. For class C that is the
current state: **distance without identity** (§7, §12.2).

---

## 2. The three settled facts the design rests on

### 2.1 It cannot be done from a single read

Measured, from several directions, and this is why no stage of the system emits a per-read
verdict.

- **No read-level rule gets near one false positive per million.** Peptide absence gives
  160,774 false flags per million at the strictest expressible threshold; read-level
  homology gives 8–10% sensitivity at 93,478 flags per million at bitscore ≥ 40 and 0%
  sensitivity at 1,489 flags at ≥ 80. ORF-likeness, k-mer composition, codon usage and GC
  have no novelty polarity at all, because a wastewater background is coding too
  (`results/readlevel_features/`, `results/peptide_absence_noise/`).
- **Thirteen axes, ~10.5 of them independent, and none discriminating**: the largest
  standardised difference between residue and identified reads is 0.40 once abundance is
  controlled (`results/residue_fields/`). The joint space beats the best single axis by
  **0.013 AUC** and shows no cluster count beating 1 in 5/5 seeds
  (`results/residue_structure/`).
- **At real distance a read recognises its own relative about one time in sixteen.** On the
  holdout at 30.3% nucleotide divergence, 5.25–7.00% of 145 nt reads clear E < 1e-5 at
  nucleotide level and 2.45–2.80% at protein, against 56.65–59.30% for referenced sequence,
  with a clean 0.00% shuffled null (`results/holdout_run/`).
- **And a per-read identity readout is biased towards *familiar*.** The same reads report
  82.56–83.33% identity where the truth is 69.72% — flattering by 12.8–13.6 points at 61–64%
  coverage, falling to +7.1 above 90% coverage. Coverage is the control that removes most of
  it and never all. Contig-level identity is essentially unbiased by comparison: 66.5%
  measured against 66.45% genome truth.

**Nucleotide beats protein at read level at every distance**, the reverse of the
contig-level result, because a 145 nt blastx alignment is ~48 residues. The protein tier's
contig-level advantage must not be carried down to reads.

So reads are used for three things only: subtraction, depth, and the measured boundary of
§7. Characterisation happens on contigs.

### 2.2 Novelty shows up as different genes at different distances

This is the discriminating shape the project found, and it is what class A ranks on.

- **Whole-contig identity is not a property of the contig.** A whole-contig best hit reports
  whichever gene won on bitscore. The median spread between a divergent contig's own genes
  is **17.4 points against a near-reference control arm's 2.35**, and it runs both ways:
  one contig reads 88.6% whole-contig with its polymerase at 39.8%, another 41.8%
  whole-contig with its polymerase at 58.6% (`results/divergent_contigs/`).
- By role, divergent against control: polymerase 66.8% / 97.4%, maturation 63.5% / 97.0%,
  coat 77.9% / 98.2%. **So the ranking is on the conserved gene — the polymerase — with the
  spread shown**, where there is one (103 of 149 contigs).
- **The holdout confirms it at a distance nobody chose**: three scorable genes spanning
  **24.3 points** inside one genome, maturation 53.5%, replicase 66.5%, coat 77.8% — the
  same role ordering, ground truth, in a genome withheld from every index.

- **The ordering is paired, and the pairing is what makes it readable.** On a pile an order
  of magnitude larger, compared role against role on contigs carrying both rather than as
  two medians over two sets, **maturation sits below polymerase on 122 of 131 contigs**,
  median difference −6.9 points, and below capsid on 28 of 31; capsid against polymerase is
  a coin toss. Counted by position instead, maturation is the most diverged role on 118 of
  203 divergent contigs and the least diverged on 7, with polymerase running the other way
  at 32 against 124 (`results/conserved_positions/`).
- **Both artefact objections were tested and neither explains it.** Reference abundance:
  distinct subject proteins matched per ORF is the same for the two roles, polymerase 0.604
  against maturation 0.598. Alignment length: inside matched length bands the gap holds, 7
  to 8 points at 150–300 aa and at 300–600 aa. Two roles carried by one polyprotein ORF are
  excluded, because those pairs are identical by construction.
- **The spread on the larger pile reads 12.8 points against a near-reference 2.45**, and it
  survives matching on the number of scorable genes — 5.6 points at two genes, 16.9 at
  three, 22.2 at four, against a control arm's 1.3, 2.65 and 3.4 (`results/divergent_known/`).
  The unmatched comparison is biased towards the divergent arm, whose contigs are longer and
  carry more genes, so the matched bands are the reading.

### 2.3 Novel strains of existing viruses are first class

Imo's call of 2026-11-22 supersedes every earlier framing: there is no separate strain
product. There is one continuum over everything that got a significant hit, ranked by
distance with coverage beside it, and **the tier a contig ended up in is a column on that
table**. §6 holds the constraints that are specific to the near end — chiefly that below
~3.8% divergence every assembler fuses lineages, so the system reports a fused lineage as
fused.

---

## 3. Why this is not one step

Imo, 2026-11-23: *"why isn't this one step? we just assemble everything, then see what reads
are left that don't match any assembly?"*

It is one **assembly**. It is not one **question**, and the difference is measured four
times over.

**(a) The subtrahend is identification, not assembly.** Subtracting every contig, including
the ones nothing identifies, removes **85.9%** of the library against **77.0%** for the
identified contigs alone, on the same block. The extra 8.9 points are not background: they
are the assembled-but-unnameable pile, which is the candidate material. Subtracting it
swallows class C's contig half whole. The run therefore records `"subtracted": false`
against that set rather than relying on the reader knowing the rule
(`results/holdout_cascade/`, and the defect that found it is in
[`hapax.py`](hapax.py)'s history).

**(b) "Reads that match no assembly" is not a stable quantity.** Ten whole-library
assemblies of the same library under different path settings move the residue by **8.79
points — 33.84% to 42.63%** — against a paired same-setting null of 0.01, while genome
recovery barely moves (20 of 25 genomes at ≥ 90% breadth for the default and four arms;
nothing beats the default). Turning off read rescue alone adds **+7.22 points, 8.24 of them
into tier 3**: reads that lose their contig rather than gaining a name
(`results/path_settings/`). So the leftover pile is substantially a property of the
assembler's path choices, and it ships with the flags that produced it or it is not a
number. *This is also why re-assembly of the residue is dead (§11): residue membership is
decided by identification, so residue-only assembly bakes in the previous run's paths and
cannot see a novel read that laxer settings absorb into a known contig.*

**(c) Assembly is not the bottleneck, and the leftover pile is not where the answer was.**
The withheld genome came out of 4.29 M reads as one 3,990 nt contig covering 96% of itself —
and was then **removed as known** at 66.5% protein identity, bitscore 812 against the tier's
threshold of 50, with 99.82% of identified contigs scoring at or below it. It never reached
the residue at all. **The thing we are looking for assembled perfectly and the pipeline
deleted it.** A design that only inspects leftovers would have missed the one positive
control whose distance nobody chose.

**(d) The most efficient index is the sample's own assembly, and it finds nothing novel.**
1.4 MB of the sample's own ribosomal contigs removes 38.85% — more than all 272 MB of SILVA
— and 15 MB of its identified non-rRNA contigs adds 30.33 marginal points where 16 GB of
Kraken2 adds 0.42. But the per-genome reanalysis of the spike shows it removes **0.0% of a
genuinely novel virus at every rung**, because a virus absent from the sample cannot be in
the sample's assembly (`results/spike_holdout/`). It is an **efficiency** instrument, never
a sensitivity one, and the 20.8% that had been read as the sample's own data catching novel
viruses was entirely an artefact of spiking resident material.

**So the shape is: assemble once, identify, and then run the distance continuum over what
was removed as well as over what was not.** The leftover reads are an output in their own
right (§7), reported as a measured boundary at median depth 1 — never as the candidate pile.

---

## 4. The pipeline

### 4.1 Resolve the per-library parameters, refusing rather than defaulting

Imo, 2026-11-23, on the mate-2 prefix: *"you mean trim primers?"* — close, and the
difference is why this is a stage of its own.

It is **not** adapter or primer trimming. Mate 2 of SRR23973193 carries ~14 nt of
**non-template 5' sequence from the kit itself**: measured by mapping untrimmed reads to
SILVA with the other end of the same read as the control, the R2 5' clip has **mode exactly
14 in 5/5 seeds on both RNA libraries**, 95% of anchored R2 reads clip ≥ 10 while 94.6% of
the *same* reads sit flush at their 3' end, and mate 1's mode is 0
(`results/nontemplate_prefix/`).

Three rules follow, and all three are refusals:

1. **It is a distribution, not a constant** (14 at 51%, 13 at 16%, with a tail), so there
   are **two operating points and a run uses one**: 14 is right for anything that reads
   sequence, because those bases are real template; **16 is the better Kraken2 setting**,
   because Kraken2's confidence is the *fraction* of a read's minimizers in the called clade,
   and two extra variable bases lower every read's confidence. Re-running at 14 moved the
   identified fraction **down** 1.03–1.11 points on this library and 1.60–1.84 on the 2022
   one. The trim ships as a curve (`trim_sensitivity_*.json`), and the run manifest states
   which point it used and what that costs.
2. **It is per library and the system refuses to guess.** `--trim-r2` is required in
   [`hapax.py`](hapax.py); the old default of 14 silently cut 14 real bases off every other
   library's mate 2. A refusing resolver, not a default — which is the general rule of §8.
3. **Mate 1 is never trimmed by a fixed amount.** It has ~2 bases of 5' excess but **no
   mode**, which is exactly what separates priming bias from a non-template prefix.

The same stage resolves the block stride (`--block`), which is not optional in practice: it
is a property of the library, and left at the Berlin default two of five seeds of a smaller
library read nothing while the directory still looked like a five-seed result.

### 4.2 Merge the mates before anything clever runs

A pair is one molecule of ~194 nt median insert, and **35.1% of this library's pairs merge
into a single read**. Measured on the spike ladder against mate 1 *of the same molecule*,
with the shuffled null searched at the merged length (`results/profile_reach/`):

| divergence | merged (~140 nt) | mate 1 (~105 nt) |
|---|---|---|
| 10% | 75.5 (72.0–78.9) | 59.6 (55.7–64.3) |
| 20% | 38.4 (33.2–43.5) | 10.5 (8.9–11.6) |
| 30% | 7.1 (5.1–9.2) | 0.0 (0.0–0.2) |
| real library reads | 16.0 (14.1–17.7) | 8.9 (7.2–10.6) |
| shuffled null, both lengths | 0.0 | 0.0 |

**So part of the 10–20% recognition wall is a read-length wall, not an index wall** — the one
wall this project has found that moves without new sequencing, new indexes or new money.

And it is a **test**, not a longer ruler: observed over length-predicted bitscore runs 0.987
on identical sequence and 0.745 on shuffles, so a merged read gains *less* bitscore than its
length predicts and least on noise. An E-value penalises a longer query, so only a true
alignment clears the penalty — which is why conversion of weak hits to significant runs
34.6% at 20% divergence and **0.0% of 387 shuffled weak hits in 5/5 seeds**
(`results/weak_hits/`). Two limits travel with it: no error-correction credit is claimed for
the overlap, and the ladder is one picornavirus with a single long polyprotein ORF, which
favours a translated search.

Merging runs **before** the cascade. Every read-level number is reported for the merged read
where one exists and per mate where it does not — never pooled, because the mates are not
interchangeable (§8).

### 4.3 Assemble the whole library

Default settings, and the flags travel with the residue number (§3b). Nothing beats the
default for recovery, and the levers that sound like levers are not: `--min-count` and `k`
were killed by holding them fixed across ten arms, and zero absorbed genomes, zero
inter-molecule joins and zero false joins appear in any arm against 53 zero-read panel
genomes (`results/path_settings/`).

**Measure chimaeras before claiming contiguity**, and read the chimaera number beside its
own false-positive rate: split-mapping fires on 46.3% of the divergent arm and 39.3% of a
near-reference control, so it describes how the assembler joins genes rather than divergence,
and nothing is dropped on it. **Coverage evenness is the field that does separate.** The
assembler's coverage header understates mapped depth about **fourfold** (37.7× header against
170.7× over 17,462,448 mates mapped back), so depth claims come from mapping reads back,
never from the header.

### 4.4 Identify — the cascade

Cheapest index first. Measured on 5 × 100,000 pairs; percentages are of read **pairs**, the
occurrence-weighted unit.

| tier | index GB | alone | marginal | cumulative |
|---|---|---|---|---|
| artefact / low complexity | 0.000 | 0.31 | 0.31 | 0.31 |
| SILVA rRNA (nucleotide) | 0.272 | 38.02 | 37.93 | 38.23 |
| Kraken2 viral | 0.667 | 5.33 | 5.21 | 43.44 |
| Kraken2 Standard-8 | 8.005 | 7.52 | 0.27 | 43.71 |
| Kraken2 PlusPF-8 | 8.005 | 6.78 | 0.00 | 43.71 |
| Kraken2 PlusPF-16 | 16.005 | 21.67 | 0.42 | 44.13 |
| **the sample's own rRNA contigs** | 0.001 | 38.85 | 1.00 | 45.13 |
| **the sample's own identified non-rRNA contigs** | 0.015 | 67.22 | 30.33 | 75.47 |
| DIAMOND RefSeq viral proteins | 0.239 | 1.43 | 1.43 | 76.90 |
| DIAMOND Swiss-Prot | 0.294 | 0.31 | 0.20 | 77.10 |

**Union 77.10% (77.02–77.28 across seeds). Residue 22.90%.**

Four design requirements come out of this table.

**Ordering is a requirement, not an accident.** rRNA is subtracted with a *nucleotide* index
before any protein index is consulted. Against a large protein index, ribosomal reads stop
looking novel faster than anything else, because of automatically annotated public entries
that are really translated rRNA — **63 entries out of 4.12 M** carry most of it, amplified by
a conservation multiplier of 388 per occurrence against 3.49 per distinct k-mer. A protein
index fails on rRNA in the *dangerous* direction: not "no match" but a spurious match, so the
read looks known.

**Thresholds dominate database size.** PlusPF-16 removes 42.82% at confidence 0, 21.67% at
0.1 and 3.52% at 0.2 — a factor of twelve inside one database against a factor of 1.01
between 8 GB and 16 GB. Every reported fraction carries its operating point; a "% unknown"
without a stated threshold is unreadable and the tool will not emit one. The same effect is
length-dependent in the direction nobody assumes: on contigs, confidence 0.1 collapses the
viral calls from 1,021 to 11 and **the survivors are the short ones**
(`results/public_health_view/`).

**Tier order is a property of the library, not of the tool.** On the DNA library from the
same water on the same morning — a virome, 0.01% ribosomal, 99.1% unidentified — 0.24 GB of
viral protein identifies **8.5×** what the 16 GB nucleotide index does, while on the RNA
library the 16 GB index wins sixfold. The order is chosen per sample by the calibration run
(§8).

**Both index kinds are needed because the length reverses the verdict.** On contigs,
protein-only identifies **26.30%** and nucleotide-only **5.73%** (34.78% of contigs, 49.13%
of bases overall); at read level nucleotide wins at every distance (§2.1). One index does not
serve both stages.

### 4.5 Score the distance on everything identified

The cascade already computes identity, alignment length, E-value, bitscore and coverage for
every contig it removes **and used to keep only a yes/no**. That discard is the defect, not
the removal: two sequences thirty identity points apart received the same single word
`identified`. `classify_contigs.py` therefore writes the removal evidence per contig, and
`removal_evidence.tsv` is a first-class output (§13). Nothing about *what* is subtracted
changes; what changes is that the subtraction is no longer silent about how far away the
thing it removed was.

### 4.6 Report the residue

§7, and with the two calibrations of §8 beside every field.

---

## 5. Output class A — divergent but related contigs

The largest and best-characterised of the three, and the one the goal of §0 actually runs
through.

### 5.1 What the system reports

Per contig, as a row on one continuum with no threshold on it:

| field | why it is there |
|---|---|
| best subject, percent identity, alignment length, E-value, bitscore | an identity without its E-value is unreadable, and identity is quantised at short alignment lengths |
| per-subject union coverage on the query | identity alone flatters towards *familiar*; coverage is the control that removes most of it |
| **per-gene identity with the role, and the spread across genes** | whole-contig identity is not a property of the contig: 17.4 points median spread on Berlin, 24.3 on the holdout (§2.2) |
| nucleotide **and** protein identity, both columns | three holdout candidates sit at 10.4–11.7% nucleotide divergence while 97.1–97.7% identical at protein level, so one column alone mis-ranks (`results/holdout_candidates/`) |
| which tier removed it, at which operating point | the same contig moves tier when a threshold moves |
| mapped depth, breadth, evenness, chimaera number | from reads mapped back, with the test's own false-positive rate beside it (§4.3) |
| the index the distance was measured against | a narrow index makes every distance an **upper bound** |

**Ranked by the conserved gene — the polymerase — with the per-gene spread shown.** Ranking
on whole-contig identity ranks on which gene happened to win.

**And the shortlist is ranked by how much of the contig its best hit explains, not by
identity.** This is the one mechanism measured to catch the positive control. Ordering the
removed pile by identity puts the withheld genome at the **93.46th percentile**, with 7,104
contigs ahead of it — identity ordering is not merely uncalibrated, it is *actively wrong*,
and `results/holdout_cascade/` also kills "flag the low-identity removals" directly: removed
contigs with a significant protein hit sit at a median 40% identity, so the holdout at 66.5%
is in the top decile of them. Its best hit covers **0.4624** of it, and identity-plus-coverage
gives a cell of **69 contigs, 0.91% of the removed pile**, with **not one** of the leak arm's
34 known contigs from the same water in it. Against core_nt, 319 GB, about a sixth of a
comparable panel collapses to strain level while the rest holds, so the cell shrinks by
roughly a sixth rather than dissolving.

**The candidate is not the top row of its own cell.** So the defensible output is **a short
ranked list with a count on it** — never a flag, never a score.

**Band labels are readings, not calls**, printed beside the number they came from and never
substituted for it: *reference or near-reference* ≥ 97%, *divergent strain of the named
taxon* ≥ 90%, *same genus or thereabouts* ≥ 70%, *family-level or looser — the subject is a
nearest neighbour, not an identification* ≥ 40%, and below that *below any band a name
survives*. A reading at 30–45% amino-acid identity is a **new species or genus borrowing a
named relative's name** — a bigger claim than a strain, not a smaller one.

**Every row ships with its control arm, drawn and scored by the same code in the same pass**:
near-reference material from the same library.

### 5.2 What this class refuses to claim

- **No host and no threat.** Host confidence over 41,829 subjects maxes at 0.40; a host
  colour labels a database record. And **what a record IS and what it INFECTS are two axes**
  — a phage record has a bacterial host and is not a bacterium; conflating them turned a
  *Bacillus subtilis* transporter into a bacteriophage in a live view. A cellular record
  reports host as *not applicable — it IS the organism*. The phage bucket is a **name test**,
  so its complement reads *virus, nothing says phage*, never *not a phage*; phage against
  bacterial is irreducible at ~50% identity over ~34 aa, and the reason is biological —
  prophages mean one protein can truthfully be both.
- **No public-health reading from architecture alone.** The seven near-complete divergent
  genomes on the Berlin panel are single-stranded RNA bacteriophage, so the reading is
  indirect and the system says so rather than implying a human target. geNomad calls 98.3%
  of contigs > 3 kb viral, so *"is it viral"* carries almost nothing.
- **No importance.** Six or seven of these are near-finished genomes, 1.8–3.9 kb at 171–467×
  with even coverage end to end, at family-level distance from anything published: the
  material worth real effort. But **importance will not come out of sequence.** Sequence
  gives distance and gene-level shape; occurrence across dates and plants gives the reason to
  care, and that axis is §12.6.

---

## 6. Output class B — novel strains of known viruses

Same table, same ranking, near end of the continuum, tier as a column. What is specific to
this end is a set of measured limits, and all of them are refusals.

- **Below about 3.8% divergence every assembler fuses lineages at any depth**; MEGAHIT
  separates from roughly 5% up (2,205 assemblies, 21 real genome pairs, 5 seeds, 0 false
  mosaics in 630 controls). **So below ~4% the system reports a fused lineage as fused**, and
  a reading at the near end ships with the divergence at which it was made
  (`results/pair_ladder/`).
- **The reads carry enough information; path selection is what is missing.** PenguiN's raw
  output holds clean multi-kb pieces of both lineages at 2% divergence, buried among ~1,000
  mostly mosaic contigs. That is a statement about tooling, not about the library.
- **Coverage pooling across lineages is dead**: it fabricates mosaics.
- **Pair scores cannot separate two hosts' swarms** — two HBV hosts share 709 exact 150-mers.
  Only population structure can.
- **A strain is 90–99% identity**, and the per-read identity bias of §2.1 lands squarely
  here: reads report ~83% where the truth is 69.72%, always towards *familiar*. Strain
  readings therefore come from contigs with coverage beside them, never from per-read
  identity.
- **Both identity columns ship** (§5.1), because a nucleotide band does not transfer to
  protein: 10.4–11.7% nucleotide divergence at 97.1–97.7% protein identity is the measured
  counter-example, and above 5% divergence the protein tier does essentially all the work.
- The lineage labeller's strand bug **voids every mosaic count published before its fix**
  (`results/labeller_audit/`): 0% false mosaics on forward-strand contigs, 65.8% on reverse.
  Any new mosaic number states that it was produced after the fix.

---

## 7. Output class C — the unplaceable residue

**This is a measured boundary, not a discovery pile**, and the system's language says so.
What lands here is mostly a function of read length and database reach: hypervariable
regions of textbook viruses land here too.

### 7.1 What the system reports

| measurement | value | control in the same pass |
|---|---|---|
| residue after the cascade | 22.90% of pairs | union 77.10% (77.02–77.28, 5 seeds) |
| of which tier 2 (assembled, unidentifiable) | ~8.4 points, 46,157 contigs, median depth 6 | 3.87% clear 500 nt and 2× — ~1,784 contigs are the characterisable set |
| of which tier 3 (no contig at all) | ~11.1 points, **median depth 1** | 86.6 for identified reads, ~3,000 for ribosomal |
| a significant named relative in the largest nucleotide database there is | 7.82% of all unplaced pairs | **10.00% of already-placed material** — a 2.2-point gap, not an order of magnitude |
| artefact signatures (hairpins, same-strand mates, poly(A/T), vector, PhiX) | under 1%, and the two fields that point anywhere point at the *identified* pile | shuffles, same pass |
| absent from the plant's own other 40 dates, at matched within-run depth | **21.5%** against **0.9%** for identified reads (± 0.6, 5 seeds) | own library 1.0000, every shuffle 0.0000 in 0 of 115 |
| both mates unwitnessed over exactly one | 2.03 residue, 0.29 identified, 0.00 rRNA, ~1,000 on shuffles | so absence is **missing molecules, not misread ones** |

**Depth decides, and the pile label almost nothing.** A union of 512 other libraries removes
87.5% of a held-out library and has **not saturated**; read down the depth column, not across
the pile row — at depth 100+ it takes everything including our own abundant residue, and at
depth 1 **identified reads survive at 0.7050 against the residue's 0.7900**. That 8.5-point
gap is the whole known-versus-unknown difference at matched depth, so a cross-plant union is
blind to *anything present once*, not to novelty. The two operating points disagree by a
factor of four and neither travels alone (`results/crossplant/`).

**At contig level the same instrument reverses and gets stronger**: tier 2 contigs cut into
90 nt windows separate from identified contigs in every length band, and **17.6% of tier 2
contigs have no exact 31-mer anywhere in four continents of wastewater against 9.0% of
identified ones**, with every shuffled contig untouched. The residue is **local**: its own
plant's other dates witness it 13.4× more than the median foreign plant.

**And the residue's two halves point in opposite directions, which is itself the finding.**
The 2,832 genuinely new contigs assembled from unplaced pairs carry **11.7% RdRp against
0.70%** for a length-matched control, dominated by Leviviridae replicase — they look viral.
The 200 highest-mass unidentifiable contigs are **less coding** than length-matched ones
(0.604 against 0.818) and lead every pile on Rfam — the deepest unnameable sequence is
non-coding and not a hidden virus. *"Not ribosomal" is scoped to a threshold, not
established*: those 200 read a median 0 matching bases against SILVA where true rRNA contigs
all clear the line, and that instrument cannot see divergent rRNA.

**One row per surviving pair, un-thresholded, each field with a shuffled null and a
known-positive comparison beside it** (`closest_relative`, `frame_evidence`, `gc`/`entropy3`,
`repetitive`, `artefact`, `chimaera` with its 26% false-positive rate on known background,
`motif` at the families' gathering thresholds, `mean_error_prob`, `recurrence`). **Nothing is
combined into a score, and `cannot discriminate here` is a legitimate row** — it is what
honest calibration of a 107 nt read produces about half the time.

The joint space ships as **a density with reference clouds and never as a clustering**
(§2.1). Cross-sample material is **context, never a subtrahend**: subtracting the other
libraries' contigs would make recurrence zero by construction.

### 7.2 What this class refuses to claim

- **No per-read verdict, ever**, and no aggregate called novelty.
- **No "matches nothing".** It always meant *matches nothing in an index that fits on this
  disk* — 86.8% of residue reads get *some* alignment, and the median one has E = 0.24 and is
  a chance match. The system says the former with an E-value and never the latter.
- **No depth claim.** A doubling of sequencing buys a second observation, not assembly depth:
  44.0% of fully unwitnessed reads gain a witness at full library and **0.3% reach depth ≥ 5**
  (41.0% and 0.3% measured on the unplaced pile itself). Below ~1 in 100,000 anything needing
  two observations of a molecule has nothing to work with.
- **No identity.** This is the honest gap of §1: for the pile sequence cannot name, we have
  distance and no identity. Structure is the candidate instrument and it is in flight in
  another thread (§12.2).

---

## 8. Calibration, controls and operating points

**Every score is referenced to this sample's own background.** The shuffled nulls, the
subtracted-read comparisons, the chimaera false-positive rate and the error profile are all
recomputed per sample. A threshold imported from another run is a threshold measuring the
other run.

This stopped being a precaution and became a measurement. Two RNA libraries from **the same
plant, the same group, the same naming series, eight months apart**: rRNA **38.0% → 67.1%**,
Kraken2 PlusPF-16 at confidence 0.1 **22.7% → 40.7%**, virus-identifiable **12.0% → 6.5%** as
ribosomal sequence crowds it out. **So every removal percentage in this repo is a statement
about one sample on one day** and is read as "on SRR23973193", not as a property of the
method. Drift is the biggest caveat on every figure.

**A null is not a bias statement, and a number needs both.** A null says what a field reads
on meaningless input; a bias statement says what it reads on *real* input for the wrong
reason, and only the second can be removed by a control. Every field states the **direction**
of its bias and the control that removes it, or the direction alone where no control exists.
Read the control first: three publishable-looking results on one day were wrong until their
bias was named, and the tell in all three was the control moving with the test.

**The bugs to watch are the ones that fake a result**, and the common shape is **a default
where a refusal belongs**:

- **the shuffle fallback** — a dinucleotide shuffle whose Eulerian walk dead-ends returned
  its input, so real sequence entered the null. It can only move the null *towards* the test,
  which makes "null scored zero" conservative and **"test and null are indistinguishable"
  the dangerous shape**: two kills were withdrawn. Dead ends happen in repetitive sequence,
  which is also what throws chance alignments, so 2.85% contamination supplied 100% of one
  null's hits. One shared module now (`shuffle_control.py`), raising rather than falling
  back, with a tripwire test on any new shuffle;
- a chimaera test that never ran and wrote `0.0`;
- a stage that silently skipped a missing index;
- an unexposed block stride that left two of five seeds reading nothing.

So: **controls, search-output parsing and provenance come from shared modules, imported,
never re-typed** (`shuffle_control.py`, `search_hits.py`, `provenance.py`), each raising where
a copy used to fall back. A stage whose reference data is absent writes
`status: "skipped"` with the reason and the run is marked incomplete; it never falls back to
a cheaper index and reports the number as if it were the real one, and it never guesses a
sample's reference.

**Sign flips, found twice, so every field states its direction measured on known material:**
"the unmatched half returns to the same subject" fires 8× more on *identified* contigs, and
the coverage-shape family has **familiarity polarity** — 0.78 AUC that means *already seen*
(a within-read depth step marks sequence shared between molecules, +0.577 with the number of
distinct contigs a read aligns to, not a junction).

**Five seeds give an honest range and a useless variance.** The min–max is an observed
spread, not an estimate of dispersion, and it reverses sign off a run's own seeds: the
inflation reading 1.44 and 1.61 over the whole library's 87 windows reads 0.65 and 0.44 off
five of them. **Any correction factor comes from the full series, never from the run's own
seeds.** The five cascade blocks are **contiguous blocks of one flow cell**: location is fine
(within 0.14 points of the whole library on every class) and dispersion is not
(**overdispersion 1.57× on identified, 1.43× on rRNA**, adjacent blocks correlated at +0.28),
so those Wilson intervals are that much too narrow and the design effect is applied and shown
beside the interval. The residue classes are nearly clean at 1.13–1.15×.

**The sample's own assembly is versioned.** It is the most powerful thing the tool subtracts
and it is produced at run time, so it gets an identifier, a digest and a manifest entry, and
every removal number names the assembly version that produced it.

**The mates are not interchangeable** and any pair-symmetric field inherits that: 107 against
93 nt is 76 scorable 31-mers against 62, and the mate direction asymmetry in linkage is
**entirely length** — matching the lengths turns 2.311/2.653 into 4.093/4.087 and nearly
doubles recovery. Report per mate, or show the length-matched control.

Plus the repo-wide rules: `--seed`, a run manifest beside every result (seed, git SHA,
**`script_sha256`** taken when the manifest is written, parameters, timings), release builds
for timings, five seeds minimum for a headline, spread reported rather than averaged away,
percentages of read pairs, figures generated from committed files and never hand-written.

---

## 9. Declared blind spots

Measured, not cautious, and the tool states them in its own output
([`blind_spots.json`](blind_spots.json) travels into every run directory).

**Four walls.** Nothing below ~1 in a million. No recognition past ~10–20% divergence, and
the band is stated at **protein** level because a nucleotide band does not transfer. Reads
find kin only to ~1 in 100,000. The unplaceable residue sits at **median depth 1**.

**Detection below about one in a million is not available at all.** At that abundance a virus
contributes ~8.7 read pairs in an 8.7 M-pair run. Detection is abundance-limited;
**characterisation at depth is where the 1–5% divergence band sits**, and that is the honest
headline.

**Dimensionality does not move either wall.** The recognition limit is a property of the
indexes, not of the decision rule, so adding axes cannot shift it; and at one in a million
there is nothing to aggregate.

**Recognition stops between 10% and 20% divergence, partly because of read length.** The
holdout spike — three real genomes simulated into this library's own measured insert
distribution and error rates, walked away from themselves, pushed through every index:

| divergence | spike pairs surviving tier 1 |
|---|---|
| 0% | 0.0% (0.0–0.0) |
| 5% | 0.6% (0.0–2.0) |
| 10% | 12.7% (6.4–19.2) |
| 20% | 82.2% (73.8–89.8) |
| 30% | 99.5% (98.2–100.0) |

Which tier does the recognising, as removal falls away:

| divergence | sample's contigs | Kraken2 viral | PlusPF-16 | viral proteins | Swiss-Prot |
|---|---|---|---|---|---|
| 0% | 65.5 | 100.0 | 64.0 | 97.5 | 97.5 |
| 10% | 21.8 | 28.1 | 0.4 | 82.9 | 83.0 |
| 20% | 0.2 | 0.0 | 0.0 | 17.4 | 17.6 |
| 30% | 0.0 | 0.0 | 0.0 | 0.5 | 0.5 |

**The nucleotide side dies first and the protein index is the last thing standing** — the
argument for keeping a protein tier in a cascade whose cheapest tiers are nucleotide. Read
the rungs as **ordinal**: uniform substitution is a crude model of real divergence. And on the
clean holdout the sample's own assembly removes 0.0% at every rung (§3d), so the per-genome
table is the one to read.

**Removal is a lower bound on "known".** An index cannot distinguish sequence with no
relative anywhere from sequence whose relative sits below the matching threshold.

**Subtraction is not free.** SURPI+ measured pre-processing removing 0.43% of true viral
reads, and at a one-in-a-million target a small false-removal rate matters more than the
headline removal.

**The cost curve has not flattened, so the far end is not quoted.** Each tenfold step of
protein index still buys ~20 points at the largest size reachable here; the nucleotide
reference ladder's last doubling still buys 15–26 points, but that instrument is weak and we
say so — reads concentrate over few taxa (top taxon 49.0% of classified reads, top ten
86.4%), so a step is largely one dominant taxon being in or out. Measured points only.

**Background is a baseline too.** Anything rising well above the panel's own level gets
reported, not removed — and the growth detector as built **points the wrong way**: the 100×
tail is 14× commoner in sequence most of the panel also carries, so it finds a common
organism blooming locally. The pooled contig growth figure is an artefact and must not be
quoted (93.2% of tier 2 contigs have zero panel depth under them).

**Provenance gaps, stated rather than tidied.** No cascade manifest carries `script_sha256`
and `subtract_known/`'s five seeds sit under three SHAs, two `-dirty`; the 77.10% headline is
read from `subtract_recut/`, whose five seeds share one clean SHA. `linkage_census/`
manifests name SILVA though the census never read it.

---

## 10. What is measured — the provenance index

Load-bearing numbers and where each comes from. Nothing else in this document is a claim.

| claim | where |
|---|---|
| Identification removes 77.10% (77.02–77.28, 5 seeds); residue 22.90% | `results/subtract_known/`, `results/subtract_recut/` |
| The subtrahend is identification, not assembly — 85.9% against 77.0% | `results/holdout_cascade/` |
| The sample's own assembly is the most efficient index (+30.33 marginal) and removes 0.0% of a genuinely novel virus | `results/subtract_known/`, `results/spike_holdout/` |
| Residue size spans 8.79 points on assembler path flags; recovery does not move | `results/path_settings/` |
| The holdout assembles at 96% breadth and is removed as known at 66.5% / bitscore 812 | `results/holdout_cascade/` |
| Identity ordering is actively wrong (93.46th percentile); coverage gives a 69-contig cell | `results/holdout_cascade/` |
| Whole-contig identity is not a property of the contig — 17.4 points median gene spread, 24.3 on the holdout | `results/divergent_contigs/`, `results/holdout_cascade/` |
| The role ordering is paired and holds on 122 of 131 contigs; neither reference abundance nor alignment length explains it | `results/conserved_positions/` |
| "Identified" spans eighty points: 4.27% of 10,172 contigs at 90%+ identity, 90.2% below 70% | `results/divergent_known/` |
| The gene spread is 12.8 against 2.45 and survives matching on gene count | `results/divergent_known/` |
| Existence of a protein alignment is the calibrated readout: 0 of 11,361 dinucleotide shuffles, 89.53% of real contigs | `results/divergent_known/` |
| The coverage-explained gap has no shuffled null; its rate on identified background is median 0.38, p99 46.2 | `results/divergent_known/` |
| Nucleotide beats protein at read level at every distance; the reverse holds on contigs | `results/holdout_run/`, `results/subtract_known/` |
| Per-read identity flatters by 12.8–13.6 points; coverage is the control; contig-level identity is unbiased | `results/holdout_run/`, `results/holdout_cascade/` |
| Merging the mates takes 20%-divergence recognition from 10.5% to 38.4%; the null stays at 0.0% | `results/profile_reach/`, `results/weak_hits/` |
| Nothing beats the production search — 11 arms, identical to the decimal; the constraint is the significance statistics | `results/reach/` |
| Structure reads fold, not distance (32.6% against 31.2% on near-reference) | `results/reach/` |
| No single-read axis discriminates (largest standardised difference 0.40); joint space +0.013 AUC; no clusters | `results/residue_fields/`, `results/residue_structure/` |
| Mate-2 prefix mode exactly 14 in 5/5 seeds, two operating points, 1.7 points apart | `results/nontemplate_prefix/` |
| Depth decides, not the pile label: 8.5 points at depth 1 after a 512-library union | `results/crossplant/` |
| 21.5% against 0.9% absent from the plant's other dates at matched depth | `results/recurrence/` |
| Residue is thin, recurring, not artefact: median depth 1, artefact under 1% | `results/residue_census/` |
| Host confidence maxes at 0.40 over 41,829 subjects | `results/full_residue_search/` |
| Below ~3.8% divergence every assembler fuses lineages; MEGAHIT separates from ~5% | `results/pair_ladder/` |
| Blocks are contiguous: overdispersion 1.57×/1.43×, design effect applied | `results/block_position_bias/` |
| Which committed claims rest on a null that did not shuffle, and which way it points | `results/fail_loudly_audit/SHUFFLES.md` |

---

## 11. What is dead — do not re-propose

Each of these was measured and failed, and the record is in
[`results/README.md`](results/README.md). Re-proposing one is the specific failure mode this
section exists to prevent.

1. **Any novelty score, any novel-or-not threshold, any flag** — score ordering tracks
   database proximity. And **any threat score**: host confidence cannot carry one.
2. **The weak-hit headroom.** "92.3% of shuffled reads hit" was taken without an E-value
   floor; with one it is 0.00%. Bigger databases buy hits, not identifications.
3. **Composition in any form** — k-mer composition, codon usage, GC, ORF-likeness: no novelty
   polarity, because a wastewater background is coding too.
4. **Peptide k-mer absence** as a detection rule: 160,774 false flags per million.
5. **Read-level Pfam / profile reach.** Pfam-A is weaker than translated search at every rung,
   and the "DIAMOND reads zero past 20%" belief was a cascade removal figure misread as a
   search's reach.
6. **Coverage-shape fields as novelty axes** — 0.78 AUC with *familiarity* polarity.
7. **Read-level chimaera and partial-match tests.** Once a 107 nt read matches ~74% of
   itself the leftover is ~10 aa, under any search's floor; the chance band is identical
   between residue and identified reads. This question belongs after assembly.
8. **Clustering the joint space** — no cluster count beats 1 in 5/5 seeds.
9. **Read-level nucleotide BLAST hit rates without an E-value**, and identity quoted without
   one at any level.
10. **Coverage pooling across lineages** — fabricates mosaics. And **`k` / `--min-count` as
    the assembly levers that matter.**
11. **"Absence in wastewater means absence"**, and *"93% unclassified"* as a statement about
    the water: it is one index's reach. Most of what we remove as known, we remove by matching
    sequence nobody has named — one index names 6.88%, all name-returning tiers 7.49%, and the
    69.54% between is recognised-but-unnameable (`results/public_health_view/`, of which
    37.90% is a SILVA match and most of the rest a contig from this sample).
12. **Flagging the low-identity removals**, and **ranking the removed pile by identity** — the
    holdout sits in the top decile of significant-hit removals and at the 93.46th percentile
    overall.
13. **Re-assembling only the residue** — residue membership is decided by identification, so
    it bakes in the previous run's path choices (Imo's correction; vary settings on the whole
    library).
14. **"Most distant relative recovered" as a statistic** — it reads 20.7–22.2% on every pile
    including a companion at 99.8% identity, because at a fixed E-value the furthest surviving
    alignment is set by database size and query length.
15. **The "all subjects explain more than the best single one" gap as a second axis** — it is
    protein-only and is coverage read twice (−0.8 at nucleotide). The field itself is untested
    rather than killed (§12.7) and nothing rests on it.
16. **Seeding as the binding constraint.** DIAMOND at three sensitivities, BLASTP and full
    Smith-Waterman over 1,210,551 subjects return the same nearest relative on all 16 contigs.
    **The constraint is the significance statistics.** Codon awareness saturates above ~30%
    protein divergence: right tool, wrong distance (11 of 22 divergent genes have no
    measurable synonymous divergence left).

---

## 12. What is open — with the experiment that would settle it

### 12.1 Paired search mode — the one live candidate (untested)

Imo, 2026-11-23: *"even when mates aren't overlapping, many tools have a paired mode that
lets them work together to find homology. have we tried that?"*

**Honest answer: not as a search mode.** What exists today:

- The Kraken2 tiers **do** run paired (`-P` in `experiments/subtract_known.py`), so a pair's
  minimizers are already pooled in the k-mer tiers — **and that has never been ablated
  against per-mate classification**, so we do not know what it is worth.
- Every **homology** tier searches mate by mate. No tier anywhere combines two
  non-overlapping mates into one significance statement.
- §4.2's merge covers only the 35.1% of pairs that **overlap**. The other ~65% carry two
  reads of one molecule and are searched as if they were two molecules.
- Mate linkage has been measured at the **geometry** level — a read's own partner is in its
  10 nearest neighbours 23× more often than chance with no shared 31-mer, null on chance —
  but that **did not move the 1-in-100,000 kin wall** (`results/mate_friends/`). It bounds
  what pairing can buy as a *locator*; it says nothing about what it buys as a *significance
  statistic*.
- Prior art makes the gap sharper, not softer: Viromancer/Lyra-Virus **has pairs and never
  uses them** (`results/readlevel_features/`).

**The experiment.** Arms, all in one pass on the holdout ladder and the Berlin spike, scored
at the same operating point:

1. best single mate (today's behaviour);
2. merged read where the pair overlaps (§4.2, as the known-positive comparator);
3. **joint significance over the two mates**: both mates against the same subject, requiring
   insert-compatible coordinates and orientation, with the two E-values combined explicitly
   (Fisher) rather than by taking the better one;
4. Kraken2 `-P` against per-mate classification — the ablation of what we already ship.

**Controls, and arm 3 is worthless without the first of them:**

- **swapped partners** — mates paired with random partners from the same library. Any joint
  statistic gains on *any* pair; if the swap control gains as much, the statistic is reading
  query length and database size, not linkage. This is the control that killed the "most
  distant relative recovered" statistic (§11.14) and it is the same shape here;
- dinucleotide shuffles at matched lengths, through `shuffle_control.py` (raising, not
  falling back — §8);
- reads the cascade already identifies, as the known-positive arm.

**Metric**: recognition at a fixed false-pair rate per million, per divergence rung, **and
whether the nearest relative changes** — not merely whether a hit appears. An arm that finds
more hits without improving the identity or the named subject has bought nothing (§5.1).

**The rule it must respect**: paired mode may improve **identification and distance**. It may
never become a novelty score. No pair rule reaches one false positive per million
(`results/mate_width/`), and that kill stands.

**Caveat on the instrument**: the strongest holdout libraries are single-end, so §4.2's
mate merge is untestable there and this experiment needs the Berlin spike ladder plus a paired
holdout library.

### 12.2 Identity for the pile sequence cannot name

Class C currently reports distance with no identity (§1, §7.2). Structure is the candidate:
against all 178,660 chains of every solved virus structure, the divergent placements are the
right proteins and the holdout's polymerase lands at E = 1.1e-22 with a shuffled null of 0 of
80 — but **the field reads fold and not distance** (32.6% median 3Di identity on the divergent
pile against 31.2% on near-reference contigs 95–99% identical to published sequence), and on
the reach panel it placed 8 of 24 divergent genes against sequence's 22, the 8 a subset. That
is an identity instrument with no distance half, which under §1 is a **row in the report, not
a failed experiment** — and the place to point it is exactly where sequence gives no identity
at all. A databank-wide structure library is being measured in a separate thread; until it
lands, class C says *no identity*.

### 12.3 The full protein index — this is disk, not science

Requirement (iii) (evidence that no-match sequence really sits below ~80% identity) needs a
translated search against full UniRef50 or nr: **~160 GB**, against ~15 GB of writable disk
after the nucleotide indexes. The streaming recipe already does the hard part (41,243,524
records, verified against UniProt's own count). One machine closes the requirement and puts
the measured far end on the cost curve instead of an extrapolation. **Still the single
highest-value piece of infrastructure in this spec.**

### 12.4 The Berlin arm of the path-settings sweep

`results/path_settings/` ran on SRR23421544 because the Berlin arm is ~8 h (Imo's call). So
how much *the discriminating dataset's* residue moves under path settings is unmeasured,
and the 8.79-point swing is quoted as a property of that library.

### 12.5 One half of the gene-spread result is untestable

`results/divergent_known/` and `results/conserved_positions/` are on `main` and §2.2 cites
them. What is open is not the merge but the claim: **the within-genome spread as a novelty
signal appears unclaimed in the literature** (the closest precedent, Friedman et al. 2009, has maturation most divergent
*between* lineages with replicase closest — our ordering; capsid-first is *within* a
lineage), so it is the most interesting unclaimed thing here.

Conservation **inside** a protein is **cannot tell, not killed**: there are no shuffled
alignments to compare against, because all 148 shuffled ORFs give depth-1 alignments and
`hhmake` refuses them. That applies to any profile-based novelty claim, not just this panel.
Domain *presence* does have a clean null on `main` — **0 of 932 shuffled ORFs** and 0 of
11,361 shuffled contigs produce a significant protein alignment under the fixed walk
(`results/fail_loudly_audit/SHUFFLES.md`) — so what is missing is the comparison *conditional
on* an alignment, not the existence test.

### 12.6 Occurrence as the importance axis

**Importance will not come out of sequence.** Sequence gives distance and gene-level shape;
occurrence across dates and plants gives the reason to care. What is measured: 21.5% of
once-seen residue reads appear on none of the plant's other 40 dates against 0.9% of
identified reads at the same depth, and transience is ruled out for whatever recurs at all.
What is **not** measured: **"appears, then climbs across dates"** — the shape that would
actually matter — and the growth detector as built points the wrong way (§9). Report
occupancy over runs, waters **and** dates, because PRJNA948850's 215 RNA runs are **85
distinct waters over 41 dates** (116 paired `RNA_WW_*` plus 99 single-end `xh_*` runs are the
same libraries through a hybrid-capture panel).

### 12.7 The coverage-explained gap as a field

The any-subject-union against per-subject coverage gap (+36.6 points on the holdout) fires at
the 94.2nd percentile of real identified background — a few hundred contigs per library — and
correct nulls on two libraries read 0 of 22,110 and 0 of 11,361, which **bounds chance at ~1
in 11,000–22,000 rather than establishing zero**. **Untested, not killed, and nothing rests on
it.** As a *second axis* it is dead (§11.15).

---

## 13. Interface

### 13.1 The command line is the contract

```
hapax run SAMPLE.R1.fq.gz SAMPLE.R2.fq.gz --seed N --trim-r2 K --out RUNDIR
```

`RUNDIR` holds the manifest, the tier table, the per-read rows as JSONL, the per-contig rows,
**`removal_evidence.tsv`** — one row per route per removed contig, carrying the identity,
alignment length, E-value, bitscore and query coverage that justified the removal (§4.5,
§5.1) — the calibration block, the assembly version and `blind_spots.json`. Everything
downstream, including the web view, reads those committed files and nothing else.

`removal_evidence.tsv` is emitted in input order with no score and no ranking on it, and
`route` is a column because a protein row carries an identity and a coverage, a Kraken2 row a
taxid and a confidence, and a ribosomal row matching bases: three kinds of number that must
not be pooled.

**Implemented in [`hapax.py`](hapax.py)**, deliberately as an orchestrator: every number in
the run directory comes from the `experiments/` script that already measures it, with that
script's own seed and manifest, because a parallel implementation would be a second
instrument and this project does not compare two. Stages: `cascade`, `contigs`, `reads`,
`linkage`, `census`, `calibration`, `intervals`, `layers`, `strains`, `divergent` — where
`divergent` reads the continuum of §5 over the removed pile.

**The driver refuses rather than defaults**, which is §8's rule as code: `--trim-r2` is
required, `--block` arithmetic is checked against the library's counted reads before the run
starts, a per-read field whose reference or tool is absent stops the run unless it is dropped
*on purpose* with `--report-without`, and a stage whose inputs are missing writes
`status: "skipped"` with a reason and marks the run incomplete. It separates *searched, zero
hits* from *not searched*, and it never guesses a sample's reference.

### 13.2 Figures and the web view are generators, not pages

`experiments/make_infographic.py` and `experiments/make_web.py` read a `RUNDIR` and write
their output; a number typed into a figure is a number without provenance, and this project
has withdrawn claims whose prose drifted from its outputs. Missing runs are left out rather
than filled in.

Panels: the breakdown as one stacked bar with intervals and operating points; the cascade
table with the threshold sweep live on a slider, so the reader sees the factor of twelve
inside one database; the cost curve with measured points only and the unflattened far end
shown as **absent** rather than extrapolated; the divergence ladder labelled ordinal; **the
class A continuum, ranked by polymerase identity with the per-gene spread drawn and the
coverage-explained shortlist shown as a cell with a count**; the residue table with every
field's null and known-positive beside it; the joint space as a density with its two
reference clouds; blind spots as a fixed panel, not a footnote.

**Three rules, enforced in the generator and not in review**: no figure without its operating
point, no pile without its interval, no score without the calibration run behind it. A panel
whose operating point is missing renders the words `operating point absent`; a fraction with
no interval is drawn as a bare tick; `cannot discriminate here` renders as a value.

Imo's reading preferences are part of the contract: summary first then detail, pictures over
words, rank and hide the tail, names not accessions, and **phage, virus and bacteria shown
separately on the two axes of §5.2**. Every figure ships its control row, drawn by the same
code.

### 13.3 Every sensitivity control the reader can honestly move, moves

Imo, 2026-11-23: *"all sensitivity controls within reason adjustable at output time. we want
this page to feel like a playground, but uncluttered and with sensible defaults."*

This is cheap because of a decision already taken: **the report is un-thresholded by
construction**. Every per-read and per-contig row already carries the continuous numbers —
identity, alignment length, E-value, bitscore, coverage, depth, occupancy — and the cascade
already writes the per-tier per-read outcome. So moving a control is **re-filtering and
re-sorting what was computed**, not recomputing it, and the page is a reader over committed
files either way (§13.1).

**Four rules, enforced in the generator.** They are what stops a playground from
manufacturing the numbers this project has had to withdraw.

1. **A control moves its null, its known-positive control arm and its census with it.** A
   threshold quoted without its own null is the exact shape behind the withdrawn readings of
   §8, and the withdrawn "92.3% of shuffled reads hit" was a figure taken with no E-value
   floor at all. So the shuffled rows and the identified-material rows are filtered by the
   same predicate in the same pass, and the panel redraws all three or none. The sliders
   are **bounded to the range where the null and the control exist**; past that edge the
   panel renders `null not measured at this setting` rather than a clean-looking number.
2. **Discrete where the measurement is discrete.** Kraken2 confidence was *run* at
   0 / 0.1 / 0.2, so that control **snaps to the measured points and never interpolates** —
   a factor of twelve lives between them and a smooth slider would invent the middle. Same
   for the index-size axis: measured points only, and the unflattened far end stays absent
   (§9).
3. **Every reading carries its setting, and the shipped default is marked.** Defaults are
   the operating points this spec measured; a moved control is shown as moved, every export
   and copied figure carries the setting string, and the page state is in the URL so a
   reading is citable. One click returns to the shipped defaults.
4. **A control may change what is shown, never what was computed**, and where moving it
   re-labels material the label re-renders with it. Thresholds move contigs between tiers, so
   a contig reads *tier 1 at your setting, tier 2 as shipped* rather than silently changing
   pile (§4.4).

**Adjustable, because the numbers already exist at every setting:** the significance floor
(E-value or bitscore) on any search panel; identity and coverage bands, including the two
axes that define class A's shortlist cell, with its count and its control count recomputing
together; the minimum contig length and depth for the characterisable set (shipped at 500 nt
and 2×); ranking by the conserved gene against whole-contig best hit, with the spread always
drawn; which cascade tiers are in the union, with the **marginals recomputed rather than
reused**, because marginal contribution depends on order; merged-read against per-mate view;
the recurrence depth stratum and occupancy cut; and the residue table's filters, with
`cannot discriminate here` always a visible value and never filtered away by default.

**Not adjustable, and the page says why rather than leaving it off silently:**

| control | why it cannot move at output time |
|---|---|
| database size or choice beyond what was run | a measurement at a size we did not run cannot be synthesised; the cost curve shows measured points and marks the far end absent |
| the mate-2 trim | it changes the reads, so it changes every number downstream. It ships as **two measured runs** (14 and 16, 1.7 points apart) with the curve, not a slider |
| assembly path settings | they move the residue by 8.79 points and need re-assembly and re-identification; they ship as arms with their flags (§3b, §4.3) |
| search sensitivity and seeding | DIAMOND at three sensitivities, BLASTP and full Smith-Waterman return the same nearest relative to the decimal, so there is nothing to tune; the constraint is the significance statistics (§11.16) |
| the interval and the design effect | not a control. The Wilson interval is widened by the measured overdispersion (1.57× / 1.43×) and cannot be switched off (§8) |
| the two cross-plant operating points | they disagree by a factor of four and **neither travels alone**, so both are shown rather than one being selectable (§7.1) |
| anything that would produce a novelty call | there is no threshold to expose, because there is no score (§0, §11.1) |

**Uncluttered means a handful visible and the rest behind one disclosure**, in the same
summary-then-detail shape as the panels: the controls that change a headline sit with that
headline, and the rest live in a Detail list. A control the reader never moves should be
invisible, and a control that cannot move honestly should be visible as a stated constraint.

---

## 14. What this is not

It is not a novel-virus detector and it emits no novel-or-not call. It is an instrument that
says, for one sample: how much is known and by what, how much assembles, what sits at real
distance from anything published and **what that distant thing is most like**, how much
cannot be characterised at all, and what each answer cost. The deliverable is the three
output classes on one distance continuum, the cost-versus-completeness curve, and a
multidimensional measurement of signal against noise in which every axis carries its own null
and its own per-sample calibration.

The residue is reported as residue. The known is reported with the distance that justified
removing it. And nothing in the output says *threat*.

---

## Appendix — where the old section numbers went

This revision renumbers the spec around the three output classes. Committed scripts and
reports cite the previous numbering, and **they are deliberately not edited**: a comment-only
change would put their committed `script_sha256` out of step with the manifests they shipped,
which is the exact condition that hash exists to rule out. They get corrected the next time
those scripts run. Until then, read citations through this table.

| cited as | now |
|---|---|
| §1, §1.1 (architecture; why a measurement, not a call) | §2, §3 |
| §2 (tier 1, the cascade table) | §4.4 |
| §2.1(i) (divergence of anything not a close homolog) | §4.5, §5.1 |
| §2.1(ii) (part aligns, rest does not) | §11.7 — withdrawn at read level, lives after assembly |
| §2.1(iii) (evidence that no-match sits below ~80%) | §12.3 |
| §3 (tier 2, assembled but unidentifiable) | §7 |
| §4, §4.1 (tier 3; the per-read evidence report) | §7, §7.1 |
| §4.2 (the axes of the measurement; item 5 withdrawn) | §7.1, §2.1, §11.5 |
| §4.3 (merge the mates first) | §4.2 |
| §5, §5.1 (the distance continuum; the four residues) | §5, §6, §7 |
| §6 (calibration and versioning) | §8 |
| §7 (declared blind spots) | §9 |
| §8.1 / §8.2 / §8.3 (CLI, figures, web view) | §13.1 / §13.2 / §13.2 |
| §9 (what this is not) | §14 |
