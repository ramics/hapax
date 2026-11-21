# hapax — production tool spec

One sample in, one breakdown out. What is known, what assembled, what did not, what codes
and what does not — each pile quantified, with its own confidence interval and its own
declared blind spot.

**The target is not a novel-or-not call. It is a multidimensional measurement of signal
against noise.** That is a deliberate change of goal, and §1.1 gives the measured reason
for it.

Every number in this spec is measured on Berlin wastewater SRR23973193, five seeds,
and is reproducible from a committed script plus a seed. Nothing here is a projection.
Provenance: `results/subtract_recut/`, `results/spike_holdout/`, `results/nucleotide_ladder/`,
`results/blast_check/`, `results/subtract_known/REPORT.md`.

---

## 1. The architecture

Three tiers, in Imo's words:

1. **Subtract everything identifiable at all.**
2. **Check what assembled but could not be identified** for novelty.
3. **Characterise the residual reads** that neither assembled nor identified.

One rule governs the whole thing: **searching only ever removes background. It is never used
to score novelty.** An index can answer *is this already known*; it cannot answer *is this
new*, because score ordering tracks database proximity, not novelty (measured, PR #8/#11:
only 11.6% of spike reads hit at all, and background outranked spike at every threshold).
The subtrahend is therefore defined by **identification**, not by assembly — otherwise tier 2
is swallowed whole by tier 1.

The output is a **breakdown plus a cost-versus-completeness curve**. It is not a discovery
claim. The residue is the honest remainder, not a finding.

### 1.1 Why the target is a measurement and not a call

Every feature this project has killed, it killed **as a threshold on one axis**: peptide
absence, read-level translated homology, ORF-likeness, k-mer composition, codon usage, GC,
and every paired-end rule. The joint distribution of those axes has never been measured.
So the kills do not say there is no signal in the residue; they say **no single axis carries
a call**, which is a weaker and much more useful statement.

The machinery was already built for this. The per-read report is deliberately
un-thresholded and already ships, per field, a shuffled null and a comparison against reads
the cascade did subtract. Changing the target from a call to a measurement does not change
the machinery; it changes what the kills mean. They stop being failures and become the
calibration.

**And the noise here is mostly not instrument noise.** It is index incompleteness — which is
what made 86% of this library look novel until the sample's own assembly was indexed —
annotation artefacts, where a handful of automatically annotated public entries make
ribosomal reads look known through a conservation multiplier of 388, and library artefacts
of the kind a match-nothing pile collects by construction. Those have different signatures
from biology. That is the reason a multidimensional separation is plausible where a
threshold was not, and it is the reason every noise source below is measured per sample
rather than assumed.

---

## 2. Tier 1 — subtract the identifiable

Cheapest index first. Measured on 5 × 100,000 pairs; percentages are of read pairs, which is
the occurrence-weighted unit (per distinct k-mer the same quantity differs by two orders —
mean reach 3.49 vs 388 on curated SSU — and a read sees the occurrence-weighted one).

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

Three things the production tool must inherit from this table.

**The sample's own assembly is the most efficient index by orders of magnitude.** 1.4 MB of
ribosomal contigs assembled from the sample itself removes 38.85% alone — slightly more than
all 272 MB of SILVA. 15 MB of its identified non-rRNA contigs adds 30.33 points marginally,
where 16 GB of Kraken2 adds 0.42. *Known to a database* and *known here* are different
questions and the cheap one is the second. The tool therefore assembles the sample and
subtracts against itself as a first-class tier, not as an afterthought.

**Thresholds dominate database size.** PlusPF-16 removes 42.82% at confidence 0, 21.67% at
0.1, 3.52% at 0.2 — a factor of twelve inside one database, against a factor of 1.01 between
8 GB and 16 GB at fixed threshold. **Every reported fraction carries its operating point.**
A "% unknown" without a stated threshold is unreadable and the tool will not emit one.

**Ordering is a design requirement, not an accident.** rRNA must be subtracted with a
*nucleotide* index before any protein index is consulted. Against a large protein index,
ribosomal reads stop looking novel faster than anything else, caused by automatically
annotated public entries that are really translated rRNA — 63 entries out of 4.12 M carry
most of it. A protein index fails on rRNA in the *dangerous* direction: not "no match" but a
spurious match, so the read looks known. Without this ordering, a third of the library's
score measures public annotation quality rather than what is in the water.

### 2.1 The three requirements on the subtraction index

Imo's, and each is a reported field rather than a threshold:

**(i) Divergence of anything that is not a close homolog.** Closest relative is reported as a
continuous number — bitscore, percent identity, alignment length — never thresholded into a
call.

**(ii) Reads where part aligns well and the rest does not.** This was written as the sharpest
open question in the design, on the grounds that for the residue it is the *normal case* — a
ribosomal read aligning end to end at 100% identity against a residue read aligning over
about a third of itself at 89%.

**That premise is withdrawn, and it was withdrawn because nobody had read the E-value.**
BLAST reports to E=10 by default, so the "third of the read at 89%" was a description of
chance matches: the median E-value in that pile is 0.24, and a 107 nt nucleotide query finds
40 bases at 88% by chance in 318.8 GB. Demanding significance inverts the shape
(`results/subtract_known/REPORT.md` §3.4):

| readout, at E < 1e-5 | residue reads against core_nt |
|---|---|
| reads with a significant hit | 73 of 500 (14.6%) |
| median fraction of the read aligned, among those | **0.96** |
| median identity, among those | 90.1% |
| significant *and* partial (< 0.95 of the read) | 28 of 500 (5.6%), aligning 87 of 107 bases |

**A real relative aligns end to end. Partial alignment is the signature of a chance match
rather than of divergence** — so at read level the shape this requirement was written about
has almost no population, and what little it has aligns 87 of 107 bases rather than half.

A second instrument reaches the same place with a control
(`results/partial_match_subjects/REPORT.md`): among reads whose best hit covers part of them,
the chance band is **identical between the residue and reads the cascade already identified**
— median E 0.82 against 0.81, identity 87.8% against 87.2%, 38 matched bases against 39. A
separator that reads the same on the residue and on known material is describing the
database. 218 of 231 such hits cannot even be reproduced by a short-read aligner against the
reference BLAST named.

So the two competing explanations are no longer symmetric:

- a **chimaera** — amplicon chimaera rates exceed 70% for the rarest taxa, and they form
  reproducibly across replicates so replicate agreement cannot filter them. But the
  split-mapping test fires on **0.04% of unidentifiable contigs against 0.15% of identified
  ones**, indistinguishable from its own 26% false-positive rate on known background, and at
  read level the question cannot be asked at all: once a 107 nt read matches ~74% of itself
  the leftover is ~10 aa, below any search's floor;
- a **genuinely divergent molecule** where only a conserved domain is still recognisable —
  which is what the significant pile looks like, holding its identity into reference sequence
  past the reported match.

**Consequences for the tool.** This question moves from the read stage to **after assembly**,
where the separator is continuous and clean: unidentifiable contigs align over 31.5% of
themselves at 35.0% identity against 56.7% at 40.6% for identified ones, with no overlap. At
read level the tool reports coverage and identity **with the E-value beside them** and makes
no chimaera claim. The field name stays; the threshold never existed and now the premise
does not either.

**(iii) Evidence that no-match sequence really sits below ~80% identity.** This is the job of
a translated search against nr, and partly discharged already by ground truth against NCBI
core_nt (133 M sequences, 1.1 Tbase — two orders of magnitude beyond any local index):

The per-class table lives in `results/blast_check/TABLES.md` and is generated from the run
rather than typed here. **Read it at the strong-hit column, not the any-hit column**: an
"any hit" rate against 318.8 GB counts chance matches, which is the error §(ii) above
records. At 95% identity over 80% of the read the classes separate — rRNA 90.0%, public
nucleotide 94.4%, identified by the sample's own assembly 25.2%, residue 9.3% and 5.6%,
tier-2 contigs 3.9% — where at "any hit" they all read 85–100% and say nothing.

This still corrected our own language in the direction that matters: **most residue reads do
hit something, and "matches nothing" always meant "matches nothing in an index that fits on
this disk".** The tool says the former and never the latter. But it says it with an E-value,
because without one the sentence is true and empty.

Two things sharpen (iii) rather than soften it. **14.6% of residue reads have a significant
relative in core_nt and it is ~90% identical over their whole length** — a divergent strain,
which is §5's business and not novelty — while **~85% have no significant relative in the
largest nucleotide database there is**. And bigger databases buy hits, not identifications:
across four cumulative protein indexes the residue's "any hit" runs 77.7 → 92.4 → 97.6 →
98.9%, flat only against a ceiling, while at bitscore ≥ 50 it runs 0.00 → 0.00 → 0.36 →
1.15% and is still tripling. **The headroom is in turning weak hits into identifications, not
in finding more weak hits.**

**What the minimal production index still needs, and it is disk.** A DIAMOND database of full
UniRef50 is ~30 GB against ~15 GB of writable disk left here after the nucleotide indexes, so
UniRef50 currently enters as streamed nested 2% and 10% subsets. A run against full UniRef50
or nr needs roughly **160 GB**. The streaming recipe in `experiments/subsample_uniref.py`
already does the hard part — 41,243,524 records, verified complete against UniProt's own
count — so only the disk is missing. That one machine closes requirement (iii) properly and
puts the measured far end on the cost curve instead of an extrapolation. It is the single
highest-value piece of infrastructure in this spec.

---

## 3. Tier 2 — assembled but unidentifiable

Contig-level identification reverses the read-level picture, and the tool needs both indexes
for that reason: **protein-only identifies 26.30% of contigs, nucleotide-only 5.73%** — at
read length the nucleotide index does the work, at contig length the protein index does, by a
factor of 4.6. Overall 34.78% of contigs and 49.13% of bases identified.

**Tier 2 is real but thin, and the tool reports it that way.** 46,157 contigs cannot be
identified at all (14.3 Mb, median 254 nt, median coverage 1.12×). Only 3.87% — about **1,784
contigs** — clear both 500 nt and 2× coverage, which is where assembly and classification
begin to work at all (PR #5). Those 1,784 are the right size for per-contig characterisation;
the rest are reported as a count and a length distribution, not as candidates.

---

## 4. Tier 3 — the residual reads

The 22.90% residue, measured by mapping back to assemblies rather than to public databases:

| pile | % of library | across seeds |
|---|---|---|
| abundant but unplaceable (tier 2's read-level population) | 8.42 | 8.33–8.51 |
| maps a contig from one of the other 115 samples | 9.97 | 9.86–10.04 |
| genuinely missing — nothing public, nothing here, nothing in the other 115 | 10.90 | 10.76–11.07 |

The first two overlap, so they do not sum.

**Cross-sample material is context, never a subtrahend.** Subtracting the other 115 samples'
contigs would make cross-sample recurrence zero by construction and destroy a field the
per-read report needs. It is reported as a reference axis and nothing else.

Subtracting the sample's *own* data is what turns the residue from "merely undeposited" into
"genuinely rare". That is the narrow replacement for the claim withdrawn in PR #13.

### 4.1 The per-read evidence report

One row per surviving pair. Every field ships with two calibrations measured on the same
data: a **null** (the field on composition-matched shuffles of the same reads) and a
**known-positive** (the field on reads the cascade did subtract). A field that reads the same
on both carries no information about the residue, and that is its honest calibration rather
than a reason to hide it.

| field | what it says | calibration that must travel with it |
|---|---|---|
| `closest_relative` | best translated hit, however far: bitscore, identity, alignment length | never thresholded; score tracks database proximity |
| `frame_evidence` | stop-free frames, longest ORF per frame, **and how many frames are candidates** | at 107 nt the true frame is the only stop-free frame in just 0.528 of coding reads |
| `gc`, `entropy3` | composition | shuffled null + subtracted-read comparison |
| `repetitive` | recurring 16-mers, tandem period | separates "no information" from "information we cannot place" |
| `artefact` | adapter, primer dimer, homopolymer, N content | residue is enriched for these by construction |
| `chimaera` | split-mapping evidence | **26% false positives on known background** — reported, never used for a call |
| `motif` | Pfam-A profile hits, six frames, per-domain E-value | the one field that can reach a relative too divergent for DIAMOND |
| `mean_error_prob` | from the quality string | lets an unexplainable read read as "probably instrument" |
| `recurrence` | occurrences across the run and the other 115 samples | sketch only overstates, so "rare" is safe |

**Nothing is combined into a score.** `cannot discriminate here` is a legitimate row, not a
gap — it is what honest calibration of a 107 nt read produces about half the time.

### 4.2 The axes of the signal-versus-noise measurement

Five dimensions, each either already measured or cheap, chosen because each is continuous
and each has a noise source whose signature differs from biology's.

**1. Coverage × identity × significance, never hit versus no-hit.** Two numbers per index,
never collapsed to a boolean — but **three**, because coverage × identity on its own is what
produced the withdrawn reading in §2.1. The "residue read aligning over a third of itself at
89%" that this axis was built on was a chance match, and without the E-value beside them
coverage and identity are a pile of noise that looks like a gradient. A ribosomal read
aligns end to end at 100% identity at E ≈ 0; a residue read with a *significant* hit aligns
over 0.96 of itself at 90.1%; a residue read's *median* hit is 40 bases at 88% at E = 0.24
and means nothing. All three numbers travel together or none of them is reportable.

**2. A read's trajectory across the index-size ladder.** Something recognised only by the
largest index is further away than something caught at 0.67 GB. That is a continuous
divergence estimate, and the per-tier per-read outcome is already computed and written to
the result files — it is currently unused. Costs nothing new.

**3. Recurrence shape across the other 115 samples.** Two-sided, which is rare here: a real
organism recurs with a pattern, random error recurs never, a systematic artefact recurs
everywhere. The project has only ever used recurrence as a yes or no. Shape — how many
samples, how evenly, alongside what else — is the unused half. Still a reference axis and
never a subtrahend.

**4. Mate disagreement as a per-sample error-floor estimator.** Instrument error is
independent between mates; organism-level absence is shared. That orthogonality is dead as
a detection rule (no pair rule approaches one per million) but it measures *this sample's*
noise floor for free, which is exactly what the calibration requirement in §6 needs.

**5. Profile reach.** Pfam-A profiles are the only axis still alive past 20% divergence,
where DIAMOND and every nucleotide index read zero (§7). The gap between "no sequence hit"
and "yes profile hit" is itself a distance measurement.

**That first experiment has been run, and the answer is cloud.** It asked whether the residue
is *structured* or *diffuse* in the joint space — whether clusters appear that behave like a
genome. **No cluster count beats 1 in 5/5 seeds** against a null carrying the residue's own
mean and covariance, and multidimensionality buys **0.013 AUC** over the best single axis
(library coverage, itself partly circular), with composition and quality fields together
reaching 0.686.

So the joint-space view ships as a **density with reference clouds, and never as a
clustering**. The axes above stay, because a per-read report with its calibration beside it
is the deliverable; what does not ship is any claim that the residue separates into
populations. Nothing downstream should be built on structure that was looked for and is not
there.

---

## 5. For known material: new strains of existing viruses

Where a read or contig *is* identified, the tool reports the strain question rather than
stopping at the taxon. Measured constraints on what that can deliver:

- Below about **3.8% divergence every assembler fuses lineages at any depth**; MEGAHIT
  separates from roughly 5% up (PR #12, 2,205 assemblies, 21 pairs, 5 seeds, 0 false mosaics
  in 630 controls).
- The reads carry enough information: PenguiN's raw output holds clean multi-kb pieces of
  both lineages at 2% divergence, buried among ~1,000 mostly mosaic contigs. What is missing
  is **choosing the right paths**, not signal.
- Coverage pooling across lineages is dead: it fabricates mosaics (PR #5).
- Pair scores cannot separate two hosts' swarms — two HBV hosts share 709 exact 150-mers.
  Only population structure can.

So strain calls ship with the divergence at which they were made, and below ~4% the tool
reports a fused lineage as fused rather than as a strain.

---

## 6. Calibration and versioning — the two run-time requirements

**Per-sample calibration run.** Every score is referenced to *this sample's own background*:
the shuffled nulls, the subtracted-read comparisons, the chimaera false-positive rate and the
error profile are all recomputed per sample. A threshold imported from another run is a
threshold measuring the other run.

This stopped being a precaution and became a measurement. Two RNA libraries from **the same
plant, the same group and the same naming series, eight months apart** differ like this
(`results/dna_rna_pair/`): rRNA **38.0% → 67.1%**, Kraken2 PlusPF-16 at confidence 0.1
**22.7% → 40.7%**, virus-identifiable **12.0% → 6.5%** as ribosomal sequence crowds it out.
The ribosomal fraction nearly doubled and the identifiable fraction with it. **So every
removal percentage in this repo is a statement about one sample on one day**, and should be
read as "on SRR23973193" rather than as a property of the method — including the holdout
spike's survival, which rides on a quantity that moved 29 points.

A second constraint follows from the same comparison: **tier order is a property of the
library, not of the tool.** On the matched DNA library — which is a virome, 0.01% ribosomal
and 99.1% unidentified — 0.24 GB of viral protein identifies **8.5×** what the 16 GB
nucleotide index does, while on the RNA library of the same water the 16 GB index wins
six-fold. A cascade tuned on one and pointed at the other spends its largest index to remove
under one percent, so the tier order is chosen per sample by the calibration run and is not
baked in.

**Operating points include the ones that look like plumbing.** Mate 2's non-template trim
moves the Kraken2 removal fraction by **1.7 points** between 14 and 16 bases — the size of
effects reported here as findings — because the prefix length is a distribution and Kraken2's
confidence is a fraction of a read's minimizers (`results/nontemplate_prefix/`). The trim
ships as a declared parameter with its curve, like every threshold.

**The sample's own assembly is versioned.** It is the most powerful thing the tool subtracts
(30.33 marginal points from 15 MB) and it is produced at run time, so two runs of the same
sample would otherwise disagree about what is known. The assembly gets an identifier, a
digest and a manifest entry, and every removal number names the assembly version that
produced it.

Both are in addition to the repo-wide rules: `--seed`, a run manifest beside every result
(seed, git SHA, parameters, timings), release builds for timings, five seeds minimum for any
headline, and spread reported rather than averaged away.

---

## 7. Declared blind spots

These are measured, not cautious. The tool states them in its own output.

**Dimensionality does not move either wall, and the spec must not imply it does.** The
recognition limit is a property of the indexes, not of the decision rule, so adding axes
cannot shift it. And at one in a million there is nothing to aggregate — the molecules are
not there. What §4.2 buys is a deliverable that is honest at depth, which is where the
1–5% divergence band this project cares about actually sits. It does not buy abundance.

**No read-level rule gets near one false positive per million.** Peptide absence gives
160,774 false flags per million at the strictest expressible threshold (PR #11/#13).
Read-level homology gives 8–10% sensitivity at 93,478 flags per million at bitscore ≥40, and
0% sensitivity at 1,489 flags at ≥80. Single-read features — ORF-likeness, k-mer composition,
codon usage, GC — all lack novelty polarity (PR #8). The pair is a multiplier on aggregation,
not a rule: no pair rule reaches one per million either.

**Recognition stops between 10% and 20% divergence.** The holdout spike — three real genomes
simulated into this library's own measured insert distribution (median 194, p95 363) and error
rates (0.00089 R1, 0.00117 R2), walked away from themselves, pushed through every index:

| divergence | spike pairs surviving tier 1 |
|---|---|
| 0% | 0.0% (0.0–0.0) |
| 5% | 0.6% (0.0–2.0) |
| 10% | 12.7% (6.4–19.2) |
| 20% | 82.2% (73.8–89.8) |
| 30% | 99.5% (98.2–100.0) |

That band is this pipeline's measured boundary between *already known* and *new*. Which tier
does the recognising, as removal falls away:

| divergence | sample's contigs | Kraken2 viral | PlusPF-16 | viral proteins | Swiss-Prot |
|---|---|---|---|---|---|
| 0% | 65.5 | 100.0 | 64.0 | 97.5 | 97.5 |
| 10% | 21.8 | 28.1 | 0.4 | 82.9 | 83.0 |
| 20% | 0.2 | 0.0 | 0.0 | 17.4 | 17.6 |
| 30% | 0.0 | 0.0 | 0.0 | 0.5 | 0.5 |

The nucleotide side dies first; the protein index is the last thing standing. That is the
argument for keeping a protein tier in a cascade whose cheapest and most effective tiers are
nucleotide.

Two caveats that ship with the ladder. Uniform substitution is a crude model of real
divergence, so read the rungs as **ordinal**, not as family distances. And two of the three
spiked genomes are already in the sample — PMMoV at 285× and human astrovirus as a 6,823 nt
contig at 11.7× — so part of the 0% removal is the sample recognising itself. Human
rhinovirus 1 is the only clean holdout. The ladder's *shape* is unaffected because the
per-tier table separates own contigs from public indexes.

**Detection below about one in a million is not available at all.** At that abundance a virus
contributes ~8.7 read pairs in an 8.7 M-pair run; the survival ladder implies roughly one
reaches the per-read report at 10% divergence and seven at 20%. That is arithmetic on top of
per-read behaviour measured at a resolvable spike rate — a sample holding 8.7 pairs cannot
resolve a rate, and pretending otherwise was the error of PR #11. Detection is
abundance-limited; **characterisation at depth is where the 1–5% divergence band sits**, and
that is the honest headline.

**Removal is a lower bound on "known".** An index cannot distinguish sequence with no
relative anywhere from sequence whose relative sits below the matching threshold. The sweeps
are the honest form of that limit.

**Subtraction is not free.** SURPI+ measured pre-processing removing 0.43% of true viral
reads. Anything this cascade removes could be a divergent member of a known family, and at a
one-in-a-million target a small false-removal rate matters more than the headline removal.

**The cost curve has not flattened, so the far end is not quoted.** Each tenfold step of
protein index still buys ~20 points at the largest size reachable here. On the nucleotide
side a reference-content ladder shows the last doubling still buying 15–26 points, but the
instrument is weak and we say so: reads concentrate over few taxa (the top taxon carries
49.0% of classified reads, the top ten 86.4%), so that step is largely one dominant taxon
being in or out of the subset. The curve's shape is measured; its value at 100% is not, and
the tool reports measured points only.

---

## 8. Interface

### 8.1 Command line is the contract

```
hapax run SAMPLE.R1.fq.gz SAMPLE.R2.fq.gz --seed N --out RUNDIR
```

`RUNDIR` holds the manifest, the tier table, the per-read rows as JSONL, the per-contig
rows, the calibration block and the assembly version. Everything downstream — including the
web interface — reads those committed files and nothing else.

### 8.2 Visual results reuse the existing generator

`experiments/make_infographic.py` already builds a one-page figure from committed result
JSONs, deliberately as a generator rather than a hand-written page: a figure with numbers
typed into it is a number without provenance, and this project has withdrawn claims whose
prose drifted from its outputs. Missing runs are left out rather than filled in.

The production figure is that generator pointed at a `RUNDIR`, not a parallel build. New
panels are added to it.

### 8.3 Web interface

A thin reader over `RUNDIR`. It adds navigation and nothing else — no number appears in it
that is not in a committed file.

- **The breakdown** as one stacked bar: identified / assembled-unidentifiable / residual,
  split coding and non-coding, each segment labelled with its interval and its operating
  point.
- **The cascade** as the tier table, with the threshold sweep live on a slider, so the reader
  sees the factor of twelve inside one database rather than being told about it.
- **The cost curve**: removal against index bytes, measured points only, with the
  unflattened far end shown as absent rather than extrapolated.
- **The divergence ladder** from this run's calibration spike, labelled ordinal.
- **The residue table**, filterable, one row per surviving pair, with every field's null and
  known-positive beside it and `cannot discriminate here` rendered as a value.
- **The joint space** of §4.2's axes, as the one view that is the actual deliverable: each
  surviving read a point, each axis selectable, the shuffled null and the subtracted-read
  population drawn underneath as the two reference clouds. This is where signal against
  noise is read off, and it replaces any panel that would have shown a novelty score.
- **Blind spots** as a fixed panel, not a footnote.

Three rules for the interface, all inherited from how numbers died in this project: no
figure without its operating point; no pile without its interval; no score without the
calibration run behind it.

---

## 9. What this is not

It is not a novel-virus detector, and it does not emit a novel-or-not call. It is an
instrument that says, for one sample, how much is known and by what, how much assembles,
how much can be characterised and how much cannot — and prices each answer. The deliverable
is the breakdown, the cost-versus-completeness curve, and a multidimensional measurement of
signal against noise in which every axis carries its own null and its own per-sample
calibration. The residue is reported as residue.
