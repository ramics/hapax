# The full cascade over a genome nobody knew: it is removed as background

**Run:** `experiments/holdout_cascade_run.sh`, which calls `experiments/classify_contigs.py`,
`hapax.py run` (seeds 1-5) and `experiments/holdout_cascade.py`. Tables in this file are
printed by `experiments/report_holdout_cascade.py` from the committed outputs; none is typed.

## The question

`results/holdout_run/` (PR #29) measured the read-level search path on this library: a 145 nt
read from a 30.3%-divergent virus clears E < 1e-5 about one time in sixteen, and the identity
it reports flatters the truth by 13 points unless coverage is read beside it.

That is not what a person running the tool would see. Imo asked the question that is: **run
the committed cascade blind over this library, let it assemble, and say whether the pipeline
would have surfaced this genome as a divergent relative worth reporting** - ideally as a
prominent flag at the top of the report. This answers it, and the answer is no, for a reason
worse than insensitivity.

## The answer in one line

The withheld genome assembles perfectly well and the cascade **classifies its contig as
identified**, which means it is subtracted in tier 1 as already-known background. It never
reaches tier 2, never reaches the residue, and the `strains` stage only ever reads the
residue. **Nothing in the report mentions it at all.**

### What the pipeline reported instead

The two outputs a reader would look at for "is there something new here" are `strains.json`,
which groups the residue by the relative it is nearest to, and `layers.json`, which ranks
residue reads by how interesting they are. Both ran and both are intact. Neither mentions the
virus, and what they do contain is worth stating because it is the chance band and nothing
else.

`strains` names **71 subjects** over 1,949 searched residue reads. **Every one of them is
supported by exactly one read** - no subject reaches two, let alone three - at a median
alignment of **42 aa** and identities from 37% to 100%. Those are the numbers a 145 nt read
produces against a protein index when it is matching by chance, and the field this project has
already retracted a claim in. The subject the withheld genome's own contig matched,
`YP_010771404.1`, is **absent**; it is searched for by name in the generator rather than
assumed missing, because "we did not see it" and "we did not look" is precisely the bug this
experiment found three of.

`layers` sorts 2,000 residue reads into `placed` 3, `distant_relative` 68,
`recurs_and_unplaceable` 822, `undetermined` 1,107. Its known-positive arm, the same 2,000
fields computed on already-identified reads, reads `placed` 280 and `distant_relative` 583 - so
the layer assignment is working and pointing the right way. It is simply looking in a pile the
virus was removed from two stages earlier.

So the failure is not that the signal was too weak to see. **The pipeline reported a clean,
well-calibrated, correct answer about a pile that no longer contained the thing it was asked
about.**

## How it happens

Assembly is not the bottleneck. 4,290,429 single-end reads assemble (MEGAHIT, 13 min) into
35,328 contigs over 22,037,613 nt, and the 4,147 nt withheld genome comes out as **one 3,990
nt contig covering 96% of it**.

Against the indexes rolled back to 2022-05-10, the day the water was collected, that contig
gets a strongly significant hit at both levels:

| level | identity | aligned | union coverage of contig | E-value | bitscore | subject |
|---|---|---|---|---|---|---|
| held-out RefSeq viral protein | **66.500%** | 620 aa | 0.4624 | 5.44e-282 | 812 | `YP_010771404.1` RdRp |
| held-out RefSeq viral nucleotide | **69.718%** | 2,873 nt | 0.7180 | 0.0 | 1,141 | `NC_074562.1` |
| today's unfiltered RefSeq viral protein | 66.500% | 620 aa | 0.4624 | 5.57e-282 | 812 | `YP_010771404.1` |
| today's unfiltered RefSeq viral nucleotide | 69.718% | 2,873 nt | 0.7180 | 0.0 | 1,141 | `NC_074562.1` |

**The last two rows are the drift control, and they are the strongest thing in this report.**
Searched against today's index with nothing held out, this contig returns **the same subject at
the same identity over the same alignment**. So the 30.3% distance is not an artefact of how
the holdout was built: this genome's nearest relative in RefSeq predates the cutoff anyway,
the rollback changed nothing about it, and the finding below does not rest on the index
reconstruction that `holdout_indexes.py` could only approximate to nine records (caveat 2).

Swiss-Prot, which is undated and is also a cascade tier, independently removes it: 36.6% over
456 aa at bitscore 246, still five times the threshold. **The removal is doubly determined** -
dropping the viral protein tier entirely would not save it.

Genome-level truth is 66.45% protein and 69.72% nucleotide, so **at contig length the distance
readout is essentially unbiased** - the +13-point flattery PR #29 measured on 145 nt reads is
gone by 4 kb. That is the one cheerful number in this report, and it is what makes a
contig-level distance worth reporting at all.

It is also why the contig is deleted. `classify_contigs.py` calls a contig identified if
DIAMOND finds a hit at bitscore >= 50. This scores **812**. The removal is not marginal and
**no tightening of that threshold could spare it**: 99.82% of identified contigs score at or
below 812, so a bar high enough to save this genome would stop the tier working at all.

### The control says this is real

The comparator genome from the same library - the one whose nearest relative was already in
the reference and long predates the cutoff - reads **99.9% identity over 6,349 nt at 99.5%
coverage**. So the holdout is not over-removing, the leak arm behaves, and nothing leaked.

Two genomes 30 identity points apart, and the cascade records the same single word for both:
`identified`. That is the whole finding. **The subtrahend is a boolean, and the distance that
produced it is computed and then discarded.**

The project's standing rule - *the subtrahend is IDENTIFICATION, not assembly* - does not help
here, and it is worth being explicit about why. That rule protects a virus **nothing can
recognise**, by routing it to tier 2 instead of deleting it. This virus **is** recognisable,
at 66.5% over most of its length. Tier 2 catches the unrecognisable. Nothing catches the
recognisable-but-distant.

## Would a distance-with-coverage readout have surfaced it? Measured: on a 70-row shortlist, but not at the top

The obvious fix is to keep the distance: report what each removal was removed *at*, rather
than only that it was. Whether that would have worked is a measurement over the pile nobody
looks at, and it was made before being recommended. The answer has two halves that point
opposite ways, and which half you read decides whether this is a usable output.

**Identity alone fails, in the direction nobody expects.** Of 22,260 identified background
contigs, 7,601 have a protein hit significant at E < 1e-5 against the held-out index. Their
median identity to the nearest known protein is **40.0%** and three quarters are below 50%:

| protein identity | contigs | % of significant |
|---|---|---|
| 0-50% | 5,766 | 75.86 |
| 50-60% | 1,015 | 13.35 |
| **60-70%** | **443** | **5.83** |
| 70-80% | 169 | 2.22 |
| 80-90% | 82 | 1.08 |
| 90-95% | 48 | 0.63 |
| 95-99% | 53 | 0.70 |
| 99-100% | 25 | 0.33 |

The withheld genome is at 66.5%, the **93.46th percentile**: it is *more* similar to a known
protein than 93% of what the cascade removed. A rule of "flag the low-identity removals" reaches
**7,104 other contigs first**. Ordering a shortlist by identity is not merely undefensible on
this project's standing rule - it puts the one genome we know is novel near the bottom.

**Adding coverage changes the answer.** Joint counts over the same pile, union coverage of the
contig against identity, printed by `report_holdout_cascade.py` section 4b:

| union coverage \ protein identity | 0-50 | 50-60 | 60-70 | 70-80 | 80-90 | 90-100 |
|---|---|---|---|---|---|---|
| 0.70-1.00 | 3,571 | 656 | 282 | 102 | 49 | 87 |
| 0.50-0.70 | 1,095 | 173 | 64 | 26 | 17 | 19 |
| **0.25-0.50** | 871 | 149 | **70** | 32 | 12 | 18 |
| 0.00-0.25 | 229 | 37 | 27 | 9 | 4 | 2 |

The candidate's protein hit covers **0.4624** of the contig at 66.5%, so it lands in the
`coverage 0.25-0.50, identity 60-70%` cell: **70 contigs, 0.92% of the removed pile.** At
nucleotide level it is coverage 0.7180 at 69.7%, a cell holding **139 contigs, 6.74%**. Either
way this is a shortlist a person can read in an afternoon, not a lottery.

**And the direction is fixed on known material from the same library.** This is what the leak
arm is for, and it is unambiguous: of its **34 contigs, not one** sits below 80% identity at
either level. At protein level 29 of 34 are in `coverage >= 0.70, identity 90-100`, with the
rest at 90-100% identity on lower coverage and a single contig at 80-90%; at nucleotide level
32 of 34 are top-right and 2 at 80-90%. **Known material from this water does not enter the
candidate's band at all.** The crowding in that cell is other *background* contigs, not
recognisable ones.

**What it does not do is rank the candidate first.** Ordered by bitscore within its own protein
cell, the top row is a 4,156 nt contig at 68.0% over 0.458 coverage scoring 881; the candidate's
812 is below it, and below three others at nucleotide level including an 8,190 nt contig at
69.6% over 0.937 coverage scoring 2,181. So this is a shortlist that **contains** the withheld
genome at a readable length. It is not a flag that points at it.

### The one field whose direction is pinned: the coverage gap

Chasing a discrepancy between two of my own coverage figures turned up the best-calibrated
signal in this experiment, so it is reported with the accident intact.

"Union coverage" has two definitions and they disagree by 36.6 points on this contig:

- **per-subject** - union of significant HSPs against the single best-covering subject:
  **0.4624**. Answers *does one known protein explain this contig*.
- **any-subject** - union against every significant subject: **0.8286**. Answers *is any of
  this contig explained by known proteins at all*.

The cause is visible in the hit table. The candidate's RdRp hits `YP_010771404.1` at 66.5%
over 620 aa covering contig positions 1-1850; a second gene hits a *different* record,
`YP_010771401.1`, at 53.5% over 493 aa covering 2478-3938. **A divergent multi-gene genome's
genes each drift onto a different nearest relative**, so no single known protein covers much
of it while known proteins between them cover most of it. Both figures are right; they are
answers to different questions, and the candidate sits at the 16.82nd coverage percentile of
the removed pile under the first and the 50.18th under the second.

**The gap between them is a field, and it has a direction fixed on known material:**

| | gap, percentage points |
|---|---|
| the withheld genome, 3,990 nt | **+36.6** |
| 193 identified background contigs of 2,793-5,187 nt | median 0.1, mean 4.5, max 54.2 |
| **the leak arm: 34 contigs of already-referenced genomes, same library** | **median 0.0, mean 0.0, max 0.0** |

Not one known contig in this water has any gap at all, and the candidate is at the **95.85th
percentile** of its length band. **It is not a gene-count artefact**: binned by number of
significant subjects, the median gap is 0.0 in every bin, including the 5-subject bin that
holds 4,399 of the 7,601 contigs and has a mean of 1.9 - and the candidate has 5 subjects.
A gap requires several genes but several genes do not produce a gap.

**Where it gets suspicious, and the limit is mine not the field's.** Four cuts - gap > 30
points, identity 50-70%, any-subject coverage >= 0.70, length >= 1500 nt - select **7
background contigs, 0.092% of the removed pile, 0 of 34 leak contigs**, and the candidate
would be the eighth row. That is a seven-row shortlist containing the withheld genome. **But
all four thresholds were chosen after seeing where the candidate sits**, on a positive set of
size one, so eight rows is what favourable tuning can achieve and not a validated operating
point. The two numbers that were *not* tuned are the ones to carry: the leak arm's three
zeroes, and the 95.85th percentile.

#### The shuffle control was run, and the field fails it

The control named above as missing - a union is the single statistic a pile of weak chance
HSPs inflates without limit, because each one contributes interval - was run by the
gene-level thread on Berlin within the hour, and **the field does not survive it. The null
is not merely comparable; it is higher.**

| gap, percentage points | real (10,172 contigs) | dinucleotide-shuffled (304 clearing the same bar) |
|---|---|---|
| median | 0.38 | **1.43** |
| p90 | 21.5 | **35.9** |
| p99 | 46.2 | 49.6 |
| above 10 points | 14.51% | **35.53%** |
| above 30 points | 5.81% | **14.47%** |

**A shuffled contig is more than twice as likely as a real one to show a gap above 30
points**, by exactly the mechanism named when the control was requested: several
weak-but-significant HSPs land on unrelated records in different places, none covers much,
and the gap follows. The candidate's +36.6 therefore sits around the **null's p90** and
cannot carry weight on its own. My gene-count control does not catch this and theirs
reproduces mine exactly - median gap 0.0 in every subject-count bin on *both* arms - so that
control was answering a question that was not the one that mattered.

**This is the project's own stated tell, and I missed it in the right direction but not the
right magnitude.** The field was reported with its missing control named, which is the rule;
what the rule did not save me from is that naming a gap in the controls is not the same as
knowing which way it points. Recorded as a kill, not a caveat.

**What survives, and it is not nothing.** A gap near *zero* means one subject explains the
contig, and that direction is solid - so the leak arm reading a maximum of 0.0 across 34
contigs remains a real statement about known material. The failure is asymmetric: a large gap
is **necessary** for multi-gene architecture and **not sufficient** for it.

**The one refinement that could survive** is a gap built only from long HSPs. The shuffled
gaps are assembled from short scattered alignments; the candidate's three genes align over
620, 493 and 135 aa - and note that the third is short enough that a length floor high
enough to exclude chance HSPs would also discard a real gene, which is the tension the sweep
below is measuring rather than assuming away.
`experiments/coverage_gap_null.py` tests this on *this* library, where the positive is a
real genome at a known distance rather than a population statistic, sweeping the length floor
rather than choosing one - because a floor picked after seeing where a single candidate sits
is the same mistake this field already made once. Its result is below.

**And the leak arm's three zeroes rest on 34 contigs from a handful of
single-ORF-dominated tobamoviruses**, so "known material has no gap" may be partly a
statement about which viruses happened to be already known in this water.

### Corroboration and a calibration, from the gene-level thread

Two results from `results/divergent_known/` (the `Subtract the known, score the unknown`
thread), which measured the same pile with a different instrument. They are quoted here
because one of them qualifies the table above rather than just agreeing with it.

**Same direction, different instrument.** Over 11,361 identified non-ribosomal contigs of at
least 300 nt, 10,170 clear E < 1e-5 at 60 aa with a **median identity of 45.1%**, and only
**4.3% sit at 90% or above**. That is consistent with the 40.0% median here once their
Swiss-Prot arm and 300 nt floor are accounted for. So "the removed pile is mostly distant
relatives" is now measured twice, independently.

**The qualification, and it matters for which axis to trust.** They shuffled each contig
dinucleotide-wise and ran it through the same search in the same pass. **Whether a significant
alignment exists at all** separates real from shuffled completely - 89.52% against 2.68%, a
factor of 33. But **conditional on a hit existing, almost nothing about the hit separates
them**: median identity 45.1% real against 43.5% shuffled. Scored as the share of real hits
above the shuffled p90, where 10% sit above by construction: identity 11.98%, alignment length
7.10%, subject coverage 7.52%, bitscore 7.63% - all at or near chance. **The one field
carrying real information is the fraction of the contig covered, at 26.34%.**

That control is about telling a real hit from a chance one, not about the coverage-by-identity
table above, which is on real material throughout. But the two measurements converge on the
same axis from opposite directions, which is the reason to trust it: **their shuffle says
coverage is the only field carrying information above chance, and our ground truth says
coverage is the field that turns a 93rd-percentile identity into a 0.92% cell.** Identity
conditional on a hit is weaker than it looks in both. A shortlist ordered on coverage with
identity beside it is better supported than the reverse, and that is now two independent
reasons rather than one preference.

## So what is defensible

**A ranked shortlist over the removed pile, carrying the distance and the coverage that caused
each removal - ordered on coverage, never scored.** That is the same shape
`results/residue_census/` settled on for the residue - rank it, ship the list with its
annotation, do not score it - applied to a pile nobody currently inspects. The gain is not
subtle: today the number of divergent relatives this pipeline reports is **zero**, because they
are all inside tier 1. With this, the withheld genome goes from silently deleted to a row on a
**70-row** list with a distance and a coverage beside it, in a band that **no recognisable
contig from the same water enters**.

**Still not a flag, and still not a score.** Imo asked for a red box reading "novel virus" at
the top of the report, and the measurement says what the box could honestly be: not a verdict,
and not the top row. Within its own cell the candidate ranks fourth on bitscore at nucleotide
level and second at protein, so a box that fires on the top row would have missed it. What is
defensible is the **cell**, not the ranking inside it: a section that says *70 contigs were
removed at family-level distance over a quarter to a half of their length, and no
already-known contig in this library looks like that*. The red box, honestly built, is a
shortlist with a count on it.

**The gap between those two statements is the whole remaining question**, and it is not
rhetorical: a 70-row list is useful and a 70-row list ordered wrongly is not, so what orders
it is the next measurement rather than this one.

**What should not order it is the coverage gap, which was the best candidate here and is
dead.** It read zero on all 34 known contigs from the same water and +36.6 points on the one
we know is novel, and it was the only field in this experiment that was not a restatement of
database proximity. Its shuffle control, run the same day, put the null *above* the test. So
the honest answer to "what orders the list" is **still open**, and the only axis with
calibration behind it from either side remains the plain fraction of the contig covered:
26.34% of real hits above the shuffled p90 against 7-12% for identity, alignment length,
subject coverage and bitscore.

**The honest size of that list is not yet known, and the direction is one-sided.** These
distances are measured against RefSeq viral plus Swiss-Prot only, which is a narrow index, so
every identity here is a **lower bound on recognisability and an upper bound on divergence**.
Some of the 70 will collapse to "known" against nr or core_nt - and the candidate will not,
because its distance is identical against today's unfiltered index. That pass belongs to the
`Subtract the known, score the unknown` thread, which owns the gene-level instrument, and its
result decides whether the shortlist is 70 rows or 20. The joint table above was sent to it.

**This recommendation has since been implemented, by that thread rather than here.** `hapax.py`
on `main` now carries a `divergent` stage over `experiments/divergent_known.py`, whose docstring
cites this experiment's numbers as the reason it exists. What this report adds that that stage
cannot measure on the Berlin library is the ground truth: **whether the output contains the one
contig we know to be novel, and where in it.** That is the measurement above, and it is the
only validation of that stage anyone can currently make.

**A single whole-contig identity is the wrong unit, and this contig proves it.** Put through
the gene-level thread's `conserved_positions.py`, the contig has **three** scorable genes,
each with a Pfam domain, not the two my own hit table made visible:

| gene | length | Pfam | identity |
|---|---|---|---|
| maturation | 487 aa | `Phage_mat-A` | **53.5%** over 493 aa |
| replicase | 616 aa | `RNA_replicase_B` | **66.5%** over 620 aa |
| coat | 135 aa | `Levi_coat` | **77.8%** over 135 aa |

So the whole-contig figure of 66.5% lands on the genome-level truth of 66.45%, but it lands
there because the replicase won on bitscore - **the spread across its genes is 24.3 points**,
and the architecture is Leviviricetes. A single identity is right here for the wrong reason.

**That ordering is a ground-truth confirmation of their asymmetry, which is worth more than
my gap ever was.** On Berlin, across 189 divergent contigs with two or more roles, maturation
is the most diverged role on 118 and the least diverged on 8, with polymerase running the
other way. This genome - at a distance nobody chose, in a different library, withheld from the
index - lands exactly on that pattern: maturation most diverged at 53.5%, coat least at 77.8%,
replicase between. The matched control behaves too: the leak contig has four genes spanning
99.6-100.0%, a spread of 0.4 points, so which role is "most diverged" there is a coin toss -
which is precisely why a near-reference arm cannot test the asymmetry and this holdout can.

**Rank by the conserved gene with the per-gene spread shown**, which is the project's
standing rule, is therefore the one recommendation in this report that ground truth supports
rather than undermines.

## The cascade, five seeds

See `report_holdout_cascade.py` section 1 for the generated table. Five seeds here are five
disjoint 200,000-read blocks of **one read pool**, so the spread between them is subsampling
noise: it is not library variation, not biological variation, and narrower than either by an
unknown factor.

## Three bugs found on the way, which matter beyond this experiment

All three are one failure mode: **an optional input is absent and the code either dies at the
end or returns a number that reads as a measurement.**

1. **`hapax.py` validated `--k2-dbs`, `--silva` and `--dmnd` and then never passed them.** The
   cascade fell back to `subtract_known.py`'s own default paths whatever the run was given. On
   a machine where those defaults exist the stage *succeeds* and the manifest records the
   indexes asked for rather than the ones consulted - the index named in the provenance not
   being the index used, which is a failure this repo has already withdrawn a figure for. It
   surfaces as a crash only on a container where the defaults are absent, which is how it was
   found. **The existence check made the bug harder to see, not easier, because it proved the
   file was there.** Where a stage validates a path, assert that path appears in the argv it
   builds.
2. **The `reads` stage never passed `--dmnd` either**, so `residue_report.py` used its
   defaults. Where those are absent every `closest_*` field is missing, and the `strains` stage
   reads exactly those fields - so it writes an empty file that reads as *no divergent
   relatives* when the truth is *nothing was searched*. Second time in this project a stage has
   quietly not searched and had its silence read as a negative result.
3. **`residue_report.py` bound `pf_known_rate` only inside its Pfam block.** With no profile
   database on disk, the entire 26-minute per-read pass completes and then dies on the last
   line of the summary with nothing written. Every other optional field in that function
   already uses an `x = None` guard. `motif_rate` also returned `0.0` when nothing had been
   searched; both now say `not measured`.

**Audit of what this could have affected: nothing committed.** No directory under `results/`
contains a manifest from `hapax.py` or a `hapax run` directory, and there is no committed
`strains.json`, `layers.json` or `reads.jsonl` anywhere. Every committed cascade number comes
from `subtract_known.py` invoked directly, whose manifests name the indexes it genuinely used.
`results/holdout_cascade/` is the first committed `hapax run`, and it carries the fixes.

## What this run could not measure

1. **No second mate.** This library is filed `PAIRED` and the archive holds one file.
   Everything SPEC 4.3 gains from the merge - which moved recognition from 10.5% to 38.4% at
   20% divergence, the most valuable free step this project has - is **untestable here**, as
   are mate disagreement and mate linkage. The cascade is run with r1 as both mates and
   `--trim-r2 0`, which is exact for every tier there because each is a min, max or union over
   the two; the per-read report drops mate 2 rather than duplicating it (see `mates()`), and
   the two mate-linkage stages **refuse** rather than report a fabricated positive. So this
   measures a **weaker instrument than the production tool**, in the direction that makes the
   withheld genome harder to find rather than easier.
2. **The held-out index is not byte-identical to PR #29's.** `holdout_indexes.py` rebuilds it
   from the committed `dropped_accessions.json` and checks six counts against
   `holdout_manifest.json`. Protein reproduces exactly (703,330 kept, 19,493 dropped, 4
   unmapped kept); nucleotide reads 19,202 kept against 19,211 recorded, because **upstream
   RefSeq viral dropped 9 records between 2026-11-21 and 2026-11-22** (today's release serves
   19,641 where PR #29 saw 19,650). The script refuses to build on a mismatch and was run with
   `--allow-mismatch`, recorded in its manifest. Direction: nine fewer pre-cutoff relatives can
   only make a candidate look *more* novel, the dangerous way - which is why every distance is
   measured against today's unfiltered index in the same pass as the drift control.
3. **Kraken2 and SILVA are not date-filtered to the day.** The Kraken2 viral database is the
   **2021-05-17** build, genuinely pre-cutoff; the 8 GB Standard and PlusPF builds are
   **2022-06-07**, 28 days late, because genome-idx publishes nothing earlier at that size. A
   2022-06-07 index is built from RefSeq as of that date and the withheld genome was not
   deposited until 2023-07-31, so it cannot contain it; the 28 days can only affect other
   sequence, and only by removing more background. SILVA is 138.1 (2020), pre-cutoff, and no
   ribosomal reference can contain this virus. Swiss-Prot is **today's** release and is not
   date-filtered - UniProt does not serve one that old here - so it is carried as its own
   column rather than merged. One-sided throughout: a post-cutoff index can only make a contig
   look more identifiable, never less, so none of this can manufacture novelty.
4. **One candidate, one library, one distance.** 30.3% is a single rung; the curve between
   0.08% and 30.3% is not measured. A genuinely paired candidate exists at ~30 GB of staging
   per library and would test the SPEC 4.3 parts this one cannot.
5. **The per-read report ran at 2,000 pairs, not hapax's default 20,000.** Measured here: the
   `fold` field alone costs **290 ms per read** (ViennaRNA, 20 dinucleotide shuffles per read
   for its own null), so 20,000 pairs is ~97 min per seed and ~8 h over five. 2,000 is the
   per-arm sample `results/holdout_run/` used. The cost sits in one field, not in the search.

## Cost

Zero. Staging, assembly and all searches in-container: 13 min of assembly, ~37 min per
`blastn -task blastn` pass over 22 Mb of contigs (megablast would be far quicker and would miss
exactly the divergent matches this experiment is about), minutes for each DIAMOND pass, ~6 min
of cascade and ~26 min of per-read report per seed.
