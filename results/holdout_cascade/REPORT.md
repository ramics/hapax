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

`strains` names **54-71 subjects** across the five seeds, over ~1,940 searched residue reads
each. **Across all five seeds exactly one subject anywhere reaches two reads** - every other
one of the 322 is supported by a single read, and none reaches three - at a median alignment
of **41-46 aa** and identities from 35.4% to 100%. Those are the numbers a 145 nt read
produces against a protein index when it is matching by chance, and the field this project has
already retracted a claim in. The subject the withheld genome's own contig matched,
`YP_010771404.1`, is **absent from all five**; it is searched for by name in the generator
rather than assumed missing, because "we did not see it" and "we did not look" is precisely
the bug this experiment found three of.

`layers` sorts 2,000 residue reads into `placed` **1-3**, `distant_relative` **53-69**,
`recurs_and_unplaceable` **716-822**, `undetermined` **1,107-1,230**. Its known-positive arm,
the same 2,000 fields computed on already-identified reads, reads `placed` **259-313** and
`distant_relative` **542-583** - so the layer assignment is working and pointing the right way,
by two orders of magnitude on `placed`. It is simply looking in a pile the virus was removed
from two stages earlier.

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
the reference and long predates the cutoff - reads **99.905% identity over 6,349 nt at 99.47%
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
`coverage 0.25-0.50, identity 60-70%` cell: **69 contigs, 0.91% of the removed pile.** At
nucleotide level it is coverage 0.7103 at 69.7%, a cell holding **137 contigs, 6.64%**. Either
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
percentile of its length band** - or the **98.68th percentile of the whole removed pile**,
the two populations being different and both reported because quoting whichever is larger is
how a field gets flattered. **It is not a gene-count artefact**: binned by number of
significant subjects, the median gap is 0.0 in every bin, including the 5-subject bin that
holds 4,399 of the 7,601 contigs and has a mean of 1.9 - and the candidate has 5 subjects.
A gap requires several genes but several genes do not produce a gap.

**The field exists at protein level only, and that is a constraint on it rather than a
caveat.** On the same contig against the nucleotide index the gap is **exactly +0.0 points**:
per-subject and any-subject union coverage are the *same number*, 0.7103. The mechanism says
why - at nucleotide level the candidate has one best genome and every significant HSP lands on
that genome, so the any-subject union has literally nothing to add. Across its length band the
nucleotide gap has a median of 0.0 and only **2 of 94** contigs exceed 10 points at all. The
gap measures *genes drifting onto different protein records*, which is a protein-space
phenomenon by construction. **So it is not a second axis on top of coverage - it is the
protein tier's coverage read twice**, and a pipeline carrying it must compute it where the
divergence is scored and not where the reads are.

**That figure is a correction, and the first version of it was impossible.** This section
first reported **-0.8 points at the 22.34nd percentile**, which cannot happen: any-subject
coverage *contains* per-subject coverage, so their difference is non-negative by construction.
The cause, found by the consolidation thread, is that `union_coverage` applied no E-value
floor while `union_coverage_any` applied one, so the two sides of the subtraction were not two
definitions but two **filters** - and the gap came out negative on 394 of 7,734 protein rows
and 508 of 2,269 nucleotide rows, to -0.50. **A sign that cannot occur is the cheapest
possible tell and it was sitting in the committed file**; the reason it was not caught here is
that the negative rows are a few percent of a column nothing was plotting.

With the floor applied to both sides and re-derived from the raw search output: **the protein
figure does not move at all** - 0.4624 against 0.8286, +36.6 points, 95.85th percentile of its
band - because every HSP in the candidate's protein hits already cleared 1e-5. The nucleotide
figure moves from an impossible -0.8 to a clean 0.0. **So the protein-only claim survives and
is now stated in its stronger form**: not "smaller at nucleotide level" but *zero* there.

Two things about the percentile that would otherwise mislead. The candidate sits at the
**64.89th percentile** of the nucleotide band rather than the 22nd, and that number is close to
meaningless because the band's median gap is 0.0 - it is a tie, and where a tie falls is an
artefact of sort order, which is a bias this repo has already had to withdraw a result over.
And `coverage_gap_null.py` was never affected, because it floors both sides deliberately: it
independently produced the identical 0.4624 and 0.8286, which is why the +36.6 was safe to
quote throughout.

**The same invariant then failed a second time, on one row, and that is the more useful
failure.** After the floor was added, 41,270 of 41,271 gap cells were non-negative and **one**
read -0.0453. The cause was `>` against `>=`: the floor was written twice and the two
spellings disagree on an HSP at **exactly** E = 1e-5, so one side kept it and the other
dropped it. **A filter written twice will differ on the boundary even when both authors mean
the same thing** - which is the argument for one shared function taking one `evalue_below`,
not for being more careful.

**So the check is now asserted where the number is computed, and it raises.** A negative gap
means the two unions disagree about which HSPs exist, which makes everything derived from
either unreadable, so it is not a row to drop or a count to report - it stops the run. Worth
recording that neither failure was visible to any tell this project relies on: *the control
moving with the test* and *an explicit census* both work by reading a distribution, and
nothing was plotting this column. **An impossible sign needs an invariant, not a reader.**

**The scoring script now reads its search output through `experiments/search_hits.py`, and
the switch found two more impossible values in this same committed table.** The two unions,
the best hit and the gap are built from one read of one row set, under one `passes_floor`,
and `coverage_gap` raises before any fraction is formed, so the filter cannot be written
twice again. What the shared parser's units then exposed:

- **Nucleotide best-HSP query coverage was above 1.0** - 1.0009, 1.0011, 1.0013, 1.0021,
  1.0026 and 1.0046 on six leak contigs - because `length` counts gaps while `qlen` does
  not. A fraction of a contig cannot exceed the contig. Taken from the alignment's own
  `qstart`/`qend` it cannot.
- **The protein version of the same field was residues over nucleotides**, so it capped at
  exactly 1/3 and read 0.4662 where the coordinate span is 0.4624. That was already known
  and fixed here; the nucleotide side was not, and it is the same bug with the opposite
  sign.

Neither field is quoted anywhere in this report, which is exactly why they survived: the
numbers it quotes are the union coverages, which were always built from query coordinates.
**Three impossible values in one table, all in columns no figure used.**

**A third thing was arbitrary rather than wrong.** `attributed_to` names whichever reference
covers most of a contig, and **38 of 35,328 contigs are covered equally by both**, so the
label was whichever the search listed first. Changing to the shared union relabelled 12 of
them. All 38 sit far below the 0.50 attribution threshold - **0 ties above it** - so no arm
changes and nothing in this report moves, but a label decided by file order is the
array-order bias this repo has withdrawn a result over, so the row now carries
`attributed_tie` and the manifest carries the count.

**The re-run the hash change forced reproduced the published table exactly.** Every figure
in `generated_tables_s1to5.txt` is byte-identical after a fresh 40-minute run of all five
searches: the cell of 69 contigs (0.91%), the +36.6 points at the 95.85th percentile, the
nucleotide 0.7103 on both sides, every `gap > 10/20/30` count. The only differences in
`contig_rows_s1.tsv` are the 12 tie labels and 2,072 gap cells moving by 1 in the fourth
decimal, because the gap is now rounded once instead of being the difference of two rounded
numbers.

**Where it gets suspicious, and the limit is mine not the field's.** Four cuts - gap > 30
points, identity 50-70%, any-subject coverage >= 0.70, length >= 1500 nt - select **7
background contigs, 0.092% of the removed pile, 0 of 34 leak contigs**, and the candidate
would be the eighth row. That is a seven-row shortlist containing the withheld genome. **But
all four thresholds were chosen after seeing where the candidate sits**, on a positive set of
size one, so eight rows is what favourable tuning can achieve and not a validated operating
point. The two numbers that were *not* tuned are the ones to carry: the leak arm's three
zeroes, and the 95.85th percentile.

#### The shuffle control: one broken arm, one clean one, and the field stands uncontrolled

This field's control history is the most instructive thing in the report, so it is given in
the order it happened rather than tidied into a conclusion.

**The field was reported with its missing control named** - a union is the single statistic a
pile of weak chance HSPs inflates without limit, because each one contributes interval. The
gene-level thread ran that arm on Berlin within the hour and reported a kill: a *shuffled*
contig looked more than twice as likely as a real one to show a gap above 30 points, which
would have put the candidate's +36.6 at the null's p90. I recorded that as a kill.

**That kill has been withdrawn, and the reason is a bug this repo has already had once.**
Their `dinuc_shuffle` built successor lists, shuffled them, walked, and **on a dead end
returned the input unchanged** after twenty retries. That fired on 324 of 11,361 contigs, and
of the 304 "shuffled" contigs carrying a significant hit, **all 304 came from that set. Zero
genuinely shuffled contigs cleared significance.** The null was a copy of the test, so what
looked like a null exceeding the test was real material compared against itself.

**My arm on this library agrees with their corrected answer, and it is the one that stands:**

| arm | contigs | HSPs | with a significant hit | best E-value |
|---|---|---|---|---|
| real | 22,107 | 148,419 | **22,084** | 0.0 |
| dinucleotide-shuffled | 22,110 | 9,608 over 5,335 contigs | **0** | 5.5e-4 |

0 walk failures, 0 returned unchanged, verified length-preserving and
dinucleotide-count-exact before the run. **So no shuffled contig produces any significant
hit, and therefore none produces any gap at all.**

**Their corrected arm, re-run on Berlin with the fixed shuffle, agrees and is the independent
replication:** **0 of 11,361** shuffled contigs carry a significant hit against **89.53%** of
real ones, and the whole shuffled library produces a single 47 aa alignment at **E = 8.45e-4**,
which clears nothing. Two libraries, two indexes, two scripts, same answer.

**But zero is a bound, not a rate, and their framing is the one to carry: 0 over 11,361 bounds
chance at roughly 1 in 11,000 rather than establishing zero** - mine at roughly 1 in 22,000 on
22,110. A field that fires on a few hundred contigs per library needs a null measured at a
resolution finer than that bound before "not chance" means "not chance at this shortlist's
size", and neither arm reaches it. What both arms do establish is that the gap is not produced
by composition alone.

**This is the same failure mode, in the same repo, for the second time in two days, and a
reader of both reports should not mistake them for independent findings.** `results/
holdout_run/` (PR #29) discarded a complete 80-job run over exactly this: naive Eulerian
walk, dead end, fallback `return seq`, 70.4% of the shuffled arm unshuffled, and a null that
read like the test. That was found, fixed, documented and committed two days ago - and then
**reimplemented from scratch in another thread, where it regressed to the same fallback.**

**Why mine was right is the whole lesson, and it is not that I was careful.**
`coverage_gap_null.py` *imports* `dinuc_shuffle` from `holdout_run.py` instead of writing its
own. That is the only reason this arm is correct. A control's implementation belongs in one
place and is imported, never rewritten, because the second writer does not know which
non-obvious condition the first one was fixing. The shuffle is now a single
`experiments/shuffle_control.py` that raises rather than returning anything, which is the
structural fix; "be careful" is not one. One detail recorded there because it is the kind of
thing a reimplementation loses: the last outgoing edge must be chosen among a vertex's
*distinct* successors, or a homopolymer vertex picks itself and poly-A contigs exhaust every
retry.

**And the explicit census is what distinguishes the two cases.** An arm with no significant
hits and an arm that silently copied the test both produce a table that *looks* fine - theirs
looked populated and was wrong in the direction that manufactures a kill. Counting contigs
with any hit, contigs with a significant hit, and the best E-value separates them, which is
why those three numbers now sit beside the distribution rather than behind it.

#### So where the field actually stands: not refuted, not established

| | |
|---|---|
| the withheld genome | **+36.6 points**, the **98.68th percentile** of the real arm |
| shuffled null | **empty** - no significant hit on any of 22,110 contigs, so no gap |
| leak arm, 34 known contigs, same library | median, mean **and max 0.0** |
| real contigs at or above the candidate | **292**, 1.3% of the removed pile |

**And the rate on real identified background, from their corrected Berlin pass, sizes the
shortlist rather than the effect.** Over the identified pile the gap has a **median of 0.38
points**, **14.51%** above 10 and **5.81%** above 30 - so a 30-point threshold sits at the
**94.2nd percentile**, and the candidate's +36.6 falls between that and the 99th. That agrees
with my own length-matched 95.85th percentile on this library, measured independently. The
operating consequence is blunt: **a threshold there fires on a few hundred identified contigs
per library**, which is the same order as the ~292 the field yields here. The gap is a
shortlist-sized filter, not a discriminator, and it should never be quoted without that number
beside it.

That is a clean negative arm on both sides and a positive set of **one**. The honest word is
**uncontrolled from the other side rather than cleared**: what has been shown is that the
gap is not chance, which was the specific objection, and what has *not* been shown is that
it is a property of novelty rather than of multi-gene architecture. A large gap is
**necessary** for multi-gene architecture and the leak arm says known multi-gene material in
this water does not produce one - but those 34 contigs are a handful of single-ORF-dominated
tobamoviruses, so that arm is thin in exactly the way that matters.

**The refinement does not help, and that result is unaffected by the retraction**, because it
was measured on my own arm throughout. Conditioning the gap on a minimum HSP length moves the
candidate from the 98.68th to the 99.27th percentile while discarding 96% of the real pile,
and above 150 aa it excludes the candidate's own 135 aa coat gene, so the gap would then be
built from a subset of its real genes:

| HSP floor | real contigs still scorable | real p90 gap | candidate percentile | real at or above |
|---|---|---|---|---|
| 0 aa | 22,084 | 1.65 | 98.68 | 292 |
| 100 aa | 18,243 | 0.91 | 99.16 | 155 |
| 200 aa | 5,397 | 0.88 | 99.15 | 47 |
| 400 aa | 953 | 0.68 | 99.27 | 8 |

A floor is shrinking the thing measured, not sharpening a contrast. **The shortlist the field
yields on this library is ~292 contigs, not the 8 that favourable tuning produced**, and 292
is the number to carry.

**And the leak arm's three zeroes rest on 34 contigs from a handful of
single-ORF-dominated tobamoviruses**, so "known material has no gap" may be partly a
statement about which viruses happened to be already known in this water.

### Corroboration from the gene-level thread, and what of it survives its own retraction

`results/divergent_known/` (the `Subtract the known, score the unknown` thread) measured the
same kind of pile on Berlin with a per-gene instrument. Its shuffled arm carried the bug
described above, so everything of theirs that rested on that arm is withdrawn and is **not**
quoted here. Listed rather than silently dropped, because an absent corroboration that a
reader remembers seeing is worse than a named one:

- *"existence of a significant alignment separates real from shuffled 89.52% against 2.68%, a
  factor of 33"* - withdrawn. The true null is zero, so the separation is total and the
  factor is undefined rather than 33.
- *"conditional on a hit, almost nothing about the alignment separates real from shuffled"* -
  withdrawn, and **not replaceable**: with the fixed shuffle no shuffled contig carries a
  significant hit at all, so there is nothing to condition on. That comparison was real
  material against real material.
- *"the fraction of the contig covered is the one field above chance, at 26.34%"* - withdrawn,
  and with it the one external support this report had for ordering a shortlist on coverage.
- their per-column conserved-core result, which used a copy of the same shuffle - withdrawn.

**What survives is the half that never touched the shuffle.** Over 11,361 identified
non-ribosomal contigs of at least 300 nt, 10,170 clear E < 1e-5 at 60 aa with a **median
identity of 45.1%**, and only **4.3% sit at 90% or above**. That is real-arm-only and
consistent with the 40.0% median here once their Swiss-Prot arm and 300 nt floor are
accounted for, so **"the removed pile is mostly distant relatives" is still measured twice,
independently.** Their per-gene role result also stands - it is an ordering within real
contigs - and it is the ground-truth confirmation reported further down.

**The consequence for this report is that I now have less external support, not more.** The
coverage axis was the one recommendation here with a second instrument behind it, and that
instrument has withdrawn its reading. Coverage is still what the ground truth in this
experiment points at - it is what turns a 93.46th-percentile identity into a 0.91% cell - but
that is now **one measurement on one library with one positive**, not two converging ones.
Stated plainly because the earlier version of this section claimed the convergence.

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
defensible is the **cell**, not the ranking inside it: a section that says *69 contigs were
removed at family-level distance over a quarter to a half of their length, and no
already-known contig in this library looks like that*. The red box, honestly built, is a
shortlist with a count on it.

**The gap between those two statements is the whole remaining question**, and it is not
rhetorical: a 70-row list is useful and a 70-row list ordered wrongly is not, so what orders
it is the next measurement rather than this one.

**The coverage gap is not recommended for it, and that is a separate statement from the
withdrawn kill.** Its status is **untested**, not cleared: the kill against it was withdrawn
because the null that produced it was broken, which removes an argument against the field
without supplying one for it. This library cannot supply one either - with no significant hit
on any of 22,110 shuffled contigs, its null has nothing in it to compare against, so "the gap
is not chance here" is as much as these data can say. The corrected Berlin null has since
landed and says the same thing on a second library - 0 of 11,361 - so the chance objection is
answered twice and **the field is still untested in the way that matters**: what is missing is
a second *positive*, and a null at a resolution finer than the few hundred contigs per library
a 30-point threshold fires on. **No ordering should be built on it until then**, and nothing in
this report's recommendations depends on it.

So **what should order the shortlist is still open.** The coverage axis is what the ground
truth here points at, and it now has one measurement behind it rather than two.

**The shortlist survives the widest index available, and that was the open question.** These
distances are measured against RefSeq viral plus Swiss-Prot only, which is a narrow index, so
every identity here is a **lower bound on recognisability and an upper bound on divergence**,
and the obvious objection was that a shortlist built on it is a shortlist of records missing
from a small database. The `Subtract the known, score the unknown` thread has now run the
candidate's divergent panel-mates against **core_nt, 319 GB**, and the attrition is partial
rather than total:

| against core_nt | panel-mates |
|---|---|
| collapse to strain level - were only distant because the index was narrow | **25** |
| stay under 80% identity **or** under half covered | **83** |
| near-reference controls, same pass | **135**, median **97.55%** over **98.8%** of their length |

So roughly a sixth of the panel is index artefact and the majority is not, while the control
arm returns essentially complete matches in the same pass - which is what says the 83 are a
property of the sequence and not of a search that was failing for everyone. **The shortlist
shrinks by about a sixth against the largest nucleotide index there is; it does not dissolve.**
The candidate itself does not move at all, because its distance is identical against today's
unfiltered index. The joint table above was sent to that thread, and this is its answer.

**This recommendation has since been implemented, by that thread rather than here.** `hapax.py`
on `main` now carries a `divergent` stage over `experiments/divergent_known.py`, whose docstring
cites this experiment's numbers as the reason it exists. What this report adds that that stage
cannot measure on the Berlin library is the ground truth: **whether the output contains the one
contig we know to be novel, and where in it.** That is the measurement above, and it is the
only validation of that stage anyone can currently make.

#### So it was run, and the answer is yes — with the position stated, not just the presence

`divergent_known_s1.json`, that stage run blind over this library's 22,110 identified
non-ribosomal contigs against the rolled-back indexes. **The withheld genome's contig is in the
output**, and its row carries exactly the fields SPEC §5 now requires:

| field | value |
|---|---|
| identity, alignment, E-value | **66.5%** over **620 aa**, E = **5.44e-282** |
| nearest subject | `YP_010771404.1`, an RdRp, at **0.984** of the subject's length |
| coverage of the contig | **0.466** |
| divergence class | *family-level or looser; the subject is a nearest neighbour, not an identification* |
| genes, with Pfam | polymerase 66.5% `RNA_replicase_B`, coat 77.8% `Levi_coat`, maturation 53.5% `Phage_mat-A` |
| spread across genes | **24.3 points** |
| read support | mean depth 30.9, breadth **1.0**, evenness 0.718, max 93 |

**Three things make this a validation rather than a sighting.**

**It got into the gene-level arm by contig LENGTH, not by identity.** The arm is capped at 600
contigs per side and the cap is applied longest-first, so at 3,990 nt the candidate is in
comfortably. That matters because this report measured identity ordering to be *actively
wrong* on this very pile — 93.46th percentile, 7,104 contigs ahead of it — and a length cap
cannot inherit that bias. Had the cap been applied by identity the candidate would have been
at risk of falling outside it.

**Its own control arm is empty where it sits.** Of the 600 gene-level rows on the divergent
side, **321** have contig coverage below 0.5 and identity below 70%; of the 600 on the
`control_near_reference` side, **zero** do. So the cell is not discriminating *within* the
divergent arm — and it is perfectly discriminating against known material, which is the
direction that was in question.

**And it is not near the top, which is the honest half.** Sorted as that table sorts, by
identity ascending, the candidate is **431 of 600**; ranked by polymerase identity among the
120 rows that have a polymerase, **100 of 120**. The 321-row cell is *not* comparable to this
report's 70-contig cell either: 70 is out of 7,601 contigs over the whole removed pile, while
321 is out of 600 already selected for being below 90%. **Two denominators, and quoting the
smaller fraction against the larger population is exactly how a shortlist gets flattered.**
The stage's value is that the genome is *present and characterised*, not that it is ranked
first.

**The null is clean and is reported as a bound.** 0 of 22,110 dinucleotide-shuffled contigs
carry a significant hit against **96.71%** of real ones, which bounds chance at roughly 1 in
22,000 rather than establishing zero. That number is a correction: the first run of this
stage reported **175 of 22,110, 0.79%**, and it is withdrawn, because `divergent_known.py`
carried a **third local copy** of the dinucleotide shuffle whose fallback returns its input
unchanged. It now imports `experiments/shuffle_control.py`, and the output carries a null
census — offered, shuffled, dropped-as-invariant — so the denominator the percentage is taken
against is visible rather than assumed. Validated before the run: 22,110 of 22,110 shuffled,
zero length violations, zero dinucleotide-count violations, zero identical to input.

**Provenance note, since this run straddles a commit.** Its manifest records `tree_dirty:
true` with a `script_sha256` of `37aed98930b5d236`, which matches `divergent_known.py` at
committed HEAD byte for byte; the dirt is the result files beside it. That is the case the
content hash exists for.

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

## The cascade: five seeds, and the three-seed first pass that was a defect not a result

`generated_tables_s1to5.txt` holds the generated table; the seeds are disjoint 200,000-read
blocks of **one read pool**, so the spread between them is subsampling noise: not library
variation, not biological variation, and narrower than either by an unknown factor.

| | five seeds, min-max |
|---|---|
| identified (tier 1) | **66.72-67.57%** |
| residue | 32.43-33.28% |
| tier 2, assembled but unidentifiable | 14.32-14.68% |
| tier 3, unplaced by anything | 18.11-18.65% |

**Where the removal comes from is the familiar shape and is the context for the finding
above.** The sample's own identified contigs are again the largest single subtrahend at a
**25.42-26.07% marginal** contribution, with ribosomal sequence at 9.33-11.53%, the three
Kraken2 databases at 1.54-1.61%, 12.77-13.13% and 0.01-0.02%, and the two protein indexes at
0.76-0.78% and 3.84-3.93%. Each figure is that tier's marginal contribution in the order the
cascade consults them, so they are not interchangeable with the alone-column. **The withheld
genome was removed by that first term** - this library's own assembly, consulted against the
rolled-back indexes - which is why the efficiency finding and this experiment's finding are
the same fact seen twice: an index built from the sample explains the sample, including the
parts of it nobody has named.

**The sample's own UNIDENTIFIED contigs are not part of that, and the output says so rather
than relying on the reader knowing the rule.** `tiers.json` records the unidentified observe
set with `"subtracted": false`, and the tier table carries the note *"tier 2 is observed and
never subtracted"*; the 66.72-67.57% union excludes it. Worth stating because the distinction
is easy to lose when both piles are "the sample's own contigs", and because the scale of what
it protects is large: **subtracted, that set would remove a further 21.06-21.41% of the
library on its own** across the five seeds, which is where tier 2's 14.32-14.68% comes from
once tier 1 has run first.

**And that is exactly why protecting it is not sufficient.** The line between the two piles is
drawn by whether a database recognises the contig, so the withheld genome - recognisable at
66.5%, novel by 30 points - lands on the *identified* side and is removed. The material this
experiment lost was never in the protected pile. **Nothing in the design protects the
recognisable-but-distant**, which is the whole finding above, restated from the other end.

The eight-month drift caveat applies with its usual force, and more so here: this is a
different plant on a different continent, so **read every percentage as "on SRR23421544"**.

**The first pass produced three seeds, not five, and the cause is a fifth defect.**
`subtract_known.py`'s `--block` stride - how far seed N skips into the FASTQ - defaults to
**1,500,000**, which is the stride for SRR23973193's 8,731,224 pairs. This library has
**4,290,429** reads, so seeds 4 and 5 start at 4.5 M and 6.0 M, read nothing, and the cascade
exits `no pairs read`. Seeds 1-3 removed 66.31-66.94% and seeds 4-5 failed.

Each seed failed loudly on its own terms, which is the system working. **What did not work is
that the run directory still looked like a five-seed result** - five `run_s*/` directories,
five manifests, and the absence visible only in `tiers.json` not existing for two of them. A
reader taking the directory at face value, or a figure script globbing it, would have
reported the spread of three seeds as the spread of five.

**And `hapax run` did not expose `--block` at all**, so the orchestrator could only ever use a
stride tuned to one library - a parameter that decides *which reads are measured*, behind a
default. Both are fixed: `hapax run --block` now exists, and
`experiments/holdout_cascade_run.sh` sets it to 1,000,000 (blocks at 0, 1 M, 2 M, 3 M, 4 M,
the last ending at 4,200,000 of 4,290,429) and **refuses to start** if the arithmetic does not
fit the library's actual read count, which it now counts rather than assumes.

All five seeds have been re-run under the corrected stride. **The three-seed numbers are not
carried forward**, because changing the stride moves every seed's block, and the convention
does not let one directory hold results from two parameterisations.

## Five defects found on the way, which matter beyond this experiment

The first three are one failure mode: **an optional input is absent and the code either dies
at the end or returns a number that reads as a measurement.** The fourth is its sibling and
the one with the longest reach: **a provenance field the conventions require is computed and
then not written**, so the run cannot say what code produced it. The fifth is the same shape
again from the other end: **a parameter that decides which reads are measured was not exposed
at all**, so it kept a default tuned to a different library - described under "The cascade"
below, because it is what produced three seeds out of five.

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
   line of the summary with nothing written. `motif_rate` also returned `0.0` when nothing
   had been searched; both now say `not measured`.

   **Correction, from the `Make hapax fail loudly` thread: my claim that every other optional
   field in that function already used an `x = None` guard was wrong, and my own committed
   runs show it.** `chimaera_split_rate` is a fraction over all rows, so with no chimaera
   reference on disk no row carries `split_disagree`, the numerator is zero, and the field
   reports **`0.0`** - while its own control, `chimaera_split_rate_on_subtracted_reads`,
   correctly reports `None`. `run_s1-3/fields.json` carry exactly that pair.

   **That combination is the inverse of this project's own heuristic and worth naming.** The
   rule is *read the control first*, on the grounds that a broken measurement shows up as the
   control moving with the test. Here the control is the honest field and the test is the
   lying one: read the control and you see `None` and know nothing ran; read the test and you
   see a rate. So "read the control first" catches a control that mimics the test and misses a
   control that is absent beside a test that defaults - and only counting what was searched
   catches both. Which is the same conclusion the shuffled-arm census reached today from the
   other direction. `crispr`, `rfam` and `cross_sample_recurrence` are `None` rather than
   `0.0`, so they are absent rather than false.

   No number in this report reads any of those fields, and the full list is in
   `results/fail_loudly_audit/`.

4. **`subtract_known.py` defines the content hash the conventions require and never calls
   it.** `provenance()` returns `git_sha`, `tree_dirty` and `script_sha256` and carries a
   docstring explaining that the hash is the only one of the three that pins the code when
   the tree is dirty - citing a real incident where a fix landed between two seeds of one
   multi-seed run. **It is never called.** The manifest instead records `"script":
   "experiments/subtract_known.py"`, which is a path, and a bare `git_sha`. `hapax.py` reads
   the stage's manifest for the hash, does not find one, and honestly records `"script_sha256":
   "absent from this stage's own manifest"` - so the orchestrator obeys the convention and the
   stage does not supply it. Same shape as bug 1: a function that exists, is documented as
   load-bearing, and is not wired up.

   **It bit exactly the case its own docstring predicts, and the final run shows both halves
   of why the hash exists.** Because report and data files were committed while the seeds ran,
   the five run-level manifests record **two different git SHAs** - `453ebcb` for seed 1 and
   `674c62a` for seeds 2-5 - which under a SHA-only convention is indistinguishable from a
   code fix landing mid-run, the exact case the convention forbids reporting as a five-seed
   result. **The run-level `script_sha256` settles it: all five read
   `950f3c05f696a6dc…`, byte-identical.** The two commits touched `results/` and nothing
   executable, which is confirmable rather than asserted: `hapax.py`,
   `subtract_known.py` and `residue_report.py` hash the same at `453ebcb`, at `674c62a` and in
   the working tree. **So the hash does its job and the SHA does not** - this directory is one
   code version on five disjoint blocks, and the field that proves it is the content hash.

   The `cascade` stage's own manifest is the half that still fails: it carries
   `git_sha` and no hash, so at stage level the question cannot be answered from the file at
   all, and `hapax.py` records `"script_sha256": "absent from this stage's own manifest"`
   rather than inventing one. The `reads` stage does propagate a hash and it is identical
   across all five seeds (`11b925110f3324ce`). Not fixed here, because seeds were in flight
   and editing a script mid-run is the one thing the conventions forbid outright; handed to the
   `Make hapax fail loudly` thread, which has it in PR #30.

**Audit, bugs 1-3: nothing committed is affected.** No directory under `results/` contains a
manifest from `hapax.py` or a `hapax run` directory, and there is no committed `strains.json`,
`layers.json` or `reads.jsonl` anywhere. Every committed cascade number comes from
`subtract_known.py` invoked directly, whose manifests name the indexes it genuinely used.
`results/holdout_cascade/` is the first committed `hapax run`, and it carries the fixes.

**Audit, bug 4: it reaches most of the archive.** Of **125 committed manifests** under
`results/`, **29 carry a `script_sha256` and 96 do not**. 36 record a `-dirty` git SHA, and
**20 are both dirty and hashless, which means their provenance cannot be recovered from the
file at all** - the SHA names a working tree that was not the code, and nothing else pins it.
Those 20 are in `nucleotide_ladder/` (5 seeds), `spike_abundance/` (4 of 5),
`partial_alignment/`, `partial_alignment_contigs/`, `peptide_absence_noise/` and a handful
more. The largest hashless groups are `novelty_blocks/` (20), `mate_rho/` (11),
`insert_profile/` (7) and `peptide_absence/` (7).

**Audit, bug 5: no committed seed group elsewhere shows the same failure, and that was
checked rather than assumed.** The `Make hapax fail loudly` thread ran
`experiments/seed_degeneracy_audit.py` over all **150 multi-seed file groups** under
`results/`, testing each for a zero-byte seed, a seed under half its group's median size, and
byte-identical seeds. **Nothing committed is empty or short anywhere.** Two groups are
byte-identical across all five seeds and both are correct by design and already stated in
their own reports: `block_position_bias/block_fractions_s*.npz`, where the block statistics
are exact over the whole library and the seed varies only the random-sample control in the
per-seed JSON, and `holdout_cascade/run_s*/blind_spots.json`, a fixed panel `hapax` copies
into every run directory. **What that check cannot see is a seed that read the right number
of reads from the wrong region** - which is the failure here, since seeds 4 and 5 did not
read a short block, they read none - and nothing closes that from committed output alone. The
`--block` guard is what closes it going forward: it refuses to start unless the arithmetic
fits the library's counted reads.

This is a **list, not a re-run**, and it is not a claim that any of those numbers are wrong -
only that the convention's own test cannot be applied to them. It is also not caused by
`subtract_known.py` alone; the field is simply absent from most writers, so **the convention
is honoured by 29 of 125 manifests**. Worth knowing before the next multi-seed campaign,
because the convention exists for precisely the situation this run hit by accident.

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
