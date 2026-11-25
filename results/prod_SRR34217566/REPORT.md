# Production from one install and one command, over an unseen library: SRR34217566

**Municipal wastewater, Szeged, Hungary, collected 2024 (the archive gives the year only and
does not say influent or effluent), RNA, 11,589,124 pairs (PRJNA1282049). The archive states a
150 nt read; the file is 151, and the run measures it rather than taking either on trust.**

The Yeongdok pass (`results/prod_SRR23702504/REPORT.md`) asked whether the pipeline would run
end to end on water nobody had measured. It did, after five things were edited in code. This
pass asks a narrower and harder question: whether the **documented install and the documented
command** reach a finished report.

```
conda env create -f environment.yml && pip install .
hapax fetch-db --all
hapax run --config run.toml --seed N --report
```

`experiments/prod_oneshot.sh` runs those three and nothing else between them. It writes its own
`run.toml`, appending the `[databases]` block that `fetch-db` prints, so the handover between
the commands is the tool's own output rather than a path typed twice. What this report
records is the honest ledger of what a person still had to supply, and a refusal that was
correct is recorded as a success.

*One water on one date. The same plant eight months apart moves its identifiable fraction from
22.7% to 40.7% (`results/dna_rna_pair/`), so nothing below is a statement about Szeged or about
Hungary.*

## Why this sample

The constraint list is the Yeongdok list with one addition, and the addition is what was hard
to satisfy.

- Not the discriminating dataset (SRR23973193), whose every per-library fact is committed, so
  running it proves only that its own defaults survive.
- Not the holdout (SRR23421544), whose value is a withheld genome.
- Not the New York correctness fixture.
- Not a member of `results/crossplant/panel.tsv`: a library must not be subtracted against a
  panel it is inside, because its own reads are abundant in that panel by construction, giving
  54% false growth against a 6% swap control (`results/crossplant/`).
- **Not a member of any of the six studies the panel was drawn from either**, which
  SRR23702504 was not required to satisfy. Its plant was in the panel while it was not, which
  is what kept the paired within-plant stratum alive. Here there is no such stratum, so the
  cross-plant stage is absent with its reason rather than run on a pairing that would compare
  continents while reporting a plant.

It also adds something no library here has had: **2 x 150 nt**. Berlin is 107 and 93 after its
trim, and the holdout is single-end despite its archive layout, so SPEC 4.3's mate merge, the
largest measured gain in recognition at distance at 10.5% to 38.4% at 20% divergence, has never
been testable on real reads in this repo.

## What the run refused, and what had to be supplied

| # | Where | What it refused | What was supplied | Was the refusal right? |
|---|---|---|---|---|
| 1 | `library_facts.resolve_args` | every stage, on an accession with no entry | a committed entry in `experiments/libraries.json`: the archive's own facts, the read paths, and the two panel pointers, each with its basis | Yes, and it is the one refusal a new library cannot avoid. |
| 2 | `stage_refs.resolve` | the truth-panel, positive-control and cross-plant stages | `truth_panel: null`, `crossplant_panel: null` and `positive_control_absent_because`, each with a reason | Yes. Absent with a reason is not zero recovery, and the cross-plant absence is the price of the condition that made this library worth running. |
| 3 | `experiments/trim_from_prefix.py` | naming a mate-2 trim, on an ample measurement | `--trim-r2 0` with the measurement as its written basis | Yes, and it is the most interesting result of the pass. See below. |
| 4 | the classify preflight | the whole run, 3 minutes in | nothing: the preflight itself was wrong | The ordering was right and the check was wrong. See below. |
| 5 | `nontemplate_prefix_map.py` | the prefix measurement, 29 minutes in | `--libraries` in place of `--r1/--r2` | Yes. That script resolves read paths through the registry precisely so a caller cannot point it at one library's reads while it reads another's facts. |
| 6 | nothing | nothing: the per-gene workup reported `0 divergent contigs, 0 near-reference controls` and DIAMOND then exited on a blank query file | the taxa stage's output in place of `strains.json` | No, and this is the one fault in the pass that did not announce itself. See below. |
| 7 | `hapax run`, after every stage of seed 1 | writing a page, on a second committed file the renderer reads | `results/web/controls.json` in the payload | Yes, and it is the same refusal as launch 5's, one import further out. |

### The trim is still the thing that cannot be inherited, and this time the measurement was good

The one fact production cannot carry across libraries is how much non-template 5' sequence to
cut from mate 2. Berlin's is 14 nt, measured; a silent default of 14 on a library with no
prefix deletes fourteen real bases from every mate, in the direction of less recognition.

On Yeongdok the instrument could not be used at all: treated effluent anchored 0.13-0.14% of
pairs to SILVA, under `trim_from_prefix.py`'s floor of 2,000, so the trim was stated rather
than measured. That is not what happened here. Over 5 x 100,000 pairs this library anchored
**17,930 to 18,192 mate-2 reads per seed**, nine times the floor and two orders of magnitude
more support than Yeongdok, and read a mate-2 5' clip **mode of 0 in 5 of 5 seeds**.

It refused anyway, on two of its four conditions:

- the mode holds only **21.4% to 21.9%** of anchored reads, under the floor of 25%
- mate 1's mode is **0 as well**, where the rule requires the two to differ

Both are the same observation: the clip histogram is **flat**. Mate-2 5'-minus-3' mean runs
+4.96 to +5.06 against mate 1's -3.08 to -3.25, an excess with no spike under it. A fixed
length of non-template sequence makes a spike; priming and end quality make a tail
(`results/nontemplate_prefix/REPORT.md`).

So the refusal is correct and the reading is clear: **this library has no non-template prefix,
and "no prefix" is not a prefix of length zero.** The rule is written to detect a prefix and
cannot certify its absence, which is why a person has to state the 0. The bias runs the safe
way: adapter read-through clips 3' ends and inflates the within-read control, so every
failure of this test errs towards refusing to trim, and not trimming leaves kit sequence in
the reads rather than deleting template.

`--trim-r2-basis` now exists so that the sentence travels with the value into the run config
and onto the page. Before it, a stated trim recorded "stated on the launcher", which says that
a person chose but not what they chose on, and a stated trim is the one library fact here that
no instrument can certify.

This is the live gap in one command on unseen water. Everything else the instance did, it
did by itself.

### The classify preflight, which cost three minutes and is worth recording

The first launch refused itself after about three minutes and $0.03, on its own check that
Kraken2's `classify` helper is findable. The environment was correct: Python 3.13.18, MEGAHIT
1.2.9, diamond 2.1.10, Kraken 2.17.2, minimap2 2.31 all resolved and logged. The check was
wrong. It ran before `pip install .` and with the working directory at `/opt/hapax`, where
`import hapax` reaches the repository's root `hapax.py`, the legacy driver, which has no
`tools`. A correct environment failed on an import that could never have worked.

The package and the driver share a name, which is why `src/hapax/cli.py` loads the driver by
path and never by name. A plain `python` launched from the checkout is the same trap one layer
out. The check now runs after the install, from `/tmp`, under `python -I`.

Recorded because the ordering worked as designed: the refusal cost three minutes rather than
four hours, and three of this pass's four launches ended the same way.

### A flag name is the cheapest thing in a launch to check and was among the last

The second launch died 29 minutes in, having completed the install, all five reference sets
(39 GB) and 1.8 GB of reads, because the template passed `--r1` and `--r2` to a script that
takes neither. `tests/test_oneshot_flags.py` now reads the long flags each EC2 template passes
against the long options the command itself declares, through `--help` for the experiments
scripts and through the installed entry point for the `hapax` subcommands. It needs no instance
and no network, and it reproduces that failure.

### Two files in one run directory answer to the word strains

`divergent_contig_workup.py --strains` documents `known_strain_divergence_s*.json` and
reads it under a `taxa` key. The package passed it `strains.json`, which is the `strains`
stage's own output and a read-level table of subjects with no such key. Both sides are a
`Path`, so nothing a type checker reads was wrong, and nothing raised: selection found no
contig in either arm, the ORF file came out empty, and DIAMOND exited 1 on a blank query
after eleven stages. The log line was **`0 divergent contigs, 0 near-reference controls`**,
which is a statement about the water. With the file the script documents there are 15
divergent contigs and 30 near-reference controls in this library.

It is the fault in this pass worth keeping, because it is the only one that produced a
plausible number rather than a traceback, and because the repo's standing rule covers it
exactly: a plausible zero and an absent measurement look identical, so an empty search is
recorded as a census. `divergent_contig_workup.py` now writes that census and runs no
search rather than handing DIAMOND an empty file; `tests/test_stage_inputs.py` compares
the filename the package passes against the filename each script's own `--help` declares,
for every file-valued flag in the glue.

## What reading the report found

The stages exited zero on every seed and the generator refused nothing. Everything below
was found by reading the page the run produced and the tables beside it, and not one of them
would have been found by reading an exit code. Each is now either covered by a test,
asserted where the quantity is computed, or made impossible by resolving a declaration
against the run instead of asserting it.

They fall into three shapes, and the shapes matter more than the count. **A declared
absence that was false**: the spectrum input, the controls that said their stage had not
run, the ring legend for marks that were not drawn, the field that printed its own
measurement as a dash, four table sections whose files had moved. **A number computed and
discarded**: the chimaera flag, the borrowed-correction flag, the closest row's own
alignment length. **A sentence true of another library**: six of them, each measured
somewhere and typed here. A thirteenth was only typography, and is recorded because it
came from the same place as the rest: the incomplete-run note appended a full stop to each
stage's reason, and a reason written as a finished sentence printed two.

### Reading the page found four claims it had no right to make

The run's numbers all carry a seed, a script hash and an operating point. The sentences
around them did not, and four of them were facts about the discriminating dataset written
where a parameter belongs.

- The heading read `<accession>, untreated wastewater` and the line under it `Untargeted
  RNA sequencing of untreated wastewater`. The archive does not say whether this run is
  influent or effluent, and this library's own committed basis says so in as many words,
  because the ribosomal fraction differs by two orders of magnitude between them. Both
  now come from the water's own facts.
- The watchlist panel read `No row reaches 1 pair in 100,000` directly under its own top
  row of 1 in 44,000. True of the library the sentence was written on. Now counted.
- The composition ring read `0% of that is viral` above a table listing eight named virus
  taxa and 448 reads. A share that is present and below the printed precision now reads
  as a bound.
- The host paragraph read `Those 0 reads are the enteric viruses sewage is expected to
  carry, matched over about 34 amino acids`: an identification the host group does not
  support, an expectation this run did not measure, and an alignment length measured on
  another library, in a sentence whose own number was zero.

The convention this repo arrived at for numbers, and then for the vocabulary file, had
not reached the renderer's own prose. These four were found by reading the page rather
than by any test, which is the argument for reading it.

### The headline figure was empty, and the page gave a reason that was not true

The distance spectrum is the figure the whole class-A output is, one mark per assembled
virus on a continuous identity axis. On the finished seed it drew nothing, and three of
its controls read "the spectrum stage was not run on this library: the divergent-contig
workup is not yet a stage of `hapax run`". The run directory next to that page held
1.7 MB of that stage's output.

The page descriptor is assembled by `hapax report`, and the line declaring this input was
a sentence typed when the statement was true. `inputs.spectrum` is required to be exactly
one of a path or `{"not_run": reason}`, so that nobody-said and does-not-apply cannot
render the same, and the hardcoded `not_run` satisfied that contract while being false.
It now resolves by looking in the run directory for the stage's own output, which is the
only way the answer is about the run in hand.

This is the harder of the two directions to notice. An input declared present but absent
is a refusal, and the generator raises on it. An input declared absent but present is a
quiet omission, and a declared absence renders as scrupulousness: the page said what it
could not show and why, in the register the rest of it uses for genuine limits.

### A third of the Methods table was blank, and nine rows were blank for the wrong reason

Every headline cites a Methods entry of six fixed rows, and a row holding a number
measured on one library is stored under that library's accession so it cannot print on
another page. On a library nobody has written about, 28 of the 90 rows therefore read
"Not measured on this library.", including the control arm and bias direction of figures
that had just rendered.

Nine of those rows state no number in the first place. Which control arm a figure carries
and which way its bias points are properties of the procedure, and they hold whatever
water is in the tube. They were blank because the sentence saying so had been written
inside the per-library text rather than beside it. Each such row now carries a committed
design sentence, printed on every page, and the per-library number is appended where one
is recorded; where none is, the row says that the number is not recorded for this library
rather than that nothing was measured. The two are different claims and the page had been
making the stronger one.

### The same shape twice more, in the tables and in a caption

The production tables printed `mate-2 trim used: unmeasured` and `block stride: -`. Both
are read from the library's registry entry, which deliberately records neither, because
both are measured on the instance and that is the point of the one-command path. The run's
own manifest held 0 and 2,500,000 and the page printed the first in its operating point
two sections down. The registry answers what was supplied; only the run answers what ran,
and the table now says which of the two answered.

The headline of the divergence section counts assembled viruses, 34 of them, and the
figure under it marked 15 with no denominator anywhere near it. The two counts come from
two piles: the gene-level workup selects its own contigs, so the page now states both side
by side and does not call the second a subset of the first, which is an overlap nothing
here has measured.

### The chimaera number was computed, carried into the page, and read by nothing

`CLAUDE.md` requires a chimaera number beside every contiguity claim, because a longer
contig that is wrong is worse than a short contig that is right. The divergence section
counts complete genomes carried at over 100x depth, and the per-contig
subject-disagreement flag that answers it was computed by the workup stage for both arms,
kept in the page's payload, and read by no line of the renderer. The Methods entry
describing it was the one entry on the page that nothing cited, which the build reports
and no test failed on.

Surfacing it exposed a second fault in the number itself. The flag can only fire on a
contig carrying two or more named genes, and the summary divided by every contig in the
arm, so it read 40.0% against a 36.7% control. Over the contigs that could be flagged the
same rows read 100% against 84.6%. The first invites the reading that the flag separates
the two arms a little; the second says what the entry's own bias row already said, that it
tracks how the assembler joined genes and not how far they sit from a reference. Both
denominators are now written out, the conditional one is `None` rather than `0.0` where
nothing was testable, and the numerator is asserted against its own denominator where it
is computed.

### A directory is not a seed

Fetching the run while its last seed was still going left four manifests and a fifth
directory holding one file. The stage census read "0 stages failed across 5 seeds" and
"0 of 5 complete" off that, while the cascade table two sections above it, which reads
numbers rather than counting directories, correctly said four. The same run described
over two denominators, with the larger one in the sentence that looks like the summary.
Counts are now over the seeds that produced a manifest, the partial ones are named, and a
results directory holding no manifest at all is a refusal rather than a table of empty
cells.

### Four whole table sections were missing, and the file ended tidily

The generated tables are this repo's defence against prose drifting from its outputs, so a
generated file that is silently short is worse than no file. Three of the records the
tables read sit at the results root on the per-library path and inside each seed's own
directory on the one-command path, and one of them is spelled differently there. Each
section was guarded by a bare "if the file loaded", so this run's tables lost the
whole-library assembly, the contig split, and the comparison of what was measured off the
reads against what the registry states, and still ended with a stage census and a
provenance line.

With them back, the comparison they carry earns its place on its own: the archive states
a 150 nt read and the reads are 151, the registry's stated value is what ran because
stating it is the override, and the table now prints both sides and which one ran rather
than reconciling them. Paths tried are recorded whether or not one hits, and a section
with no file anywhere is a line naming the file it needed.

### The named band lost every Kraken2 tier to a database release date

`hapax fetch-db` names each Kraken2 directory after the release it downloaded, so this
run's cascade tiers arrive as `k2_viral_20261022.d` and `k2_pluspf_16_GB_20261022.d`. Two
panels held a literal map of which tiers can return a name, keyed on the undated names the
per-library scripts use. Nothing matched. Every Kraken2 read moved out of the named band
and into recognised-but-unnamed: the band read 2.24% against the cascade's own 4.07%, and
the reconciliation panel went on printing "agrees with the row above" beside two figures
fifteen points apart.

A tier is now classified by what it is, with an unrecognised one reported as unknown
rather than defaulting into the larger band, and whether the two named rows agree is read
off the gap instead of asserted.

Two faults in the same figure came out with it. The sub-label under the unnamed band read
"mostly own-sample contigs (2.26%), rest adapter and low complexity" over a band of
14.18%, where the 2.26% named only the non-ribosomal self-assembly tier and left out the
9.73% of own-sample ribosomal contigs, and the artefact filters it blamed for the
remainder had removed 0.61%. And the ribosomal band carried that index's own reach while
the band beside it was computed from marginals, so the stacked parts summed to more than
the stack. Every band is a marginal, the sub-label lists each tier with its own share, and
the five are asserted to sum to the library where they are built, because a stacked bar
normalises by its own total and will draw any set of numbers as a full width.

### Six more sentences were about another library

The page's numbers are read from the run. A set of sentences around them were not, and
each was true of the sample it was written on.

- The cannot-move table said the mate-2 trim "ships as two measured runs, 14 and 16, 1.7
  points apart" on a page whose own trim is 0, said assembly flags "move the residue by
  8.79 points", and reported an agreement between five search methods measured over
  sixteen contigs of another library. Each reason is general and now prints everywhere;
  each number is stored under the accession it was measured on, as the Methods rows
  already were.
- The distance spectrum's reading note said the genes of a divergent contig "disagree by a
  median 12.8 points", which is the figure from an 11,361-contig pile elsewhere. This run
  measured 8.6 against a near-reference control of 0.95, and had printed both two sections
  above.
- The same note taught a three-colour code for rings marking contigs checked against a
  319 GB nucleotide index, over a count of zero. No larger-index block ran, so the figure
  says that and describes no colours.
- The interval widening read "Measured over all 8.7 million pairs" and put its source in
  the final clause. The producer records that the correction is borrowed and applied here
  unmeasured, and the page printed neither. Borrowed now comes before the numbers.
- The linkage stage's own bias statement ended "mate 2 is 93 nt against mate 1's 107 and
  carries 62 scorable 31-mers against 76". Every library's output carries that string, so
  a library of two 151 nt mates was told its mates were 93 and 107. The stage measures
  them now; the page states the run's measured lengths beside the sentence, because runs
  produced before the fix still carry it.

The asymmetry those lengths describe is real and load-bearing where it exists, which is
why the sentence was written. It is a property of one sequencing run, so it is measured
like one.

### A measurement rendered as its own absence

One field shipped a bare number where the calibration table expected a mean and a standard
deviation. The table's fallback for an unreadable value is the same dash it uses for not
measured, so `frame_identifiable_frac` printed empty cells on both arms, at 0.1615 on the
residue against 0.0950 on the null, with "calibrated" in the state column beside them. The
bare number is read now, and a shape the function cannot read at all is loud rather than
blank: a dash is indistinguishable from the honest case, which is what let this sit.

## What the run produced

Every number here is in `TABLES.md`, generated from the committed run directories by
`experiments/report_prod.py`. Spread is the observed range over the seeds and is not an
estimate of dispersion: five seeds are contiguous blocks of one flow cell, and off five of
them this repo's own dispersion statistic reverses its sign.

Identification removes most of the library, and the sample's own assembly is again where the
efficiency is. The whole-library assembly is 45,180 contigs over 25.8 Mb, N50 574,
with 215 contigs at or above 3 kb; identification names 45.9% of those contigs and 50.1%
of their bases. Against the cascade, the ribosomal index takes the largest single bite and
the sample's own ribosomal contigs, 7.2 MB of them, take the most per byte of any reference
in the stack: more than the 16 GB Kraken2 database by four orders of magnitude on that
ratio. This is the third library to say so.

The residue is small here, which is a property of this water rather than of the tool. The
union of every tier removes 94.10% of read pairs over five seeds, 92.40 to 96.54, so the
remainder is 5.90%, against 22.9% on the discriminating dataset and 21.6% on Yeongdok. The ribosomal fraction is the reason: it is
about three quarters of this library, where effluent had almost none. Read every one of
these fractions as "on this sample".

Nothing here is pointed at a person, and the page says so rather than scoring it. Of
the leftover mates that carry a named relative at all, none points at a record whose host
is human-associated, and host confidence never exceeds 0.30 across the matched records, so
no relative is close enough for its host label to transfer. The closest resemblance in the
residue is 96.3% identity over 31 amino acids, which is a ratio of small integers and is
printed with that length beside it.

The watchlist fires low, with its control beside it. The highest row is Leptospira at about
one pair in 44,000, nine of forty-four taxa are detected at the working confidence, and the
same panel over composition-matched shuffled reads detects twenty-three of forty-four at
confidence zero and none at the working setting. The shuffled arm is why the detected rows
are readable at all.

At the contig level the gene-by-gene asymmetry reproduces. Fourteen of thirty-four
assembled viruses sit below 90% protein identity to their closest published protein. Within
one contig the genes disagree by a median 8.6 points against 0.95 on contigs that match
their reference, so a single identity per contig reports whichever gene won on bitscore.
By role, divergent against near-reference: polymerase 80.9 against 98.6, maturation 83.0
against 96.3, coat 81.5 against 98.7. That is the same ordering, on a third library, with
the control arm drawn by the same code.

The chimaera flag carries almost no information here either. Over the contigs it could
fire on it reads 100% on the divergent arm against 84.6% on the near-reference control,
which is what its own bias row says: it tracks how the assembler joined genes, not how far
they sit from a reference. Nothing is dropped on it, and the contiguity claims above now sit
beside it on the page.

Two stages are absent with reasons rather than estimated. The cross-plant stage refuses
because Szeged is not one of the thirty-two plants in the panel and no library from this
plant is, so there is no within-plant stratum; searching it against the existing 530 would
measure how far Hungary is from four other continents, which is a real quantity and not the
one the stage reports. The calibration ladder was not requested, so no sensitivity-at-
divergence figure may be read off this page. Both appear on the report as stated absences.

## What a person still had to supply, end to end

The claim under test is one install and one command on a library nobody here had looked
at. What the pass needed from a person, in order:

1. **The accession.** `experiments/prod_oneshot.sh` stages the reads from SRA in region;
   nothing was downloaded by hand and nothing was copied between libraries.
2. **An entry in `experiments/libraries.json`.** `library_facts` raises for a library it
   does not know rather than falling back, so there is no run without one. Every sample
   fact in it came from this study's ENA `filereport`, with the basis recorded beside it:
   the archive gives the collection year and no month, so the page says 2024 and invents
   no date, and it does not say influent or effluent, so the entry does not either.
3. **The reasons the panels are absent.** Four of them, hand-written: no depositor truth
   panel exists for this water, Szeged is not one of the 32 plants in the cross-plant
   panel, the spiked ladder needs a reference assembly this pass does not stage, and the
   mate-2 trim is deliberately not stated. `stage_refs` refuses a null panel with no
   reason, so these are the price of a stage being absent rather than silently skipped.

What a person did **not** have to supply is the part that had been written into about
twenty scripts as a default. The mate-2 trim was measured on the instance against SILVA
before the expensive half of the pass and came back 0, which is the honest answer for a
library with no non-template prefix and the answer a Berlin default would have got
wrong by fourteen bases in the direction of less recognition. The seed stride was
measured off the read count in the file at 2,500,000 pairs, with 200,000 pairs read per
seed, where Berlin's committed stride would have crowded all five seeds into part of this
library. Two quantities, two names, and the production table prints each with the value the
run used beside the registry's: the archive states a 150 nt read here and the file is 151,
and the table says which of the two answered.

What remains hand-written is the registry entry itself, and it is the last committed file
an unseen library needs a line in. Most of it is sample facts that the archive will
answer for and panel pointers that are absences with reasons. The run configuration
already carries sample facts as an override the page marks as supplied; letting it carry
the panel pointers too would leave an accession and a pair of FASTQs as the whole input,
which is what the claim says. That is not done here, and the Methods numbers for this
library are still unwritten: twelve of the fifteen entries print their design and say
that no number is recorded for this water.

## Cost and wall clock

`LEDGER.md` carries the table, generated by `experiments/run_ledger.py` from the committed
launch manifests and each instance's own writes to S3. Every figure there is a **lower
bound**: a terminated instance leaves `describe_instances` within hours and every one of
them had already gone, so the gap between an instance's last upload and its terminate call
is not recoverable and is not estimated. The previous library's report carries the same
table typed by hand with two rows marked as approximations, which is why this one is a
script.

The pass that produced the committed result is the last row, on `r6i.2xlarge` on-demand
with a 300 GB volume, at **$3.82 for five seeds** against its launcher's own projection
of $6.39 and the 720-minute dead-man timer's worst case of $6.36. It is inside the $10 that
is pre-approved for a single run. The five seeds took 65.2, 65.9, 66.1, 68.7 and 68.3
minutes, so a seed is about an hour on eight cores and the pass is about five and a half
hours of compute on top of staging.

The eight launches before it cost **$5.76** between them and every one of them bought
a defect in staging rather than in the pipeline: the two that wrote nothing to S3 at all
died before their first upload. That is the cost of a one-command path being genuinely one
command: each failure was a line in a shell script that only a clean machine could expose.

A second library on this instance type should cost about what the last row did, with no
further staging.

## Provenance

- `launch_i-*.json`, one per launch, keyed by the commit the payload was built from. The
  launcher refuses to ship a dirty tree, so the commit identifies the code that ran.
- `LEDGER.md` and `ledger.json`, generated, with the rate table and the date it was taken
  written into the script rather than queried, so regenerating the table later cannot
  re-price this run.
- `TABLES.md`, generated by `experiments/report_prod.py` from the run directories. No
  number in this file is typed by hand that is not in there.
- `nontemplate_prefix/nontemplate_prefix_s{1..5}.json`: the trim measurement, which
  refused, with the rule it applied and the direction of its bias.
- `run_s{1..5}/library_measured.json`: what the launcher measured off the FASTQ before any
  stage ran, which is the only answer to what the reads are as against what the archive
  says they are.
- `pip_freeze.txt`, `run.toml`: the environment and the configuration the pass ran on.
- `run_s{1..5}/manifest.json`: per seed, the git SHA, the script content hash, the
  parameters and the per-stage status and timings. Read the per-seed manifest rather than
  the directory: the pages were rebuilt in-container at a later commit than the one the
  stages ran at, and each page records both.

The pages in this directory were rebuilt in-container rather than taken from the instance.
The run executed at the commit in its manifests; every fault in the list above was found by
reading what that build produced and fixed after it. Rebuilding is what CLAUDE.md's rule
about a provenance change forcing a re-run asks for, and it is free: the page generator
reads committed result files and takes about two seconds a seed.
