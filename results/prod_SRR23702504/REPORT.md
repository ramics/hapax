# Production over a library nobody had measured: SRR23702504

**Yeongdok plant effluent, South Korea, collected 2020-09-11, RNA, 2 x 151 nt, 32,125,614 pairs
(PRJNA941081).**

The claim being tested is not that the pipeline produces good numbers on this water. It is that
**the pipeline runs end to end on a sample nobody has touched without anyone editing a flag**, and
the deliverable is the honest list of what had to be supplied by hand. A refusal that was correct
is recorded here as a success, not a defect.

*Drift is the largest caveat on every figure in this repo: the same plant, the same assay, eight
months apart moves the identifiable fraction from 22.7% to 40.7% (`results/dna_rna_pair/`). This is
one water on one date. Nothing here is a statement about Yeongdok, about Korea, or about effluent
in general.*

## Why this sample

Three libraries in this repo cannot test genericness, by construction:

- **SRR23973193** (Berlin influent) is the discriminating dataset. Every per-library fact in
  `experiments/libraries.json` was measured on it, so running it proves only that the defaults it
  set are still there.
- **SRR23421544** (the holdout) is the positive control. Its value is the withheld genome, and its
  facts are likewise committed.
- **PRJNA746354** (New York) is a correctness fixture.

What was wanted was a library from a different plant and ideally a different continent, with no
entry in `experiments/libraries.json`, **and not a member of the cross-plant panel** — a library
must not be subtracted against a panel it is inside, because its own reads are abundant in that
panel by construction (54% false growth against a 6% swap control, `results/crossplant/`); the
stage refuses that case and the refusal is correct.

SRR23702504 satisfies all of it and keeps one thing the panel needs: it is **not** in
`results/crossplant/panel.tsv`, while **its plant is** — 12 panel libraries from Yeongdok, three
influent from the same date — so the paired within-plant stratum the cross-plant measurement
requires still exists. `tests/test_prod_run.py` asserts both halves of that, because the
not-a-member half is the one refusal that would invent recurrence.

It is also, deliberately, **effluent** rather than influent, and that turned out to be the single
most consequential difference between it and Berlin (see the trim, below).

## What the run produced

Five seeds, 200,000 pairs each, disjoint blocks at the registry's stride. Every figure is in
`TABLES.md`, generated from the committed run directories; the three that say whether the pass
worked are:

- **Nothing failed.** Zero failed stages across 5 seeds x 13 stages. Three stages are skipped with
  their committed reason because SPEC 13.4 keeps them out of production, and three more declined
  with a reason rather than returning a number (the trim measurement, the truth panel, and the
  four-nickname prefix map before it was fixed).
- **Identifiable material removes 78.42% of read pairs** (observed range over the five seeds
  78.31 to 78.59), leaving a residue of 21.58% of which 11.74% is unplaceable at the read level.
  Berlin reads 77.10% and about 11.1% on the same cascade, so a library nobody had measured lands
  within about 1.3 points of the one every threshold was set on. That is a statement about these
  two waters and not an accuracy claim.
- **The sample's own assembly removes 71.92% alone where the 16 GB Kraken2 database removes 0.66%**,
  the sharpest version of this repo's indexing result so far (Berlin: +30.3 points for 1.4 MB of
  self-contigs). Efficiency, not sensitivity: self-assembly removes 0.0% of a genuinely novel
  virus, measured in `results/crossplant/`.

All five seed manifests carry the same `script_sha256`, so one body of code produced all of them.
They carry no commit, which is bug 4 below.

## What the run refused, and what had to be supplied

| # | Where | What it refused | What I supplied | Was the refusal right? |
|---|---|---|---|---|
| 1 | `library_facts.resolve_args` | every stage, on an accession with no entry | a committed entry in `experiments/libraries.json`: pairs, read length, block stride, read paths, each with its basis | Yes. This is the whole point of the resolver raising, and it is the one refusal that is unavoidable for a new library. |
| 2 | `stage_refs.resolve` | the truth-panel and positive-control stages | `truth_panel: null` with `truth_panel_absent_because` (no depositor genome set for this water) and the same for the spiked control | Yes. **Absent with a reason is not zero recovery**, and the earlier design would have reported a zero. |
| 3 | `experiments/trim_from_prefix.py` | naming a mate-2 trim at all | the value `0` stated on the launcher's `--trim-r2`, with a paragraph of basis and `trim_r2_unverified: true` | Yes, and it is the most interesting result of the pass. See below. |
| 4 | `hapax.py run` preflight (local, against dummy inputs) | the whole run, twice, before any money was spent | the Rust binaries (`linkage`, `census`) built on the instance, and the relabelled assembly written under the file name `crossplant_stage.py` requires | Yes. Both would have been discovered eleven hours and six dollars later. |

### The trim is the thing that could not be inherited

The one fact production cannot carry across libraries is how much non-template 5' sequence to cut
from mate 2. Berlin's is 14 nt, measured; a silent default of 14 on a library with no prefix
deletes fourteen real bases from every mate, in the direction of less recognition.

The committed instrument (`experiments/nontemplate_prefix_map.py`) maps untrimmed reads against
SILVA and reads the soft-clip at the 5' end, with mate 1 and the 3' end of the same mate as
controls. It therefore needs ribosomal content. On this library, over 5 x 100,000 pairs:

- mate-2 5' clip **mode 0 in 5 of 5 seeds**, holding 47.1–52.5% of anchored mates
- mate 1 also 0, at 73.3–78.5%
- resting on **127–144 anchored mates per seed: 0.13–0.14% of pairs**, against roughly 32,000 and
  ~32% on Berlin
- non-zero clips spread flat over 3–73 nt at 1–3% each; 3' clip mean 11.0–13.6 nt

The shape is a library with no prefix. The support is 139 reads. `trim_from_prefix.py` refused on
its own anchored floor (2,000), which is the correct call: at that depth the histogram describes
SILVA's reach into treated effluent, not the library's chemistry. So the trim is **stated, not
measured** — SPEC 4.1 allows stating it because stating it is the override, and a stated value is
never re-checked. `0` is chosen because it is the only value that cannot delete real template; the
residual risk, if this kit does carry a prefix, is kit sequence left in the reads, which costs
recognition rather than inventing it.

The by-product is a third rRNA instrument, committed as `rrna_anchored_fraction`: this water is two
to three orders of magnitude less ribosomal than the Berlin influent.

**The better instrument is obvious and untested**: clip the reads against the sample's *own*
assembled contigs, which needs no SILVA at all. It needs an untrimmed assembly to measure against,
so it is a chicken and egg this pass did not break.

## What I edited, and why each was code rather than data

The measure of genericness is how much had to change in code. Four things did:

1. **`experiments/nontemplate_prefix_map.py`** — knew three nicknames and nothing else, so a new
   accession could not be measured at all. It now resolves any accession through
   `library_facts`. This is the one edit that was unavoidable: a script that measures a library had
   no way to be pointed at one.
2. **`experiments/provenance.py`** — raises rather than writing "unknown", which is what it is for,
   but an instance has a tarball and no `git`, so every script importing it raised there. It now
   takes the commit from `HAPAX_GIT_SHA`, which the launcher sets, **only when git cannot answer**,
   so a machine with a checkout cannot dress a dirty tree as clean.
3. **Three new scripts rather than edits**: `trim_from_prefix.py` (the five seed files to an integer
   or a refusal), `trim_mate2.py` (a trimmed copy for the assembler, since MEGAHIT reads files and
   `path_settings.py` explicitly refuses to trim), `megahit_rename.py` (MEGAHIT headers to the
   SPAdes-style names the classification stages parse).
4. **`experiments/prod_run.sh`** — the ordering from `holdout_cascade_run.sh` with every
   sample-specific value required from the environment (`${X:?...}`) instead of defaulted. Nothing
   in it names a library.
5. **`experiments/measure_library.py`** (new) — the pair count, the modal read length per mate, the
   layout and the seed stride, measured from the FASTQ instead of asked of a person. Imo's question,
   on being shown the ledger above: *"why are we not discovering the facts from the sample itself?"*
   Three of the four things item 1 required are in the file. The registry stays what it is good for:
   an override with a basis, and the facts that genuinely are not in the file.

Nothing in any scoring stage was touched. Two things in the driver were, and both are bugs rather
than adaptations, so they are recorded as such.

## Bugs this pass found, which are bugs and not refusals

1. **The tier-2 stage hashed its input once per contig** (`hapax.py`, fixed in `50c2e97`). Every row
   of `contigs.jsonl` carries the digest of the fasta it was read from, and the digest was taken
   inside the per-contig writer — so the stage re-read its whole input once for each contig. On this
   library that is 194,703 contigs against 104 MB, about **20 TB of SHA-256** for a stage whose own
   arithmetic is a length and a GC fraction. Measured: it had run **1h45 and was roughly two thirds
   through**, and it is a per-seed cost, so five seeds would have spent most of a day on it. **The
   output is byte-identical either way**, which is why no run directory ever looked wrong — only a
   clock tells, and nobody had run this stage over a third of a million contigs before. The test
   counts digests rather than timing anything, because a count is deterministic;
   `tests/test_hapax_stages.py::test_the_source_digest_is_read_once_not_once_per_contig`, verified
   to fail on the old code.
2. **The stage list was typed in three places and the launcher's copy is the one that ran.** `panel`
   and `calibration` came out of production (SPEC 13.4) in `hapax.py` and in `prod_run.sh`, and
   `ec2_prod.py`'s own `--stages` default still carried `panel`. So `i-05594b6c0bc02e744` ran a validation stage
   in a production pass. It refused — correctly, with the committed reason, which is the design
   working — but the disagreement should not have been possible. The launcher now imports
   `hapax.PRODUCTION_STAGES`, and a test asserts the two cannot differ. Same family as a default
   where a refusal belongs: a literal in a second place.
3. **`hapax run --library X` told two of its thirteen stages which library they were running.**
   The other stage scripts take the library through `library_facts.add_argument`, whose default is
   the discriminating dataset so that every committed command reproduces unchanged — so a stage
   that was not told asked for **Berlin's** facts while reading this library's reads. It surfaced
   because `linkage` needs a block stride, was not given one, and refused: the read path named
   SRR23702504 while `--library` said SRR23973193, which is `resolve_args`' path check catching the
   hole one layer in. Under a read path that names no run — a work directory, which is the common
   case — there is no disagreement to catch and the stage would have used Berlin's stride in
   silence. `Run.call` now appends `--library` for every script that takes one, and the registry is
   checked against the scripts' own source rather than retyped
   (`tests/test_hapax_stages.py::test_every_library_aware_script_is_told_which_library`). Fixed in
   `e7ab28f`; this is the defect that made a fourth production launch necessary.
4. **`hapax.py` kept a private `git()` that returns `""`, so every seed manifest of this run
   records an empty commit and asserts a clean tree it never checked.** `provenance.py` was fixed
   for exactly this machine — an instance has a payload tarball and no `git`, so the launcher
   passes the commit in `HAPAX_GIT_SHA` — and the driver never imported it. It has its own helper
   that swallows the exception and returns the empty string, which `git_sha()` raises on. Two
   consequences, and the second is worse than the first. `git_sha` reads `""`, which a reader can
   at least see is nothing. But `tree_dirty` is `bool(git("status", "--porcelain"))`, and an empty
   string is `False`, so the manifest states positively that the tree was clean when no check ran.
   A refusal would have been a third launch; a silent `False` is a provenance claim nobody made.
   What pins this run's code is therefore `script_sha256` alone, identical across all five seeds —
   which is the field CLAUDE.md says earns the manifest its place, working as intended, but working
   alone. **The fix is deliberately not in the same commit as these results**: changing `hapax.py`
   changes its hash, and the five manifests pin `78dcb6cb…`, so a fix landing beside them would
   leave a directory whose pin matches no committed file. It lands in the commit after, where the
   history shows the results pinned to the code that produced them and the fix following.

## What the tool could not show while it ran

The assembly is the longest step in the pass and the only one with no visible progress.
`path_settings.py` captures MEGAHIT's stdout, so from outside the instance the log shows the
command and then nothing until a contig file exists. Over an hour-plus step that is the moment a
progress line is worth most, and its absence is indistinguishable from a hang. It is the same
family as the refusals above — an empty output and a step that never ran look identical — one
layer out: the run says nothing rather than saying it is working. Noted here rather than fixed,
because fixing it mid-run would have stamped a new script hash onto this run's output.

## The dead-man timer sizes the run, so five seeds needed two instances

The instance is held by `shutdown -h +N` and terminate-on-shutdown, which is what stops a forgotten
run from billing for days. It also cannot be moved once the instance is up: the timer is set in
UserData, this account's instances carry no SSM agent (`describe-instance-information` returns an
empty list for them), and there is no other way in. So the timer has to be right at launch, and
sizing it needs a per-seed wall clock that only the first seed of a new library can give.

Measured here: seed 1 of this library is **67 min** of stages after a **38 min** setup — fetching
the assembly from S3, trimming, classifying the contigs — which puts five seeds at about 5h45 past
a 5h30 timer. Rather than accept a run directory with a hole in it, or relaunch and pay for the
setup a fifth time, the remaining seeds were split across two instances working the same seed list
from opposite ends: `i-034b76bf8e2a55734` upward from seed 1, `i-01073d986644e28cc` downward from
seed 5. The two produce interchangeable seeds because everything that decides what a seed reads is
stated rather than discovered — the same commit, the same payload tarball, the same assembly object
in S3, and `block = 6,000,000` from `libraries.json` — and each seed's own manifest carries the
commit, the script hashes and the seed, so which instance ran it is recorded and not load-bearing.

Two notes fall out of that arithmetic rather than being designed:

- **The launcher's stride is advisory and the registry's is the one that runs.** `ec2_prod.py`
  preflights with the `block` in `libraries.json` (6,000,000) and hands it to the instance;
  `measure_library.py` independently derives 7,100,000 from the 32,125,614 pairs it counts in the
  file. Both fit five disjoint blocks, the stated value wins because stating it is the override,
  and the measurement ships beside the run unused. That is the intended behaviour, but the block in
  the registry is a value **I supplied**, and the file could have supplied it: it belongs in the
  ledger above as hand-supplied, not as measured.
- Per-seed cost on this instance type is about **$0.56**, so a seed is cheap and the setup is not.
  `--assembly-from` is what makes a relaunch cost one seed rather than a whole pass.

## What the generator could not build

`experiments/report_prod.py` turns a run directory into `TABLES.md`, and the point of the
convention that a figure is generated from committed files is that prose cannot drift from the
outputs. Run against a directory nobody had run it against, it produced a page that rendered
cleanly and was missing most of the run. Five things, four of them in the generator:

1. **The cascade curve was a column of dashes, on this library and on the holdout.** The per-tier
   removal a reader wants — alone, marginal, cumulative, and the index size that bought it — is one
   row per tier in `tiers.json`'s `curve`, with the sizes in a top-level `index_bytes` dict. The
   generator read those four fields out of `tiers["tiers"][name]`, which holds each tier's internal
   sweep instead. **A table of dashes renders as cleanly as a table of numbers**, which is this
   repo's recurring failure one layer out: an empty cell and a measurement that never happened look
   identical. Fixed, and the table it now writes is the most informative thing in the file — on this
   water the sample's own assembly removes 71.9% alone where the 16 GB Kraken2 database removes
   0.65%.
2. **A tier that is measured and never subtracted vanished from the table rather than being named.**
   `observe_sample_contigs_unidentified` has no row in the curve by design, because SPEC's
   subtrahend is identification and matching an unnamed contig does not remove the read. Its absence
   now carries a sentence, since subtracting that tier is exactly the error that inflated removal to
   85.9% against the true 77.0%.
3. **The assembly section was empty, as a consequence of `--assembly-from`.** The assembly ran on
   one instance and the seeds on another, so its manifest was not in the directory the generator
   reads. The resume path that makes a relaunch cost one seed instead of a whole pass also splits
   the run's provenance across prefixes, and nothing noticed. Fetched by hand into the committed
   directory; the generator still has no way to follow a run that spans instances.
4. **The file's own measurements were not in the report at all.** `measured_library.json` is
   written beside every run and nothing read it, so the one comparison a reader of a new library
   wants — what the FASTQ says against what the registry states — could not be made. It is now a
   table, and it shows a disagreement worth seeing: the registry's stride is 6,000,000, the file's
   own arithmetic gives 7,100,000, and the stated one is what ran.
5. **Fractions printed to sixteen digits** at n = 323,430, and five seed-level refusals printed as
   five near-identical sentences. Both are cosmetic and both are the kind of thing that makes a
   generated page read as a dump rather than as a result.

A sixth belongs here although it is mine rather than the generator's. This run's unidentifiable
pile is 194,703 contigs, so `contigs.jsonl` is 79 MB a seed and 395 MB over five, in a repo whose
whole `results/` tree is 583 MB. Committing it compressed is the right call; compressing
`reads.jsonl` and `layers.jsonl` with it was not. **This was the first run directory in the repo to
carry a `.jsonl.gz`**, the page builder that reads a run directory looks for the plain name, and it
drew an empty residue table and two empty search arms for a 32 M pair library without saying
anything. Both spellings are read now, and a directory carrying both refuses rather than choosing;
here the two small files are committed plain, under the name every other run directory uses, and
only the 79 MB one is compressed. A reader that cannot find a file nobody told it to expect reports
a zero, which is this repo's oldest failure mode wearing a new hat.

None of these was a wrong number. All five are a generated page that looked finished and was not,
which is the failure mode the convention exists to prevent, so they belong here rather than in a
commit message.

## Cost and wall clock

| Launch | Commit | Result | Wall | Cost |
|---|---|---|---|---|
| i-0c35504358c8f69e3 | ed01fc3 | `pip install` exited 0, `import mappy` failed under `-I`; `--quiet` hid the reason | ~3 min | ~$0.03 |
| i-0d20d59fe31946d78 | c76d6a5 | Kraken2 smoke test piped through `head` under `pipefail`: SIGPIPE read as a broken binary | ~90 s | ~$0.02 |
| i-00007d94fdc55f541 | bb0f48b | provenance raised with no `git`, and the shell reported a refusal that had not happened — plus the real prefix measurement | ~13 min | ~$0.07 |
| i-0e44a2d889f41d119 | 90060fa | the assembly: 2h07 of MEGAHIT over 32.1 M pairs, kept in S3 and reused by every launch after it | ~2h53 | ~$1.45 |
| i-05594b6c0bc02e744 | 0bd82f8 | terminated: the tier-2 stage was 1h45 into a per-seed hash that should take seconds (bug 1) | ~3h35 | ~$1.80 |
| i-0b62260310974e946 | 50c2e97 | terminated: `linkage` refused because the driver never passed `--library` (bug 3) | ~1h45 | ~$0.90 |
| i-034b76bf8e2a55734 | e7ab28f | the production pass, seeds 1 upward: seeds 1, 2 and 3, then terminated mid-seed-4 once the other instance had it | 4h35 | $2.31 |
| i-01073d986644e28cc | e7ab28f | the same pass, seeds 5 downward: seeds 5 and 4, which completed the set, then terminated | 2h47 | $1.41 |

The two production instances cost **$3.72** together, against the $4.71 their timers would have
billed if both had run to the end — the saving is only that seed 4 landed on one of them while the
other was still working on it, so that one could be killed early. Running the same five seeds on one
instance would have cost about **$3.15**: the duplication is the price of a timer that cannot be
moved, and it is about 18% of the pass.

The three dead staging launches cost about **$0.12** in total. Every one of them was a failure in
staging, not in the pipeline, and each is now a test or an assertion. The two terminated production
launches cost **$2.70** and each bought one defect that only a library nobody had measured could
show: a stage whose cost is quadratic in contig count, and a driver that handed the wrong library's
facts to eleven of its thirteen stages.

**The whole day is $7.99**, of which the pass that produced the committed result is $3.72 and the
assembly it rests on is $1.45. The remaining $2.82 bought the three bugs above and three staging
failures, so a second library should cost about **$5.20** on this instance type: the assembly plus
five seeds, with no further staging. All of it sits inside the $10 pre-approved for a single run.

## Provenance

- Launch manifests: `launch_i-*.json`, one per launch, keyed by the commit the payload was built
  from. The launcher refuses to ship a dirty tree, so the commit identifies the code.
- `pip_freeze.txt`, `tool_versions.txt`: the environment the pass ran in.
- `nontemplate_prefix/nontemplate_prefix_s{1..5}.json`: the trim measurement, which refused.
- Tables in `TABLES.md` are generated by `experiments/report_prod.py` from the committed run
  directories. No number in this file is typed by hand that is not in those files.
