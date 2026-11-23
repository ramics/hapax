# Operating brief

Read [`SPEC.md`](SPEC.md) first — it is the production tool spec and the source of truth for
intent. Then [`results/README.md`](results/README.md), the experiment archive, which says in
one line each what has been measured and what has been killed. This file covers how to work
in the repo.

**Current direction (Imo, 2026-11-20, recut 2026-11-23): the production system in
`SPEC.md`** — subtract everything identifiable, check what assembled but could not be
identified, characterise the residual reads, and report the whole thing as a
multidimensional measurement of signal against noise rather than a novel-or-not call.

The spec is now organised around **three output classes on one distance continuum** —
divergent-but-related contigs, novel strains of known viruses, and the unplaceable residue —
each with what it reports and what it refuses to claim, plus three sections that are the
spec's working memory: **what is measured** (§10), **what is dead and must not be
re-proposed** (§11), and **what is open, with the experiment that would settle it** (§12).
Read §11 before proposing anything. The goal the whole thing serves, in Imo's words, is
flagging a novel virus that might be pointed at humans, and the combined signal is **"new
AND like X"**: distance without a nearest-relative identity is half an answer.

Two earlier directions are archived rather than deleted, because what they killed is why the
current one is shaped as it is: the codon-aware assembler under
`experiments/codon_assembler/` (its spec is `experiments/codon_assembler/METHOD.md`), and
strain separation during assembly, whose live result is `results/pair_ladder/REPORT.md`.
Read the archive index before starting anything new.

## Non-negotiables

**Searching only ever removes background. It is never used to score novelty.** An index can
answer *is this already known*; it cannot answer *is this new*, because score ordering tracks
database proximity. Measured, repeatedly. Any design that scores novelty with a search has
reintroduced the failure this project is built around.

**Nucleotides are never discarded.** Codon or protein space is where sequence is *scored*.
Anything carried, laid out and emitted is nucleotide. A protein-only representation has lost
linkage.

**Every number ships with its operating point.** Thresholds dominate database size by an
order of magnitude more than size does — one Kraken2 database removes 42.8%, 21.7% or 3.5%
of the same library depending only on its confidence setting. A reported fraction without its
threshold is unreadable, including the ones this project has quoted at itself.

**Every field ships with its calibration.** A shuffled null says what the field reads on
sequence carrying no information; a known-positive comparison against subtracted reads says
whether it carries anything about the residue. A field that reads the same on both is
reported with that fact, not hidden.

**A shuffled null comes from `experiments/shuffle_control.py`, imported, never
reimplemented, and it raises rather than falling back.** Two scripts each carried their own
naive dinucleotide shuffle, and both handled a dead-ended walk by retrying and then
returning the input unchanged. That fired for 324 of 11,361 contigs, 2.85%, and all 304
"shuffled" contigs that produced a significant protein hit were among them: zero genuinely
shuffled contigs produced one. The bad null published an existence ratio, a claim that
identity cannot separate real sequence from noise, and two kills of other people's fields,
all withdrawn. It was the same bug twice in two days, already fixed once in PR #29, and the
only arm that was right imported its shuffle instead of writing one. Two consequences
beyond the module: **when two instruments disagree about a CONTROL, chase that before any
result** - 0.00% against 2.68% on the same nominal null is what exposed this - and **record
"no significant hits" as an explicit census**, because an empty table and a search that
never ran look identical.

**A null is not a bias statement, and a number needs both.** A null says what a field reads
on meaningless input. A bias statement says what it reads on *real* input for the wrong
reason — and only the second can be removed by a control. So every number ships with the
direction of its bias and the control that removes it, or, where no control exists, with the
direction alone. Three measurements on one day were publishable-looking and wrong until
their bias was named: a nearest-neighbour test whose ties were broken by array order
inflated its own control; a pair-symmetric field that averaged two opposite leans caused by
the mates' unequal length; and a chimaera test that read 100% on both piles because each
query sat inside the reference it was searched against. The tell in all three was the
control moving with the test, so read the control first.

**A provenance change forces a re-run, and that is how reproduction gets done here.** The
convention below gives reproducibility as a property of the design — a committed script plus
a seed. Having actually re-run a result from clean and got the same numbers is a stronger
claim, and almost nothing in this repo can make it. Nobody schedules a reproduction
campaign, so the rule is the cheap version: any change to how provenance is recorded means
re-running what it describes, and the reproduction falls out free. It is how
`results/residue_structure/` came to carry one (§6c) — not by intending to.

**Measure chimaeras before claiming contiguity.** Longer contigs that are wrong are worse
than short contigs that are right. Every contiguity number ships beside a chimaera number
from the same run — and beside that test's own false-positive rate, which for split-mapping
is 26% on known background.

**Negative results ship.** Most of what this repo knows, it knows from things that failed.
Record them and say so rather than tuning a simulator until it looks good.

## Conventions

- Every experiment is a script in `experiments/`, takes `--seed`, and writes a run manifest
  (seed, git SHA, parameters, timings) beside its results. **The git SHA is the field that
  earns the manifest its place**: a fix landing midway through a multi-seed run leaves a
  directory where some seeds are right and some are wrong and nothing visible says which. A
  mixture is not a five-seed result — re-run every seed after a fix, and until they are all
  re-run, read the per-seed manifest rather than the directory. **A `-dirty` SHA is
  unresolved provenance, not provenance**: it names HEAD, not the code that ran, so two runs
  can carry the same SHA and different code, which is exactly the case this rule comes from.
  Every manifest therefore carries a content hash of the script beside the SHA
  (`script_sha256`), which pins the code whatever the tree looked like. File timestamps
  answer the same question only while the container lives, which is not long enough.
  **The hash is taken when the manifest is written, not when the run starts**, so editing a
  script while a run of it is in flight stamps the new hash onto the old code's output and
  destroys exactly what the hash exists to protect. Finish the run or kill it; this has now
  cost two restarts.
- Any reported number is reproducible from a committed script plus a seed.
- **Controls, search-output parsing and provenance come from shared modules, imported,
  never re-typed**: `shuffle_control.py`, `search_hits.py`, `provenance.py`. Each raises
  where a copy used to fall back. The same shuffle bug shipped twice in two days because
  the second writer did not know which condition the first was fixing. Existing copies
  are grandfathered, with their measured behaviour, in `experiments/SHARED.md`. A tripwire
  test fails on a new shuffle.
- **Five seeds give an honest range and a useless variance.** The min-max this repo reports
  is an observed spread, not an estimate of dispersion: the sampling error of a five-sample
  standard deviation swamps the effect it would be measuring. Measured, and it reverses the
  sign - the inflation that reads 1.44 and 1.61 over the whole library's 87 windows reads
  0.65 and 0.44 off the same five seeds, which would have been reported as *narrower* than
  random sampling (`results/block_position_bias/`). Any correction factor comes from the
  full series, never from the run's own seeds.
- Headline results run across at least 5 seeds, reporting spread. A single-seed conclusion
  from simulated data is noise.
- **The two mates are not interchangeable, and any field computed symmetrically over a pair
  inherits that.** Mate 2 of SRR23973193 carries ~14 nt of non-template 5' sequence which
  every experiment trims, so R2 is 93 nt against R1's 107 — 62 scorable 31-mers against 76 —
  and it maps at a different rate. Measured consequences so far: asymmetric absence leans to
  R2 (64:36) while a fixed-width error signature leans to R1 (58:42), opposite directions
  from one cause (`results/residue_structure/REPORT.md` §6b); mate novelty and mapping rates
  differ likewise (`results/residue_structure/`, `results/subtract_known/`). Report a
  pair-symmetric field per mate, or show the length-matched control. A pooled number here
  averages two biases that point opposite ways and cannot be unpicked afterwards.
- Percentages of a library are of read pairs, the occurrence-weighted unit. Never average
  over distinct k-mers: for conserved sequence the two differ by two orders of magnitude, and
  a read sees the occurrence-weighted one.
- Figures are generated from committed result files, never hand-written. A number typed into
  a figure is a number without provenance, and this repo has had to withdraw claims whose
  prose drifted from its outputs.
- Benchmark numbers come from release builds, never debug.
- Docstrings explain *why* a step exists and what it is calibrated against. That reasoning is
  the research.

## Language

**Rust for a hot path that exists, Python for experiments.** Nothing in the current direction
is compute-bound in Python yet; the indexes and aligners it drives are already compiled. The
archived crate (`experiments/codon_assembler/Cargo.toml`) has rayon and a release profile and
modules can be lifted out of it. Add PyO3 bindings when there is something worth binding.

**The bigger lever is indexing, not the language.** The measured lesson from the cascade is
sharper than that: 1.4 MB of this sample's own contigs removed more of the library than 272
MB of SILVA, and 15 MB of them removed 30 points where 16 GB of Kraken2 removed 0.4. Which
index you consult beats how fast you consult it, by orders of magnitude.

## Repo layout

```
SPEC.md              the production tool; source of truth for intent
CLAUDE.md            this file
README.md            front door
results/README.md    experiment archive index — read this before adding an experiment
requirements.txt     Python dependencies; pinned where a manifest records the version
TOOLS.md             external binaries and their pins (foldseek, DIAMOND, BLAST, ...)
experiments/*.py     current experiments, one per results/ directory
experiments/codon_assembler/
                     the first direction, moved here unchanged with its results so
                     provenance survives. Run its scripts from that directory: it is
                     laid out as the repo root once was, and its data/ symlinks to the
                     top-level data/.
data/MANIFEST.md     accessions and download dates (data/ itself is gitignored)
```

## Data

`data/` is gitignored and fetched locally; this sandbox cannot reach NCBI. Record accessions
and download dates in `data/MANIFEST.md`.

The discriminating dataset is Berlin wastewater SRR23973193. The New York cryptic-lineage
runs are a correctness fixture only.

## Sample facts are data, not defaults

**A value that is a property of one library never has a default.** A trim length, a block
size and a read path are facts about a sequencing run, and they were written into roughly
twenty scripts as `default=14`, into five more as a bare `14` inside the call, and into two
as a module-level `TRIM_R2`. On a library with no non-template prefix that deletes fourteen
real bases from every mate, silently, in the direction of less recognition; a block sized
to Berlin's 8.73 M pairs leaves seeds 4 and 5 of a 4.29 M-read library reading nothing, and
the directory still looks like a five-seed result. They live in
`experiments/libraries.json`, keyed by accession, each with how it was established, and are
read through `experiments/library_facts.py`, which **raises for a library it does not know
rather than falling back**. Same rule as the shuffle: the fallback is the bug.

**Resolving the value is only half of it: `--library` has to be checked against the
reads.** It defaults to the discriminating dataset so that every committed command
reproduces unchanged, which means a caller who points `--r1` at another library's FASTQ and
says nothing would still have got Berlin's trim — the same silence, one layer in.
`resolve_args` therefore refuses when a read path names a *different* run, and a path
naming no run at all (`data/berlin`, a work directory) cannot disagree. Stating the value
on the command line always wins and is never checked, because stating it is the override.
A second operating point of the same quantity gets its own key rather than a reconciliation:
`trim_r2_for_mapback` is 16 because that is what `results/divergent_known/` was produced
with, and the entry says so. `tests/test_library_facts.py` asserts the values rather than
assuming them, the negatives included.

**The same applies to words.** Every plain-English claim the web view makes about biology
is in `experiments/virus_vocabulary.json` with its basis, applied by
`experiments/architecture_call.py`, which records the evidence that fired and declines
with a reason rather than ever recording "not a phage". The renderer assembles sentences
and has no vocabulary of its own. This does not make the English stop being hand-written;
it makes it reviewable, attributable, and wrong in exactly one place. The convention that
a figure is generated from committed files was being applied to the numbers and not to
the claims beside them, and the claims are usually larger.
