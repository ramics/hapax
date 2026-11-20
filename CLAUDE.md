# Operating brief

Read [`SPEC.md`](SPEC.md) first — it is the production tool spec and the source of truth for
intent. Then [`results/README.md`](results/README.md), the experiment archive, which says in
one line each what has been measured and what has been killed. This file covers how to work
in the repo.

**Current direction (Imo, 2026-11-20): the production tool in `SPEC.md`** — subtract
everything identifiable, check what assembled but could not be identified, characterise the
residual reads, and report the whole thing as a multidimensional measurement of signal
against noise rather than a novel-or-not call.

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

**Measure chimaeras before claiming contiguity.** Longer contigs that are wrong are worse
than short contigs that are right. Every contiguity number ships beside a chimaera number
from the same run — and beside that test's own false-positive rate, which for split-mapping
is 26% on known background.

**Negative results ship.** Most of what this repo knows, it knows from things that failed.
Record them and say so rather than tuning a simulator until it looks good.

## Conventions

- Every experiment is a script in `experiments/`, takes `--seed`, and writes a run manifest
  (seed, git SHA, parameters, timings) beside its results.
- Any reported number is reproducible from a committed script plus a seed.
- Headline results run across at least 5 seeds, reporting spread. A single-seed conclusion
  from simulated data is noise.
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
