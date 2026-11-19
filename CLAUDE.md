# Operating brief

Read `METHOD.md` first. This file covers how to work in the repo.

## What this is

A codon-aware assembler for viral populations. Overlaps are detected in codon space
so synonymous variation collapses; nucleotides are kept throughout so linkage and
gene order survive.

## Non-negotiables

**Nucleotides are never discarded.** Codon space is where overlaps are *scored*.
The sequence carried, laid out and emitted is nucleotide. Any design that drops to
protein-only representation has lost linkage and is wrong for this tool.

**Reads are joined only when they share nucleotide identity.** Divergence tolerance
is for comparing a read to a reference, which this tool never does. A permissive
overlap criterion is a chimaera generator.

**Measure chimaeras before claiming contiguity.** Longer contigs that are wrong are
worse than short contigs that are right. Every contiguity number ships beside a
chimaera number from the same run.

**Negative results ship.** The kill criteria in METHOD.md §5 are real. If the
synonymous collapse buys nothing over MEGAHIT at realistic diversity, record it and
say so rather than tuning the simulator until it looks good.

## Order of work

1. **Pairwise overlap scoring.** The core primitive: two reads, three candidate
   offsets, best score wins. No global frame is elected — viral genomes read in
   several frames at once, and assembly only needs the offset between a pair.
2. **Two-tier dispatch.** Nucleotide seeding, then the codon HMM only where seeding
   is ambiguous. Instrument what fraction falls through to the HMM — if it is most
   of them, the design has failed and you should find out early.
3. **Per-region scoring mode.** Codon scoring where coding, nucleotide scoring
   through untranslated regions and dual-coding regions.
4. **Layout and variant-aware consensus.**
5. **The EvolveAGene diversity ladder**, against MEGAHIT, metaSPAdes and PLASS.
6. Real data: HIV *env*, then *gag/pol* as the dual-coding stress test, then public
   wastewater.

## Language

**Rust for the core, Python for experiments.** This is an assembler; the hot path
compiles. `Cargo.toml` is scaffolded with rayon and a release profile; `src/` is
yours to create. Add PyO3 bindings when there is something worth binding, not before.

**The bigger lever is indexing, not the language.** All-against-all overlap is
quadratic in any language. minimap2 is fast because of minimizer indexing, not
because it is C. Build the seeding index first — minimizers, sketching, bucketing —
and treat the Rust rewrite as the constant factor on top. A fast quadratic assembler
is still a quadratic assembler.

## Repo layout

```
METHOD.md            build spec; source of truth for intent
CLAUDE.md            this file
Cargo.toml           scaffolded; src/ to be created
src/                 Rust core: index, overlap, codon scoring, layout
python/
  simulate.py        sequence evolution + surrogates  [seeded, repurpose]
experiments/
results/             committed: tables, figures, run manifests
tests/
```

`simulate.py` is a working generator for coding sequences under substitution and
indel pressure, with reading frames preserved. It was written for a different
experiment; repurpose it to simulate within-population swarms rather than separate
families. Its Markov surrogate sampler is the fiddly part and it works — keep it
under test if you refactor. It stays Python; it is nowhere near the hot loop.

## Conventions

- Every experiment is a script in `experiments/`, takes `--seed`, writes a run
  manifest (seed, git SHA, parameters, timings) beside its results.
- Benchmark numbers come from release builds, never debug.
- Any reported number is reproducible from a committed script plus a seed.
- Headline results run across at least 5 seeds, reporting spread. A single-seed
  conclusion from simulated data is noise.
- Docstrings explain *why* a step exists and what it is calibrated against. That
  reasoning is the research.

## Data

`data/` is gitignored and fetched locally; this sandbox cannot reach NCBI. Record
accessions and download dates in `data/MANIFEST.md`.
