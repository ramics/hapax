# Operating brief

This file covers how to work in the repo.

**Current direction (Imo, 2026-11-20): separating closely related strains during
assembly.** The codon-aware assembler that came first is archived, verbatim and still
runnable, under `experiments/codon_assembler/` (its spec is
`experiments/codon_assembler/METHOD.md`, its findings
`experiments/codon_assembler/results/`). Read those before starting anything new: they
record what has been measured and what was killed by measurement. A spec for the new
direction has not been written yet.

## What this was

The archived first design: a codon-aware assembler for viral populations. Overlaps are detected in codon space
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

**Negative results ship.** The kill criteria in the archived METHOD.md §5 are real. If the
synonymous collapse buys nothing over MEGAHIT at realistic diversity, record it and
say so rather than tuning the simulator until it looks good.

## Order of work (the archived codon-assembler plan)

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
compiles. The archived crate (`experiments/codon_assembler/Cargo.toml`) has rayon and a release
profile; new core code may start a fresh crate or lift modules out of it. Add PyO3 bindings when there is something worth binding, not before.

**The bigger lever is indexing, not the language.** All-against-all overlap is
quadratic in any language. minimap2 is fast because of minimizer indexing, not
because it is C. Build the seeding index first — minimizers, sketching, bucketing —
and treat the Rust rewrite as the constant factor on top. A fast quadratic assembler
is still a quadratic assembler.

## Repo layout

```
CLAUDE.md            this file
data/MANIFEST.md     accessions and download dates (data/ itself is gitignored)
experiments/
  codon_assembler/   the first assembler and every experiment run on it, moved here
                     unchanged with its results so provenance survives. Run its
                     scripts from this directory: it is laid out as the repo root was
                     (Cargo.toml, src/, tests/, python/, experiments/, results/), and
                     its data/ is a symlink to the top-level data/.
```

`simulate.py` (now `experiments/codon_assembler/python/simulate.py`) is a working
generator for coding sequences under substitution and indel pressure, with reading
frames preserved. Its Markov surrogate sampler is the fiddly part and it works — keep
it under test if you refactor or lift it out. It stays Python.

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
