# hapax

Take a wastewater sample, return the full breakdown: what is known and by what, what
assembled, what did not, what codes and what does not — each pile quantified, and each
answer priced.

**Start here: [`SPEC.md`](SPEC.md).** It is the production tool spec and the source of truth
for intent. It is short by design and every number in it is measured on real data with a
committed script and a seed behind it.

Then: [`CLAUDE.md`](CLAUDE.md) for how to work in the repo, and
**[`results/README.md`](results/README.md) for the experiment archive** — one line per
experiment saying what it measured, what it killed, and where its report is.

## What this is not

It does not emit a novel-or-not call. The target is a multidimensional measurement of signal
against noise, because every feature this project has tried as a novelty score has been
measured to fail as a threshold on one axis (`SPEC.md` §1.1). The declared blind spots are
part of the deliverable: no read-level rule approaches one false positive per million,
recognition of a known virus falls away between 10% and 20% divergence, and detection below
about one in a million is not available at all (`SPEC.md` §7).

## History

Two earlier directions are archived rather than deleted, because what they killed is the
reason the current one is shaped as it is.

- **The codon-aware assembler**, the first attempt, is under
  [`experiments/codon_assembler/`](experiments/codon_assembler/) with its own spec and all
  its results, laid out as the repo root once was so its scripts still run from there. It
  ships as a negative: scoring overlaps in codon space gave a small gain over nucleotide
  scoring that showed only at high simulated diversity, and on real multi-host swarms it
  lost to metaSPAdes and MEGAHIT by merging sister swarms.
- **Strain separation during assembly** produced the measured benchmark in
  [`results/pair_ladder/`](results/pair_ladder/REPORT.md): below about 3.8% divergence every
  assembler fuses lineages at any depth, and the reads carry enough information to separate
  them, so what is missing is path selection rather than signal. That result is live and is
  cited by the spec.
