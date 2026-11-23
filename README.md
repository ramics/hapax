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

## Running it

```
python3 hapax.py run SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz --seed 1 --out RUNDIR \
    --trim-r2 14 \
    --contigs IDENTIFIED_RRNA.fa IDENTIFIED_OTHER.fa --unidentified-contigs UNIDENTIFIED.fa \
    --identified-contigs IDENTIFIED_OTHER.fa --rrna-contigs IDENTIFIED_RRNA.fa \
    --assembly ASSEMBLY.fa --cross-sample-ref OTHER_SAMPLES.fa --spacers SPACERS.fa
python3 experiments/make_web.py results/web/pages/SAMPLE.json   # the run viewer
```

The web view is built by one command for every sample. What a page says about its sample
that is not in the run directory - its title, one sentence on what the water is, and every
result file outside the run that a panel reads, each as a path or as `not_run` with a
reason - is in that sample's descriptor under `results/web/pages/`. Methods entries live in
`results/web/methods.json`, with any row measured on one library stored under its
accession; on another sample's page that row reads "Not measured on this library".

Everything that belongs to one sample - its trim, its own contigs, its assembly - has to
be given: there are no defaults for those, because a default is some other sample's file
and the run would use it without complaint. A per-read field whose reference you do not
have is left out by name (`--report-without cross_sample crispr ...`) and recorded as not
searched. The run checks every input, index and tool before it starts and refuses, exit
status 2, if anything is missing.

`hapax.py` is an orchestrator: every number in `RUNDIR` comes from the `experiments/`
script that already measures it, with that script's own seed and manifest. The web view is
a thin reader over `RUNDIR` and contains no number that is not in a committed file. Add
`--plan` to print the stages without running them.

## How far signal-to-noise gets today, and where it stops

**It gets to:** a three-way breakdown of the library with intervals and operating points
(77.0% identified, 11.9% assembled-but-unidentified, 11.1% residual reads on this sample);
per read, the distance to its closest relative with the alignment's coverage and E-value,
what kind of sequence it is, whether it recurs in the library, and what it physically sits
next to; per subject, how many residue reads point at one known relative from how far,
which is where a new strain of a known virus appears; and a divergence ladder from a spiked
holdout genome that says what the whole cascade would have missed at 10% and 20%.

**It stops at:** any per-read verdict, and grouping. A mate finds its own partner in its ten
nearest neighbours 23 times above chance with no shared sequence, which says the geometry
carries real molecule-level information — but the rate falls as the space grows, so at
library scale a read cannot find even its own mate, let alone a stranger from the same
genome ([`results/mate_friends/`](results/mate_friends/REPORT.md)). Nothing here approaches
one false positive per million.

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
