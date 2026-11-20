# hapax

Assembly of viral populations that keeps closely related strains apart.

The current direction is separating near-relative lineages (roughly 1-10% divergent)
during assembly, where incumbent assemblers either collapse the minor lineage or fuse
the two into mosaic contigs. No code for it exists yet.

The first attempt, a codon-aware assembler, is archived with all its experiments and
results under [`experiments/codon_assembler/`](experiments/codon_assembler/). It ships
as a negative: scoring overlaps in codon space gave a small gain over nucleotide
scoring that showed only at high simulated diversity, and on real multi-host swarms
the assembler lost to metaSPAdes and MEGAHIT by merging sister swarms. Its results tables are the measured baseline this
direction starts from. See `CLAUDE.md` for how to work in the repo.
