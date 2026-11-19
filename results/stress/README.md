# Can the incumbents be broken?

Imo's redirect (2026-11-19): before building more of hapax, try to break the
published viral haplotype assemblers on the hardest realistic case. If they
hold up, hapax is not needed. If they break, that is the justification and the
test set.

This directory holds that attempt. Two arms, both scored with the same rule
CLAUDE.md imposes on hapax itself: every contiguity number ships beside a
chimaera number from the same run.

## Arm 1: within-species resolution (`python/bench_panel.py`)

Truth is real sequence, not simulated evolution: published single-genome
amplification env sequence sets from the Zambia transmission cohort, one set
per host, so the panel is a sample of haplotypes that genuinely co-existed in
one person. Length variation in the hypervariable loops comes with the data
(2412-2604 nt for the same gene within one host) rather than being imposed.
Reads are simulated on top, since that part has to be controlled.

Panels record their own minimum pairwise divergence, because no assembler can
be asked to separate haplotypes closer than its error model, and identical
sequences are dropped: a panel holding the same sequence twice makes recovery
unreachable by construction.

Scored by `python/score_haplotypes.py`: haplotypes recovered at 99% identity
over 95% of their length, haplotype switches within a contig, and frameshifts
per kb against the best-matching truth haplotype. Frameshifts are the number
the codon-space argument stands on: the truth haplotypes are real coding
sequences whose own indels are codon multiples, so a non-multiple-of-three
indel in a contig is the assembler's.

## Arm 2: novelty in a mixed community (`python/mock_community.py`)

Three bacterial genomes as background so the viruses are a minority of the
reads, phages for their mosaic structure, thirteen known viruses, and seven
of those paired with a divergent relative also present in the sample:

| pair | divergence |
|---|---|
| SARS-CoV-2 / RaTG13 | 0.040 |
| HBV genotype D / A | 0.105 |
| HIV-1 subtype B / C | 0.163 |
| SARS-CoV-2 / bat CoV BM48-31 | 0.240 |
| phage lambda / HK97 | 0.458 |
| HPV16 / HPV11 | 0.506 |
| HHV-5 / HHV-7 | 0.546 |

The env swarm goes into the same sample, so between-species novelty and
within-species swarm are present at once.

Scored by `python/score_community.py`. The headline is not contig count but
**collapse**: whether the divergent relative comes back as its own genome, or
its reads were absorbed into the known neighbour. A relative that collapses is
a novel virus that the pipeline would never have reported. Against that sits
**false joins**: a contig whose windows match two different genomes, which is
how a novel virus gets invented rather than missed.

## Results

### The community, Illumina and nanopore (seed 1)

Panel: 33 truth sequences, 23 scored genomes (three bacterial background
chromosomes plus twenty viral) and ten env swarm haplotypes. Wall clock and
peak memory on 4 cores: MEGAHIT 1228 s / 142 MB, metaSPAdes 708 s / 1078 MB,
Haploflow 203 s, metaFlye 125 s / 303 MB (nanopore, 27 Mbp at ~8x mean).

See `table.md`. In short: nobody collapsed more than one relative, nobody
recovered a single swarm haplotype, and the two tools with the cleanest
novelty scores are the two plain metagenome assemblers.

- **MEGAHIT** 18 of 23 genomes, 5 of 7 relatives separate, no collapse, no
  real cross-genome chimaera. Its three chimaeras are all joins between env
  swarm haplotypes.
- **metaSPAdes** deleted RaTG13 outright: zero contigs assigned to it, 0%
  covered, while SARS-CoV-2 came back complete. It also emitted a 13.8 kb
  contig that runs HIV-1 subtype B from 1275 to 9275, wraps through B's
  start, then continues as subtype C, every 500 nt window matching its part
  at 100% identity. A fabricated inter-subtype recombinant.
- **Haploflow** 1 of 7 relatives, bat CoV BM48-31 collapsed, and 7
  cross-genome chimaeras. It is built for two-strain mixtures and this is not
  one.
- **metaFlye** on nanopore collapsed nothing but joined bacteria to viruses:
  E. coli to HHV-5, and E. coli to bat CoV BM48-31 over 36 kb.

One join is not an error and is excluded: phage lambda and E. coli K-12 share
3.5 kb at 98% identity, so a contig running between them is correct sequence.

### Where the collapse threshold sits

`ladder.md`: the same two-member mixture at 100x and 25x, MEGAHIT, across real
sequence pairs. The minor member is **erased at or below 2.5% divergence**
(zero contigs, zero coverage, and at the two closest rungs the single output
contig is a mosaic of both), half recovered at 4%, and cleanly separate from
10.5% upwards.

That is the number this exercise was for. Between-species novelty is solved:
a relative 10% away comes back as its own genome from a 2010s de Bruijn
assembler. Below about 4%, a virus that is really there disappears into its
neighbour. Intra-host swarms, and the novel viruses that are near-relatives of
known ones, live entirely inside that blind spot.

### The haplotype panels

`within_*`, `ctrl*` and `p*` JSONs: on real within-host env sets, no
short-read tool recovered more than 2 of 20 co-existing haplotypes, and on
5-haplotype panels at even abundance no tool exceeded 2 of 5. These runs are
single-seed and the tools are built for whole genomes rather than a 2.5 kb
amplicon, so they are a weaker claim than the community results and are kept
as supporting evidence, not a headline.
