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

(filled in by the runs; each table carries the panel's divergence and the
tool's wall clock and peak memory from the same run)
