# Can the incumbents be broken?

Imo's redirect (2026-11-19): before building more of hapax, try to break the
published viral haplotype assemblers on the hardest realistic case. If they
hold up, hapax is not needed. If they break, that is the justification and the
test set.

This directory holds that attempt. Two arms, both scored with the same rule
CLAUDE.md imposes on hapax itself: every contiguity number ships beside a
chimaera number from the same run.

## Arm 1: within-species resolution (`experiments/codon_assembler/python/bench_panel.py`)

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

Scored by `experiments/codon_assembler/python/score_haplotypes.py`: haplotypes recovered at 99% identity
over 95% of their length, haplotype switches within a contig, and frameshifts
per kb against the best-matching truth haplotype. Frameshifts are the number
the codon-space argument stands on: the truth haplotypes are real coding
sequences whose own indels are codon multiples, so a non-multiple-of-three
indel in a contig is the assembler's.

## Arm 2: novelty in a mixed community (`experiments/codon_assembler/python/mock_community.py`)

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

Scored by `experiments/codon_assembler/python/score_community.py`. The headline is not contig count but
**collapse**: whether the divergent relative comes back as its own genome, or
its reads were absorbed into the known neighbour. A relative that collapses is
a novel virus that the pipeline would never have reported. Against that sits
**false joins**: a contig whose windows match two different genomes, which is
how a novel virus gets invented rather than missed.

## Results

Two community seeds (`meta1`, `meta2`), four assemblers plus viralFlye, and a
twelve-rung divergence ladder. Tables in `table.md` and `ladder.md`, per-run
JSON beside them. Wall clock on 4 cores: MEGAHIT 615-1228 s, metaSPAdes
419-708 s, Haploflow 182-203 s, metaFlye 95-125 s on 27-55 Mbp of nanopore,
viralFlye 261 s on top of metaFlye.

### Untargeted novelty is not the gap

A relative 10% or more away from its known neighbour comes back as its own
complete genome. MEGAHIT recovered 5 of 7 relatives in both seeds with
**nothing collapsed and no cross-genome chimaera**, including a bat
coronavirus 24% from SARS-CoV-2 and HBV genotype A 10.5% from genotype D.
If the claim for a new assembler is "it finds novel viruses in a mixed
sample without a reference", a 2015 de Bruijn metagenome assembler already
does that, and there is no point building for it.

### The near-relative is where everything fails

RaTG13 sits 4% from SARS-CoV-2, at 25x against its 100x. In seed 2
metaSPAdes, Haploflow and metaFlye **all three collapsed it entirely**: zero
contigs, zero coverage, while SARS-CoV-2 came back complete. MEGAHIT got 55%
of it in fragments. **No tool recovered it in any run.** A novel sarbecovirus
next to a known one is exactly the discovery case that matters, and the
field's standard pipelines delete it.

Failures go the other way too. metaSPAdes (seed 1) emitted a 13.8 kb contig
running HIV-1 subtype B from 1275 to 9275, wrapping through B's start, then
continuing as subtype C, every 500 nt window matching its part at 100%
identity: a fabricated inter-subtype recombinant. Haploflow (seed 2) joined
RaTG13 to SARS-CoV-2 and West Nile lineage 1 to lineage 2. metaFlye joined
HHV-5 to a bacterial chromosome in both seeds. viralFlye, which is the
discovery pipeline rather than just an assembler, filtered metaFlye's
assembly down to 10 viral contigs and in doing so lost three relatives that
were present in the assembly it was given: the classification step throws
away the novel ones.

### And nobody resolves a swarm

The ten env haplotypes in the same sample: 0 or 1 of 10 recovered by every
tool in every run, with the chimaeras that MEGAHIT does produce being
exclusively joins between swarm haplotypes. The swarm is invisible to all of
them while the between-species part of the same sample is handled well.

### The threshold

`ladder.md`, MEGAHIT, two members at 100x and 25x, twelve rungs of real
sequence pairs from 0.15% to 46% divergence:

* below about **3%**: the minor member is erased, output is one mosaic contig;
* **3% to 5.5%**: 27% to 77% recovered, in fragments, and two rungs 0.08%
  apart straddle total collapse, so no threshold can honestly be quoted in
  this band;
* **10.5% and up**: one contig each, both over 99% covered, no chimaeras.

Intra-host swarms live at 0.1-3%. Near-relative novel viruses live at 1-5%.
Both sit inside the band where the incumbents either erase a genome or return
a coin toss. That is the gap, it is measured, and it is what hapax has to beat
to be worth building.

### What would falsify this

* Parameter tuning. Every tool ran at defaults except Strainline (seed length
  lowered for 2.5 kb amplicons). Haploflow in particular is built for
  two-strain mixtures and its `--strict` and `--error-rate` settings were not
  swept; its poor showing here is partly a statement about defaults.
* One seed per ladder rung, two per community. The 3-5.5% band is exactly
  where more seeds would change the picture, and it is the band that matters.
* The reference-guided long-read tools (RVHaplo, HaploDMF) were run only on
  the haplotype panels, not the community, and Strainline did not finish
  inside the sandbox's limits. A tool given the right reference may well do
  better on the near-relative case than these de novo assemblers.
* Simulated reads throughout. ART and Badread are what the tools' own papers
  use, but a real library has biases neither models.