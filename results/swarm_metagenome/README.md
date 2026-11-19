# Several real viral swarms in one metagenome

What happens when a sample holds more than one within-host swarm, which is the
case hapax exists for and the case the stress-test work identified as the real
gap (between-species novelty is already handled by a plain de Bruijn
assembler).

Reproduce with:

```
python experiments/fetch_swarm_panels.py                                  # real panels, once
python experiments/swarm_metagenome.py --seed 1 --n-seeds 5 --out results/swarm_metagenome
python experiments/report_swarm_metagenome.py results/swarm_metagenome
```

`work/seed*/manifest.json` holds the git SHA, parameters and timings per seed.
Short version: **on this benchmark hapax is behind metaSPAdes and MEGAHIT on
every aggregate measure except contigs per swarm, because it merges swarms
from different hosts and loses three of eight.** The single thing it does
better is read the variant structure of the most diverse swarm.

## The sample

Eight swarms, each one real single-genome-amplification or clonal sequence data
from one host (`experiments/fetch_swarm_panels.py`, accessions in
`data/MANIFEST.md`). Nothing about the sequences is simulated; only the reads
are.

| swarm | virus, region | haplotypes | length |
|---|---|---|---|
| `hiv_env_ZM247F` | HIV-1 C *env* | 72 | ~3.2 kb |
| `hiv_env_ZM246F` | HIV-1 C *env* | 61 | ~3.1 kb |
| `hiv_gagpol_2139` | HIV-1 *gag-pol-vif-vpr* | 19 | ~5.4 kb |
| `hcv_Pt3P` | HCV core..NS3, primary strain | 49 | ~4.9 kb |
| `hcv_Pt3S` | HCV core..NS3, superinfecting strain in the same host | 39 | ~4.9 kb |
| `hcv_Pt7` | HCV core..NS3 | 54 | ~4.9 kb |
| `hbv_p03` | HBV complete genome | 24 | 3215 nt |
| `hbv_p10` | HBV complete genome | 24 | 3215 nt |

Background: *E. coli* K-12, *Bacteroides fragilis*, *Staphylococcus
epidermidis* at 5x, phages lambda, T4, crAssphage and phiX174 at 30x, and one
HIV-1 HXB2 genome at 30x, so the swarms are a minority of the reads. Each swarm
gets 300x in total, split across its haplotypes with log-normal abundance
(sigma = 1), redrawn per seed, so a few haplotypes dominate and many sit near
1x. Reads are ART HiSeq 2500, 2x150, ~265k pairs. hapax is single-end and gets
both mates as independent reads; the others get pairs.

The sister pairs are the point of the design. By consensus identity: the two
*env* hosts are 86% identical, the two HBV hosts 91%, HXB2 is 89% identical to
the *gag-pol* swarm, and the two HCV strains in host Pt3 share no exact 31-mers
at all. A tool that cannot tell two hosts' swarms apart will show it here.

## Headline: 5 seeds, mean ± sd

| tool | swarms found /8 | contigs per found swarm | completeness | cross-swarm chimaeras | haplotypes resolved | variant recall | variant precision | wall (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|
| hapax codon | 4.80 ± 0.45 | 1.25 ± 0.09 | 0.57 ± 0.05 | 3.60 ± 0.89 | 0.00 | 0.47 ± 0.05 | 0.55 ± 0.06 | 284 ± 105 | 2.29 |
| hapax codon, tier-1 groups | 5.00 ± 0.00 | 1.00 ± 0.00 | 0.62 ± 0.00 | 3.80 ± 0.45 | 1.00 | 0.48 ± 0.03 | 0.52 ± 0.03 | 284 ± 109 | 2.31 |
| hapax nucleotide (ablation) | 4.80 ± 0.45 | 1.25 ± 0.09 | 0.57 ± 0.05 | 3.60 ± 0.89 | 0.00 | 0.47 ± 0.05 | 0.55 ± 0.06 | 281 ± 107 | 2.33 |
| MEGAHIT | 7.60 ± 0.89 | 3.28 ± 0.43 | 0.73 ± 0.04 | 1.80 ± 0.84 | 1.00 | 0.53 ± 0.08 | 0.86 ± 0.02 | 301 ± 121 | 0.39 |
| metaSPAdes | 7.00 ± 0.00 | 1.63 ± 0.22 | 0.81 ± 0.04 | 2.00 ± 0.71 | 1.40 ± 0.55 | 0.56 ± 0.04 | 0.82 ± 0.02 | 374 ± 13 | 1.31 |
| Haploflow | 8.00 ± 0.00 | 3.95 ± 0.78 | 0.87 ± 0.03 | 12.40 ± 3.13 | 1.00 | 0.65 ± 0.06 | 0.75 ± 0.04 | 126 ± 4 | 3.46 |

Definitions are in `experiments/score_swarm_metagenome.py`; all of them are
contig-only, so every tool is judged by the same rules. In brief: a contig
belongs to the swarm (or background genome) that owns most of its 250-nt
windows; a **cross-swarm chimaera** is a contig where a second swarm owns at
least 10% of its windows; **completeness** is the longest single contig's
aligned span on its best haplotype; **haplotypes resolved** counts haplotypes
recovered at better than both 1% divergence and half their distance to their
nearest neighbour, over 95% of their length; **variant recall and precision**
are against the abundance-weighted minor-allele frequencies of the swarm's own
alignment, from reads mapped back to each tool's contigs (the standard
assemble-then-map pipeline, identical for every tool). metaSPAdes runs with
`--only-assembler`: with read correction on it did not finish in an hour.

## Findings

**1. hapax merges swarms from different hosts, and that costs it three swarms
of the eight.** `hbv_p03`, `hcv_Pt3S` and `hiv_env_ZM246F` come back empty or
nearly so on most seeds. The mechanism is not that reads are discarded: on the
isolated *env* pair, the single merged contig holds 5,522 reads from one host
and 4,243 from the other, and its consensus is 97% identical to the more
abundant host's haplotypes. The minority host is still present, but only as
variant calls on its sister's consensus. One contig per population is doing
exactly what it says, with the population boundary drawn in the wrong place.
Both other short-read assemblers draw it correctly: MEGAHIT and metaSPAdes
each make about two cross-swarm chimaeras, against hapax's 3.6 to 3.8.

**2. Where hapax does hold a swarm, it holds it in one piece, and that is the
only place it wins.** On the swarms it keeps (`hiv_gagpol_2139`, `hbv_p10`,
`hcv_Pt3P`) it returns 1.00 contigs at 1.00 completeness on every seed, where
MEGAHIT needs 2.2 to 3.8 contigs for the same swarms. But one contig is not by
itself worth anything: metaSPAdes reaches 1.00 contigs at 1.00 completeness on
`hbv_p10` and `hcv_Pt3P` too, and beats hapax's variant recall there (0.89
against 0.81 on `hbv_p10`). The exception is the most diverse swarm,
`hiv_env_ZM247F`, where hapax with tier-1 grouping recovers 0.49 of the true
variant sites against metaSPAdes' 0.24 and MEGAHIT's 0.21 — roughly double,
when it keeps that swarm whole, which is on four of five seeds.

**3. hapax's variant precision is its worst number: 0.52 to 0.55, against 0.82
to 0.86 for the de Bruijn assemblers.** A merged sister contig is the likely
cause: reads from a host 14% divergent from the consensus they are aligned to
produce columns of apparent variation that is really misalignment. The one
swarm where precision is high (`hcv_Pt3P`, 1.00) is the one whose sister shares
no exact 31-mers with it, so no merge is possible.

**4. Codon scoring and nucleotide scoring are indistinguishable here, and part
of the reason is a bug in hapax's own region caller.** The two modes agree to
three decimal places on every aggregate number. Dual-coding regions fall back
to nucleotide scoring by design, and on real genomes the ORF caller
(`find_orfs`, `src/consensus.rs`) overcalls them badly: it accepts every
stop-free stretch of 300 nt in all six frames, with no start codon and no
scoring, so the HCV polyprotein contigs — one long ORF, plus a short ARF — come
back 60% dual-coding, and HIV *gag-pol* 21%. HBV's 71% is about right, since
HBV genuinely is overlapping ORFs almost throughout. So codon scoring is
switched off over most of the coding sequence in this sample, and no
codon-versus-nucleotide number from real data should be believed until the ORF
caller is fixed. (Found independently by the wastewater work on a real
astrovirus genome.)

**5. Haploflow is the cautionary comparison.** It finds all eight swarms at the
highest completeness of any tool and the best variant recall, and makes 12.4
cross-swarm chimaeras doing it, on a sample where MEGAHIT makes 1.8. Chasing
haplotype-level output without a chimaera budget produces exactly what
CLAUDE.md warns about.

**6. hapax also emits about 7,400 contigs of bacterial background it has no
business assembling**, against MEGAHIT's 4,170 and metaSPAdes' 3,370 for the
same reads, at 2.3 GB against MEGAHIT's 0.39 GB. Total contig counts across
tools are therefore not comparable, and a real pipeline would put a viral read
filter in front of hapax rather than feed it a whole metagenome.

## Three fixes for the merge, and why each failed

`experiments/sister_pairs.py` isolates one sister pair at a time, with no
background, at 30x and 300x, so the only variable is how the grouping or
splitting rule behaves. `results/sister_pairs/summary.md` has the table.

**Thresholds cannot do it.** Raising `--bridge-frac` to 0.6, raising
`--min-bridge` to 50, or using tier-1 components as the first-stage clusters
all leave the pair merged, at 30x exactly as at 300x. The reason is that the
joining overlaps are not errors: over a 150-nt window two hosts' sequences
often differ at one or two positions, which is inside sequencing error, so
those really are nucleotide-identical overlaps, and CLAUDE.md's own join
criterion accepts them. The two *env* hosts share no exact 150-mers at all and
still merge on near-identical overlaps; the two HBV hosts share 709 exact
150-mers and 194 exact 300-mers, which no read-length-limited assembler can
refuse. `split-cut` does produce one contig per swarm instead of one for both,
but each is half each host (the cut runs across the genome), so it trades a
merge for two mosaics and the chimaera count rises. **No per-pair score, of any kind, can separate two hosts' swarms.** That
is the central result of this experiment, and it is a statement about the
information in a read pair, not about hapax.

**Post-consensus phasing converges to mosaics** (`--split-phase`,
`src/phase.rs`, off by default). Fitting two populations to the contig's
variant sites — a two-component Bernoulli mixture by EM, then rounds of
sequence-level reassignment — does separate most reads, but it settles at about
70/30 per side rather than 99/1, and a 70/30 side gives a mosaic consensus:
each piece comes back 92% identical to *both* hosts instead of 98% to one. The
prior attempt, chaining neighbouring sites in strong linkage, fails outright:
only 97 of 677 neighbouring site pairs on the merged *env* contig reach
r² ≥ 0.9 and the longest unbroken chain is 5 sites, because each host is itself
a 60-haplotype swarm and the between-host sites are interleaved with far more
within-host ones. The unit test in `src/phase.rs` shows the mixture fit works
on two clean populations with within-population noise; two swarms are a harder
object than that.

**A sparse cut in the overlap graph cuts the genome, not the population**
(`--split-cut`, `layout::split_sparse`, off by default). The edge counts
suggest a cut exists: overlaps between the two *env* hosts' reads are 1.1% of
the overlaps within the sparser host at 300x, and 4.4% for the HBV pair. That
ratio is itself coverage-dependent — at 30x it is 21.5% and 11.6%, because the
within-host overlaps are the ones that thin out — so a fixed threshold on it
could not travel between samples even if the cut were the right cut. But in
the overlap graph of a linear genome the cheapest cut is always across the
middle of the genome, so label propagation splits by coordinate and each side
comes back holding both hosts in equal measure. Any cut rule needs a constraint
that both sides span the same coordinate range, which is phasing again, with
the reads anchored in layout coordinates. That is the next thing to try, and it
is not a tweak.

## Where this leaves the kill criteria (METHOD.md §5)

| criterion | status after this experiment |
|---|---|
| No gain over nucleotide assemblers | **Failing, on short reads with several swarms.** metaSPAdes beats hapax on swarms found, completeness, chimaeras, variant recall and precision, at half the memory. hapax's only edge is variant recall on the most diverse swarm, about double metaSPAdes', and only when it keeps that swarm whole. The internal codon-versus-nucleotide ablation shows no difference at all here. |
| Collapsing erases differences in overlapping genes | **Untested, and now blocked by a bug.** The region caller overcalls dual-coding on real genomes (finding 4), so the dual-coding path is being taken where it should not be. HBV, where two thirds of the genome genuinely is overlapping ORFs, is the test to run once that is fixed. |
| Too slow | **Passed.** 280 s on 530k reads at 2.3 GB, against metaSPAdes' 374 s. The HMM fraction stays near 0.01 at this coverage. |
| Chimaeras increase | **Confirmed, and worse than the simulator showed.** 3.6 to 3.8 cross-swarm chimaeras against 1.8 for MEGAHIT, every one of them a merge of two hosts' swarms. The simulator's sister populations were 0.30 subs/site apart, which is far enough to be separable; real sister swarms at 9 to 14% are not. |

The honest summary is that hapax's premise — one population, one variant-aware
contig — is not wrong about what is useful, but the tool has no mechanism for
deciding where one population ends and the next begins, and on real
multi-swarm data that decision is the whole problem. Everything in this
experiment points at the same missing piece: a position-anchored phasing step
that can divide a component into populations, with the chimaera number as its
scoreboard.
