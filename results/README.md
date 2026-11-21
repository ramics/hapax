# Experiment archive

One line per experiment: what it measured, what it settled, and where to read it. The spec
that these results justify is [`../SPEC.md`](../SPEC.md).

**Everything here is kept, including what was killed and what was withdrawn.** The negatives
are load-bearing — they are why the spec's target is a measurement rather than a call — so
nothing is deleted for tidiness. Where a result has been superseded or retracted, the row
says so and the files stay put.

Conventions: every script takes `--seed` and writes a run manifest (seed, git SHA,
parameters, timings) beside its output; headline results run across at least five seeds with
spread reported. Scripts live in [`../experiments/`](../experiments/) under the same name as
their results directory.

## The current direction — subtract the known, measure the residue

| directory | what it measured | status |
|---|---|---|
| [`readlevel_features/`](readlevel_features/REPORT.md) | Ceilings on single-read features — ORF-likeness, k-mer composition, codon usage, GC. None has novelty polarity, because a wastewater background is coding too. Includes the prior-art review and the Viromancer/Lyra-Virus review. | live; the reason no axis is thresholded |
| [`residue_structure/`](residue_structure/REPORT.md) | SPEC.md §4.2's gate experiment: is the residue structured or diffuse in the joint space? **No clusters beyond its own covariance (k=1, 5/5 seeds), the joint space beats the best single axis by 0.013 AUC, but a planted genome's reads find each other by neighbourhood down to about 1 in 100,000 and not below.** Adds exact whole-library 31-mer coverage shape as an axis set. | live; the gate is answered |
| [`profile_reach/`](profile_reach/REPORT.md) | SPEC.md's two unmeasured assertions about profile search, on the spike ladder with shuffled and known-material controls in the same pass. **A negative and a correction**: Pfam-A is *weaker than translated search at every rung* (91.5% vs 39.9% at 0% divergence, 61.3% vs 18.6% at 10%, 10.0% vs 5.8% at 20%), DIAMOND does **not** read zero past 20% (39.9% at 20%, 5.7% at 30% at bitscore >=40 — the "zero" was a cascade removal figure read as a search's reach), and the no-sequence-hit-yes-profile-hit gap is 3.8% on the holdout against 2.6% on already-identified reads with overlapping ranges. **And a positive**: merging the mates of the 35.1% of pairs that overlap, compared against mate 1 of the *same* molecule, takes recognition from 10.5% to 38.4% at 20% divergence and from 0.0% to 7.1% at 30%, with the shuffled null at 0.0% at every length — so part of the published 10-20% recognition wall is a read-length wall. | live, 5 seeds; SPEC.md §4.2 item 5 withdrawn, §4.3 added |
| [`depth_rarefaction/`](depth_rarefaction/REPORT.md) | What deeper sequencing would buy, priced from the library in hand rather than by spending: a fixed query block scored against nested subsamples at 1/8, 1/4, 1/2, 1. **A doubling buys a second observation, not assembly depth.** Of the reads fully unwitnessed at half the library, 44.0% (40.6-46.2) gain a witness at full library and 0.3% reach depth >=5; the rescue rate *rises* along the ladder (24.4, 34.0, 44.0) because each doubling adds twice the absolute reads of the one before, so the last transition is a floor for the next. Across 8x effort the below-depth-5 fraction moves only 38.5% to 24.4%. | live, 5 seeds; no 10x figure is reported and the reason is in the report |
| [`peptide_absence/`](peptide_absence/) | Peptide k-mer absence under family holdout. Saturates at k=10; graded at k=6. | live |
| [`peptide_absence_noise/`](peptide_absence_noise/REPORT.md) | The instrument error floor under peptide absence. Peptide beats nucleotide (synonymous substitutions), but one 107 nt read cannot carry a one-in-a-million call. | live |
| [`berlin_quality/`](berlin_quality/) | This library's actual error rates: 0.00089 on R1, 0.00117 on R2, modal length 107, flat profile. Every simulation uses these. | live |
| [`absence_curve/`](absence_curve/) | Absence against protein index size, to 1.73 G 10-mers. **The curve has not flattened** — each tenfold step still buys ~20 points, so the far end must not be extrapolated. | live; corrects an earlier "diminishing returns" claim |
| [`absence_residual/`](absence_residual/) | What the residue is made of: part maps the sample's own assembly and sits outside gene calls, part is whole molecules missing from it with the mate agreeing. | live |
| [`rrna_content/`](rrna_content/REPORT.md) | rRNA is about a third of this library, and the conservation multiplier (388 per occurrence against 3.49 per distinct k-mer) that makes a few badly annotated public entries fool a protein index. | live; the generalisable result |
| [`rrna_contigs/`](rrna_contigs/), [`rrna_match_source/`](rrna_match_source/) | The sample's own ribosomal contigs, and which public entries the spurious protein matches come from — 63 of 4.12 M. | live |
| [`mate_width/`](mate_width/REPORT.md) | The wide-insert question. Width is negligible for hitting coding sequence, worth a lot for novelty; no pair rule reaches one per million, so the pair is a multiplier on aggregation. Includes the insert prior-art review. | live |
| [`insert_profile/`](insert_profile/), [`mate_separation/`](mate_separation/), [`mate_rho/`](mate_rho/), [`novelty_blocks/`](novelty_blocks/) | The inputs to that verdict: the library's measured insert distribution, novelty per residue against mate separation, the correlation between mates' false flags, and the block structure of novelty in holdout genomes. | live |

The three-tier cascade itself — the re-cut subtraction table, the index-size ladders, the
holdout spike ladder and the NCBI ground truth — is the newest work and lands with the
open pull request from the `subtract the known` thread. The spec cites it as
`results/subtract_recut/`, `results/spike_holdout/`, `results/nucleotide_ladder/` and
`results/blast_check/`; those paths appear here when that branch merges.

## Strain separation

| directory | what it measured | status |
|---|---|---|
| [`pair_ladder/`](pair_ladder/REPORT.md) | The de-collinearised ladder: 2,205 assemblies, 21 real genome pairs, 5 seeds, 0 false mosaics in 630 controls. Below ~3.8% divergence every tool fuses lineages at any depth; MEGAHIT separates from ~5% up. PenguiN's raw output holds clean pieces of both lineages at 2%, buried among mosaics — so the reads carry the information and path selection is what is missing. Confirmed against real Berlin background. | live |
| [`labeller_audit/`](labeller_audit/README.md) | The strand bug in the lineage labeller: the archived caller gave 0% false mosaics on forward-strand contigs and 65.8% on reverse. | live, and it **voids every mosaic count published before the fix** |

## Archived: the codon-aware assembler

The first direction, kept verbatim and still runnable under
[`../experiments/codon_assembler/`](../experiments/codon_assembler/) — spec, crate, scripts
and its own `results/` tree, laid out as the repo root once was. Ships as a negative result,
per its own kill criteria.
