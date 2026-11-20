# Lineage labeller audit

Script: `experiments/labeller_audit.py`, seeds 1-5, 100 contigs per (rung, condition, seed),
on the seven real rabies genome pairs of the archived two-lineage ladder (2.0-9.7%
divergence). Rows in `rows.jsonl`, pooled table in `summary.md`, run manifest beside them.
Runtime 4.3 minutes.

## Finding

The per-window caller that produced every mosaic count in
`experiments/codon_assembler/results/abundance_sweep` aligned each 300 nt window to both
lineages in the forward orientation only. Assemblers emit contigs on either strand. On
a reverse-strand contig both edit distances are near random, so the call is a coin flip.

| pure lineage-A contigs, 700-3,000 nt | archived caller | strand-aware caller |
|---|---|---|
| forward strand, 0 / 0.5 / 1% added error | 0.0% / 0.0% / 0.0% called mosaic | 0.0% / 0.0% / 0.0% |
| reverse strand, 0 / 0.5 / 1% added error | 65.8% / 64.9% / 66.2% (range 52-78%) | 0.0% / 0.0% / 0.0% (max cell 1%) |

Each cell is 3,500 contigs (7 rungs x 5 seeds x 100). So the archived caller has no false
mosaics at all on correctly oriented contigs, even with 1% consensus error, and its whole
false-mosaic rate is the strand bug. With roughly half of contigs reversed, that explains
13 of 42 single-lineage controls being called mixtures.

On synthetic one-join mosaics (lineage A then lineage B at homologous coordinates, the easy
case), the strand-aware caller detects 99.9% forward and 100.0% reverse; the archived caller
detects 99.9% forward and 77.2% reverse, the latter by chance. Join position is only known
to the window: median error 150 nt for both callers, which is half the window. Base-resolution
join features need a per-diagnostic-site caller, not this one.

## What it changes

The strand-aware caller is fit to count mosaics. Every mosaic count, lineage-switch count and
join observable in the archived rows is void and cannot be rescored, because the contigs were
not kept: rescoring means re-assembly. Lengths, mass, breadth and the 1,500 nt crossing event
come from minimap2 PAF, which is strand-aware, and are unaffected.
