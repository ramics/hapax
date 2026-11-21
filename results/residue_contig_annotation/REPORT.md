# What the residue contigs are, read by profile rather than by index

`experiments/residue_contig_annotation.py`, 5 seeds, aggregated by
`experiments/report_residue_contig_annotation.py`. **Every number is in
[`TABLES.md`](TABLES.md)**, generated from the committed per-seed JSONs.

Three questions that do not need an index to name an organism: does it code (Prodigal in
metagenomic mode, scored as coding density), is it structured non-coding RNA (Rfam cmscan
at the gathering thresholds), and does any protein family recognise it, RdRp in particular
(Pfam-A hmmsearch at each family's own gathering threshold).

Two piles were asked, and **they answer in opposite directions**, which is the result.

## Pile one: the 200 highest-mass unidentifiable contigs

The head of the residue's read-mass distribution — 36.9% of the residue's contig mass,
median 660 nt at median 57x, nine above 500x. Long enough and deep enough for a profile
search, abundant enough to matter, and nothing identifies them. If there were a hidden
abundant virus in this library, this is where it would be.

| | these 200 | identified contigs | unidentifiable contigs | re-derived | shuffled |
|---|---|---|---|---|---|
| coding density | **0.604** | 0.919 | 0.818 | 0.891 | 0.203 |
| any Pfam | **7.5%** | 69.5% | 11.2% | 38.3% | 0.0% |
| an RdRp family | **2.5%** | 28.5% | 1.0% | 21.1% | 0.0% |
| any Rfam | **3.0%** | 1.0% | 1.0% | 0.0% | 0.0% |

Read the controls first. The shuffle null reads 0.203 coding density and fires zero Pfam
and zero Rfam, so all three fields are reading real sequence rather than composition. The
identified control is where every protein field is highest, as it must be.

**The abundant unidentifiable material is less protein-coding than ordinary unidentifiable
material, and more structured.** Coding density 0.604 against 0.818 for contigs of the same
length drawn from the same unidentifiable set; 61% of these contigs have half their length
inside a predicted CDS against 90% of the control. Rfam is the only field where this pile
leads every other pile, at 3.0%, and the hits are large- and small-subunit rRNA.

That is a negative for the obvious hypothesis and a positive for a duller one. The most
abundant thing in this library that nothing can name is not a hidden coding virus; it
looks like non-coding and structural RNA, some of it ribosomal sequence SILVA did not
catch. **A profile search on the deepest residue contigs does not find a virus, and the
reason is not that the contigs were too short.**

## Pile two: the 2,832 genuinely new contigs from unplaced pairs

Contigs assembled from pairs that neither SILVA nor the sample's own assembly touched, and
that are then under 10% covered by that assembly, so they are new sequence rather than a
re-derivation. The split is the residue thread's (`experiments/residue_contig_overlap.py`);
without it, 650 re-derived contigs would have sat inside their own control.

| | new contigs | identified | unidentifiable | re-derived | shuffled |
|---|---|---|---|---|---|
| coding density | **0.790** | 0.916 | 0.803 | 0.885 | 0.282 |
| any Pfam | **23.1%** | 50.7 (49.9-51.5) | 6.6 (6.2-7.0) | 30.2 (29.8-30.6) | 0.0 |
| an RdRp family | **11.7%** | 16.2 (15.4-16.7) | 0.70 (0.56-1.02) | 19.5 (19.3-19.8) | 0.0 |
| any Rfam | **2.1%** | 0.97 (0.71-1.06) | 0.94 (0.74-1.20) | 2.5 | 0.0 |

Here the direction reverses. These contigs carry an RdRp **seventeen times as often** as
length-matched contigs from the sample's own unidentifiable set, and about three quarters
as often as contigs we did identify. The dominant family is `RNA_replicase_B` — Leviviridae
replicase, ssRNA phage — 1,615 hits pooled over seeds against 83 in the unidentifiable
control.

**State the bias before reading further.** The pairs these contigs were built from were
selected by SILVA and own-contig alignment alone, with no read-level Kraken2 and no
DIAMOND, so the input is a superset of the cascade's residue and some of it is material
the cascade would have identified. That bias inflates the protein numbers here; it cannot
deflate them. The comparison that survives it is against the unidentifiable control, which
was built at the cascade's own 40-base operating point: whatever this pile is, it is not
the same population as the contigs we already call unidentifiable, and it is not an
assembly artefact, because the shuffle null is flat at zero.

The re-derived arm is the useful third control and it reads close to the identified pile
(30.2% Pfam, 19.5% RdRp). Same assembler, same read class, already in our assembly. So
assembling unplaced pairs produces recognisable phage sequence whether or not we had it
before; what the new pile adds is that a quarter of it is recognisable and was not in the
assembly at all.

## What this changes

1. **Depth and novelty point in opposite directions here.** The deepest unidentifiable
   contigs are the least coding; the newly assembled ones are the most protein-rich. A
   pipeline that ranks residue by abundance and runs profiles on the top of the list will
   find rRNA. Both lists are worth keeping and they are not the same list.
2. **There is a protein tier's worth of signal in the unplaced pairs.** 23% of new contigs
   carry a Pfam family at the curators' own threshold, from reads that four nucleotide and
   two protein indexes did not place. That is consistent with the holdout ladder, where
   nucleotide identification dies at 20% divergence and protein is last standing.
3. **A small annotation gap is measurable.** Rfam fires on 2.1% of new contigs against
   0.94% of the control, mostly large-subunit rRNA, so roughly one point of this pile is
   ribosomal sequence the nucleotide tier missed. Small, and it is the first time this
   project has put a number on it.
4. **Rfam is otherwise close to silent on the residue**, which closes part of the
   non-coding ladder left open in `results/README.md`: at these contig lengths, Rfam adds
   about one point over its own control and nothing like a new family.

## Operating point

Prodigal (pyrodigal 3.7.1) in metagenomic mode, because 1.3 Mb of unrelated fragments
cannot train a single-genome model. Pfam-A (pyhmmer 0.12.3 hmmsearch, models streamed from
the flatfile) at each family's **gathering threshold**, which is the curators' operating
point rather than one chosen here; RdRp counted by accession, not by name match, and
reverse transcriptase (PF00078) reported separately and never folded in. Rfam via
Infernal cmscan `--rfam --cut_ga --nohmmonly` with `-Z` fixed across piles so E-values are
comparable. Controls drawn from the sample's own contigs in log length bins of ratio 1.25
from 200 nt, without replacement, with every query name removed from every pool first. The
shuffle null is dinucleotide-preserving (Altschul-Erickson), because stop-codon density is
a function of dinucleotide composition and a mononucleotide shuffle would hand Prodigal
longer open frames than the real sequence.

## What it cannot say

A Pfam or Rfam hit on a residue contig means the read-level indexes missed it, not that
the contig is novel — that is what tier 2 is. Absence of a hit is a statement about Pfam
and Rfam at their gathering thresholds and about a median of 319 nt, which is short for
both tools. The `control_resasm_rederived` arm could not fill the upper length bins (2,194
contigs short over 6 bins on the new pile), so it is matched only below the top of the
length range and leans short; the other two controls are matched within three contigs.
