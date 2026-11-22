# The contigs that sit far from any published protein

`experiments/divergent_contig_workup.py`, seed 1, in container, no spend, 342 s on four
cores. Output: `workup_s1.json`.

Imo, 2026-11-22: *"those assembled viruses at weak identity are going to turn out to be
absolutely key. what can we do to enrich those data?"* This is the first answer, and it
begins by correcting the question.

## 1. "Divergent strain" is the wrong reading of the divergent end

Strain-level variation lives at 90-99% protein identity. The rows leading the page sit at
30-45%, where the named subject is the nearest neighbour in a sparse database rather than
an identification. `divergence_class` now states the band on every row and the word
"strain" appears only above 90%. The 92.2% astrovirus is a divergent strain; the 30.6%
*Sholavirus* is not a strain of anything.

## 2. The number we were ranking by is not a property of the contig

`known_strain_divergence.py` takes one DIAMOND hit per contig. That hit is whichever gene
won on bitscore, and on this material that is as often the fastest-evolving gene as the
slowest.

| | divergent arm | near-reference control |
|---|---|---|
| contigs | 149 | 135 |
| with two or more scorable genes | 67 | 54 |
| **median spread between a contig's own genes** | **17.4 points** | **2.35 points** |

The control is the calibration that makes the test readable: on contigs that match their
reference, a contig's genes agree with each other to within about two points. On the
divergent arm they disagree by seventeen, and it runs in both directions.

- *Lirnavirus borborovicinum*: 88.6% whole-contig, polymerase **39.8%**, spread 48.8.
- *Mihkrovirus pelocola*: 41.8% whole-contig, polymerase **58.6%**, spread 37.7.

By gene role, divergent arm against control:

| gene | divergent median | control median |
|---|---|---|
| polymerase | 66.8% (n=105) | 97.4% (n=63) |
| maturation protein | 63.5% (n=77) | 97.0% (n=37) |
| coat protein | 77.9% (n=63) | 98.2% (n=48) |

That ordering is the expected one: the polymerase places a virus in a family, the
maturation protein moves fastest. **The ranking is now on the polymerase where the contig
has one** (103 of 149), with the whole-contig hit as the stated fallback, and every row
carries the spread so one figure can never stand in for a contig again.

## 3. The assembler's coverage header understates read support about fourfold

Mapping all 17,462,448 mates back (mate 2 trimmed 16 nt, as everywhere in this repo):
*Rodtovirus pelovivens* reads 37.7x in the SPAdes header and **170.7x** mapped. Median
mapped depth on the divergent arm is 18.4x, maximum 467x, and median breadth covered is
99.3% with nothing below 90%.

**Seven contigs are near-full-length, deep and divergent at once** — 3 to 4.5 kb is a
complete *Fiersviridae* genome:

| taxon named by the nearest subject | nt | depth | polymerase | evenness |
|---|---|---|---|---|
| Kahshuvirus borborenecus | 3,484 | 467x | 44.9% | 0.84 |
| Fiersviridae | 3,452 | 402x | 48.6% | 0.82 |
| Koteshevirus lutadaptatum | 3,058 | 297x | 46.6% | **0.51** |
| Sholavirus sp. 'caenicola' | 3,725 | 282x | 83.8% | 0.77 |
| Rodtovirus pelovivens | 3,930 | 171x | 45.6% | 0.81 |
| Dugnivirus solivicinum | 3,996 | 142x | 47.8% | 0.96 |
| Chimpavirus luticola | 4,526 | 116x | 84.8% | 0.81 |

Five of the seven sit below 50% polymerase identity. These, not the thin single-contig
rows at the top of the old ranking, are the material worth effort.

## 4. The chimaera flag is honest and unexciting; evenness is the field that separates

Subject disagreement across a contig's ORFs fires on 46.3% of the divergent arm and 39.3%
of the control. **Firing almost equally on both means it is describing how the assembler
joins genes, not divergence**, so it is reported and nothing is dropped on it. Coverage
evenness does separate: three divergent contigs have under half their length within 2x of
their own median depth, *Koteshevirus lutadaptatum* among the deep seven, and those are
the rows to distrust.

Pfam reach is the same on both arms (85.9% against 80.7%), which is the expected direction
and confirms the contig side of the profile line: profiles are a measured negative at read
level (SPEC §4.2 item 5, withdrawn) and work here because an ORF runs to hundreds of
residues.

## 5. What this is

Coat plus maturation protein plus polymerase, 3 to 4.5 kb, is the architecture of
single-stranded RNA bacteriophage. These infect bacteria. That is consistent with the
standing negative that nothing in this sample is human-associated, and it means the
public-health reading of this panel is indirect.

## 6. The gene that picks the host is the one furthest from the database

Every one of the seven deep genomes carries the complete *Fiersviridae* gene set, and the
domains confirm it rather than inferring it: `RNA_replicase_B`, `Levi_coat`, `Phage_mat-A`.
The alignments run the full length of each ORF (polymerase 460-785 aa, coat 130-160 aa,
maturation protein 290-515 aa), which is what distinguishes a complete genome from a
fragment that happens to be long.

In this class the maturation protein is the receptor-binding protein: it attaches to the
host pilus and it is what decides which bacterium the phage can infect. It is also the
most divergent gene on these contigs.

| | divergent arm | near-reference control |
|---|---|---|
| contigs with both a polymerase and a maturation protein | 50 | 19 |
| maturation protein below the polymerase | **86%** | 63% |
| median polymerase minus maturation protein | **7.7 points** | 0.4 points |

The control is what makes this readable: on contigs that match their reference the two
genes agree to within half a point, so the 7.7-point lean on the divergent arm is not a
property of the two proteins' general evolutionary rates as this instrument sees them. Six
of the seven deep genomes follow it; *Chimpavirus luticola* is the exception.

**The direction of the bias is stated and not controlled**: the maturation protein is the
longest and the fastest-evolving gene in this family regardless of host, so a lean of this
sign is expected even with no host difference at all. What the control establishes is that
the lean is specific to the divergent arm. Whether these phages infect different bacteria
cannot be answered from sequence identity, and the route that could answer it is
co-occurrence with bacterial abundance across the plant's 115 other libraries.

Two contigs carry a short ORF with no subject and no domain at all - 74 aa on *Rodtovirus
pelovivens*, 121 aa on *Dugnivirus solivicinum*. Position and length fit the lysis gene,
which is small, overlapping and poorly conserved in this family, but nothing here
establishes that. They are the clearest candidates for a structural prediction.

## 7. A 319 GB nucleotide index collapses a sixth of the divergent panel

Every identity in sections 1 to 6 is against `viral.dmnd` plus `sprot.dmnd`, a few hundred
megabytes. The objection that answers is that the distance is a property of the index. So
the same 284 contigs, both arms, went to NCBI core_nt (319 GB, db version 2026-10-22) on
one EC2 instance: `blastn -task dc-megablast`, 8 threads, 4,427 s, $0.42 of instance time.
Harvest in [`core_nt_s1.json`](core_nt_s1.json), summarised by
`experiments/bigindex_harvest.py` from the committed row file.

A hit counts at E < 1e-5, and identity is read beside the fraction of the contig the
alignment covers, never alone.

| | near-reference arm | divergent arm |
|---|---|---|
| contigs queried | 135 | 149 |
| significant hit | **135** | 147 |
| no significant hit | 0 | 2 |
| median identity | 97.55% | **81.56%** |
| median alignment | 1,332 nt, 98.8% of contig | 846 nt, 82.8% of contig |
| **≥90% over ≥half the contig** | **131 (97.0%)** | **25 (16.8%)** |
| **<80% identity, or under half covered** | **0** | **83 (55.7%)** |

The near-reference arm is the calibration, and it passes: core_nt returns all 135 at a
median 97.55% across 98.8% of the contig, and the four that miss the strain cut miss it
narrowly (85.3% and 83.1% identity, 89.6% twice). Nothing about the search moved.

On the divergent arm the big index does real work and does not finish the job. Twenty-five
contigs, a sixth of the panel, collapse to strain-level identifications: their protein
distance was an artefact of a narrow protein index. Eighty-three stay distant on at least
one axis, 32 of them on both, and the floor is 65.0% nucleotide identity. The two axes
separate cases a single number merges: 49 contigs have their alignment confined to under
half their length, and 5 of those sit above 90% identity inside it — a conserved fragment
in an otherwise unmatched contig, which is the shape a recombinant or a multi-segment
assembly makes, not an identification. The worst of them covers 2.2% of the contig at
93.5% over 77 nt.

So the published numbers in this report are upper bounds on distance, as section 6's
blind spot said, and the measured size of that bound is now one in six at nucleotide
level. `blastx` against nr is running on the same instance and will put the same bound on
the protein numbers, which are the ones section 6's ordering rests on.

## Blind spots

- An ORF with no subject and no domain is unexplained, not non-coding. Past roughly 20-30%
  identity sequence search is at the end of its reach and the instrument that still works
  is structure, which is not attempted here.
- pyrodigal in meta mode calls ORFs on non-coding loci; an ORF is a prediction.
- The contigs are the study authors' assembly, so each is already a consensus and hides
  within-sample strain mixtures below about 3.8% divergence (`results/pair_ladder/`).
- The control arm shares the assembler and the library with the test arm but not its
  divergence, so it cannot control for what the assembler does to divergent sequence
  specifically.
- Every name in sections 1 to 6 came from RefSeq viral proteins plus Swiss-Prot, which are
  narrow, so every identity there is an upper bound on distance from what is known.
  Section 7 measures that bound against core_nt and finds it is one contig in six at
  nucleotide level; the matching protein pass against nr is not in yet, so the gene-level
  ordering in section 6 still carries the unmeasured version of this caveat.
- Section 7 ran `blastn -task dc-megablast`, which is more sensitive than megablast and
  less sensitive than the `blastn` task. The direction is fixed: a more sensitive search
  can only collapse more of the divergent arm, never less, so 16.8% is a floor on how much
  of the panel a full-size nucleotide index explains, not an estimate of it. The
  near-reference arm shows the task is not what is limiting there, since it recovers all
  135 contigs, but a near-reference arm cannot show that for divergent sequence.
- Section 7 keeps one best alignment per contig, so a contig whose segments match
  different subjects reads as the single best of them.

## Single seed

One seed, stated. The inputs are this sample's whole assembly and its whole read set, both
deterministic, so the seed selects nothing here but the pilot draws in the query builder.
Nothing in this report is a fraction of a sampled library.
