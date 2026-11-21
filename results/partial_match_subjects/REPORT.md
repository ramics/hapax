# What the residue's partial matches actually match

`experiments/partial_match_subjects.py`, one run over 1,133 core_nt rows (the 500-read pilot
plus three 211-read blocks from `i-019ad3b008c57a5d1`), 567 of which match ≥30 bases over
less than 80% of the read. Mate 1 only, because the query sets are mate 1 only.

## The question, and why it is answered here without a paid run

Imo asked for pairwise alignments of residue reads that match part of a reference but were
not identified — the whole read, the reference's matching nucleotides, and enough reference
past the match to equal the read's length — to eyeball what the shape is. His stated
suspicion: *"you're comparing apples with pairs and what you're seeing is just the
untranscribed portions of the blast."*

The rows this project holds cannot be rendered directly, because the pilot's BLAST payload
carries no subject coordinates. A re-run with the alignment fields is a paid launch. But two
things turned out to be free, and between them they answer the question:

1. **The subject inventory needs no alignment at all.** The accession is in the row, and one
   batched esummary call turns 543 accessions into 543 organism names. What a read matched is
   often the entire answer.
2. **The alignment can be rebuilt locally where the reference is small enough to fetch.**
   That worked for 13 of them — and the rate at which it *failed* is itself the result.

## 1. His concern was legitimate, and the control says so

The block files score residue reads beside reads our cascade already identified, so the same
partial-match question can be asked of material we know is known. Median figures per pile:

| pile | band | n | median E | median identity | median matched bases | median fraction of read |
|---|---|---|---|---|---|---|
| residue | significant (E<1e-5) | 20 | 4.9e-10 | 85.8% | 74.0 | 0.724 |
| identified | significant | 8 | 4.7e-08 | 81.9% | 73.5 | 0.707 |
| **residue** | **chance (E≥1e-5)** | **490** | **0.82** | **87.8%** | **38.0** | **0.3645** |
| **identified** | **chance** | **49** | **0.81** | **87.2%** | **39.0** | **0.3645** |

**In the chance band the two piles are identical to three figures.** Same E-value, same
percent identity, same match length, same fraction of the read. 95% of this population
(539 of 567) sits in that band. A separator that reads the same on the residue and on known
material is describing the database, not the reads.

That the medians converge is expected *if and only if* these are chance matches: both piles
take the best hit per read from the same database at the same cutoff, so if neither pile has
real relatives, both medians land on the chance distribution's median by construction. The
convergence is therefore confirmatory rather than surprising — and it is the shape of the
failure this repo keeps finding, where the control moves with the test.

The identical `0.3645` was checked rather than trusted: it is 39/107, the modal read length
with the modal match length, and the values either side of each median are spread
(0.3619, 0.3645, 0.3689, 0.3738), so it is a genuine coincidence of two distributions with
the same mode and not a duplicated column.

## 2. The subject names alone are damning

In the chance band the most common subjects of a Berlin wastewater read's "partial match"
include: `Caudoviricetes sp.` (11), `Acidobacteriota bacterium` (9), `uncultured bacterium`
(7), `Verrucomicrobiia bacterium` (5), `Actinomycetes bacterium` (5), `Homo sapiens` (4) —
and, in the pilot's own list, *Timema monikensis* (a Californian stick insect), *Apteryx
mantelli* (kiwi), *Pan troglodytes*, *Callithrix jacchus*, *Gallus gallus*, *Mus musculus*,
*Oryza sativa* and *Scomber australasicus* (mackerel), the last several as BAC clones and
whole chromosomes. A 107 nt wastewater read sharing 39 bases at 88% identity with a kiwi BAC
at E=0.8 needs no picture.

## 3. 94% of these hits cannot be reproduced by a short-read aligner

Of the 231 partial matches whose reference is at or under 1 Mb and could therefore be
fetched, **218 could not be aligned locally at all** — not by minimap2's short-read preset,
not at k=11, not at k=9, against the very reference BLAST named. Only 13 produced an
alignment, and the seeding tier they needed is a ladder of weakness: 1 at the default
preset, 5 at k=11, 7 only at k=9.

Across the padded window — a window the width of the read, centred on the reported match —
the aligned 13 hold a median identity of **0.617**, and only **7.7%** exceed 0.8.

"BLAST reports it and a standard short-read aligner cannot seed it" is the sharpest available
statement that a hit is not a relative.

## 4. The significant band is a different animal, and it is the interesting one

The 28 significant partial matches do not look like the chance band at all. They match over
~74 of ~107 bases, and **the identity holds into the padding**:

```
SRR23973193.80670/1 [residue] E=2.21e-14 pident=84.3 blast_aln=83/106
    | BK028722.1 Caudoviricetes sp. (124155 nt)
  read TACTGTCTGAAGGGTCTCAACTATTTTGTAGGTCCTAATGGTGCAGGGAAATCTACTGTTCTACAAGCTGTT...
       |||   || ||   | | || ||||||||||||||||| |||||||| |||||||||||  | || ||  ||
  ref  TACAAACTCAATAATATGAATTATTTTGTAGGTCCTAACGGTGCAGGTAAATCTACTGTAATGCAGGCAATT...
  identity over the padded window 0.764
```

That is not a chimaera and not noise: it is a phage protein-coding region about 16% diverged
from the nearest thing in core_nt, matching along its whole length at a rate a nucleotide
search can barely hold. Likewise `PP215199.1 Penaeus vannamei levi-like virus 2` at 87.5%
over 72 bases with padded identity 0.738.

**So requirement 1.ii's "part of the read matches and the rest does not" splits in two.** At
E<1e-5 it is not really happening — the match covers 72% of the read and keeps going into
the padding, so these are divergent relatives, which is requirement 2's material rather than
1.ii's. Below significance it is happening constantly and means nothing. This agrees with
`results/subtract_known/` §3.4, reached there from the opposite direction, and with the
chimaera result in `results/residue_structure/`: the shape has no population at read level.

## 5. One genuinely odd case worth keeping

```
SRR23973193.5923/1 [identified] E=2.1e-08 pident=97.6 blast_aln=42/107
    | EU195987.1 Acrobeloides maximus (3483 nt)
  identity over the padded window 0.617, local match at read positions 56-98
```
42 bases at 97.6% identity to a nematode record, sitting in the middle of the read with
gibberish either side. This is the one shape in the file that looks like 1.ii as originally
imagined, it is in the **identified** pile, and there is one of it. Worth remembering when
the artefact thread opens, which `results/README.md` already flags as the weakest
load-bearing piece.

## What this cannot say

- **Mate 1 only.** The query sets are mate 1 only, so no pair statement is available and per
  `CLAUDE.md` none may be inferred.
- **One run, not five seeds.** The 1,133 rows are one block series from one instance; the
  pilot half is a genuine second staging of the same pinned database and reproduces it, but
  this is not a five-seed result and the medians above should not be quoted as one.
- **The aligned subset leans viral.** The 1 Mb fetch cap correlates with organism, because
  large genomes are cellular, so §3's 13 alignments are not a random sample of §1's 567.
- **Organism names are NCBI's annotation** and travel through unchecked. A mislabelled
  record would arrive here looking like a finding, which is a reason to read the alignments
  and not only the names.
- **The 49 identified reads in the chance band are a thin control.** The medians agree to
  three figures, which is strong, but the pile is small; more blocks from the running pass
  will thicken it.

## Files

- `inventory_s1.tsv` — all 567, with organism, reference length, identity, matched bases,
  E-value and band. Sorted significant first.
- `partial_matches_s1.aln.txt` — the 13 renderings, read against reference plus padding.
- `partial_matches_s1.fa` — the same as FASTA pairs.
