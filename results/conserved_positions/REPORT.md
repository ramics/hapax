# Which bits of a divergent protein are conserved

`experiments/conserved_positions.py`, seed 1, 300 divergent and 300 near-reference contigs
from SRR23973193, arms taken from `results/divergent_known/`. Result file:
[`conserved_positions_s1.json`](conserved_positions_s1.json). 1,317 ORFs called, 1,063 with
a significant protein alignment, 790 of those carrying a Pfam domain. Two contigs from
the holdout library ride along as labelled calibration rows. Every alignment column
walked from the DIAMOND traceback and assigned to a query position, then split by whether
that position falls inside a Pfam envelope at the curators' gathering cutoffs.

Two questions, two different answers.

## 1. Inside a protein: the conserved-core reading FAILS its own null

The expectation was that the columns inside a domain envelope hold while the rest drifts,
and the raw numbers look like that:

| arm | identity inside | identity outside | paired difference | share inside higher |
|---|---|---|---|---|
| divergent, 358 ORFs | 37.95% | 29.34% | **+7.60 points** | 83.5% |
| near-reference, 159 ORFs | 98.72% | 98.17% | +0.00 points | 49.7% |
| **dinucleotide-shuffled, 27 ORFs** | **45.39%** | **34.10%** | **+6.35 points** | **74.1%** |

The shuffled arm carries no information by construction, and it reproduces 84% of the point
difference and 89% of the directional share. Expressed as the ratio of mismatch rates
outside to inside, which has no ceiling on either arm, the three are indistinguishable:
**1.12 divergent, 1.12 near-reference, 1.14 shuffled.**

The mechanism is that both the envelope and the alignment are chosen to maximise a score.
hmmer places the envelope where the protein best matches its family and DIAMOND places the
alignment where the query best matches its subject, so the two agree on the same stretch
before any biology is involved. A shuffled contig that happens to produce an ORF with a
domain call and an alignment shows the same agreement.

**So no per-column or per-domain claim can be made from this instrument.** The ceiling
argument does not rescue it: the near-reference arm sits at 98% where there are two points
of room, which is why its difference in points reads 0.00, but its mismatch-rate ratio is
the same 1.12 as the divergent arm.

What does separate is, again, existence. **74.3% of real ORFs with a significant hit carry
a Pfam domain, against 5.1% of all shuffled ORFs** (34 of 43 shuffled ORFs that clear the
hit bar, out of 905 called). Whether a domain is there is calibrated; where it is, and how
well it matches, is not.

## 2. Between proteins on one contig: the polymerase holds and the maturation protein carries the distance

This one holds. Over 203 divergent contigs carrying two or more domain roles:

| role | most diverged on | least diverged on |
|---|---|---|
| maturation (receptor binding) | **118** | 7 |
| polymerase | 32 | **124** |
| capsid | 17 | 38 |
| other | 19 | 25 |
| protease | 8 | 4 |

**Maturation is the most diverged role on 58.1% of divergent contigs and the least diverged
on 3.4%, a 17-fold asymmetry. Polymerase runs the other way, 32 against 124.** Median
within-contig spread is 5.6 points.

By role, median whole-ORF identity on the divergent arm: maturation 28.75% over a median
300 aa, protease 28.8% over 1,237 aa, helicase 28.8% over 1,330 aa, polymerase 35.5% over
492 aa, capsid 37.25% over 146 aa. This reproduces `results/divergent_contigs/` on a pile
three times larger, where maturation read 63.5% against coat 77.9%.

**The symmetry argument is what makes it readable, because the near-reference control cannot
speak here.** Near-reference contigs agree across their roles to a median 0.8 points, so
which role comes out "most diverged" there is close to a coin toss: 12 other, 12 polymerase,
8 capsid, 7 movement, 6 maturation over 54 contigs. The obvious confound is alignment
length, since maturation aligns over a median 300 aa against polymerase's 492 and a shorter
alignment gives a noisier identity, which makes extremes more likely. That confound is
symmetric. It would put maturation at both ends, and maturation sits at one end 17 times
more often than the other.

## 3. The holdout genome lands on the same ordering

Two contigs from the holdout library were carried through the identical chain as labelled
rows. They are a different water, so they are not pooled into either arm.

`candidate_withheld`, 3,990 nt, covering 96% of a genome withheld from every index, true
distance 30.3% at nucleotide level:

| gene | Pfam | identity | alignment |
|---|---|---|---|
| maturation | Phage_mat-A | **53.5%** | 493 aa |
| polymerase | RNA_replicase_B | 66.5% | 620 aa |
| capsid | Levi_coat | **77.8%** | 135 aa |

Most diverged maturation, least diverged capsid, polymerase between them, spread 24.3
points. A genome at a distance nobody chose, in another library, falls on the Berlin
pattern. The coat protein is a third scorable gene the holdout thread's own table did not
carry.

`leak_already_referenced`, 6,383 nt from an already-referenced genome in the same library,
is the matched control and behaves: five roles between 99.6% and 100.0%, spread 0.4 points.
Which role is "most diverged" there is a coin toss, which is the same reason the Berlin
near-reference arm cannot test the asymmetry.

## 4. Reading

The functional answer to "which bits are conserved" is at the level of the gene, not the
column. On these contigs the replication machinery stays closest to what is published and
the protein that binds the host sits furthest away, which is the ordering host-range
pressure predicts, it is confirmed on the one genome whose true distance is known, and it is
the ordering a surveillance reader should be shown. The per-column version of the same
question does not survive the control that was run beside it, and is recorded here as killed
rather than left for someone to try again.

## 5. Blind spots

- A Pfam domain is a statement about a curated family, so the inside arm is biased towards
  well-studied families and the outside arm towards everything else. The direction is fixed
  and no control here removes it.
- The shuffled arm is 27 scorable ORFs. It is enough to kill a claim, since it reproduces
  most of the effect, and not enough to put an interval on how much.
- Role assignment is keyword matching on Pfam family names and descriptions. The first full
  run put 263 RNA_replicase_B domains, the Leviviricetes polymerase, into "other", which
  made the largest bucket on the divergent arm a mislabelled polymerase bucket. The word is
  in the list now; any role table from this instrument should be read against the family
  names behind it.
- The holdout rows are one contig each. They are a calibration point, not a sample.
- One best hit per ORF, so every identity is to the nearest relative the index holds, and is
  an upper bound on identity.
- 300 contigs per arm taken in file order from the divergent stage output, not sampled.
- One library, one date, plus two contigs from a second library.

## Cost

434 s on 4 cores: DIAMOND blastp across both arms and the null, Pfam at gathering cutoffs,
ORF calling. Two bugs were caught before anything was believed. DIAMOND emits `qseq` and
`sseq` ungapped, so zipping them walks out of register at the first indel and reported 5%
per-column identity on ORFs whose whole-ORF identity was 21%; the walk uses the `btop`
traceback instead. And every Pfam envelope and every alignment is now asserted to lie inside
its translated ORF, so a minus-strand gene cannot place its columns in the wrong frame.
