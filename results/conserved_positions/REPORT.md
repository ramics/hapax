# Which bits of a divergent protein are conserved

`experiments/conserved_positions.py`, seed 1, 300 divergent and 300 near-reference contigs
from SRR23973193, arms taken from `results/divergent_known/`. Result file:
[`conserved_positions_s1.json`](conserved_positions_s1.json). 1,317 ORFs called, 1,063 with
a significant protein alignment, 811 carrying a Pfam domain. Two contigs from the holdout
library ride along as labelled calibration rows. Every alignment column is walked from the
DIAMOND `btop` traceback and assigned to a query position, then split by whether that
position falls inside a Pfam envelope at the curators' gathering cutoffs.

**Operating point:** E < 1e-5 and at least 60 aa aligned, against RefSeq viral proteins plus
Swiss-Prot; Pfam-A at gathering cutoffs. The null is a dinucleotide-preserving shuffle of the
same contigs through the identical chain, from `experiments/shuffle_control.py`.

Two questions, two different answers, and the first of them is an honest *cannot tell*.

## 1. Inside a protein: the per-column reading is untested, not killed and not confirmed

The expectation was that columns inside a domain envelope hold while the rest drifts, and
the raw numbers look like that:

| arm | identity inside | identity outside | paired difference | share inside higher |
|---|---|---|---|---|
| divergent, 358 ORFs | 37.95% | 29.34% | **+7.60 points** | 83.5% |
| near-reference, 159 ORFs | 98.72% | 98.17% | +0.00 points | 49.7% |
| **dinucleotide-shuffled, 0 ORFs** | — | — | — | — |

**The shuffled arm is empty.** Of 932 ORFs called on shuffled contigs, 0 produce a
significant alignment, so there is no shuffled column distribution of any kind. An earlier
version of this section killed the per-column reading on a shuffled arm that reproduced 84%
of the effect. That arm was not shuffled: the shuffle returned its input on a dead-ended
walk, and every one of those 27 "shuffled" ORFs came from a contig that was never shuffled.
**The kill is withdrawn.** It is not replaced by a pass, because an empty control supplies no
argument in either direction.

What is calibrated is the step before the alignment. **Pfam runs on every ORF in both arms
and does not need a BLAST hit: 811 of 1,317 real ORFs carry a domain at gathering cutoffs,
against 0 of 932 shuffled ORFs.** Whether a domain is there is therefore a calibrated
statement about real sequence. Where it is, and how well the columns inside it match, is not.

The near-reference arm says a little more than it did. Read as the ratio of mismatch rates
outside to inside, which has no ceiling on either arm, the divergent arm reads a median 1.12
with **83.2% of its 358 ORFs above 1**, and the near-reference arm a median 1.12 with
**52.7% of its 112 scorable ORFs above 1** — a coin toss. So the mechanism the kill rested
on, that hmmer places an envelope and DIAMOND places an alignment by maximising the same
kind of score and therefore agree before any biology, does not by itself produce the
asymmetry: it operates on the near-reference arm too, where the asymmetry is absent.

That is suggestive and it is not enough. At 98% identity a near-reference ORF has one or two
mismatches on each side of the envelope boundary, so its ratio is a coin toss for want of
mismatches rather than for want of an effect, and the already-referenced holdout contig makes
the instability visible by reading 4.85 off a single mismatch. **So this instrument records
"cannot discriminate here".** The per-column claim is open, it needs a control that produces
alignments — a length-matched real-but-unrelated subject, or a reference-free structural
reading — and nothing in this run licenses a per-column or per-domain conclusion either way.

## 2. Between proteins on one contig: the polymerase holds and the maturation protein carries the distance

This one holds, and it survived two faults found while re-reading it.

**Paired, role against role, on contigs carrying both.** A median taken per role compares
different sets of contigs, because the role counts here are lopsided. Paired on the divergent
arm, excluding contigs where one ORF carries both roles:

| pair | contigs | first role more conserved | median difference |
|---|---|---|---|
| maturation vs polymerase | **131** | 6.9% | **−6.9 points** |
| capsid vs maturation | 31 | 90.3% | +10.7 points |
| capsid vs polymerase | 63 | 58.7% | +2.6 points |
| maturation vs other | 17 | 5.9% | −8.2 points |
| capsid vs other | 13 | 84.6% | +7.6 points |

So maturation sits below polymerase on 122 of 131 contigs that carry both, and below capsid
on 28 of 31. Capsid against polymerase is close to a coin toss. Counted by position rather
than paired, the same thing: maturation is the most diverged role on 118 of 203 divergent
contigs and the least diverged on 7, a 17-fold asymmetry, with polymerase running the other
way at 32 against 124. Median within-contig spread 5.6 points on the divergent arm against
0.8 on the near-reference arm.

**Two roles on one ORF are not two genes.** A picorna-like polyprotein carries protease,
helicase and RdRp in a single ORF, and this script credits that ORF's identity to every role
its domains name, so those pairs are identical by construction. That is why
helicase-against-polymerase read a median difference of exactly 0.0 points and a
share-higher of 0.000 over 31 contigs: one number compared with itself, on 28 of the 31.
Those pairs are now counted and excluded, which removes four of the fourteen pairs entirely
and leaves the maturation-polymerase comparison untouched at 131 of 131 distinct-ORF
contigs.

**The two confounds the literature predicted were tested and neither explains it.**

- *Reference abundance.* Detection models for this family are keyed on the conserved
  replicase, so the polymerase could look conserved because the database holds more of them.
  Measured as distinct subject proteins matched per ORF, the two roles are the same:
  polymerase 0.604 over 275 ORFs, maturation 0.598 over 132.
- *Alignment length.* Maturation aligns over a shorter stretch than polymerase, and a
  shorter alignment gives a noisier identity. Inside matched length bands the gap holds:
  at 150-300 aa, polymerase 37.2% against maturation 29.1%; at 300-600 aa, 35.6% against
  28.3%. Seven to eight points in both bands.

**The near-reference arm cannot test any of this, and that is itself the calibration.** Its
roles agree to a median 0.8 points and its paired comparisons sit near 0.5 — maturation
against polymerase 0.273 over 11 contigs, median −0.4 points — so on material that matches
its reference, which role comes out most diverged is close to a coin toss.

By role, median whole-ORF identity on the divergent arm: maturation 28.75% over a median
300.5 aa, protease 28.8% over 1,237 aa, helicase 28.8% over 1,329.5 aa, polymerase 35.5% over
492 aa, capsid 37.25% over 146 aa. **Capsid cannot be placed in the ordering.** Half its
ORFs align over under 150 aa, where identity is quantised and optimistic, and in the one
length band where it can be compared directly it falls the other way, reading 27.0% at
300-600 aa over 7 ORFs. The claim here is the maturation-against-polymerase end.

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
pattern at the maturation end. Its capsid is the most conserved of the three, which the
Berlin panel does not reproduce and which the capsid's 135 aa alignment is enough to
explain either way.

`leak_already_referenced`, 6,383 nt from an already-referenced genome in the same library,
is the matched control and behaves: five roles between 99.6% and 100.0%, spread 0.4 points.

## 4. The ordering is published; this use of it is not

Friedman, Cooper, Bollback et al., *J. Virol.* 83(21):11233-43, 2009
([PMC2772794](https://pmc.ncbi.nlm.nih.gov/articles/PMC2772794/)), mapped 30 complete
Leviviridae genomes. Quotations below are from the full text, checked against the page
rather than a search snippet; the per-gene numbers appear only in running prose, and the
paper marks several as "data not shown".

**Two different quantities, and they give two different rankings in the same paper.**

- *Within a genogroup*, amino acid distance by ClustalW and neighbour joining, the capsid is
  most conserved in all four: genogroup I "The capsid protein was the most conserved
  (distance of 0.0000 to 0.0411)", genogroup II "the capsid protein was the most conserved,
  with a genetic distance of 0.0135 to 0.1116", and likewise III and IV. The greatest
  distance falls on the maturation protein in genogroups III and IV — "The greatest genetic
  distance was observed in the maturation protein (0.0607 to 0.5646)" — but not in I or II.
- *Across genogroups*, by nucleotide similarity plot, "the replicase genes were most similar
  (approximately 0.5 to 0.6, or 50 to 60% similarity)" while "the assembly or maturation
  genes were the most dissimilar (<10%)".

Our quantity is the distance from a contig to the nearest reference the database holds,
which is a between-lineage comparison, and the paper's between-lineage readout is the second
one. **It matches ours on both ends: replicase closest, maturation furthest.** The
capsid-most-conserved ranking is variability among known close relatives, measured at a
distance where our divergent arm has no members, and it is the ranking our one holdout
genome reproduces. The paper also notes the capsid's length is tightly constrained at "130
to 133 amino acids", which is the same fact that makes its alignment short in our data.

So the ordering is precedented and should be cited rather than presented as new. What
appears not to be precedented is the use: novelty in this literature is scored on a single
marker or a whole-genome best hit, and within-genome gene-wise divergence *spread* has not
been proposed as a ranking signal. The closest adjacent work is gene-tree discordance as
information and per-protein novelty fractions in large-scale virus surveys. That split —
ordering known, use new — is the stronger claim and the one to make.

## 5. Reading

The functional answer to "which bits are conserved" is at the level of the gene, not the
column. On these contigs the replication machinery stays closest to what is published and
the protein that binds the host sits furthest away. That is the ordering host-range pressure
predicts, it is what the 2009 between-genogroup comparison reports, it is confirmed on the
one genome whose true distance is known, it survives matching on alignment length and on
reference abundance, and it is the ordering a surveillance reader should be shown. The
per-column version of the same question is recorded as untested: its kill is withdrawn and
its null is empty.

## 6. Blind spots

- A Pfam domain is a statement about a curated family, so the inside arm is biased towards
  well-studied families and the outside arm towards everything else. The direction is fixed
  and no control here removes it.
- The shuffled arm produces no alignments at all, which calibrates existence and nothing
  else. Every field conditional on an alignment is uncalibrated in this run, including the
  per-column ones, and that is stated rather than converted into a kill.
- Role assignment is keyword matching on Pfam family names and descriptions. The first full
  run put 263 RNA_replicase_B domains, the Leviviricetes polymerase, into "other", which
  made the largest bucket on the divergent arm a mislabelled polymerase bucket. The word is
  in the list now; any role table from this instrument should be read against the family
  names behind it.
- The maturation-polymerase result is 131 contigs from one library, and this family is
  single-stranded RNA bacteriophage, so it does not generalise to other virus architectures
  without being measured there.
- The holdout rows are one contig each. They are a calibration point, not a sample.
- One best hit per ORF, so every identity is to the nearest relative the index holds, and is
  an upper bound on identity. `results/divergent_contigs/` §7 bounds that at nucleotide level
  against core_nt; the protein equivalent against nr has not landed.
- 300 contigs per arm taken in file order from the divergent stage output, not sampled.
- One library, one date, plus two contigs from a second library.

## Cost

471 s on 4 cores: DIAMOND blastp across both arms and the null, Pfam at gathering cutoffs,
ORF calling. Four faults were caught before anything was believed, three of them after a
version of this report had already been written. DIAMOND emits `qseq` and `sseq` ungapped,
so zipping them walks out of register at the first indel and reported 5% per-column identity
on ORFs whose whole-ORF identity was 21%; the walk uses the `btop` traceback instead. Every
Pfam envelope and every alignment is asserted to lie inside its translated ORF, so a
minus-strand gene cannot place its columns in the wrong frame. The shuffle returned its input
on a dead-ended walk, which is what produced the withdrawn kill in section 1. And the
mismatch-rate ratio took its median over the ORFs where it is defined and its share over all
of them, which made the near-reference arm read a median of 1.12 with 37.1% above 1 —
arithmetically impossible, and the specific number the kill leaned on.
