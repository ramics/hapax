# How far each method reaches, and what the shortcut costs

`experiments/reach_queries.py` builds the panel, `experiments/reach_arms.py` runs the arms on
one EC2 instance, `experiments/report_reach.py` generates [`TABLES.md`](TABLES.md) from the
rows. Seed 1 for the panel, five seeds for the null. Two runs: `i-0672a2ca4a782aa21`, which
ran five arms and found three bugs, and `i-0499bf15b019eba9f`, which ran all nine. Both on
c6i.2xlarge on demand.

Imo, 2026-11-22: *"sounds good, let's try it. structure comparison also sounds fun."*

## The question, and the one thing it needs held fixed

For a contig that sits far from anything published, which search method recovers the most
distant real relative, and at what cost?

Reach is a property of an instrument **and** of the set it searches, and in this project the
set dominates by orders of magnitude: 1.4 MB of this sample's own contigs removed more of the
library than 272 MB of SILVA. So a comparison of instruments is only a comparison of
instruments if the subjects are held fixed, and where that is possible here, they are:

| alphabet | subjects | arms |
|---|---|---|
| protein sequence | Swiss-Prot + RefSeq viral proteins, 1,210,551 sequences | DIAMOND x3, BLASTP, exhaustive Smith-Waterman, jackhmmer |
| protein structure | AlphaFold/Swiss-Prot | Foldseek |
| profile library | Pfam-A 39.0, 53,275 families | hmmsearch, HHsearch |
| nucleotide genomes | RefSeq viral genomic release | tblastn, codon scoring |

**AlphaFold's Swiss-Prot set is Swiss-Prot folded.** That is what makes the structure arm a
controlled comparison rather than a different database: the same subjects, in a different
alphabet. Pfam, PDB and the viral genomes are separate universes, so no difference against
the sequence baseline is computed for them and their rows say so.

## The panel, and the two controls that decide every row

| pile | what | n |
|---|---|---|
| test | the near-complete divergent *Fiersviridae* genomes of `results/divergent_contigs/` §3 — 3.0-4.5 kb, 116-467x mapped depth, polymerase 44.9-84.8% | 7 |
| ground truth | the withheld-genome contig of `results/holdout_cascade/` — 3,990 nt at ~452x, true nucleotide divergence 30.28% | 1 |
| zero distance | its already-referenced companion from the same library, 6,383 nt at 1,188x | 1 |
| negative | near-reference contigs from the same assembly, same family, same length band, 95.2-99.3% identity | 7 |
| null | every record above, dinucleotide-shuffled, five seeds | 80 |

**The negative arm is matched on family as well as on length**, because a plant virus at 98%
identity has a different gene architecture and would make the control unreadable. Every
control contig carries two or more of `RNA_replicase_B`, `Levi_coat` and `Phage_mat-A`, as
every test contig does.

The null is `holdout_run.dinuc_shuffle`, the fixed Altschul-Erikson walk, **imported and not
reimplemented**: two hand-rolled walks in this repo answered a dead end with `return seq`,
put real sequence into the null pile and manufactured two kills in two days
(`results/fail_loudly_audit/`). The fallback rate is a property of the input as much as of
the code, so it is measured here too: **0 of 16 records at every one of the five seeds.**

A dinucleotide shuffle preserves stop-codon density, so the shuffled records yield few and
short ORFs — 36 against the real pile's 54 at one seed. The null therefore calibrates
**existence** of a significant alignment and nothing conditional on one. That is the same wall
`results/divergent_known/` §2 reached from the other side.

## The answer to the first two arms: nothing beats the production search

| arm | contigs where it found a nearer relative than the baseline |
|---|---|
| DIAMOND `--sensitive` | 0 of 16 |
| DIAMOND `--very-sensitive` | 0 of 16 |
| BLASTP, composition-based statistics | 0 of 16 |
| **Smith-Waterman, exhaustive, no seeding** | **0 of 16** |

Three sensitivity settings, an aligner with exact Karlin-Altschul statistics, and full
Smith-Waterman against all 1,210,551 subjects with no seed at all return **the same nearest
relative for every query on the panel**. The per-contig table is [`TABLES.md`](TABLES.md) §8.

**So the binding constraint is the significance statistics, not the seeding.** That was the
expectation this arm was built to test cheaply, and it is confirmed: removing the heuristic
entirely moves nothing. The exhaustive arm cost 550 s against DIAMOND's 94 s for three modes
and BLASTP's 40 s — affordable, and it buys nothing on this material.

### What the exhaustive arm does produce, and why it is not reach

It reports 595 query-subject pairs the baseline did not, which looks like reach until the
control is read: **273 of them are on the divergent contigs and 260 on the near-reference
contigs.** BLASTP's 283 extras split 124 to 128 the same way. A contig that already matches
its reference at 96% picks up just as many, so these are the tail of marginal alignments
rather than a longer reach, and none of them is the nearest relative for its contig.

## The statistic in the brief does not work, and the control is how we know

"The most distant relative recovered" was the headline column this comparison was designed
around. At E < 1e-5 over at least 60 aligned residues it reads:

| pile | median most-distant identity |
|---|---|
| divergent, polymerase 44.9-84.8% | 21.8% |
| near-reference, 95-99% | 22.2% |
| the withheld genome, 30.28% diverged | 21.1% |
| its already-referenced companion | 20.7% |

**A contig matching its reference at 99.8% has the same most-distant relative as one at
44.9%.** At a fixed E-value the furthest alignment a 1.2-million-sequence database will
return is set by that database's size and the query's length — a longer query clears a fixed
E-value on a worse alignment — so the column describes the operating point and not the query.

That is the control moving with the test, which is this repo's standing tell for a field
describing the instrument rather than the material (CLAUDE.md, "the tell in all three was the
control moving with the test"). The column stays in the tables because its flatness across
piles is the result. **Nothing is ranked on it.** What separates piles is the *nearest*
relative's identity with its coverage beside it; what would separate arms is finding a
significant relative where another arm found none, and on this panel no arm does.

## The panel reproduced exactly, which almost nothing here can claim

CLAUDE.md: *"Having actually re-run a result from clean and got the same numbers is a
stronger claim, and almost nothing in this repo can make it."*

Different code, a fresh container, databases staged fresh from the pinned NCBI snapshot and
from EBI:

| contig | published | this run |
|---|---|---|
| Kahshuvirus borborenecus | 44.9% | 44.9% |
| Fiersviridae | 48.6% | 48.6% |
| Koteshevirus lutadaptatum | 46.6% | 46.6% |
| Sholavirus sp. 'caenicola' | 83.8% | 83.8% |
| Rodtovirus pelovivens | 45.6% | 45.6% |
| Dugnivirus solivicinum | 47.8% | 47.8% |
| Chimpavirus luticola | 84.8% | 84.8% |
| the withheld genome | 66.5% over 620 aa | 66.5% over 620 aa |

Seven polymerase identities from `results/divergent_contigs/` and the holdout's protein
identity from `results/holdout_cascade/`, to the decimal place.

## What the first run bought, which was three bugs

The first run was a shakedown and it is reported because it is where most of the correctness
came from. Five arms ran; the two structure arms failed; three errors were found that would
each have produced a plausible wrong table.

**RefSeq names its genomic release files differently from its protein files** —
`viral.N.1.genomic.fna.gz` against `viral.N.protein.faa.gz`. Eight 404s dropped the codon arm
from that run with a named gap rather than a wrong answer, which is what the graded staging
is for.

**HMMER's two coordinate pairs are not the same way round for both programs.** `--domtblout`
puts "hmm coord" on the query and "ali coord" on the target, and which of those is our ORF
depends on the program: hmmsearch searches a profile library against our sequences, so the
ORF is the target; jackhmmer searches our sequence against a database, so the ORF is the
query. One mapping for both put the Pfam profile's length into the ORF's `qlen` column. In
57,917 profile rows **the only visible symptom was a single query coverage of 1.156.**
Coverage above 1.0 is now an assertion, both programs go through one parser that takes which
side the ORF is on, and `tests/test_reach_domtbl.py` holds that run's own output as its
fixture — including a 1,132 residue ORF against a 173-position profile, where every quantity
differs so a swap cannot pass by coincidence.

**Foldseek's `--prostt5-model` wants the directory, not the weights file inside it.** It
accepted the option, echoed it back and exited 1, which cost both structure arms.

The codon walk was tested before it ever ran, in `tests/test_reach_codon.py`, for the same
reason: it reconstructs a codon alignment from a gapped protein alignment, a CDS in
nucleotides and subject coordinates that can run backwards, and an off-by-one in any of them
gives dN/dS values that are wrong and entirely plausible. The four cases are identical
sequences, third-position changes only (which a frame off by one turns non-synonymous), the
same alignment with the genome flipped, and a saturated pair that must come back undefined
rather than clipped.

## Blind spots

- **One panel, sixteen contigs, one family.** Every test contig is an ssRNA bacteriophage
  genome with a polymerase, a coat protein and a maturation protein. "No arm beats the
  baseline" is measured on material where the baseline already finds a full-length polymerase
  alignment at 45-85% identity. It says nothing about a contig with no recognisable gene at
  all, which is the tier-2 pile `results/residue_contig_annotation/` describes.
- **The subject sets are narrow by choice.** Swiss-Prot plus RefSeq viral is what made the
  arms comparable; it is not what a production search would use. Every identity here is an
  upper bound on identity and every divergence a lower bound on divergence, and
  `results/divergent_contigs/` §7 has already measured that a broader nucleotide index
  collapses 16.8% of the divergent panel to strain level.
- **One seed on the real arms.** The real contigs are fixed sequence, so there is nothing to
  re-draw; the five seeds are on the null, which is the only stochastic element. A five-seed
  range on the test arm would be a range over nothing.
- **The null cannot calibrate anything conditional on a hit**, because it contains no
  significant alignments at all. Identity, coverage and aligned length are reported for what
  they are.
- **The exhaustive arm is exhaustive in alignment, not in scoring.** BLOSUM62 at 11/1 with
  published lambda and K: a different matrix is a different instrument and none was tried.
- **Timings are from one instance type on one day.** Fixed and marginal cost are not
  separated here as `results/bigindex/` separates them, because every arm's query set was the
  same size.
