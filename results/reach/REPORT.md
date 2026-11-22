# How far each method reaches, and what the shortcut costs

`experiments/reach_queries.py` builds the panel, `experiments/reach_arms.py` runs the eight
arms that need an instance, `experiments/reach_structure.py` runs the structure arm in the
container, and `experiments/report_reach.py` generates [`TABLES.md`](TABLES.md) from both
sets of rows. Seed 1 for the panel, five seeds for the null.

**The table is merged from two sources and [`TABLES.md`](TABLES.md) prints each one's
provenance separately**, because one script hash covering work from two hosts is exactly what
this repo's `script_sha256` convention exists to prevent. The eight paid arms ship from
`i-055ea9ae3eaf1809c`; the structure arm runs in the container against 178,660 solved virus
chains, for nothing. They agree on the genes by construction and the manifests prove it: the
structure arm records the sha256 of the `orfs.faa` it called, and it is byte-identical to the
instance's.

Seven runs on c6i.2xlarge on demand, about $3 in total. Five were shakedowns, they found ten
faults between them, and they are described below rather than discarded: five of the ten were
in readouts rather than in tools, and every one of those produced a plausible wrong table.
The largest single line in that $3 is not a shakedown but a staging download that stalled
dead and took a run's whole dead-man timer with it, which has its own entry.

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
| test | the near-complete divergent *Fiersviridae* genomes of `results/divergent_contigs/` §3 — 3,058-4,526 nt, 116-467x mapped depth, polymerase 44.9-84.8% | 7 |
| ground truth | the withheld-genome contig of `results/holdout_cascade/` — 3,990 nt at ~452x, true nucleotide divergence 30.28% | 1 |
| zero distance | its already-referenced companion from the same library, 6,383 nt at 1,188x | 1 |
| negative | near-reference contigs from the same assembly, same family, same length band (3,118-4,047 nt), polymerase 95.2-99.2% | 7 |
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
short ORFs — 21, 29, 30, 31 and 37 across the five seeds, against the real pile's 54. The
null therefore calibrates **existence** of a significant alignment and nothing conditional
on one. That is the same wall `results/divergent_known/` §2 reached from the other side.

## The answer to the first two arms: nothing beats the production search

The baseline is DIAMOND `--ultra-sensitive`, the most sensitive setting the production search
offers; every other arm is scored against it.

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
entirely moves nothing. The exhaustive arm cost 902 s against DIAMOND's 93 s for all three
modes and BLASTP's 53 s — affordable at this panel's size, and it buys nothing on this
material.

### What the exhaustive arm does produce, and why it is not reach

It reports 595 query-subject pairs the baseline did not, which looks like reach until the
control is read: **273 of them are on the divergent contigs and 260 on the near-reference
contigs.** BLASTP's 283 extras split 124 to 128 the same way. A contig that already matches
its reference at 96% picks up just as many, so these are the tail of marginal alignments
rather than a longer reach, and none of them is the nearest relative for its contig. The
furthest of the 595 is 19.3% identity over 450 residues — on the floor described in the next
section, which is where a database of this size puts its last significant alignment
regardless of the query.

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

**And the table reproduces itself.** Runs 2 and 3 were separate instances, separate
containers and separately staged databases, a little under an hour apart. Regenerating
[`TABLES.md`](TABLES.md) from each gives files that differ in exactly two places: the
provenance line naming the run, and the wall-clock column. Every identity, every count,
every control, every saturation figure is identical. The timings are the part that does not
reproduce — jackhmmer ran 1,008 s and then 1,150 s on the same input — which is the ordinary
reason this repo's convention says a benchmark number needs its instance type and its day.

## Profiles: the literature's claim reproduces, but only profile against profile

Two profile arms, and they disagree in a way that matters.

**hmmsearch**, a curated profile against our sequence, places a domain on every contig in
every pile, at the gathering cutoffs and at E ≤ 10 alike, and puts nothing on the null. It
agrees with the sequence arms about which contigs are recognisable, which is a cross-check
worth having and is not reach: there is no contig profiles find and sequence search misses.

**HHsearch**, our profile against Pfam's profiles, places about three times as many models —
and that is where the family literature's claim shows up:

| arm | pile | ORFs with a model | ORFs carrying 2+ models | most on one ORF | distinct models |
|---|---|---|---|---|---|
| hmmsearch at GA | divergent | 17 | **0** | 1 | 3 |
| hmmsearch at GA | near_reference | 17 | 1 | 2 | 4 |
| hmmsearch at GA | its companion | 4 | 1 | 4 | 7 |
| HHsearch | divergent | 19 | **12** | 6 | 9 |
| HHsearch | near_reference | 20 | 13 | 6 | 9 |
| HHsearch | the withheld genome | 3 | 2 | 6 | 9 |
| HHsearch | its companion | 4 | 3 | **10** | 22 |

Callanan 2020 and Rūmnieks 2026 report that several profiles per protein are needed to
recover this family's divergent members. **Profile against sequence does not reproduce that
at all** — not one divergent ORF carries two models. **Profile against profile does**: the
median divergent ORF carries two and the maximum is six, and the models it adds are the
other viral polymerase families, PF00680, PF00978, PF00998, PF02123 and PF05919. So the
claim is about the comparison being profile-to-profile, not about profiles in general, and
reading it off hmmsearch would have produced a flat contradiction of the literature.

**It is still not reach, and the control is again what says so.** HHsearch is very nearly a
superset of hmmsearch — zero hmmsearch-only contig-model pairs on three of the four piles,
one on the fourth — and the pairs it adds number **34 on the divergent pile and 34 on the
95–99% control.** Identical. A contig matching its reference at 96% picks up exactly as many
extra models as one at 45%, so the multiplicity is a property of Pfam's overlapping RdRp
models rather than evidence about the query. And HHsearch's best hit on the withheld genome
is PF00978 at E = 3.9e-06, against hmmsearch's PF03863 at E = 1e-25: on the one query whose
distance is known, the more sensitive instrument returns the *weaker* statistic.

The one row that does stand apart is the already-referenced companion, at a median of eight
models per ORF and a maximum of ten against everything else's two and six. That is the
familiarity polarity `cov_step_conservation/` found in the coverage-shape family, one level
up: **a count of matched profiles says "already described".** It cannot be used the other
way round.

**The null here is zero for a reason that is not a clean null.** HHsearch reads 0 of 80
shuffled records, but no shuffled record could have scored: the query profile is built from
jackhmmer's alignment, every one of the 148 shuffled ORFs yields an alignment of depth 1,
and `hhmake` refuses a profile built from one sequence. The arm's null is empty by
construction, not by calibration, so its only working control is the near-reference pile.
That limitation belongs to any profile-based novelty claim, not just to this panel.

## Codon space: the saturation claim, measured on both bands

The brief asked for the codon arm on the divergent and the near-reference band both, so that
saturation would be measured rather than asserted. Nei-Gojobori counting with a
Jukes-Cantor correction that returns undefined instead of clipping, over tblastn alignments
to the RefSeq viral genomic release, taking each ORF's best-scoring alignment:

| pile | median %id | dS undefined | median dS | median dN | median dN/dS |
|---|---|---|---|---|---|
| near_reference | 96.1 | **0 of 22** | 0.189 | 0.019 | 0.104 |
| divergent | 52.3 | **11 of 22** | 1.164 | 0.466 | 0.109 |
| the withheld genome | 66.9 | 0 of 3 | 1.333 | 0.271 | 0.204 |
| its companion | 99.8 | 0 of 4 | 0.002 | 0.001 | 0.000 |

**dN/dS is the same on both bands — 0.104 against 0.109 — and that is not the result.** The
result is the column beside it: half the divergent pile's genes have no measurable
synonymous divergence left at all, and the half that do sit at dS ≈ 1.16, six times the
control's and close enough to saturation that the ratio's denominator is barely carrying
information. So codon-level selection scoring works on the 95–99% band, degrades through
the 67% band, and has lost half its measurements by the time the nearest relative is at
52%. It orders the panel correctly — the companion at 99.8% reads dS 0.002 — and it answers
nothing a protein alignment did not already answer, because it needs the protein alignment
to exist first.

**This arm's pooled readout was wrong before it was right, and the control is what showed
it.** Reported over every significant alignment instead of each ORF's best, the divergent
pile read 35.8% median identity with 73% of alignments saturated and the 95–99% control read
36.9% with 71% — indistinguishable. A pile's alignments are dominated by its long tail of
distant subjects rather than by its own nearest relative, so the pooled number describes the
subject set. Both blocks are in [`TABLES.md`](TABLES.md) §7 with the pooled one marked as
the one not to read.

## Structure: it finds the right thing, it finds nothing new, and it cannot tell the piles apart

This is the arm Imo asked for by name, and it took the longest to get running for reasons
that were entirely mechanical — they are in the faults section. What it measures is worth
stating carefully, because the obvious reading of it is wrong in two directions at once.

**The target changed, and the measurement is why.** The plan was AlphaFold/Swiss-Prot,
chosen because it is Swiss-Prot folded: Foldseek against it would be arms 1–2's own subjects
in a different alphabet, which is the only way to compare instruments rather than databases
in a project whose standing confound is that the subject set dominates reach. That
reasoning had a hole, and [`TABLES.md`](TABLES.md) §5b counts it: the shared sequence
database is Swiss-Prot **plus** the RefSeq viral proteins, and **47 of the 51 ORFs with a
significant baseline hit find it in the RefSeq viral half, which is folded nowhere**. A null
from that target would have been a statement about Swiss-Prot's coverage of ssRNA phage
proteins, and it would have read as a statement about structure.

So the arm runs against **every experimentally determined PDB entry whose source lineage
includes Viruses** — 27,439 entries, 178,660 chains — selected by a query to the RCSB search
API with its date recorded, not by a list anybody typed. That universe is chosen by what has
been crystallised, which is independent of anything the sequence arms found, and it is where
a remote structural relative of a phage replicase or capsid actually sits.
`experiments/reach_structure.py` runs it in the container for nothing: ProstT5 writes the
3Di alphabet for all 202 ORFs in 652 s, the subject database builds in 442 s, and the
exhaustive search and readout bring the whole arm to 1,249 s on 4 vCPU. It calls
`reach_arms.call_orfs`, and its manifest carries the sha256 of the `orfs.faa` it produced,
which is **byte-identical to the paid run's** — so this is the same 202 genes, checkable
rather than asserted.

**What it found.** Every real contig gets a significant structural alignment, and the
assignments are correct biology arrived at independently:

[`TABLES.md`](TABLES.md) §5c counts it per gene, because §1 counts contigs and a contig
counts as placed if any one of its genes is — which for this arm hides the whole shape of the
result. Every structural placement is listed there with the entry it landed on.

| pile | contigs placed (§1) | genes placed, of (§5c) | genes structure placed that sequence did not | median 3Di identity of the best placement |
|---|---|---|---|---|
| divergent | 7 of 7 | 8 of 24 | **0** | 32.6% |
| near-reference | 7 of 7 | 8 of 23 | **0** | 31.2% |
| holdout | 1 of 1 | 1 of 3 | **0** | 31.5% |
| already-referenced companion | 1 of 1 | 4 of 4 | **0** | 47.1% |
| shuffled, 5 seeds | 0 of 80 | 0 of 148 | **0** | — |

Seven of the eight divergent placements are the **polymerase on a Qβ replicase structure**
(3MMP, 3AVX, 3AVT, 4R71 — the RNA-dependent RNA polymerase of the type species of this very
family), and the eighth is a **coat protein on a bacteriophage virus-like particle** (6YFM).
The holdout's polymerase lands on Qβ replicase at E = 1.1e-22, the strongest statistic any
test contig gets from this arm. The already-referenced companion is the one that breaks the
pattern, placing all four of its genes and topping out on a tobamovirus helicase at 81.7% —
a different molecule, which is what it is.

**Three readings, and only the third is the result.**

It is tempting to read "0 structure-only placements" as structure failing. It is not that:
the sequence arms place 22 of 24 divergent genes and structure places 8, against a subject
set two orders of magnitude smaller in which most of these genes have no solved relative at
all. One gene per contig is what it looks like when the PDB holds this family's replicase
fold and little else of it. Coverage is not comparable between the two arms and nothing here
makes it so.

It is equally tempting to read the correct assignments as structure working. That is the
weaker half of it. **The number that decides is the control, and the control moves with the
test**: 32.6% median 3Di identity on the divergent pile against 31.2% on the near-reference
pile, whose contigs are 95–99% identical to published sequence. The same subjects, the same
E-value range, the same identities, on two piles 48 sequence-identity points apart. The one
place it does separate is the already-referenced companion at 47.1% — which is the
familiarity polarity this project keeps finding, not a novelty signal.

So the result is the third reading: **against the structural universe that exists, predicted
3Di comparison recovers the right fold for the conserved gene and carries no information
about how diverged that gene is.** The shuffled null at 0 of 80 says the arm is not reading
noise; the control at 31.2% says what it does read is not distance. That is a field reported
with the fact that it reads the same on both piles, which is this repo's rule, and it means
the arm belongs in the identification half of the pipeline rather than the novelty half.

**What a stronger version would cost, priced rather than waved at.** Two numbers neither of
which came from a guess. No TM-score or LDDT is reportable from this arm at all: both are
superposition scores needing C-alpha on both sides, and a ProstT5 query has no coordinates —
foldseek says so itself. Folding the 54 real ORFs with ESMFold is about $4 on a GPU instance
and buys those two numbers. Making the comparison genuinely controlled is the expensive one:
it needs the *subject* side encoded, and the RefSeq viral half alone is **175,064,339
residues**, which at this arm's own measured **66.3 residues/s on 8 vCPU** is thirty days of
CPU. That is a GPU job whose rate we have not measured, so it is named as the next step and
not costed to a figure this report would have to stand behind.

## The one null leak, and what it was

Every arm but one put nothing above significance on any of the 400 shuffled records — the
structure arm included, at 0 of 80 real-contig-matched shuffles per seed. One arm did not,
and it is worth the space because the fix generalises.

Floored on the full-sequence E-value, jackhmmer returned seven significant alignments on
1 of 80 null records at one seed. All seven were the same 23 residues of one shuffled ORF
against one Swiss-Prot protein, at subject positions 318, 346, 374, 402, 430, 458 and 486 —
**a 28-residue tandem repeat, matched once per copy.** Every individual domain sat at an
E-value between 2,000 and 2,900, which is nothing. The full-sequence value reached 1.1e-08
by summing seven of them.

A full-sequence E-value answers *does this query have a relative*; a table of alignments
needs the value that answers *is this alignment real*, which is the domain's own. Flooring
on that removes the leak and costs the real piles nothing — 6 of 107 hmmsearch rows and 2.7%
of jackhmmer's go, every query on every real pile keeps its hit, and the null is 0 of 80 at
every seed for every arm. **So this is a repeat artefact named, not a sensitivity traded**,
and it is the same aggregate-against-per-unit shape as the three failures of 2026-11-20.

It also puts a number on what iterative profile search actually bought. jackhmmer returns
45,860 significant query-subject pairs against the baseline's 11,611 — **3.9 times as many,
and a strict superset: it misses none of the baseline's.** Of the 34,595 pairs it adds,
12,804 are on the divergent contigs and 13,350 on the 95–99% control. Equal again, so they
are the marginal tail, like the exhaustive arm's. Five iterations of profile building
quadruple the hit list, find no nearer relative for any of the sixteen contigs, and are the
only arm on the panel that puts anything above significance on a shuffle.

## What the runs bought, which was ten faults

They are reported because it is where most of this experiment's correctness came from, and
because the split between the two kinds matters. **A tool that fails exits non-zero; a
readout that fails returns a number — and the worst kind exits zero and returns nothing.**
Five of the ten were in readouts and every one produced a plausible table. Five were in
tools or in the services behind them; three of those announced themselves, and two did not.

The two that did not are the same fault twice, and it cost more than any other single thing
here: **`foldseek databases` exits 0 when its download fails.** Every prebuilt Foldseek
database is served by one Cloudflare worker. That worker answered HTTP 429 with a 17-byte
body from three separate networks across several hours on 2026-11-22; foldseek untarred the
error body, printed `gzip: stdin: not in gzip format`, moved nothing, and returned success.
It bit once on the ProstT5 weights, where a file test caught it, and then again on the
subject databases, where nothing was testing — so a run added both structure arms on that
exit code and spent 500 s per arm failing on `Input .../afdb_sprot does not exist`. **A
whole arm was reported absent for a reason that was not the real one**, which is worse than
reporting it absent, because the recorded reason was a search failure rather than a staging
failure. Every structure database is now decided on the files and the arm refuses before
searching if its subjects are missing.

The follow-on is worth recording too, because it is the one failure here that no amount of
checking would have fixed. The only bulk-servable alternative to that worker, EBI's 28.6 GB
`swissprot_pdb_v6.tar`, came down at 6.7 MB/s and then **stopped advancing entirely** at
33 GB — no error, no exit, and no socket timeout on the curl, so it consumed a whole run's
dead-man timer at $0.34/h before anything ran. The lesson is the cheap one: a long download
in staging needs `--speed-limit`, and a target that takes an hour to stage does not belong on
a paid instance when the alternative runs in the container for nothing.

### The one that mattered: a control had replaced the test, and the run looked clean

The profile-to-profile arm splits jackhmmer's concatenated Stockholm into one alignment per
query and names each file after that query. It found the name by substring — the first ORF
name appearing anywhere in the block. **Every shuffled record on this panel is named
`sh<seed>_<qid>`, so a real ORF's name is a substring of its own null's name.** Block
`sh5_ctl1__orf1` matched `ctl1__orf1`, was written to that filename, and overwrote the real
query's alignment with its own control. 202 ORFs produced 59 files, nearly all of them
holding a shuffle.

Then the second half: a shuffled ORF has no homologs, so its alignment holds one sequence;
`hhmake` refuses a profile built from one sequence; so the arm finished with **zero rows and
`ok: true`**, which reads as "profile search found no relative" and is a different claim
entirely. Nothing in the run said otherwise.

The name is now parsed from jackhmmer's own `#=GF ID <query>-i<N>` line and matched exactly,
with unknown names, absent names and collisions all raising, and an arm that builds no
profile at all now fails rather than returning an empty table. Validated against that run's
real 182 MB alignment: 202 of 202 blocks named, no collisions, real alignments at median
depth 951 sequences against the nulls' 1. `tests/test_reach_hhsearch_split.py` pins it.

### The codon readout pooled its piles, and the control moved with the test

Described in full above. The fix is one line — take each ORF's best-scoring alignment rather
than every alignment — and it moves the divergent pile's saturation from "indistinguishable
from the 95–99% control" to "half its genes, against none of the control's".

### Three in the tools and the services behind them

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

**hhmake does not read Stockholm, and said so as "alignment too shallow".** It reads A3M,
A2M and aligned FASTA; handed jackhmmer's output it exits non-zero, and the arm's one
`except` branch labelled every failure a shallow alignment — which is the correct label for
a depth-1 null and the wrong one for 51 real alignments at a median depth of 951. The
zero-profiles guard added after the collision is what caught it, by failing the arm instead
of reporting an empty table. The blocks are now converted here, concatenating each
sequence's wrapped rows rather than taking its first occurrence, and a refusal of an
alignment with two or more sequences raises with hhmake's own words.

**Foldseek's segfault, which four runs read as an exit 1.** `--prostt5-model` takes a
directory and foldseek looks inside it for a weights file *by name*, and the accepted names
changed between releases. Pinned release 9-427df8a accepts `model.safetensors` or
`model.gguf`, and the only gguf its own downloader ships is a q4_0 quantisation; handed the
f16 gguf this account holds — the right model, the right format, 2.4 GB, in the directory it
was pointed at — it crashes inside the model load. Release 10-941cd33 looks for
`prostt5-f16.gguf`, which is exactly what that file is, and reads it.

Three things made this expensive and all three are worth keeping in mind. The arm used
`check=True` without capturing output, so foldseek's own message was discarded for three
runs. `foldseek databases ProstT5` is not a usable fallback: its weights URL answers
Cloudflare 1027, daily quota exhausted, from two separate networks, and foldseek untars the
17-byte error body, fails, and **exits 0** — which is why run 3's staging reported success
with no weights on disk. And the whole thing was diagnosable in the container for nothing:
the binary is a 50 MB download and the weights come out of S3 in twenty seconds, so the
segfault reproduced in two seconds locally once anyone tried. Running the arm's full command
sequence locally then caught a further fault before it cost anything — `convertalis` cannot
report `alntmscore` or `lddt` from a ProstT5 database, because both need C-alpha coordinates
and a predicted 3Di string has none.

The codon walk was tested before it ever ran, in `tests/test_reach_codon.py`, for the same
reason: it reconstructs a codon alignment from a gapped protein alignment, a CDS in
nucleotides and subject coordinates that can run backwards, and an off-by-one in any of them
gives dN/dS values that are wrong and entirely plausible. The four cases are identical
sequences, third-position changes only (which a frame off by one turns non-synonymous), the
same alignment with the genome flipped, and a saturated pair that must come back undefined
rather than clipped.

## Blind spots

- **The control is matched on family and length but not on depth.** The test contigs sit at
  116-467x mapped depth and the near-reference controls at 27-112x, with no overlap. Depth
  should not reach a protein search on an assembled contig, but it does reach the contig's
  own accuracy, so a residual difference between the piles could in principle be assembly
  quality rather than divergence. The three arms that would be most sensitive to that - the
  ones scoring per-residue identity - put the two piles 48 identity points apart, which is
  far more than any plausible assembly-error term, so this is recorded rather than
  controlled.
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
- **The profile-to-profile arm has no shuffled control, and cannot have one.** Its query
  profile comes from jackhmmer's alignment, and all 148 shuffled ORFs produce an alignment
  of depth 1, which `hhmake` refuses. Profile building needs homologs and shuffled sequence
  has none, so that arm's only control is the near-reference pile. This is a property of
  the method rather than of this panel, and it applies to any profile-based novelty claim:
  the usual null is unavailable by construction.
- **The query profile is shallower than the method's own recipe.** The standard build
  iterates against UniRef30; this one iterates against Swiss-Prot plus RefSeq viral, to keep
  the query side's subject universe identical to every other arm's. So this arm's reach is
  a lower bound on HHsearch rather than a measurement of it, and the trade was made for
  comparability on purpose.
- **The structure arm's universe is not comparable with the sequence arms', and no table
  here pretends otherwise.** 178,660 solved chains against 1,210,551 sequences, selected by
  taxonomy rather than by sampling fold space. Its coverage number is therefore not a
  sensitivity and its "0 structure-only placements" is not a kill: both are properties of
  what has been crystallised. The claim it supports is the calibration one — same identity
  on both piles, null at zero.
- **A predicted 3Di string is not a structure.** ProstT5 was trained on sequences with known
  folds, so on a protein with no structural relative in its training set the alphabet it
  writes is an informed guess, and its errors point towards folds the model has seen. That
  is the same direction as the bias this arm would need to be free of, so a negative from it
  is weaker than a negative from folded queries, and the $4 that fixes it is named above.
- **3Di identity and sequence identity are different scales and are never compared here.**
  The near-reference pile reads 31.2% structural identity at 95-99% sequence identity. Any
  reading across the two columns is meaningless; only the within-column comparison between
  piles is used.
- **The exhaustive arm is exhaustive in alignment, not in scoring.** BLOSUM62 at 11/1 with
  published lambda and K: a different matrix is a different instrument and none was tried.
- **Timings are from one instance type on one day.** Fixed and marginal cost are not
  separated here as `results/bigindex/` separates them, because every arm's query set was the
  same size.
