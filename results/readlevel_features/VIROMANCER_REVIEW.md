# Viromancer / Lyra-Virus: full-paper review

**Date:** 2026-11-20
**Paper:** *Generalizable and Agent-Contextualized Metagenomic Viral Detection with Viromancer*.
NeurIPS 2026 workshop oral, "Agentic AI for Biological Discovery".
Public listing (neurips.cc/virtual/2026/172554) names Feldman, Arze & Sabeti; the PDF supplied
to this thread is the anonymised submission ("Anonymous Author(s)"), 39 pages including a full
Methods appendix.
**Status of the read:** the complete paper, methods and all, was read. This supersedes the
abstract-only entry in `PRIOR_ART.md`, which flagged calibration and pair aggregation as unknown.
Both are now answered.

This file exists because `PRIOR_ART.md` named this paper the one live threat to the per-read
novelty direction and said to read it before claiming the ground was open. The verdict is that it
takes two things from us and leaves the core claim standing, with a new obligation attached.

---

## 1. What the system actually is

**Lyra-Virus** is a binary, read-level, **nucleotide** classifier: raw read in, one scalar
probability that the read is viral, no taxonomy. Architecture is the Lyra structured state-space
model — six-symbol token embedding (A/C/G/T/N/pad, d = 384), two Projected Gated Convolution
layers, three stacked S4D diagonal state-space layers at O(L log L) via FFT, learned positional
pooling, three-layer MLP to a sigmoid. About **2.1 M parameters**, trains to convergence in
**under two hours on two RTX 4090s**.

Separate checkpoints are trained per read length (75, 100, 150, 250, 300 bp); the pipeline infers
the dominant read length from the first 1,000 alignments and picks the nearest checkpoint. Reads
under 50 bp are excluded.

Training data: ~430,000 curated sequences, 199,000 viral and 233,000 non-viral. Viral from NCBI
Virus, 116 families across six ICTV realms, per-family download cap of 20,000 (5,000 for
pandemic-inflated families such as *Coronaviridae*), dereplicated at 99% with MMseqs2 linclust.
Non-viral from 12 categories chosen to mimic *clinical* background: gut commensals, environmental
bacteria, respiratory pathogens and flora, oral and skin flora, archaea, reagent "kitome"
contaminants, clinical fungi, human genomic fragments, eukaryotic model organisms, plasmids.
No-template-control reads and ERCC spike-in constructs are injected at training time as hard
negatives.

Holdout: entire viral families assigned to five bins by realm-stratified greedy bin-packing
(their Algorithm 1), four folds evaluated, ~60/20/20 by family. Genuinely conservative.

**Viromancer** wraps it: a read branch (Lyra-Virus + Kraken2 + Centrifuger, then BLASTn on
Lyra-Virus calls) and a contig branch (MEGAHIT → geNomad + Foldseek), merged by a concordance
layer, with four Claude Sonnet 4.6 agents (read, contig, adversarial skeptic, orchestrator)
adjudicating disagreements.

## 2. The numbers that matter to us

**Family-level generalisation.** 150 bp model, mean test **AUROC 0.88** across four folds; above
0.89 in folds 0–2, falling to 0.82 in fold 3 (which withheld *Flaviviridae* and *Picornaviridae*).
AUROC rises monotonically with read length from 75 to 300 bp.

**Operating point.** Deployment threshold t = 0.8. The F1-maximising threshold was t ≈ 0.075–0.10,
which they reject because "even a low false-positive rate can produce an unmanageable number of
candidates when applied to tens or hundreds of millions of reads". The checkpoint-selection
criterion during training is **TPR@FPR10 — sensitivity at a 10% false-positive budget**.

**Versus a k-mer baseline** (logistic regression on k-mer frequencies, k ∈ {3,4,5,6}, same
family-level holdout): AUROC 0.929 vs 0.81; at t = 0.8, precision 96.7% vs 83.5%; and the one they
lead on, **at a fixed 1% false-positive rate, Lyra-Virus recovers 34.8% of viral reads against the
baseline's 7.5%**. They note themselves that the baseline is linear on raw frequencies and so is a
floor, not a ceiling, for what composition can do.

**How much of it is composition.** Trinucleotide frequencies reconstruct about **29% of the
variance** in Lyra-Virus scores across held-out families (declining with k, negative by k = 6).
The S4D layers add ~0.18 AUROC over the convolutional embedding (0.579 → 0.760 across the three
layers). Learned features include reverse-complement-consistent tetramer coefficients, suppression
of homopolymers (AAAA, TTTT) and of CpG-containing tetramers.

**The Evo2 benchmark** — the closest thing in the paper to a novelty test. Evo2 generates
continuations from short primers drawn from known genomes: a gut phage corpus (1,493 contigs,
73.8% with no detectable homology to source) and a *Microviridae* corpus from a fine-tuned model
(3,716 contigs, 80.6% high identity), plus 1,074 bacterial/human control contigs. Reads are
stratified by MMseqs2 homology to the source corpus into primer / detectable / novel (<30% identity
or no alignment).

| | Microviridae reads | gut phage reads | control FPR |
|---|---|---|---|
| Lyra-Virus, t = 0.5, all reads | 54.3% | 21.9% | 9.9% |
| Lyra-Virus, t = 0.5, **novel regions only** | 48.3% | 21.6% | — |
| Lyra-Virus, t = 0.8 | 33.0% | — | 2.9% |
| Centrifuger, novel regions | 0.2% | — | (99.8% specificity) |
| Kraken2, novel regions | ~0 | ~0 | — |

Their Figure 3E gives the false-positive rate on the negative control corpus falling monotonically
from 37.6% at t = 0.1 to 1.4% at t = 0.9.

**The cleanest background number in the paper.** Five no-template controls from the same
respiratory cohort, 24.7 M reads, no biological input. At t = 0.8 Lyra-Virus flagged 2,363 reads,
**0.010% of control reads — about 100 per million** — ranging 0.00008% to 0.044% (roughly 0.8/M to
440/M) across individual controls. Kraken2 flagged 41 (1.7/M), Centrifuger none. No flagged read
had multi-classifier support and no viral contig assembled from any control.

**Prophage as a false-positive source.** 11.5% of reads from annotated prophage regions exceed
t = 0.5 against 5.2% of reads from the surrounding chromosome, with the prophages left unremoved
and labelled non-viral in training. Free lytic phage scores much higher (median 0.81 for T4, 0.60
for T7).

**Clinical deployment.** 297 nasopharyngeal respiratory metagenomes (283 with a reported viral
infection, 14 bacterial-only), Illumina NovaSeq 150 bp paired-end. Plausible viral signal in
261/283; the panel-reported pathogen recovered in 219/283 (77.4%); at least one reported pathogen
in 228/283 (80.6%). Per-pathogen recovery ranges from influenza 79/80 (98.8%) down to adenovirus
1/10 (10.0%). Ten stool spike-in samples with ERCC, MS2, Twist SARS-CoV-2 and Twist measles at
defined copy numbers.

**Grouping and triage.** The read agent does not see individual reads. Reads are clustered by top
BLAST accession, or — when they have no BLAST hit but Lyra-Virus and Kraken2/Centrifuger agree —
by species-rank taxid. Each group is represented by up to five reads and carries its **reads per
million (RPM)**. A stratified RPM filter gates agent review: **≥0.5 RPM for a curated list of
high-consequence human pathogens, ≥20 RPM for commensal, environmental and bacteriophage groups.**

**Not in the paper at all:** any mate-pair or read-pair aggregation (R1/R2 are extracted from the
host-depleted BAM and scored independently); any environmental or wastewater data; any protein-space
or absence-based scoring; any code, weights or data availability statement. One cited reference by
an author is a position paper on restricted release of advanced biological models, so released
weights should not be assumed.

---

## 3. What this takes from us

**(a) Per-read scoring under family-level holdout is theirs.** They claim it explicitly — "to our
knowledge, Viromancer is the first metagenomic viral discovery pipeline to demonstrate this
generalization directly at the level of individual short reads" — and the evidence supports it.
We cite this as the per-read precedent and stop implying the read is an unoccupied unit.

**(b) "Nobody combines a per-read score with read grouping into one calibrated call" must be
narrowed.** They do combine them: per-read score → grouping → an RPM threshold → a call. What they
do not have is *calibration*. The RPM thresholds are hand-set and asymmetric by a curated pathogen
list, chosen for clinical triage cost rather than fitted to any error rate, and there is no curve
relating group-level sensitivity to group-level false-positive rate. The surviving form of our
claim is narrower and should be written that way.

There is also a structural hole in their grouping that is directly our target. A group is formed
by a shared BLAST accession, or by a shared species taxid where Lyra-Virus and a reference
classifier agree. **A read from something genuinely novel has no BLAST accession and no concordant
taxid, so it is never grouped and never reaches an agent.** Their own numbers say this is the
common case: at the deployment threshold only a minority of divergent reads survive, and the whole
BLAST-adjudication layer presumes a hit exists.

## 4. What survives, and why

**(a) The polarity is still ours.** Their question is viral versus non-viral throughout. The
family holdout tests whether *viralness* generalises to unseen lineages, not whether a read is new.
The Evo2 benchmark is the nearest approach and they are careful about it themselves: it "measures
detection of viral-conditioned sequence across decreasing homology, not the biological sensitivity
of each tool against a collection of known functional viruses". A read from a well-known virus and
a read from a novel one both score high by design. That is the axis our own work already found
uninformative in a virome — geNomad calls 98.3% of Berlin contigs over 3 kb viral — and their
background is human-dominated clinical respiratory, where viral-versus-not is a real contrast.

**(b) The operating point is still ours, and the gap is larger than we thought.** Every number in
the paper lives at 1–10% false positives, i.e. 10⁴–10⁵ per million. Their training criterion is
sensitivity at FPR 10%. Their headline comparison is at FPR 1%. Their best control-corpus figure is
1.4%. The one genuinely low number, ~100 per million on no-template controls, is measured on
*empty* libraries rather than complex background, and it is a background flag rate with no
sensitivity attached at the same point on the same data. So the claim stands and sharpens:
**nobody reports sensitivity at a fixed, low false-positive rate per million on real complex
background.** Sensitivity and false positives continue to come from different experiments.

**(c) The representation is untouched, and they state the objection to ours in one sentence.**
Lyra-Virus is nucleotide end to end. There is no protein space anywhere in the read branch, and
Foldseek is confined to contigs with this justification: *"This analysis is restricted to contigs
because a 150-bp read can encode at most approximately 50 amino acids, too little information for
reliable protein-level comparison."* That is about homology and structure comparison, which needs
an alignment good enough to trust. Absence does not: asking whether a peptide 10-mer appears in any
known proteome needs an index lookup, not a reliable alignment. The distinction is the whole
proposition, and a serious lab has now put the general objection in print without testing the
absence form of it.

**(d) Paired reads are untouched.** They have paired-end data and never use the pairing. The
sibling thread's axis is unaffected.

**(e) Environmental background is untouched.** Clinical respiratory and stool only.

## 5. The new obligation, and a warning

**Lyra-Virus is now the baseline any per-read score has to beat.** At 2.1 M parameters and under
two hours on two consumer GPUs, it is cheap for us to be asked for and expensive for us to refuse.
Their own k-mer baseline is linear on raw frequencies and they concede it is a floor, so "we beat a
k-mer baseline" will not substitute.

But there is a reason to expect its headline not to transfer to our setting, and it is worth saying
plainly because it is testable. **A substantial part of Lyra-Virus's signal is compositional** —
trinucleotide frequencies reconstruct ~29% of its score variance, and the features it learned are
tetramer, CpG and homopolymer patterns. Our own feature map killed composition as a *novelty*
feature on polarity grounds: in a wastewater virome the background is coding viral sequence too, so
the viral-versus-bacterial compositional contrast that Lyra-Virus exploits largely disappears.
Their prophage result is the same effect seen from the other side — viral sequence sitting inside a
bacterial genome scores viral and becomes a false positive.

The honest prediction is therefore that Lyra-Virus transfers poorly to an environmental virome, and
that its 0.88 AUROC is partly an artefact of a background chosen to be non-viral. **We should not
assert that; we should measure it**, and we cannot yet, because no weights are released.

## 6. Consequences for this thread

1. Rewrite the Viromancer entry in `PRIOR_ART.md` against this file. Done in the same commit.
2. Restate the surviving claim as: *sensitivity at a fixed low per-million false-positive rate on
   real complex background, for a NEW-versus-KNOWN score, with a calibrated grouping layer.* Three
   qualifiers, each of which this paper leaves open.
3. Any build proposal must name Lyra-Virus as the baseline and say what we do when the weights are
   not available — most likely retrain the published architecture, which is fully specified in
   their appendix and cheap, rather than claim it as unavailable.
4. The compositional-transfer prediction in §5 is a cheap, falsifiable experiment on data already
   on disk once a model exists, and it is the kind of continuous readout the project's process
   rules ask for. It belongs early.

## 7. Caveats on this review

The PDF was supplied by Imo; it is the anonymised submission and carries no DOI. Figure panels were
read as rendered text and numbers quoted from figure captions (notably the 37.6% → 1.4% false-
positive curve) come from caption prose rather than from underlying data. Percentages converted to
per-million in §2 and §4 are my arithmetic on their percentages, not figures the paper states. The
workshop is non-archival and the paper is under review, so numbers may move.
