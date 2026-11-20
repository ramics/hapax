# Prior art for the two live read-level directions

Reviewed 2026-11-20, on Imo's request, before running the `rho_fp` measurement that
`REPORT.md` §3.1 recommends. Two strands: (1) peptide k-mer absence as a per-read novelty
score, (2) paired-read evidence for novelty and false-positive suppression.

Three findings change what `REPORT.md` says. They are applied there as well as recorded
here, and listed first because two of them cost us claims.

---

## 1. Corrections to our own claims

### 1.1 The per-million claim must be narrowed. PACIFIC already uses the unit.

Standing claim in project memory: *nobody publishes sensitivity at a fixed false-positive
rate per million background reads; published rates are in percent, four orders of magnitude
from what one-in-a-million needs.*

PACIFIC (Acera Mateos et al., *Sci Rep* 11:3368, 2021,
[PMC7864945](https://pmc.ncbi.nlm.nih.gov/articles/PMC7864945/)) reports false positives in
exactly our unit. Table 2 is captioned "Average false positive reads across 100 experiments
for each viral class", with columns headed **`FP/1 M reads (Balanced)`**,
**`FP/1 M reads (Unbalanced)`** and `FP% threshold`:

| Class | FP/1M (balanced) | FP/1M (unbalanced) | FP% threshold |
|---|---|---|---|
| Coronaviridae | 81.6 | 66.6 | 0.0405 |
| Influenza | 1.9 | 6 | 0.000807 |
| Metapneumovirus | 2 | 1 | 0.000154 |
| Rhinovirus | 35.5 | 53.9 | 0.0418 |
| SARS-CoV-2 | 2 | 1 | 0.000213 |

So the unit is not ours and never was. What survives, checked against the full text:

- **No sensitivity-at-a-fixed-FP operating point.** Sensitivity (Table 1, recall
  0.997–0.999) and false positives (Table 2) come from different experiments, at one fixed
  per-read rule, with no curve and no FP column beside the sensitivity. You cannot read
  sensitivity off at a chosen per-million rate.
- **All simulated.** Both FP tables use ART-simulated reads. The 63 real SRA runs are
  reported as read counts with no ground truth, so no real-background FP rate is computed.
- **Thresholds are sample-level, not calibrated per read**: the 0.99 quantile of per-sample
  FP% from a moment-matched Beta fit, with a separate posterior ≥0.95 per-read rule.
- **No mate-pair aggregation.** PACIFIC's "double prediction" is a read and its own reverse
  complement, not the two mates of a pair. Benchmarks are single-end.

**Revised safe claim:** *no published read-level detector reports sensitivity at a fixed
false-positive rate per million background reads, measured on real background.* Do not
write "nobody expresses false positives per million" — a 2021 paper has it in a table.

### 1.2 Merged mates (§3.3) is dead, and not for the reason I flagged.

I flagged §3.3 as resting on an unverified insert size. The mechanism is worse off than
that: it has no headroom in this substrate, and the one group that analysed it concluded
the opposite of our hypothesis.

- The NAO's own read-length analysis argues on cost-per-informative-base grounds that **a
  positive gap beats any overlap, and more overlap is worse** — a 2×150 run on an 80 bp
  fragment yields only 80 bp of unique sequence, and short inserts read through into
  adapters.
- Wastewater RNA at 2×150 is **already mostly overlapping**: NAO reports viral read pairs
  not collapsed at 42% for Yang et al. 2020 and 24% for Spurbeck et al. 2023, i.e. roughly
  58% and 76% already merge. There is nothing left to engineer by shortening inserts.
- Merging tools (BBMerge *PLoS ONE* 2017, FLASH 2011, PEAR 2014) are framed throughout as
  preprocessing for assembly and amplicon work. No paper argues insert length as a
  detection-sensitivity lever. Published levers are enrichment and depth.

**Action: drop mechanism C.** Unoccupied because it does not work, which is the least
useful kind of unoccupied.

### 1.3 The mate conjunction has been tried, informally, and rejected.

Our §3.1 framing — that nobody frames the conjunction this way — is too strong.

The Nucleic Acid Observatory's open lab notebook (Bradshaw, entries
[2024-02-08](https://naobservatory.github.io/wills-public-notebook/notebooks/2024-02-08_crits-christoph-2.html)
and
[2024-02-15](https://naobservatory.github.io/wills-public-notebook/notebooks/2024-02-15_crits-christoph-3.html),
on Crits-Christoph wastewater data) implements both rules — both mates above a threshold,
or either mate above it — across thresholds 15–45 on length-adjusted Bowtie2 scores, and
**chooses the disjunction**: "pass a pair if either the forward or reverse read had an
adjusted alignment score of at least 20". Their whole FP-reduction trajectory
(992 → 549 → ~100 → 40) is computed under the either-mate rule. They plot forward against
reverse scores and never compute a correlation; no per-rule sensitivity/FP numbers are
printed. Unpublished, and the most likely thing to be raised as prior art.

This is better for us than it looks, and it sharpens the experiment rather than killing it.
Our model says the conjunction's value swings four orders of magnitude on `rho_fp`, so a
group finding the disjunction better on their data is **a data point about their `rho_fp`,
not a refutation**. The defensible claim is now narrower and more interesting:

> The mate-score correlation has never been measured. taxMaps shows mate concordance tracks
> precision at ρ = 0.994 without explaining why. Our model predicts both why NAO's
> disjunction won on their data and the conditions under which a conjunction wins instead.

Also relevant, and now the thing to calibrate against: **taxMaps** (Corvelo et al.,
*Genome Res* 28:751, 2018) classifies each mate independently and uses mate concordance as
a precision proxy across seven real metagenomes and five classifiers, validating it at
Spearman ρ = 0.994 against true precision on simulated data. It never uses concordance as a
filter and never measures mate-vs-mate correlation.

---

## 2. Strand 1: peptide k-mer absence — not occupied, and the parts exist

No paper scores individual short reads for novelty by six-frame translation plus counting
peptide k-mers absent from a protein reference. Every near neighbour differs on polarity
(presence, not absence), unit (proteomes, genes, contigs, samples — not reads), or purpose
(classification or annotation, not novelty).

**Use rather than rebuild:**

- **Rarity index** (*NAR Genom Bioinform* 2024,
  [PMC10993293](https://pmc.ncbi.nlm.nih.gov/articles/PMC10993293/)) — a graded per-k-mer
  absence score over 45,785 genomes and 21,871 proteomes, including per-superkingdom
  indices for viruses. **This is the scoring primitive, already published**, applied to
  k-mers as objects of study and never to reads as queries. Prior art for the score, not
  for the detector; take the weighting from here.
- **MetaProFi** (*Bioinformatics* 2023) — chunked Bloom filter over amino-acid k-mers,
  queryable with nucleotide input via six-frame translation, **zero false-negative rate**.
  The index this method needs, already built, with exactly the property an absence score
  requires.
- **k calibration** (Pierce-Ward & Brown 2022, sourmash group) — at k_aa = 7 over 80% of
  peptide k-mers are shared across phyla; k_aa ≈ 10 approaches nucleotide-level
  specificity. So k = 7 is too promiscuous to support an absence signal; start at 10.

**Adjacent, correct polarity, wrong unit:** nullomers/nullpeptides (*Genome Biol* 2021,
absence is under negative selection, so it carries signal); WhatsGNU (*Genome Biol* 2020,
per-protein novelty from exact matches, GNU score 0 = first known allele); DE-kupl
(*Genome Biol* 2017) and Kevlar (2019), which retain reads by nucleotide k-mer absence from
a reference or control — the closest precedent for absence-filtered reads as an object.

**Adjacent infrastructure, opposite polarity:** Kaiju (*Nat Commun* 2016) already does
six-frame translation with maximum exact amino-acid matches — variable-length MEMs, not
fixed 12-mers — and nothing scores its unclassified fraction as novelty. sourmash supports
protein/Dayhoff/hp encodings and `gather --output-unassigned`, but reports unassigned
fraction **per sample** and states it does not classify individual reads. Metabuli
(*Nat Methods* 2024) is the closest per-read protein scoring in production, presence-based.
KrakenUniq's unique-k-mer count is used in the opposite direction, as a confidence filter.

**The cleaner gap than the one we claimed:** nobody quantifies unknown-ness at **read**
granularity. AGNOSTOS (*eLife* 2022) is the canonical framework for the unknown coding
space — 416M genes, known/unknown gene clusters — and requires assembly and gene
prediction, the step this method skips. Viral dark-matter reviews put 40–90% of viral genes
without known homologs, always per-gene or per-sample. (Note: AGNOSTOS's own unknown-fraction
percentages are inconsistent between its abstract, digest and gene counts; check Methods
before quoting one.)

**The rival framing looks weak at our length.** "Beyond Nativeness: Viral Proteins in
Protein Language Models" (ICML 2026) finds a dominant nativeness axis aligned with masked
perplexity, ordering cellular → viral — but **linear probes beat perplexity-only
classifiers**, and all of it is on full proteins. Nothing evaluates ~35–50 aa fragments,
where perplexity should degrade badly. This argues for peptide k-mer absence ahead of a
PLM, consistent with `REPORT.md`'s A12-before-A13 ordering.

**The failure mode to measure first.** Error-driven absence is characterised in nucleotide
k-mer work ("a sequencing error generates, with high probability, k-mers that occur only
once") and **never for protein k-mers** — where it is worse: one substitution corrupts a
whole codon and therefore every peptide k-mer overlapping it, roughly k consecutive k-mers
per error, and a frameshifting indel destroys the remainder of the frame. No paper
quantifies this amplification. Under-sequenced known organisms producing false novelty is
documented and **uncorrected** (Chorlton, *Front Bioinform* 2024, recommends broad database
inclusion and offers no statistical correction); the FracMinHash authors hit 45% unmatched
k-mers in a mock metagenome and explicitly declined to apportion it between error and
novelty. Low-complexity masking is standard, with no protein-k-mer-specific treatment.

---

## 3. Strand 2: paired-read mechanisms

| Mechanism | Status | Key prior art |
|---|---|---|
| **A** — mate conjunction as FP suppression | **Open in published literature**, tried informally and rejected by NAO; the correlation it depends on measured nowhere | NAO notebooks 2024-02; taxMaps 2018 (ρ=0.994); ROVER 2014 |
| **B** — one-sided anchor | **Method fully occupied since 2010**; novel only as a re-application, and its confound dominates | NovelSeq 2010; ViFi 2018; InserTag 2020; dengue/Nextera chimaera study |
| **C** — merged mates / short inserts | **Open but wrong** — drop | NAO read-length analysis; BBMerge/FLASH/PEAR |

**How classifiers actually combine mates — all summation, none conjunction.** Kraken 2
concatenates mates with an `N` and runs LCA once over the union of hits, so a strong R1 with
an empty R2 still classifies; the paper makes no claim that pairing reduces false positives.
Kaiju emits one line per pair and never states in Methods how mate scores combine, offering
only a **disjunction** rationale ("the chance of one mate overlapping with a protein-coding
gene is higher than for single-end"). Centrifuge sums hits over the combined pair length.
DIAMOND has no paired support at all, so protein-space pipelines built on it have no pair
logic to inherit.

**The pair as two quasi-independent observations: the middle ground is empty.** The variant-
calling literature treats mates as *fully dependent* and engineers the double-count away —
Strelka's documented false-positive inflation from double-counting overlaps, GATK/MuTect2
capping combined base quality and treating an overlap mismatch as one observation of one
molecule, fgbio `ClipBam` clipping template overlap precisely so callers do not double-count.
ROVER (*Source Code Biol Med* 2014) is the one published instance of our mechanism A in any
application — "only when both reads of a read-pair concur can the read-pair contribute" —
used as a design choice with **no quantification of the benefit and no correlation model**.
So the literature is binary: fully dependent, or silently fully independent. **The partial-
correlation middle our model occupies is genuinely empty, and that is the strongest claim
in this review.** Our ≪1x regime is also the low-depth corner where the independence
approximation is least bad.

**B is weaker than I wrote it.** "One-end-anchored" (OEA) pairs are a named primitive since
NovelSeq (*Bioinformatics* 2010) and are standard in insertion discovery (InserTag 2020,
1,696 non-reference insertions) and the SV-caller family. ViFi (*NAR* 2018) already does
anchor-plus-divergence-tolerant-scoring: it keeps a pair if one end maps to human and the
other maps to a viral genome **or scores well against profile HMMs for viral families** —
structurally almost exactly our B, so "it's only used for known viruses' integration sites"
undersells it. Applying it to free-living novel viruses in a metagenome is not found, but it
is a re-application, not a new mechanism. Binpairs (*PLoS ONE* 2014) reconciles mates that
both got assignments and pointedly does **not** transfer a classified mate's taxonomy to an
unassigned one — the exact move B needs, left undone.

And the confound is worse than I assumed. A dengue/mosquito Nextera study found **most, if
not all, detected host–virus chimaeras were artificial**, with control RNA from an unrelated
species showing about the same rate, mechanism identified as intra-molecular template
switching during PCR; the same signature is reported in Ebola and herpesvirus Nextera data.
No shotgun-specific chimaera rate is published; amplicon figures run 5–45%, and **>70% for
less abundant species** (Haas et al., *Genome Res* 2011) — our target is by construction the
rarest thing in the library — and chimaeras **form reproducibly across independent
amplifications**, so replicate agreement will not filter them. If we keep B, the chimaera
rate is the first experiment, not a caveat.

**What published methods actually use at ~1 read per million.** One read, remote protein
homology, a Bernoulli detection probability, and a hand-tuned per-study threshold. Ebinger
et al. (*CSBJ* 2021) invert the Bernoulli formula to predict the dataset size needed to see
**one** virus read at 99% probability, validated on 30 datasets — the direct intellectual
competitor to any pair-level detection statistic, with no pair structure. Grimm et al.
(medRxiv 2023) model sensitivity from relative abundance and depth, again with no pair layer.
vFams (*PLoS ONE* 2014) and SURPI make the case for translated remote homology at low
abundance, noting identity below ~30% defeats BLAST. Serratus calls novelty by a flat 90%
amino-acid identity cut on palmprints. **No published threshold in reads per million**, and
calling criteria are per-study.

**geNomad's score calibration is the machinery to extend**: it states plainly that raw
classifier scores are not probabilities because they ignore sample composition, and
recalibrates against the sample's empirical composition to give probabilities and FDRs.
Sequence-level, no pair layer. A reviewer will ask why we did not build on it.

**The live threat, now read in full — see `VIROMANCER_REVIEW.md`.** *Viromancer* /
*Lyra-Virus* (Feldman, Arze & Sabeti, NeurIPS 2026 **workshop** oral, "Agentic AI for
Biological Discovery"). The complete paper and its methods appendix have now been read, so
the "unknown calibration, unknown pair aggregation" caveat above this line is retired.

Verdict in short. Lyra-Virus is a 2.1 M-parameter **nucleotide** S4D state-space model that
scores single reads for **viral vs non-viral** — never new vs known — at AUROC 0.88 under
genuine family-level holdout, trained in under two hours on two consumer GPUs. It takes two
things from us: **per-read scoring under family holdout is now theirs to cite**, and our
"nobody combines a per-read score with grouping into one call" must narrow, because they do
combine them (group by BLAST accession or species taxid, gate by an RPM threshold) — what
they lack is *calibration*, those thresholds being hand-set and asymmetric by a curated
pathogen list.

What survives. Polarity: viral-vs-not is the axis our own geNomad reframe found
uninformative in a virome, and their background is clinical respiratory, where it is a real
contrast. Operating point: everything in the paper sits at 1–10% false positives (10⁴–10⁵
per million; their training criterion is literally TPR@FPR10), so §1.1's narrowed claim
stands and sharpens. Representation: nucleotide throughout, no protein space, no absence
scoring — and they confine Foldseek to contigs because *"a 150-bp read can encode at most
approximately 50 amino acids, too little information for reliable protein-level
comparison"*, which is an objection to homology comparison, not to an index lookup for
absence. Paired reads: they have them and never use the pairing. Environmental background:
absent.

The new obligation: Lyra-Virus is the baseline any per-read score must now beat, and it is
cheap enough that refusing to run it will not pass. No weights are released. One prediction
worth measuring rather than asserting: ~29% of its score variance is reconstructible from
trinucleotide frequencies, so much of its signal is compositional and may not transfer to a
virome background that is itself viral.

---

## 4. Net effect on the plan

Unchanged: **measure `rho_fp` first**. It is now better motivated, not worse — it explains a
published-adjacent group's choice rather than merely proposing a new knob, and the
partial-correlation treatment is the clearest empty space either strand found.

Changed:

1. **Drop mechanism C** (merged mates). No headroom in wastewater RNA; NAO's analysis argues
   against; §3.3 of `REPORT.md` amended.
2. **Demote mechanism B** (one-sided anchor) below A12. The algorithm is 15 years old, ViFi
   already does the divergence-tolerant version, and artificial chimaeras dominated real
   junctions in the one controlled study. Keeping it means measuring the chimaera rate first.
3. **Narrow the per-million claim** to sensitivity at a fixed FP rate *on real background*,
   per §1.1.
4. **Build A12 on existing parts**: rarity index for the score, MetaProFi for the index,
   k_aa ≈ 10 as the floor, and quantify peptide-k-mer error amplification before anything
   else, since it is unmeasured and structurally worse than the nucleotide case.
5. ~~**Read the Viromancer OpenReview paper** before claiming per-read novelty is
   unoccupied.~~ **Done, 2026-11-20** — full paper in `VIROMANCER_REVIEW.md`. Per-read
   *scoring* is occupied and must be cited; per-read *novelty* is not. Any build proposal
   now has to name Lyra-Virus as its baseline.

Nothing found makes the peptide-absence direction redundant. The honest statement of the
contribution is that the components are published and the composition is not — with the
correlation parameter that decides whether the composition is worth anything never having
been measured.

### Provenance and caveats

This review was conducted by literature search with primary-source verification on the two
claims that cost us something: PACIFIC's Table 2 and threshold rule were read in the PMC
full text, and the SecureBio engineered-virus post was read directly. Items not verified at
source, and flagged accordingly: AGNOSTOS's internally inconsistent unknown-fraction
percentages (the Viromancer paper has since been read in full, see `VIROMANCER_REVIEW.md`);
"kmdiff" could not be
confirmed as a real tool name. Several cited 2026 items are at or past the current date, so
treat their publication status as provisional. The SecureBio figure sometimes quoted as
"~0.4 alerts per billion read pairs" is arithmetic on 14 flags / 35 billion read pairs, is
not a figure that post states, counts pre-review flagged read collections rather than
reviewed alerts, excludes a second dataset (64 flags from 25 billion pairs), and its one
true positive is described there as a benign lab vector (a lentiviral vector, most likely
pHR'-CMVLacZ), not a threat.
