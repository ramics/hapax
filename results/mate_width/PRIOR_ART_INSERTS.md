# Can the wet lab deliver wide inserts from wastewater RNA?

Literature pass, 2026-11-20, done before committing to the wide-insert recommendation
in `REPORT.md` §5. The question: is SRR23973193's measured insert distribution (median
208 nt, p95 387, p99 ~515, 38% of pairs overlapping) a property of the **material** or
of the **protocol**?

**Answer: below ~1 kb it is the protocol.** Every total-RNA kit in wastewater use is
random-primed and tops out near 200–350 nt even on intact RNA. Above ~1 kb both the
material and Illumina clustering bite. Sources below, separated into measured numbers,
vendor specifications, and inference.

## 1. The kit that produced this run

- Wyler et al., *Comprehensive profiling of wastewater viromes by genomic sequencing*
  (bioRxiv 2022.12.16.520800; peer-reviewed as Environ Int 190:108875, 2024) state the
  library kit: **SMARTer Stranded Total RNA-Seq Kit v3 – Pico Input Mammalian** (Takara
  634485), 2×109 bp on a NovaSeq 6000.
  https://www.biorxiv.org/content/10.1101/2022.12.16.520800v1.full.pdf
- That kit's manual (Takara doc 111126): successful libraries show "a distinct curve
  spanning 200–1,000 bp, with a local maximum at ~300–400 bp", and "library preparation
  adds 150 bp to the size of the original RNA molecules" → **implied insert mode
  150–250 nt**. The manual also warns the kit "is good at retaining small inserts".
  *Vendor specification.* Our measured median of 208 nt is the kit's own stated output.
- **Not reported** by the authors: which fragmentation option was used (Option 1, 2–4
  min at 94 °C, versus Option 2, no fragmentation), and no RIN or DV200 for the
  extracts. So we cannot tell whether the RNA was deliberately sheared at all.

## 2. The decisive table: fragmentation time is not the lever

Illumina TruSeq Stranded Total RNA reference guide (doc 1000000040499), Table 12, with
inserts "determined after clustering, and sequencing with a paired-end sequencing run",
i.e. measured inserts rather than Bioanalyzer traces:

| 94 °C fragmentation | insert range (bp) | median insert (bp) |
|---|---|---|
| **none** | 130–350 | **200** |
| 1 min | 130–310 | 190 |
| 2 min | 130–290 | 185 |
| 4 min | 120–225 | 160 |
| 8 min | 120–210 | 155 |
| 12 min | 115–180 | 140 |

With the explanation: "the resulting cDNA fragments are smaller than the mRNA due to
internal priming by the random hexamers". *Vendor-measured.*

**Skipping fragmentation entirely on intact high-RIN RNA still gives median 200 bp
inserts.** Since every kit used in this field is random-primed (Illumina Stranded Total
RNA with Ribo-Zero Plus: fixed 94 °C/2 min, "expected insert size ~190 bp"; NEBNext
Ultra II Directional: times "calibrated for creating insert sizes ~200 nt", and >300 bp
inserts require both shorter fragmentation *and* a longer 42 °C RT), **400 nt+ inserts
are not reachable by tuning fragmentation in any standard total-RNA kit.** Reaching them
needs full-length RT (template-switching/oligo-dT) plus a DNA-side size selection.

Illumina also notes the discrepancy between Bioanalyzer and post-clustering insert size
is "due to the bias towards clustering smaller fragments", and that hitting a target
size requires a gel size-selection step after adapter ligation. No wastewater study
found reports a size-selection window, a fragmentation time, or an insert distribution.

## 3. Is the material too fragmented? Partly, and it is per-virus

Evidence that it is fragmented:

- Within this very dataset: RNA-level hybridisation capture recovered SARS-CoV-2 only as
  coverage "islands" inside probe regions while a plant virus (TBRV) was covered evenly,
  from which the authors conclude "the SARS-CoV-2 RNA is already fragmented in the
  wastewater … In contrast, the TBRV RNA would be intact" (Wyler et al., above). A
  coverage argument, not a length measurement.
- Nanopore wastewater work comparing tiled schemes concludes "shorter PCR amplification
  is more suitable for wastewater samples due to viral genome fragmentation", with the
  1,200 bp scheme giving less complete coverage (Water Research 2024,
  doi 10.1016/j.watres.2024.121623).

Evidence against, which matters because it keeps the recommendation alive:

- Capsid state, not nucleic-acid chemistry, decides survival: encapsidated viral RNA
  T90 ≈ 8–12 days in untreated wastewater near room temperature, 2.74–26 days across
  (+)ssRNA viruses, while unprotected RNA is lost fast (Wigginton group, *The Protective
  Effect of Virus Capsids on RNA and DNA Virus Genomes in Wastewater*, Environ Sci
  Technol 2023, 10.1021/acs.est.3c03814).
- Purified SARS-CoV-2 RNA spiked into raw wastewater "was strongly and quickly
  degraded" while encapsidated vRNA was quantified without significant loss for up to
  7 days; but **"less than 30% of the viral genomes on the average were under a
  protected form"** (Wurtzer et al., Water Research 198:117183, 2021).
- 1,200 nt tiled amplicons do work on municipal wastewater, with improvements from the
  400 bp scheme described as "marginal" (mSystems 2021, PMC8525555) - contiguous 1.2 kb
  templates survive in sewage.
- Our own measured tail: p99 ≈ 515 nt, max 1,158 nt. Molecules of that length are
  present and did sequence.

**Not reported anywhere found:** RIN or DV200 for any wastewater RNA extract; any
fragment-length distribution (electropherogram, N50, DV200) for wastewater RNA as
opposed to for a finished library; and any insert distribution from a capsid-enriched or
nuclease-protected sewage fraction. The last one is the experiment this recommendation
would need, and it does not exist.

## 4. The Illumina ceiling

- Illumina DNA PCR-Free Prep tunes median insert by SPRI ratio to measured 339 / 457 /
  551 bp - **550 bp is the top of Illumina's own documented tunable range** (tech note
  770-2020-004).
- Core-facility limits: "the DNA insert size(s) should not exceed 700 bp"; for NovaSeq X,
  "the average insert size should ideally be 350 bp and the 'tail' of longer fragments
  should not exceed 550 bp" (UC Davis DNA Technologies Core).
- Mechanism: "larger fragments do not cluster very well and therefore interfere with
  accurate cluster density prediction" (*Improved Protocols for Illumina Sequencing*,
  PMC3849550), plus Illumina's own statement of a bias towards clustering smaller
  fragments. A wide library is preferentially *un*-sequenced.
- ≥1 kb: Illumina patent prose states surface-mediated bridge PCR "performs poorly for
  inserts ≥1 kb". No authoritative statement was found for ExAmp on patterned
  flowcells; claims that patterned cells are worse for long libraries come only from
  forum posts and should not be cited.

## 5. Long-insert libraries for virus discovery

- **Not reported for RNA viromes.** No mate-pair, jumping, or long-insert library
  applied to a virome or to degraded environmental RNA was found. The mate-pair
  literature is bacterial/mammalian genome scaffolding from high-molecular-weight DNA
  (3–25 kb inserts), mechanically incompatible with cDNA at these lengths.
- Large inserts *are* used for virus discovery, but on the DNA side via fosmid cloning
  (inserts to 40 kb; Mizuno et al., PLoS Genetics 2013, 208 new marine phage genomes).
- The published long-range route for environmental RNA viromes is long reads rather than
  long inserts, and benchmarks favour Illumina or hybrid assembly over nanopore alone
  (*The long and short of it: benchmarking viromics using Illumina, Nanopore and
  PacBio*, Microbial Genomics 2024).

## 6. What this changes in our own claims

1. **Withdraw** "wastewater RNA is too fragmented to give wide templates" as the reason
   for the measured insert distribution. The kit explains the median; fragmentation of
   the free (>70%) fraction explains the absence of a long tail beyond ~1 kb.
2. **Keep** the wide-insert idea, but as a bespoke-prep recommendation with a ceiling of
   ~550 nt and an untested dependence on the protected fraction, not as something
   obtainable from public runs.
3. **New, cheap, unoccupied measurement**: nobody has published an insert-length
   distribution for a wastewater RNA virome. Ours (`results/insert_profile/`, five
   seeds, two estimators) may be the first, and the same script run across several SRA
   runs and kits would be a contribution in its own right.
