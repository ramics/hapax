# The contig half of the residue

Companion to [`REPORT.md`](REPORT.md), which established that the residue is two piles an
order of magnitude apart in depth: about 11 points of the library unplaced at median depth 1,
and about 8.4 points on contigs the cascade could not identify at median depth 6. The depth-1
half is a separate thread's. **This is the depth-6 half — the part contig-level work can
reach — and the one place in this project where profile annotation has length to work with.**

Three results, each with the control that makes it readable. The annotation itself is the
companion thread's and lives in [`../residue_contig_annotation/`](../residue_contig_annotation/);
cited here because the two halves only mean something together.

## 1. 77% of the residue-only contigs are sequence the assembly does not hold — and 18% is

`experiments/residue_contig_overlap.py`, `contig_overlap_s{1..5}.json`.

The 3,695 contigs of `REPORT.md` §2 were assembled from pairs *selected* for not matching any
of this sample's own contigs at 40 bases, so they ought to be new sequence. Ought to is not
measured: a read can miss a contig at that threshold and still come from the same molecule —
the assembly can stop short of it, or the read can sit over a variant. Aligned back against
all 82,947 own contigs, scoring the **union of alignment intervals** so a contig re-derived
in two pieces still counts as covered:

| | contigs | bases | |
|---|---|---|---|
| covered < 10% | **2,832** | 950,937 | new sequence |
| covered 10–50% | 213 | 85,288 | partial |
| covered ≥ 50% | **650** | 270,770 | re-derivation; 10.5% of all 3,695 are ≥90% covered |

| control, same pass | mean covered | ≥0.5 | <0.1 |
|---|---|---|---|
| residue-only contigs (the test) | 0.172 | 17.6% | 76.6% |
| an own-unidentified contig vs the rest of the same assembly, itself excluded | 0.001–0.004 | 0.1–0.4% | 99.6–99.9% |
| the same contigs shuffled | 0.000 | 0.0% | 100.0% |

The control is the load-bearing part. Contigs of one assembly share repeats and conserved
regions with each other, so an overlap figure alone is unreadable — and it reads 0.1–0.4%
across five seeds, far below the test, so **the 18% is genuine re-derivation rather than
incidental sharing.** The split is written out as `resasm_{new,partial,rederived}.fa` with the
continuous value per contig in `covered_fraction.tsv`.

**It mattered downstream twice.** The companion thread's length-matched control pool is the
sample's unidentifiable contigs, which is where the 650 re-derived ones live, so annotating
the unsplit 3,695 would have put them inside their own control — the circularity this repo
keeps catching. And the re-derived pile then became a **third control arm** that nobody had:
same assembler, same read class, already in our own assembly, which is what separates
*assembled from shallow unplaced material* from *actually new*. Both of the annotation's
headline conclusions turn on that arm.

## 2. The residue is LESS concentrated than known material, but its head is deep

`experiments/residue_contig_mass.py`, `contig_mass_s1.json`. Deterministic over all 8,731,224
pairs — every mate counted, nothing sampled — so there is no spread to report and the script
says so rather than shipping five copies of one number.

| | residue on unidentifiable contigs | identified reads on identified contigs (control) |
|---|---|---|
| mate observations | 1,478,246 | 5,768,467 |
| contigs carrying any | 46,157 | 22,247 |
| contigs for 50% of the mass | **762** | **112** |
| contigs for 90% | 22,646 | 2,026 |
| Gini | 0.729 | 0.922 |
| top 100 share | 29.8% | 48.5% |

Concentration is a property of any metagenome, so this number without that control says
nothing. With it, the answer runs against expectation: **known material in this library is
dominated by a handful of very abundant things and the residue is not.** There is no small
set that is most of the residue.

The head is nonetheless where the length and depth are. The 200 highest-mass unidentifiable
contigs (`top_unidentifiable_contigs.tsv`) carry **36.9%** of the residue's contig mass at
median 660 nt and median 57× SPAdes coverage; 59 are above 100×, 9 above 500×, and the leader
is 1,004 nt at 748× with 33,432 residue mate observations. The top twenty alone carry 1.37%
of every mate in the library.

And the reach is far better than the per-contig figure suggests. SPEC.md §3 reports only
3.87% of unidentifiable contigs clearing 500 nt and 2×, which read as *tier 2 is thin*; but
read mass concentrates on the long deep ones, so **42.6% of the residue's contig mass sits on
contigs ≥500 nt and ≥2×**, and 18.6% on contigs ≥1 kb. An order of magnitude more of the
residue is annotatable than the contig count implies.

## 3. Depth and recognisability are anti-correlated — and it is not ribosomal leakage

The annotation (companion thread, 5 seeds) found the two piles pointing opposite ways:

- **`resasm_new` looks viral**: 23.06% carry a Pfam family at its gathering threshold and
  11.69% an RdRp family, against 0.70% for length-matched unidentifiable contigs, dominated
  by RNA_replicase_B (1,615 hits pooled against 83). The re-derived arm reads 30.16% and
  19.47%, close to identified material, which is what makes the new pile legible.
- **The top 200 by mass are the reverse**: coding density 0.604 against 0.818 for
  length-matched unidentifiable contigs, 61% half-coding against 90%, Pfam 7.50% against
  11.17%. The deeper the unnameable sequence, the less it looks like a gene.

The first reading of that was ribosomal sequence SILVA had missed, on 6 Rfam rRNA hits against
2 in the control. `experiments/top_contig_rrna_check.py` tested it against the threshold that
defined the pile — a contig was called ribosomal at ≥100 matching bases on SILVA SSU+LSU, so
missed rRNA must pile up just under that line. 5 seeds, matching bases, with the contigs the
split **did** call ribosomal as the known-positive arm:

| pile | n | median | zero | ≥40 | ≥100 |
|---|---|---|---|---|---|
| top 200 by residue mass | 200 | 0 | 97.5% | 2.5% | **0.0%** |
| `resasm_new` | 2,832 | 0 | 99.1% | 0.4% | 0.2% |
| unidentified, length-matched | 194–2,832 | 0 | 98.5–100% | 0–1.0% | 0.0% |
| identified non-rRNA | 199–2,832 | 0 | 99.2–100% | 0–0.7% | 0.0% |
| **true rRNA contigs** | 194–2,832 | **262–499** | **0.0%** | **100%** | **100%** |

The known-positive clears the threshold on every contig at 262–499 median matching bases; the
query piles read 0 median with 97.5–99.1% at exactly zero and essentially nothing near 100.
Two populations, not a close call. **The ribosomal reading was withdrawn** (companion thread,
`5f0d3ef`), and 6 counts against 2 at n=200 was never a mechanism.

### What this check does NOT establish, and an arm that failed

minimap2 `sr` answers *did the split's threshold nearly catch this* and nothing more. It is
tuned for short reads at high identity, so **divergent ribosomal RNA is invisible to it** and
neither run licenses "not ribosomal at all". An `asm20` arm was run for that question and
**failed as an instrument**: the known-positive rRNA arm fell from 100% over the line to
49.5%, so the control moved with the test. It is committed as that
(`resasm_new_rrna_asm20_s1.json`) rather than dropped, because the failure is the reason the
claim is scoped the way it is.

The remaining tension closed from the other side, and as a kill. Split by model, ribosomal
CMs fire on 1.80% of `resasm_new` against 0.62% of identified non-rRNA contigs and 0.54% of
unidentifiable ones — so the false-positive rate on **real non-ribosomal sequence** at these
lengths is about 0.6%, which is the calibration the pooled column lacked, and a dinucleotide
shuffle's 0.00% is a null and not that rate. But the re-derived arm reads **2.51%, higher
than the residue itself**. The excess tracks *assembled from shallow unplaced pairs*, not
*new*, so it is not an annotation gap in the cascade. Structural reach versus model false
positives on genuinely divergent rRNA is still undistinguished and needs a sensitive
blastn-class search neither container has.

## 4. What this says about the deliverable

- **The depth-6 half is in far better shape than the depth-1 half.** 42.6% of its mass is on
  contigs a profile search can work on, and on the new pile that search returns a real signal
  against a measured control (11.69% RdRp against 0.70%).
- **Ranking the residue by abundance is worth keeping, and not for the reason first offered.**
  It does not surface rRNA leaking past an operating point. It surfaces long, deep, non-coding
  sequence that nothing can name — which is a better reason to rank.
- **The deliverable for this half is a ranked shortlist, not a classifier.** 762 contigs for
  half the mass is too many to read; the 200 in `top_unidentifiable_contigs.tsv` are 36.9% of
  it, abundant, long enough to annotate, and unidentified — which is the form Imo's recut of
  the target asks for, a shortlist for a skilled reader rather than a per-read verdict.
- **Open.** Whether that non-coding head is structured RNA, and the profile-reach question
  above. Both need an instrument neither container has.

## Provenance

| file | script |
|---|---|
| `contig_overlap_s{1..5}.json`, `covered_fraction.tsv`, `resasm_{new,partial,rederived}.fa` | `experiments/residue_contig_overlap.py` |
| `contig_mass_s1.json`, `top_unidentifiable_contigs.tsv` | `experiments/residue_contig_mass.py` |
| `top_mass_unidentifiable_rrna_s{1..5}.json`, `resasm_new_rrna_s{1..5}.json`, `resasm_new_rrna_asm20_s1.json` | `experiments/top_contig_rrna_check.py` |

Annotation numbers quoted in §3 are `../residue_contig_annotation/` (companion thread). Each
result carries its seed, git SHA, dirty flag and script content hash.
