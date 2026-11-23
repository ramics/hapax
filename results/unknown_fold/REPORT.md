# Structure search on the pile sequence cannot name

Does a structure-alphabet search place genes that sequence search leaves unidentified, and
if it does, onto what? **It places none of them**, and that answer comes with an operating
point, three controls and a measured reason.

Every number here is in [`ec2_5seed/TABLES.md`](ec2_5seed/TABLES.md), generated from the
committed run directory by `experiments/report_unknown_fold.py`. Nothing in this file is
typed by hand.

## The question, and why it is not the reach arm again

The reach arm (`results/reach/`, PR #38) asked how far structure sees on genes whose relatives are
already known: on 24 divergent genes it placed 8 where sequence placed 22, and the 8 were a
subset. That measured structure's *reach on named material*. It did not measure what
structure does where sequence returns no identity at all, which is what this project's
residue consists of. "The unknown pile differs" was the premise to test, not to inherit.

The target is tier 2 of the cascade: the 46,157 contigs that assembled from Berlin
SRR23973193 and that identification could not name. The comparator is the 25,209 identified
contigs from the same assembly of the same water, rRNA removed.

## What was done

ProstT5 writes Foldseek's 3Di structural alphabet directly from sequence, so there is no
folding and no GPU. The subject library is the whole solved databank made non-redundant by
RCSB's own 30%-identity clustering: 272,663 clusters, 40,066 distinct entries, **368,095
chains**. Virus-only would have been the wrong library, because the nameable part of this
water is 76% bacterial against 8% viral, so a virus-only library can answer "is this like a
solved virus" and cannot answer "is this like anything".

Genes come from pyrodigal in metagenome mode at a 50 aa floor. Five seeds are disjoint
blocks of 400 contigs per pile, from one fixed shuffle. **Placement means a Foldseek hit at
E < 1e-5**, and the search runs prefiltered at `-s 9.5`, so every rate here is a floor.

Because the query side is written from sequence it has no C-alpha, so **no TM-score and no
LDDT is reportable anywhere in this arm**. Similarity here is an alignment score over a
20-letter alphabet, not a superposition.

Cost: 8,178 s of encoding across 8 independent shards on one `r6i.2xlarge`, 548,959
residues at 67.1 residues/s, plus 66 s of search. $2.58 at worst case.

## Three controls, read first

*The null.* A dinucleotide-preserving shuffle of the same contigs through the identical gene
caller and the identical search, from the shared `shuffle_control` module, with its census.
**It reads 0.00% on all five seeds of both arms**, and no shuffle failed or returned an
unchanged sequence, so the null is a real null and not the fallback bug of 2026-11-22.

*The comparator.* The identified contigs place at **7.78-14.09%** per seed. Without that a
rate on the unnameable pile would be unreadable.

*The length match.* The two piles' genes are not the same length — median 78 aa against 104
— and alignment significance is a function of query length. In the comparator, placement
runs from **1.6%** in the 50-80 aa band to **35.2%** past 200 aa, a twenty-fold range. The
length effect is therefore larger than anything else in this arm, so every comparison is
reported inside length bands and the headline is the comparator's per-band rates reweighted
onto the unnameable pile's own band mix: **3.14-6.51%**.

## The result

**Zero placements.** Pooled over five seeds the unnameable pile offers **1,636 genes** and
places **0**. The length-matched comparator predicts **70.5**. Zero of 1,636 puts the
one-sided 95% upper bound on this pile's placement rate at **0.18%**, below the matched
comparator's 3.14-6.51%, so the deficit is not an absence of chances. The nearest miss in
the whole pile is a single gene at **E = 1.016e-05**, just the wrong side of the threshold.

**But zero at E < 1e-5 is not silence.** Un-thresholded (TABLES §2b), the unnameable arm
sits between both nulls and the comparator on every column:

| arm | genes | aligned at all | best E, median | best E seen |
|---|---|---|---|---|
| identified | 1,889 | 82.4% | 1.1 | 6.4e-41 |
| **unnameable** | **1,636** | **63.6%** | **2.8** | **1.0e-05** |
| identified, shuffled | 708 | 45.1% | 4.7 | 0.045 |
| unnameable, shuffled | 535 | 33.6% | 5.8 | 0.15 |

The ordering is the same on all four columns and the real arm beats even the *comparator's*
null. That is what makes the two statements consistent: this is a signal that does not clear
its own significance statistics, not the absence of one. A distance, which is what this
project means by novelty.

**The comparator's table is coherent, which is the evidence the instrument works.**
(TABLES §4) Its 205 placements onto 108 distinct chains are dominated by polymerases —
picobirnavirus and Norwalk RdRp, foot-and-mouth 3D-pol, influenza PB1, a tetramycovirus
RdRp, tobacco mosaic replicase — alongside reverse transcriptases and retroelements, phage
coat and structural proteins, and bacterial defence systems. The unnameable pile's table is
empty. Each row names the chain that aligned, its polymer type and its source organism as
three separate fields.

## Why, and what it costs to fix

Most of this answer is **gene length**. Three quarters of the unnameable pile's genes are
under 100 aa, where the comparator itself places only 1.6-4.5%, and only 37 of 1,636 are
past 200 aa, where the comparator places 35.2%. Short genes do not carry enough 3Di to
clear the statistics, exactly as 145 nt reads do not at nucleotide level
(`results/holdout_run/`). The band-matched deficit — 0 against 70.5 — is what remains after
length is accounted for, and it is real.

That makes this arm a statement about the pile's **contigs**, not about structure. The
lever is assembly, not the encoder or the library: this pile's median contig is 254 nt. The
agreed order in `SPEC.md` already points there — subtract across plants, assemble what
survives, then point structure at the contig residual — and this run is the measurement that
says structure should wait for the longer contigs rather than be pointed at these.

## What this does not say

A subject's description and source organism are a database record's words about a database
entry. They say what a placed gene is *like*. They are not a statement about what is in the
water and never a statement about what anything infects.

The library is what has been crystallised, not a sample of fold space, and the search is
prefiltered. Both make every rate a floor. **A floor of zero is still a floor** — a more
sensitive search or a larger library could move it, and this run does not bound by how much.

## A bug this run found, in its own output

The first readout resolved a Foldseek subject through its PDB *entry* and took the first
polymer entity the data API returned. A Foldseek target is a *chain*, and a solved entry is
usually several different molecules, so the names were being chosen by array order. In this
run's own comparator table a reverse transcriptase was reported as "28S DNA top strand", an
influenza polymerase as an RNA oligomer, and uL13 as L21; eight of the top fifteen rows
changed once it was fixed. Placement rates never consulted that field and are unaffected —
it was the "and like what" half that was wrong, which is the half this arm exists for.

The fix resolves per chain (`experiments/rcsb_chains.py`, tests in
`tests/test_rcsb_chains.py`) and raises where it cannot. Three traps there all failed
silently before:

- Foldseek writes the *auth* chain while the API indexes instances by mmCIF *label*. In 9TS7
  label "A" is the 26S rRNA and auth "A" is ribosomal protein S2.
- Nucleic entities had to be excluded, which turns the polymer-type census into a **check on
  the lookup**: a structure library holds protein chains only, so anything but 100% protein
  means a chain resolved to the wrong molecule. It read 33 of 1,590 onto rRNA and DNA
  records before; it now reads 3,419 of 3,419 protein.
- NMR ensembles write the chain last (`2MTV_MODEL_10_A`), so splitting on the first
  underscore gave 40 subjects a chain of "MODEL".

Same family as the shuffle fallback and the chimaera test that wrote 0.0: a default where a
refusal belongs.

## Provenance

Five seeds under one commit, `7ea6a47`, with `script_sha256` `3b22a60dcc0ce0ea` and
`code_sha256` `03ff44fb23072f7d` in `ec2_5seed/manifest.json`. The manifest carries
`git_sha_source: "launcher"`, because an EC2 instance has no checkout and no git binary, so
the commit is the one its launcher built the payload from; the content hashes pin the code
regardless.

An independent re-run of the same five seeds is **in flight in-container**, on different
hardware against a **separately built copy** of the subject library: 368,085 chains against
this run's 368,095, because two downloads of the databank days apart differ by ten chains.
Its gene census already matches this run digit for digit — 4,768 genes, 548,959 residues,
1,636/535/1,889/708 across the four arms — and its `script_sha256` and `code_sha256` are
identical, so the two ran the same code over the same material. It will be added here as
`container_5seed/` when it finishes.

Note that block membership depends on the `--seeds` *list* rather than on the seed value,
because the slice index is the position in that list. A run of a subset of seeds is
therefore not a subset of the same sample, and `blocks.json` is the cheap check before
treating two runs as comparable.
