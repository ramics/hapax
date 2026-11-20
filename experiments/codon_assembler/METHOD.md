# Codon-aware assembly for viral populations

**Status:** build spec. Narrowed from earlier drafts to one component.

---

## 1. The problem

A viral population in a sample is a swarm of related variants, not one genome.
Nucleotide assemblers handle this badly: the graph tangles, and the result is either
many short contigs from one entity or chimaeras stitched between variants. The
field's own account of the failure is that complete viral genomes are usually not
recovered — instead you get several contigs from a single entity, some matching
nothing known.

Most of the variation causing the tangle is synonymous. Third-position changes carry
no protein consequence, and they are the bulk of within-population diversity.

## 2. Thesis

**Detect overlaps in codon space, keep nucleotides throughout.**

Working in codon space makes synonymous variation invisible, so a swarm collapses
into one entity rather than fragmenting the graph. Keeping the nucleotide sequence
preserves gene order, intergenic regions and read linkage — everything protein-level
assembly throws away.

What this does **not** do: join reads that lack nucleotide identity. Reads from one
molecule match each other by definition, and joining ones that don't is the chimaera
failure. Divergence tolerance belongs in read-to-*reference* comparison, which
assembly never performs. Structural homology likewise belongs downstream, in deciding
what an assembled thing is.

**Where the benefit shrinks.** In regions coding in two frames at once, selection
acts on both, so there is far less synonymous variation to collapse. Those regions
are correspondingly more conserved, so nucleotide assembly already handles them. The
method degrades gracefully rather than breaking — but see §5 for the matching risk.

Prior art to build on rather than repeat: RAMICS already aligns in codon space with
explicit states separating codon-sized indels from frameshifts, and carries
platform-specific error probabilities. It is the author's.

## 3. Design

**Relative offset, rather than frame election.** Viral genomes are dense and read in
several frames at once — *gag/pol* in HIV, most of the hepatitis B genome — so there
is no single frame to elect. Assembly does not need one. Two reads from the same
locus share a frame *relative to each other*, whatever that frame encodes. So test
the three possible offsets between a read pair and take the best-scoring. No
reference, no global frame, and overlapping genes drop out of the mechanics.

**Overlap detection, two-tier for speed.**
- Seed with exact nucleotide k-mers. Cheap, and it catches the easy majority.
- Where seeding is ambiguous or fails, score with the codon HMM. Expensive, rare.

This is what makes it fast: the hard machinery runs only where the cheap machinery
is uncertain. GPU execution of the HMM, as RAMICS did, applies again on far better
hardware.

**Scoring mode chosen per region.** Viral genomes carry untranslated sequence —
UTRs, packaging signals, origins — where codon scoring means nothing. The assembler
needs nucleotide scoring alongside codon scoring, selected per region rather than
set once for a run.

**Output is variant-aware.** A single consensus discards the thing we care about.
Emit a contig plus its variant structure, so downstream work can ask whether this is
one entity or several. This is the output that makes closure scoring possible later.

## 4. Evaluation

**Simulated populations at known diversity.** EvolveAGene generates coding sequences
at controlled diversity; the author's thesis used 0.05 to 0.20 substitutions per
site. Run the same ladder.

Compare against MEGAHIT, metaSPAdes and PLASS on:

- **contigs per true entity** — the headline number; one population should yield one
  contig
- **chimaera rate** — joins between variants that never co-occurred
- **recovered length** at each diversity level
- **wall clock**

**Then real data, in order of difficulty.** HIV *env*, where the author has curated
alignments and ground truth. Then a public wastewater sample.

## 5. What would kill this

- **No gain over nucleotide assemblers** at realistic within-population diversity.
  The whole claim is the synonymous collapse; if assemblers already cope, there is
  nothing here.
- **Collapsing erases real differences in overlapping genes.** A change silent in one
  frame can be meaningful in another, so synonymous collapse can merge variants that
  differ in a protein encoded by the other frame. Measure this directly on HIV
  *gag/pol* and on hepatitis B, where most of the genome is dual-coding. If the
  method cannot be made safe there, it needs a flag marking dual-coding regions and
  falling back to nucleotide scoring across them.
- **Too slow.** If the codon HMM is invoked on most overlaps rather than a minority,
  the two-tier design has failed and the tool is unusable at metagenomic scale.
- **Chimaeras increase.** A more permissive overlap criterion joins more things
  wrongly. Measure this before claiming any improvement in contiguity, because
  longer contigs that are wrong are worse than short ones that are right.

## 6. Where this is heading

The longer arc is identifying novel viruses by whether a group of fragments covers a
complete replication cycle — closure as the identification rather than as quality
control afterwards. That needs assembled entities with their variation intact, which
is what this component produces. It is recorded here so the output format serves it,
and is otherwise out of scope until the assembler works.
