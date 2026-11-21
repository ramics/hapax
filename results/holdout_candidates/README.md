# A post-hoc holdout: candidate search and design

**Status: candidate search complete, run not started.** This directory is a shortlist and
a design. Nothing here says the cascade finds anything; it says which target would make
that question answerable, what it would cost, and — the part that changed during the
search — which band the target has to sit in.

Script: [`experiments/holdout_candidates.py`](../../experiments/holdout_candidates.py).
Table: [`SHORTLIST.md`](SHORTLIST.md), generated, not typed. Every number below comes from
a file in this directory.

---

## 1. Why

`results/spike_holdout/` is the only positive control this project has, and its blind is
not blind. We chose the genome, the divergence ladder and the library, so it answers
*would tier 1 remove a read at 20% divergence* and cannot answer *would we have found
something nobody had seen*. Its clean holdout arm is one rhinovirus, and at depth 1 that
is a handful of pairs.

The archive supplies a real blind for nothing. **If a viral species' first deposited
genome postdates an untargeted shotgun RNA metagenome of wastewater that contains it, then
at the moment that library was sequenced the species was unknown to every database on
earth.** Withholding that species, and everything deposited after the run's collection
date, reconstructs the state of knowledge at the time of the run exactly. The experiment
becomes a historical counterfactual: run the committed cascade over those reads with the
clock turned back and ask whether the identity-to-closest-relative readout surfaces the
thing at its true distance, at its true depth, without being told it is there.

No spike, no chosen divergence, no simulator. The deposition dates do the work.

## 2. What was searched

**Runs** — ENA inventory of paired, untargeted wastewater RNA libraries. The census is
**73,186** runs; **14,621** in **483** studies survive as usable hosts
([`runs_by_study.tsv`](runs_by_study.tsv), [`runs_summary.json`](runs_summary.json), and
[`runs.tsv.gz`](runs.tsv.gz) per run). Amplicon and probe capture are excluded by construction —
what they contain was chosen in advance, so a species unknown at the time could not have
been targeted — and so is single-agent surveillance sequencing, which passes every
structural filter because it is neither amplicon nor capture, yet an index holdout cannot
make such a library untargeted after the fact.

**Deposits** — NCBI nucleotide records of RNA viruses whose source text says wastewater or
sewage, deposited since 2022-01-01 ([`pool.tsv`](pool.tsv)), with the first-deposit date
of every taxon they name ([`first_deposit.json`](first_deposit.json)). Of **523** such
taxa, **37** have no deposit of any kind before the cutoff.

**The join that matters** is [`links.tsv`](links.tsv): each record's GenBank DBLINK gives
its BioProject and BioSample, and from those the **run that contains it**. A candidate
with a linked public run has guaranteed presence — its genome was assembled from those
very reads — and a measurable depth. A candidate without one needs an argument that the
species happened to be in some other sewershed on some other date, which is the kind of
argument this project does not accept, and 44 of the 77 candidates are flagged
`no_containing_run` for exactly that reason.

Starting from runs and asking what was deposited under the same BioProject was tried and
is worse: it found five projects and neither of the two real leads, because a run's taxon
and strategy labels are noisy (§6) while a deposit's `isolation_source` and
`collection_date` were written by whoever held the sample.

## 3. How divergence is measured, and in which direction it is wrong

Against **viral RefSeq**, nucleotide and protein, which is the index the pipeline itself
consults (SPEC §2: DIAMOND RefSeq viral proteins), with the candidate's own taxon
excluded. Nucleotide identity is taken at three sensitivities — `megablast`,
`dc-megablast`, `blastn` — because the pipeline's tiers die in that order and the gap
between them *is* the band where nucleotide identification is failing
([`hits_*.tsv`](.)).

Every figure carries its alignment coverage and E-value. This repo has had to retract an
identity claim that turned out to be chance matching over 40 bases, and identity is
quantised at short alignment lengths, so a bare percentage is unreadable.

**The reference is today's, not the run's.** That is deliberate and one-sided: a reference
holding *more* than existed at the collection date can only make a candidate look **less**
novel. A candidate measured at 27% divergence today was at least 27% divergent then. Read
every number here as a lower bound on divergence and an upper bound on recognisability.

The exact as-of-date alternative was tried and fails in the dangerous direction. **A
RefSeq record's `createdate` is when the RefSeq entry was built, not when the sequence was
first deposited**: 3,926 of 19,650 viral RefSeq records carry a 2023 `createdate` while
wrapping sequence decades older. Filtering on it would have deleted genuinely old
relatives and inflated every divergence figure. So the holdout for the real run must key
its negative lists to the **original GenBank accession's** deposit date, never to the
RefSeq record's. NCBI's remote BLAST, which accepts a `PDAT` range as an Entrez filter and
would have done this directly, returns `Problem with taxonomic lookup and/or entrez` on
that query, so it is not an option either.

## 4. The finding that reordered the shortlist

**The project's 10–20% band is a nucleotide band, and it does not transfer to real
candidates.** It comes from `results/spike_holdout/`, where a genome is mutated uniformly,
which degrades nucleotide and protein identity together. Natural divergence does not
behave that way, and the shortlist shows both failure modes:

- Three candidates sit at **10.4–11.7% nucleotide divergence** — squarely in the band —
  while being **97.1–97.7% identical at protein level**. Their distance is almost entirely
  synonymous. The protein tier, which does essentially all the work above 5% divergence
  (SPEC §4.2: 91.3 / 84.8 / 20.0 on the holdout), identifies them immediately.
- Four candidates sit at **27–30% nucleotide divergence**, past the nucleotide wall
  entirely, yet at **60.8–69.4% protein identity over 80–96% of the genome**, which is
  where the protein tier is marginal rather than blind.

So the band has to be stated at the level the tier works at. Ranking is on protein
identity first (60–75% over at least half the genome), nucleotide second, and both columns
ship. **Eight** candidates are in the nucleotide band and **four** in the protein band
with no flags — and they are not the same four.

The top four, all with public containing runs, all first deposited 2023-07-31 from
libraries collected May–August 2022, are divergent ssRNA phage and one unclassified virus
at 60.8–69.4% protein identity. They are the ones that test whether the cascade can see
something at the edge of its reach. The 10–12% nucleotide trio tests something different
and still worth knowing: whether the reported distance is *right* when the thing is
findable.

**And depth settles which one to run.** `StochSRVP_11` (OQ722337), at 66.5% protein
identity over 80% of a 4,147 nt genome and 30.3% divergent in nucleotide, carries **32,575
mapped reads at 452× mean coverage** in its own library and **71,360 reads at 866×** in a
sibling library collected three months later ([`depth.tsv`](depth.tsv)). That is four
orders of magnitude more material than the rhinovirus holdout arm it replaces, and it
removes the one weakness this search was commissioned to fix. The other three sit at
12–19× and a few hundred reads, which is thin but real; the 17% nucleotide pair reads 63×
and 1,943 reads.

A caution that travels with the table: where `same rel` is 0, the nucleotide and protein
best hits name **different organisms**, so the nucleotide-to-protein divergence ratio is
comparing two relatives rather than two levels, and must not be read as a synonymous
fraction. It is 0 for all four leading candidates.

## 5. What the measurement killed

**Renamings, not new sequence.** Several taxa whose first deposit is recent are
**98.0–99.9% identical** to viruses already in viral RefSeq over the whole genome,
including one that is 99.86% identical to a RefSeq entry bearing almost the same name.
They are new *names* for old sequence. From the deposition dates alone they were perfect
candidates. **This is why the date is a filter and never a score: a new taxon is not a new
sequence, and only the alignment can tell them apart.**

**No containing run** — 44 candidates, including a novel sewage tobamovirus at 64.7%
protein identity which would otherwise be interesting, and every one of the recently named
enterovirus taxa.

**The best date structure in the whole search, killed by molecule and depth.** Five novel
genomes deposited in February 2025 from wastewater collected in March 2015 — a ten-year
gap, better than anything else found — come from DNA libraries (WGS and MDA) of 0.85–2.3 M
pairs. A DNA library cannot reliably hold an RNA virus, and 1 M pairs cannot hold anything
at depth. Their nucleotide hits cover 0.6–5% of the query, which is a chance match, not a
relative.

## 6. Labels are evidence, and so are taxon ids

Six things went wrong during this search in ways worth carrying forward, because every one
of them produced a plausible-looking number.

1. **Taxon ids written from memory, twice.** Four wastewater taxids: three were wrong (a
   crustacean, a bacterium and a liverwort) and the inventory they defined silently
   contained bat and bird viromes. Then the protein reference was scoped to orders chosen
   the same way, and the four candidates closest to the protein wall came back with **no
   protein relative at all** — a striking finding, and an artefact of a database that did
   not contain their class. The ids in the script are now the ones the ENA taxonomy service
   returns, and the protein reference is the whole RefSeq viral release.
2. **A coverage flag that punished the right candidates.** `partial` was defined on the
   single best HSP, but one HSP cannot cover a multi-gene genome: a levivirus encodes
   maturation, coat and replicase, so its best blastx HSP covers about 45% of the genome.
   That flag demoted precisely the candidates this search exists to find. Coverage is now
   the union over HSPs at both levels, and the two levels are flagged separately, because a
   nucleotide alignment covering 65% while protein covers 88% is **the signature of a
   divergent relative, not of a fragment**.
3. **Berlin SRR23973193, our own dataset, is filed as `TRANSCRIPTOMIC` rather than
   `METATRANSCRIPTOMIC` and carries no collection date in ENA at all** — its date survives
   only inside the sample title. A query on either field alone misses it.
4. **Conversely**, the deep-sewage candidates' libraries are filed
   `METAGENOMIC / OTHER / unspecified`, yet the paper behind them reports sixteen
   full-length RNA virus genomes assembled from those same libraries. So labels are
   collected and emitted as columns (`rna_evidence`, `cutoff_source`, `run_labels`), never
   silently trusted.
5. **ENA's `filereport` takes one accession.** A comma-separated list returns either a 500
   or an empty table with a 200, and the second is the worse failure: it looks like a
   negative result. The batchable form is the `search` endpoint with OR'd
   `sample_accession`.
6. **For a whole index, take the release file, not the API.** Paging 722,823 RefSeq viral
   proteins through `efetch` was going to take about 100 minutes; the same set as a RefSeq
   release file on FTP took 6 seconds.

## 7. The run this shortlist is for

Not started. The design, so it can be argued with before it costs anything:

**Holdout construction, per index.** No rebuild is needed for most of it. `blastn` takes
`-negative_seqidlist` and `-negative_taxidlist` at search time; DIAMOND takes
`--taxon-exclude`; the protein indexes in SPEC §2 are small enough to rebuild from
filtered FASTA in minutes. **Kraken2 has no negative-list option**, so there the holdout
is a post-filter on per-read taxid calls, and it ships with the direction of its error
stated: a post-filter removes a call after the fact, so anything the species pulled in
through a shared ancestor stays removed and the residue is **under**-stated, which flatters
the cascade. Exclusion is by deposit date on the original GenBank accession (§3), applied
both to the species and to everything deposited after the run's collection date.

**Merge the mates first** (SPEC §4.3) — free, and recognition at 20% divergence goes
10.5% → 38.4%. **Five seeds**, min–max reported as an observed range, not a dispersion
estimate.

**But the leading libraries have no second mate, and that is a real cost.** All three are
filed `library_layout=PAIRED`, and the archive holds one record per spot with no `/2`
reads: mean read length 144.99 nt against `base_count/read_count` = 144.8, so the file is
the whole library, not half of it. So SPEC §4.3's merge — the single most valuable free
step this project has — **cannot be applied to this candidate at all**, and neither can
anything that uses mate disagreement or mate linkage. The 145 nt single reads are longer
than Berlin's 107/93 and recover part of it, but the merged read's gain is a *test* rather
than a longer ruler (SPEC §4.3), and a longer single read is not that test. Any comparison
to Berlin numbers must say so. The 10–12% nucleotide trio's libraries *are* genuinely
paired (2 × 150 nt, 231 M pairs, 69 Gbase), so the merge applies there — at the price of
30 GB of staging per library, which is an EC2 job and not a download.

**The leak check is the first thing to read.** A holdout that leaks shows up as
suspiciously perfect identification, so the run reports, before anything else, the identity
of the candidate's best surviving hit. If that is ~99%, the holdout failed and nothing
downstream is readable.

**What counts as which answer.** The readout is identity-to-closest-relative with coverage
beside it, not a verdict, so success is not "it was flagged". Success is: the candidate's
reads survive tier 1, assemble in tier 2, and the reported distance to their nearest
surviving relative matches the distance measured here from the finished genome — the
instrument reads the true distance from reads alone. The informative failure is the reads
surviving but the distance being unreadable, which is a statement about where the wall is.
The expensive failure is the candidate's reads being removed by a relative at 10–20%,
which would say tier 1 eats novel material at exactly the distance we care about.

**Three controls come free with this design.** The three libraries of the leading study
were collected on different dates, so each candidate genome mapped against the *other* two
runs is a presence control — recurrence or a single shedding event, measured rather than
assumed ([`depth.tsv`](depth.tsv) is per run for this reason, never pooled). It separates
them already: `StochSRVP_11` reads 32,575 / 71,360 / 13,518 across the three dates and is
resident, while `StochSRVP_12` reads 954 / 18 / 0 from the same sample and is not. That
distinction is the one `results/recurrence/` says the depth-1 wall turns on, and here it
costs nothing. The 10–12%
nucleotide trio is a known-positive arm: it *should* be identified, and if it is not, the
instrument is broken rather than blind. And the renamings of §5 are a leak detector: they
are 98–99.9% identical to RefSeq entries, so after the holdout is applied they must still
be identified, because their relatives predate everything.

**Depth is the binding constraint, not divergence.** The residue this project cares about
sits at median depth 1, so a candidate present at a handful of pairs reproduces the
weakness of the spike instead of fixing it. [`depth.tsv`](depth.tsv) is measured by
mapping, not projected, and it is the column that decides which candidate is run first.
It has already done so: the leading four live in libraries of 4.3–10.6 M reads, small
enough to pull into a container, and one of them reads 452×. The 10–12% trio lives in
libraries of 168–231 M pairs, which is a staging job for EC2 rather than a download, so
its depth is the one number in this shortlist that is still unmeasured.

## 8. Cost

The search cost nothing: ENA and NCBI metadata, a 590 MB viral RefSeq nucleotide fetch, a
107 MB protein release file, local BLAST and local mapping, all inside the container. Read
staging for the depth measurement ran at ~1.3 MB/s, which is the one real constraint on
which candidate can be measured here and which needs EC2.

The run itself needs the chosen library staged and the cascade executed, which is EC2 work
well inside the pre-approved limit; the figure goes in the message that proposes it, not
here.
