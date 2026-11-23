# The pair as the unit a search scores, and two things it found about the holdout

**Scripts:** [`experiments/pair_search.py`](../../experiments/pair_search.py) (arms, modes,
searches), [`experiments/report_pair_search.py`](../../experiments/report_pair_search.py)
(tables). **Tables, generated:** [`berlin/TABLES.md`](berlin/TABLES.md),
[`holdout/TABLES.md`](holdout/TABLES.md). **Manifests:** `*/manifest_near.json`, carrying
`script_sha256` and the git SHA.

Provenance note, because it is the kind this repo has been bitten by: the hit tables were
produced by commit `b70a689`, which committed the script before the run. The union mode
(`best_or_fisher`), the `strict` and `positional` stages and the tables came in `db2fffe`
afterwards. The search invocations are byte-identical between the two, so every cached
search is what the later code would have run; the two later stages ran their own searches
under `db2fffe`. Nothing in the tables depends on code that did not exist when the data
were made.

---

## 1. The question, and what in the record is not it

Many search tools advertise a paired mode, where two mates are handed over as one unit.
Three things in this repo look adjacent and are not that:

- **Merging** overlapping mates into one longer read (`results/profile_reach/`). Needs
  overlap, covers 35.1% of Berlin's pairs, and it works: 10.5% -> 38.4% recognition at 20%
  divergence. The pair becomes a longer *single* query.
- **The pair as evidence** (`results/crossplant/`, `results/residue_structure/`): both
  mates unwitnessed against exactly one. The search still scores one read at a time.
- **Mate width** (`results/mate_width/`): what a wider insert buys a per-read flag.

The unmeasured thing was the **non-overlapping pair as the unit whose significance the
search computes**. That is what this directory measures, and the answer is that it does not
pay.

## 2. Design

**Primary block, `berlin/`: real pairs, nothing constructed.** 1,000 non-overlapping pairs
per pile per seed from SRR23973193, five seeds as five disjoint 1,500,000-stride blocks.
Piles come from the committed whole-library label table, so `identified` (both mates on this
sample's own identified contigs) and `residue` (neither mate on SILVA or any own contig) mean
what `results/residue_census/` made them mean. Overlapping pairs are excluded with
`profile_reach.merge_real_pair`, the committed merger, so non-overlap is the complement of a
criterion already measured here; **58-65% of a pile survives**, consistent with that
directory's 35% merge rate. Mate 2 is trimmed of its 14 nt non-template prefix. Three index
levels: date-held-out viral RefSeq nucleotide, viral protein, and Swiss-Prot, the last
because a wastewater residue pile is mostly not viral and a viral index is the wrong
denominator for it.

**Secondary block, `holdout/`: the known distance.** The holdout library is single-end in
practice, so a pair cannot be taken from it; mate 1 is a real read at its real position and
only its partner is constructed, from Berlin's measured insert distribution conditioned on
the non-overlapping band and Berlin's measured per-mate error rates. **Direction of the bias,
with no control that removes it:** the constructed partner carries no library artefact and no
within-population variation, so it is cleaner than a real mate, and any gain a joint rule
shows there is an upper bound. Two uniform-sampling arms (`candidate_sim`, `rung17_sim`) are
kept because the gap between them and the real arm is a result in its own right (§5).

**Nulls.** 10,000 shuffled pairs per block, both mates shuffled independently with the shared
`shuffle_control` module, census recorded: 0 dropped identical anywhere. Shuffling both mates
while keeping the pairing is the right null for a joint rule - it leaves the rule its
structure and takes away the sequence.

**Two operating points, and why the second one is not per mode.** Every mode is read at
E < 1e-5, this repo's point, and at a `common` point: the strictest E-value at which the
*baseline's* null is still silent, applied to every mode. A per-mode matched-null threshold
is not usable here, because several joint modes' nulls never fire at all (0 of 10,000), which
leaves their threshold unbounded rather than strict - reading each mode at its own unbounded
point would hand the joint modes a looser test than the baseline got.

## 3. The answer: every joint mode is net negative

Pair rate at E < 1e-5, five seeds, min-max, from [`berlin/TABLES.md`](berlin/TABLES.md):

| level | pile | best single mate | concatenated | concatenated + spacer | Fisher per subject | both mates agree |
|---|---|---|---|---|---|---|
| nucleotide | identified | **26.40-29.70** | 26.00-29.30 | 26.00-29.20 | 19.90-23.20 | 18.50-20.40 |
| nucleotide | residue | **5.90-10.60** | 5.80-10.30 | 5.40-10.30 | 2.40-5.10 | 1.60-2.80 |
| viral protein | identified | **39.10-43.60** | 36.80-41.60 | 35.90-41.60 | 26.00-30.40 | 14.20-18.70 |
| viral protein | residue | **13.30-17.20** | 12.40-16.30 | 11.80-15.80 | 8.10-11.30 | 2.40-3.70 |
| Swiss-Prot | identified | **16.70-19.30** | 16.60-19.30 | 15.70-18.90 | 14.30-16.20 | 11.60-14.00 |
| Swiss-Prot | residue | **3.60-5.20** | 3.80-4.60 | 2.80-4.40 | 1.80-3.30 | 0.30-1.30 |

The null reads 0.00% in every cell but one: a single nucleotide-level pair in 10,000, which
bounds chance at about 1 in 10,000 and says nothing below that.

**The table that settles it is conversion, not rate.** Of the pairs the single mate misses,
each joint mode rescues **0-8 per 1,000**; of the pairs it catches, the modes lose **10 to
261 per 1,000**. Net negative in every cell of every level. A rate column alone would hide
that, because a mode can hold its rate while trading one set of pairs for another.

**`agree` cannot add sensitivity at all, by construction.** Its statistic is the *worse*
mate's E-value, so it is a strict subset of the single-mate rule. Its only possible value is
precision, and the null leaves no room for that - it is already silent. A write-up quoting a
both-mates rule as a sensitivity gain has misread it; `tests/test_pair_search.py` asserts
the subset relation so nothing here can drift back into that reading.

## 4. Why concatenation loses, measured at the join

An E-value scales with the query's length while the alignment does not, so concatenating two
mates doubles the length and leaves the HSP where it was - unless an alignment genuinely
spans the join. The junction census says whether any does. Every HSP of every concatenated
query, Berlin block, nucleotide level:

| | inside mate 1 | inside mate 2 | crossing the join |
|---|---|---|---|
| concatenated | 166,980 | 142,731 | **41,820** |
| concatenated, 10 nt N spacer | 171,942 | 148,117 | **7** |

With a separator, **7 of ~320,000 alignments cross** - so essentially nothing real spans the
join, and a concatenated pair is one mate's alignment charged for two mates' length. Without
a separator, 41,820 cross and **2,979 of them win their query**: the search invents an
alignment across the join, a chimaera of the experimenter's own making. The holdout block
reads the same way, 5,936 best HSPs crossing without a spacer and none with one.

**So "feed the tool both mates" is not a neutral convenience.** Unseparated concatenation
manufactures hits that are not in either mate, and separated concatenation is a pure
significance penalty.

## 5. The one honest positive, and its size

The Fisher combination's null never fires - 0 of 10,000 pairs, at every protein level - so
taken as an *addition* rather than a replacement (`best_or_fisher`: whichever of the baseline
and the joint statistic is better) it cannot cost sensitivity and costs nothing measurable in
false positives. What it buys:

| level | pile | gain /1,000 pairs | loss /1,000 | null |
|---|---|---|---|---|
| nucleotide | identified / residue | 0-2 / 0-2 | 0 | 0.00% |
| viral protein | identified / residue | 2-5 / **1-8** | 0 | 0.00% |
| Swiss-Prot | identified / residue | 1-3 / 0-4 | 0 | 0.00% |
| holdout, real reads at 30.3% divergence | candidate | **0** | 0 | 0.00% |

**At most 8 pairs in 1,000, and exactly zero on the blind.** That is the whole of what
pairing adds to a search, and it is the number to hold against the cost: a joint mode needs
mate identity carried through the whole search path, and requiring non-overlap discards
35-42% of the library's pairs to the merge criterion.

**Recommendation for SPEC, and now in it** (§4.3 and §5.1 item 2, this branch): keep the
pair out of the search and keep the merge. The pair
earns its place where this repo already puts it - as evidence after scoring, and as a longer
single query where the mates overlap.

## 6. What this ran into: `results/holdout_run/`'s nucleotide index was not viral

Rebuilding the date-held-out indexes for §2 produced an index that disagreed with the
committed one, and chasing the disagreement found two things. Both are in the direction of
"the published number is not measuring what it says", and both were invisible to the gate
that exists to catch exactly this.

**(a) The index contained bacterial chromosomes.** `results/holdout_run/REPORT.md` names its
nucleotide index as "RefSeq viral nucleotide (19,650 records)". Its own committed hit tables
name *Haemophilus influenzae* (40 reads on seed 1), *Candidatus Hamiltonella defensa* (21),
*Staphylococcus epidermidis* (6), *Lactobacillus johnsonii*, *Staphylococcus aureus*,
*Streptococcus pyogenes* as the best subjects of its candidate arm. Of the 122 reads that run
reports as identified at nucleotide level, **81 have a non-viral best subject**; across five
seeds the split is 67-69%. None of those records is in the RefSeq viral genomic release at
all, and the run's manifest cannot say what the database was built from: that invocation ran
`--stages reads,search,report` and its `nt_fa` argument is empty. Re-searching those 122
reads against an index rebuilt from the viral release, **71 lose their hit entirely**. The
protein side is clean - all 48 protein hits on seed 1 are viral.

**The count gate could not see this.** `holdout_indexes.py` refuses to build unless it
reproduces six recorded counts, and it reproduced the protein counts exactly while the
nucleotide count differed by 9 records. A total is not a membership check. The cheap test
that does work, and which this directory now demonstrates: take the previous run's committed
hit tables, re-search its identified queries against the rebuilt index, and compare which
queries survive.

**(b) The candidate arm is 94% a partial-match halo, not reads of the organism.** That arm is
defined as every read whose primary minimap2 `sr` alignment goes to the held-out genome, and
is described as "reads of the held-out virus, 32,575 of them at 452x". Measured from the same
mapping: those alignments span a **median of 48 reference bases out of a 150 nt read**, and
only **1,856 of 32,575 (5.7%)** span 80% or more of their read. **31,061 of them - 95.4% -
start at one position**, base 70 of 4,147, a window whose best nucleotide hit is a 31 nt exact
match to a relative's own 5' end and which has no significant protein hit. The arm is a
pile-up on one short conserved 5' element.

The control that makes this unarguable is in the same mapping pass: the `leak` genome, a real
abundant organism 0.08% divergent, has **24,290 of 45,265 reads (53.7%)** spanning 80% or
more, median span 134. A genuine organism gives end-to-end alignments. So the dilution is a
property of that arm, not of the mapper, and the genome's end-to-end depth is about 67x
rather than 452x.

### The corrected read-level numbers

Reads aligning over >=80% of their length, whole eligible pool (no subsampling, so a count
rather than a five-seed spread), both indexes the viral release minus the same 439 post-cutoff
accessions, E < 1e-5. From [`holdout/strict.json`](holdout/strict.json):

| arm | true distance | reads | nucleotide | viral protein | shuffled null |
|---|---|---|---|---|---|
| `candidate` | **30.3% nt divergent** | 1,856 | **18.91%** | **29.63%** | 0.00% / 0.00% |
| `leak` | 0.08% divergent | 24,290 | 99.75% | 88.28% | 0.01% / 0.00% |

Two published claims move:

1. **The read-level wall at 30% divergence is softer than recorded.** The published figure is
   5.25-7.00% at nucleotide level, with the gloss "a single 145 nt read survives 30%
   divergence about one time in sixteen". On reads that actually belong to the organism it is
   about **one in five** at nucleotide level and nearly **one in three** at protein. The
   published number is low because of the halo and high because of the index, and the two
   errors are not of the same size.
2. **"Nucleotide beats protein at read level at every distance" does not hold at this rung.**
   On indexes that are both viral, protein is ahead at 30.3% divergence (29.63 against 18.91)
   and nucleotide is ahead at 0.08% (99.75 against 88.28). The honest statement is a
   **crossover**, not an ordering: nucleotide wins near, protein wins far. SPEC §5.1 item 2
   now carries the crossover instead of the ordering, and §7 carries the depth note below.

3. **A depth figure on this genome is not even coverage.** The recorded 452x is aligned-base
   depth over a read set 95.4% of which starts at one position; end-to-end depth over the
   body of the genome is about 67x. The contig-level facts other directories quote from this
   library - the 3,990 nt contig covering 96% of the genome, and its protein identification -
   do not depend on the faulty nucleotide index, which `results/holdout_cascade/` never
   consulted (it uses Kraken2, SILVA and a protein tier), but any claim that reads a depth
   number as uniform support does depend on the profile in
   [`holdout/positional_profile.tsv`](holdout/positional_profile.tsv).

**What is NOT re-derived here, and should be.** The same arm carries that directory's
per-read identity bias (+12.8 to +13.6 points at nucleotide level), its coverage-binned table
and its p90 identity figures. Those were computed on the diluted arm and need re-running on
the strict one before they are quoted again. This directory does not do it, because it is not
this experiment's question; it is recorded as an open item so nobody quotes those numbers in
the belief that they have been checked.

## 7. A warning about simulated reads, from the two arms kept for it

Pairs sampled uniformly off the candidate genome identify at **65.0-68.5%** per pair at
nucleotide level. Real reads of the same genome identify at 2.4-3.5% per pair on the diluted
arm and 18.91% per read on the strict one. So uniform simulation overstates the real strict
arm by roughly **two-fold**, and overstates the arm as the published run defined it by about
twenty-fold. `holdout/positional.json` separates the two causes: 49.6% of the genome's 145 nt
windows have a significant nucleotide hit, so recognisability is concentrated, and the real
reads' start positions are not spread over the genome at all.

**So a simulated-read experiment is measuring a different, easier instrument than a real
library**, and the size of the difference depends on the pile-up structure of the arm it is
compared against. Any ladder built on simulated reads - `results/spike_holdout/` included -
carries this.

## 8. Caveats, in the order they would bite

1. **The Berlin block carries no known distance.** Nothing in that water has a true answer
   attached, so it answers the mode question (does a joint rule clear what neither mate
   clears) and not a distance question. The holdout block carries the distance and a
   constructed partner.
2. **Non-overlap is a lower bound on separation.** A pair the committed merger cannot merge
   is either genuinely non-overlapping or overlapping with more mismatches than the criterion
   allows.
3. **Five seeds are five blocks of one flow cell** for the Berlin block and five subsamples
   of one pool for the holdout block. The min-max is an observed range; nothing is estimated
   from its spread.
4. **The protein levels are DIAMOND `--very-sensitive`,** not BLASTX. `results/reach/` found
   DIAMOND at three sensitivities, BLASTP and full Smith-Waterman returning the same nearest
   relative on all 16 contigs, so the substitution is one this repo has calibrated; it is
   still a different instrument, and no comparison here crosses levels.
5. **The index is today's viral release minus the committed 439 accessions**, not a
   re-derivation of the date rule. No record in today's release derives from the study's
   accessions or carries its organism names, so the holdout is intact for this candidate, but
   a post-cutoff record newly added to the release would not be covered by that list.
6. **The strict arm's filter is a choice with a direction**: requiring an end-to-end alignment
   to the deposited consensus could discard genuine reads of a variant whose ends diverge,
   which would make the strict arm read high. A 48 nt median is far too short for that to be
   the main effect, and the leak arm's 53.7% is what the filter does to reads that genuinely
   are the organism's.

## 9. Cost

Zero. In-container throughout: local BLAST+ 2.17.0 and DIAMOND 2.1.9, the RefSeq viral
release, Swiss-Prot, and the two staged libraries. 285 searches, about 2.5 hours of wall
time. No EC2.

## 10. Re-running this

```
# indexes (committed accession list applied to today's release; --allow-mismatch is
# required and is itself part of finding (a) above)
python3 -I experiments/holdout_indexes.py --seed 1 --allow-mismatch \
  --nt-fa viral.1.1.genomic.fna.gz --prot-fa viral.1.protein.faa.gz \
  --gbff viral.1.genomic.gbff.gz --blastbin $BIN --out $S/holdout
diamond makedb --in $S/holdout/holdout_prot.fa -d $S/holdout/db_holdout_prot_dmnd
diamond makedb --in uniprot_sprot.fasta.gz -d $S/holdout/db_sprot_dmnd

# primary block
python3 -I experiments/pair_search.py --seed 1 2 3 4 5 --band near --pairs 1000 \
  --arms berlin_identified berlin_residue shuffled --levels nt prot_dmnd sprot_dmnd \
  --berlin-r1 SRR23973193_1.fastq.gz --berlin-r2 SRR23973193_2.fastq.gz \
  --blastbin $BIN --dbdir $S/holdout --out results/pair_search/berlin

# holdout block, then the two stages that chased the instrument findings
python3 -I experiments/pair_search.py --seed 1 2 3 4 5 --band near --pairs 1000 \
  --arms candidate leak candidate_sim rung17_sim shuffled --levels nt prot_dmnd \
  --fastq SRR23421544.fastq.gz --blastbin $BIN --dbdir $S/holdout \
  --out results/pair_search/holdout
python3 -I experiments/pair_search.py --stages strict,positional --seed 1 \
  --blastbin $BIN --dbdir $S/holdout --out results/pair_search/holdout

python3 -I experiments/report_pair_search.py --dir results/pair_search/berlin
python3 -I experiments/report_pair_search.py --dir results/pair_search/holdout
```

The hit tables are gzipped per block (`pairs_near.tsv.gz`, one row per pair per level with
every mode's E-value). The per-query BLAST output and the mapped pool file (which holds library
sequence) are gitignored and rebuilt by the `reads` and `strict` stages; what is committed
instead is `holdout/library_pool_geometry.tsv.gz`, the alignment geometry of all 77,840
pooled reads without their sequence, which is what both instrument findings in §6 are read
off, and `holdout/positional_profile.tsv`, the per-window recognisability and read-start
profile behind §7.
