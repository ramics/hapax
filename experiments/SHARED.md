# Shared code: what lives in one place, and what is still copied

Imo, 2026-11-22: *"writing the same code over and over is a huge smell."* It is also how
the same shuffle bug shipped twice in two days. The second writer of a walk did not know
which condition the first was fixing. This file is the index of shared modules and the
ranked list of code that is still duplicated. It is ranked by how much a copy that drifts
could distort a published number, not by how many times the code appears.

## The rule

**A control, a parser of search output, or a provenance writer lives in one module and is
imported.** The module **raises** when it cannot do its job. It never falls back to
something the caller might mistake for a result. Where a published number depends on a
behaviour, a test checks the module's own invariant by running it on real and adversarial
input. Nobody gets to pass it by reading it.

**Switching a script to a shared module changes its `script_sha256`.** CLAUDE.md ties a
results directory to the hash of the script that wrote it. So existing copies are
grandfathered. Each is listed in a test with its measured behaviour, and it switches when
its script is next touched for its own reason. Where the test proves the shared function
draws identically (the mononucleotide shuffles), the switch changes no output, only the
hash. New code imports; a tripwire test fails on a new copy.

**For the dinucleotide walks the switch is a re-run, not a no-op, and that is measured.**
`holdout_run.dinuc_shuffle` is correct — over 400 real sequences it preserves length,
dinucleotide composition and both endpoints on every one and never fails — and so is
`shuffle_control.dinuc_shuffle`. They agree on **0 of 400**, because Altschul-Erikson admits
many valid walks and the two consume the generator differently. So switching
`reach_queries.py` or `holdout_run.py` to the canonical module redraws every shuffled
sequence and moves the null in `results/reach/` and `results/holdout_run/`. Both of those
nulls read 0.00%, so the *conclusion* survives a redraw and any *count* does not.
`tests/test_shuffle_control.py::test_dinuc_copies_do_not_draw_identically` pins this, so
that a well-meant tidy-up cannot quietly become a measurement change. **Correctness is the
thing to check on a copy; interchangeability is a separate question and the answer here is
no.**

## Shared modules

| module | what | invariant test |
|---|---|---|
| `shuffle_control.py` | `dinuc_shuffle` (Altschul-Erikson walk with the tree condition, raises), `mono_shuffle`, `check_shuffle`, `shuffle_arm` (validates every output, returns a census) | `tests/test_shuffle_control.py` |
| `search_hits.py` | `read_hits` (named columns, raises on a missing file or a wrong field count), `best_per_query` (strict `>` on bitscore, as every copy does), `pair_key`, `query_coverage` and `subject_coverage` from coordinates in the query's own units, `union_query_coverage` with one E-value floor for both modes, `coverage_gap` (both unions from the same rows under one `passes_floor`, and it raises on a negative gap before any table is written), `outfmt_args` so the command and the parser share a column list | `tests/test_search_hits.py` |
| `provenance.py` | `provenance(script, git_sha_override=None)`: the same `git_sha`, `tree_dirty` and `script_sha256` as `subtract_known.provenance`, plus `modules_sha256` (every repo module loaded) and `code_sha256` over all of them. It raises rather than writing "unknown". An EC2 instance has no checkout and no git, so a launcher passes its own commit in and `git_sha_source` records which of the two answered; a `-dirty` or empty override raises | `tests/test_provenance.py` |
| `loud.py` | refuse-to-start preflight, consulted-index record, `run_tool` | `tests/test_loud.py` (on main) |
| `subtract_known.py` (de facto) | `provenance`, `kraken2`, `map_filter`, `_mm_index`, `dir_bytes`, imported by 47 scripts | none of its own; see rank 3 |

## Inventory, ranked

### 1. Controls and nulls

**Sequence shuffles: 4 dinucleotide walks and 8 mononucleotide copies** (two of the
mononucleotide copies arrive with PR #32). Measured by `tests/test_shuffle_control.py`
(`pytest -s` prints the table). Results are returned-input plus `None`, per family:

| copy | real (1,079) | dinuc repeat (400) | AT only (400) | poly-A flanked (400) |
|---|---|---|---|---|
| `shuffle_control` (canonical) | 0 | **0** | 0 | 0 |
| `holdout_run` | 0 | 21 `None`, counted by its caller | 0 | 0 |
| `residue_structure` | 0 | **78 returned input** | 0 | 1 |
| `residue_contig_annotation` | 0 | **41 returned input** | 3 | 1 |
| `divergent_known` on main | 5 | 269 | 0 | 9 |

"Real" is the 78 committed holdout-candidate genomes, whole and in 107 and 93 nt windows.
On real sequence the three copies with the tree condition are clean. That is why their
committed results stand, and it matches the shuffle audit. On a dinucleotide repeat with a
short random tail, **every local copy gives up and two hand back their input**, while the
canonical walk shuffles all 400. "Clean" is a property of the input, not of the code.
Neither committed null was run on such input, so nothing published moves. The danger is the
next person who points one of these copies at satellite or low-complexity sequence.

`divergent_known.py`'s copy is replaced by an import in PR #33.

**Why the canonical walk survives where the others give up.** It picks each vertex's
reserved last edge among *distinct* successors. The others weight by edge multiplicity, so
on `ACACAC…` they almost always reserve the dominant edge and fail the tree test until
their retries run out. The cost is stated rather than hidden: choosing among distinct
successors does not sample Eulerian walks uniformly. Composition is exact either way, so a
null built this way is still valid. It is a slightly different random draw from the one a
multiplicity-weighted walk would make, which matters only if someone compares two nulls
drawn by different walks. Don't.

**Other randomisations are already safe by construction.** These are partner and label
permutations (`mate_rho`, `cov_step_*`, `cross_library_witness`, `remote_blast_check`),
curveball swaps (`recurrence_grouping`), the Markov-1 fold resample (`residue_report.folding`),
and numpy permutations (`prodigal_on_rrna`, `rrna_peptide_reach`). None can return the test's
data under the label of a null. They are not consolidated: each permutes a different
structure, and one shared function would only hide that.

**Length-matched comparators: two `length_bin`s that are not the same function.**
`residue_contig_annotation.length_bin` uses log bins from 200 nt at a factor of 1.25.
`holdout_cascade.length_bin` uses the fixed `LENGTH_BINS` table. Both are right for their
own experiments. Matching across the two would compare different strata.

### 2. Reading search output: identity, coverage, alignment length, E-value

**About 14 copies of "best DIAMOND or BLAST hit per query by bitscore"**:
`subtract_known`, `residue_report`, `full_residue_search`, `known_strain_divergence`,
`holdout_cascade`, `divergent_known`, `profile_reach`, `classify_contigs`, `spike_holdout`,
`dna_rna_pair`, `coverage_gap_null`, `partial_alignment_anatomy`, `divergent_contig_workup`
and `bigindex_search`. They agree on the core: strict `>` on bitscore, so the first-seen
hit wins a tie. They differ in four ways that each move a number:

- **The key.** Some take the read pair (`rsplit("/")`), some the mate, some the raw query.
  A pair keyed by mate and then counted as pairs double-counts.
- **Columns by position.** Each script requests its own `--outfmt 6` column list and
  unpacks it by position. The same index means `evalue` in one script and `qstart` in
  another.
- **Whether the E-value is kept.** Three copies (`classify_contigs`, `spike_holdout`,
  `dna_rna_pair`) keep only the bitscore, or the bitscore and identity. Under this repo's
  own rule, an identity without its E-value is not quotable.
- **Coverage units, now wrong twice.** `holdout_cascade.parse_best` divided a DIAMOND
  blastx alignment length, which is in amino acids, by `qlen`, which is in nucleotides. So
  every protein label's best-HSP coverage capped at exactly 0.3333 in PR #33's
  `holdout_cascade_s1.json`. `divergent_known` computes the same quantity as
  `3 * aln_aa / contig_nt`. The holdout thread has fixed `parse_best` in PR #33. The field
  was never quoted. `search_hits.query_coverage` takes coverage from query coordinates,
  in the query's own units.
- **Two definitions of the coverage gap, and one of them went negative.** Before PR #33's
  fix, `holdout_cascade` subtracted a per-subject union over *all* HSPs from an
  any-subject union over HSPs with E < 1e-5. `coverage_gap_null` floored both. Under one
  floor the gap cannot be negative, but it was negative on 394 of 7,734 protein rows and
  508 of 2,269 nucleotide rows. The candidate's +36.6 was the same under both definitions.
  The -0.8 nucleotide figure was an artefact. Once the floors were matched, a `>` against
  `>=` difference still left one negative row in 41,271. So `search_hits` has a single
  comparison (`passes_floor`) and `coverage_gap`, which forms both unions from the same
  rows and raises on a negative gap before anything is written.

### 3. Manifest and provenance

- **`subtract_known.provenance` hashes only the calling script.** 47 scripts import
  their mapping, Kraken2 and filtering code from `subtract_known.py`. When that shared code
  changes, their `script_sha256` does not. The git SHA is the only record, and a `-dirty`
  SHA is no record. Consolidating more code makes this hole bigger, so `provenance.py`
  comes with the consolidation, not after it. New scripts use it. Existing ones keep
  `subtract_known.provenance` until they are touched, because the output keys are a
  superset of the old ones.
- About 33 scripts call `git rev-parse` themselves. There are 14 private
  `_git_sha`, `_provenance` and `manifest` helpers. All of them flag a dirty tree, but
  only some hash the script.

### 4. Statistics and small helpers (low risk, consolidate when touched)

- **`wilson`: 7 copies.** Three of them (`bigindex_report`, `dna_rna_pair`,
  `remote_blast_check`) return the interval `(0.0, 0.0)` for `n = 0`. That is a certain
  zero manufactured from no data, which is the "default where a refusal belongs" shape. The
  other three return `None`, and `block_position_bias` divides by zero. Their z values differ
  in the fifth digit (1.96 against 1.959964).
- **`read_fasta`: 18 copies.** They differ in return type (a dict or a list of pairs) and in
  whether they upper-case. A dict silently keeps the last of any duplicate names.
- `revcomp` (9), `translate` and `six_frames` (3), and the subprocess-and-time wrappers
  around DIAMOND (one per script).

### 5. Plotting

Figures are generated by about 20 `report_*.py` and `make_*.py` scripts with no shared
style or colour module. This is the lowest risk to a number and the most visible to a
reader. It is not consolidated in this pass.
