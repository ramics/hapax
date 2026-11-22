# Which committed claims rest on a null that did not shuffle

Imo, 2026-11-22: *"so the fact that our shuffling runs didn't actually shuffle — is there
anything we should rerun after discovering that?"*

The bug: a dinucleotide-preserving shuffle is an Eulerian walk over the doublet graph, and a
naive walk dead-ends. Two implementations in this repo answered a dead end with `return seq`,
which puts **real sequence into the null pile**. A null contaminated with the test can only
move one way — towards the test — so the direction of the error is fixed and readable:

- A claim of the form **"the null scores zero / far below the test, so the signal is real"**
  is *conservative* under contamination. The true separation is at least as large.
- A claim of the form **"test and null are indistinguishable, so the field is dead"** is
  *dangerous*. Contamination is exactly what makes them look alike.

Read the second kind first. There are three, and all three were already withdrawn before
this audit, by the threads that own them.

**Correction, after this file was first committed: the bug is still live on `main`.**
`experiments/divergent_known.py` defines its own naive walk with `return s` after 20 retries,
and that is the copy `main` runs; `experiments/shuffle_control.py`, the corrected module, is
on a branch and not on `main`. So a fourth claim belongs in the table below: running `main`'s
`divergent` stage over the holdout library's identified pile reported 175 shuffles of 22,110,
0.79%, carrying a significant hit. Withdrawn, same shape as the 324-of-11,361 case, and the
file records no count of how often the fallback fired, so the share of those 175 that were
never shuffled is unrecoverable from the output. Being fixed in PR #33, which imports the
canonical module rather than writing a fourth copy of it and adds the `null_census` that three
broken copies would have had to show. **The lesson this file now carries about itself**: an
audit that reads committed *results* does not see a live defect in committed *code*, and this
one did not.

## 1. Measured, not read off the code

`experiments/shuffle_fallback_audit.py --seed 1..5 --fasta <resasm contigs>`: every
implementation still in the tree, run over the same real sequence, cut into read-sized
windows (107 and 93 nt, this library's two mate lengths) and whole contigs of 300–5,000 nt.
Output `shuffle_fallback_s{1..5}.json`. The column is the fraction of inputs the function
handed back unchanged.

| implementation | dead-end answer | read windows | whole contigs |
|---|---|---|---|
| `holdout_run.dinuc_shuffle` (the fixed walk, last-edge arborescence first) | raises | **0.00%** | **0.00%** |
| `divergent_known.dinuc_shuffle` (naive walk, 20 retries) | `return s` | **0.65–1.15%** | **0.46–0.91%** |
| `residue_structure.dinuc_shuffle` (walk, 50 retries) | `return s` | 0.00–0.10% | 0.00% |
| `residue_contig_annotation.dinuc_shuffle` (Altschul–Erikson, 100 retries) | `return s` | 0.00% | 0.00% |

Five seeds, min–max. The fixed walk reading 0.00% is the calibration that makes the other
rows readable. On Berlin's identified contigs the gene-level thread measured 2.85% for
`divergent_known` — longer, more repetitive sequence dead-ends more often, so **the rate is a
property of the input as much as of the code** and the figures here are not a universal rate.

**The rate alone understates the damage, and by how much is measurable.** The sequences that
dead-end a walk are repetitive and low-complexity, which are the same sequences that then
produce a chance alignment, so contamination is concentrated in exactly the inputs that score.
On Berlin identified contigs the 2.85% that fell through supplied **304 of 304** of the null's
significant hits. So the pair of numbers to report is the return-the-input rate beside the
share of the null's signal that came from it: 2.85% of the pile, 100% of the hits. A rate of
1% is not a 1% perturbation of a null; here it was the whole of it.

Re-run with the fixed walk, from the thread that owns that result: **0 of 11,361** shuffled
Berlin contigs produce a significant protein alignment and **0 of 932** shuffled ORFs do, the
whole library yielding one 47 aa alignment at E=8.45e-4, which clears neither floor.

### The defect is the missing tree condition, not the bare `return s`

A grep for the fallback across `experiments/` finds five local copies, and **two of them are
false positives for a reader greping on that string**. `residue_structure.py` and
`residue_contig_annotation.py` both end with a silent `return s` after their retry loop, and
both carry the Altschul–Erikson connectivity condition — each vertex's chosen *last* edge has
to keep a path to the sequence's final base — which is what makes the walk succeed. The two
broken copies lacked that condition. **The fallback made the failure silent; it did not cause
it.** So the thing to grep for is a walk with no last-edge or tree condition, and the thing to
check afterwards is the measured rate.

Measured over 7,390 read-length windows of real sequence (one RNG seed per call, the way each
script calls it):

| implementation | has the tree condition | returned input |
|---|---|---|
| `residue_contig_annotation` | yes | **0 of 7,390** |
| `residue_structure` | yes | **2 of 7,390, 0.027%** |

`residue_structure` is **not** zero, and the two failures are ordinary-looking 93 and 107 nt
reads rather than degenerate sequence. Its 50 retries cannot fix a structural failure, so for
those sequences the fallback fires on every retry. The gene-level thread measured 0 of 5,000
for the same function; both readings are consistent with a rate near 0.03%, which neither
sample size pins. **No committed number from either script moves at this magnitude**, and the
direction is the conservative one, but the file should not record it as zero, so it does not.

**The hygiene fix is a decision about re-runs, not a free tidy-up.** Making either script raise
instead of falling back, and report its own fallback count, changes its `script_sha256`, which
under this repo's provenance rule forces re-running `results/residue_structure/` and
`results/residue_contig_annotation/`. Folding `holdout_run.py`'s own correct copy into the
shared module forces `results/holdout_run/` too, 92 minutes. **Recommendation, recorded rather
than acted on: do not re-run for this.** The measured rates say nothing published moves, so the
fix should ride along the next time each of those scripts is touched for a reason of its own,
and the rates above are what makes that deferral readable instead of a shrug.

**A mononucleotide shuffle (`rng.shuffle` of a list of bases) has no walk and no dead end, so
it cannot hit this bug.** Neither can a permutation of labels, partners, columns or rows. Most
of this repo's nulls are one of those, which is why the blast radius is small.

## 2. Every committed claim that rests on a null

| claim | null it used | affected | bias direction | verdict |
|---|---|---|---|---|
| `divergent_known/` §2 — existence of a significant alignment separates real from shuffled, 89.53% against 2.68%, a factor of 33 | `divergent_known` dinuc | **yes**, 324 of 11,361 contigs returned unchanged, and they supplied all 304 of the null's hits | null inflated towards the test; the true separation is larger | **already withdrawn** — re-run gives 0 of 11,361, correction on `claude/project-thread-do0tzd`-adjacent branch `…eswrkj`, commit 457628e |
| `divergent_known/` §2 — conditional on a hit, nothing about the alignment separates real from shuffled (identity, aligned length, coverage, bitscore) | same | **yes**, and this is the dangerous shape: it compared real sequence against real sequence | the kill was manufactured by the contamination | **already withdrawn**; with a correct shuffle the control pile is empty, so these fields have no null at all, neither calibrated nor killed |
| `holdout_cascade/` — the coverage-gap field is dead, the shuffled arm reading *higher* than the real one (14.47% against 5.81%) | same | **yes** | same | **already withdrawn** (6c1fc37); `experiments/coverage_gap_null.py` now imports the fixed walk and has published nothing |
| `holdout_run/` — 5.25–7.00% of 145 nt reads find a relative at 30.3% divergence, shuffled null **0.00%**, five seeds | `holdout_run` dinuc, the fixed walk | no | n/a | **safe**, and it is the calibration the table above leans on |
| `residue_contig_annotation/` — the 2,832 new contigs carry an RdRp at 11.7% against 0.70%, shuffle null flat at zero coding density and zero Pfam | `residue_contig_annotation` dinuc | no, 0.00% measured | would have raised the null towards the test | **safe** |
| `residue_structure/` — the shuffled null is degenerate on the coverage axes (`zero_frac` 0.9999) and calibrates composition and nothing else | `residue_structure` dinuc | ≤0.10% of null reads | an unshuffled null read keeps its own 31-mers, so contamination lowers `zero_frac` towards the test | **safe at this magnitude**; and the report's one kill (`mate_disagree` is an error-floor estimator) deliberately calibrates against a known floor rather than this null |
| `residue_census/` — the artefact signatures together are under 1% and shuffled reads score zero | `residue_artefacts`, **mononucleotide** | no, no dead end exists | n/a | **safe by construction** |
| `weak_hits/`, `profile_reach/` — the weak-hit band is indistinguishable from composition-matched shuffles; with an E-value floor the shuffled hit rate is 0.00% | `profile_reach.shuffle_seq`, **mononucleotide** | no | separate caveat below | **safe from this bug** |
| `readlevel_features/`, `recurrence/` queries, `public_health_view/` 44-taxon panel | mononucleotide | no | n/a | **safe by construction** |
| `rrna_content/` (`prodigal_on_rrna`, `rrna_peptide_reach`) | numpy permutation | no | n/a | **safe by construction** |
| `mate_width/`, `peptide_absence_noise/`, `mate_friends/` | partner permutation | no | n/a | **safe by construction** |
| `recurrence/` grouping | curveball row/column swap | no | n/a | **safe by construction** |
| per-read `fold` field | Markov-1 resampling, which cannot copy its input | no | n/a | **safe by construction** |

**The separate caveat on the weak-hit kill, which is not this bug.** A mononucleotide shuffle
destroys dinucleotide composition, so it inflates stop-codon density and deflates ORFs
relative to real sequence. That makes the null *easier* to distinguish, not harder — so a kill
of the form "indistinguishable from shuffles" survives it, and the direction is stated rather
than controlled.

## 3. Priority order, nothing re-run

1. **One thing is outstanding and it is code, not a number**: `main`'s own
   `divergent_known.py` still carries the naive walk, so any run of the `divergent` stage from
   `main` produces a contaminated null today. PR #33 replaces it with the canonical
   `shuffle_control` module. Until that lands, treat any null from that stage as withdrawn on
   sight.
2. **No committed claim is outstanding.** Every claim this bug touched sits in
   `divergent_known/` and `holdout_cascade/`, and both were re-run and corrected by the
   threads that own them. `divergent_known/`'s re-run reproduced every number to the digit
   from a materially different code state (§5b there), which is the stronger reproduction
   claim this repo rarely gets to make.
3. **`divergent_known.py` is edited by the fail-loudly pass**, and the owning thread has
   agreed the order: this branch merges first, that thread takes these hunks rather than
   re-deriving them, and the single re-run happens afterwards. Its side is small and does not
   conflict semantically — the local `dinuc_shuffle` deleted for an import of
   `shuffle_control`, plus `--max-target-seqs 25` and a wider outfmt on the whole-contig
   search. (required `--dmnd-sprot`, r1/r2
   and the diamond binary; three `if not db.exists(): continue` removed). Under the
   provenance rule that forces a re-run of its committed result, and the owning thread's
   `experiments/shuffle_control.py` replaces the shuffle in the same file, so the re-run
   should happen **once, after both land**, not twice.
4. **`residue_structure/`'s null arm** is a ≤0.10% contamination pointing the conservative
   way inside a field already reported as non-discriminating. Re-running it changes a number
   nobody quotes. Worth doing when that experiment is next touched, not on its own.
5. **The one measurement nobody has made**: whether the dead-end rate is higher on the
   sequence classes this project cares about most — low-complexity and repetitive residue
   contigs are exactly where a walk dead-ends, and that is also where a chance alignment comes
   from. The 304-of-304 concentration above establishes that the two coincide on one pile; what
   is not measured is whether a *correct* shuffle's residual chance rate is itself concentrated
   in the same classes, which is the one way this bug's reading could survive the corrections.
   The control for it is the aligner's own compositional masking, still not run.

## 4. Five-seed directories: are any seeds empty, identical or short?

Same failure family, asked by the holdout thread after finding that an unexposed block stride
left two of five seeds reading nothing while the directory still looked complete. Checked
across all 150 multi-seed file groups in `results/` for a zero-byte seed, a seed under half
its group's median size, and byte-identical seeds.

**No seed is empty or short anywhere.** Two groups are byte-identical across all five seeds,
and both are correct by design, stated in their own reports:

- `block_position_bias/block_fractions_s*.npz` — the block statistics are exact over the whole
  library and do not vary; the seed varies only the random-sample control, which is in the
  per-seed JSON beside it.
- `holdout_cascade/run_s*/blind_spots.json` — a fixed declared panel copied into every run
  directory by `hapax.py`.

What this check cannot see: a seed that read the *right* number of reads from the *wrong*
region, and a seed whose output differs only in a field nobody reads. The stride fix is what
closes the first; nothing closes the second from committed output alone.
