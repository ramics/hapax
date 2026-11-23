# Provenance debts: the committed record's own weak points

Every entry here is already measured and recorded somewhere in the tree. What was missing
was one place that lists them, so that a reader of a headline can tell in one look whether
its provenance is clean. The convention they are read against is in
[`CLAUDE.md`](../CLAUDE.md): a manifest carries the seed, the git SHA, the parameters and a
content hash of the script, a `-dirty` SHA is unresolved provenance rather than provenance,
and a provenance change forces a re-run.

**None of these debts moves a number.** Where one could have, that is said explicitly and
the direction is given. A withdrawn number stays in its table with its item closed; nothing
here is a quiet deletion.

| # | debt | state | moves a number? |
|---|---|---|---|
| 1 | Manifests name SILVA in `parameters` for stages that never consulted it | recorded, open | no — the census never read it |
| 2 | `subtract_known/`'s five seeds sit under two commits, three of the five `-dirty` | **closed by a clean re-run**: `subtract_recut/` | no — the headline is the clean directory's |
| 3 | 47 scripts import shared code whose hash their manifests never recorded | fixed forward in `provenance.py`, not retrofitted | no |
| 4 | `divergent_known` and `divergent_contig_workup` mapped back at a mate-2 trim of 16 where the measured prefix is 14 | pinned with the discrepancy stated; **decision below** | no load-bearing number; see 4 |
| 5 | Three scripts default their read paths to `data/ww_berlin/`, which neither `libraries.json` nor `data/MANIFEST.md` documents | **found 2026-11-23, open, not fixed here** | no — it fails loudly |
| 6 | A *correct* dinucleotide shuffle is not interchangeable with the canonical one, so "switch it to the shared module" is a re-run | measured and pinned by a test | no — but it would if switched blind |

## 1. SILVA is named in manifests that did not read it

`results/fail_loudly_audit/audit.json` records the verdict in the stage's own words:
*"numbers unaffected (the census never read SILVA); these manifests name an index in
`parameters` that was not consulted"*. It is a manifest-writing defect, not a measurement
defect: `parameters` was written from the argument namespace, so an index the run was
*offered* appears whether or not a tier reached it.

It is open on purpose rather than retrofitted. Rewriting a committed manifest to remove a
key would change the file that pins a result without re-running the result, which is the
one thing the convention exists to prevent. The forward fix is the consulted-index record
that PR #30 added — since then a run states what it actually read, and a stage that cannot
reach a named index refuses to start rather than skipping it silently. So a manifest written
before PR #30 is read as "these indexes were available", and one written after it as "these
indexes were consulted".

**What a reader should do with it:** take the index list in a pre-PR-#30 manifest as an upper
bound on what was consulted, and read the tier table in the same directory for what actually
removed anything. A tier that removed nothing and a tier that never ran look identical in
`parameters` and are distinguishable in the curve.

## 2. The 77.10% headline comes from the clean directory, and it is the recut

This one is closed, and the closing is worth stating because the two directories are easy to
confuse.

| directory | git SHA per seed | `script_sha256` |
|---|---|---|
| `subtract_known/` | `c7dfe7f9…-dirty` (s1, s2), `74ea377c…` (s3), `74ea377c…-dirty` (s4, s5) | absent |
| `subtract_recut/` | `1278233f…` for all five seeds | absent |

So `subtract_known/` is exactly the mixture `CLAUDE.md` says is not a five-seed result: two
commits, three trees that were not committed, and no script hash to pin the code whatever
the tree looked like. `subtract_recut/` is five seeds under one clean commit.

**The 77.10% (77.02–77.28) union is the recut's number**, and `subtract_known/REPORT.md`
§table quotes it from `subtract_recut/TABLES.md` rather than from its own run. SPEC §10 names
both directories on that row because the report and the numbers live apart; the number itself
has one source.

Neither directory carries `script_sha256`, because both predate the field. That is debt 3.

## 3. The import-hash gap is fixed forward, not retrofitted

`subtract_known.provenance` hashed only the calling script, and 47 scripts import `provenance`,
`kraken2`, `map_filter`, `_mm_index` and `dir_bytes` from it. A change to any of that shared
code was therefore invisible to the caller's `script_sha256`: two runs could carry the same
hash and different code.

`experiments/provenance.py` (PR #35) closes it, adding `modules_sha256` — every module loaded
from this repo, hashed — and `code_sha256` over the script and those modules together. Two
runs with equal `code_sha256` ran the same Python from this repo.

**Nothing in `results/` carries `modules_sha256` yet**, because switching a script to the
shared module changes its `script_sha256`, and under the re-run convention that means
re-running the directory it wrote. The ranked list of what is still copied, with each copy's
measured behaviour, is `experiments/SHARED.md`; each script switches when it is next touched
for its own reason. New code imports, and a tripwire test fails on a new copy.

**What a reader should do with it:** a `script_sha256` on a pre-PR-#35 directory pins the
entry script and not the shared code underneath it. Where that matters, the git SHA is the
stronger field, and a clean one pins everything.

## 4. The mate-2 trim of 16 in the two divergent stages: the decision is not to re-run

`divergent_known.py` and `divergent_contig_workup.py` trimmed 16 nt off mate 2 when mapping
reads back onto contigs, where the measured non-template prefix is 14
(`results/nontemplate_prefix/`, mode exactly 14 in 5/5 seeds on both RNA libraries). It is
pinned as `trim_r2_for_mapback` in `experiments/libraries.json` with the discrepancy stated,
rather than reconciled behind the results' back.

**Asked to decide rather than leave it pinned forever: do not re-run, and retire the pin when
either stage next re-runs for its own reason.** The reasoning, so the decision can be
overturned on evidence rather than re-argued:

- **What it can touch is depth and evenness only.** The trim applies to the map-back stage.
  Identity, aligned length, subject and query coverage, the gene-level spread, the role
  ordering and both artefact controls come from the protein search over contig sequence,
  which never sees a trimmed read.
- **Nothing load-bearing is a depth number.** Of the rows SPEC §10 draws from these two
  directories, none is a depth or an evenness figure. Depth appears in the description of the
  six or seven near-finished genomes (171–467x) and in the evenness field of
  `divergent_contigs/` §4; both are read as bands, and a 2-base over-trim on one mate of a
  107/93 pair is a ~2% uniform reduction in mapped bases. It cannot reorder contigs, because
  it applies to every read in the arm equally.
- **It is the conservative direction.** Over-trimming removes real template, so every depth
  and breadth figure in those two directories is a *floor*. A re-run at 14 would move them up
  slightly, which is the direction that makes the near-finished genomes look better, not
  worse.
- **The cost is not the objection.** `divergent_known` is 951 s on 4 cores and fits
  in-container; the re-run is affordable. The objection is that it would move committed
  numbers for no gain in what those numbers are used for, and a re-run that changes a figure
  buys a reconciliation debt of its own in every report and figure that quotes it.

What *does* change here is one line of prose: `divergent_contigs/REPORT.md` said mate 2 was
trimmed 16 nt "as everywhere in this repo", which was true of neither the repo nor the measured
prefix. Corrected in place, with the discrepancy named.

## 5. Three scripts defaulted to a data path the manifest does not document — resolved, and the premise corrected

`absence_curve.py`, `absence_residual.py` and `rrna_match_source.py` default `--r1`/`--r2`
to `data/ww_berlin/SRR23973193_{1,2}.fastq.gz`. Every other script, `libraries.json`'s `r1`
and `r2`, and every Berlin entry in `data/MANIFEST.md` say `data/berlin/`. Both spellings
appear in committed manifests — 62 occurrences of `data/berlin/` against 18 of
`data/ww_berlin/` — so the directory was renamed at some point and three defaults were left
behind.

**It cannot corrupt a number.** A wrong path is a missing file, so the run fails rather than
reading the wrong library; this is the loud kind of defect, which is why it has survived. The
three scripts do resolve `trim_r2` through `library_facts`, and because the stale path still
carries the right accession it satisfies `paths_name_library`, so the refusal that catches a
*different* library's reads correctly does not fire here.

**Resolved as recommended, 2026-11-23.** A read path is a fact about a sequencing run, so
`library_facts.resolve_args` now fills `r1`/`r2` from `libraries.json` where they appear in
its `keys`, and the three scripts carry `default=None`. Nothing points at a path this file
chose: Berlin's entry already recorded `data/berlin/SRR23973193_{1,2}.fastq.gz`, which is
the spelling 62 committed manifests use. Two refusals come with it, both asserted in
`tests/test_library_facts.py`: **agreement is decided before any path is filled**, so giving
mate 1 of one run and letting mate 2 default refuses rather than handing over the named
library's mate 2; and a library whose entry has no `r1` refuses rather than borrowing
Berlin's, which is the holdout library's case. A resolved path is made absolute against the
repo root, because the same command is run from scratch directories and EC2 mounts.

**And the premise of this entry was half wrong, which is worth more than the fix.** The
directory was not simply renamed. The READS are documented under `data/berlin/` (62 against
18), but every committed manifest that names an ASSEMBLY subset says
`data/ww_berlin/spades_sample*.fa` — 11 to 13 occurrences, and **none** under
`data/berlin/`. So `absence_residual.py --ref` and `rrna_contigs.py --contigs` are not
stale: they name the spelling their own results were produced under, and "correcting" them
to match the reads would have pointed them at a path nothing documents. They are still
sample-specific defaults and belong here as an `assembly` fact, which needs somebody who
can see `data/` on disk — the one thing this container cannot check is which spelling
exists. **Open, scoped, and smaller than it was.**

**What moved as a consequence.** The `script_sha256` of `absence_curve.py`,
`absence_residual.py`, `rrna_match_source.py` and `library_facts.py` changed, so
`results/absence_curve/`, `results/absence_residual/` and `results/rrna_match_source/` now
carry a hash that does not match the committed code. No number moved — a read path cannot
corrupt one, because a wrong path is a missing file — and this is recorded rather than
regenerated, the same decision as §3 and the same precedent as the trim-and-stride pass.

## 6. Not a debt, but the trap next to one: a correct copy is still not interchangeable

Worth a section because a tidy-up is exactly when this gets broken, and because it was nearly
broken during this one.

`holdout_run.dinuc_shuffle` is **correct**. Over 400 real sequences from this repo's corpus it
preserves length, dinucleotide composition and both endpoints on every one, and never fails.
So does `shuffle_control.dinuc_shuffle`. They agree on **0 of 400**, because Altschul-Erikson
admits many valid walks and the two consume the generator differently.

So the instinct that follows from `CLAUDE.md`'s rule — a shuffle lives in one module and is
imported, never re-typed — has a measured exception for the *existing* dinucleotide copies:
switching `reach_queries.py` or `holdout_run.py` to the canonical module **redraws every
shuffled sequence** and so moves the null in `results/reach/` and `results/holdout_run/`.
Both of those nulls read 0.00%, so the **conclusion** survives a redraw and any **count**
does not. That makes it a re-run, not a switch.

This is the opposite of the mononucleotide case, where
`test_mono_copies_draw_identically` proves the switch changes nothing but the hash — which is
why `experiments/SHARED.md` could say a switch is free there and cannot say it here.

**Correctness is what to check on a copy; interchangeability is a separate question, and here
the answer is no.** Pinned by
`tests/test_shuffle_control.py::test_dinuc_copies_do_not_draw_identically`, which fails if a
walk is ever changed in a way that makes them agree — at which point the committed nulls that
used the old one need re-running.
