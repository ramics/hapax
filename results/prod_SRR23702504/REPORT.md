# Production over a library nobody had measured: SRR23702504

**Yeongdok plant effluent, South Korea, collected 2020-09-11, RNA, 2 x 151 nt, 32,125,614 pairs
(PRJNA941081).**

The claim being tested is not that the pipeline produces good numbers on this water. It is that
**the pipeline runs end to end on a sample nobody has touched without anyone editing a flag**, and
the deliverable is the honest list of what had to be supplied by hand. A refusal that was correct
is recorded here as a success, not a defect.

*Drift is the largest caveat on every figure in this repo: the same plant, the same assay, eight
months apart moves the identifiable fraction from 22.7% to 40.7% (`results/dna_rna_pair/`). This is
one water on one date. Nothing here is a statement about Yeongdok, about Korea, or about effluent
in general.*

## Why this sample

Three libraries in this repo cannot test genericness, by construction:

- **SRR23973193** (Berlin influent) is the discriminating dataset. Every per-library fact in
  `experiments/libraries.json` was measured on it, so running it proves only that the defaults it
  set are still there.
- **SRR23421544** (the holdout) is the positive control. Its value is the withheld genome, and its
  facts are likewise committed.
- **PRJNA746354** (New York) is a correctness fixture.

What was wanted was a library from a different plant and ideally a different continent, with no
entry in `experiments/libraries.json`, **and not a member of the cross-plant panel** — a library
must not be subtracted against a panel it is inside, because its own reads are abundant in that
panel by construction (54% false growth against a 6% swap control, `results/crossplant/`); the
stage refuses that case and the refusal is correct.

SRR23702504 satisfies all of it and keeps one thing the panel needs: it is **not** in
`results/crossplant/panel.tsv`, while **its plant is** — 12 panel libraries from Yeongdok, three
influent from the same date — so the paired within-plant stratum the cross-plant measurement
requires still exists. `tests/test_prod_run.py` asserts both halves of that, because the
not-a-member half is the one refusal that would invent recurrence.

It is also, deliberately, **effluent** rather than influent, and that turned out to be the single
most consequential difference between it and Berlin (see the trim, below).

## What the run refused, and what had to be supplied

| # | Where | What it refused | What I supplied | Was the refusal right? |
|---|---|---|---|---|
| 1 | `library_facts.resolve_args` | every stage, on an accession with no entry | a committed entry in `experiments/libraries.json`: pairs, read length, block stride, read paths, each with its basis | Yes. This is the whole point of the resolver raising, and it is the one refusal that is unavoidable for a new library. |
| 2 | `stage_refs.resolve` | the truth-panel and positive-control stages | `truth_panel: null` with `truth_panel_absent_because` (no depositor genome set for this water) and the same for the spiked control | Yes. **Absent with a reason is not zero recovery**, and the earlier design would have reported a zero. |
| 3 | `experiments/trim_from_prefix.py` | naming a mate-2 trim at all | the value `0` stated on the launcher's `--trim-r2`, with a paragraph of basis and `trim_r2_unverified: true` | Yes, and it is the most interesting result of the pass. See below. |
| 4 | `hapax.py run` preflight (local, against dummy inputs) | the whole run, twice, before any money was spent | the Rust binaries (`linkage`, `census`) built on the instance, and the relabelled assembly written under the file name `crossplant_stage.py` requires | Yes. Both would have been discovered eleven hours and six dollars later. |

### The trim is the thing that could not be inherited

The one fact production cannot carry across libraries is how much non-template 5' sequence to cut
from mate 2. Berlin's is 14 nt, measured; a silent default of 14 on a library with no prefix
deletes fourteen real bases from every mate, in the direction of less recognition.

The committed instrument (`experiments/nontemplate_prefix_map.py`) maps untrimmed reads against
SILVA and reads the soft-clip at the 5' end, with mate 1 and the 3' end of the same mate as
controls. It therefore needs ribosomal content. On this library, over 5 x 100,000 pairs:

- mate-2 5' clip **mode 0 in 5 of 5 seeds**, holding 47.1–52.5% of anchored mates
- mate 1 also 0, at 73.3–78.5%
- resting on **127–144 anchored mates per seed: 0.13–0.14% of pairs**, against roughly 32,000 and
  ~32% on Berlin
- non-zero clips spread flat over 3–73 nt at 1–3% each; 3' clip mean 11.0–13.6 nt

The shape is a library with no prefix. The support is 139 reads. `trim_from_prefix.py` refused on
its own anchored floor (2,000), which is the correct call: at that depth the histogram describes
SILVA's reach into treated effluent, not the library's chemistry. So the trim is **stated, not
measured** — SPEC 4.1 allows stating it because stating it is the override, and a stated value is
never re-checked. `0` is chosen because it is the only value that cannot delete real template; the
residual risk, if this kit does carry a prefix, is kit sequence left in the reads, which costs
recognition rather than inventing it.

The by-product is a third rRNA instrument, committed as `rrna_anchored_fraction`: this water is two
to three orders of magnitude less ribosomal than the Berlin influent.

**The better instrument is obvious and untested**: clip the reads against the sample's *own*
assembled contigs, which needs no SILVA at all. It needs an untrimmed assembly to measure against,
so it is a chicken and egg this pass did not break.

## What I edited, and why each was code rather than data

The measure of genericness is how much had to change in code. Four things did:

1. **`experiments/nontemplate_prefix_map.py`** — knew three nicknames and nothing else, so a new
   accession could not be measured at all. It now resolves any accession through
   `library_facts`. This is the one edit that was unavoidable: a script that measures a library had
   no way to be pointed at one.
2. **`experiments/provenance.py`** — raises rather than writing "unknown", which is what it is for,
   but an instance has a tarball and no `git`, so every script importing it raised there. It now
   takes the commit from `HAPAX_GIT_SHA`, which the launcher sets, **only when git cannot answer**,
   so a machine with a checkout cannot dress a dirty tree as clean.
3. **Three new scripts rather than edits**: `trim_from_prefix.py` (the five seed files to an integer
   or a refusal), `trim_mate2.py` (a trimmed copy for the assembler, since MEGAHIT reads files and
   `path_settings.py` explicitly refuses to trim), `megahit_rename.py` (MEGAHIT headers to the
   SPAdes-style names the classification stages parse).
4. **`experiments/prod_run.sh`** — the ordering from `holdout_cascade_run.sh` with every
   sample-specific value required from the environment (`${X:?...}`) instead of defaulted. Nothing
   in it names a library.

Nothing in `hapax.py` or in any scoring stage was touched.

## What the tool could not show while it ran

The assembly is the longest step in the pass and the only one with no visible progress.
`path_settings.py` captures MEGAHIT's stdout, so from outside the instance the log shows the
command and then nothing until a contig file exists. Over an hour-plus step that is the moment a
progress line is worth most, and its absence is indistinguishable from a hang. It is the same
family as the refusals above — an empty output and a step that never ran look identical — one
layer out: the run says nothing rather than saying it is working. Noted here rather than fixed,
because fixing it mid-run would have stamped a new script hash onto this run's output.

## What the generator could not build

(to be filled once the run directory is committed)

## Cost and wall clock

| Launch | Commit | Result | Wall | Cost |
|---|---|---|---|---|
| i-0c35504358c8f69e3 | ed01fc3 | `pip install` exited 0, `import mappy` failed under `-I`; `--quiet` hid the reason | ~3 min | ~$0.03 |
| i-0d20d59fe31946d78 | c76d6a5 | Kraken2 smoke test piped through `head` under `pipefail`: SIGPIPE read as a broken binary | ~90 s | ~$0.02 |
| i-00007d94fdc55f541 | bb0f48b | provenance raised with no `git`, and the shell reported a refusal that had not happened — plus the real prefix measurement | ~13 min | ~$0.07 |
| i-0e44a2d889f41d119 | 90060fa | the production pass | (pending) | projected worst case $5.86 |

The three dead launches cost about **$0.12** in total. Every one of them was a failure in staging,
not in the pipeline, and each is now a test or an assertion.

## Provenance

- Launch manifests: `launch_i-*.json`, one per launch, keyed by the commit the payload was built
  from. The launcher refuses to ship a dirty tree, so the commit identifies the code.
- `pip_freeze.txt`, `tool_versions.txt`: the environment the pass ran in.
- `nontemplate_prefix/nontemplate_prefix_s{1..5}.json`: the trim measurement, which refused.
- Tables in `TABLES.md` are generated by `experiments/report_prod.py` from the committed run
  directories. No number in this file is typed by hand that is not in those files.
