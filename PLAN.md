# Plan: from research archive to a reusable tool

## Summary

hapax today is a research archive with a thin driver on top. `hapax run` starts after the
hardest steps (assembly, splitting the assembly into identified and unidentified piles,
staging databases), which a user must do by hand, and the report renderer accepts a run
directory that is missing most of its inputs. The Yeongdok run (PR #55) finished with zero
failed stages and produced a mostly empty page, which is the expected outcome of the current
design, since nothing tests the path from raw reads to the finished report.

The fix is structural. Build one installable package whose single command goes from reads
to report, make every report section the output of a declared stage, make a missing section
an error, and gate every merge on an end-to-end test over a small fixture. Experiments move
out of the product's import path and stay as an archive.

Nothing below changes a pipeline number. Each step is one PR.

## 1. Gaps, worst first

1. No test runs reads to report. CI (`.github/workflows/tests.yml`) runs `pytest` with no
   external tools or indexes, so every test stops at a unit or a mocked stage.
2. The driver starts after assembly. `hapax.py` never calls an assembler; `cascade`
   demands `--contigs`, `--unidentified-contigs` and `--cross-sample-ref` from the user
   (`hapax.py:313`, `hapax.py:534`, `hapax.py:1533`, `hapax.py:1549`). Assembly and
   splitting live in shell and EC2 UserData (`experiments/*_userdata.sh`,
   `experiments/holdout_cascade_run.sh`).
3. Report sections with no producing stage. `make_web.py:623-638` reads
   `public_health_profile_s*.json`, `known_strain_divergence_s*.json` and `ranking.json`;
   no stage on `main` writes them (`STAGE_ORDER`, `hapax.py:1432`).
4. The renderer treats a missing input as valid. `load()` returns `None` for an absent file
   (`make_web.py:258`) and `_one()` does the same by design (`make_web.py:125`); the page
   then prints "absent" text instead of failing the build.
5. Success is measured by stage exit codes, not by output. PR #55 changes the run summary
   so a run reads as fine when stages were left out (diff of `hapax.py` around "complete
   runs: 0 of 5").
6. Fixes are split across unmerged branches that overlap. PRs #46 and #55 both edit
   `hapax.py`, `experiments/concern_ranking.py`, `experiments/libraries.json` and
   `tests/test_hapax_stages.py`, and have never been run together.
7. The product is 10 scripts inside a directory of 146. `hapax.py` shells out to
   `subtract_known.py`, `residue_report.py`, `linkage_census.py`, `spike_holdout.py`,
   `divergent_known.py`, `mate_friends.py`, `removal_shortlist.py`, `truth_panel.py`,
   `crossplant_stage.py` and itself; their dependencies on the other scripts are implicit.
8. Not installable. No `pyproject.toml`; the entry point is `python hapax.py`; the
   renderer is a separate 4,113-line script run by hand.
9. Configuration is spread across CLI flags (about 30 on `hapax run`), JSON fact files
   (`libraries.json`, `stage_references.json`) and paths typed into scripts. A new input
   needs code or JSON edits before it runs.
10. The script-hash provenance rule (CLAUDE.md, `script_sha256`) discourages editing the
    scripts that are, in practice, the product, so fixes get layered beside them.

## 2. Target architecture

- **One package.** `src/hapax/` with `pyproject.toml`, installed by `pip install .`, entry
  point `hapax`. Product code lives only there. `experiments/` keeps the research scripts
  and their results and is never imported by the package.
- **One command.** `hapax run --config run.toml READS... --out DIR` performs, in order:
  input check and trim measurement, assembly, identification and splitting of contigs,
  cascade, every analysis stage, and the report. `hapax report DIR` rebuilds the page only.
- **Stages as a registry.** Each stage is a Python function with declared inputs, declared
  output files and a JSON schema for each output. The runner builds the order from the
  declarations. A report section names the stage output it renders, and a test fails if
  any section names an output no stage declares.
- **Missing means failed.** A stage that cannot run because a database is absent is
  recorded as failed with the reason, and the report build exits non-zero unless the user
  passed `--allow-partial`, in which case the page carries a banner listing what is absent.
  "Run complete" means every declared section rendered from a validated file.
- **Config, not edited constants.** One TOML file per run: inputs, database paths, threads,
  seeds, operating points. Library facts that must be measured (trim, layout) are measured
  by a stage and written into the run directory; a stated value in config overrides it.
  `libraries.json` becomes a store of known values for the regression tests.
- **Declared dependencies and data.** Python deps in `pyproject.toml`; external binaries
  (MEGAHIT, DIAMOND, Kraken2, minimap2, BLAST) pinned in a conda `environment.yml` plus a
  container image; databases fetched by `hapax fetch-db NAME` into a configured directory
  with checksums, so EC2 UserData shrinks to "install, fetch, run".
- **Provenance stays, scoped to the package.** Manifests record the package version and git
  SHA. The per-script hash rule applies to `experiments/`; product code is versioned and
  edited normally, and a number in `results/` names the package version that produced it.
- **Archive and delete.** Move the shell runners and `ec2_*.py` launchers either into the
  package as one `hapax cloud` helper or into `experiments/archive/`. Delete driver code
  paths that exist only to accept hand-made piles once the assembly stage lands.

## 3. Steps, one PR each

1. **Reconcile the open branches.** Merge PR #46 and PR #55 into one branch, resolve the
   four shared files, revert the "complete" counter change, and run the existing suite.
   Close the other branch.
2. **End-to-end fixture test, failing on purpose.** Add a tiny fixture (a few thousand
   simulated read pairs from two short public virus genomes plus background, and toy
   databases built from them, all under a few MB in `tests/fixtures/`). Add a CI job that
   installs MEGAHIT, DIAMOND, Kraken2 and minimap2 from conda and runs the whole command,
   then asserts that every report section rendered. It fails until step 6, which makes the
   remaining work visible.
3. **Package skeleton.** `pyproject.toml`, `src/hapax/`, the `hapax` entry point, the
   stage registry and the TOML config, wrapping the existing scripts unchanged.
4. **Strict report.** Section-to-output contract, JSON schemas, missing file is a build
   error, `--allow-partial` banner. Move `make_web.py` into the package as `hapax.report`.
5. **Assembly and splitting stages.** Trim measurement, MEGAHIT, and identified/unidentified
   splitting as stages, removing the hand-supplied `--contigs`,
   `--unidentified-contigs` and `--cross-sample-ref` inputs.
6. **Adopt the three missing producers** (`profile`, `taxa`, `ranking`) as stages with
   schemas. The fixture test from step 2 goes green here.
7. **Lift script code into the package** stage by stage, replacing subprocess calls with
   imports, and add `hapax fetch-db`. Archive what the package no longer calls.
   Finished as agreed on 2026-11-24 (Imo): the driver is `hapax.driver`, run in-process by
   `hapax run`, with `hapax.py` kept as a shim for committed command lines; the package's
   glue is lifted with byte-equality tests against the scripts (PR #62); and the science
   stages stay byte-identical child processes started only through `hapax.scripts` (PR
   #64), because the process boundary scopes each stage's `code_sha256` to its own code.
   Nothing was archived: every script with no other caller is the producer of a committed
   number.

## 4. Definition of done

Judged on the report a user sees, not on stage exit codes:

- On a fresh machine, `conda env create -f environment.yml && pip install .` is the whole
  install.
- `hapax fetch-db --all` then `hapax run --config run.toml R1.fq.gz R2.fq.gz --out DIR` on an
  input the repo has never seen produces a report where every section is filled from a file
  the run wrote, with no file supplied or edited by hand.
- The same command on the CI fixture runs on every pull request, and the job fails if any
  section is missing, empty, or rendered from a file that does not match its schema.
- Re-running the Berlin and Yeongdok inputs reproduces their committed headline numbers,
  checked by a regression test against stored values.
