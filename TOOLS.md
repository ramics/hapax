# External tools

The binaries hapax drives. None is pip-installable, so none is in `requirements.txt`; the
Python side is there. **The test suite needs none of these**: `tests/` replaces Kraken2 and
DIAMOND with stubs written into a temp directory, because what it tests is whether a stage
refuses and what it records, not the numbers. A fresh container runs the suite green with
`pip install -r requirements.txt` alone.

A run needs them. `hapax.py run` checks every input, index and tool before it starts and
refuses with exit status 2 if one is missing, so a missing binary is found in seconds, not
at the stage that needs it.

**Where they come from, since 2026-11-23.** Kraken2, MEGAHIT and minimap2 used to exist
only as bare paths into somebody's staging directory (`$S/k2pkg/root/share/...` appears in
nine scripts) with nothing in the repo that built them, so a production run on a new
machine could not be started by reading this repo. `experiments/prod_userdata.sh` stages
all three from pinned sources and proves each one runs before anything expensive starts.

## Pinned, and why

A version is **pinned** where a committed number depends on it or where an earlier version
has failed. Moving a pin is a provenance change: re-run what it describes.

| Tool | Version | Pinned in | Used by | Why this version |
|---|---|---|---|---|
| foldseek | **10-941cd33 or later** | `ec2_reach.py` `FOLDSEEK_VERSION` | `reach_arms.py` structure arm | **9-427df8a segfaults** loading the f16 ProstT5 gguf (it reads only `model.gguf` / `model.safetensors`); 10-941cd33 reads `prostt5-f16.gguf` and writes 3Di on CPU. Cost four EC2 runs. See the comment at the pin. |
| DIAMOND | 2.1.10 | `ec2_reach.py` `DIAMOND_VERSION` | the protein tier of `subtract_known.py` and 20 other scripts | Every timing in `results/diamond_scaling/` was measured with it, and it is the production search. Five older manifests record 2.1.9. |
| BLAST+ | 2.17.0 | `ec2_reach.py` `BLAST_VERSION` | `blastn`: `holdout_*`, `bigindex_search.py`; `blastp`: `reach_arms.py`, `divergent_contig_workup.py`; `makeblastdb`: `holdout_indexes.py` | Recorded in 14 manifests. |
| HMMER | 3.4 | `ec2_reach.py` `HMMER_VERSION` | `hmmsearch`, `jackhmmer` arms of `reach_arms.py` | Recorded in 9 manifests. Separate from `pyhmmer` in `requirements.txt`. |
| HH-suite | 3.3.0 | `ec2_reach.py` `HHSUITE_VERSION` | `hhsearch` arm of `reach_arms.py` | Recorded in 8 manifests; searched against `pfamA_35.0` (`HH_PFAM_URL`). |
| Kraken2 | 2.17.2 | default `--classify-bin` path in `subtract_known.py` and 8 other scripts | Kraken2 tiers of the cascade | The cascade calls the `classify` binary directly, not the `kraken2` wrapper, so `--classify-bin` points inside the package. Database is part of the operating point: manifests record `k2_viral_20210517.d`, `k2_pluspf08`, `k2_pluspf16`. |
| Kraken2 | 2.17.2 | `ec2_prod.py` `KRAKEN2_URL` | staged by `prod_userdata.sh` | The same version, now with a source: `conda.anaconda.org/bioconda/linux-64/kraken2-2.17.2-pl5321h3be2455_0.conda`, extracted under `k2pkg/root` so the path matches the one nine scripts already name. The `classify` binary links only against system libraries on Amazon Linux 2023, and the instance proves it runs before staging 33 GB of databases. |
| MEGAHIT | 1.2.9 | `ec2_prod.py` `MEGAHIT_URL` | `path_settings.py`, `residue_assembly.py`, `pair_ladder.py` | `bioconda/linux-64/megahit-1.2.9-haf24da9_8.tar.bz2`. Moved from "recorded, not pinned" to pinned when it became a committed download: every assembly in this repo records v1.2.9 and an assembly is not re-run cheaply. |
| minimap2 | 2.31 | `ec2_prod.py` `MINIMAP2_URL` | `path_settings.py` `--minimap2`, `truth_panel.py` | `bioconda/linux-64/minimap2-2.31-h118bc1c_0.conda`, the same release as `mappy` 2.31 in `requirements.txt`, so the binary and the library side of the cascade are one version. |
| Node.js | 22 | `environment.yml`, `.github/workflows/end_to_end.yml` | `experiments/make_web.py`'s render check | The page is one HTML file with its data inlined, and the check that it DREW is a browser loading it. A section can resolve from its file and render nothing, which no file test can see. |
| Playwright (node) | 1.49.1 | `.github/workflows/end_to_end.yml` | the same render check | Installed with its chromium in the end-to-end job. Without the module the check returns nothing and the build records that it did not run, so the job would pass while asserting less than it appears to. Separate from `playwright` in `requirements.txt`, which nothing in the report path uses. |

## Recorded, not pinned

Versions the committed manifests record but no code enforces. A run with another version
writes that version into its own manifest; whether it changes a number is unmeasured.

| Tool | Recorded | Used by |
|---|---|---|
| minimap2 | 2.31-r1302 (3 manifests), 2.28-r1209 (1) | `pair_ladder.py`, `path_settings_absorption.py`. The cascade uses `mappy` 2.31 instead, the same release. Pinned above since `ec2_prod.py` stages it. |
| MEGAHIT | v1.2.9 | `residue_assembly.py`, `pair_ladder.py`. Pinned above since `ec2_prod.py` stages it. |
| SPAdes / metaSPAdes | v4.3.0 | `pair_ladder.py` |
| PenguiN | 5.cf8933-MPI | `pair_ladder.py` |
| Haploflow | not recorded | `pair_ladder.py` |
| Infernal (`cmscan`) | not recorded | `residue_report.py`'s `rfam` field, and nothing else. **It is in neither `environment.yml` nor any staging script, so the documented install cannot run that field** — found when a production run's first seed finished with `--report-without rfam` and a tripwire asked why nothing passed `--cmscan-bin`. The field is therefore reported as not searched by name rather than as a zero, and `hapax run` resolves the binary only when the field is requested, so asking for it on a machine without Infernal is a refusal before the run starts instead of a container path that exists on one machine. Adding it to the install means pinning a version and re-running what the `rfam` field appears in. |

`experiments/codon_assembler/` is the archived first direction and keeps its own
dependencies (it also drives `seqkit`, `samtools` and `genomad`). Its Rust crate builds
with `cargo test --release` from that directory and its Python tests run with
`python3 -m pytest python/tests` from there; neither is part of the root suite.
