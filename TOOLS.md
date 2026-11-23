# External tools

The binaries hapax drives. None is pip-installable, so none is in `requirements.txt`; the
Python side is there. **The test suite needs none of these**: `tests/` replaces Kraken2 and
DIAMOND with stubs written into a temp directory, because what it tests is whether a stage
refuses and what it records, not the numbers. A fresh container runs the suite green with
`pip install -r requirements.txt` alone.

A run needs them. `hapax.py run` checks every input, index and tool before it starts and
refuses with exit status 2 if one is missing, so a missing binary is found in seconds, not
at the stage that needs it.

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

## Recorded, not pinned

Versions the committed manifests record but no code enforces. A run with another version
writes that version into its own manifest; whether it changes a number is unmeasured.

| Tool | Recorded | Used by |
|---|---|---|
| minimap2 | 2.31-r1302 (3 manifests), 2.28-r1209 (1) | `pair_ladder.py`, `path_settings_absorption.py`. The cascade uses `mappy` 2.31 instead, the same release. |
| MEGAHIT | v1.2.9 | `residue_assembly.py`, `pair_ladder.py` |
| SPAdes / metaSPAdes | v4.3.0 | `pair_ladder.py` |
| PenguiN | 5.cf8933-MPI | `pair_ladder.py` |
| Haploflow | not recorded | `pair_ladder.py` |

`experiments/codon_assembler/` is the archived first direction and keeps its own
dependencies (it also drives `seqkit`, `samtools` and `genomad`). Its Rust crate builds
with `cargo test --release` from that directory and its Python tests run with
`python3 -m pytest python/tests` from there; neither is part of the root suite.
