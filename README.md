# hapax

Codon-aware assembly for viral populations. Overlaps scored in codon space so
synonymous variation collapses; nucleotides kept throughout so linkage survives.

See `METHOD.md` for the spec, `CLAUDE.md` for how to work in the repo, and
`results/README.md` for what the first experiments found (short version: the
codon gain is real but small, and only shows at high diversity).

## Build and run

```
cargo build --release
./target/release/hapax assemble -i reads.fa -o out/            # codon mode
./target/release/hapax assemble -i reads.fa -o out/ --seed spaced
./target/release/hapax assemble -i reads.fa -o out/ --mode nucleotide   # ablation
./target/release/hapax assemble -i reads.fa -o out/ --split-phase 1    # experimental, see below
```

All options are listed at the top of `src/main.rs`.

Outputs in `out/`:

| file | contents |
|---|---|
| `contigs.fa` | one consensus per read group (target: one per population) |
| `variants.tsv` | per-site allele counts, frequency, region type, syn/nonsyn effect per covering ORF frame |
| `linkage.tsv` | within-read haplotype counts and r² between consecutive variant sites |
| `regions.tsv` | consensus annotated noncoding / coding / dual-coding (from ORFs in six frames) |
| `layout.tsv` | which read went where, on which strand |
| `stats.json` | dispatch counters (incl. `hmm_fraction`, the METHOD.md §5 speed number) and stage timings |

## How it works

| stage | file | what and why |
|---|---|---|
| seeding | `src/index.rs` | (w,k) minimizers, canonical across strands; optional codon-spaced seed (mask `110`) that ignores wobble positions in whichever phase lines up with the frame |
| candidates | `src/overlap.rs` | anchors grouped by read pair and diagonal; top partners per read kept |
| tier 1 | `src/overlap.rs` | ungapped nucleotide check on the seed diagonal: clearly-same pairs accepted, clearly-different rejected |
| tier 2 | `src/hmm.rs` | RAMICS-style pair HMM with read A as the codon grid, run at the three relative offsets; in-codon frameshift-error states, codon indels, and a nucleotide-mode state so scoring mode is chosen per region; dual-coding signal falls back to nucleotide scoring |
| scoring model | `src/model.rs` | GY94 codon log-odds (ω = 0.23, κ = 2) and K80 nucleotide log-odds on the same background, both with sequencing error convolved in |
| population splitting | `src/phase.rs`, `layout::split_sparse` | **experimental and off by default.** Two swarms of one virus from different hosts share nucleotide-identical overlaps through conserved windows, so no per-pair score can keep them apart, and hapax merges them. `--split-phase` fits two populations to a contig's variant sites; `--split-cut` looks for a thin cut in the overlap graph. Both fail on real sister swarms, in the ways `results/swarm_metagenome/README.md` records |
| layout | `src/layout.rs` | best-k clusters merged by bridge counts (stops haplotype splitting; it does not stop two hosts' swarms merging, see `results/swarm_metagenome/README.md`), then a maximum-spanning-tree embedding |
| consensus | `src/consensus.rs` | iterative realignment to the draft, pileup, variants, linkage, ORF/region annotation |

## Python

```
python -m pytest python/tests            # simulator + surrogate tests
python python/swarm.py --seed 1 --diversity 0.1 --out /tmp/world
python python/evaluate.py /tmp/world out/             # membership-based scoring
python python/evaluate.py /tmp/world out/ --kmer      # contig-only scoring (any assembler)
python experiments/diversity_ladder.py --seed 1 --n-seeds 5 --out results/ladder
python experiments/report_ladder.py results/ladder
python experiments/scaling.py --seed 1 --out results/scaling
python experiments/fetch_swarm_panels.py                  # real within-host panels
python experiments/swarm_metagenome.py --seed 1 --n-seeds 5 --out results/swarm_metagenome
python experiments/report_swarm_metagenome.py results/swarm_metagenome
python experiments/sister_pairs.py --seed 1 --n-seeds 3 --out results/sister_pairs
```
