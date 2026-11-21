# The depth ladder, per pile

Generated from the committed per-level TSVs by `experiments/report_depth_piles.py`. 5 seeds, 4,000 query pairs per pile per seed, drawn at random across the whole library and excluded from every library at every level.

Piles are pair-level, from the residue thread's shared per-mate label table; the pair rule is in `experiments/depth_rarefaction.py` and `unplaced` reproduces `results/residue_census/` §2 to the pair (1,372,021 pairs, 15.71%).

## 1. Fully unwitnessed fraction, by pile and effort (mate 1)

| pile | 0.125x | 0.25x | 0.5x | 1x |
|---|---|---|---|---|
| ribosomal | 0.5 (0.4–0.5) | 0.3 (0.2–0.4) | 0.2 (0.1–0.2) | 0.1 (0.1–0.1) |
| identified (non-rRNA contig) | 4.8 (4.1–5.3) | 2.5 (2.2–2.7) | 1.0 (0.9–1.1) | 0.3 (0.2–0.4) |
| contig outside the three-way split | 9.0 (8.8–9.2) | 6.8 (6.5–7.3) | 4.8 (4.5–5.0) | 2.7 (2.3–3.0) |
| residue on an unidentifiable contig | 25.8 (25.0–26.5) | 15.2 (14.5–16.0) | 7.4 (7.0–8.1) | 2.5 (2.4–2.8) |
| **residue, unplaced by anything** | 75.6 (75.0–76.7) | 62.4 (61.8–63.0) | 45.2 (44.5–46.6) | 26.7 (25.4–28.2) |

## 2. Median depth, by pile and effort (mate 1)

| pile | 0.125x | 0.25x | 0.5x | 1x |
|---|---|---|---|---|
| ribosomal | 381 (368–406) | 759 (738–810) | 1519 (1483–1600) | 3041 (2958–3214) |
| identified (non-rRNA contig) | 13 (13–14) | 28 (26–31) | 56 (52–61) | 112 (105–120) |
| contig outside the three-way split | 411 (397–421) | 812 (802–822) | 1614 (1586–1633) | 3230 (3158–3303) |
| residue on an unidentifiable contig | 1 (1–1) | 2 (2–2) | 3 (3–3) | 6 (6–7) |
| **residue, unplaced by anything** | 0 (0–0) | 0 (0–0) | 0 (0–0) | 1 (1–1) |

## 3. The transition: what one doubling rescues, by pile

Of the reads fully unwitnessed at one level, the fraction with any witness one doubling later. The same named reads scored against both libraries, joined on the name; no fit and no estimator.

**Mate 1.**

| pile | 0.125x -> 0.25x | 0.25x -> 0.5x | 0.5x -> 1x | n unwitnessed at 0.5x |
|---|---|---|---|---|
| ribosomal | 37.8 (19.0–57.9) | 47.3 (33.3–52.9) | 35.5 (0.0–62.5) | 6 |
| identified (non-rRNA contig) | 48.3 (45.1–52.3) | 58.2 (53.7–67.6) | 68.0 (65.7–74.4) | 41 |
| contig outside the three-way split | 24.5 (20.5–27.7) | 30.2 (26.7–32.6) | 42.7 (39.2–47.8) | 190 |
| residue on an unidentifiable contig | 41.0 (38.2–42.6) | 51.6 (47.2–53.6) | 65.6 (59.7–68.2) | 294 |
| **residue, unplaced by anything** | 17.5 (16.8–18.0) | 27.5 (26.1–28.5) | 41.0 (39.5–42.8) | 1,797 |

**Mate 2.**

| pile | 0.125x -> 0.25x | 0.25x -> 0.5x | 0.5x -> 1x | n unwitnessed at 0.5x |
|---|---|---|---|---|
| ribosomal | 34.1 (28.9–41.9) | 39.3 (17.6–54.5) | 34.1 (22.2–50.0) | 14 |
| identified (non-rRNA contig) | 45.7 (43.1–46.8) | 54.8 (49.3–60.5) | 57.7 (55.2–60.6) | 66 |
| contig outside the three-way split | 23.0 (18.5–26.9) | 28.1 (25.7–30.9) | 40.9 (36.4–46.6) | 193 |
| residue on an unidentifiable contig | 37.6 (36.6–38.5) | 49.0 (45.8–51.7) | 61.2 (53.7–65.0) | 365 |
| **residue, unplaced by anything** | 16.4 (15.7–16.8) | 26.8 (25.8–27.7) | 40.5 (39.0–42.8) | 1,856 |

## 4. A witness is not a depth

Of the same fully-unwitnessed reads, the fraction reaching median depth 5 or more one doubling later — the depth where assembly and classification begin to work at all (`results/spike_abundance/`). Mate 1.

| pile | gained a witness | reached depth >= 5 |
|---|---|---|
| ribosomal | 35.5 (0.0–62.5) | 0.0 (0.0–0.0) |
| identified (non-rRNA contig) | 68.0 (65.7–74.4) | 0.4 (0.0–2.2) |
| contig outside the three-way split | 42.7 (39.2–47.8) | 2.2 (1.5–2.8) |
| residue on an unidentifiable contig | 65.6 (59.7–68.2) | 0.7 (0.0–1.5) |
| **residue, unplaced by anything** | 41.0 (39.5–42.8) | 0.3 (0.2–0.5) |

## 5. The number this run exists for

**The unplaced residue's own rescue rate over the last measured doubling is 41.0 (39.5–42.8)% (mate 1, 5 seeds).** `results/depth_rarefaction/` measured 44.0 (40.6–46.2) on a random block of the library, and flagged that a random block is not the residue. This replaces that figure for the pile the question was about; the library figure stands for the library.

The library-wide run is [`../depth_rarefaction/TABLES.md`](../depth_rarefaction/TABLES.md), same instrument, same levels, same seeds.

## Provenance

- generator `experiments/depth_rarefaction.py`, content hash `2163ca17d074b4fc`
- this report `experiments/report_depth_piles.py`, content hash `e90b35c10c0b5b49`
- git SHA `9c8736ce44d6b7f78382c4adac3d2a0c2949abd3-dirty`
- label table `/mnt/project-files/shared/resasm/library_label_classes.npz`
- k = 31, levels [0.125, 0.25, 0.5, 1.0], mate 2 trimmed 14 nt
- library sizes at each level: 0.125x = 1,088,906, 0.25x = 2,176,908, 0.5x = 4,354,839, 1x = 8,711,224
- total runtime 1253s

