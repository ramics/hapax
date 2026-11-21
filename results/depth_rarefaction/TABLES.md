# What deeper sequencing of this library would buy

Generated from the committed per-level TSVs by `experiments/report_depth_rarefaction.py`. 5 seeds, each drawing its own query block of 10,000 real pairs AND its own nested subsample, so the spread covers both. Query pairs are excluded from every library at every level, so a depth is always occurrences elsewhere.

## 1. The curve

| effort | library pairs | mate | median depth | fully unwitnessed | depth < 5 |
|---|---|---|---|---|---|
| 0.125x | 1,090,240 | 1 | 22 (21–24) | 16.2 (16.0–16.6) | 38.5 (38.0–39.4) |
| 0.125x | 1,090,240 | 2 | 17 (15–19) | 17.3 (17.2–17.7) | 40.3 (39.8–41.1) |
| 0.25x | 2,179,512 | 1 | 46 (42–48) | 12.2 (12.0–12.6) | 34.0 (33.6–34.7) |
| 0.25x | 2,179,512 | 2 | 34 (30–37) | 13.3 (13.0–13.5) | 35.6 (35.1–36.3) |
| 0.5x | 4,359,953 | 1 | 90 (81–96) | 8.1 (8.0–8.2) | 29.4 (29.0–30.1) |
| 0.5x | 4,359,953 | 2 | 68 (62–75) | 8.9 (8.6–9.2) | 31.1 (30.6–31.7) |
| 1x | 8,721,224 | 1 | 180 (163–190) | 4.5 (4.4–4.8) | 24.4 (23.8–25.3) |
| 1x | 8,721,224 | 2 | 137 (122–146) | 5.0 (4.8–5.3) | 26.0 (25.6–26.4) |

'Fully unwitnessed' is every scorable 31-mer of the read having no witness anywhere else in the library — the strict version, not 'median depth zero'. A read with an internal zero run of about 31 in otherwise covered sequence is a substitution error, not a rare molecule, and the strict definition excludes it by construction.

## 2. The transition: what one doubling actually rescues

Of the reads fully unwitnessed at one level, the fraction with any witness one doubling later. The same named reads scored against both libraries; no fit and no estimator.

| doubling | mate | n unwitnessed before | gained a witness | reached depth >= 5 |
|---|---|---|---|---|
| 0.125x -> 0.25x | 1 | 1,616 | 24.4 (23.7–25.4) | 0.0 (0.0–0.1) |
| 0.125x -> 0.25x | 2 | 1,730 | 23.5 (22.3–24.9) | 0.0 (0.0–0.1) |
| 0.25x -> 0.5x | 1 | 1,222 | 34.0 (33.1–34.7) | 0.1 (0.1–0.2) |
| 0.25x -> 0.5x | 2 | 1,324 | 32.9 (30.8–34.4) | 0.1 (0.0–0.3) |
| 0.5x -> 1x | 1 | 807 | 44.0 (40.6–46.2) | 0.3 (0.1–0.5) |
| 0.5x -> 1x | 2 | 888 | 43.8 (41.1–45.8) | 0.4 (0.1–0.5) |

The rate RISES along the ladder, because each doubling adds twice the absolute number of reads of the one before it while witnessing depends on absolute depth. So the last measured transition — 44.0 (40.6–46.2)% of mate 1's fully unwitnessed reads gaining a witness when the library goes from half to full — is a FLOOR for what the next doubling would do, not a ceiling.

No 10x figure is reported. Interpolation inside the measured 8x range is measurement; one doubling is a short extrapolation anchored on three measured transitions; ten times is neither, and this repo has already withdrawn one claim for reading a curve past where it was measured.

## Stated biases, with direction

- **Subsampling a finished library models re-sequencing the SAME library deeper, not sequencing new water.** The same plant's rRNA fraction went 38.0% to 67.1% in eight months, so new water is a different population. Direction: optimistic about new water, roughly right about deeper sequencing of this library.
- **Duplicates read as depth.** A rescued read whose witness is a PCR duplicate of its own molecule has gained confirmation that it is not an error, but no new sequence. Direction: inflates the rescue rate as a measure of assembly gain, while remaining honest as a measure of 'seen more than once'.
- **A random query block is not the residue.** These are real library reads of every kind, so the curve describes the library. The unplaced residue sits at the bottom of this distribution by definition and its own transition rate will be lower than the library average; the ladder here bounds the shape, not that pile's value.

## Provenance

- generator `experiments/depth_rarefaction.py`, content hash `8d3ceb328e3eddd3`
- this report `experiments/report_depth_rarefaction.py`, content hash `a1dae993272f3905`
- git SHA `6a941d66fc07d453aff4654fb5636c1af4004da7-dirty`
- k = 31, levels [0.125, 0.25, 0.5, 1.0], 10,000 query pairs per seed, mate 2 trimmed 14 nt
- total runtime 1325s; it shared the machine with another job, so the timings in `manifest.json` are operational, not benchmarks

