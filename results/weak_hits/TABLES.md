# Is there headroom in weak protein hits?

Generated from `bands.json` and `conversion.json` by `experiments/report_weak_hits.py`. 5 seeds, 100,912 queries, DIAMOND `--ultra-sensitive --max-target-seqs 1` reported down to E <= 100 so that a weak-hit population exists at all, with E <= 0.001 as the call line. Mean over seeds with the min–max range.

## 1. The weak band reads the same on real sequence as on shuffles

Mate 1. 'Any hit' is at the reporting cutoff, which is why it is near 100% everywhere **including the null** — that column is the chance band, not evidence.

| pile | any hit (E<=100) | significant (E<=1e-3) | median E of hits | median bits | median %id |
|---|---|---|---|---|---|
| holdout, 0% divergence | 99.9 (99.8–100.0) | 92.2 (91.2–92.9) | 4.2e-17 (4.2e-17–4.2e-17) | 74.7 (74.7–74.7) | 100.0 (100.0–100.0) |
| holdout, 5% divergence | 98.4 (97.6–99.1) | 88.9 (87.7–90.2) | 2.3e-13 (6.2e-14–4.9e-13) | 64.5 (63.5–65.9) | 88.6 (88.6–88.6) |
| holdout, 10% divergence | 98.0 (96.5–98.6) | 81.8 (79.8–84.4) | 8.8e-10 (2.5e-10–2e-09) | 54.6 (53.5–55.8) | 79.4 (78.8–80.0) |
| holdout, 20% divergence | 95.4 (93.6–97.6) | 47.8 (42.8–54.0) | 0.0014 (0.00018–0.0031) | 37.5 (36.2–39.3) | 61.0 (60.6–61.5) |
| holdout, 30% divergence | 92.7 (90.7–94.2) | 10.0 (7.5–14.8) | 1.7 (1.3–2.4) | 28.6 (28.1–28.9) | 50.7 (50.0–52.0) |
| real library reads (known material) | 95.9 (95.2–96.4) | 20.0 (18.6–22.2) | 1.7 (1.7–1.8) | 28.4 (28.1–28.5) | 50.0 (50.0–50.0) |
| **composition-matched shuffle (null)** | 92.3 (91.2–93.3) | 0.0 (0.0–0.0) | 3.3 (3.3–3.3) | 27.7 (27.7–27.7) | 47.7 (47.6–47.8) |

The weak hits themselves, with the null beside them:

| pile | weak hits (1e-3 < E <= 100) | median bits of those | median %id of those |
|---|---|---|---|
| holdout, 0% divergence | 7.7 (7.0–8.8) | 27.3 (27.3–27.3) | 47.2 (46.2–49.0) |
| holdout, 5% divergence | 9.5 (8.7–10.5) | 27.7 (27.3–28.1) | 48.0 (44.6–50.0) |
| holdout, 10% divergence | 16.3 (14.2–17.1) | 28.3 (27.7–28.5) | 51.4 (47.8–54.2) |
| holdout, 20% divergence | 47.5 (41.8–50.8) | 29.0 (28.5–29.6) | 53.9 (52.6–55.6) |
| holdout, 30% divergence | 82.7 (75.8–85.2) | 28.1 (28.1–28.1) | 49.6 (48.1–50.0) |
| real library reads (known material) | 75.8 (73.4–77.8) | 27.7 (27.7–27.7) | 47.2 (46.4–47.8) |
| **composition-matched shuffle (null)** | 92.3 (91.2–93.3) | 27.7 (27.7–27.7) | 47.7 (47.6–47.8) |

## 2. Merging converts weak hits, and converts none of the noise

Among pairs whose mate 1 produced only a weak hit, the fraction whose merged read produces a significant one. The control column is pairs already significant on mate 1, which merging should keep. The last column is the length-arithmetic null: a merged read's bitscore divided by its mate's scaled by the length ratio, where 1.0 means a longer query simply accumulated more of the same score.

| pile | n weak on mate 1 | converted | kept (already significant) | observed/length-predicted bits |
|---|---|---|---|---|
| holdout, 0% divergence | 40 | 3.9 (0.0–6.2) | 100.0 (100.0–100.0) | 0.987 (0.983–0.990) |
| holdout, 5% divergence | 51 | 23.5 (15.9–40.9) | 100.0 (99.8–100.0) | 0.995 (0.990–0.998) |
| holdout, 10% divergence | 87 | 35.6 (27.5–45.0) | 99.9 (99.5–100.0) | 0.983 (0.978–0.989) |
| holdout, 20% divergence | 246 | 34.6 (29.9–40.5) | 99.0 (98.4–100.0) | 0.952 (0.943–0.959) |
| holdout, 30% divergence | 413 | 18.2 (16.2–21.5) | 98.7 (96.3–100.0) | 0.873 (0.849–0.894) |
| real library reads (known material) | 330 | 6.6 (5.4–8.0) | 99.7 (98.5–100.0) | 0.794 (0.783–0.803) |
| **composition-matched shuffle (null)** | 387 | 0.0 (0.0–0.0) | n/a | 0.745 (0.735–0.760) |

## 3. The merged gain at a fixed E-value, not a bitscore floor

`../profile_reach/` §2 measured the merged advantage at a bitscore floor, where a longer query has an arithmetic edge. An E-value penalises query length, so if the gain were arithmetic it would shrink here. Same overlapping pairs, merged against mate 1 of the same molecule.

| pile | E<=1e-3 merged / mate 1 | E<=1e-5 merged / mate 1 | E<=1e-10 merged / mate 1 |
|---|---|---|---|
| holdout, 0% divergence | 92.2 (90.8–93.2) / 91.8 (90.2–92.8) | 92.0 (90.8–93.0) / 91.7 (90.2–92.6) | 89.8 (88.8–91.4) / 89.5 (88.4–91.4) |
| holdout, 5% divergence | 91.3 (88.2–93.4) / 88.2 (85.7–90.2) | 90.1 (87.5–92.3) / 85.5 (83.4–87.6) | 81.5 (79.2–83.6) / 70.9 (68.7–72.5) |
| holdout, 10% divergence | 87.8 (85.1–91.4) / 80.8 (77.6–83.7) | 82.9 (80.3–85.7) / 73.3 (69.2–75.7) | 64.5 (61.4–67.8) / 38.4 (35.0–42.5) |
| holdout, 20% divergence | 64.5 (61.8–72.4) / 47.4 (40.8–54.9) | 52.7 (49.0–60.6) / 29.5 (25.5–37.6) | 23.2 (19.2–26.9) / 2.5 (1.6–3.8) |
| holdout, 30% divergence | 25.4 (21.7–30.0) / 10.1 (6.9–15.8) | 13.9 (11.2–16.1) / 2.6 (1.0–4.7) | 2.1 (0.8–3.3) / 0.0 (0.0–0.0) |
| real library reads (known material) | 22.1 (20.2–22.9) / 16.8 (14.8–18.1) | 19.1 (16.9–20.8) / 12.6 (11.3–13.1) | 12.8 (11.1–15.6) / 6.3 (4.8–9.0) |
| **composition-matched shuffle (null)** | 0.0 (0.0–0.0) / 0.0 (0.0–0.0) | 0.0 (0.0–0.0) / 0.0 (0.0–0.0) | 0.0 (0.0–0.0) / 0.0 (0.0–0.0) |

The advantage survives at every threshold and **grows** at the strictest, which is the opposite of what a length artefact does.

## Provenance

- generator `experiments/weak_hit_conversion.py`, content hash `36b494243c468ff7`
- this report `experiments/report_weak_hits.py`, content hash `a34122435bc3af10`
- git SHA `c712f233259462807fed7c4f627e481b525dcceb-dirty`
- the query set is byte-identical to results/profile_reach/'s; only the reporting E-value cutoff differs, so the two runs are the same instrument at two operating points
- DIAMOND time 0s (0 means the search was reused from an earlier identical run)

