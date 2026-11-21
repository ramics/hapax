# Profile reach, and the pair as one molecule

Generated from `cells.json` and `paired_merged.json` by `experiments/report_profile_reach.py`. 5 seeds (1200 simulated pairs per seed per rung, 1200 real pairs per seed). Every cell is mean over seeds with the min–max range.

Operating points, stated because thresholds dominate everything else here: DIAMOND `--ultra-sensitive` against viral.1.protein.faa.gz + uniprot_sprot.fasta.gz, reported at a bitscore floor; Pfam-A searched with no score cutoff and reported at per-domain E-value floors, with the families' own curated gathering thresholds (`GA`) as one point on that curve.

## 1. Reach along the divergence ladder (mate 1, 107 nt)

| pile | DIAMOND bs>=50 | DIAMOND bs>=40 | Pfam E<=1e-2 | Pfam E<=1e-5 | Pfam at GA |
|---|---|---|---|---|---|
| holdout, 0% divergence | 91.5 (90.6–91.9) | 92.1 (91.0–92.8) | 39.9 (38.2–41.1) | 26.6 (24.1–28.1) | 37.7 (36.4–38.6) |
| holdout, 5% divergence | 81.8 (78.7–83.8) | 88.2 (86.8–89.8) | 27.6 (25.5–29.1) | 14.4 (13.8–15.2) | 25.0 (22.8–26.6) |
| holdout, 10% divergence | 61.3 (57.2–64.3) | 78.9 (76.6–81.9) | 18.6 (16.6–21.2) | 7.6 (6.2–10.0) | 16.7 (14.8–18.4) |
| holdout, 20% divergence | 10.0 (9.1–10.6) | 39.9 (35.9–46.7) | 5.8 (4.3–7.0) | 1.4 (0.1–3.2) | 4.6 (2.7–6.7) |
| holdout, 30% divergence | 0.1 (0.0–0.2) | 5.7 (2.1–8.9) | 0.5 (0.0–1.1) | 0.0 (0.0–0.0) | 1.1 (0.0–3.1) |
| real library reads (known material) | 11.1 (10.2–12.2) | 18.0 (16.8–20.6) | 5.9 (5.4–6.6) | 2.9 (2.6–3.2) | 4.6 (4.2–5.0) |
| composition-matched shuffle (null) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.4 (0.1–0.7) |

## 2. The gap: no translated hit, yes profile hit

The reportable quantity, because a profile hit rate on a coding background is a familiarity signal. Read the background and shuffle rows FIRST: if the gap reads the same on known material, it describes the database rather than the query.

| pile | gap (bs<50, E<=1e-2) | gap (bs<50, E<=1e-5) | gap (bs<40, E<=1e-2) | gap (bs<50, GA) |
|---|---|---|---|---|
| holdout, 0% divergence | 0.6 (0.4–0.8) | 0.3 (0.1–0.4) | 0.2 (0.2–0.3) | 0.6 (0.4–0.8) |
| holdout, 5% divergence | 1.3 (0.5–1.8) | 0.4 (0.1–0.9) | 0.3 (0.2–0.6) | 1.1 (0.5–1.7) |
| holdout, 10% divergence | 1.9 (0.8–2.2) | 0.8 (0.2–1.4) | 0.3 (0.0–0.6) | 1.7 (0.8–2.2) |
| holdout, 20% divergence | 3.8 (2.2–4.8) | 0.9 (0.0–2.2) | 1.4 (1.1–2.0) | 3.2 (1.8–4.6) |
| holdout, 30% divergence | 0.5 (0.0–1.1) | 0.0 (0.0–0.0) | 0.3 (0.0–1.0) | 1.1 (0.0–3.1) |
| real library reads (known material) | 2.6 (2.2–3.1) | 1.1 (0.7–1.4) | 1.1 (0.8–1.8) | 2.1 (1.3–2.7) |
| composition-matched shuffle (null) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.4 (0.1–0.7) |

## 3. Per mate, never pooled

Mate 2 is 93 nt against mate 1's 107 after the non-template prefix trim, so a pooled figure averages two instruments.

| pile | mate | mean query nt | peptide segments | DIAMOND bs>=50 | Pfam E<=1e-2 |
|---|---|---|---|---|---|
| holdout 0% | m1 | 106.18 (106.05–106.29) | 7.36 (7.31–7.38) | 91.5 (90.6–91.9) | 39.9 (38.2–41.1) |
| holdout 0% | m2 | 92.74 (92.64–92.78) | 6.61 (6.56–6.67) | 91.3 (90.5–92.0) | 33.4 (31.8–35.9) |
| holdout 5% | m1 | 106.16 (105.98–106.27) | 7.45 (7.36–7.50) | 81.8 (78.7–83.8) | 27.6 (25.5–29.1) |
| holdout 5% | m2 | 92.75 (92.69–92.77) | 6.71 (6.65–6.77) | 75.2 (72.0–78.5) | 21.9 (20.8–23.0) |
| holdout 10% | m1 | 106.15 (106.04–106.23) | 7.59 (7.55–7.64) | 61.3 (57.2–64.3) | 18.6 (16.6–21.2) |
| holdout 10% | m2 | 92.72 (92.64–92.79) | 6.81 (6.76–6.87) | 44.1 (40.7–49.8) | 14.0 (11.4–17.1) |
| holdout 20% | m1 | 106.14 (105.89–106.33) | 7.77 (7.71–7.83) | 10.0 (9.1–10.6) | 5.8 (4.3–7.0) |
| holdout 20% | m2 | 92.71 (92.60–92.77) | 6.96 (6.92–7.01) | 3.0 (1.8–4.4) | 3.8 (2.8–5.3) |
| holdout 30% | m1 | 106.24 (106.17–106.28) | 7.84 (7.59–7.98) | 0.1 (0.0–0.2) | 0.5 (0.0–1.1) |
| holdout 30% | m2 | 92.76 (92.74–92.78) | 7.03 (6.88–7.12) | 0.0 (0.0–0.0) | 0.2 (0.1–0.7) |
| real library | m1 | 106.13 (105.97–106.38) | 7.73 (7.68–7.79) | 11.1 (10.2–12.2) | 5.9 (5.4–6.6) |
| real library | m2 | 92.11 (91.95–92.38) | 7.01 (6.93–7.06) | 9.4 (7.8–10.7) | 4.7 (3.7–5.0) |

## 4. The pair as one molecule: merged against mate 1 of the SAME pair

Restricted to pairs whose mates overlap, so the comparison is within a molecule rather than between insert lengths. The merged read carries no error-correction credit for its overlap (stated bias 2), so these are lower bounds.

| pile | mean merged nt | mean mate-1 nt | DIAMOND bs>=50 merged | same, mate 1 | Pfam E<=1e-2 merged | same, mate 1 |
|---|---|---|---|---|---|---|
| holdout, 0% divergence | 139.71 (138.50–142.00) | 105.05 (104.83–105.27) | 91.1 (89.8–92.0) | 90.9 (89.6–92.0) | 50.8 (49.7–54.1) | 39.3 (37.1–41.5) |
| holdout, 5% divergence | 139.18 (138.39–140.83) | 105.06 (104.56–105.36) | 87.2 (85.5–88.9) | 80.1 (79.2–81.5) | 37.2 (35.5–39.0) | 26.2 (23.9–29.2) |
| holdout, 10% divergence | 139.99 (138.49–141.91) | 105.02 (104.79–105.16) | 75.5 (72.0–78.9) | 59.6 (55.7–64.3) | 27.9 (24.3–33.0) | 18.0 (15.7–22.5) |
| holdout, 20% divergence | 140.06 (138.38–141.18) | 104.98 (104.48–105.39) | 38.4 (33.2–43.5) | 10.5 (8.9–11.6) | 10.5 (7.9–14.9) | 5.8 (4.9–6.4) |
| holdout, 30% divergence | 139.84 (137.87–141.55) | 105.17 (105.07–105.30) | 7.1 (5.1–9.2) | 0.0 (0.0–0.2) | 2.2 (1.2–3.4) | 0.5 (0.0–1.2) |
| real library reads (known material) | 148.95 (146.84–150.73) | 106.56 (106.54–106.60) | 16.0 (14.1–17.7) | 8.9 (7.2–10.6) | 8.6 (7.9–9.4) | 5.5 (4.5–6.8) |
| composition-matched shuffle (null) | 148.95 (146.84–150.73) | 106.56 (106.54–106.60) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) | 0.0 (0.0–0.0) |

## Provenance

- generator `experiments/profile_reach.py`, content hash `c5ee9180fa12d098`
- this report `experiments/report_profile_reach.py`, content hash `851e014e2e3e737e`
- git SHA `6a941d66fc07d453aff4654fb5636c1af4004da7-dirty`
- holdout: D00239 (human rhinovirus A1B), ENA fasta, 7133 nt
- library: SRR23973193 (ENA fastq), mate 2 trimmed 14 nt
- timings: DIAMOND 530s, Pfam 3340s, total 3879s

