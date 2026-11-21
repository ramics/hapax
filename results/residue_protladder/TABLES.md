# Protein index-size ladder on the residue

5 seeds; 8000 residue pairs per seed, 8000 identified pairs as the control. Cumulative, smallest index first. DIAMOND blastx, --sensitive.

| cumulative index | GB | residue any hit | marginal | residue bitscore >=50 | identified any hit |
|---|---|---|---|---|---|
| + viral | 0.24 | 77.70 (77.18-78.24) | base | 0.00 (0.00-0.00) | not searched |
| + sprot | 0.53 | 92.44 (92.07-92.79) | +14.73 | 0.00 (0.00-0.00) | not searched |
| + uniref50_02pct | 0.90 | 97.50 (97.12-98.10) | +5.06 | 0.36 (0.33-0.40) | not searched |
| + uniref50_10pct | 2.72 | 98.90 (98.75-99.16) | +1.40 | 1.18 (0.96-1.30) | not searched |

## What the curve does

- Any hit rises from 77.7% to 98.9% across the ladder. The last rung, a fivefold step in index size, buys +1.40 points of it, so this readout IS flattening - but it is flattening against a ceiling at 100%, which is the uninformative way to flatten. Read the bitscore column instead.
- At a bitscore floor of 50 the same ladder reaches 1.18%. That gap is the result: a bigger protein index finds *something* for almost every read in the residue and identifies almost none of them. Database size buys weak hits.
- **The control is missing, and that is a limitation rather than a result.** The identified reads are written into these rows for their composition fields, but the protein search ran on the residue only, so there is no identified-read hit column to compare against and the table says 'not searched' rather than 0%. The identified-versus-residue contrast that this ladder cannot make is made instead by the BLAST ground truth in the main report, section 3.4, on a smaller sample against a far larger database.
