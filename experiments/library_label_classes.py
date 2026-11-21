"""Derive a per-mate class from `library_labels.py`'s table, using the cascade's own contig split.

Classes, at the cascade's 40-matching-base operating point, in this order of precedence:
  1 rrna        SILVA >= 40, or best own contig (>= 40) is in sample_contigs_rrna.fa
  2 identified  best own contig (>= 40) is in sample_contigs_identified_norrna.fa
  3 unident     best own contig (>= 40) is in sample_contigs_unidentified.fa
  4 other_ctg   best own contig (>= 40) is in none of the three (7,944 contigs the split omits)
  0 none        nothing at >= 40
The contig split's own operating points (from the pipeline thread, experiments/classify_contigs.py):
rRNA contig = >= 100 matching bases on SILVA; identified = Kraken2 PlusPF-16 at confidence 0.1 OR
DIAMOND bitscore >= 50 against Swiss-Prot / UniRef50-10%.

STATED BIAS: no read-level Kraken2 or DIAMOND, so "none" over-counts against the cascade.
On the five seed blocks the table's rRNA fraction (37.96-38.13%) matches the cascade's rrna tier
(37.95-38.13%) and no residue pair reaches 40 on SILVA.
"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
from subtract_known import provenance

work = ROOT / "data/work"
sets = ROOT / "data/geo/sets"
d = np.load(work / "library_labels.npz")
names = d["contig_names"]
code = {}
for f, c in (("sample_contigs_unidentified.fa", 3), ("sample_contigs_identified_norrna.fa", 2),
             ("sample_contigs_rrna.fa", 1)):
    for line in open(sets / f):
        if line[0] == ">":
            code[line[1:].split()[0]] = c
ctg_class = np.array([code.get(n, 4) for n in names], dtype=np.int8)
out = {"spot": d["spot"]}
summary = {}
for m in (1, 2):
    ctg, cm, sv = d[f"ctg{m}"], d[f"cm{m}"], d[f"silva{m}"]
    cl = np.zeros(len(ctg), dtype=np.int8)
    ok = (cm >= 40) & (ctg >= 0)
    cl[ok] = ctg_class[ctg[ok]]
    cl[sv >= 40] = 1
    out[f"class{m}"] = cl
    summary[f"mate{m}"] = {k: float((cl == v).mean()) for k, v in
                           (("none", 0), ("rrna", 1), ("identified", 2), ("unident", 3), ("other_ctg", 4))}
np.savez_compressed(work / "library_label_classes.npz", **out,
                    class_codes=np.array(["none", "rrna", "identified", "unident", "other_ctg"]))
man = {"script": "experiments/library_label_classes.py", **provenance(__file__),
       "labels_manifest": json.loads((work / "library_labels_manifest.json").read_text()),
       "class_fractions_of_library_mates": summary}
(work / "library_label_classes_manifest.json").write_text(json.dumps(man, indent=2))
print(json.dumps(summary, indent=1))
