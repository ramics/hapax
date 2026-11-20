"""Fill the placeholders in results/peptide_absence_noise/REPORT.md from the JSON.

Kept as a script rather than done by hand so that the prose in the report is generated
from the committed outputs, not transcribed from them. The project has lost claims to
report prose that drifted from the numbers beside it.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "peptide_absence_noise"
d = json.loads((RES / "peptide_absence_noise.json").read_text())
rep = (RES / "REPORT.md").read_text()

agg = d["aggregate"]
s0 = d["per_seed"][0]["synonymous"]["synonymous_fraction_by_codon_position"]


def pct(x):
    return f"{100 * x:.1f}%"


def spread(sp):
    return f"{100 * sp['mean']:.2f}% ({100 * sp['min']:.2f}–{100 * sp['max']:.2f})"


rep = rep.replace("PLACEHOLDER_SYN", f"**{pct(agg['synonymous_fraction']['mean'])}** "
                                     f"({pct(agg['synonymous_fraction']['min'])}–"
                                     f"{pct(agg['synonymous_fraction']['max'])} across seeds)")
rep = rep.replace("PLACEHOLDER_P0", pct(s0["0"]))
rep = rep.replace("PLACEHOLDER_P1", pct(s0["1"]))
rep = rep.replace("PLACEHOLDER_P2", pct(s0["2"]))

rows = ["| substitution rate | peptide corrupted | nucleotide corrupted | difference | reads with none corrupted |",
        "|---|---|---|---|---|"]
for r in agg["corruption"]:
    if r["profile"] != "flat":
        continue
    rows.append(f"| {r['error_rate']} | {spread(r['peptide_corrupted_fraction'])} | "
                f"{spread(r['nucleotide_corrupted_fraction'])} | "
                f"{100 * r['peptide_minus_nucleotide']['mean']:+.2f} pp | "
                f"{100 * r['reads_with_zero_corrupted']['mean']:.1f}% |")
rep = rep.replace("PLACEHOLDER_CORRUPTION_TABLE", "\n".join(rows))

rows = ["| indel rate | peptide corrupted | nucleotide corrupted |", "|---|---|---|"]
for r in agg["indels"]:
    rows.append(f"| {r['indel_rate']} | {spread(r['peptide_corrupted_fraction'])} | "
                f"{spread(r['nucleotide_corrupted_fraction'])} |")
rep = rep.replace("PLACEHOLDER_INDEL_TABLE", "\n".join(rows))

# headroom: aggregate the tail across seeds, reporting the worst (highest) threshold seen
budgets = [b for b in d["per_seed"][0]["tail"]["rows"][0]["thresholds"]]
rows = ["| substitution rate | " + " | ".join(f"≤{float(b):,.0f}/M" for b in budgets) + " |",
        "|---" * (len(budgets) + 1) + "|"]
for i, ref in enumerate(d["per_seed"][0]["tail"]["rows"]):
    cells = []
    for b in budgets:
        ts, hs, res = [], [], []
        for p in d["per_seed"]:
            v = p["tail"]["rows"][i]["thresholds"][b]
            if v["min_threshold"] is not None:
                ts.append(v["min_threshold"])
                hs.append(v["headroom_kmers"])
            res.append(v["resolved"])
        if not ts:
            cells.append("not reachable")
        else:
            t = max(ts)
            h = min(hs)
            flag = "" if all(res) else " \\*"
            cells.append(f"T ≥ {t}, headroom {h}{flag}" if h >= 0 else f"not reachable{flag}")
    rows.append(f"| {ref['error_rate']} | " + " | ".join(cells) + " |")
rows.append("")
rows.append(f"Worst threshold and smallest headroom across {len(d['per_seed'])} seeds, "
            f"{d['per_seed'][0]['tail']['rows'][0]['n_samples']:,} samples per seed "
            f"(resolution {d['per_seed'][0]['tail']['rows'][0]['resolution_per_million']:.2f} per million). "
            "An entry marked \\* hit a zero count in at least one seed, so it is a bound rather "
            "than an estimate.")
rep = rep.replace("PLACEHOLDER_HEADROOM_TABLE", "\n".join(rows))

if "PLACEHOLDER" in rep:
    sys.exit("placeholders left unfilled: " +
             ", ".join(w for w in rep.split() if "PLACEHOLDER" in w))
(RES / "REPORT.md").write_text(rep)
print("filled")
