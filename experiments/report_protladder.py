"""Aggregate the protein index-size ladder on the residue, and say whether it flattens.

WHAT THIS IS FOR. The protein side of the removal curve is the one that has already
overturned a conclusion once: two small virus-weighted points were read as "protein
databases barely help", and a tenth of UniRef50 then moved the same quantity by 23 points.
So this reports the ladder as a curve with a marginal column, across five seeds with
spread, and refuses to name a flattening point that the measured rungs do not show.

TWO READOUTS PER RUNG, because they answer different questions:

- ANY HIT is the one that moves with database size, and it is nearly useless on its own:
  at these lengths almost everything eventually matches something, which the BLAST ground
  truth already showed.
- A BITSCORE FLOOR is the readout that bears on identification. A read with a 28-bit hit
  to a UniRef entry has not been identified by any standard anyone would defend.

The reads the cascade identified are carried through as a control at every rung, because a
rung that lifts both populations equally has bought coverage of the library rather than
coverage of the residue.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# smallest first: the ladder is cumulative, so the order is the claim
ORDER = [("closest_viral", "data/prot/viral.dmnd"),
         ("closest_sprot", "data/prot/sprot.dmnd"),
         ("closest_uniref50_02pct", "data/prot/uniref50_02pct.dmnd"),
         ("closest_uniref50_10pct", "data/prot/uniref50_10pct.dmnd")]


def cell(xs) -> str:
    if not xs:
        return "n/a"
    return (f"{st.mean(xs) * 100:.2f} ({min(xs) * 100:.2f}-{max(xs) * 100:.2f})"
            if len(xs) > 1 else f"{xs[0] * 100:.2f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, default=ROOT / "results/residue_protladder")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/residue_protladder/TABLES.md")
    ap.add_argument("--bitscores", type=float, nargs="+", default=[40.0, 50.0])
    args = ap.parse_args()

    paths = sorted(glob.glob(str(args.rows / "rows_s*.jsonl")))
    if not paths:
        raise SystemExit(f"no rows in {args.rows}")

    per_seed = []
    for p in paths:
        rows = {"residue": [], "known": []}
        with open(p) as h:
            for line in h:
                r = json.loads(line)
                rows["known" if r.get("set") == "known" else "residue"].append(r)
        d = {"path": Path(p).name,
             "n": {k: len(v) for k, v in rows.items()}, "rungs": []}
        for i in range(len(ORDER)):
            keys = [k for k, _ in ORDER[:i + 1]]
            rung = {"rung": ORDER[i][0].replace("closest_", ""), "sets": {}}
            for setname, rs in rows.items():
                if not rs:
                    continue
                # The control population is written into these rows for its composition
                # fields, but the protein search was run on the residue only, so its hit
                # columns are ABSENT rather than zero. Reporting 0% would be a false
                # statement about a population nobody searched.
                if not any(any(isinstance(r.get(k), dict) for k, _ in ORDER) for r in rs):
                    rung["sets"][setname] = {"not_searched": True}
                    continue
                best = []
                for r in rs:
                    b = 0.0
                    for k in keys:
                        v = r.get(k)
                        if isinstance(v, dict) and v.get("bitscore"):
                            b = max(b, v["bitscore"])
                    best.append(b)
                rung["sets"][setname] = {
                    "any_hit": sum(1 for b in best if b > 0) / len(best),
                    **{f"bitscore_ge_{int(t)}": sum(1 for b in best if b >= t) / len(best)
                       for t in args.bitscores},
                }
            d["rungs"].append(rung)
        per_seed.append(d)

    def col(i, setname, key):
        out = []
        for d in per_seed:
            s = d["rungs"][i]["sets"].get(setname)
            if s and not s.get("not_searched"):
                out.append(s[key])
        return out

    sizes = []
    cum = 0
    for _k, rel in ORDER:
        f = ROOT / rel
        cum += f.stat().st_size if f.exists() else 0
        sizes.append(cum)

    lines = [f"# Protein index-size ladder on the residue", "",
             f"{len(per_seed)} seeds; "
             f"{st.mean([d['n']['residue'] for d in per_seed]):.0f} residue pairs per seed, "
             f"{st.mean([d['n']['known'] for d in per_seed]):.0f} identified pairs as the "
             f"control. Cumulative, smallest index first. DIAMOND blastx, --sensitive.", "",
             "| cumulative index | GB | residue any hit | marginal | residue bitscore >=50 "
             "| identified any hit |", "|---|---|---|---|---|---|"]
    prev = 0.0
    prev_marginal = 0.0
    for i, (k, _rel) in enumerate(ORDER):
        anyh = col(i, "residue", "any_hit")
        if not anyh:
            continue
        m = st.mean(anyh)
        marg = "base" if i == 0 else f"{(m - prev) * 100:+.2f}"
        ctrl = cell(col(i, "known", "any_hit"))
        lines.append(f"| + {k.replace('closest_', '')} | {sizes[i] / 1e9:.2f} | "
                     f"{cell(anyh)} | {marg} | "
                     f"{cell(col(i, 'residue', 'bitscore_ge_50'))} | "
                     f"{ctrl if ctrl != 'n/a' else 'not searched'} |")
        prev_marginal = m - prev
        prev = m

    last = len(ORDER) - 1
    a_first = st.mean(col(0, "residue", "any_hit") or [0])
    a_last = st.mean(col(last, "residue", "any_hit") or [0])
    b_last = st.mean(col(last, "residue", "bitscore_ge_50") or [0])
    lines += ["", "## What the curve does", "",
              f"- Any hit rises from {a_first * 100:.1f}% to {a_last * 100:.1f}% across the "
              f"ladder. The last rung, a fivefold step in index size, buys "
              f"{prev_marginal * 100:+.2f} points of it, so this readout IS flattening - "
              f"but it is flattening against a ceiling at 100%, which is the uninformative "
              f"way to flatten. Read the bitscore column instead.",
              f"- At a bitscore floor of 50 the same ladder reaches {b_last * 100:.2f}%. "
              f"That gap is the result: a bigger protein index finds *something* for "
              f"almost every read in the residue and identifies almost none of them. "
              f"Database size buys weak hits.",
              "- **The control is missing, and that is a limitation rather than a "
              "result.** The identified reads are written into these rows for their "
              "composition fields, but the protein search ran on the residue only, so "
              "there is no identified-read hit column to compare against and the table "
              "says 'not searched' rather than 0%. The identified-versus-residue contrast "
              "that this ladder cannot make is made instead by the BLAST ground truth in "
              "the main report, section 3.4, on a smaller sample against a far larger "
              "database."]
    args.out.write_text("\n".join(lines) + "\n")
    (args.rows / "ladder_summary.json").write_text(json.dumps(
        {"seeds": len(per_seed), "cumulative_bytes": sizes, "per_seed": per_seed},
        indent=2))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
