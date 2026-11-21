"""How many measurements are there really? Correlation structure of the per-read fields.

WHY THIS IS A RESULT AND NOT A DIAGNOSTIC. The deliverable Imo asked for is a
multidimensional measurement of signal against noise, so the DIMENSIONALITY is part of the
claim. A dozen columns that all track the same underlying thing is one measurement wearing
a hat, and reporting it as twelve would overstate what the report can tell you about a
read. This script measures how much of the variation is shared.

Three readouts, in increasing order of how much they assume:

1. Spearman rank correlation between every pair of fields, on the residue reads. Rank
   rather than Pearson because several fields are bounded, skewed or zero-inflated (fold
   z-scores, duplicate fractions, bitscores with a floor at no-hit).
2. Groups of fields at |rho| >= 0.5, which is the practical statement: within a group,
   knowing one field tells you most of another.
3. The eigenvalue spectrum of the rank-correlation matrix, and its participation ratio
   (exp of the entropy of the normalised eigenvalues). That number is the honest answer to
   "how many independent things are you measuring": 1 would mean a single axis wearing
   twelve hats, and the field count itself would mean no redundancy at all.

Every field is also compared between residue reads and the reads the cascade identified,
because a field that does not separate those two populations contributes nothing to a
signal-against-noise statement no matter how independent it is. Independence and
usefulness are different properties and the table reports both.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (column name, how to pull it out of a row). Kept explicit rather than walking the JSON,
# because a field that means nothing at read length should not quietly join the matrix.
FIELDS = {
    "gc": lambda r: r.get("gc"),
    "entropy3": lambda r: r.get("entropy3"),
    "homopolymer": lambda r: r.get("homopolymer"),
    "dup_kmer16": lambda r: _g(r, "repetitive", "kmer16_duplicate_frac"),
    "longest_orf_aa": lambda r: _g(r, "frame", "longest_orf_aa"),
    "stop_free_frames": lambda r: _g(r, "frame", "stop_free_frames"),
    "CpG_oe": lambda r: _g(r, "dinuc", "CpG_oe"),
    "UpA_oe": lambda r: _g(r, "dinuc", "UpA_oe"),
    "fold_z": lambda r: _g(r, "fold", "z"),
    "fold_mfe": lambda r: _g(r, "fold", "mfe"),
    "library_abundance": lambda r: _g(r, "library_abundance", "median_31mer_count"),
    "cross_sample_samples": lambda r: _g(r, "cross_sample", "distinct_samples"),
    "closest_bitscore": None,  # filled below, once _best_bits is defined
    "mean_error_prob": lambda r: r.get("mean_error_prob"),
    "frac_bases_under_q20": lambda r: r.get("frac_bases_under_q20"),
}


def _g(row, *path):
    cur = row
    for p in path:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(p)
    return cur


def _best_bits(row):
    best = None
    for k, v in row.items():
        if k.startswith("closest_") and isinstance(v, dict):
            b = v.get("bitscore")
            if b is not None and (best is None or b > best):
                best = b
    return best if best is not None else 0.0


FIELDS["closest_bitscore"] = _best_bits


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    return num / (da * db) if da and db else 0.0


def load_rows(paths):
    rows = []
    for p in paths:
        with open(p) as h:
            for line in h:
                rows.append(json.loads(line))
    return rows


def matrix(rows, names):
    cols = {}
    for n in names:
        vals = [FIELDS[n](r) for r in rows]
        if sum(1 for v in vals if isinstance(v, (int, float))) < 0.5 * len(rows):
            continue
        cols[n] = [float(v) if isinstance(v, (int, float)) else 0.0 for v in vals]
    return cols


def eigenvalues(m):
    """Symmetric eigenvalues by Jacobi rotation: the matrix is at most ~15x15 here, so a
    dependency on numpy's linear algebra is not worth adding for it."""
    n = len(m)
    a = [row[:] for row in m]
    for _sweep in range(100):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-12:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-15:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return sorted((a[i][i] for i in range(n)), reverse=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", type=Path, default=ROOT / "results/residue_fields")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/residue_fields/INDEPENDENCE.md")
    ap.add_argument("--json-out", type=Path,
                    default=ROOT / "results/residue_fields/independence.json")
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    paths = sorted(glob.glob(str(args.fields / "rows_s*.jsonl")))
    if not paths:
        raise SystemExit(f"no per-read rows in {args.fields}")
    rows = load_rows(paths)
    residue = [r for r in rows if r.get("set", "residue") != "known"]
    known = [r for r in rows if r.get("set") == "known"]

    cols = matrix(residue, list(FIELDS))
    names = list(cols)
    ranks = {n: rank(cols[n]) for n in names}
    rho = {a: {b: (1.0 if a == b else pearson(ranks[a], ranks[b])) for b in names}
           for a in names}

    # groups at the threshold, by single linkage
    groups, seen = [], set()
    for a in names:
        if a in seen:
            continue
        grp = {a}
        changed = True
        while changed:
            changed = False
            for b in names:
                if b in grp:
                    continue
                if any(abs(rho[x][b]) >= args.threshold for x in grp):
                    grp.add(b)
                    changed = True
        seen |= grp
        groups.append(sorted(grp))

    m = [[rho[a][b] for b in names] for a in names]
    ev = [max(v, 0.0) for v in eigenvalues(m)]
    tot = sum(ev) or 1.0
    p = [v / tot for v in ev]
    entropy = -sum(x * math.log(x) for x in p if x > 0)
    participation = math.exp(entropy)

    # separation from the known control, where both populations exist
    separation = {}
    if known:
        kcols = matrix(known, names)
        for n in names:
            if n not in kcols:
                continue
            a, b = cols[n], kcols[n]
            # A FIELD THAT IS IDENTICALLY ZERO IN THE CONTROL WAS NOT MEASURED THERE, and
            # saying so is not pedantry: `closest_bitscore` read 0 for every identified
            # read and produced a standardised difference of +9.5, the largest in the
            # table, purely because the protein search was run on the residue only. The
            # protein ladder had to re-label a whole column for the same reason. A constant
            # zero beside a varying column is absence, not a measurement of zero, so it is
            # named rather than scored.
            if st.pstdev(b) == 0 and st.mean(b) == 0 and st.pstdev(a) > 0:
                separation[n] = {"residue_mean": st.mean(a), "known_mean": None,
                                 "standardised_difference": None,
                                 "not_measured_in_control": True}
                continue
            pooled = (st.pstdev(a) + st.pstdev(b)) / 2 or 1.0
            separation[n] = {"residue_mean": st.mean(a), "known_mean": st.mean(b),
                             "standardised_difference": (st.mean(a) - st.mean(b)) / pooled}

    out = {"reads": len(residue), "fields": names, "spearman": rho, "groups": groups,
           "eigenvalues": ev, "participation_ratio": participation,
           "separation_from_known": separation, "threshold": args.threshold,
           "seeds": len(paths)}
    args.json_out.write_text(json.dumps(out, indent=2))

    lines = [f"# How many measurements are in the per-read report?", "",
             f"{len(residue):,} residue reads over {len(paths)} seeds, {len(names)} "
             f"numeric fields.", "",
             f"**Participation ratio {participation:.1f} of {len(names)} fields.** That "
             f"is the effective number of independent axes: the report is measuring "
             f"about {participation:.0f} things, not {len(names)}.", "",
             "## Groups of fields that move together (|rho| >= "
             f"{args.threshold})", ""]
    for g in groups:
        lines.append(f"- {', '.join(g)}" + ("" if len(g) > 1 else "  (independent)"))
    lines += ["", "## Does the field separate residue from the reads the cascade "
              "identified?", "",
              "| field | residue mean | known mean | standardised difference |",
              "|---|---|---|---|"]
    for n, s in sorted(separation.items(),
                       key=lambda kv: -abs(kv[1]["standardised_difference"] or 0)):
        if s.get("not_measured_in_control"):
            lines.append(f"| {n} | {s['residue_mean']:.4g} | *not searched* | "
                         f"*cannot be computed* |")
            continue
        lines.append(f"| {n} | {s['residue_mean']:.4g} | {s['known_mean']:.4g} | "
                     f"{s['standardised_difference']:+.2f} |")
    nm = [n for n, s in separation.items() if s.get("not_measured_in_control")]
    if nm:
        lines += ["", "**" + ", ".join(f"`{n}`" for n in nm) + "** read identically zero "
                  "across the control, which means the measurement was never made there - "
                  "the protein search ran on the residue only. Reported as *not searched* "
                  "rather than as a difference: scored naively, `closest_bitscore` showed "
                  "the largest separation in this table (+9.5 standardised) and it was an "
                  "artefact of a population nobody searched."]
    lines += ["", "Independence and usefulness are different properties: a field can be "
              "its own axis and still say nothing about whether a read is interesting, "
              "and both columns are needed to read the report honestly.", "",
              "Two results from the gate experiment (PR #17) belong beside this table. "
              "First, the composition and quality fields were tested JOINTLY on residue "
              "against subtracted reads and reach AUC 0.686, where the full space reaches "
              "0.905 and the coverage axes alone reach 0.904: the axes PR #8 killed one at "
              "a time also fail together, so a high participation ratio among them is "
              "dimensionality without discrimination. Second, the residue is a continuum - "
              "no cluster count beats k=1 against a covariance-matched reference over five "
              "seeds - so nothing here should be read as carving groups.", "",
              "A measured caveat on `library_abundance`: it comes from a count-min sketch "
              "that only ever overstates and that includes the read's own occurrences, so "
              "it cannot resolve coverage SHAPE. PR #17's exact 31-mer counter "
              "(`experiments/kmer_coverage/`, whole library, own occurrences subtracted) "
              "does, at 0.76-0.78 AUC per shape field. Those fields belong in THIS matrix "
              "when they are adopted rather than being assumed independent of abundance."]
    args.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
