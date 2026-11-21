"""Aggregate the residue-contig annotation seeds into one table. Tables are generated.

The layout puts the controls to the right of the measurement and the shuffle null last,
because the reading order this project keeps getting wrong is the one where the test is
read first. Three of the four columns here are controls, and two of them move with the
test in ways that change what the measurement means.

Refusals match the other aggregators: differing or missing `script_sha256` is a mixture
and is refused; a dirty tree is recorded rather than refused, because the content hash is
what pins the code.
"""

from __future__ import annotations

import argparse
import glob
import json
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

RATE_FIELDS = [
    ("any_pfam", "any Pfam family at its gathering threshold"),
    ("rdrp", "an RNA-dependent RNA polymerase family"),
    ("any_rfam", "any Rfam family at its gathering threshold"),
    ("coding_density_above_half", "at least half the contig inside a predicted CDS"),
]


def pcs(vals) -> str:
    xs = [100 * v for v in vals if v is not None]
    if not xs:
        return "not measured"
    if len(xs) == 1 or max(xs) - min(xs) < 1e-9:
        return f"{xs[0]:.2f}"
    return f"{st.mean(xs):.2f} ({min(xs):.2f}-{max(xs):.2f})"


def nums(vals) -> str:
    xs = [v for v in vals if v is not None]
    if not xs:
        return "not measured"
    if len(xs) == 1 or max(xs) - min(xs) < 1e-9:
        return f"{xs[0]:.3f}"
    return f"{st.mean(xs):.3f} ({min(xs):.3f}-{max(xs):.3f})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path,
                    default=ROOT / "results/residue_contig_annotation")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.results / "residue_contig_annotation_s*.json")))
    top = sorted(glob.glob(str(args.results / "top200" /
                                "residue_contig_annotation_s*.json")))
    if not files:
        raise SystemExit(f"no seeds in {args.results}")
    ds = [json.loads(Path(f).read_text()) for f in files]
    hashes = {d.get("script_sha256") for d in ds} | {
        json.loads(Path(f).read_text()).get("script_sha256") for f in top}
    if None in hashes or "" in hashes:
        raise SystemExit("a seed carries no script_sha256")
    if len(hashes) > 1:
        raise SystemExit(f"these runs came from different scripts: "
                         f"{sorted(h[:12] for h in hashes)}; re-run every seed")
    dirty = sum(1 for d in ds if d.get("tree_dirty"))
    piles = list(ds[0]["piles"].keys())
    par = ds[0]["parameters"]
    reused = sum(1 for d in ds if d.get("residue_pile_reused"))

    def block(dss: list[dict], title: str, note: str) -> list[str]:
        ps = list(dss[0]["piles"].keys())
        L = [f"## {title}", "", note, "",
             "| | " + " | ".join(ps) + " |", "|---" * (len(ps) + 1) + "|"]
        L.append("| contigs scored | " +
                 " | ".join(str(dss[0]["piles"][p]["contigs"]) for p in ps) + " |")
        L.append("| median length nt | " +
                 " | ".join(nums([d["piles"][p]["length"]["median"] for d in dss])
                            for p in ps) + " |")
        L.append("| coding density, mean | " +
                 " | ".join(nums([d["piles"][p]["coding_density"]["mean"] for d in dss])
                            for p in ps) + " |")
        for key, label in RATE_FIELDS:
            L.append(f"| {label}, % | " +
                     " | ".join(pcs([(d["piles"][p].get(key) or {}).get("point")
                                     for d in dss]) for p in ps) + " |")
        L += ["", "Wilson 95% intervals on the single-seed counts are in the JSONs; the "
                  "spread above is across seeds, which is the stronger statement where "
                  "there is more than one.", ""]
        fam = {}
        for p in ps:
            c: dict[str, int] = {}
            for d in dss:
                for name, k in d["piles"][p]["top_pfam"]:
                    c[name] = c.get(name, 0) + k
            fam[p] = sorted(c.items(), key=lambda kv: -kv[1])[:6]
        L += ["Most frequent Pfam families, pooled over seeds:", ""]
        for p in ps:
            L.append(f"- **{p}**: " + (", ".join(f"{n} x{k}" for n, k in fam[p])
                                       or "none"))
        rf = {}
        for p in ps:
            c = {}
            for d in dss:
                for name, k in d["piles"][p]["top_rfam"]:
                    c[name] = c.get(name, 0) + k
            rf[p] = sorted(c.items(), key=lambda kv: -kv[1])[:6]
        L += ["", "Rfam families:", ""]
        for p in ps:
            L.append(f"- **{p}**: " + (", ".join(f"{n} x{k}" for n, k in rf[p])
                                       or "none"))
        short = dss[0].get("length_match_shortfall") or {}
        if short:
            L += ["", "**Length matching did not fully succeed for every control.** A "
                      "control pool with fewer long contigs than the query pile cannot "
                      "fill its upper bins, so that control is matched below the top of "
                      "the length range only:"]
            for tag, bins in short.items():
                w = sum(b["wanted"] for b in bins.values())
                a = sum(b["available"] for b in bins.values())
                L.append(f"- `{tag}`: {w - a} contigs short across {len(bins)} length "
                         f"bins, so it is smaller than the query pile and leans short.")
        exc = dss[0].get("query_contigs_removed_from_control_pool") or {}
        if exc:
            L += ["", "Query contigs removed from their own control pool before "
                      "sampling: " + ", ".join(f"{k} ({v})" for k, v in exc.items()) +
                  ". Without this the query would sit inside its own control."]
        return L + [""]

    L = ["# Annotating the residue contigs", "",
         f"Generated by `experiments/report_residue_contig_annotation.py`. Script "
         f"`{list(hashes)[0][:12]}`. Prodigal in metagenomic mode; Pfam-A at each family's "
         f"own gathering threshold; {par['rfam_settings']}. Controls are "
         f"{par['length_matching']}."
         + (f" {dirty} of {len(ds)} seeds ran on a dirty tree, so their commit names HEAD; "
            f"the script hash pins them." if dirty else "")
         + (f" The residue pile is deterministic - every contig is annotated and nothing "
            f"is sampled - so it is computed once and reused by {reused} of {len(ds)} "
            f"seeds, which vary only the controls and the shuffle." if reused else ""),
         ""]

    if top:
        L += block([json.loads(Path(f).read_text()) for f in top],
                   "The 200 highest-mass unidentifiable contigs",
                   "The head of the residue's read-mass distribution: long enough and deep "
                   "enough for a profile search to have a chance, and the place any "
                   "annotation signal should appear first.")
    L += block(ds, "The residue-only assembly, new contigs",
               "Contigs assembled from unplaced pairs that are under 10% covered by the "
               "sample's own assembly, so not a re-derivation of sequence already in it.")

    L += ["", "---", "", f"Provenance: {json.dumps(provenance(__file__))}", ""]
    out = args.out or (args.results / "TABLES.md")
    out.write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwritten {out}")


if __name__ == "__main__":
    main()
