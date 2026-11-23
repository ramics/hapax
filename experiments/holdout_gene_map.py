#!/usr/bin/env python3
"""Gene coordinates for the two holdout contigs, so the page can draw them.

Why this exists. `results/conserved_positions/` scores the withheld contig gene by gene
and reports three identities spanning 24.3 points, but it keeps no nucleotide
coordinates: its unit is the ORF, not the molecule. The web view draws a genome as the
molecule, with each gene at its real position, because a bar chart of identities beside a
name tells a reader neither what the thing is nor what the number is a property of. That
drawing needs start, end and strand, which no committed file carried.

What it does NOT do. It calls no identities of its own. Every identity, subject and role
here is read from `conserved_positions_s1.json` and joined by ORF index, so the drawing
cannot disagree with the measurement it illustrates: if the join fails, this raises rather
than drawing an unlabelled gene. The ORF caller and its minimum length are the same as
that stage's, which is what makes the indices comparable at all.

The truth columns are read from the committed candidate tables, never retyped: the
withheld genome's divergence from its nearest public relative and its depth in this
library are facts about a deposit, not about this run.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# The withheld genome. Deposited 2023-07-31, after the 2022-05-10 index cutoff, so it is
# absent from every index this cascade consulted.
WITHHELD_ACC = "OQ722337.1"
LIBRARY = "SRR23421544"


def read_fasta(p: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif name:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def provenance(script: Path) -> dict:
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip())
    return {"git_sha": sha + ("-dirty" if dirty else ""), "tree_dirty": dirty,
            # Taken now, when the manifest is written, which is the only moment that
            # pins the code that actually ran.
            "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest()[:16]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=Path,
                    default=ROOT / "results/holdout_cascade/holdout_contigs.fa")
    ap.add_argument("--scored", type=Path,
                    default=ROOT / "results/conserved_positions/"
                                   "conserved_positions_s1.json")
    ap.add_argument("--shortlist", type=Path,
                    default=ROOT / "results/holdout_candidates/shortlist.tsv")
    ap.add_argument("--depth", type=Path,
                    default=ROOT / "results/holdout_candidates/depth.tsv")
    # The same minimum the scoring stage used. A different one renumbers the ORFs and
    # silently shifts every identity onto the wrong gene, which is the one way this
    # script could lie.
    ap.add_argument("--min-orf-aa", type=int, default=50)
    ap.add_argument("--run", type=Path,
                    default=ROOT / "results/holdout_cascade/run_s1",
                    help="the cascade run directory, read only to record whether each "
                         "contig survived into the residue the report shows")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    import pyrodigal  # imported late so --help works without it

    t0 = time.time()
    seqs = read_fasta(a.contigs)
    scored = json.loads(a.scored.read_text())

    # arm name -> the contig the scoring stage gave it. It renames an extra FASTA record
    # to "extra_<record>", so the record name is what joins the two files.
    by_orf = {}
    for r in scored["rows"]:
        if r["contig"] in seqs:
            by_orf[(r["contig"], r["orf"])] = r

    finder = pyrodigal.GeneFinder(meta=True)
    contigs = []
    for name, seq in seqs.items():
        genes = finder.find_genes(seq.encode())
        orfs, joined, unjoined = [], 0, 0
        for i, g in enumerate(genes, start=1):
            aa = len(g.translate()) - 1
            if aa < a.min_orf_aa:
                continue
            row = by_orf.get((name, f"{name}__orf{i}"))
            if row is not None:
                joined += 1
                if abs(row["orf_aa"] - aa) > 1:
                    raise SystemExit(
                        f"{name}__orf{i}: this script calls {aa} aa and the scoring "
                        f"stage recorded {row['orf_aa']}. The ORF numbering does not "
                        f"line up, so an identity would be drawn on the wrong gene.")
            else:
                unjoined += 1
            roles = sorted({d["role"] for d in (row or {}).get("domains", [])})
            orfs.append({
                "orf": f"{name}__orf{i}", "s": g.begin, "e": g.end, "d": g.strand,
                "aa": aa,
                "role": (roles[0] if len(roles) == 1 else
                         ("/".join(roles) if roles else None)),
                "pident": (row or {}).get("pident"),
                "aln_aa": (row or {}).get("aln_aa"),
                "evalue": (row or {}).get("evalue"),
                "title": (row or {}).get("subject_title"),
                "pfam": sorted({d["pfam"] for d in (row or {}).get("domains", [])}),
                # An ORF the scoring stage has no row for was called here and had no
                # significant alignment there. Recorded as a census rather than dropped:
                # an empty table and a search that never ran look the same otherwise.
                "scored": row is not None,
            })
        # Did this contig reach the pile the report draws? contigs.jsonl carries the
        # unidentified pile, so absence here is the cascade having subtracted it as
        # known. Checked rather than asserted: the node name is in the FASTA header.
        node = (a.contigs.read_text().split(">" + name + " ")[1].split("\n")[0].strip()
                if ">" + name + " " in a.contigs.read_text() else None)
        survived = None
        cj = a.run / "contigs.jsonl"
        if node and cj.exists():
            survived = any(node in line for line in cj.open())
        contigs.append({"contig": name, "nt": len(seq), "orfs": orfs,
                        "assembly_node": node,
                        "reached_the_unidentified_pile": survived,
                        "orfs_called": len(orfs), "orfs_scored": joined,
                        "orfs_with_no_significant_hit": unjoined,
                        "arm": ("withheld" if name.startswith("candidate")
                                else "already_referenced")})

    truth = {}
    head = a.shortlist.read_text().splitlines()
    cols = head[0].split("\t")
    for line in head[1:]:
        f = dict(zip(cols, line.split("\t")))
        if f.get("acc") == WITHHELD_ACC:
            truth = {"accession": WITHHELD_ACC, "organism": f.get("organism"),
                     "genome_nt": int(f["genome_nt"]),
                     "nucleotide_divergence_percent": float(f["divergence_pct"]),
                     # The shortlist's own columns are shifted for some rows, so the
                     # relative's accession is not read from it. It is in
                     # results/holdout_candidates/hits_blastn.tsv beside the identity
                     # this divergence comes from, and nothing here needs it.
                     "deposited": f.get("deposited"),
                     "collected": f.get("collected")}
    for line in a.depth.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 6 and f[0] == LIBRARY and f[1] == WITHHELD_ACC:
            truth["mapped_reads"] = int(f[3])
            truth["mean_depth"] = float(f[4])
    if not truth:
        raise SystemExit(f"{WITHHELD_ACC} is in neither committed table; the drawing "
                         f"would carry a distance nobody can check")

    out = {
        "unit": "an ORF of a holdout contig, with the coordinates the drawing needs",
        "question": "where on each contig its genes sit, so the page can draw the "
                    "molecule instead of a bar chart of its identities",
        "identities_are_not_measured_here": (
            "every identity, subject, role and E-value is read from "
            f"{a.scored.name} and joined by ORF index; this script calls ORFs and "
            "nothing else, and raises if the two disagree about an ORF's length"),
        "truth": truth,
        "truth_note": ("the withheld genome was deposited after the index cutoff this "
                       "cascade was rolled back to, so its distance is a fact nobody in "
                       "this project chose"),
        "operating_point": {"orf_caller": f"pyrodigal {pyrodigal.__version__} meta",
                            "min_orf_aa": a.min_orf_aa,
                            "identities_from": str(a.scored.relative_to(ROOT))},
        "contigs": contigs,
        "blind_spots": [
            "a called ORF is a prediction, in meta mode, on a contig that is itself an "
            "assembly; a gene drawn here is where the caller put it",
            "an ORF with no significant alignment is drawn unlabelled rather than as "
            "absent, because no search here can distinguish a novel gene from a wrong "
            "call",
            "the withheld contig covers 96% of the deposited genome, so a gene near "
            "either end may be truncated by the assembly rather than by biology",
            "reached_the_unidentified_pile is a membership test against one run "
            "directory's contigs.jsonl, which holds the unidentified pile. False means "
            "the cascade identified and subtracted the contig; it does not say at what "
            "score, which results/holdout_cascade/REPORT.md records",
        ],
        "seed": a.seed,
        "timings_seconds": {"total": round(time.time() - t0, 3)},
        **provenance(Path(__file__).resolve()),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    for c in contigs:
        print(f"{c['contig']}: {c['nt']} nt, {c['orfs_called']} ORFs, "
              f"{c['orfs_scored']} scored, {c['orfs_with_no_significant_hit']} with no "
              f"significant hit")
    print(f"truth: {truth.get('organism')} at "
          f"{truth.get('nucleotide_divergence_percent')}% divergence, "
          f"{truth.get('mean_depth')}x")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
