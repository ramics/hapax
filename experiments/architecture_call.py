#!/usr/bin/env python3
"""Apply the committed vocabulary to scored contigs, and record what it could not call.

Why this exists, and it is a correction rather than a feature. The web view was saying
"a single-stranded RNA bacteriophage: a virus that infects bacteria, not people" from a
regular expression written inside the renderer. Every number on that page comes from a
committed result file with a seed and a script hash; that sentence came from nowhere, and
it is a larger claim than most of the numbers beside it. This repo's own convention is
that a figure is generated from committed files because prose drifts from outputs - the
convention was being applied to the numbers and not to the words.

So the rule now lives in `experiments/virus_vocabulary.json`, where it can be read and
corrected as data, and this script applies it and writes the answer out with the evidence
that fired. The renderer gets a field, not a judgement. It still cannot make the English
stop being hand-written; it makes it visible, attributable and wrong in one place.

Two things it deliberately does NOT do. It calls no architecture from length or from the
absence of evidence: a contig with no matching domain and no matching subject is recorded
as `null` with a reason, never as "not a phage". And it records the per-contig role
ordering - which role is furthest from the reference and which is closest - as a fact
rather than letting the page assert an agreement with `results/conserved_positions/`.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def provenance(script: Path) -> dict:
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip())
    return {"git_sha": sha + ("-dirty" if dirty else ""), "tree_dirty": dirty,
            "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest()[:16],
            "vocabulary_sha256": hashlib.sha256(
                (ROOT / "experiments/virus_vocabulary.json").read_bytes()
            ).hexdigest()[:16]}


def genes_of(row: dict) -> list[dict]:
    """One shape out of two input files, because they were written for different jobs.

    The divergent workup carries `orf_rows`; the holdout gene map carries `orfs`. Both
    name a role, a percent identity, a subject title and a Pfam list, which is everything
    the vocabulary's rule reads.
    """
    out = []
    for o in row.get("orf_rows") or row.get("orfs") or []:
        pf = o.get("pfam") or []
        out.append({
            "role": o.get("gene_role") or o.get("role"),
            "pident": o.get("pident"),
            "title": o.get("best_title") or o.get("title") or "",
            "pfam": [p["id"] if isinstance(p, dict) else p for p in pf],
            "significant": o.get("significant", o.get("pident") is not None),
        })
    return out


def call_one(genes: list[dict], vocab: dict) -> dict:
    archs = vocab["architectures"]
    for key, a in archs.items():
        rule = a["rule"]
        fired_pfam, fired_title = [], []
        for g in genes:
            for name in rule.get("any_pfam_name", []):
                if name in g["pfam"]:
                    fired_pfam.append({"role": g["role"], "pfam": name})
            for sub in rule.get("any_subject_title_substring", []):
                if sub.lower() in (g["title"] or "").lower():
                    fired_title.append({"role": g["role"], "substring": sub})
        n = len({(f["role"], f["pfam"]) for f in fired_pfam}) + len(fired_title)
        if n >= rule.get("minimum_genes_firing", 1):
            return {"architecture": key,
                    "evidence_pfam": fired_pfam, "evidence_title": fired_title,
                    "genes_firing": n}
    return {"architecture": None,
            "not_called_because": ("no gene carries a domain or a subject title that any "
                                   "architecture rule in the vocabulary names. This is "
                                   "an absence of evidence and not a call."),
            "evidence_pfam": [], "evidence_title": [], "genes_firing": 0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scored", type=Path, required=True,
                    help="a divergent_contig_workup_s*.json or a gene_map_s*.json")
    ap.add_argument("--vocabulary", type=Path,
                    default=ROOT / "experiments/virus_vocabulary.json")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    t0 = time.time()
    vocab = json.loads(a.vocabulary.read_text())
    src = json.loads(a.scored.read_text())
    rows = src.get("rows") or src.get("contigs") or []

    out_rows, counts, uncalled = [], {}, 0
    for r in rows:
        key = r.get("contig")
        genes = genes_of(r)
        c = call_one(genes, vocab)
        scored = [g for g in genes if g["significant"] and g["pident"] is not None
                  and g["role"]]
        scored.sort(key=lambda g: g["pident"])
        roles_present = sorted({g["role"] for g in scored})
        c.update({
            "contig": key,
            "arm": r.get("arm"),
            "genes_called": len(genes),
            "genes_with_a_role": len(scored),
            "roles_present": roles_present,
            "has_every_expected_role": (
                bool(c["architecture"]) and
                set(vocab["architectures"][c["architecture"]]["expected_roles"])
                <= set("coat" if x == "capsid" else x for x in roles_present)),
            # stated rather than asserted: the page should not claim this genome agrees
            # with the Berlin ordering, it should show which gene is where and cite it
            "most_diverged_role": scored[0]["role"] if scored else None,
            "least_diverged_role": scored[-1]["role"] if len(scored) > 1 else None,
            "identity_spread_points": (round(scored[-1]["pident"] - scored[0]["pident"], 1)
                                       if len(scored) > 1 else None),
            "nearest_subject_provenance": subject_provenance(genes, vocab),
        })
        if c["architecture"] is None:
            uncalled += 1
        counts[c["architecture"] or "not called"] = (
            counts.get(c["architecture"] or "not called", 0) + 1)
        out_rows.append(c)

    out = {
        "unit": "a scored contig",
        "question": "which architecture the committed vocabulary's rule calls for each "
                    "contig, and on what evidence",
        "nothing_here_is_measured": (
            "the rule and every word it attaches are hand-written background knowledge "
            "in experiments/virus_vocabulary.json. This script only applies them and "
            "records what fired. A contig with no matching evidence is null, never "
            "'not a phage'."),
        "source": str(a.scored),
        "vocabulary": str(a.vocabulary.relative_to(ROOT)),
        "census": {"contigs": len(out_rows), "by_architecture": counts,
                   "not_called": uncalled},
        "blind_spots": [
            "the rule fires on a single matching domain or subject title, so one "
            "mislabelled gene calls the whole contig; the evidence that fired is "
            "recorded per contig so a reader can see how thin a call is",
            "there is one architecture in the vocabulary, so this cannot distinguish "
            "between architectures, only recognise the one it knows or decline",
            "a subject title is a label on a database record. 'ssRNA phage' in a title "
            "is a claim made by whoever deposited it",
        ],
        "rows": out_rows,
        "seed": a.seed,
        "timings_seconds": {"total": round(time.time() - t0, 3)},
        **provenance(Path(__file__).resolve()),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"{len(out_rows)} contigs: " +
          ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    print(f"wrote {a.out}")


def subject_provenance(genes, vocab) -> str | None:
    sp = vocab.get("subject_provenance", {}).get("environmental_record")
    if not sp:
        return None
    sig = [g for g in genes if g["significant"] and g["title"]]
    if not sig:
        return None
    if all(any(s in g["title"] for s in sp["title_substrings"]) for g in sig):
        return "environmental_record"
    return None


if __name__ == "__main__":
    main()
