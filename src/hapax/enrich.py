"""The three enrichment stages, run by the package. PLAN.md step 6.

WHAT THESE ARE, AND WHY THEY ARE NOT THE OTHER THREE. Writing the registry's declarations
down found six files the report reads that no stage wrote. Three of them -
`public_health_profile_s*.json`, `known_strain_divergence_s*.json`, `ranking.json` - were
sections rendering a census of an absent file, and are ordinary stages now: on real water
that page told its reader something was missing from their sample when what was missing
was a stage of ours.

These three are different and the difference decides what the report does. The sections
that read them render WITHOUT them and say what that costs: the distance panel ranks on
the whole-contig best hit rather than on the conserved gene, and the subject table shows
an accession where it could show an organism. So they are declared optional, with the
consequence of each absence stated, and a run that cannot produce them produces a report
that says so rather than refusing.

Calling an enrichment a hole would make every run without a Pfam library refuse. Calling a
hole an enrichment is what the page was doing when a run with three empty sections built
cleanly. Both are wrong, and only a declaration can tell them apart.

WHY A STAGE REFUSES SOFTLY HERE AND NOWHERE ELSE. Each of these needs a database the
cascade does not: the per-gene workup needs a Pfam HMM library and a Swiss-Prot protein
index, and the name map needs the FASTA the protein index was built from rather than the
index itself. A run given none of those is a legitimate run of a tool whose job is the
cascade, so the absence is RECORDED - what was missing, what the page shows instead -
rather than raised. The record is the difference between this and a fallback.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from .config import Config


def _attempt(cmd: list[str], cwd: Path, stage: str, needs: list[str]) -> dict:
    """Run one enrichment, recording the outcome either way."""
    absent = [n for n in needs if not Path(n).exists()]
    if absent:
        return {"stage": stage, "status": "skipped",
                "reason": f"these are not staged on this machine: {', '.join(absent)}"}
    proc = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        return {"stage": stage, "status": "failed",
                "reason": (proc.stderr or proc.stdout)[-1500:],
                "argv": [str(c) for c in cmd]}
    return {"stage": stage, "status": "ok", "argv": [str(c) for c in cmd]}


def enrich(cfg: Config, out: Path, repo: Path, seed: int) -> list[dict]:
    """Run the three, in their declared order, and write what happened.

    The order is the registry's: the architecture call reads the per-gene workup, so it is
    absent wherever the workup is, and recording that explicitly is why it is attempted
    rather than skipped on inspection.
    """
    records: list[dict] = []
    py = sys.executable
    strains = out / "strains.json"
    workup = out / f"divergent_contig_workup_s{seed}.json"

    pfam = cfg.databases.get("pfam")
    sprot = cfg.databases.get("sprot")
    contigs = None
    prepared = out / "prepared.json"
    if prepared.exists():
        contigs = json.loads(prepared.read_text()).get("assembly")

    if not (pfam and sprot):
        records.append({
            "stage": "workup", "status": "skipped",
            "reason": "the per-gene workup needs a Pfam HMM library and a Swiss-Prot "
                      "protein index; neither pfam nor sprot is named under [databases]"})
    else:
        records.append(_attempt(
            [py, "-I", repo / "experiments/divergent_contig_workup.py", "--seed", seed,
             "--strains", strains, "--contigs", contigs or "",
             "--pfam", pfam, "--dmnd-sprot", sprot,
             "--dmnd", cfg.database("protein"), "--threads", cfg.threads,
             # A file, and a work directory beside it: this script writes the path it is
             # given rather than a name of its own inside it, and `--work` is required.
             "--work", out / "work" / "workup",
             "--out", out / f"divergent_contig_workup_s{seed}.json"],
            repo, "workup", [str(strains), str(pfam), str(sprot)]))

    if workup.exists():
        records.append(_attempt(
            [py, "-I", repo / "experiments/architecture_call.py", "--seed", seed,
             "--scored", workup, "--out", out / f"architecture_s{seed}.json"],
            repo, "architecture", [str(workup)]))
    else:
        records.append({"stage": "architecture", "status": "skipped",
                        "reason": "it reads the per-gene workup, which this run has not "
                                  "produced"})

    # The name map reads the FASTA the protein index was BUILT FROM, not the index: a
    # second source for the same names could disagree with the index the search used, and
    # a label that does not come from the searched database is a label nobody can check.
    sources = cfg.databases.get("protein_fasta")
    if not sources:
        records.append({"stage": "names", "status": "skipped",
                        "reason": "no protein_fasta under [databases]: the name map is "
                                  "read from the FASTA the protein index was built from, "
                                  "and a run given only the index cannot have it"})
    else:
        records.append(_attempt(
            [py, "-I", repo / "experiments/subject_names.py", out,
             "--sources", sources, "--out", out / "subject_names.json"],
            repo, "names", [str(sources)]))

    (out / "enrichments.json").write_text(json.dumps({
        "what_this_is": "the enrichment stages, which the sections that read them render "
                        "without. An absence here is recorded with what the page shows "
                        "instead; it is not a hole in the report.",
        "seed": seed,
        "stages": records,
    }, indent=1) + "\n")
    return records
