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
from pathlib import Path

from . import scripts, tools
from .config import Config


def _attempt(repo: Path, name: str, args: list, stage: str, needs: list[str]) -> dict:
    """Run one enrichment, recording the outcome either way."""
    absent = [n for n in needs if not Path(n).exists()]
    if absent:
        return {"stage": stage, "status": "skipped",
                "reason": f"these are not staged on this machine: {', '.join(absent)}"}
    proc, code = scripts.run(repo, name, args)
    if proc.returncode != 0:
        return {"stage": stage, "status": "failed",
                "reason": (proc.stderr or proc.stdout)[-1500:], **code}
    return {"stage": stage, "status": "ok", **code}


def enrich(cfg: Config, out: Path, repo: Path, seed: int) -> list[dict]:
    """Run the three, in their declared order, and write what happened.

    The order is the registry's: the architecture call reads the per-gene workup, so it is
    absent wherever the workup is, and recording that explicitly is why it is attempted
    rather than skipped on inspection.
    """
    records: list[dict] = []
    # THE TAXA STAGE'S OUTPUT, not `strains.json`. Two files in one run directory answer
    # to the word strains: the `strains` stage writes a read-level table of SUBJECTS, and
    # `known_strain_divergence.py` writes the per-CONTIG table of identities this script's
    # `--strains` documents and its `pick_contigs` reads under a `taxa` key. The read-level
    # table has no such key, so selection returned nothing, the ORF file was empty, and
    # DIAMOND exited 1 on a blank query 11 stages into a paid run - reported as `0
    # divergent contigs, 0 near-reference controls`, which reads like a property of the
    # water and was a property of the filename.
    strains = out / f"known_strain_divergence_s{seed}.json"
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
            repo, "divergent_contig_workup.py", ["--seed", seed,
             "--strains", strains, "--contigs", contigs or "",
             "--pfam", pfam, "--dmnd-sprot", sprot,
             "--dmnd", cfg.database("protein"), "--threads", cfg.threads,
             # Always passed, never left to the script's default, which is the layout of
             # one container written where a parameter goes: `/home/user/tools/diamond`.
             # On an instance that is not there, and the workup failed with a
             # FileNotFoundError 11 stages into a paid run. Same rule as prepare.split.
             "--diamond-bin", tools.resolve(cfg, "diamond"),
             # THE READS AND THE TRIM, because this script maps the selected contigs back
             # to the library and its own defaults are the discriminating dataset's FASTQ
             # paths. The driver guards its own stages against exactly this with two
             # tables (LIBRARY_AWARE and SAMPLE_SPECIFIC) asserted against the scripts'
             # source; the package runs this script and inherited neither, so seed 1 of
             # an unseen library called 76 genes, searched them, and then opened
             # `data/berlin/SRR23973193_1.fastq.gz`. `--library` goes with them so that
             # library_facts.resolve_args can check the paths against the accession
             # rather than filling a fact from whichever library its default names.
             # A trim that was NOT stated is left out rather than passed as the word
             # None: the script then resolves it from the library registry, which raises
             # for a library it does not know. Stating it is the override, and the
             # absence is the question, never a value.
             "--library", cfg.label,
             "--r1", cfg.inputs["r1"], "--r2", cfg.inputs["r2"],
             *(["--trim-r2", cfg.fact("trim_r2")]
               if cfg.fact("trim_r2") is not None else []),
             # HOW MANY MATES THE MAP-BACK STREAMS, stated rather than left at the
             # script's 0, which is the whole library. Every other stage of the run reads
             # one seed's block of this many pairs, and a mappy pass over 11.6 M pairs is
             # about half an hour per seed against about a minute here. The script records
             # what it streamed under `mates_actually_streamed`, so the operating point
             # travels with the depth. It streams from the head of each file rather than
             # from the seed's block, so this read support is one per-contig depth and not
             # a per-seed estimate.
             "--map-limit", cfg.run.get("pairs", 200_000),
             # A file, and a work directory beside it: this script writes the path it is
             # given rather than a name of its own inside it, and `--work` is required.
             "--work", out / "work" / "workup",
             "--out", out / f"divergent_contig_workup_s{seed}.json"],
            "workup", [str(strains), str(pfam), str(sprot)]))

    if workup.exists():
        records.append(_attempt(
            repo, "architecture_call.py", ["--seed", seed,
             "--scored", workup, "--out", out / f"architecture_s{seed}.json"],
            "architecture", [str(workup)]))
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
            repo, "subject_names.py", [out,
             "--sources", sources, "--out", out / "subject_names.json"],
            "names", [str(sources)]))

    (out / "enrichments.json").write_text(json.dumps({
        "what_this_is": "the enrichment stages, which the sections that read them render "
                        "without. An absence here is recorded with what the page shows "
                        "instead; it is not a hole in the report.",
        "seed": seed,
        "stages": records,
    }, indent=1) + "\n")
    return records
