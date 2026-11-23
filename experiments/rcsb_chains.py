#!/usr/bin/env python3
"""Resolve a Foldseek subject name to the molecule that was actually aligned to.

WHY THIS IS A MODULE AND NOT A FEW LINES IN A REPORT. A Foldseek target is a CHAIN
(`8GH6_A`), a solved PDB entry is usually several different molecules, and the entry-level
answer names whichever entity the data API happened to return first. That is the
"ties broken by array order" shape this repo has already been bitten by, and it lands in the
one field whose entire purpose is the name. Measured on this experiment's first output,
before the fix: `8GH6_A` was reported as "28S DNA top strand" and is "Reverse
transcriptase-like protein"; `7Z4O_A` was reported as an RNA oligomer and is "Polymerase
acidic protein"; `10QM_R` was reported as L21 and is uL13. Placement RATES never consult
this and are unaffected - it is the "and like what" half that was wrong.

THREE TRAPS, ALL LIVE, ALL FOUND BY THE REFUSAL RATHER THAN BY A TABLE THAT LOOKED FINE.

- *The chain id convention.* GraphQL instance ids are `ENTRY.asym_id`, the mmCIF label,
  and Foldseek writes the AUTH chain. In 9TS7 asym "A" is the 26S rRNA while auth "A" is
  ribosomal protein S2, which is the chain Foldseek held. Resolving against asym_id first
  therefore names a different molecule with no error anywhere. Auth is preferred here.

- *Protein only, as a tripwire and not as a filter.* Every chain in a Foldseek structure
  library is a protein chain by construction, because 3Di comes from backbone geometry it
  only reads for protein. So a target that resolves to an RNA or DNA entity has not found a
  nucleic subject - it has resolved to the wrong molecule. Registering nucleic entities let
  a label-id collision answer in place of the real chain: it put 33 of 1,590 chains onto
  rRNA and DNA records. With them excluded the census reads 100% protein, which is what the
  library is, so the census is now a check on this module rather than a result.

- *The NMR ensemble.* An NMR entry is deposited as many models, so Foldseek writes
  `2MTV_MODEL_10_A` and the chain lands at the END of the name. Splitting on the first
  underscore gave 40 of them a chain of "MODEL".

A chain that still does not resolve RAISES. A wrong name is worse than a missing one, and
the fallback is what produced every error above.

TWO AXES, NEVER ONE. `description` is what the record IS - a polymerase, a capsid - and
`source_organisms` is where that record came from, which is not what a query infects and
not what is in the water.
"""
from __future__ import annotations

import collections
import json
import re
import subprocess
from pathlib import Path

RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

QUERY = ('{{entries(entry_ids:[{ids}]){{rcsb_id struct{{title}} '
         'polymer_entities{{rcsb_id '
         'rcsb_polymer_entity{{pdbx_description}} '
         'entity_poly{{rcsb_entity_polymer_type}} '
         'rcsb_entity_source_organism{{ncbi_scientific_name}} '
         'polymer_entity_instances{{rcsb_id '
         'rcsb_polymer_entity_instance_container_identifiers{{asym_id auth_asym_id}}'
         '}}}}}}}}')


def split_target(target: str) -> tuple[str, str]:
    """`8GH6_A` -> ("8GH6", "A"), and the NMR shape `2MTV_MODEL_10_A` -> ("2MTV", "A").

    An NMR entry is deposited as an ensemble, so foldseek writes one subject per model and
    the chain lands at the END of the name rather than straight after the entry. Parsing on
    the first underscore gave those a chain of "MODEL", which resolved to nothing - which is
    how the 40 of them in the first run were found, by the refusal below rather than by a
    table that looked fine.
    """
    # The filename segment is REMOVED, not truncated at: a Foldseek library built from
    # files rather than ids writes `8GH6.cif.gz_A`, and cutting the name at ".cif" takes
    # the chain with it. That drops silently into the chainless path below, which is a
    # different and wrong answer.
    stem = re.sub(r"\.(cif|pdb|ent|mmcif)(\.gz)?", "", target, flags=re.I)
    stem = re.sub(r"_MODEL_\d+", "", stem)
    entry, _, chain = stem.partition("_")
    return entry.upper(), chain


def annotate_subjects(targets: list[str], out: Path) -> dict:
    """Per-CHAIN description, polymer type and source organism for the chains actually hit.

    Keyed by the Foldseek target name, so a caller asks about the chain that aligned and
    cannot be handed a different molecule from the same entry. Fetched in batches of 100
    entries, so the cost is bounded by the result rather than by the library.
    """
    want = {t: split_target(t) for t in sorted(set(targets))}
    entries = sorted({e for e, _ in want.values()})
    by_chain: dict[tuple[str, str], dict] = {}
    by_entity: dict[str, list[dict]] = {}
    titles: dict[str, str] = {}
    for i in range(0, len(entries), 100):
        batch = entries[i:i + 100]
        ids = ", ".join(f'"{e}"' for e in batch)
        payload = out / f"gql_{i}.json"
        payload.write_text(json.dumps({"query": QUERY.format(ids=ids)}))
        done = subprocess.run(
            ["curl", "-sS", "--max-time", "180", "-X", "POST", RCSB_GRAPHQL,
             "-H", "Content-Type: application/json", "--data", f"@{payload}"],
            capture_output=True, text=True)
        if done.returncode != 0 or not done.stdout.strip():
            raise RuntimeError(
                f"the RCSB data API returned nothing for entries {i}+. Unannotated "
                f"placements would read as placements onto nothing, which is a different "
                f"claim. curl said: {done.stderr[-400:]}")
        body = json.loads(done.stdout)
        if body.get("errors"):
            raise RuntimeError(f"the RCSB data API refused the query: {body['errors'][:2]}")
        for entry in (body.get("data", {}).get("entries") or []):
            if not entry:
                continue
            eid = entry["rcsb_id"].upper()
            titles[eid] = (entry.get("struct") or {}).get("title")
            for pe in (entry.get("polymer_entities") or []):
                desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")
                ptype = (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type")
                orgs = sorted({so["ncbi_scientific_name"]
                               for so in (pe.get("rcsb_entity_source_organism") or [])
                               if so.get("ncbi_scientific_name")})
                row = dict(description=desc, polymer_type=ptype,
                           source_organisms=orgs[:6], entity=pe.get("rcsb_id"),
                           entry_title=titles[eid])
                # PROTEIN ONLY, and that is a tripwire rather than a filter. Every chain
                # in this subject library is a protein chain by construction: foldseek
                # builds 3Di from backbone geometry it only reads for protein. So a target
                # that resolves to an RNA or DNA entity has not found a nucleic subject, it
                # has resolved to the wrong molecule - which is exactly what happened on
                # 9TS7, where asym "A" is the 26S rRNA and auth "A" is ribosomal protein
                # S2, the chain foldseek actually held. Registering the nucleic entity would
                # have let a label-id collision answer in place of the real chain, silently.
                if ptype != "Protein":
                    continue
                by_entity.setdefault(eid, []).append(row)
                for inst in (pe.get("polymer_entity_instances") or []):
                    ids_ = (inst.get(
                        "rcsb_polymer_entity_instance_container_identifiers") or {})
                    # Auth first: a foldseek target name carries the auth chain, and the
                    # label id of a different entity can collide with it.
                    for key in ("auth_asym_id", "asym_id"):
                        ch = ids_.get(key)
                        if ch and (eid, ch) not in by_chain:
                            by_chain[(eid, ch)] = row
        payload.unlink(missing_ok=True)

    rows: dict[str, dict] = {}
    missing = []
    for target, (eid, chain) in want.items():
        if chain:
            hit = by_chain.get((eid, chain))
        else:
            # No chain in the target name: unambiguous only for a single-entity entry.
            ents = by_entity.get(eid) or []
            hit = ents[0] if len(ents) == 1 else None
        if hit is None:
            missing.append(target)
            continue
        rows[target] = hit
    if missing:
        raise RuntimeError(
            f"{len(missing)} of {len(want)} hit chains did not resolve to a polymer entity, "
            f"e.g. {missing[:5]}. The old code answered these with the entry's first entity, "
            f"which is how a reverse transcriptase came to be reported as '28S DNA top "
            f"strand'. A wrong name is worse than a missing one here.")
    census = collections.Counter(r["polymer_type"] for r in rows.values())
    return {"by_target": rows, "polymer_type_census": dict(census),
            "entry_titles": titles, "chains": len(rows)}
