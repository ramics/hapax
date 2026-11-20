"""Build the protein database the read-level and contig-level arms search, with
the spike's whole family removed.

Why a family and not just the genome: a database that still holds a close
relative turns "detect novel coding structure" into "align to a known protein",
which is the thing every incumbent already does well. Removing the whole family
of the spike (Rhabdoviridae for the rabies spike) is the weakest honest version
of novelty: the detector may still use any other virus, any bacterium and any
host protein, but nothing from the spike's own family.

Removal is by NCBI taxonomy rather than by name matching, because RefSeq
organism strings include unclassified and strain-level names that no hand-built
keyword list covers. The script verifies afterwards that no protein of the spike
accession itself survives, which is the one check that must never fail.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

ORG = re.compile(r"\[([^\[\]]+)\]\s*$")


def descendants(nodes_dmp: Path, root: int) -> set[int]:
    children: dict[int, list[int]] = {}
    with nodes_dmp.open() as fh:
        for line in fh:
            f = line.split("\t|\t")
            tid, parent = int(f[0]), int(f[1])
            if tid != parent:
                children.setdefault(parent, []).append(tid)
    out, stack = {root}, [root]
    while stack:
        for c in children.get(stack.pop(), ()):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def names_for(names_dmp: Path, taxids: set[int]) -> set[str]:
    out = set()
    with names_dmp.open() as fh:
        for line in fh:
            f = line.split("\t|\t")
            if int(f[0]) in taxids:
                out.add(f[1].strip().lower())
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", type=Path, required=True, help="RefSeq viral protein faa.gz")
    ap.add_argument("--tax", type=Path, required=True, help="directory with nodes.dmp and names.dmp")
    ap.add_argument("--exclude-taxid", type=int, default=11270, help="11270 = Rhabdoviridae")
    # Protein accessions, not the genome accession: RefSeq protein headers
    # carry NP_/YP_ ids and never the genome's NC_ id, so checking the latter
    # passes vacuously. The default five are NC_001542's own proteins.
    ap.add_argument("--spike-prot", default="NP_056793,NP_056794,NP_056795,NP_056796,NP_056797",
                    help="protein accessions of the spike; none may survive the filter")
    ap.add_argument("--family-name", default="Lyssavirus",
                    help="a name string that must no longer appear in any header")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    taxids = descendants(args.tax / "nodes.dmp", args.exclude_taxid)
    banned = names_for(args.tax / "names.dmp", taxids)
    print(f"excluding {len(taxids)} taxa, {len(banned)} name strings", flush=True)

    kept = dropped = 0
    keep_this = False
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.faa, "rt") as fh, args.out.open("w") as out:
        for line in fh:
            if line.startswith(">"):
                m = ORG.search(line.rstrip())
                org = m.group(1).strip().lower() if m else ""
                keep_this = org not in banned
                kept += keep_this
                dropped += not keep_this
            if keep_this:
                out.write(line)
    print(f"kept {kept} proteins, dropped {dropped}", flush=True)

    # The checks that must not fail: none of the spike's own proteins survive,
    # and the family's name appears nowhere. Both are grep over headers, which
    # is cheap compared with the cost of silently searching a database that
    # still holds the answer.
    headers = [ln for ln in args.out.open() if ln.startswith(">")]
    bad = [a for a in args.spike_prot.split(",") if any(a in h for h in headers)]
    fam = sum(1 for h in headers if args.family_name.lower() in h.lower())
    print(f"spike proteins surviving: {len(bad)} {bad}; "
          f"headers still naming {args.family_name}: {fam}")
    if bad or fam:
        sys.exit("FATAL: the spike's family survived the filter")


if __name__ == "__main__":
    main()
