"""Turn the accessions a run reports into the names a person can read.

WHY THIS IS A SCRIPT AND NOT A LOOKUP AT READ TIME. An accession is the right key for a
pipeline and the wrong label for a reader: `YP_009552769.1` says nothing, and the report
this project is building is meant to be read by someone deciding whether a thing is worth
a second look. The names come from the same local FASTA files the DIAMOND databases were
built from, so the label on the page and the sequence that was searched are the same
release by construction. Nothing is fetched over the network, which also means the mapping
is reproducible from the committed data manifest rather than from whatever NCBI returns
today.

WHAT IT DOES NOT DO. It does not taxonomically classify anything. The organism string is
whatever the record's own header says, which for a viral RefSeq protein is the source
organism in brackets and for Swiss-Prot is the OS= field. A header that says "Escherichia
phage" is a claim by the record, not by this tool, and a read matching that record at 60%
identity over 34 amino acids is not that organism - which is why the page prints the
identity and the E-value next to the name rather than the name alone.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_SOURCES = [ROOT / "data/prot/viral.1.protein.faa.gz",
                   ROOT / "data/prot/uniprot_sprot.fasta.gz"]

#: RefSeq protein headers put the source organism in brackets at the end.
BRACKET = re.compile(r"\[([^\[\]]+)\]\s*$")
#: Swiss-Prot headers are `sp|ACC|ENTRY Description OS=Organism OX=... GN=...`.
SPROT_OS = re.compile(r"\bOS=(.+?)(?=\s+[A-Z][A-Z]=|$)")


def parse_header(h: str) -> tuple[str, dict]:
    """Return (accession as the search reports it, {name, organism})."""
    h = h.lstrip(">").rstrip()
    if h.startswith("sp|") or h.startswith("tr|"):
        key = "|".join(h.split()[0].split("|")[:3])
        rest = h.split(None, 1)[1] if " " in h else ""
        org = SPROT_OS.search(rest)
        name = rest.split(" OS=")[0].strip()
        return key, {"name": name, "organism": org.group(1).strip() if org else None}
    first = h.split(None, 1)
    key = first[0]
    rest = first[1] if len(first) > 1 else ""
    m = BRACKET.search(rest)
    org = m.group(1).strip() if m else None
    name = BRACKET.sub("", rest).strip() if m else rest.strip()
    return key, {"name": name or None, "organism": org}


def wanted_from_run(rundir: Path) -> set[str]:
    want: set[str] = set()
    s = rundir / "strains.json"
    if s.exists():
        for row in json.loads(s.read_text()).get("subjects", []):
            if row.get("subject"):
                want.add(row["subject"])
    # The whole-residue search names subjects the block sample never saw - it is the
    # reason this run has subjects at all beyond the 6,000-mate sample - so its output is
    # a source of accessions here too, or the page would print bare accessions for exactly
    # the rows the search was run to find.
    for f in sorted(rundir.glob("full_residue_search_s*.json")):
        for row in (json.loads(f.read_text()).get("piles", {})
                    .get("residue", {}) or {}).get("subjects", []) or []:
            if row.get("subject"):
                want.add(row["subject"])
    r = rundir / "reads.jsonl"
    if r.exists():
        with open(r) as fh:
            for line in fh:
                try:
                    row = json.loads(line)
                except ValueError:
                    continue
                for k in ("closest_viral", "closest_sprot"):
                    h = row.get(k)
                    if h and h.get("subject"):
                        want.add(h["subject"])
    return want


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--sources", type=Path, nargs="+", default=DEFAULT_SOURCES)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    want = wanted_from_run(args.rundir)
    if not want:
        raise SystemExit(f"no subject accessions in {args.rundir}")
    found: dict[str, dict] = {}
    scanned = {}
    for src in args.sources:
        if not Path(src).exists():
            scanned[str(src)] = "absent"
            continue
        n = 0
        op = gzip.open if str(src).endswith(".gz") else open
        with op(src, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                n += 1
                key, rec = parse_header(line)
                if key in want and key not in found:
                    rec["source"] = Path(src).name
                    found[key] = rec
        scanned[str(src)] = f"{n} records"

    out = {"names": found, "wanted": len(want), "resolved": len(found),
           "sources": scanned,
           "note": "the organism string is the source record's own header, not a "
                   "classification of the read; identity and E-value belong beside it",
           **provenance(__file__)}
    dest = args.out or (args.rundir / "subject_names.json")
    dest.write_text(json.dumps(out, indent=1))
    print(f"{len(found)} of {len(want)} accessions named -> {dest}")


if __name__ == "__main__":
    main()
