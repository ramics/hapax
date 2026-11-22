"""Build the query sets for a big-index pass over the divergent contigs.

WHY A SEPARATE QUERY BUILDER AND NOT A SEPARATE LAUNCHER. `ec2_bigindex.py` and
`bigindex_userdata.sh` already stage core_nt and nr, measure their own fixed scan cost
from a probe before committing to a pass, chunk the work against a dead-man timer and
terminate in a `finally`. None of that is specific to reads. The only thing that was
specific to reads is the query files, so this writes query files of the same names from a
different source and the launcher is reused unchanged.

WHAT IS BEING ASKED. Every name on the divergence panel came from `data/prot/viral.dmnd`
plus Swiss-Prot: RefSeq viral proteins and a curated protein set. Those are narrow. A
contig at 35% identity to the closest RefSeq viral protein may sit at 90% to something in
nr or core_nt that nobody has curated into a viral reference set, and the two readings
mean completely different things. This pass is what separates them.

THE CONTROL SHIPS WITH IT, as everywhere else in this repo: the near-reference contigs
from the same assembly go through the same search in the same chunks. Without them a
column of identities has nothing to be read against, and this project has already had to
relabel such a column.

COST SHAPE, which is the reason this is worth running at all. `results/subtract_known/`
measured BLAST against core_nt at 1,082 s fixed plus 1.97 s per 107 nt read. The fixed
part is a scan of the database and does not care how many queries follow it. A few
hundred contigs is therefore nearly free once the machine is up, where the same pass over
8.7 million reads is impossible. A fixed cost masquerading as a rate is a trap this repo
has recorded; here it works in our favour and the pilot measures it again anyway.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_fasta(path: Path) -> dict:
    seqs, name, buf = {}, None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write(path: Path, names: list, seqs: dict, arms: dict) -> int:
    # The arm travels in the name, as `<arm>|<contig>`, so that the instance-side script
    # stays a search wrapper with no reference data and the arms can be separated again
    # from the results alone.
    path.write_text("".join(f">{arms[n]}|{n}\n{seqs[n]}\n" for n in names))
    return len(names)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workup", type=Path, required=True,
                    help="divergent_contig_workup output JSON")
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--out", type=Path, required=True, help="query directory to write")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--probe", type=int, default=2,
                    help="contigs in scan_probe.fa, which measures the fixed scan cost "
                         "alone; keep it tiny or it measures the pass instead")
    ap.add_argument("--pilot", type=int, default=12)
    ap.add_argument("--eyeball", type=int, default=6)
    a = ap.parse_args()

    w = json.loads(a.workup.read_text())
    arms = {r["contig"]: ("divergent" if r["arm"] == "divergent" else "control")
            for r in w["rows"]}
    seqs = read_fasta(a.contigs)
    names = [n for n in arms if n in seqs]
    missing = [n for n in arms if n not in seqs]
    if missing:
        print(f"{len(missing)} contigs not in {a.contigs}", file=sys.stderr)

    div = sorted(n for n in names if arms[n] == "divergent")
    con = sorted(n for n in names if arms[n] == "control")
    rng = random.Random(a.seed)

    a.out.mkdir(parents=True, exist_ok=True)
    # The probe and pilot are drawn from both arms, because a pilot that is all divergent
    # sequence measures the wrong per-query cost: a query with no hits returns sooner.
    def mix(k):
        pick = rng.sample(div, min(k - k // 2, len(div)))
        pick += rng.sample(con, min(k // 2, len(con)))
        return sorted(pick)

    counts = {
        "scan_probe.fa": write(a.out / "scan_probe.fa", mix(a.probe), seqs, arms),
        "pilot.fa": write(a.out / "pilot.fa", mix(a.pilot), seqs, arms),
        "eyeball.fa": write(a.out / "eyeball.fa", mix(a.eyeball), seqs, arms),
        "block_001.fa": write(a.out / "block_001.fa", div + con, seqs, arms),
        "nr_sample.fa": write(a.out / "nr_sample.fa", div + con, seqs, arms),
    }
    manifest = {
        "built_from": str(a.workup),
        "contigs": {"divergent": len(div), "control_near_reference": len(con)},
        "files": counts,
        "name_format": "<arm>|<contig>, arm in {divergent, control}",
        "seed": a.seed,
        "contigs_not_found": missing,
        "sha256_of_block": hashlib.sha256(
            (a.out / "block_001.fa").read_bytes()).hexdigest()[:16],
    }
    (a.out / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    print(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
