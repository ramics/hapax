"""Find real genome pairs at target divergences, so a ladder rung is never simulated.

The two-lineage rung (section 5) has a hole: the four pairs reachable from the
first 40 deposited rabies genomes sit at 2.0, 3.8, 7.3 and 9.7% divergence, and
the 4-7% band - exactly where the pooling effect turns erratic - was empty. The
sister thread proposed filling it by mutating one genome to order. CLAUDE.md
forbids that: truth is real published sequence, and a rung whose divergence was
dialled in cannot falsify a claim about real lineages, because the substitution
process that made it is the same process the assembler's k-mer matching sees.

The deposited data is not the limit here. There are ~2,840 complete rabies
genomes; the empty band was an artefact of having downloaded 40 of them. This
script downloads a larger sample, computes the full pairwise divergence matrix
with edlib, and reports the closest real pair to each target. It is a selection
helper, not an experiment: it has no seed because it is deterministic given the
accession list, and the accession list it writes is what makes the rung
reproducible.
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import time
from pathlib import Path

import edlib

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(term: str, retmax: int, out: Path) -> Path:
    """Accessions then sequences, in batches NCBI tolerates."""
    if out.exists():
        return out
    es = subprocess.run(["curl", "-sS", "-G", ESEARCH,
                         "--data-urlencode", "db=nuccore",
                         "--data-urlencode", f"term={term}",
                         "--data-urlencode", f"retmax={retmax}",
                         "--data-urlencode", "retmode=json"],
                        capture_output=True, text=True, check=True)
    ids = json.loads(es.stdout)["esearchresult"]["idlist"]
    print(f"{len(ids)} accessions", flush=True)
    with out.open("w") as fh:
        for i in range(0, len(ids), 100):
            batch = ",".join(ids[i:i + 100])
            r = subprocess.run(["curl", "-sS", "-G", EFETCH,
                                "--data-urlencode", "db=nuccore",
                                "--data-urlencode", f"id={batch}",
                                "--data-urlencode", "rettype=fasta",
                                "--data-urlencode", "retmode=text"],
                               capture_output=True, text=True, check=True)
            fh.write(r.stdout)
            time.sleep(0.5)
    return out


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.open():
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--term", default="Lyssavirus rabies[Organism] AND complete genome[Title] "
                                      "AND 11000:12500[SLEN]")
    ap.add_argument("--retmax", type=int, default=400)
    ap.add_argument("--targets", default="0.045,0.055,0.065")
    ap.add_argument("--min-len", type=int, default=11000)
    ap.add_argument("--max-n", type=float, default=0.01, help="max fraction of ambiguous bases")
    ap.add_argument("--fasta", type=Path, default=Path("/home/user/sweepdata/ref/rabies_wide.fa"))
    ap.add_argument("--out", type=Path, default=Path("results/abundance_sweep/divergence_pairs.json"))
    args = ap.parse_args()

    fa = fetch(args.term, args.retmax, args.fasta)
    seqs = {k: v for k, v in read_fasta(fa).items()
            if len(v) >= args.min_len
            and sum(v.count(c) for c in "NRYKMSWBDHV") / len(v) <= args.max_n}
    print(f"{len(seqs)} usable genomes", flush=True)

    keys = sorted(seqs)
    t0 = time.perf_counter()
    pairs = []
    for a, b in itertools.combinations(keys, 2):
        d = edlib.align(seqs[a], seqs[b], task="distance")["editDistance"]
        pairs.append((d / max(len(seqs[a]), len(seqs[b])), a, b))
    pairs.sort()
    print(f"{len(pairs)} pairs in {time.perf_counter() - t0:.0f}s; "
          f"range {pairs[0][0]:.4f} to {pairs[-1][0]:.4f}", flush=True)

    # Histogram of what the deposited data actually offers, so an empty band is
    # reported as empty rather than silently filled.
    edges = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10, 0.15, 1.0]
    hist = [{"from": lo, "to": hi,
             "pairs": sum(1 for d, _, _ in pairs if lo <= d < hi)}
            for lo, hi in zip(edges, edges[1:])]

    chosen = []
    for t in [float(x) for x in args.targets.split(",")]:
        d, a, b = min(pairs, key=lambda p: abs(p[0] - t))
        chosen.append({"target": t, "divergence": round(d, 4), "a": a, "b": b,
                       "len_a": len(seqs[a]), "len_b": len(seqs[b])})
        print(json.dumps(chosen[-1]), flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "term": args.term, "retmax": args.retmax,
        "genomes_usable": len(seqs), "pairs_compared": len(pairs),
        "histogram": hist, "chosen": chosen,
        "accessions": keys,
    }, indent=2))


if __name__ == "__main__":
    main()
