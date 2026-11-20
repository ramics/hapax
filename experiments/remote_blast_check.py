"""Ground-truth check: BLAST a stratified sample against the real NCBI nt and nr.

WHY. Imo, 2026-11-20: "why are we not running an actual blast at least once to confirm?"
Every number in `results/subtract_known/REPORT.md` rests on local indexes this project
openly calls thin - a capped Kraken2 hash, Swiss-Prot, a 10% UniRef50 subset - so the word
"unidentifiable" currently means "not identified by what fits on this disk". One honest
check against the full public databases decides whether that word is doing any work at
all, and it calibrates every other number in the report.

DESIGN. A stratified random sample across the classes the cascade assigns, because the
interesting question is per class rather than overall:

  rrna                 removed by SILVA
  public_nucleotide    removed by Kraken2 and not by SILVA
  selfasm_identified   removed by this sample's own identified contigs only
  residue_assembled    survived, but maps to one of this sample's UNIDENTIFIED contigs:
                       present here, often abundant, unplaceable - tier 2's read-level
                       population
  residue_missing      survived and maps to no contig of this sample at all - tier 3, by
                       construction rare AND unsequenced
  tier2_contigs        unidentifiable contigs of >=500 nt and >=2x coverage: the subset
                       where classifiers actually operate, and the class that should
                       discriminate best because a contig carries far more signal than
                       107 nt

Labels are recomputed here from the same references and the same thresholds as the
cascade, so this script is self-contained and reproducible from a seed rather than
depending on another run's bookkeeping.

WHAT IT CAN AND CANNOT SETTLE. A few hundred sequences gives percentage-level resolution
and no better, so every proportion is reported with a Wilson 95% interval and the
report should quote the interval, not the point estimate. BLAST against nt/nr is also not
a definition of "known": nt is a nucleotide database with its own gaps, megablast is
tuned for close matches, and a true novel family would be expected to miss nt and hit nr
only weakly - which is exactly the pattern this script is looking for, and the reason
both programs are run on the same sequences.

ETIQUETTE AND PRIVACY. Submissions go in modest batches with polling waits rather than a
flood. The sequences are public SRA data, so nothing sensitive leaves the container; no
contact details are sent, only a generic tool name.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs
from subtract_known import _mm_index, dir_bytes, git_sha, kraken2, map_filter

ROOT = Path(__file__).resolve().parents[1]
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def wilson(k: int, n: int) -> tuple[float, float]:
    """95% interval for a proportion. Wilson rather than normal because these counts are
    small and often near 0 or 1, where the normal interval leaves the unit interval."""
    if n == 0:
        return (0.0, 0.0)
    z = 1.959964
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - s) / d, (c + s) / d)


def blast_submit(seqs: dict[str, str], program: str, database: str,
                 extra: dict | None = None) -> str:
    query = "".join(f">{k}\n{v}\n" for k, v in seqs.items())
    data = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": query,
            "HITLIST_SIZE": "5", "FORMAT_TYPE": "JSON2_S", "TOOL": "hapax"}
    if program == "blastn":
        data["MEGABLAST"] = "on" if (extra or {}).get("megablast", True) else "off"
    data.update({k: v for k, v in (extra or {}).items() if k != "megablast"})
    req = urllib.request.Request(BLAST_URL,
                                 data=urllib.parse.urlencode(data).encode(),
                                 headers={"User-Agent": "hapax-subtraction-check"})
    with urllib.request.urlopen(req, timeout=180) as r:
        body = r.read().decode("utf-8", "replace")
    m = re.search(r"RID = (\S+)", body)
    if not m:
        raise RuntimeError("no RID returned; first 300 chars: " + body[:300])
    return m.group(1)


def blast_fetch(rid: str, poll: int, max_wait: int) -> str:
    waited = 0
    while waited < max_wait:
        time.sleep(poll)
        waited += poll
        url = (f"{BLAST_URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}&TOOL=hapax")
        with urllib.request.urlopen(url, timeout=120) as r:
            info = r.read().decode("utf-8", "replace")
        if "Status=READY" in info:
            break
        if "Status=UNKNOWN" in info or "Status=FAILED" in info:
            raise RuntimeError(f"RID {rid} failed: {info[:200]}")
    # JSON2_S rather than Tabular: the URL API returns an empty body for Tabular on a
    # finished job (verified 2026-11-20), and the JSON carries the subject's taxid and
    # scientific name, which is what makes "what are the hits actually" answerable.
    url = f"{BLAST_URL}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}&TOOL=hapax"
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read().decode("utf-8", "replace")


def parse_json2(text: str) -> dict[str, dict]:
    """Best hit per query from the JSON2_S report, with identity and query coverage.

    Percent identity is computed from identities over alignment length rather than taken
    from a field, and coverage is the aligned query span over the query length: an
    alignment covering 30 of 107 nt at 95% identity is a different statement from one
    covering 100, and a report that printed only identity would hide that.
    """
    best: dict[str, dict] = {}
    data = json.loads(text[text.index("{"):])
    for block in data.get("BlastOutput2", []):
        res = block.get("report", {}).get("results", {}).get("search", {})
        q = res.get("query_title") or res.get("query_id")
        qlen = res.get("query_len") or 0
        for hit in res.get("hits", []):
            d = (hit.get("description") or [{}])[0]
            for hsp in hit.get("hsps", []):
                alen = hsp.get("align_len") or 0
                rec = {"subject": d.get("accession"), "title": d.get("title"),
                       "taxid": d.get("taxid"), "sciname": d.get("sciname"),
                       "pident": round(100 * hsp.get("identity", 0) / alen, 2) if alen else None,
                       "aln_len": alen,
                       "query_cov": (round(abs(hsp.get("query_to", 0) - hsp.get("query_from", 0) + 1)
                                           / qlen, 4) if qlen else None),
                       "evalue": hsp.get("evalue"), "bitscore": hsp.get("bit_score")}
                if q not in best or (rec["bitscore"] or 0) > (best[q]["bitscore"] or 0):
                    best[q] = rec
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=20_000,
                    help="pairs read to draw the stratified sample from")
    ap.add_argument("--per-class", type=int, default=60)
    ap.add_argument("--batch", type=int, default=15, help="sequences per BLAST job")
    ap.add_argument("--poll", type=int, default=20)
    ap.add_argument("--max-wait", type=int, default=1200)
    ap.add_argument("--between-batches", type=int, default=10)
    ap.add_argument("--programs", nargs="+", default=["blastn", "blastx"])
    ap.add_argument("--megablast", choices=["on", "off"], default="off",
                    help="off runs the slower, far more sensitive blastn task. The "
                         "question here is whether a relative exists 'however far', and "
                         "megablast is tuned for close matches, so it would answer a "
                         "different question and flatter the cascade")
    ap.add_argument("--only-classes", nargs="*", default=None,
                    help="restrict BLAST to these classes. Used to run the expensive "
                         "blastx-vs-nr pass on the residue and tier-2 classes only: the "
                         "removed classes are already answered by blastn, and a protein "
                         "search on them would spend hours of a public service to "
                         "confirm what nt already says")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--silva-mmi", type=Path, default=Path("/tmp/claude-0/wk/subtract/silva.mmi"))
    ap.add_argument("--silva", type=Path, default=ROOT / "data/silva/silva_dna.fa")
    ap.add_argument("--strong-pident", type=float, default=95.0)
    ap.add_argument("--strong-qcov", type=float, default=0.8)
    ap.add_argument("--identified", type=Path,
                    default=ROOT / "data/geo/sample_contigs_identified.fa")
    ap.add_argument("--unidentified", type=Path,
                    default=ROOT / "data/geo/sample_contigs_unidentified.fa")
    ap.add_argument("--k2-db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--k2-confidence", type=float, default=0.1)
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--contig-min-len", type=int, default=500)
    ap.add_argument("--contig-min-cov", type=float, default=2.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/blastchk"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/blast_check")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    t_start = time.perf_counter()

    # ---- label a block of reads with the same tiers and thresholds as the cascade
    rows = list(read_pairs(args.r1, args.r2, args.pairs, 14, (args.seed - 1) * 1_500_000))
    seqs = {r[0]: (r[1], r[2]) for r in rows}
    fq1, fq2 = args.work / "s1.fq", args.work / "s2.fq"
    for path, i in ((fq1, 1), (fq2, 2)):
        with open(path, "w") as h:
            for name, s1, s2, q1, q2 in rows:
                s, q = (s1, q1) if i == 1 else (s2, q2)
                h.write(f"@{name}\n{s}\n+\n{q}\n")

    silva = set(map_filter(_mm_index(args.silva, args.silva_mmi, args.threads),
                           seqs, args.min_mlen))
    k2_calls, _w, _r = kraken2(args.k2_db, fq1, fq2, args.classify_bin, args.threads,
                               args.work / "k2.out", args.k2_confidence)
    k2 = {k for k, tx in k2_calls.items() if tx != 0}
    idx_id = _mm_index(args.identified, args.work / "identified.mmi", args.threads)
    ident = set(map_filter(idx_id, seqs, args.min_mlen))
    del idx_id
    idx_un = _mm_index(args.unidentified, args.work / "unidentified.mmi", args.threads)
    unident = set(map_filter(idx_un, seqs, args.min_mlen))
    del idx_un

    classes: dict[str, list[str]] = {k: [] for k in
                                     ("rrna", "public_nucleotide", "selfasm_identified",
                                      "residue_assembled", "residue_missing")}
    for name in seqs:
        if name in silva:
            classes["rrna"].append(name)
        elif name in k2:
            classes["public_nucleotide"].append(name)
        elif name in ident:
            classes["selfasm_identified"].append(name)
        elif name in unident:
            classes["residue_assembled"].append(name)
        else:
            classes["residue_missing"].append(name)

    sample: dict[str, dict[str, str]] = {}
    for cls, names in classes.items():
        rng.shuffle(names)
        # Mate 1 only: a pair would double the BLAST load for a near-duplicate question,
        # and the per-class proportions are per molecule either way.
        sample[cls] = {f"{cls}|{n}": seqs[n][0] for n in names[:args.per_class]}

    # ---- tier 2 contigs where classifiers can operate
    tier2 = {}
    name = None
    buf: list[str] = []
    cand = []
    with open(args.unidentified) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    cand.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        cand.append((name, "".join(buf)))
    keep = [(n, s) for n, s in cand if len(s) >= args.contig_min_len
            and _cov(n) >= args.contig_min_cov]
    rng.shuffle(keep)
    for n, s in keep[:args.per_class]:
        tier2[f"tier2_contigs|{n}"] = s
    sample["tier2_contigs"] = tier2

    # ---- BLAST, in modest batches
    results: dict[str, dict] = {}
    query_rows: list[dict] = []
    for program in args.programs:
        db = "nt" if program == "blastn" else "nr"
        got: dict[str, dict] = {}
        queued: list[tuple[str, str]] = []
        for cls, d in sample.items():
            if args.only_classes and cls not in args.only_classes:
                continue
            queued += list(d.items())
        for i in range(0, len(queued), args.batch):
            chunk = dict(queued[i:i + args.batch])
            try:
                rid = blast_submit(chunk, program, db,
                                   {"megablast": args.megablast == "on"})
                text = blast_fetch(rid, args.poll, args.max_wait)
                got.update(parse_json2(text))
            except Exception as e:  # a failed batch is recorded, not fatal
                results.setdefault("errors", []).append(
                    f"{program} batch {i // args.batch}: {type(e).__name__}: {e}")
            time.sleep(args.between_batches)
        per_class = {}
        for cls, d in sample.items():
            if args.only_classes and cls not in args.only_classes:
                continue
            n = len(d)
            hits = [got[k] for k in d if k in got]
            lo, hi = wilson(len(hits), n)
            per_class[cls] = {
                "n": n,
                "any_hit": len(hits),
                "any_hit_frac": len(hits) / n if n else None,
                "any_hit_ci95": [round(lo, 4), round(hi, 4)],
                "pident": _dist(h["pident"] for h in hits if h["pident"] is not None),
                "query_cov": _dist(h["query_cov"] for h in hits
                                   if h["query_cov"] is not None),
                "top_names": Counter(h.get("sciname") for h in hits).most_common(8),
                "aln_len": _dist(h["aln_len"] for h in hits),
                "bitscore": _dist(h["bitscore"] for h in hits),
                "top_titles": Counter((h.get("title") or "")[:90]
                                      for h in hits).most_common(6),
                # Hit-or-no-hit is the wrong readout against a 1.1-Tbase database: almost
                # everything hits something. What separates the classes is how much of the
                # read is covered at how high an identity, so the joint criterion is
                # reported alongside, and the product identity x coverage as a continuous
                # version of it.
                "strong_hit_frac": (sum(
                    1 for h in hits
                    if (h["pident"] or 0) >= args.strong_pident
                    and (h["query_cov"] or 0) >= args.strong_qcov) / n) if n else None,
                "pident_x_qcov": _dist((h["pident"] or 0) / 100 * (h["query_cov"] or 0)
                                       for h in hits),
            }
            for k in d:
                h = got.get(k, {})
                query_rows.append({
                    "program": program, "class": cls, "query": k.split("|", 1)[1],
                    "hit": k in got, "pident": h.get("pident"),
                    "query_cov": h.get("query_cov"), "bitscore": h.get("bitscore"),
                    "aln_len": h.get("aln_len"), "sciname": h.get("sciname"),
                    "title": (h.get("title") or "")[:120],
                })
        results[program] = {"database": db, "per_class": per_class}

    out = {
        "seed": args.seed,
        "pairs_labelled": len(seqs),
        "class_sizes_in_block": {k: len(v) for k, v in classes.items()},
        "class_fracs_in_block": {k: len(v) / len(seqs) for k, v in classes.items()},
        "tier2_contig_candidates": len(keep),
        "sampled_per_class": {k: len(v) for k, v in sample.items()},
        "blast": results,
        "wall_seconds": time.perf_counter() - t_start,
    }
    (args.out / f"blast_check_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    with open(args.out / f"blast_rows_s{args.seed}.jsonl", "w") as h:
        for r in query_rows:
            h.write(json.dumps(r) + "\n")
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/remote_blast_check.py",
        "git_sha": git_sha(),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "note": "queries are public SRA sequence; no contact details are sent to NCBI",
    }, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k != "blast"}, indent=2))
    for program, r in results.items():
        if program == "errors":
            print("ERRORS:", r)
            continue
        print(f"\n== {program} vs {r['database']}")
        for cls, c in r["per_class"].items():
            print(f"{cls:>20}  n={c['n']:>3}  any hit {c['any_hit_frac'] * 100:5.1f}% "
                  f"[{c['any_hit_ci95'][0] * 100:.1f}-{c['any_hit_ci95'][1] * 100:.1f}]  "
                  f"pident median {(c['pident'] or {}).get('median')}  "
                  f"aln median {(c['aln_len'] or {}).get('median')}")


def _cov(header: str) -> float:
    m = re.search(r"cov_([0-9.]+)", header)
    return float(m.group(1)) if m else 0.0


def _dist(xs):
    xs = sorted(xs)
    if not xs:
        return None
    return {"n": len(xs), "min": xs[0], "p25": xs[len(xs) // 4],
            "median": xs[len(xs) // 2], "p75": xs[3 * len(xs) // 4], "max": xs[-1]}


if __name__ == "__main__":
    main()
