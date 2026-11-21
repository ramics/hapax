"""Render partially-matching residue reads against the reference they matched, to be read by eye.

WHY THIS EXISTS, AND WHOSE QUESTION IT ANSWERS. The cascade reports that a residue read's
best relative in core_nt is typically ~40 of 107 bases at ~88% identity, and calls that
distant homology. Imo's objection on 2026-11-21 is that it might be nothing of the kind:
if the read runs off the end of what the reference record actually covers - a spliced
transcript against a genomic subject, a partial record, an operon boundary - then the
alignment stops for a reason that has nothing to do with divergence, and "partial match"
is an artefact of the comparison rather than a property of the sequence.

That objection has a signature and this script is built to expose it. **The boundary story
predicts near-perfect identity INSIDE the match**, because nothing is diverging: the same
organism's sequence simply runs out of comparable reference. Distant homology predicts the
opposite - identity in the 80-95% band, degrading smoothly, which is what the spike ladder
measures at 10-20% divergence. So the discriminating field is identity within the aligned
block, and the three fields that say which boundary was hit are carried beside it:

- `to_record_end`: bases from the end of the match to the end of the subject record. Zero or
  near it means the reference ran out, which is the artefact.
- `molecule`: whether the subject is a genomic record or an mRNA/cDNA one, parsed from its
  title. Splicing can only truncate an alignment against a genomic subject, so a read whose
  match stops mid-record against an mRNA cannot be explained that way.
- `q_stop`: where in the read the alignment stops. A population of artefacts clusters at a
  consistent read position; biology does not.

WHAT THE FILE SHOWS AND HOW IT IS ANCHORED. For each read the reference segment is cut so
that **position i of the reference record corresponds to position i of the read**: the match
itself, extended backwards by however many read bases precede it and forwards by however
many follow, so both records are the read's length and can be compared column by column by
eye. The padding past the match is the point of the exercise - it is the sequence the read
was supposed to resemble and does not - so it is taken from the subject record itself
(`efetch` with a coordinate range, reverse-complemented when the match is on the minus
strand), not from the alignment, which by definition stops at the match.

SELECTION, AND THE CENSORING IT HAS TO UNDO. The reads are drawn from those the cascade did
NOT identify, which means the operating point has already removed everything at >=90%
identity over >=40 bases - exactly the population the objection is about. A sample drawn
naively from the unidentified pile therefore under-represents the artefact by construction.
So `--high-identity-share` deliberately reserves part of the sample for partial matches at
>=95% identity, and the manifest records how the sample was composed. This file is for
looking at, not for a rate: every proportion quoted from it would be a proportion of a
stratified sample.
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

#: Words in a subject title that say the record is a transcript rather than genomic DNA.
#: Only a genomic subject can truncate an alignment at a splice junction, so this is what
#: makes that explanation falsifiable per row rather than plausible in general.
MRNA_WORDS = ("mrna", "cdna", "transcript", "complete cds", "partial cds", "est")
GENOMIC_WORDS = ("genome", "chromosome", "genomic", "scaffold", "contig", "plasmid",
                 "complete sequence", "dna")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def molecule_of(title: str) -> str:
    t = title.lower()
    for w in MRNA_WORDS:
        if w in t:
            return "mrna"
    for w in GENOMIC_WORDS:
        if w in t:
            return "genomic"
    return "unknown"


def blast_submit(seqs: dict[str, str], database: str) -> str:
    fasta = "".join(f">{n}\n{s}\n" for n, s in seqs.items())
    data = {"CMD": "Put", "PROGRAM": "blastn", "TASK": "blastn", "DATABASE": database,
            "QUERY": fasta, "HITLIST_SIZE": "1", "FORMAT_TYPE": "JSON2_S",
            "TOOL": "hapax"}
    req = urllib.request.Request(BLAST_URL,
                                 data=urllib.parse.urlencode(data).encode(),
                                 headers={"User-Agent": "hapax"})
    with urllib.request.urlopen(req, timeout=180) as r:
        body = r.read().decode(errors="replace")
    m = re.search(r"RID = (\S+)", body)
    if not m:
        raise RuntimeError("no RID in submission reply")
    return m.group(1)


def blast_fetch(rid: str, poll: int, max_wait: int) -> str:
    waited = 0
    while waited < max_wait:
        time.sleep(poll)
        waited += poll
        url = f"{BLAST_URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}&TOOL=hapax"
        with urllib.request.urlopen(url, timeout=120) as r:
            info = r.read().decode(errors="replace")
        if "Status=READY" in info:
            break
        if "Status=FAILED" in info or "Status=UNKNOWN" in info:
            raise RuntimeError(f"RID {rid} failed at NCBI")
    else:
        raise TimeoutError(f"RID {rid} not ready after {max_wait}s")
    url = f"{BLAST_URL}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}&TOOL=hapax"
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read().decode(errors="replace")


def parse_json2(text: str) -> dict[str, dict]:
    """Best HSP per query, with everything needed to cut the reference and to judge it."""
    doc = json.loads(text)
    out: dict[str, dict] = {}
    reports = doc.get("BlastOutput2", [])
    if isinstance(reports, dict):
        reports = [reports]
    for rep in reports:
        res = rep.get("report", {}).get("results", {})
        search = res.get("search") or {}
        qid = search.get("query_title") or search.get("query_id")
        hits = search.get("hits") or []
        if not qid or not hits:
            if qid:
                out[qid] = {}
            continue
        hit = hits[0]
        hsp = sorted(hit.get("hsps", []), key=lambda h: -h.get("bit_score", 0))[0]
        desc = (hit.get("description") or [{}])[0]
        aln = hsp.get("align_len") or 0
        ident = hsp.get("identity") or 0
        out[qid] = {
            "accession": desc.get("accession") or desc.get("id", ""),
            "title": desc.get("title", ""),
            "sciname": desc.get("sciname", ""),
            "taxid": desc.get("taxid"),
            "subject_len": hit.get("len"),
            "bitscore": hsp.get("bit_score"), "evalue": hsp.get("evalue"),
            "align_len": aln, "identities": ident,
            "pident": (100.0 * ident / aln) if aln else 0.0,
            "q_start": hsp.get("query_from"), "q_stop": hsp.get("query_to"),
            "s_start": hsp.get("hit_from"), "s_stop": hsp.get("hit_to"),
            "q_strand": hsp.get("query_strand"), "s_strand": hsp.get("hit_strand"),
            "qseq": hsp.get("qseq", ""), "hseq": hsp.get("hseq", ""),
            "midline": hsp.get("midline", ""),
        }
    return out


def efetch_range(acc: str, start: int, stop: int, pause: float) -> str:
    """One subject slice, 1-based inclusive, forward orientation.

    Ranges rather than whole records because a subject in core_nt can be a chromosome: the
    point is 107 bases of it, and fetching the record to take 107 bases of it would move
    gigabytes for a file meant to be read by eye.
    """
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "fasta",
                                "retmode": "text", "seq_start": max(1, start),
                                "seq_stop": max(1, stop), "tool": "hapax"})
    time.sleep(pause)
    with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=120) as r:
        text = r.read().decode(errors="replace")
    lines = [l.strip() for l in text.splitlines()]
    return "".join(l for l in lines if l and not l.startswith(">")).upper()


def reference_window(hit: dict, qlen: int, pause: float) -> dict:
    """Cut the subject so that column i of the reference is column i of the read.

    The match anchors it: whatever precedes the match in the READ is taken from the subject
    before the match, and whatever follows it is taken after, which is the sequence the read
    failed to resemble and the whole reason for the file. On the minus strand the subject's
    coordinates run the other way, so the window is mirrored and reverse-complemented -
    getting that wrong would silently produce a reference that cannot match anything, which
    is exactly how a strand bug fooled this repo once already.
    """
    pre = (hit["q_start"] or 1) - 1
    post = qlen - (hit["q_stop"] or qlen)
    minus = (hit.get("s_strand") == "minus") or (hit["s_start"] > hit["s_stop"])
    slen = hit.get("subject_len") or 0
    if minus:
        hi, lo = max(hit["s_start"], hit["s_stop"]), min(hit["s_start"], hit["s_stop"])
        start, stop = lo - post, hi + pre
        seq = revcomp(efetch_range(hit["accession"], start, stop, pause))
        to_record_end = lo - 1  # the match's far edge, in the read's direction of travel
    else:
        start, stop = hit["s_start"] - pre, hit["s_stop"] + post
        seq = efetch_range(hit["accession"], start, stop, pause)
        to_record_end = (slen - hit["s_stop"]) if slen else None
    clipped = start < 1 or (slen and stop > slen)
    return {"ref_seq": seq, "ref_from": max(1, start), "ref_to": stop,
            "strand": "minus" if minus else "plus", "to_record_end": to_record_end,
            "window_clipped_by_record": bool(clipped),
            "pad_before": pre, "pad_after": post}


#: The strata, by the only field that turned out to decide what a partial match IS.
#: Measured on the pilot's 500 residue reads: the median partial-but-unidentified match has
#: E = 0.24 and 70% of them sit above 0.05, so most of what the cascade was calling a distant
#: relative is the best chance match a 107 nt query finds in 319 GB. Demand significance and
#: the shape inverts - real hits align over 96% of the read. So the file is stratified by
#: significance first and identity second, and the three strata are what have to be compared
#: by eye: the real partials (the only place the splice/boundary question can live), the
#: full-length divergent strains, and the chance matches that outnumber both.
STRATA = ("significant_partial", "significant_full", "chance")


def stratum_of(row: dict, sig: float, chance: float) -> str | None:
    if row["evalue"] < sig:
        return "significant_partial" if row.get("aln_frac", 1) < 0.95 \
            else "significant_full"
    if row["evalue"] >= chance:
        return "chance"
    return None  # the marginal band is deliberately not sampled; nothing to see by eye


def pick(rows: list[dict], reads: dict[str, str], n: int, min_match: int, sig: float,
         chance: float, rng: random.Random) -> tuple[list[str], dict[str, int]]:
    """Equal quotas per stratum, not a random draw from the pile.

    A proportional sample would be ~90% chance matches, which is the true composition and
    the wrong file: the populations worth looking at are the small ones. Quotas make the
    comparison possible and make the file useless as a rate, which the docstring says and
    the manifest records.
    """
    cand: dict[str, list[dict]] = {s: [] for s in STRATA}
    for r in rows:
        if not (r.get("hit") and r.get("aln_len", 0) >= min_match
                and r.get("read") in reads):
            continue
        s = stratum_of(r, sig, chance)
        if s:
            cand[s].append(r)
    for s in STRATA:
        rng.shuffle(cand[s])
    per = max(1, n // len(STRATA))
    chosen: list[dict] = []
    for s in STRATA:  # a stratum short of its quota gives its places to the others
        chosen += cand[s][:per]
    spare = n - len(chosen)
    if spare > 0:
        pool = [r for s in STRATA for r in cand[s][per:]]
        rng.shuffle(pool)
        chosen += pool[:spare]
    rng.shuffle(chosen)
    available = {s: len(cand[s]) for s in STRATA}
    return [r["read"] for r in chosen], available


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reads", type=Path,
                    default=Path("/tmp/claude-0/wk/bigindex_queries/pilot.fa"),
                    help="the FASTA the pilot searched; names are <class>|<read>")
    ap.add_argument("--rows", type=Path,
                    default=ROOT / "results/bigindex/rows_pilot_core_nt_i-0b42267a.jsonl",
                    help="one row per query from the staged core_nt pilot. Used only to "
                         "choose which reads are worth a remote search - the alignment "
                         "itself is re-derived here because those rows carry no coordinates")
    ap.add_argument("--n", type=int, default=45)
    ap.add_argument("--min-match", type=int, default=30,
                    help="bases of alignment to qualify, per Imo's ask")
    ap.add_argument("--significant", type=float, default=1e-5,
                    help="E-value below which a match is treated as real. 1e-5 rather than "
                         "BLAST's reporting default of 10, because the default is what let "
                         "a 27-base chance match be quoted as the residue's closest relative")
    ap.add_argument("--chance-above", type=float, default=0.05,
                    help="E-value at or above which a match is treated as chance")
    ap.add_argument("--database", default="core_nt")
    ap.add_argument("--batch", type=int, default=15)
    ap.add_argument("--between-batches", type=int, default=10)
    ap.add_argument("--poll", type=int, default=20)
    ap.add_argument("--max-wait", type=int, default=1800)
    ap.add_argument("--efetch-pause", type=float, default=0.4,
                    help="NCBI asks for no more than three requests a second without a key")
    ap.add_argument("--out", type=Path, default=ROOT / "results/partial_eyeball")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.out.mkdir(parents=True, exist_ok=True)
    reads_raw = read_fasta(args.reads)
    reads = {n.split("|", 1)[1] if "|" in n else n: s for n, s in reads_raw.items()}
    # The reads have to be the ones the rows were computed from. The instance regenerated
    # its query set between launches, so a local copy of an earlier set shares only a
    # handful of names with these rows - which is how the first run of this script picked
    # four reads instead of forty-five, silently. Hence the explicit check.
    def _check(rows_: list[dict]) -> None:
        if not rows_:
            return
        hit = sum(1 for r in rows_ if r.get("read") in reads)
        if hit < 0.9 * len(rows_):
            raise SystemExit(
                f"only {hit} of {len(rows_)} rows name a read present in {args.reads}. "
                f"These rows were computed from a different query set - fetch the one the "
                f"instance actually searched rather than sampling whatever is local.")
    rows = [json.loads(l) for l in open(args.rows)] if args.rows.exists() else []
    rng = random.Random(args.seed)
    _check(rows)
    names, available = pick(rows, reads, args.n, args.min_match, args.significant,
                            args.chance_above, rng)
    if not names:
        raise SystemExit("no qualifying reads; check --rows and --reads")
    print(f"{len(names)} reads selected of {len(rows)} pilot rows; available per "
          f"stratum {available}", flush=True)

    hits: dict[str, dict] = {}
    errors: list[str] = []
    for i in range(0, len(names), args.batch):
        chunk = {n: reads[n] for n in names[i:i + args.batch]}
        try:
            rid = blast_submit(chunk, args.database)
            hits.update(parse_json2(blast_fetch(rid, args.poll, args.max_wait)))
            print(f"batch {i // args.batch}: {len(chunk)} queries, "
                  f"{sum(1 for n in chunk if hits.get(n))} with a hit", flush=True)
        except Exception as e:  # a lost batch is dropped and named, never scored as a miss
            errors.append(f"batch {i // args.batch}: {e!r}")
            print(f"batch {i // args.batch} FAILED: {e!r}", flush=True)
        time.sleep(args.between_batches)

    records, table = [], []
    for name in names:
        h = hits.get(name)
        if not h or not h.get("accession"):
            continue
        read = reads[name]
        try:
            win = reference_window(h, len(read), args.efetch_pause)
        except Exception as e:
            errors.append(f"{name} efetch: {e!r}")
            continue
        pilot_row = next((r for r in rows if r["read"] == name), {})
        row = {"read": name, "read_len": len(read),
               "stratum": stratum_of({**pilot_row, "evalue": h["evalue"],
                                      "aln_frac": (h["align_len"] / len(read))},
                                     args.significant, args.chance_above) or "marginal",
               "molecule": molecule_of(h["title"]), **{k: v for k, v in h.items()
                                                       if k not in ("qseq", "hseq",
                                                                    "midline")},
               **{k: v for k, v in win.items() if k != "ref_seq"}}
        table.append(row)
        sci = h["sciname"] or "unknown organism"
        records.append((
            f">{name} read len={len(read)} match={h['align_len']}nt "
            f"pident={h['pident']:.1f} qpos={h['q_start']}-{h['q_stop']} "
            f"stratum={row['stratum']}\n{read}\n"
            f">{h['accession']}:{win['ref_from']}-{win['ref_to']} {win['strand']} "
            f"ref_for={name} [{sci}] molecule={row['molecule']} "
            f"to_record_end={win['to_record_end']} "
            f"pad={win['pad_before']}+{win['pad_after']} | {h['title'][:120]}\n"
            f"{win['ref_seq']}\n"), )
        records[-1] = records[-1][0]

    fa = args.out / f"partial_eyeball_s{args.seed}.fa"
    fa.write_text("".join(records))

    # The rendered alignment is the thing that can actually be eyeballed: query over
    # midline over subject, exactly as BLAST computed it, with the padded reference beneath
    # so the part that does NOT match is visible rather than merely absent.
    aln_lines = []
    for row in table:
        h = hits[row["read"]]
        win = {k: row[k] for k in ("ref_from", "ref_to", "strand", "to_record_end")}
        aln_lines += [
            f"=== {row['read']}  ({row['stratum']})",
            f"    subject {h['accession']}  [{h['sciname'] or 'unknown organism'}]  "
            f"{row['molecule']}",
            f"    {h['title'][:140]}",
            f"    match {h['align_len']} nt at {h['pident']:.1f}%, read {h['q_start']}-"
            f"{h['q_stop']} of {row['read_len']}, subject {h['s_start']}-{h['s_stop']} "
            f"({win['strand']}), {win['to_record_end']} nt of record past the match",
            f"  read {h['qseq']}",
            f"       {h['midline']}",
            f"  subj {h['hseq']}",
            "",
        ]
    (args.out / f"partial_eyeball_s{args.seed}.aln.txt").write_text("\n".join(aln_lines)
                                                                    + "\n")

    with open(args.out / f"partial_eyeball_s{args.seed}.tsv", "w") as h:
        cols = ["read", "stratum", "sciname", "accession", "molecule", "pident",
                "align_len", "q_start", "q_stop", "read_len", "s_start", "s_stop",
                "strand", "subject_len", "to_record_end", "window_clipped_by_record",
                "evalue", "bitscore", "title"]
        h.write("\t".join(cols) + "\n")
        for row in table:
            h.write("\t".join(str(row.get(c, "")).replace("\t", " ") for c in cols) + "\n")

    manifest = {"script": "experiments/partial_alignment_eyeball.py",
                **provenance(__file__),
                "parameters": {k: (str(v) if isinstance(v, Path) else v)
                               for k, v in vars(args).items()},
                "selected": len(names), "rendered": len(table),
                "available_per_stratum_in_pilot": available,
                "not_a_rate": "quotas per stratum, so no proportion in this file is an "
                              "estimate of anything in the library",
                "dropped": len(names) - len(table), "errors": errors,
                "strata": {s: sum(1 for r in table if r["stratum"] == s)
                           for s in (*STRATA, "marginal")},
                "wall_seconds": time.perf_counter() - t0,
                "python": platform.python_version(), "platform": platform.platform()}
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps(manifest, indent=2))
    print(f"{len(table)} alignments written to {fa}", flush=True)


if __name__ == "__main__":
    main()
