"""What do the residue's partial matches actually match? Inventory first, alignments second.

THE ASK AND WHY IT CHANGED SHAPE. Imo asked for pairwise alignments of residue reads that
match part of a reference but were not identified, with the scientific name of what they
matched, the whole read, the reference's matching nucleotides and enough reference past the
match to equal the read's length - to eyeball what the shape really is. His stated suspicion
was that we are looking at nothing: "you're comparing apples with pairs".

The rows this project already holds cannot be rendered as alignments, because the pilot's
BLAST payload carries no subject coordinates (qseqid, sseqid, pident, length, qlen, evalue,
bitscore and nothing else), so there is no way to say WHERE in a reference the match sat. A
re-run with the alignment fields is a paid launch and is on a decision card. This script is
what can be answered without it, and it turns out to be most of the question:

  1. THE SUBJECT INVENTORY NEEDS NO ALIGNMENT AT ALL. The accession is in the row, and one
     batched esummary call turns 300 accessions into 300 organism names and lengths. What a
     read matched is often the whole answer: a 107 nt wastewater read sharing 47 bases at
     83% with a chimpanzee BAC clone needs no picture.
  2. THE ALIGNMENT CAN BE REBUILT LOCALLY FOR SMALL REFERENCES. Fetch the reference, align
     the read to it here, and the coordinates BLAST did not report are recovered. This is
     capped by reference length, because fetching whole bacterial or fish chromosomes to
     place one 107 nt read is both slow and rude to NCBI. Over the cap the row keeps its
     organism and length and says plainly that it has no alignment.

EVERY ROW CARRIES ITS E-VALUE, AND THAT IS THE POINT. BLAST reports down to E=10 by default,
so "matched part of a reference" in this population mostly means a chance match: the median
E-value here is a fraction of one. This repo has already had to withdraw a reading of these
same rows ("the residue aligns ~36-48 nt at ~89% identity") for exactly that reason. So the
inventory is sorted and summarised BY E-value band, and no figure in it is quoted without
one. A rendering that let a reader forget the E-value would reproduce the error it exists to
show.

THE CONTROL. Chance matches and real relatives are not told apart by eye, so the inventory
reports the significant band (E < 1e-5) beside the chance band from the same run, with the
same columns. If the organism mix is the same in both, the subjects are telling us about the
database's composition rather than about the reads - which is itself the answer to Imo's
question, in the direction he suspected.

WHAT IS NOT CONTROLLED. These are mate 1 only, because the pilot query set is mate 1 only; a
pair-level statement cannot be made from it, and per CLAUDE.md must not be inferred. The
reference lengths and organism names are NCBI's own annotation, and a mislabelled record
would travel straight through - which is a reason to read the alignments and not only the
names.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import statistics as st
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: E-value bands. The boundary is this project's own significance cutoff, used so the two
#: bands here mean the same thing as the ones in results/subtract_known/REPORT.md.
SIGNIFICANT = 1e-5


def accession(sseqid: str) -> str:
    """BLAST's `gi|2452632002|ref|NC_071921.1|` down to `NC_071921.1`."""
    parts = sseqid.split("|")
    if sseqid.startswith("gi|") and len(parts) >= 4:
        return parts[3]
    return parts[-1] if parts[-1] else parts[-2] if len(parts) > 1 else sseqid


def eutils(endpoint: str, pause: float, **params) -> bytes:
    """One e-utilities call, with the courtesy pause NCBI asks of unkeyed clients.

    No contact details are sent: the queries are public SRA sequence and public accessions,
    and this project's user identity is not NCBI's business.
    """
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    time.sleep(pause)
    with urllib.request.urlopen(url, timeout=120) as h:
        return h.read()


def summaries(accs: list[str], pause: float, batch: int = 150) -> dict[str, dict]:
    """{accession: {organism, length, title}} from batched esummary calls."""
    out: dict[str, dict] = {}
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        raw = eutils("esummary.fcgi", pause, db="nuccore", id=",".join(chunk),
                     retmode="json")
        res = json.loads(raw)["result"]
        for uid in res.get("uids", []):
            r = res[uid]
            cap = r.get("caption", "")
            rec = {"organism": r.get("organism", ""), "length": r.get("slen"),
                   "title": r.get("title", ""), "uid": uid}
            out[cap] = rec
            if r.get("accessionversion"):
                out[r["accessionversion"]] = rec
        print(f"  esummary {min(i + batch, len(accs))}/{len(accs)}", flush=True)
    return out


def fetch_fasta(acc: str, pause: float) -> str | None:
    raw = eutils("efetch.fcgi", pause, db="nuccore", id=acc, rettype="fasta",
                 retmode="text")
    lines = raw.decode("utf-8", "replace").splitlines()
    seq = "".join(l.strip() for l in lines if l and not l.startswith(">"))
    return seq.upper() or None


def load_reads(names: set[str], r1: Path) -> dict[str, str]:
    """The named mate-1 reads, straight out of the FASTQ.

    Taken by name rather than by regenerating the query set, because the query builder is not
    stable across versions - a local regeneration once shared 6 of 500 read names with the
    rows it was supposed to match, and that mismatch is the kind that produces a plausible
    file about the wrong reads.
    """
    want = {n.rsplit("/", 1)[0] for n in names}
    out: dict[str, str] = {}
    with gzip.open(r1, "rt") as h:
        for i, line in enumerate(h):
            if i % 4 == 0:
                cur = line[1:].split()[0]
            elif i % 4 == 1 and cur in want:
                out[cur] = line.strip().upper()
                if len(out) == len(want):
                    break
    return out


def render(read: str, ref: str, r_st: int, r_en: int, strand: int, q_st: int,
           q_en: int) -> tuple[str, str, str]:
    """(read line, reference line, match bar) anchored so column i faces column i.

    The reference is extended PAST the match to the read's length, which is what Imo asked
    for and is the whole point: a real relative keeps matching into the padding and a chance
    match stops dead at the edge of the reported alignment. The padding is therefore the
    control built into the picture.
    """
    import mappy
    ref_seq = mappy.revcomp(ref) if strand == -1 else ref
    if strand == -1:
        r_st, r_en = len(ref) - r_en, len(ref) - r_st
    # Anchor the reference so its matched block sits under the read's matched block, then
    # pad both ends out to the read's length.
    start = r_st - q_st
    seg = ""
    for p in range(start, start + len(read)):
        seg += ref_seq[p] if 0 <= p < len(ref_seq) else "-"
    bar = "".join("|" if a == b and a != "-" else " " for a, b in zip(read, seg))
    return read, seg, bar


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rows", type=Path, action="append", required=True,
                    help="BLAST row JSONL; give it more than once to pool blocks")
    ap.add_argument("--reads-r1", type=Path,
                    default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--min-aln", type=int, default=30, help="Imo's floor: 30 matched bases")
    ap.add_argument("--max-aln-frac", type=float, default=0.8,
                    help="partial means the match does not cover the read")
    ap.add_argument("--max-ref-length", type=int, default=1_000_000,
                    help="references longer than this are inventoried but not fetched; "
                         "placing one 107 nt read in a 38 Mb fish genome is not worth the "
                         "download to NCBI or to us")
    ap.add_argument("--max-alignments", type=int, default=400)
    ap.add_argument("--pause", type=float, default=0.40,
                    help="NCBI asks unkeyed clients for no more than three requests a second")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/partialsubj"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/partial_match_subjects")
    args = ap.parse_args()

    t0 = time.perf_counter()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in args.rows:
        if not p.exists():
            raise SystemExit(f"{p} is missing")
        with open(p) as h:
            rows += [json.loads(l) for l in h if l.strip()]
    part = [r for r in rows if r.get("hit") and r.get("aln_len", 0) >= args.min_aln
            and r.get("aln_frac", 1.0) < args.max_aln_frac]
    if not part:
        raise SystemExit(f"no rows pass aln_len>={args.min_aln} and "
                         f"aln_frac<{args.max_aln_frac} in {len(rows)} rows")
    mates = Counter(r["read"].rsplit("/", 1)[-1] for r in part)
    print(f"{len(part)} partial matches of {len(rows)} rows; mates {dict(mates)}", flush=True)

    accs = sorted({accession(r["sseqid"]) for r in part})
    print(f"{len(accs)} distinct subjects; esummary", flush=True)
    meta = summaries(accs, args.pause)

    inventory = []
    for r in part:
        a = accession(r["sseqid"])
        m = meta.get(a, {})
        inventory.append({
            "read": r["read"], "class": r.get("class", "?"),
            "accession": a, "organism": m.get("organism", "?"),
            "ref_length": m.get("length"), "title": m.get("title", ""),
            "pident": r["pident"], "aln_len": r["aln_len"], "qlen": r["qlen"],
            "aln_frac": round(r["aln_frac"], 4), "evalue": r["evalue"],
            "band": "significant" if r["evalue"] < SIGNIFICANT else "chance"})

    def band_summary(band: str, cls: str | None = None) -> dict:
        v = [x for x in inventory if x["band"] == band
             and (cls is None or x["class"] == cls)]
        if not v:
            return {"n": 0}
        return {"n": len(v),
                "median_evalue": st.median([x["evalue"] for x in v]),
                "median_pident": round(st.median([x["pident"] for x in v]), 2),
                "median_aln_len": st.median([x["aln_len"] for x in v]),
                "median_aln_frac": round(st.median([x["aln_frac"] for x in v]), 4),
                "top_organisms": Counter(x["organism"] for x in v).most_common(12)}

    bands = {b: band_summary(b) for b in ("significant", "chance")}
    # THE CONTROL, where the rows carry one. The block files score residue reads beside
    # identified ones, so the same partial-match question can be asked of material we KNOW
    # is known. If the two classes give the same organism mix and the same median E-value,
    # the subjects describe the database rather than the reads.
    classes = sorted({x["class"] for x in inventory})
    by_class = {c: {b: band_summary(b, c) for b in ("significant", "chance")}
                for c in classes} if len(classes) > 1 else None
    for b, d in bands.items():
        print(f"{b}: n={d.get('n')} median E={d.get('median_evalue')} "
              f"pident={d.get('median_pident')} aln={d.get('median_aln_len')}", flush=True)
        for org, n in d.get("top_organisms", [])[:6]:
            print(f"    {n:4d}  {org}", flush=True)

    # ---------------------------------------------------- the alignments, where affordable
    reads = load_reads({x["read"] for x in inventory}, args.reads_r1)
    missing = [x["read"] for x in inventory if x["read"].rsplit("/", 1)[0] not in reads]
    if len(missing) > 0.1 * len(inventory):
        raise SystemExit(
            f"{len(missing)} of {len(inventory)} reads are not in {args.reads_r1}. The rows "
            f"and the FASTQ are not the same library or the names were parsed wrong; a file "
            f"rendered from this would be about the wrong reads.")

    import mappy
    small = [x for x in inventory
             if x["ref_length"] and x["ref_length"] <= args.max_ref_length]
    small.sort(key=lambda x: (x["band"] != "significant", x["evalue"]))
    small = small[:args.max_alignments]
    print(f"{len(small)} of {len(inventory)} have a reference at or under "
          f"{args.max_ref_length} nt and will be aligned", flush=True)

    cache: dict[str, str] = {}
    aligned, failed = [], []
    fa = open(args.out / f"partial_matches_s{args.seed}.fa", "w")
    aln = open(args.out / f"partial_matches_s{args.seed}.aln.txt", "w")
    aln.write(
        "Residue reads whose best core_nt hit covers part of the read, rendered against the\n"
        "reference they hit plus reference padding out to the read's length.\n"
        "READ THE E-VALUE FIRST: BLAST reports to E=10 by default and most of this\n"
        "population is a chance match. The padding is the control - a real relative keeps\n"
        "matching past the reported block, a chance match stops at its edge.\n"
        "Mate 1 only (the query set is mate 1 only), so nothing here is a pair statement.\n\n")
    for x in small:
        seq = reads.get(x["read"].rsplit("/", 1)[0])
        if not seq:
            continue
        if x["accession"] not in cache:
            try:
                s = fetch_fasta(x["accession"], args.pause)
            except Exception as exc:                      # noqa: BLE001 - network, logged
                failed.append({"accession": x["accession"], "error": str(exc)[:120]})
                cache[x["accession"]] = ""
            else:
                cache[x["accession"]] = s or ""
        ref = cache[x["accession"]]
        if not ref:
            continue
        # A SENSITIVITY LADDER, AND IT DOUBLES AS A READOUT. minimap2's short-read preset
        # seeds on exact 21-mers, and a 39-base match at 87% identity often contains none,
        # so the default preset reproduced only 1 of 25 of BLAST's hits here. Rather than
        # record that as a failure, the tiers are walked from the default down to a very
        # sensitive seeding, and the tier that first succeeds is kept beside the row. It is
        # a continuous measure of how weak the match is - "BLAST sees it and a standard
        # short-read aligner cannot" is a statement about the hit, not about the tool.
        best, tier = None, None
        for label, kw in (("sr", {"preset": "sr"}),
                          ("sr_k11", {"preset": "sr", "k": 11, "w": 5}),
                          ("sr_k9", {"preset": "sr", "k": 9, "w": 3})):
            a = mappy.Aligner(seq=ref, **kw)
            for h in a.map(seq):
                if best is None or h.mlen > best.mlen:
                    best = h
            if best is not None:
                tier = label
                break
        if best is None:
            failed.append({"read": x["read"], "accession": x["accession"],
                           "error": "no local alignment even at k=9; BLAST's hit is below "
                                    "what a short-read aligner can seed, which is itself "
                                    "the informative result"})
            continue
        rd, seg, bar = render(seq, ref, best.r_st, best.r_en, best.strand,
                              best.q_st, best.q_en)
        rec = dict(x, local_mlen=best.mlen, local_q_st=best.q_st, local_q_en=best.q_en,
                   local_r_st=best.r_st, local_r_en=best.r_en, local_strand=best.strand,
                   read_len_fastq=len(seq), seed_tier=tier,
                   padded_identity=round(sum(1 for c in bar if c == "|") / len(bar), 4))
        aligned.append(rec)
        head = (f"{x['read']} [{x['class']}] E={x['evalue']:g} {x['band']} "
                f"pident={x['pident']} "
                f"blast_aln={x['aln_len']}/{x['qlen']} | {x['accession']} "
                f"{x['organism']} ({x['ref_length']} nt)")
        fa.write(f">{head} READ\n{rd}\n>{head} REF_PADDED\n{seg}\n")
        aln.write(f"{head}\n")
        aln.write(f"  local match {best.q_st}-{best.q_en} of read, strand "
                  f"{'+' if best.strand == 1 else '-'}, seeded at {tier}, "
                  f"identity over the padded window {rec['padded_identity']:.3f}\n")
        aln.write(f"  read {rd}\n       {bar}\n  ref  {seg}\n\n")
    fa.close()
    aln.close()

    with open(args.out / f"inventory_s{args.seed}.tsv", "w") as h:
        cols = ["read", "class", "accession", "organism", "ref_length", "pident",
                "aln_len", "qlen", "aln_frac", "evalue", "band", "title"]
        h.write("\t".join(cols) + "\n")
        for x in sorted(inventory, key=lambda y: (y["band"] != "significant", y["evalue"])):
            h.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")

    pad = [x["padded_identity"] for x in aligned]
    out = {"script": "experiments/partial_match_subjects.py", **provenance(__file__),
           "parameters": {k: ([str(x) for x in v] if isinstance(v, list)
                              else str(v) if isinstance(v, Path) else v)
                          for k, v in vars(args).items()},
           "rows_read": len(rows), "partial_matches": len(part),
           "mates": dict(mates), "distinct_subjects": len(accs),
           "bands": bands, "by_class": by_class,
           "aligned": len(aligned), "fetch_failures": failed[:40],
           "unreproducible_locally": len(failed),
           "seed_tier": dict(Counter(x["seed_tier"] for x in aligned)),
           "padded_identity": ({"median": round(st.median(pad), 4),
                                "mean": round(st.mean(pad), 4),
                                "frac_above_0.8": round(
                                    sum(1 for p in pad if p > 0.8) / len(pad), 4)}
                               if pad else None),
           "readout": "identity across a window the width of the read, centred on the "
                      "reported match. A real relative holds its identity into the "
                      "padding; a chance match decays to background outside its block",
           "bias": "mate 1 only, so no pair statement is available. Organism names are "
                   "NCBI's annotation and travel through unchecked. References over the "
                   "length cap are inventoried without an alignment, and the cap "
                   "correlates with organism - large genomes are cellular - so the "
                   "aligned subset leans viral relative to the inventory",
           "wall_seconds": time.perf_counter() - t0,
           "python": platform.python_version(), "platform": platform.platform()}
    (args.out / f"partial_match_subjects_s{args.seed}.json").write_text(
        json.dumps(out, indent=2))
    print(f"written {args.out}/ ({len(aligned)} alignments, "
          f"{len(inventory)} inventory rows)", flush=True)


if __name__ == "__main__":
    main()
