"""Split this sample's own assembly into what can be identified and what cannot.

WHY THIS EXISTS - a correction, not an addition. The first cut of the cascade
(`experiments/subtract_known.py`) subtracted the sample's own assembly wholesale, and that
was the most efficient filter measured by a wide margin: 35 MB of contigs removed 83.6% of
pairs alone. But Imo's three-tier architecture (2026-11-20) makes clear that this is the
wrong subtrahend:

  1. remove everything we can assemble AND identify - or identify at all, assembled or not;
  2. check everything we can assemble but CANNOT identify, for novelty;
  3. check the residual reads for genuinely new biology.

Subtracting the assembly wholesale swallows tier 2 entirely: a novel virus abundant enough
to assemble would be deleted as background before anyone looked at it. So the subtrahend
has to be defined by IDENTIFICATION, not by assembly. This script does that split, and the
cascade is re-run against the identified part only, with the unidentified contigs routed
to tier 2 rather than into the bin.

HOW A CONTIG COUNTS AS IDENTIFIED. Two independent routes, because they fail differently
and a contig only needs one of them:

  - nucleotide: Kraken2 assigns it a taxon at the stated confidence. Cheap, and the only
    route that can place non-coding sequence such as rRNA;
  - protein: DIAMOND blastx finds a hit at or above the stated bitscore. Reaches further
    down the divergence axis, which is the whole reason the project kept a protein tier.

Both are swept rather than fixed, and the partition is written at the operating point so
the downstream cascade is reproducible. A contig is "unidentifiable" only in the sense of
"not identified by these indexes at these thresholds", which is the same one-sided
statement the rest of this work makes about absence; the sweep is how that limit is
reported rather than hidden.

The tier-2 population is also described here rather than just counted: length and SPAdes
k-mer coverage from the contig headers, because tier 2's own success measure is
sensitivity on contigs where existing classifiers actually operate, and this project
measured that boundary at ~500 nt and ~2x coverage (PR #5, the abundance sweep).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import dir_bytes, git_sha, kraken2

ROOT = Path(__file__).resolve().parents[1]

# SPAdes contig names carry length and k-mer coverage: NODE_1_length_12710_cov_727.441011
HEADER = re.compile(r"length_(\d+)_cov_([0-9.]+)")


def read_fasta(path: Path):
    name, buf = None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="unused: the input is a fixed assembly, but the repo requires it")
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--k2-db", type=Path, default=ROOT / "data/k2_pluspf16")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--dmnd", type=Path, nargs="+",
                    default=[ROOT / "data/prot/sprot.dmnd",
                             ROOT / "data/prot/uniref50_10pct.dmnd"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--k2-confidence", type=float, nargs="+", default=[0.0, 0.05, 0.1, 0.2])
    ap.add_argument("--operating-confidence", type=float, default=0.1)
    ap.add_argument("--operating-bitscore", type=float, default=50.0)
    ap.add_argument("--silva-mmi", type=Path,
                    default=Path("/tmp/claude-0/wk/subtract/silva.mmi"),
                    help="minimap2 index of SILVA SSU+LSU built by subtract_known.py")
    ap.add_argument("--rrna-min-mlen", type=int, default=100,
                    help="matching bases against SILVA for a CONTIG to count as "
                         "ribosomal. Stricter than the read-level floor because a contig "
                         "is long: a 40 nt match inside a 3 kb contig is a shared motif, "
                         "not a ribosomal contig")
    ap.add_argument("--min-len", type=int, default=200,
                    help="contigs shorter than this are not routed to tier 2: no "
                         "classifier in the 2024 benchmark operates there, so calling "
                         "them unidentifiable would overstate the tier-2 population")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/contigs"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/subtract_known")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()

    contigs = dict(read_fasta(args.contigs))
    lens = {k: len(v) for k, v in contigs.items()}
    cov = {}
    for k in contigs:
        m = HEADER.search(k)
        if m:
            cov[k] = float(m.group(2))

    # ---- ribosomal identification, kept separate from everything else
    #
    # Why its own route. The biggest step on the subtraction curve is this sample's own
    # assembly, and PR #14 found evidence that ribosomal sequence is leaking into the
    # protein side of the same trick (ribosomal reads came out LESS peptide-absent than
    # other reads, which cannot happen against an index of real proteins). On the
    # nucleotide side the equivalent worry is that the assembly's big step is partly
    # ribosomal RNA that SILVA missed - legitimate removal, but it belongs to the rRNA
    # tier, not to "sequence unknown to databases". Splitting the contigs this way lets
    # the self-assembly step be reported with ribosomal contigs excluded as well as
    # included, which is the check that decides whether that step survives.
    timings: dict[str, float] = {}
    rrna_hit: set[str] = set()
    rrna_evidence: dict[str, dict] = {}
    if args.silva_mmi and Path(args.silva_mmi).exists():
        import mappy
        t0 = time.perf_counter()
        sidx = mappy.Aligner(str(args.silva_mmi), preset="sr", n_threads=args.threads)
        for k, v in contigs.items():
            best, bestref = 0, None
            for a in sidx.map(v):
                if a.mlen > best:
                    best, bestref = a.mlen, a.ctg
            if best >= args.rrna_min_mlen:
                rrna_hit.add(k)
                # The ribosomal route has no identity to report, so it reports what it
                # actually used: matching bases and the reference. Writing a percent
                # identity here would invent a number of a different kind.
                rrna_evidence[k] = {"matching_bases": best, "reference": bestref,
                                    "of_contig": round(best / max(1, len(v)), 4)}
        timings["silva_contigs"] = time.perf_counter() - t0
        del sidx

    # ---- nucleotide identification
    fa = args.work / "contigs.fa"
    with open(fa, "w") as h:
        for k, v in contigs.items():
            h.write(f">{k}\n{v}\n")
    k2_sweep, k2_hit, k2_evidence = {}, set(), {}
    for conf in args.k2_confidence:
        out = args.work / f"k2_{conf}.out"
        t0 = time.perf_counter()
        cmd = [str(args.classify_bin), "-H", str(args.k2_db / "hash.k2d"),
               "-t", str(args.k2_db / "taxo.k2d"), "-o", str(args.k2_db / "opts.k2d"),
               "-M", "-p", str(args.threads), "-T", str(conf), "-O", str(out), str(fa)]
        subprocess.run(cmd, check=True, capture_output=True)
        timings[f"kraken2_conf{conf}"] = time.perf_counter() - t0
        hit = set()
        with open(out) as h:
            for line in h:
                _st, name, taxid, _l, _hits = line.rstrip("\n").split("\t")
                if taxid != "0":
                    hit.add(name)
                    if conf == args.operating_confidence:
                        # Kraken2 gives a taxon and a confidence, never an identity, and
                        # the two are not interchangeable: a taxid is a name, a percent
                        # identity is a distance. Recorded as what it is.
                        k2_evidence[name] = {"taxid": taxid, "confidence": conf,
                                             "db": args.k2_db.name}
        k2_sweep[f"{conf:.2f}"] = len(hit) / len(contigs)
        if conf == args.operating_confidence:
            k2_hit = hit
        out.unlink(missing_ok=True)

    # ---- protein identification
    # evidence[contig] is the removal's own justification, carried out of this script
    # so that nothing downstream has to re-derive it from a boolean. SPEC 2.1(i), 5.
    prot_sweep, prot_hit, evidence = {}, {}, {}
    for db in args.dmnd:
        if not db.exists():
            continue
        tsv = args.work / f"{db.stem}.tsv"
        t0 = time.perf_counter()
        # SPEC 2.1(i) IS A REQUIREMENT ON THE OUTPUT, AND THIS CALL IS WHERE IT WAS FAILED.
        # These columns were already being computed and then thrown away: the loop below
        # kept the maximum bitscore and nothing else, so a contig removed at 66.5% identity
        # and one removed at 99.9% were indistinguishable downstream. results/holdout_cascade/
        # is that failure with a known answer in it - a genome withheld from every index,
        # removed in tier 1 at 66.5% over 620 aa against this tier's threshold of 50, with no
        # trace of the distance in any output. The extra columns (alignment length, E-value,
        # query start and end, query length) cost nothing: DIAMOND has them in hand. qstart
        # and qend are here because identity alone flatters towards "familiar" and coverage
        # is the control that removes most of that bias, so a distance without a coverage is
        # not reportable. The identification decision itself is UNCHANGED - still the maximum
        # bitscore against --operating-bitscore - because widening what is recorded must not
        # quietly widen what is subtracted.
        subprocess.run([str(args.diamond_bin), "blastx", "--db", str(db), "--query", str(fa),
                        "--out", str(tsv), "--threads", str(args.threads), "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "bitscore", "length", "evalue",
                        "qstart", "qend", "qlen", "--max-target-seqs", "1",
                        "--evalue", "10", "--sensitive", "--quiet"], check=True)
        timings[f"diamond_{db.stem}"] = time.perf_counter() - t0
        best: dict[str, float] = {}
        with open(tsv) as h:
            for line in h:
                f = line.rstrip("\n").split("\t")
                q, subj, pid, bits = f[0], f[1], float(f[2]), float(f[3])
                aln, ev = int(f[4]), float(f[5])
                qs, qe, qlen = int(f[6]), int(f[7]), int(f[8])
                if bits > best.get(q, 0):
                    best[q] = bits
                # The evidence table is keyed on the same comparison the decision uses, so
                # the row kept is the row that caused (or would have caused) the removal.
                prev = evidence.get(q)
                if prev is None or bits > prev["bitscore"]:
                    evidence[q] = {
                        "contig": q, "contig_nt": qlen, "index": db.stem, "level": "protein",
                        "subject": subj, "pident": pid, "aln_aa": aln, "evalue": ev,
                        "bitscore": bits,
                        # One HSP per subject at --max-target-seqs 1, so this is the best
                        # subject's coverage of the query and not a union over subjects. The
                        # any-subject union is a different field with a different meaning and
                        # is measured in results/holdout_cascade/, not conflated here.
                        "query_coverage": (round((abs(qe - qs) + 1) / qlen, 4)
                                           if qlen else None),
                    }
        prot_sweep[db.stem] = {str(b): sum(1 for v in best.values() if v >= b) / len(contigs)
                               for b in (30, 40, 50, 60, 80, 100)}
        prot_hit[db.stem] = {q for q, v in best.items() if v >= args.operating_bitscore}
        tsv.unlink(missing_ok=True)

    identified = set(k2_hit)
    for s in prot_hit.values():
        identified |= s
    identified |= rrna_hit
    identified_norrna = identified - rrna_hit
    unidentified = [k for k in contigs if k not in identified]
    tier2 = [k for k in unidentified if lens[k] >= args.min_len]

    id_fa = args.contigs.with_name(args.contigs.stem + "_identified.fa")
    idnr_fa = args.contigs.with_name(args.contigs.stem + "_identified_norrna.fa")
    rr_fa = args.contigs.with_name(args.contigs.stem + "_rrna.fa")
    un_fa = args.contigs.with_name(args.contigs.stem + "_unidentified.fa")
    with open(id_fa, "w") as h:
        for k in identified:
            h.write(f">{k}\n{contigs[k]}\n")
    with open(idnr_fa, "w") as h:
        for k in identified_norrna:
            h.write(f">{k}\n{contigs[k]}\n")
    with open(rr_fa, "w") as h:
        for k in rrna_hit:
            h.write(f">{k}\n{contigs[k]}\n")
    with open(un_fa, "w") as h:
        for k in tier2:
            h.write(f">{k}\n{contigs[k]}\n")

    def desc(keys):
        ls = [lens[k] for k in keys]
        cs = [cov[k] for k in keys if k in cov]
        return {
            "contigs": len(keys),
            "bases": sum(ls),
            "len_median": st.median(ls) if ls else None,
            "len_p90": sorted(ls)[int(0.9 * len(ls))] if ls else None,
            "len_max": max(ls) if ls else None,
            "cov_median": st.median(cs) if cs else None,
            "frac_over_500nt": sum(1 for x in ls if x >= 500) / len(ls) if ls else None,
            "frac_over_2x": sum(1 for x in cs if x >= 2) / len(cs) if cs else None,
            "frac_over_500nt_and_2x": (sum(1 for k in keys
                                           if lens[k] >= 500 and cov.get(k, 0) >= 2)
                                       / len(keys)) if keys else None,
        }

    # ---- the removal's own justification, one row per identified contig ------------------
    #
    # WHY THIS FILE EXISTS AND WHY IT IS A FILE. SPEC 5: there is one distance continuum and
    # the tier is a column on it, not a separate table. A contig removed here at family-level
    # distance is still reportable material, so the distance that justified the removal has to
    # leave this script. It did not before: the summary JSON held fractions and the splits held
    # sequence, and nothing held "this one was removed at 66.5% over 620 aa".
    #
    # ROUTES ARE NOT POOLED. A protein row carries an identity and a coverage; a Kraken2 row
    # carries a taxid and a confidence; a ribosomal row carries matching bases. Those are
    # three kinds of number and averaging or mixing them would be comparing two instruments,
    # so `route` is a column and a reader can group on it.
    #
    # THERE IS NO SCORE AND NO ORDERING IN THIS FILE. Ordering by distance is ordering by
    # database proximity, which is the one thing this project does not do: measured on the
    # holdout, identity ordering puts the known-novel genome at the 93.46th percentile of the
    # removed pile. The file is emitted in input order and the consumer joins coverage to it.
    ev_path = args.out / "removal_evidence.tsv"
    cols = ["contig", "contig_nt", "route", "index", "subject_or_taxon", "pident",
            "aln_aa", "evalue", "bitscore", "query_coverage", "matching_bases",
            "confidence", "tier"]
    n_rows = 0
    with open(ev_path, "w") as h:
        h.write("\t".join(cols) + "\n")
        for k in contigs:
            if k not in identified:
                continue
            tier = "rrna" if k in rrna_hit else "identified_norrna"
            rows = []
            e = evidence.get(k)
            if e is not None and e["bitscore"] >= args.operating_bitscore:
                rows.append({"route": "protein", "index": e["index"],
                             "subject_or_taxon": e["subject"], "pident": e["pident"],
                             "aln_aa": e["aln_aa"], "evalue": e["evalue"],
                             "bitscore": e["bitscore"],
                             "query_coverage": e["query_coverage"]})
            if k in k2_evidence:
                g = k2_evidence[k]
                rows.append({"route": "kraken2", "index": g["db"],
                             "subject_or_taxon": g["taxid"],
                             "confidence": g["confidence"]})
            if k in rrna_evidence:
                r = rrna_evidence[k]
                rows.append({"route": "rrna", "index": "silva",
                             "subject_or_taxon": r["reference"],
                             "matching_bases": r["matching_bases"],
                             "query_coverage": r["of_contig"]})
            if not rows:
                # Identified by a route that recorded nothing is a bug, not an empty row, so
                # it is written as an explicit unknown rather than omitted. An absent row and
                # a row that says nothing was recorded are different facts.
                rows.append({"route": "unrecorded", "index": "",
                             "subject_or_taxon": "ROUTE DID NOT RECORD ITS EVIDENCE"})
            for r in rows:
                r.update({"contig": k, "contig_nt": lens[k], "tier": tier})
                h.write("\t".join("" if r.get(c) is None else str(r.get(c, ""))
                                   for c in cols) + "\n")
                n_rows += 1

    result = {
        "contigs_file": str(args.contigs),
        "contigs": len(contigs),
        "bases": sum(lens.values()),
        "operating_point": {"kraken2_confidence": args.operating_confidence,
                            "diamond_bitscore": args.operating_bitscore,
                            "kraken2_db": args.k2_db.name,
                            "protein_dbs": [d.stem for d in args.dmnd if d.exists()]},
        "kraken2_confidence_sweep": k2_sweep,
        "protein_bitscore_sweep": prot_sweep,
        "identified_frac": len(identified) / len(contigs),
        "identified_bases_frac": sum(lens[k] for k in identified) / sum(lens.values()),
        "route_nucleotide_only": len(k2_hit - set().union(*prot_hit.values())) / len(contigs)
        if prot_hit else None,
        "route_protein_only": (len(set().union(*prot_hit.values()) - k2_hit) / len(contigs)
                               if prot_hit else None),
        "all": desc(list(contigs)),
        "identified": desc(list(identified)),
        "rrna_contigs": desc(list(rrna_hit)),
        "rrna_contig_frac": len(rrna_hit) / len(contigs),
        "rrna_bases_frac": (sum(lens[k] for k in rrna_hit) / sum(lens.values())
                            if rrna_hit else 0.0),
        "identified_excluding_rrna": desc(list(identified_norrna)),
        "unidentified": desc(unidentified),
        "tier2_population": desc(tier2),
        "tier2_min_len": args.min_len,
        "removal_evidence": {
            "file": str(ev_path),
            "rows": n_rows,
            "a_row_is": ("one route's justification for removing one contig. There is no "
                         "score on this table and no ordering in this file"),
            "routes_are_not_pooled": ("protein carries identity and coverage, kraken2 a "
                                      "taxid and a confidence, rrna matching bases; three "
                                      "kinds of number, so group on `route`"),
            "upper_bound_warning": ("every identity is measured against the indexes this "
                                    "run was given, so it is an upper bound on distance"),
        },
        "files": {"identified": str(id_fa),
                  "identified_excluding_rrna": str(idnr_fa),
                  "rrna_contigs": str(rr_fa),
                  "unidentified_tier2": str(un_fa)},
        "index_bytes": {"kraken2": dir_bytes(args.k2_db / "hash.k2d"),
                        **{d.stem: dir_bytes(d) for d in args.dmnd if d.exists()}},
        "timings_seconds": timings,
    }
    (args.out / "contig_identification.json").write_text(json.dumps(result, indent=2))
    (args.out / "contig_identification_manifest.json").write_text(json.dumps({
        "script": "experiments/classify_contigs.py",
        "git_sha": git_sha(),
        # The convention's own test, which this manifest was failing. A -dirty SHA names
        # HEAD and not the code that ran, so two runs can carry the same SHA and different
        # code; the content hash pins the bytes whatever the tree looked like. Taken HERE,
        # when the manifest is written, not at startup - see CLAUDE.md.
        "script_sha256": hashlib.sha256(
            Path(__file__).resolve().read_bytes()).hexdigest(),
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }, indent=2))
    fa.unlink(missing_ok=True)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("files", "index_bytes", "timings_seconds")}, indent=2))


if __name__ == "__main__":
    main()
