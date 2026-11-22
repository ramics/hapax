"""Gene-level workup of the contigs that sit far from any published protein.

WHY THIS EXISTS. `known_strain_divergence.py` gives one verdict per contig: the best
DIAMOND hit over the whole thing. That is the right instrument for the question it asks
(is this the reference organism or a drifted version of it) and the wrong one for the
question that follows, which is what the divergent rows actually ARE. A single best hit
over a 3.9 kb contig is one gene's worth of evidence standing in for a genome.

THE CORRECTION THIS SCRIPT CARRIES, and it changes how the output must be read. A row at
30-45% amino acid identity is NOT a divergent strain of the taxon it is named after.
Strain-level variation lives at 90-99% protein identity. At 30-45% the named subject is
the nearest neighbour in a sparse database, which is a statement about the database, not
a taxonomy. So every row here reports `divergence_class` from the identity itself, and
the word "strain" appears only above 90%.

WHAT IT MEASURES, each with the control that makes it readable.

  1. ORFs, and a per-ORF search rather than a per-contig one. Three things fall out that
     a whole-contig hit cannot give: which gene carries the divergence, whether different
     ORFs point at unrelated subjects (the signature of a mis-joined contig, not of a
     virus), and what the proteins are for. The control is the same pass over contigs the
     same run placed at >=95% identity, length-matched, where the ORFs must agree with
     each other and with one subject.

  2. Pfam-A domains over those ORFs. Profile reach is a MEASURED NEGATIVE at read level
     (results/residue_structure/, SPEC §4.2 item 5 withdrawn: weaker than DIAMOND at
     every rung) and a measured positive on contigs, where an ORF runs to hundreds of
     residues: 23.1% of new residue contigs carry a Pfam domain against a length-matched
     control (results/residue_contig_annotation/). This script is on the contig side of
     that line, which is why the profile arm is here at all.

  3. True depth from mapping every read back, not the assembler's `cov_` header. The
     header is a k-mer coverage estimate from the assembly graph; it is not read support
     and the two differ. Evenness is reported beside the mean, because a contig whose
     coverage is a spike over 10% of its length is a repeat or a mis-join wearing a
     genome's length.

  4. A chimaera reading, BECAUSE NO CONTIGUITY NUMBER SHIPS WITHOUT ONE (CLAUDE.md). The
     evidence is disagreement between ORFs over which subject they point at, plus the
     coverage discontinuity at the junction. Its false-positive direction is stated: ORF
     disagreement also fires on a genuine virus whose genes have different nearest
     neighbours, which at this distance is common, so the field is reported as a flag to
     read beside the coverage, never as a verdict.

WHAT IT CANNOT DO. It cannot say what an ORF does when nothing in Pfam or the protein
databases is close enough to it. Past roughly 20-30% identity, sequence search is at the
end of its reach, and the honest instrument is structure, which needs a GPU and is not
attempted here. An ORF that gets no domain and no subject is reported as exactly that.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The identity bands, stated once and used everywhere, so that the word "strain" cannot
# drift onto a row that does not earn it. The boundaries are conventions, not
# measurements, and are printed in the output beside every use.
BANDS = [
    (97.0, "reference or near-reference"),
    (90.0, "divergent strain of the named taxon"),
    (70.0, "same genus or thereabouts, not a strain"),
    (40.0, "family-level or looser; the subject is a nearest neighbour, not an identification"),
    (0.0, "below any band a name survives"),
]

# A DIAMOND alignment shorter than this says nothing about identity at these distances:
# this repo has measured the chance band for short queries and withdrawn a claim to it.
MIN_ALN_AA = 60
SIGNIFICANT_EVALUE = 1e-5


def band(pident: float) -> str:
    for lo, name in BANDS:
        if pident >= lo:
            return name
    return BANDS[-1][1]


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def read_fasta(path: Path) -> dict:
    """Minimal FASTA reader; the assembly is plain text and this avoids a dependency."""
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


def pick_contigs(strain_json: Path, lo: float, hi: float) -> dict:
    """Contigs whose own alignment identity falls in [lo, hi), with their taxon row.

    Selection is on the CONTIG's identity, not the taxon's median: a taxon can carry one
    divergent contig and four that match the reference, and pooling those would average
    away the row worth looking at.
    """
    d = json.loads(strain_json.read_text())
    out = {}
    for tax in d.get("taxa", []):
        for r in tax.get("rows", []):
            if r.get("aln_aa", 0) < MIN_ALN_AA:
                continue
            p = r.get("pident")
            if p is None or not (lo <= p < hi):
                continue
            out[r["contig"]] = {
                "taxon": tax["name"],
                "on_the_panel": tax.get("on_the_panel", False),
                "contig_nt": r["nt"],
                "whole_contig_subject": r["subject"],
                "whole_contig_pident": p,
                "whole_contig_aln_aa": r["aln_aa"],
                "whole_contig_evalue": r["evalue"],
                "divergence_class": band(p),
            }
    return out


def call_orfs(seqs: dict, min_aa: int) -> list:
    """pyrodigal in meta mode, the mode a metagenomic assembly requires.

    meta=True because a single-genome training pass needs a megabase of one organism's
    sequence and these contigs are kilobases of many. The same call produced the proteins
    indexed elsewhere in this repo, kept identical so the two are comparable.
    """
    import pyrodigal

    gf = pyrodigal.GeneFinder(meta=True)
    orfs = []
    for name, seq in seqs.items():
        for i, g in enumerate(gf.find_genes(seq.encode()), start=1):
            aa = g.translate().rstrip("*")
            if len(aa) < min_aa:
                continue
            orfs.append({
                "orf": f"{name}__orf{i}",
                "contig": name,
                "start": int(g.begin),
                "end": int(g.end),
                "strand": int(g.strand),
                "aa_len": len(aa),
                "partial": bool(g.partial_begin or g.partial_end),
                "aa": aa,
            })
    return orfs


def diamond_blastp(orfs: list, db: Path, work: Path, binpath: Path, threads: int,
                   mode: str, label: str) -> dict:
    """Best hit per ORF against one protein index. Returns {orf: hit}."""
    fa = work / f"orfs_{label}.faa"
    fa.write_text("".join(f">{o['orf']}\n{o['aa']}\n" for o in orfs))
    tsv = work / f"orfs_{label}.tsv"
    cmd = [str(binpath), "blastp", "-d", str(db), "-q", str(fa), "-o", str(tsv),
           "-p", str(threads), mode, "--quiet",
           "--evalue", "1e-3", "--max-target-seqs", "5", "--outfmt", "6",
           "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    best, per_orf = {}, collections.defaultdict(list)
    with tsv.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            row = {"subject": f[1], "pident": float(f[2]), "aln_aa": int(f[3]),
                   "evalue": float(f[4]), "bitscore": float(f[5]),
                   "title": f[6] if len(f) > 6 else "", "index": label}
            per_orf[f[0]].append(row)
    for q, rows in per_orf.items():
        rows.sort(key=lambda r: -r["bitscore"])
        best[q] = rows[0]
        best[q]["other_subjects"] = [r["subject"] for r in rows[1:4]]
    return {"best": best, "seconds": round(time.time() - t0, 1)}


def pfam_scan(orfs: list, hmm: Path, work: Path, threads: int) -> dict:
    """Pfam-A over the ORFs with pyhmmer, at the curators' own gathering thresholds.

    bit_cutoffs="gathering" rather than an E-value of our choosing: the GA cutoff is the
    family curator's statement of what belongs in that family, and substituting a
    threshold of our own would make every domain call this script's opinion.
    """
    import pyhmmer

    alpha = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=o["orf"].encode(),
                                       sequence=o["aa"]).digitize(alpha)
            for o in orfs]
    hits = collections.defaultdict(list)
    t0 = time.time()

    def txt(v) -> str:
        """pyhmmer hands these back as bytes on some versions and str on others."""
        if isinstance(v, (bytes, bytearray)):
            return v.decode()
        return "" if v is None else str(v)

    with pyhmmer.plan7.HMMFile(str(hmm)) as fh:
        for top in pyhmmer.hmmsearch(fh, seqs, cpus=threads, bit_cutoffs="gathering"):
            fam = txt(top.query.name)
            desc = txt(top.query.description)
            for h in top:
                if h.included:
                    hits[txt(h.name)].append(
                        {"pfam": fam, "description": desc,
                         "bitscore": round(h.score, 1), "evalue": h.evalue})
    return {"hits": dict(hits), "seconds": round(time.time() - t0, 1)}


def map_back(contig_fa: Path, r1: Path, r2: Path, trim_r2: int, work: Path,
             limit: int) -> dict:
    """Read support per contig: mean depth, evenness, and the breadth actually covered.

    The assembler's `cov_` header is a k-mer coverage estimate off the assembly graph. It
    is not read support, and a contig can carry a high header coverage with reads piled on
    a fraction of its length. Evenness is reported as the fraction of positions within a
    factor of two of the contig's own median, so a spike cannot read as a genome.

    Mate 2 is trimmed by the same amount every experiment in this repo trims it: it
    carries a non-template 5' prefix, so the two mates are 107 and 93 nt and are not
    interchangeable.
    """
    import gzip
    import mappy

    aln = mappy.Aligner(str(contig_fa), preset="sr", n_threads=4)
    if not aln:
        raise SystemExit("mappy could not index the contig file")
    lens = {n: len(aln.seq(n)) for n in aln.seq_names}
    depth = {n: [0] * lens[n] for n in lens}
    seen = 0

    def stream(path: Path, trim: int):
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    yield line.strip()[trim:]

    for path, trim in ((r1, 0), (r2, trim_r2)):
        for n, seq in enumerate(stream(path, trim)):
            if limit and n >= limit:
                break
            seen += 1
            for h in aln.map(seq):
                if h.ctg in depth:
                    for p in range(h.r_st, min(h.r_en, lens[h.ctg])):
                        depth[h.ctg][p] += 1
    out = {}
    for n, d in depth.items():
        covered = [x for x in d if x > 0]
        med = statistics.median(d) if d else 0
        even = (sum(1 for x in d if med and med / 2 <= x <= med * 2) / len(d)) if d else 0
        out[n] = {"mapped_mean_depth": round(sum(d) / len(d), 2) if d else 0,
                  "mapped_median_depth": med,
                  "breadth_covered": round(len(covered) / len(d), 4) if d else 0,
                  "evenness_within_2x_of_median": round(even, 4),
                  "max_depth": max(d) if d else 0}
    return {"per_contig": out, "mates_streamed": seen}


def workup(contigs: dict, seqs: dict, orfs: list, hits: dict, pfam: dict,
           depth: dict) -> list:
    """One row per contig, with the per-ORF evidence that produced it."""
    by_contig = collections.defaultdict(list)
    for o in orfs:
        by_contig[o["contig"]].append(o)
    rows = []
    for name, meta in contigs.items():
        mine = sorted(by_contig.get(name, []), key=lambda o: o["start"])
        orf_rows, subjects, named = [], [], 0
        for o in mine:
            h = hits.get(o["orf"])
            dom = pfam.get(o["orf"], [])
            if h:
                named += 1
                subjects.append(h["subject"])
            orf_rows.append({
                "orf": o["orf"], "start": o["start"], "end": o["end"],
                "strand": o["strand"], "aa_len": o["aa_len"], "partial": o["partial"],
                "best_subject": h["subject"] if h else None,
                "best_title": (h or {}).get("title", "")[:160],
                "pident": (h or {}).get("pident"),
                "aln_aa": (h or {}).get("aln_aa"),
                "evalue": (h or {}).get("evalue"),
                "significant": bool(h and h["evalue"] < SIGNIFICANT_EVALUE
                                    and h["aln_aa"] >= MIN_ALN_AA),
                "divergence_class": band(h["pident"]) if h else None,
                "pfam": [{"id": d["pfam"], "description": d["description"]}
                         for d in dom],
            })
        # Chimaera reading. ORFs pointing at unrelated subjects is the signature of a
        # mis-join AND of a genuine virus whose genes have different nearest neighbours,
        # which at this distance is ordinary. So it is a flag with its direction stated,
        # read beside coverage evenness, and never a verdict on its own.
        distinct = len({s for s in subjects})
        rows.append({
            **meta,
            "contig": name,
            "orfs": len(mine),
            "orfs_with_a_subject": named,
            "orfs_with_a_pfam_domain": sum(1 for r in orf_rows if r["pfam"]),
            "distinct_subjects_across_orfs": distinct,
            "coding_density": round(
                sum(abs(o["end"] - o["start"]) for o in mine) / max(1, len(seqs[name])), 3),
            "read_support": depth.get(name, {}),
            "subject_disagreement_flag": distinct > 1 and len(mine) > 1,
            "orf_rows": orf_rows,
        })
    rows.sort(key=lambda r: (r["whole_contig_pident"], -r["contig_nt"]))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--strains", type=Path, required=True,
                    help="known_strain_divergence_s*.json from the run directory")
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--trim-r2", type=int, default=16)
    ap.add_argument("--map-limit", type=int, default=0,
                    help="mates per file to stream, 0 for the whole library")
    ap.add_argument("--diverged-below", type=float, default=90.0)
    # The control arm. Contigs the same run placed at or above this identity are the
    # known-positive: their ORFs should agree with each other and with one subject, and
    # if the test's fields fire equally on them the fields are describing assembly rather
    # than divergence.
    ap.add_argument("--control-above", type=float, default=95.0)
    ap.add_argument("--min-orf-aa", type=int, default=50)
    ap.add_argument("--dmnd", type=Path, default=ROOT / "data/prot/viral.dmnd")
    ap.add_argument("--dmnd-sprot", type=Path, default=ROOT / "data/prot/sprot.dmnd")
    ap.add_argument("--pfam", type=Path, default=ROOT / "data/pfam/Pfam-A.hmm")
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-mode", default="--sensitive")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    for p in (a.strains, a.contigs, a.dmnd, a.pfam):
        if not p.exists():
            raise SystemExit(f"missing input: {p}")
    a.work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    test = pick_contigs(a.strains, 0.0, a.diverged_below)
    ctrl = pick_contigs(a.strains, a.control_above, 101.0)
    for v in test.values():
        v["arm"] = "divergent"
    for v in ctrl.values():
        v["arm"] = "control_near_reference"
    both = {**test, **ctrl}
    print(f"{len(test)} divergent contigs, {len(ctrl)} near-reference controls",
          file=sys.stderr)

    allseqs = read_fasta(a.contigs)
    seqs = {n: allseqs[n] for n in both if n in allseqs}
    missing = [n for n in both if n not in allseqs]
    sub = a.work / "workup_contigs.fa"
    sub.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))

    orfs = call_orfs(seqs, a.min_orf_aa)
    print(f"{len(orfs)} ORFs >= {a.min_orf_aa} aa", file=sys.stderr)

    viral = diamond_blastp(orfs, a.dmnd, a.work, a.diamond_bin, a.threads,
                           a.diamond_mode, "viral")
    hits = dict(viral["best"])
    sprot_sec = 0.0
    if a.dmnd_sprot.exists():
        sp = diamond_blastp(orfs, a.dmnd_sprot, a.work, a.diamond_bin, a.threads,
                            a.diamond_mode, "sprot")
        sprot_sec = sp["seconds"]
        # Best across the two indexes by bitscore, ties to the viral index, which is the
        # specific one. Recorded per hit as "index" so a row says where its name came from.
        for q, r in sp["best"].items():
            if q not in hits or r["bitscore"] > hits[q]["bitscore"]:
                hits[q] = r

    pf = pfam_scan(orfs, a.pfam, a.work, a.threads)
    depth = map_back(sub, a.r1, a.r2, a.trim_r2, a.work, a.map_limit)

    rows = workup(both, seqs, orfs, hits, pf["hits"], depth["per_contig"])
    div = [r for r in rows if r["arm"] == "divergent"]
    con = [r for r in rows if r["arm"] == "control_near_reference"]

    def summarise(rs):
        if not rs:
            return {}
        return {
            "contigs": len(rs),
            "median_orfs": statistics.median([r["orfs"] for r in rs]),
            "fraction_with_any_orf_subject":
                round(sum(1 for r in rs if r["orfs_with_a_subject"]) / len(rs), 4),
            "fraction_with_any_pfam_domain":
                round(sum(1 for r in rs if r["orfs_with_a_pfam_domain"]) / len(rs), 4),
            "fraction_flagged_subject_disagreement":
                round(sum(1 for r in rs if r["subject_disagreement_flag"]) / len(rs), 4),
            "median_coding_density":
                round(statistics.median([r["coding_density"] for r in rs]), 3),
            "median_mapped_depth": statistics.median(
                [r["read_support"].get("mapped_median_depth", 0) for r in rs]),
            "median_breadth_covered": round(statistics.median(
                [r["read_support"].get("breadth_covered", 0) for r in rs]), 4),
        }

    out = {
        "unit": "one row per contig; orf_rows is one row per predicted protein",
        "question": "what are the contigs that sit far from any published protein",
        "the_word_strain": (
            "reserved for 90-97% identity. A row at 30-45% is not a divergent strain of "
            "the taxon it is named after; the subject is the nearest neighbour in a "
            "sparse database. divergence_class states the band on every row."),
        "bands": [{"at_least_percent": lo, "reading": name} for lo, name in BANDS],
        "operating_point": {
            "min_alignment_aa": MIN_ALN_AA,
            "significant_evalue": SIGNIFICANT_EVALUE,
            "diamond_mode": a.diamond_mode,
            "min_orf_aa": a.min_orf_aa,
            "pfam_threshold": "curators' gathering cutoffs, not an E-value of ours",
            "divergent_arm_is_below_percent": a.diverged_below,
            "control_arm_is_at_or_above_percent": a.control_above,
            "mate_2_trimmed_nt": a.trim_r2,
            "mates_streamed_per_file": a.map_limit or "whole library",
        },
        "summary": {"divergent": summarise(div), "control_near_reference": summarise(con)},
        "how_to_read": (
            "Compare the two arms column by column. A field that reads the same on both "
            "is describing assembly, not divergence. subject_disagreement_flag fires on "
            "a mis-joined contig AND on a genuine virus whose genes have different "
            "nearest neighbours, so read it beside evenness_within_2x_of_median and "
            "never alone."),
        "blind_spots": [
            "An ORF with no subject and no domain is unexplained, not non-coding: past "
            "about 20-30% identity sequence search is at the end of its reach and the "
            "instrument that still works is structure, which is not attempted here",
            "pyrodigal in meta mode calls ORFs on ribosomal and other non-coding loci; "
            "an ORF is a prediction, not an observed protein",
            "the contigs are the study authors' assembly, so each is already a consensus "
            "and hides within-sample strain mixtures below about 3.8% divergence",
            "the control arm is near-reference contigs from this same assembly, which "
            "shares the assembler and the library with the test arm but not its "
            "divergence; it cannot control for anything the assembler does to divergent "
            "sequence specifically",
        ],
        "contigs_not_found_in_the_assembly": missing,
        "rows": rows,
        "seed": a.seed,
        "timings_seconds": {"diamond_viral": viral["seconds"], "diamond_sprot": sprot_sec,
                            "pfam": pf["seconds"], "total": round(time.time() - t0, 1)},
    }
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                             capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                    check=True, capture_output=True,
                                    text=True).stdout.strip())
    except Exception:
        sha, dirty = "unknown", True
    out["git_sha"] = sha + ("-dirty" if dirty else "")
    out["tree_dirty"] = dirty
    # Taken when the manifest is written, which pins the code that ran whatever the tree
    # looked like. Editing this script while a run of it is in flight stamps the new hash
    # onto the old code's output; that has cost this repo two restarts.
    out["script_sha256"] = sha256_of(Path(__file__).resolve())
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"wrote {a.out} in {out['timings_seconds']['total']}s", file=sys.stderr)


if __name__ == "__main__":
    main()
