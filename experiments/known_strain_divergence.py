"""How far this sample's own known viruses have drifted from the reference.

THE QUESTION, which is SPEC.md's and a public-health reader's at once: when we find a
virus we already know, is it the one in the database or a divergent strain of it? That is
the question wastewater surveillance is actually used for - a new variant shows up in the
sewer before it shows up in a clinic - and nothing in this repo answered it, because
everything here is built to ask what the indexes CANNOT name.

WHY CONTIGS AND NOT READS, which is the whole design. A 107 nt read is about 35 amino
acids once translated, and strain-level variation within a known virus is mostly
synonymous. At 35 aa a protein identity cannot separate a reference strain from a
divergent one: this repo has measured the chance band for exactly that query length and
found 36-48 nt matching at ~89% on reads that are pure background. So this script works
on the sample's OWN SPAdes contigs, where an alignment runs to hundreds of amino acids
and a percentage point of identity is a measurement rather than a coin flip. The read
level is reported as the blind spot it is, not quietly skipped.

WHAT IT DOES. Kraken2 over this sample's contigs at the cascade's operating point, then
DIAMOND blastx of every contig the classifier places under a watchlist taxon or anywhere
in Viruses, against the same RefSeq viral protein database the cascade uses. Per taxon it
reports the contigs, their length, and the identity distribution WITH the aligned length
beside it, because this repo has had to withdraw an identity quoted without one.

HOW TO READ A ROW, and the direction of every bias in it.

  - High identity over a long alignment: this is the reference organism, or close to it.
  - LOW identity over a LONG alignment: a divergent strain of something known. This is
    the row a public-health reader wants, and it is the only row this instrument can
    honestly produce.
  - Low identity over a SHORT alignment: says nothing. Reported, and marked.

  Protein identity understates nucleotide divergence, because synonymous changes are
  invisible to it. So a row saying 97% protein identity is NOT a claim of 97% nucleotide
  identity - the nucleotide figure is lower and this instrument cannot see how much
  lower. Direction known, magnitude not: a nucleotide comparison needs reference genomes,
  which this sandbox cannot fetch. That is the stated blind spot rather than an estimate.

  The contigs are the study authors' assembly of this sample, so a contig is already a
  consensus over reads and hides within-sample strain mixtures. Assembly fuses lineages
  below about 3.8% divergence at any depth (results/pair_ladder/), which is exactly the
  range a "new strain of a known virus" lives in. So this measurement finds divergence
  from the DATABASE and is blind to divergence WITHIN the sample. Both halves are
  printed.

  A taxon absent from this table is absent from the contigs, which is a higher bar than
  absent from the reads: a virus present at low depth never assembles. The read-level
  panel in public_health_profile.py is the sensitive instrument and this is the specific
  one; they are reported side by side and must not be merged.
"""

from __future__ import annotations

import argparse
import collections
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from public_health_profile import load_taxonomy, lineage_maps, say  # noqa: E402
from provenance import rel  # noqa: E402
from subtract_known import provenance  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

#: An alignment at least this long in amino acids before an identity is read as a
#: measurement rather than noise. 60 aa is ~180 nt, well past the chance band this repo
#: measured for short queries, and short enough that a single viral gene fragment
#: qualifies.
MIN_ALN_AA = 60


def read_fasta(path: Path):
    name, buf = None, []
    with open(path) as h:
        for line in h:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def classify_contigs(classify_bin: Path, db: Path, fa: Path, out: Path, threads: int,
                     confidence: float) -> dict:
    cmd = [str(classify_bin), "-H", str(db / "hash.k2d"), "-t", str(db / "taxo.k2d"),
           "-o", str(db / "opts.k2d"), "-p", str(threads),
           "-T", str(confidence), "-O", str(out), str(fa)]
    subprocess.run(cmd, check=True, capture_output=True)
    calls = {}
    with open(out) as h:
        for line in h:
            f = line.split("\t", 3)
            calls[f[1]] = int(f[2])
    out.unlink(missing_ok=True)
    return calls


def diamond_best(dbin: Path, dmnd: Path, fa: Path, tsv: Path, threads: int,
                 mode: str) -> dict:
    subprocess.run([str(dbin), "blastx", "--db", str(dmnd), "--query", str(fa),
                    "--out", str(tsv), "--threads", str(threads), "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                    "--max-target-seqs", "1", "--evalue", "1e-5", mode, "--quiet"],
                   check=True)
    best: dict[str, dict] = {}
    with open(tsv) as h:
        for line in h:
            q, s, pid, ln, ev, bits = line.rstrip("\n").split("\t")
            bits = float(bits)
            if bits > best.get(q, {}).get("bitscore", -1.0):
                best[q] = {"subject": s, "pident": float(pid), "aln_aa": int(ln),
                           "evalue": float(ev), "bitscore": bits}
    tsv.unlink(missing_ok=True)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=Path, default=ROOT / "data/geo/sample_contigs.fa")
    ap.add_argument("--db", type=Path, default=ROOT / "data/k2_pluspf08")
    ap.add_argument("--dmnd", type=Path, default=ROOT / "data/prot/viral.dmnd")
    ap.add_argument("--panel", type=Path,
                    default=Path(__file__).resolve().parent / "pathogen_panel.json")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/"
                                 "kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-mode", default="--sensitive")
    # 0.0, and the reason is a measurement rather than a preference. Kraken2's
    # confidence is the fraction of a sequence's minimizers that fall in the called
    # clade, so it is length-dependent in a way the cascade's read-level operating point
    # does not carry over: at 0.1 the viral calls on these contigs collapse from 1,021 to
    # 11, and the 11 survivors are the SHORT ones - a 9 kb viral contig has thousands of
    # minimizers and a divergent virus matches few of them. The filter that replaces it
    # here is stronger and is the point of this script: a row only exists if DIAMOND
    # finds a protein alignment at E<1e-5 over at least 60 amino acids. Kraken2 proposes,
    # the alignment disposes.
    ap.add_argument("--confidence", type=float, default=0.0)
    ap.add_argument("--also-report-confidence", type=float, default=0.1,
                    help="a second pass, reported but not used, to keep the length "
                         "effect on the record beside the numbers it would change")
    ap.add_argument("--min-contig", type=int, default=500)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()

    parent, rank, name = load_taxonomy(args.db)
    sk_of, sp_of = lineage_maps(parent, rank, name)
    panel = json.loads(args.panel.read_text())
    panel_taxa = {int(e["taxid"]): e for e in panel["entries"]}

    seqs = {n: s for n, s in read_fasta(args.contigs) if len(s) >= args.min_contig}
    say(f"{len(seqs)} contigs at or above {args.min_contig} nt")

    calls = classify_contigs(args.classify_bin, args.db, args.contigs,
                             args.work / "contigs.k2", args.threads, args.confidence)
    say(f"{sum(1 for v in calls.values() if v)} of {len(calls)} contigs classified")
    other = classify_contigs(args.classify_bin, args.db, args.contigs,
                             args.work / "contigs2.k2", args.threads,
                             args.also_report_confidence)
    length_effect = {
        f"{args.confidence}": sum(1 for v in calls.values()
                                  if v and sk_of.get(v) == "Viruses"),
        f"{args.also_report_confidence}": sum(1 for v in other.values()
                                              if v and sk_of.get(v) == "Viruses"),
    }
    say(f"viral contigs by confidence: {length_effect}")

    def panel_ancestor(taxid: int):
        cur, seen = taxid, 0
        while cur and cur not in (0, 1) and seen < 64:
            if cur in panel_taxa:
                return cur
            cur = parent.get(cur, 0)
            seen += 1
        return None

    keep = {}
    for cname, seq in seqs.items():
        t = calls.get(cname, 0)
        if not t:
            continue
        pa = panel_ancestor(t)
        if pa is None and sk_of.get(t) != "Viruses":
            continue
        keep[cname] = (t, pa, seq)
    say(f"{len(keep)} contigs are viral or under a watchlist taxon")

    qfa = args.work / "query.fa"
    with open(qfa, "w") as h:
        for cname, (_t, _pa, seq) in keep.items():
            h.write(f">{cname}\n{seq}\n")
    best = diamond_best(args.diamond_bin, args.dmnd, qfa,
                        args.work / "hits.tsv", args.threads, args.diamond_mode)
    qfa.unlink(missing_ok=True)
    say(f"{len(best)} contigs have a significant protein hit")

    groups: dict[str, dict] = {}
    for cname, (t, pa, seq) in keep.items():
        key = name.get(pa) if pa else (name.get(sp_of.get(t) or t) or f"taxid {t}")
        g = groups.setdefault(key, {
            "name": key,
            "on_the_panel": pa is not None,
            "panel_group": panel_taxa[pa]["group"] if pa else None,
            "contigs": 0, "total_nt": 0, "longest_nt": 0,
            "with_protein_hit": 0, "rows": []})
        g["contigs"] += 1
        g["total_nt"] += len(seq)
        g["longest_nt"] = max(g["longest_nt"], len(seq))
        hit = best.get(cname)
        if hit:
            g["with_protein_hit"] += 1
            g["rows"].append({"contig": cname, "nt": len(seq), **hit})

    out_groups = []
    for g in groups.values():
        scorable = [r for r in g["rows"] if r["aln_aa"] >= MIN_ALN_AA]
        g["scorable_alignments"] = len(scorable)
        if scorable:
            pid = sorted(r["pident"] for r in scorable)
            aln = sorted(r["aln_aa"] for r in scorable)
            g["median_identity_percent"] = round(statistics.median(pid), 1)
            g["identity_range"] = [round(pid[0], 1), round(pid[-1], 1)]
            g["median_aligned_aa"] = int(statistics.median(aln))
            g["longest_aligned_aa"] = aln[-1]
            far = min(scorable, key=lambda r: r["pident"])
            g["most_divergent"] = far
            g["reading"] = (
                "a divergent strain of something known: the alignment is long enough for "
                "the identity to be a measurement"
                if g["median_identity_percent"] < 90 else
                "close to the reference over a long alignment")
        else:
            g["median_identity_percent"] = None
            g["reading"] = ("no alignment reached "
                            f"{MIN_ALN_AA} amino acids, so no identity here is readable")
        g["rows"] = sorted(g["rows"], key=lambda r: -r["nt"])[:5]
        out_groups.append(g)
    out_groups.sort(key=lambda g: (-g["total_nt"]))

    out = {
        "unit": "one row per taxon, over this sample's own assembled contigs",
        "question": ("when we find a virus we already know, is it the one in the database "
                     "or a divergent strain of it"),
        "operating_point": {
            "contigs": rel(args.contigs),
            "min_contig_nt": args.min_contig,
            "kraken2_confidence": args.confidence,
            "kraken2_database": args.db.name,
            "protein_database": rel(args.dmnd),
            "diamond_mode": args.diamond_mode,
            "diamond_evalue": 1e-5,
            "min_alignment_aa_for_an_identity": MIN_ALN_AA,
        },
        "confidence_is_length_dependent": {
            "viral_contigs_by_confidence": length_effect,
            "reading": ("Kraken2 confidence is the fraction of a sequence's minimizers "
                        "in the called clade, so the threshold a 107 nt read is scored "
                        "at does not transfer to a 9 kb contig: the survivors at the "
                        "higher threshold are the SHORT contigs, which is the opposite "
                        "of what a reader would assume. This script therefore takes the "
                        "assignment at the lower threshold and requires a protein "
                        "alignment before a row exists"),
        },
        "contigs_considered": len(seqs),
        "contigs_viral_or_watchlist": len(keep),
        "contigs_with_protein_hit": len(best),
        "taxa": out_groups,
        "how_to_read": ("low identity over a LONG alignment is a divergent strain and is "
                        "the row worth attention. Low identity over a short alignment "
                        "says nothing and is marked as unreadable rather than scored"),
        "blind_spots": [
            "protein identity understates nucleotide divergence, because synonymous "
            "changes are invisible to it. Direction known, magnitude not: a nucleotide "
            "comparison needs reference genomes this sandbox cannot fetch",
            "a contig is a consensus over reads, and assembly fuses lineages below about "
            "3.8% divergence at any depth (results/pair_ladder/), so this finds "
            "divergence from the DATABASE and is blind to strain mixtures WITHIN the "
            "sample",
            "a virus present at low depth never assembles, so absence from this table is "
            "a weaker statement than absence from the read-level panel",
            "at read level this measurement is impossible: 107 nt is about 35 amino "
            "acids, where strain-level variation is mostly synonymous and the identity "
            "sits in the chance band this repo measured",
        ],
        "seed": args.seed,
        "wall_seconds": round(time.perf_counter() - t_start, 1),
        **provenance(__file__),
    }
    dest = args.out / f"known_strain_divergence_s{args.seed}.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    say(f"{len(out_groups)} taxa -> {dest}")


if __name__ == "__main__":
    main()
