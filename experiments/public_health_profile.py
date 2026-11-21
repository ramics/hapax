"""What is in this water, as a public-health reader would ask it, over the whole library.

THE GAP THIS FILLS. Everything else in this repo asks what the indexes CANNOT name,
because that is the research question. A public-health reader asks the opposite question
first - what IS here, how much of it, and is any of it something we screen for - and the
run directory could not answer it. The cascade records whether a read was removed and by
which tier, never what it was removed AS, so "47% of this library is bacterial" was not a
number this pipeline could produce. This script produces it, at the cascade's own
operating point, over every pair rather than a block.

WHAT IT IS AND IS NOT.

  It is a classified-read profile: for each read pair, the taxon Kraken2 calls it, rolled
  up to superkingdom and to species. Counted in read PAIRS, the occurrence-weighted unit
  this repo reports everything in.

  It is NOT an abundance of organisms. This is a metatranscriptome, so a read count is a
  transcript count filtered through extraction, rRNA content and sequencing depth; a
  bacterium with a thousand ribosomes contributes a thousand times what its genome count
  would suggest. A viral fraction of x% is "x% of the RNA we sequenced", never "x% of what
  is in the sewer". The output says so in every pile and the page prints it.

  It is NOT a diagnosis. A species call from a 107 nt read against a k-mer index is a
  claim about shared k-mers with a reference, and the confidence sweep below shows how
  much that claim moves with one parameter.

THE THRESHOLD IS THE MEASUREMENT. This repo's standing finding is that Kraken2's
confidence moves the classified fraction by an order of magnitude more than database size
does: the same PlusPF-16 index removes 42.8%, 21.7% or 3.5% of this library at confidence
0, 0.1 or 0.2. So every number here is reported at three confidences and the page may not
show one without the others. A single classified fraction is a choice wearing the costume
of a measurement.

WHAT IS SCREENED FOR, AND WHAT WAS NOT FOUND. `pathogen_panel.json` is a committed,
auditable list of taxa a public-health reader would want asked about. Every entry is
reported whether or not it was seen, because "we looked for measles and did not see it" is
a result and an empty row is not. A not-detected row carries the rule-of-three upper bound
over the pairs actually classified, which is a limit of THIS SEQUENCING and not of the
water: untargeted shotgun is reported in the literature as insensitive for human pathogens
at wastewater concentrations, and this library is untargeted. The bound is printed with
that sentence attached.

CALIBRATION AND DIRECTION OF BIAS.

  - Unclassified is reported as its own pile and is the largest one at any useful
    confidence. It is not evidence of novelty: results/full_residue_search/ measures the
    same reads against protein indexes and finds 92% of them point at nothing nameable,
    which is a statement about the references.
  - A superkingdom rollup inherits every misassignment below it, and the direction is
    known: Kraken2 over-calls at low confidence, so the viral and bacterial fractions at
    confidence 0 are upper bounds and at 0.2 lower bounds. The sweep brackets it rather
    than estimating it.
  - Mate 2 carries a non-template prefix and is trimmed by the same amount the cascade
    uses, so the two agree. 16 is the better Kraken2 operating point and 14 is the
    cascade's; this script takes it as an argument and records it.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PANEL = Path(__file__).resolve().parent / "pathogen_panel.json"

#: The confidences this repo sweeps everywhere. Fixed in the script, not a flag, because
#: a run that reported one of them would be the failure this file exists to prevent.
CONFIDENCES = (0.0, 0.1, 0.2)



def say(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def load_taxonomy(db: Path) -> tuple[dict, dict, dict]:
    """parent, rank and scientific name per taxid, from the database's own dump files.

    From the database directory rather than a separate NCBI download, so the tree that
    interprets a call is the tree the call was made against. A taxonomy fetched today
    would quietly renumber clades that this database still places elsewhere.
    """
    parent, rank, name = {}, {}, {}
    with open(db / "nodes.dmp") as h:
        for line in h:
            f = [x.strip() for x in line.split("|")]
            parent[int(f[0])] = int(f[1])
            rank[int(f[0])] = f[2]
    with open(db / "names.dmp") as h:
        for line in h:
            f = [x.strip() for x in line.split("|")]
            if f[3] == "scientific name":
                name[int(f[0])] = f[1]
    return parent, rank, name


#: The groups a reader means by "viral", "bacterial" and so on, as TAXIDS and not as a
#: rank name. NCBI renamed the `superkingdom` rank to `domain` and files Viruses under
#: `acellular root`, so a rollup that keys on the rank string silently puts every read in
#: "above superkingdom" - which is exactly what the first run of this script did. The
#: taxids are stable; the rank names were not.
GROUP_TAXIDS = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses",
                12884: "Viroids"}


def lineage_maps(parent: dict, rank: dict, name: dict) -> tuple[dict, dict]:
    """For every taxid: which of GROUP_TAXIDS it sits under, and its species ancestor.

    Walked once per taxid with the path memoised on the way back down, because the
    alternative is walking the tree per read over 8.7 M of them. A cycle or a missing
    parent terminates the walk rather than hanging: a truncated dump should give a
    missing group, not an infinite loop.
    """
    grp: dict[int, str | None] = {}
    spp: dict[int, int | None] = {}
    for leaf in parent:
        if leaf in grp:
            continue
        chain, cur, g, sp, guard = [], leaf, None, None, 0
        while cur and cur != 1 and guard < 256:
            if cur in grp:
                g, sp = grp[cur], spp[cur]
                break
            chain.append(cur)
            if cur in GROUP_TAXIDS:
                g = GROUP_TAXIDS[cur]
                break
            nxt = parent.get(cur)
            if nxt == cur:
                break
            cur = nxt
            guard += 1
        for c in reversed(chain):
            if rank.get(c) == "species":
                sp = c
            grp[c], spp[c] = g, sp
    return grp, spp


def write_fastq(src: Path, dest: Path, trim: int, limit: int) -> int:
    """Copy (and for mate 2, trim) a FASTQ, so the classifier sees the cascade's reads.

    Kraken2 is given plain files rather than a stream because it reads the two mates in
    lockstep; the trim is applied here so that the only difference between this profile
    and the cascade's own Kraken2 tier is the taxonomy rollup.
    """
    n = 0
    with gzip.open(src, "rt") as h, open(dest, "w") as out:
        while True:
            a = h.readline()
            if not a:
                break
            s = h.readline().rstrip("\n")
            p = h.readline()
            q = h.readline().rstrip("\n")
            if trim:
                s, q = s[trim:], q[trim:]
            out.write(a)
            out.write(s + "\n")
            out.write(p)
            out.write(q + "\n")
            n += 1
            if limit and n >= limit:
                break
    return n


def shuffle_fastq(src: Path, dest: Path, seed: int) -> None:
    """A composition-matched shuffle of every read, for the null pass.

    WHAT IT CONTROLS AND WHAT IT DOES NOT, which matters more than the number it
    produces. Shuffling a read keeps its base composition and destroys every k-mer, so
    the null answers one question only: does this classifier manufacture panel hits out
    of base composition. It does NOT control the failure that actually threatens a
    low-count detection here, which is misassignment between close relatives - a gut
    Bacteroides read landing on a Shigella k-mer is a real k-mer match and a shuffle says
    nothing about it. So a near-zero null is necessary and nowhere near sufficient, and
    the output says that rather than letting a clean null stand in for a clean panel.
    """
    rng = random.Random(seed)
    with open(src) as h, open(dest, "w") as out:
        while True:
            a = h.readline()
            if not a:
                break
            s = list(h.readline().rstrip("\n"))
            p = h.readline()
            q = h.readline()
            rng.shuffle(s)
            out.write(a)
            out.write("".join(s) + "\n")
            out.write(p)
            out.write(q)


def classify(classify_bin: Path, db: Path, fq1: Path, fq2: Path, out: Path,
             threads: int, confidence: float) -> tuple[collections.Counter, float]:
    """One Kraken2 pass; returns pairs per taxid.

    The per-read output is streamed into a counter and never kept: 8.7 M lines per
    confidence is 2 GB on disk for a count that fits in a dictionary.
    """
    t0 = time.perf_counter()
    # No -M. Memory-mapping the 8 GB hash makes every minimizer lookup a random read
    # from disk: measured at 100,000 pairs, the mapped run had not finished in 9 minutes
    # while the loaded run took 27 s, of which 0.26 s was classification. Loading costs
    # 8 GB of RAM for the length of one pass, which this container has.
    cmd = [str(classify_bin), "-H", str(db / "hash.k2d"), "-t", str(db / "taxo.k2d"),
           "-o", str(db / "opts.k2d"), "-p", str(threads), "-P",
           "-T", str(confidence), "-O", str(out), str(fq1), str(fq2)]
    subprocess.run(cmd, check=True, capture_output=True)
    counts: collections.Counter = collections.Counter()
    with open(out) as h:
        for line in h:
            f = line.split("\t", 3)
            counts[int(f[2])] += 1
    out.unlink(missing_ok=True)
    return counts, time.perf_counter() - t0


def rollup(counts: collections.Counter, sk_of: dict, sp_of: dict, name: dict,
           rank_of: dict, top: int) -> dict:
    """Superkingdom shares, the top species in each, and the unclassified pile."""
    total = sum(counts.values())
    by_sk: collections.Counter = collections.Counter()
    sp_by_sk: dict[str, collections.Counter] = {}
    for taxid, n in counts.items():
        if taxid == 0:
            by_sk["unclassified"] += n
            continue
        sk = sk_of.get(taxid) or "classified, but above any of these groups"
        by_sk[sk] += n
        sp = sp_of.get(taxid)
        # A call that never reached species is labelled with the rank it DID reach. The
        # first version of this listed "Bacteria" among the top bacterial taxa, which
        # reads as an organism and is really 321,748 pairs the classifier could not place
        # below the domain. That is a fact about the reads and it should look like one.
        if sp:
            key = name.get(sp) or f"taxid {sp}"
        else:
            key = (name.get(taxid) or f"taxid {taxid}")
            r = rank_of.get(taxid)
            key += f" (no deeper than {r})" if r and r != "no rank" else " (no rank given)"
        sp_by_sk.setdefault(sk, collections.Counter())[key] += n
    return {
        "pairs": total,
        "by_superkingdom": [
            {"group": g, "pairs": n, "fraction": round(n / total, 6)}
            for g, n in by_sk.most_common()],
        "top_taxa": {
            g: [{"name": k, "pairs": n, "fraction": round(n / total, 8)}
                for k, n in c.most_common(top)]
            for g, c in sp_by_sk.items()},
    }


def panel_result(counts: collections.Counter, panel: dict, parent: dict, name: dict,
                 total: int) -> list[dict]:
    """Every watchlist entry, found or not, with the pairs under its whole clade.

    Counted over the CLADE and not the exact taxid: a read called `Enterovirus C` and a
    read called `Human poliovirus 1` are both enterovirus C to a reader, and requiring an
    exact species match would report absence for a taxon the classifier placed one rank
    up. The clade is walked from the leaf upwards per observed taxid, which is cheap
    because only observed taxids are walked.
    """
    want = {int(e["taxid"]): e for e in panel["entries"]}
    hits: dict[int, int] = {t: 0 for t in want}
    for taxid, n in counts.items():
        cur, seen = taxid, 0
        while cur and cur not in (0, 1) and seen < 64:
            if cur in hits:
                hits[cur] += n
                break
            cur = parent.get(cur, 0)
            seen += 1
    rows = []
    for e in panel["entries"]:
        k = hits[int(e["taxid"])]
        rows.append({
            "name": e["name"], "taxid": e["taxid"], "group": e["group"],
            "why_on_the_panel": e["why"],
            "pairs": k,
            "fraction": (k / total) if total else None,
            "detected": k > 0,
            # Rule of three: with zero observations in n trials the 95% upper bound on the
            # rate is 3/n. It is a bound on what this SEQUENCING could have seen, not on
            # what is in the water.
            "upper_bound_if_absent": (3.0 / total) if (k == 0 and total) else None,
        })
    rows.sort(key=lambda r: (-r["pairs"], r["name"]))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r1", type=Path, required=True)
    ap.add_argument("--r2", type=Path, required=True)
    ap.add_argument("--trim-r2", type=int, default=14)
    ap.add_argument("--db", type=Path, default=ROOT / "data/k2_pluspf08")
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/"
                                 "kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--panel", type=Path, default=PANEL)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--limit-pairs", type=int, default=0,
                    help="0 means the whole library; a small value is a pilot")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()

    fq1, fq2 = args.work / "r1.fq", args.work / "r2.fq"
    t0 = time.perf_counter()
    n1 = write_fastq(args.r1, fq1, 0, args.limit_pairs)
    n2 = write_fastq(args.r2, fq2, args.trim_r2, args.limit_pairs)
    extract_s = time.perf_counter() - t0
    if n1 != n2:
        raise SystemExit(f"mate counts differ: {n1} and {n2}")
    say(f"{n1} pairs staged in {extract_s:.0f}s")

    parent, rank, name = load_taxonomy(args.db)
    sk_of, sp_of = lineage_maps(parent, rank, name)
    panel = json.loads(args.panel.read_text())
    say(f"taxonomy: {len(parent)} nodes")

    at: dict[str, dict] = {}
    for conf in CONFIDENCES:
        counts, secs = classify(args.classify_bin, args.db, fq1, fq2,
                                args.work / f"calls_{conf}.k2", args.threads, conf)
        total = sum(counts.values())
        r = rollup(counts, sk_of, sp_of, name, rank, args.top)
        r["panel"] = panel_result(counts, panel, parent, name, total)
        r["classify_seconds"] = round(secs, 1)
        at[f"{conf}"] = r
        cl = total - counts.get(0, 0)
        say(f"confidence {conf}: {cl} of {total} classified ({100 * cl / total:.2f}%), "
            f"{secs:.0f}s")
    # The null, through the same classifier at the same thresholds, in the same run.
    n1f, n2f = args.work / "null1.fq", args.work / "null2.fq"
    shuffle_fastq(fq1, n1f, args.seed)
    shuffle_fastq(fq2, n2f, args.seed + 1)
    null = {}
    for conf in CONFIDENCES:
        counts, secs = classify(args.classify_bin, args.db, n1f, n2f,
                                args.work / f"null_{conf}.k2", args.threads, conf)
        total = sum(counts.values())
        cl = total - counts.get(0, 0)
        null[f"{conf}"] = {
            "pairs": total, "classified": cl,
            "fraction_classified": round(cl / total, 8) if total else None,
            "panel_rows_detected": sum(
                1 for r in panel_result(counts, panel, parent, name, total)
                if r["detected"]),
            "panel_pairs": sum(r["pairs"] for r in
                               panel_result(counts, panel, parent, name, total)),
        }
        say(f"null at {conf}: {cl} classified, "
            f"{null[f'{conf}']['panel_rows_detected']} panel rows")
    for f in (fq1, fq2, n1f, n2f):
        f.unlink(missing_ok=True)

    out = {
        "unit": "read PAIRS, the occurrence-weighted unit; never distinct k-mers",
        "library_pairs": n1,
        "whole_library": not args.limit_pairs,
        "operating_point": {
            "classifier": "kraken2 2.17.2",
            "database": args.db.name,
            "database_bytes": sum(f.stat().st_size for f in args.db.glob("*.k2d")),
            "confidences": list(CONFIDENCES),
            "trim_r2": args.trim_r2,
            "paired": True,
        },
        "at_confidence": at,
        "shuffled_null": null,
        "null_reading": ("the same library with every read shuffled, through the same "
                         "classifier at the same thresholds. It controls one thing: "
                         "whether base composition alone manufactures a call. It does "
                         "NOT control misassignment between close relatives, which is "
                         "the failure that threatens a low-count detection here, and no "
                         "control in this run does"),
        "reading": ("a classified-read profile of a METATRANSCRIPTOME. A fraction here is "
                    "a fraction of the RNA that was sequenced, not of the organisms in "
                    "the water: rRNA content, extraction and depth all sit between the "
                    "two, and a bacterium with many ribosomes contributes far more reads "
                    "than its cell count would suggest"),
        "bias": ("Kraken2 over-calls at low confidence, so the fractions at confidence 0 "
                 "are upper bounds and those at 0.2 lower bounds. The sweep brackets the "
                 "answer; it does not estimate it. Unclassified is not evidence of "
                 "novelty - see results/full_residue_search/, where the same reads are "
                 "searched against protein indexes"),
        "panel_note": ("every entry is reported whether or not it was seen. A "
                       "not-detected row carries the rule-of-three 95% upper bound over "
                       "the pairs classified, which bounds what THIS SEQUENCING could "
                       "have seen. Untargeted shotgun is reported as insensitive for "
                       "human pathogens at wastewater concentrations, so absence here is "
                       "much weaker evidence than absence from a targeted assay"),
        "panel_source": str(args.panel.relative_to(ROOT)),
        "seed": args.seed,
        "wall_seconds": round(time.perf_counter() - t_start, 1),
        **provenance(__file__),
    }
    dest = args.out / f"public_health_profile_s{args.seed}.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    say(f"{out['wall_seconds']}s -> {dest}")


if __name__ == "__main__":
    main()
