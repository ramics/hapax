"""Which parts of a divergent protein are conserved, and which parts carry the distance.

WHY THIS EXISTS. The gene-level pass established that one identity per contig is a
property of the best gene rather than of the contig: a contig's own genes disagree by a
median of 12.8 points on divergent material against 2.45 on near-reference material
(`results/divergent_known/`). The holdout shows the same shape inside one genome, its
polymerase at 66.5% and its second gene at 52.2% on the same subject organism. The next
rung down is inside a single gene. A protein is not uniformly diverged: a catalytic core
is held by selection while the surface drifts, so an identity averaged over a whole gene
mixes two populations of columns that mean different things.

WHAT IS MEASURED. For every ORF with both a Pfam domain and a significant protein
alignment, every aligned column is assigned to a position in the query and labelled
inside or outside a domain envelope. The readout is the identity of the columns inside
against the identity of the columns outside, per arm, per domain family, and per contig.

WHY THIS IS NOT A NOVELTY SCORE. The quantity is a shape, not a ranking. A conserved core
beside a diverged surface is what an old, real protein looks like; it is also what a
well-sampled protein family looks like when the database holds a distant member. The
measurement says which columns carry the distance. It does not say whether the molecule
is new, and ordering contigs by it would order them by how well their family is sampled,
which §1.1 of SPEC.md forbids.

THREE CONTROL ARMS, because a per-column claim is exactly the kind that fails on its own
control. (1) Near-reference contigs from the same assembly through the identical
pipeline: if inside beats outside there too, the gap is a property of alignment and not
of divergence. (2) A dinucleotide-preserving shuffle of the same contigs: it says what
the whole chain reads on sequence carrying no information. (3) Column counts beside every
figure, because identity over a short stretch is quantised and a domain envelope can be
30 columns.

Usage:
  python3 -I experiments/conserved_positions.py --work W --out out.json \
      --arms-from results/divergent_known/divergent_known_s1.json
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "experiments"))

from divergent_contig_workup import (  # noqa: E402
    MIN_ALN_AA, SIGNIFICANT_EVALUE, call_orfs, read_fasta, sha256_of,
)

# A domain envelope can be short, and identity over a short stretch is a ratio of small
# whole numbers. Columns below this are counted but never turned into a percentage.
MIN_COLUMNS = 20

# Pfam names carry the function. The map is keyword to role, applied to the family name
# and its description, longest-match first, and anything unmatched stays "other" rather
# than being forced into a bucket.
ROLE_WORDS = [
    ("rna_dep_rnapol", "polymerase"), ("rdrp", "polymerase"),
    ("rna_pol", "polymerase"), ("viral_rdrp", "polymerase"),
    ("helicase", "helicase"), ("peptidase", "protease"), ("protease", "protease"),
    ("capsid", "capsid"), ("coat", "capsid"), ("maturation", "maturation"),
    ("lysis", "lysis"), ("methyltransferase", "methyltransferase"),
    ("glycosyl", "glycosyltransferase"), ("ribosomal", "ribosomal"),
]


def role_of(name: str, desc: str) -> str:
    hay = (name + " " + desc).lower()
    for word, role in ROLE_WORDS:
        if word in hay:
            return role
    return "other"


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Eulerian-walk shuffle preserving length and dinucleotide composition.

    Dinucleotides rather than single bases because single-base shuffling destroys local
    composition and makes the null easier than the real thing for reasons that have
    nothing to do with homology.
    """
    s = seq.upper()
    if len(s) < 4:
        return s
    for _ in range(20):
        nxt = collections.defaultdict(list)
        for a, b in zip(s, s[1:]):
            nxt[a].append(b)
        for v in nxt.values():
            rng.shuffle(v)
        out, cur = [s[0]], s[0]
        ok = True
        for _ in range(len(s) - 1):
            if not nxt[cur]:
                ok = False
                break
            cur = nxt[cur].pop()
            out.append(cur)
        if ok:
            return "".join(out)
    return s


def pfam_envelopes(orfs: list, hmm: Path, threads: int) -> dict:
    """Pfam-A with pyhmmer at the curators' gathering cutoffs, keeping the ENVELOPE.

    The envelope rather than the alignment coordinates: hmmer's envelope is its statement
    of where the domain is on the query, including the ends it is less certain of, and
    using the tighter alignment bounds would shrink the inside arm towards the best-fitting
    columns, which is the direction that would manufacture the effect this script is
    testing for.
    """
    import pyhmmer

    alpha = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=o["orf"].encode(),
                                       sequence=o["aa"]).digitize(alpha) for o in orfs]

    def txt(v) -> str:
        if isinstance(v, (bytes, bytearray)):
            return v.decode()
        return "" if v is None else str(v)

    env = collections.defaultdict(list)
    t0 = time.time()
    with pyhmmer.plan7.HMMFile(str(hmm)) as fh:
        for top in pyhmmer.hmmsearch(fh, seqs, cpus=threads, bit_cutoffs="gathering"):
            fam, desc = txt(top.query.name), txt(top.query.description)
            for h in top:
                if not h.included:
                    continue
                for d in h.domains:
                    if not d.included:
                        continue
                    env[txt(h.name)].append({
                        "pfam": fam, "description": desc,
                        "role": role_of(fam, desc),
                        "env_from": int(d.env_from), "env_to": int(d.env_to),
                        "bitscore": round(d.score, 1),
                    })
    return {"env": dict(env), "seconds": round(time.time() - t0, 1)}


def blastp_with_alignment(faa: Path, db: Path, out: Path, binpath: Path, threads: int,
                          mode: str) -> dict:
    """Best hit per ORF, carrying the aligned sequences so columns can be walked.

    qseq and sseq are in the output format because the question is per-column and no
    summary statistic can answer it. One target per query: the comparison is against the
    nearest relative the index holds, which is what the cascade used to remove the contig.
    """
    cmd = [str(binpath), "blastp", "-d", str(db), "-q", str(faa), "-o", str(out),
           "-p", str(threads), mode, "--quiet", "--evalue", "1e-3",
           "--max-target-seqs", "1", "--outfmt", "6",
           "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
           "qstart", "qend", "qlen", "slen", "qseq", "sseq", "stitle"]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    best = {}
    with out.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            row = {"subject": f[1], "pident": float(f[2]), "aln_aa": int(f[3]),
                   "evalue": float(f[4]), "bitscore": float(f[5]),
                   "qstart": int(f[6]), "qend": int(f[7]), "qlen": int(f[8]),
                   "subject_aa": int(f[9]), "qseq": f[10], "sseq": f[11],
                   "title": (f[12] if len(f) > 12 else "")[:200]}
            prev = best.get(f[0])
            if prev is None or row["bitscore"] > prev["bitscore"]:
                best[f[0]] = row
    return {"best": best, "seconds": round(time.time() - t0, 1)}


def columns(hit: dict) -> list:
    """Walk the alignment and return (query position, identical) for every aligned column.

    Insertions in the query advance the query position and are columns; deletions do not
    advance it and are dropped, because a column with no query position cannot be assigned
    to a domain envelope. Both are counted in the totals so the two arms cannot differ by
    how gaps were handled.
    """
    out = []
    q = hit["qstart"]
    for a, b in zip(hit["qseq"], hit["sseq"]):
        if a == "-":
            continue
        out.append((q, a == b and a != "-"))
        q += 1
    return out


def inside_outside(hit: dict, envs: list) -> dict:
    """Split one ORF's aligned columns by whether they fall in any Pfam envelope."""
    marks = set()
    for e in envs:
        marks.update(range(e["env_from"], e["env_to"] + 1))
    ins = [ok for q, ok in columns(hit) if q in marks]
    out = [ok for q, ok in columns(hit) if q not in marks]
    def pct(xs):
        return round(100.0 * sum(xs) / len(xs), 2) if len(xs) >= MIN_COLUMNS else None
    return {"columns_inside": len(ins), "columns_outside": len(out),
            "identity_inside_percent": pct(ins), "identity_outside_percent": pct(out)}


def med(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.median(xs), 2) if xs else None


def arm_rollup(rows: list) -> dict:
    """Per arm: the inside against outside comparison, paired where both are scorable.

    Paired rather than two independent medians: the two numbers come from one alignment,
    so the difference per ORF is the quantity with a control built into it, and a
    difference of medians would mix ORFs that have only one of the two.
    """
    both = [r for r in rows if r["identity_inside_percent"] is not None
            and r["identity_outside_percent"] is not None]
    diffs = [r["identity_inside_percent"] - r["identity_outside_percent"] for r in both]
    return {
        "orfs_with_a_significant_hit": len(rows),
        "orfs_with_a_pfam_domain": sum(1 for r in rows if r["domains"]),
        "orfs_scorable_both_sides": len(both),
        "median_identity_inside_percent": med([r["identity_inside_percent"] for r in both]),
        "median_identity_outside_percent": med([r["identity_outside_percent"] for r in both]),
        "median_paired_difference_points": med(diffs),
        "share_inside_higher": (round(sum(1 for d in diffs if d > 0) / len(diffs), 4)
                                if diffs else None),
        "median_columns_inside": med([r["columns_inside"] for r in both]),
        "median_columns_outside": med([r["columns_outside"] for r in both]),
    }


def by_role(rows: list) -> dict:
    out = {}
    g = collections.defaultdict(list)
    for r in rows:
        for role in {d["role"] for d in r["domains"]} or {"no domain"}:
            g[role].append(r)
    for role, rs in sorted(g.items()):
        both = [r for r in rs if r["identity_inside_percent"] is not None
                and r["identity_outside_percent"] is not None]
        if len(rs) < 5:
            continue
        out[role] = {
            "orfs": len(rs),
            "median_whole_orf_identity_percent": med([r["pident"] for r in rs]),
            "median_aln_aa": med([r["aln_aa"] for r in rs]),
            "scorable_both_sides": len(both),
            "median_identity_inside_percent":
                med([r["identity_inside_percent"] for r in both]),
            "median_identity_outside_percent":
                med([r["identity_outside_percent"] for r in both]),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=Path,
                    default=ROOT / "data/geo/sample_contigs_identified_norrna.fa")
    ap.add_argument("--arms-from", type=Path,
                    default=ROOT / "results/divergent_known/divergent_known_s1.json",
                    help="the divergent stage output, whose rows already carry the arm "
                         "each contig is in. Reusing it keeps one definition of the two "
                         "arms rather than recomputing a second one that could differ")
    ap.add_argument("--max-per-arm", type=int, default=300)
    ap.add_argument("--min-orf-aa", type=int, default=50)
    ap.add_argument("--pfam", type=Path, default=ROOT / "data/pfam/Pfam-A.hmm")
    ap.add_argument("--dmnd", type=Path, default=ROOT / "data/prot/viral.dmnd")
    ap.add_argument("--dmnd-sprot", type=Path, default=ROOT / "data/prot/sprot.dmnd")
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-mode", default="--sensitive")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    t_all = time.time()

    arms_doc = json.loads(a.arms_from.read_text())
    arms = {}
    for r in arms_doc.get("rows", []):
        arms.setdefault(r["arm"], []).append(r["contig"])
    allseqs = read_fasta(a.contigs)

    picked = {}
    for arm, names in arms.items():
        keep = [n for n in names if n in allseqs][: a.max_per_arm]
        for n in keep:
            picked[n] = arm
    print(f"{len(picked)} contigs: " +
          ", ".join(f"{k} {sum(1 for v in picked.values() if v == k)}"
                    for k in sorted(arms)), file=sys.stderr)

    seqs = {n: allseqs[n] for n in picked}
    shuf = {n + "__shuf": dinuc_shuffle(s, rng) for n, s in seqs.items()}

    orfs = call_orfs(seqs, a.min_orf_aa)
    sorfs = call_orfs(shuf, a.min_orf_aa)
    print(f"{len(orfs)} ORFs real, {len(sorfs)} shuffled", file=sys.stderr)

    secs = {}
    results = {}
    for tag, olist in (("real", orfs), ("shuffled", sorfs)):
        faa = a.work / f"orfs_{tag}.faa"
        faa.write_text("".join(f">{o['orf']}\n{o['aa']}\n" for o in olist))
        hits = {}
        for label, db in (("viral", a.dmnd), ("sprot", a.dmnd_sprot)):
            if not db.exists():
                continue
            r = blastp_with_alignment(faa, db, a.work / f"{tag}_{label}.tsv",
                                      a.diamond_bin, a.threads, a.diamond_mode)
            secs[f"blastp_{tag}_{label}"] = r["seconds"]
            for q, row in r["best"].items():
                row["index"] = label
                if q not in hits or row["bitscore"] > hits[q]["bitscore"]:
                    hits[q] = row
        pf = pfam_envelopes(olist, a.pfam, a.threads) if a.pfam.exists() else {"env": {}}
        secs[f"pfam_{tag}"] = pf.get("seconds")
        results[tag] = (olist, hits, pf["env"])

    def rows_for(tag):
        olist, hits, env = results[tag]
        out = []
        for o in olist:
            h = hits.get(o["orf"])
            if not h or h["evalue"] >= SIGNIFICANT_EVALUE or h["aln_aa"] < MIN_ALN_AA:
                continue
            envs = env.get(o["orf"], [])
            row = {
                "orf": o["orf"],
                "contig": o["contig"].replace("__shuf", ""),
                "arm": "shuffled" if tag == "shuffled"
                       else picked.get(o["contig"], "unknown"),
                "orf_aa": o["aa_len"],
                "pident": h["pident"], "aln_aa": h["aln_aa"], "evalue": h["evalue"],
                "subject": h["subject"], "subject_title": h["title"],
                "index": h["index"],
                "domains": [{k: d[k] for k in ("pfam", "role", "env_from", "env_to",
                                               "bitscore")} for d in envs],
            }
            row.update(inside_outside(h, envs))
            out.append(row)
        return out

    real_rows = rows_for("real")
    shuf_rows = rows_for("shuffled")
    by_arm = collections.defaultdict(list)
    for r in real_rows:
        by_arm[r["arm"]].append(r)
    by_arm["shuffled_null"] = shuf_rows

    # One contig at a time: which function on it is the diverged one. A contig whose
    # polymerase is close and whose capsid is far is a different object from the reverse,
    # and neither is visible in a single identity.
    per_contig = []
    cg = collections.defaultdict(list)
    for r in real_rows:
        cg[(r["arm"], r["contig"])].append(r)
    for (arm, contig), rs in cg.items():
        roles = {}
        for r in rs:
            for d in r["domains"]:
                roles.setdefault(d["role"], []).append(r["pident"])
        if len(roles) < 2:
            continue
        flat = {k: round(statistics.median(v), 1) for k, v in roles.items()}
        lo = min(flat, key=lambda k: flat[k])
        hi = max(flat, key=lambda k: flat[k])
        per_contig.append({
            "arm": arm, "contig": contig, "by_role": flat,
            "most_diverged_role": lo, "least_diverged_role": hi,
            "spread_points": round(flat[hi] - flat[lo], 1),
        })
    per_contig.sort(key=lambda r: -r["spread_points"])

    which = collections.Counter((r["arm"], r["most_diverged_role"]) for r in per_contig)

    out = {
        "unit": "one row per ORF with a significant protein alignment",
        "question": ("within one protein, do the columns inside a Pfam domain hold while "
                     "the columns outside it drift, and does that shape differ between "
                     "divergent and near-reference material"),
        "operating_point": {
            "significant_evalue": SIGNIFICANT_EVALUE,
            "min_alignment_aa": MIN_ALN_AA,
            "min_orf_aa": a.min_orf_aa,
            "min_columns_for_a_percentage": MIN_COLUMNS,
            "pfam_threshold": "curators' gathering cutoffs, envelope coordinates",
            "diamond_mode": a.diamond_mode,
            "arms_defined_by": str(a.arms_from),
            "max_contigs_per_arm": a.max_per_arm,
            "null": "dinucleotide-preserving shuffle of the same contigs, same chain",
        },
        "scope": {
            "contigs": len(picked),
            "contigs_per_arm": {k: sum(1 for v in picked.values() if v == k)
                                for k in sorted(set(picked.values()))},
            "orfs_real": len(orfs), "orfs_shuffled": len(sorfs),
            "orfs_scored_real": len(real_rows), "orfs_scored_shuffled": len(shuf_rows),
        },
        "inside_against_outside": {k: arm_rollup(v) for k, v in sorted(by_arm.items())},
        "by_domain_role": {k: by_role(v) for k, v in sorted(by_arm.items())},
        "which_function_is_diverging": {
            "readout": ("contigs carrying two or more domain roles, by which role is the "
                        "most diverged. A count, never an ordering of contigs"),
            "counts": {f"{arm} | {role}": n for (arm, role), n in which.most_common()},
            "contigs_with_two_or_more_roles": len(per_contig),
        },
        "blind_spots": [
            "a Pfam domain is a statement about a family that has been curated, so the "
            "inside arm is biased towards well-studied protein families and the outside "
            "arm towards everything else. The direction is fixed and cannot be removed "
            "by a control here",
            "domain envelopes and alignments are both chosen to maximise a score, so some "
            "agreement between them is expected before any biology; the near-reference "
            "arm is what says how much",
            "identity over a short envelope is quantised, so no percentage is reported "
            f"below {MIN_COLUMNS} columns and every figure carries its column count",
            "one best hit per ORF, so every distance is to the nearest relative the index "
            "holds and is an upper bound on identity",
            "this says which columns carry the distance. It does not say the molecule is "
            "new, and ranking by it would rank by how well a family is sampled",
        ],
        "rows": real_rows,
        "shuffled_rows": shuf_rows,
        "per_contig_roles": per_contig,
        "seed": a.seed,
        "timings_seconds": {**secs, "total": round(time.time() - t_all, 1)},
        "script_sha256": sha256_of(Path(__file__)),
    }
    try:
        out["git_sha"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                        capture_output=True, text=True,
                                        check=True).stdout.strip()
        out["tree_dirty"] = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                                capture_output=True, text=True,
                                                check=True).stdout.strip())
    except Exception:
        out["git_sha"] = "unknown"
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"wrote {a.out} in {time.time() - t_all:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
