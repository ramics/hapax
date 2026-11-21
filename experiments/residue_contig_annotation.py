"""What the residue-only contigs ARE, read three ways, against controls that cannot be
explained by length.

THE QUESTION. Tier 3 of the cascade is reads that no index could place. Assembling only
those reads gives a small, diffuse set of contigs - 3,695 contigs, 1.31 Mb, median 331 nt
from 1,372,021 unplaced pairs, with 3.8% of those pairs mapping back. That shape already
says the unplaced fraction is not a few hidden abundant genomes. What it does not say is
what the sequence is. This script asks three questions of it that do not need an index to
name the organism:

  coding     does it look like it codes? Prodigal in metagenomic mode, scored as coding
             DENSITY (fraction of contig bases inside a predicted CDS) rather than ORF
             count, because count is length in disguise.
  Rfam       is it structured non-coding RNA? cmscan with the Rfam gathering thresholds.
             This is the one axis the whole project has never run on the residue, and it
             is the obvious blind spot of a protein-first cascade.
  Pfam       does any protein family recognise the predicted peptides, and in particular
             any RdRp family? hmmsearch at the family gathering threshold, which is the
             curators' own operating point rather than one chosen here.

WHY LENGTH-MATCHED CONTROLS, AND WHY THEY ARE THE WHOLE EXPERIMENT. All three fields rise
with contig length: a longer contig has more bases to be coding, more chance of carrying
an Rfam model, more peptide for Pfam to recognise. The residue assembly has a median of
331 nt against the sample's own contigs' very different distribution, so an unmatched
comparison would measure the assembly's length distribution and report it as biology. The
control here is therefore drawn from the sample's own contigs BINNED BY LENGTH and
sampled to match the residue's histogram bin for bin, so the two piles have the same
length distribution by construction. CLAUDE.md's rule that the control is read first
applies with force: if the control moves with the residue, the field is reading length.

AND WHY A SHUFFLE NULL BESIDE IT. The length-matched control says what the field reads on
real sequence we could identify. A dinucleotide-preserving shuffle of the residue contigs
themselves says what it reads on sequence carrying no information at all but the same
composition. Mononucleotide shuffling is not enough for coding density: stop-codon
frequency is a function of dinucleotide composition, so a mononucleotide shuffle would
hand Prodigal longer open frames than the real sequence and inflate the apparent signal.
Both controls ship with every number, which is the project rule.

WHAT THIS CANNOT SAY. A Pfam or Rfam hit on a residue contig means the cascade's read-level
indexes missed it, not that the contig is novel - that is tier 2 by definition, and the
point of running this is to size tier 2 from the residue side. Absence of a hit is a
statement about Pfam and Rfam at their gathering thresholds, which are recorded, and about
331 nt of sequence, which is short for both tools. Every rate here ships with its Wilson
interval and its operating point.
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
import time
from collections import Counter
from math import sqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

# Pfam families that are RNA-dependent RNA polymerases or their close relatives. Named by
# accession rather than by name match, because a substring search for "RdRp" also catches
# unrelated descriptions and misses the ones curated under other names.
RDRP_ACCESSIONS = {
    "PF00680",  # RdRP_1, picornavirus-like
    "PF00978",  # RdRP_2, alphavirus-like
    "PF00998",  # RdRP_3, carmovirus-like
    "PF02123",  # RdRP_4, birnavirus-like
    "PF00946",  # Mononeg_RNA_pol
    "PF00602",  # Flavi_NS5
    "PF05183",  # RdRP, mitovirus / narnavirus
    "PF04196",  # Birna_RdRp
    "PF00078",  # RVT_1, reverse transcriptase - a sister signal, reported separately
    "PF03431",  # RNA_replicase_B
    "PF08467",  # Peptidase_C37 (picornaviral polyprotein neighbour)
}
RT_ACCESSIONS = {"PF00078"}


def wilson(k: int, n: int, z: float = 1.96) -> dict | None:
    """A rate with its interval. A bare fraction of a few thousand contigs is unreadable."""
    if not n:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {"point": p, "lo": max(0.0, c - h), "hi": min(1.0, c + h), "k": k, "n": n}


def mean_ci(xs: list[float], z: float = 1.96) -> dict | None:
    """Mean with a normal interval. Reported beside the median because the distributions
    here are skewed and the two disagreeing is itself informative."""
    if not xs:
        return None
    n = len(xs)
    m = sum(xs) / n
    if n < 2:
        return {"mean": m, "lo": m, "hi": m, "median": m, "n": n}
    var = sum((x - m) ** 2 for x in xs) / (n - 1)
    h = z * sqrt(var / n)
    s = sorted(xs)
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    return {"mean": m, "lo": m - h, "hi": m + h, "median": med, "n": n}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")


def length_bin(n: int) -> int:
    """Log-spaced bins. Linear bins would put every short contig in one bucket and leave
    the long tail unmatched, which is the half of the distribution the fields care about."""
    b = 0
    x = 200
    while n >= x:
        x = int(x * 1.25)
        b += 1
    return b


def length_matched_sample(residue: list[tuple[str, str]],
                          pool: list[tuple[str, str]],
                          rng: random.Random) -> tuple[list[tuple[str, str]], dict]:
    """Draw from `pool` so the drawn pile's length histogram matches `residue`'s, bin for
    bin, without replacement. Bins that the pool cannot fill are recorded rather than
    silently short-changed: an unreported shortfall would reintroduce exactly the length
    confound the matching exists to remove.
    """
    by_bin: dict[int, list[tuple[str, str]]] = {}
    for rec in pool:
        by_bin.setdefault(length_bin(len(rec[1])), []).append(rec)
    for v in by_bin.values():
        rng.shuffle(v)
    want = Counter(length_bin(len(s)) for _, s in residue)
    out, short = [], {}
    for b, k in sorted(want.items()):
        have = by_bin.get(b, [])
        take = min(k, len(have))
        out.extend(have[:take])
        if take < k:
            short[str(b)] = {"wanted": k, "available": take}
    return out, short


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erickson dinucleotide-preserving shuffle, by Eulerian walk.

    Mononucleotide shuffling is wrong for this experiment: stop-codon density depends on
    dinucleotide composition, so a mononucleotide null hands Prodigal longer open reading
    frames than the real sequence and makes the coding signal look larger than it is.
    """
    s = seq.upper()
    if len(s) < 4:
        return s
    edges: dict[str, list[str]] = {}
    for a, b in zip(s, s[1:]):
        edges.setdefault(a, []).append(b)
    last = s[-1]
    # Pick, for each vertex but the last, a final edge forming a tree rooted at `last`;
    # retry until the choice is connected, which is what makes the walk Eulerian.
    for _ in range(100):
        lastedge = {}
        ok = True
        for v, outs in edges.items():
            if v != last:
                lastedge[v] = rng.choice(outs)
        for v in edges:
            if v == last:
                continue
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in lastedge:
                    ok = False
                    break
                seen.add(cur)
                cur = lastedge[cur]
            if not ok:
                break
        if ok:
            break
    else:
        return s
    pools = {v: [e for e in outs] for v, outs in edges.items()}
    for v, e in lastedge.items():
        pools[v].remove(e)
        rng.shuffle(pools[v])
        pools[v].append(e)
    out, cur = [s[0]], s[0]
    for _ in range(len(s) - 1):
        nxt = pools[cur].pop(0)
        out.append(nxt)
        cur = nxt
    return "".join(out)


def prodigal(recs: list[tuple[str, str]]) -> tuple[dict, list[tuple[str, str]]]:
    """Coding density per contig, plus the predicted peptides for the Pfam pass.

    Metagenomic mode, because single-genome training needs far more sequence than 1.3 Mb of
    unrelated fragments and would train on a mixture. Density is reported, not ORF count:
    count is length wearing a different hat.
    """
    import pyrodigal
    gf = pyrodigal.GeneFinder(meta=True)
    per, peptides = [], []
    for name, seq in recs:
        genes = gf.find_genes(seq.encode())
        coding = sum(abs(g.end - g.begin) + 1 for g in genes)
        longest = max((abs(g.end - g.begin) + 1 for g in genes), default=0)
        per.append({"contig": name, "length": len(seq), "n_orf": len(genes),
                    "coding_density": min(1.0, coding / len(seq)) if seq else 0.0,
                    "longest_orf": longest,
                    "longest_orf_frac": min(1.0, longest / len(seq)) if seq else 0.0})
        for i, g in enumerate(genes):
            peptides.append((f"{name}_{i + 1}", g.translate().rstrip("*")))
    return {"per_contig": per}, peptides


def _file_sha(path) -> str:
    """Content hash of an input file. The reuse check needs to know the contigs did not
    change, not merely that the filename did not."""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def _txt(v) -> str:
    """pyhmmer hands back bytes on some versions and str on others."""
    if v is None:
        return ""
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)


def pfam(peptides: list[tuple[str, str]], hmm: Path, threads: int) -> dict:
    """hmmsearch every Pfam-A model against the predicted peptides, at each family's own
    gathering threshold.

    The threshold is Pfam's, not ours, which is the point: a bitscore chosen here would be
    one more operating point to defend, and the GA cutoffs are what the curators use to
    decide a family membership. Models are streamed from the flatfile rather than pressed,
    which costs a little time and saves 4 GB of disk the container does not have.
    """
    import pyhmmer
    if not peptides:
        return {"per_peptide": [], "models_scanned": 0}
    alpha = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=n.encode(), sequence=s).digitize(alpha)
            for n, s in peptides if s]
    hits_by_pep: dict[str, dict] = {}
    n_models = 0
    with pyhmmer.plan7.HMMFile(str(hmm)) as fh:
        for top in pyhmmer.hmmer.hmmsearch(fh, seqs, cpus=threads, bit_cutoffs="gathering"):
            n_models += 1
            acc = _txt(top.query.accession).split(".")[0]
            fam = _txt(top.query.name)
            for hit in top:
                if not hit.included:
                    continue
                key = _txt(hit.name)
                cur = hits_by_pep.get(key)
                if cur is None or hit.score > cur["bitscore"]:
                    hits_by_pep[key] = {"peptide": key, "accession": acc, "family": fam,
                                        "bitscore": float(hit.score),
                                        "evalue": float(hit.evalue)}
    return {"per_peptide": list(hits_by_pep.values()), "models_scanned": n_models}


def rfam(recs: list[tuple[str, str]], cm: Path, cmscan: Path, work: Path,
         threads: int, total_mb: float) -> dict:
    """cmscan against Rfam at the gathering thresholds.

    --rfam and --nohmmonly are Rfam's own recommended cmscan settings for a genome-scale
    scan; -Z fixes the effective database size so E-values are comparable between the two
    piles, which they would not be if each pile set its own.
    """
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "rfam_in.fa"
    write_fasta(fa, recs)
    tbl = work / "rfam.tbl"
    cmd = [str(cmscan), "--cpu", str(threads), "--rfam", "--cut_ga", "--nohmmonly",
           "--oskip", "--fmt", "2", "-Z", f"{total_mb:.6f}",
           "--tblout", str(tbl), "-o", str(work / "rfam.out"), str(cm), str(fa)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return {"failed": r.stderr.strip()[-600:], "seconds": time.time() - t0}
    per: dict[str, dict] = {}
    for line in tbl.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 18:
            continue
        acc, target, ev = f[2], f[3], f[17]
        try:
            ev = float(ev)
        except ValueError:
            continue
        cur = per.get(target)
        if cur is None or ev < cur["evalue"]:
            per[target] = {"contig": target, "rfam_accession": acc,
                           "family": f[1], "evalue": ev}
    return {"per_contig": list(per.values()), "seconds": time.time() - t0,
            "command": " ".join(cmd[:-2]) + " <Rfam.cm> <contigs>"}


def summarise(tag: str, recs, pro, pf, rf) -> dict:
    n = len(recs)
    pep_owner = {}
    for h in pf.get("per_peptide", []):
        pep_owner.setdefault(h["peptide"].rsplit("_", 1)[0], []).append(h)
    rdrp = {c for c, hs in pep_owner.items()
            if any(h["accession"] in RDRP_ACCESSIONS - RT_ACCESSIONS for h in hs)}
    rt = {c for c, hs in pep_owner.items()
          if any(h["accession"] in RT_ACCESSIONS for h in hs)}
    rf_hit = {h["contig"] for h in rf.get("per_contig", [])}
    dens = [p["coding_density"] for p in pro["per_contig"]]
    lens = [p["length"] for p in pro["per_contig"]]
    fams = Counter(h["family"] for h in pf.get("per_peptide", []))
    rfams = Counter(h["family"] for h in rf.get("per_contig", []))
    return {
        "tag": tag,
        "contigs": n,
        "bases": sum(lens),
        "length": mean_ci([float(x) for x in lens]),
        "coding_density": mean_ci(dens),
        "coding_density_above_half": wilson(sum(1 for d in dens if d >= 0.5), n),
        "any_pfam": wilson(len(pep_owner), n),
        "rdrp": wilson(len(rdrp), n),
        "reverse_transcriptase": wilson(len(rt), n),
        "any_rfam": wilson(len(rf_hit), n),
        "top_pfam": fams.most_common(15),
        "top_rfam": rfams.most_common(15),
        "rdrp_contigs": sorted(rdrp),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", type=Path, required=True,
                    help="the residue-only assembly")
    ap.add_argument("--control-contigs", type=Path, nargs="+",
                    default=[ROOT / "data/geo/sample_contigs_identified_norrna.fa",
                             ROOT / "data/geo/sample_contigs_unidentified.fa"],
                    help="pool the length-matched control is drawn from")
    ap.add_argument("--pfam", type=Path, default=ROOT / "data/pfam/Pfam-A.hmm")
    ap.add_argument("--rfam", type=Path, default=ROOT / "data/rfam/Rfam.cm")
    ap.add_argument("--cmscan", type=Path, default=Path("/home/user/tools/inf/bin/cmscan"))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-contigs", type=int, default=0,
                    help="cap each pile, for a smoke run; 0 means all")
    ap.add_argument("--skip-rfam", action="store_true")
    ap.add_argument("--reuse-residue-from", type=Path, default=None,
                    help="an earlier seed's JSON whose residue pile to reuse. The residue "
                         "pile is deterministic - every contig is annotated, nothing is "
                         "sampled - so recomputing it per seed burns hours to produce the "
                         "same five numbers. Only the controls and the shuffle are random. "
                         "The reuse is refused unless the script hash and the contig file "
                         "hash both match, because a pile carried across a code change is "
                         "the mixture this repo re-runs seeds to avoid.")
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/resann"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/residue_contig_annotation")
    args = ap.parse_args()

    for p in [args.contigs, args.pfam] + ([] if args.skip_rfam else [args.rfam, args.cmscan]):
        if not Path(p).exists():
            raise SystemExit(f"missing input: {p}")

    rng = random.Random(args.seed)
    t0 = time.time()
    residue = read_fasta(args.contigs)
    if args.max_contigs:
        rng.shuffle(residue)
        residue = residue[:args.max_contigs]

    # One control per pool, not one pooled control. The sample's identified contigs and its
    # assembled-but-unidentified contigs are different populations - tier 1 and tier 2 - and
    # pooling them would average two controls that are expected to read differently, which
    # is the pair-symmetric-field mistake in another costume.
    # A control pool may be a superset of the query pile - the 200 highest-mass
    # unidentifiable contigs are drawn from the same 46,157 the control is sampled from -
    # and a query sitting inside its own control is the failure this repo has already made
    # once, with a chimaera test that read 100% on both piles for exactly that reason. So
    # every query name is removed from every pool, and how many were removed is recorded.
    query_names = {n for n, _ in residue}
    controls, shortfall, excluded = [], {}, {}
    for p in args.control_contigs:
        if not Path(p).exists():
            continue
        pool = read_fasta(Path(p))
        kept = [r for r in pool if r[0] not in query_names]
        if len(kept) != len(pool):
            excluded[Path(p).name] = len(pool) - len(kept)
        recs, short = length_matched_sample(residue, kept, rng)
        tag = "control_" + Path(p).stem.replace("sample_contigs_", "")
        controls.append((tag, recs))
        if short:
            shortfall[tag] = short
    shuffled = [(f"{n}_shuf", dinuc_shuffle(s, rng)) for n, s in residue]

    args.work.mkdir(parents=True, exist_ok=True)
    total_mb = sum(len(s) for _, s in residue) * 2e-6 or 1e-6
    piles, timings = {}, {}
    reused = None
    todo = [("residue", residue)] + controls + [("dinucleotide_shuffled_residue", shuffled)]
    if args.reuse_residue_from:
        prev = json.loads(Path(args.reuse_residue_from).read_text())
        want = (provenance(__file__)["script_sha256"], _file_sha(args.contigs))
        got = (prev.get("script_sha256"), prev.get("contigs_sha256_16"))
        if want != got:
            raise SystemExit(
                f"refusing to reuse {args.reuse_residue_from}: it was made by script "
                f"{got[0]} over contigs {got[1]}, this run is {want[0]} over {want[1]}")
        piles["residue"] = prev["piles"]["residue"]
        timings["residue"] = prev["timings"]["residue"]
        reused = {"from": str(args.reuse_residue_from), "seed": prev["parameters"]["seed"]}
        todo = [t for t in todo if t[0] != "residue"]
        print(f"[residue] reused from seed {reused['seed']}", flush=True)
    for tag, recs in todo:
        ts = time.time()
        pro, peps = prodigal(recs)
        t_pro = time.time() - ts
        ts = time.time()
        pf = pfam(peps, args.pfam, args.threads)
        t_pf = time.time() - ts
        rf = ({"per_contig": [], "skipped": True} if args.skip_rfam
              else rfam(recs, args.rfam, args.cmscan, args.work / tag, args.threads, total_mb))
        piles[tag] = summarise(tag, recs, pro, pf, rf)
        piles[tag]["peptides"] = len(peps)
        timings[tag] = {"prodigal_s": t_pro, "pfam_s": t_pf,
                        "rfam_s": rf.get("seconds")}
        print(f"[{tag}] n={len(recs)} coding_density_mean="
              f"{piles[tag]['coding_density']['mean']:.3f} "
              f"any_pfam={piles[tag]['any_pfam']['point']:.3%} "
              f"rdrp={piles[tag]['rdrp']['point']:.3%} "
              f"any_rfam={piles[tag]['any_rfam']['point']:.3%}", flush=True)

    out = {
        "parameters": {"contigs": str(args.contigs),
                       "control_contigs": [str(p) for p in args.control_contigs],
                       "seed": args.seed, "threads": args.threads,
                       "max_contigs": args.max_contigs,
                       "pfam_threshold": "per-family gathering (GA) cutoff",
                       "rfam_settings": "cmscan --rfam --cut_ga --nohmmonly, -Z over both piles",
                       "skip_rfam": args.skip_rfam,
                       "length_matching": "log bins, ratio 1.25 from 200 nt, without replacement"},
        "contigs_sha256_16": _file_sha(args.contigs),
        "query_contigs_removed_from_control_pool": excluded,
        "residue_pile_reused": reused,
        "length_match_shortfall": shortfall,
        "piles": piles,
        "timings": timings,
        "wall_seconds": time.time() - t0,
        **provenance(__file__),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / f"residue_contig_annotation_s{args.seed}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"written {dest}")


if __name__ == "__main__":
    main()
