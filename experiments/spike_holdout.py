"""The positive control: spike a held-out virus into the real background and see where it
lands in the cascade, and how it scores on every field of the ladder.

WHY THIS IS THE EXPERIMENT THAT MATTERS. Everything measured so far is background-only, so
we know what the cascade REMOVES and nothing about whether it keeps the thing we care
about. A spike whose reads are removed at tier 1 by a public database is a false negative
of the whole architecture, and that is precisely the failure nobody checks.

THE DIVERGENCE LADDER, instead of one "novel" spike. A simulated novel family is a fiction
however it is built, so this script does not pretend to produce one. It takes a real viral
genome and walks it away from itself under a substitution model at several rates, which
turns the question into a measurable one: HOW DIVERGENT MUST A VIRUS BE BEFORE THE CASCADE
STOPS IDENTIFYING IT, and what does it look like on the ladder once it does? At 0%
divergence the spike is a known virus and should be removed - that is the control that the
cascade works at all. The rate at which survival crosses from 0 to 1 is the boundary
between "discovery" and "already known" for this pipeline, measured rather than asserted.

Substitutions are applied uniformly over the nucleotide sequence. That is a crude model of
real divergence - real viral families differ by selection-shaped, codon-structured change,
with synonymous sites far more variable - so the rates here are NOT equivalent to a family
distance, and the honest reading is ordinal: more divergence, less identifiable. Stated
here rather than in a footnote because the temptation to read "20% divergent" as "a new
genus" is strong.

ABUNDANCE, and why the headline rate is not simulated directly. At 1 in a million a spike
contributes about one read pair per million, so the whole 8.7M-pair run holds ~9 pairs and
a single 100,000-pair block holds 0.1. Simulating that directly measures nothing: the
sample size cannot resolve the quantity, which is the error this project made once already
(PR #11, where 1-per-million thresholds came from zero counts in a 20,000-read run). So
the per-read behaviour is measured at a spike rate high enough to resolve it, and the
1-in-a-million statement is then arithmetic on top of per-read survival, reported as
arithmetic. The one place the real rate matters is the library-recurrence field, which
depends on how many copies of the spike are present; there the sketch is built with exactly
the number of spike pairs that 1-in-a-million implies, so that field is measured at the
real abundance rather than at the convenient one.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import statistics as st
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from berlin_io import read_pairs, revcomp
from subtract_known import (_mm_index, artefact_fields, dir_bytes, git_sha, kraken2,
                            map_filter, provenance, shannon_entropy, spread)
from loud import Preflight, require_output, run_tool
from residue_report import (LibraryAbundance, dinucleotide_bias, folding, frame_evidence,
                            repetitiveness, _sample_of)

ROOT = Path(__file__).resolve().parents[1]


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    """Uniform substitution at `rate` per site, no indels.

    No indels on purpose: an indel frameshifts everything downstream, which would make the
    coding fields collapse for reasons that have nothing to do with how far away a relative
    is. Indel effects on peptide-level features were measured separately in PR #11.
    """
    out = []
    alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    for b in seq:
        if b in alt and rng.random() < rate:
            out.append(rng.choice(alt[b]))
        else:
            out.append(b)
    return "".join(out)


def simulate_pairs(genome: str, n: int, rng: random.Random, read_len: int,
                   insert_quantiles: list[tuple[float, int]], err1: float, err2: float):
    """Paired reads off a genome with this library's measured insert and error profile.

    The insert distribution is the one measured on SRR23973193 itself (PR #10: median 194,
    p75 246, p95 363, p99 485), not a round number, because 53% of this library's pairs
    overlap their mate and that changes how much independent sequence a pair carries.
    Error rates are Berlin's own measured means, 0.00089 on R1 and 0.00117 on R2 (PR #11),
    rather than a nominal Q30.
    """
    alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}

    def add_errors(s: str, e: float) -> str:
        if e <= 0:
            return s
        out = list(s)
        for i, b in enumerate(out):
            if b in alt and rng.random() < e:
                out[i] = rng.choice(alt[b])
        return "".join(out)

    pairs = []
    for i in range(n):
        u = rng.random()
        ins = insert_quantiles[-1][1]
        for q, v in insert_quantiles:
            if u <= q:
                ins = v
                break
        ins = max(ins, read_len + 1)
        if len(genome) <= ins:
            ins = len(genome) - 1
        start = rng.randrange(0, max(len(genome) - ins, 1))
        frag = genome[start:start + ins]
        if rng.random() < 0.5:
            frag = revcomp(frag)
        r1 = add_errors(frag[:read_len], err1)
        r2 = add_errors(revcomp(frag)[:read_len], err2)
        pairs.append((f"spike{i}", r1, r2))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--genome", type=Path, nargs="+", required=True,
                    help="fasta of the held-out genome(s) to spike")
    ap.add_argument("--divergence", type=float, nargs="+",
                    default=[0.0, 0.05, 0.10, 0.20, 0.30])
    ap.add_argument("--spike-pairs", type=int, default=2000,
                    help="pairs simulated per divergence rung: enough to resolve a "
                         "survival fraction to ~1%, unlike the real 1-in-a-million rate")
    ap.add_argument("--read-len", type=int, default=107)
    ap.add_argument("--err1", type=float, default=0.00089)
    ap.add_argument("--err2", type=float, default=0.00117)
    ap.add_argument("--background-pairs", type=int, default=100_000,
                    help="real pairs used for the background score distributions")
    ap.add_argument("--library-pairs", type=int, default=400_000,
                    help="real pairs counted into the abundance sketch")
    ap.add_argument("--abundance-only", action="store_true",
                    help="run only the library-recurrence arm, which is the one arm "
                         "measured at the real 1-in-a-million rate and therefore the one "
                         "worth re-running on its own")
    ap.add_argument("--real-rate", type=float, default=1e-6,
                    help="the abundance the headline is about; sets how many spike pairs "
                         "go into the abundance sketch")
    ap.add_argument("--r1", type=Path, default=ROOT / "data/berlin/SRR23973193_1.fastq.gz")
    ap.add_argument("--r2", type=Path, default=ROOT / "data/berlin/SRR23973193_2.fastq.gz")
    ap.add_argument("--silva", type=Path, default=ROOT / "data/silva/silva_dna.fa")
    ap.add_argument("--identified", type=Path,
                    default=ROOT / "data/geo/sample_contigs_identified.fa")
    ap.add_argument("--unidentified", type=Path,
                    default=ROOT / "data/geo/sample_contigs_unidentified.fa")
    ap.add_argument("--cross-sample-ref", type=Path,
                    default=ROOT / "data/geo/all_contigs_1kb.fa")
    ap.add_argument("--k2-dbs", type=Path, nargs="+",
                    default=[ROOT / "data/k2_viral", ROOT / "data/k2_pluspf16"])
    ap.add_argument("--classify-bin", type=Path,
                    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify"))
    ap.add_argument("--k2-confidence", type=float, default=0.1)
    ap.add_argument("--dmnd", type=Path, nargs="*",
                    default=[ROOT / "data/prot/viral.dmnd", ROOT / "data/prot/sprot.dmnd"])
    ap.add_argument("--diamond-bin", type=Path, default=Path("/home/user/tools/diamond"))
    ap.add_argument("--diamond-bitscore", type=float, default=50.0)
    ap.add_argument("--min-mlen", type=int, default=40)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", type=Path, default=Path("/tmp/claude-0/wk/spike"))
    ap.add_argument("--out", type=Path, default=ROOT / "results/spike_holdout")
    ap.add_argument("--preflight", action="store_true",
                    help="check every input, index and tool, then exit: 0 if the run "
                         "could start, 2 if not")
    args = ap.parse_args()

    # The ladder is the run's only positive control, so an index that silently drops out
    # of it makes the spike survive tier 1 more often and reads as lower sensitivity of
    # the TOOL. Every index it was asked for has to be there.
    pf = Preflight("spike_holdout.py")
    for g in args.genome:
        pf.file(g, "--genome")
    pf.file(args.r1, "--r1")
    pf.file(args.r2, "--r2")
    if not args.abundance_only:
        pf.file(args.silva, "--silva")
        pf.file(args.identified, "--identified")
        pf.file(args.unidentified, "--unidentified")
        pf.file(args.cross_sample_ref, "--cross-sample-ref")
        for db in args.k2_dbs:
            pf.kraken2_db(db)
        if args.k2_dbs:
            pf.exe(args.classify_bin, "--classify-bin")
        for db in args.dmnd:
            pf.file(db, "--dmnd")
        if args.dmnd:
            pf.exe(args.diamond_bin, "--diamond-bin")
    pf.done(args.preflight)

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    t_start = time.perf_counter()

    genomes = {}
    for g in args.genome:
        name, buf = None, []
        for line in open(g):
            if line.startswith(">"):
                if name:
                    genomes[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
        if name:
            genomes[name] = "".join(buf)

    insert_q = [(0.50, 194), (0.75, 246), (0.95, 363), (0.99, 485), (1.0, 600)]

    # ---- background reads, for the distributions every spike score is read against
    bg = ([] if args.abundance_only else
          list(read_pairs(args.r1, args.r2, args.background_pairs, 14,
                          (args.seed - 1) * 1_500_000)))

    # ---- indexes, built once
    idx = {}
    if args.abundance_only:
        idx = {}
    else:
        idx["rrna"] = _mm_index(args.silva, args.work / "silva.mmi", args.threads)
    if not args.abundance_only:
        for tag, ref in (("selfasm_identified", args.identified),
                         ("selfasm_unidentified", args.unidentified),
                         ("cross_sample", args.cross_sample_ref)):
            idx[tag] = _mm_index(ref, args.work / f"{ref.stem}.mmi", args.threads)

    rungs: list[dict] = []
    if not args.abundance_only:
        # ---- every rung's reads are generated first and passed through each index ONCE.
        # Loading a 16 GB Kraken2 hash per rung dominated the runtime and did nothing for the
        # measurement: the per-read answer does not depend on what else is in the file.
        all_pairs: list[tuple[str, str, str]] = []
        key_of: dict[str, tuple[float, str]] = {}
        for div in args.divergence:
            for gname, gseq in genomes.items():
                mut = mutate(gseq, div, random.Random(args.seed * 1000 + int(div * 1000)))
                tag = f"d{int(div * 1000)}_{gname}"
                for nm, a, b in simulate_pairs(mut, args.spike_pairs, rng, args.read_len,
                                               insert_q, args.err1, args.err2):
                    key = f"{tag}|{nm}"
                    all_pairs.append((key, a, b))
                    key_of[key] = (div, gname)
        seqs = {k: (a, b) for k, a, b in all_pairs}

        hit: dict[str, set[str]] = {}
        for tag in ("rrna", "selfasm_identified", "selfasm_unidentified", "cross_sample"):
            if tag in idx:
                hit[tag] = set(map_filter(idx[tag], seqs, args.min_mlen))
        fq1, fq2 = args.work / "sp1.fq", args.work / "sp2.fq"
        for path, i in ((fq1, 1), (fq2, 2)):
            with open(path, "w") as h:
                for nm, a, b in all_pairs:
                    sq = a if i == 1 else b
                    h.write(f"@{nm}\n{sq}\n+\n{'I' * len(sq)}\n")
        for db in args.k2_dbs:
            calls, _w, _r = kraken2(db, fq1, fq2, args.classify_bin, args.threads,
                                    args.work / "k2.out", args.k2_confidence)
            if len(calls) != len(all_pairs):
                raise SystemExit(f"Kraken2 {db.name} returned {len(calls)} calls for "
                                 f"{len(all_pairs)} pairs")
            hit[db.name] = {k for k, tx in calls.items() if tx != 0}
        fa = args.work / "spike.fa"
        with open(fa, "w") as h:
            for nm, a, b in all_pairs:
                h.write(f">{nm}/1\n{a}\n>{nm}/2\n{b}\n")
        prot_best: dict[str, dict[str, tuple[float, float]]] = {}
        for db in args.dmnd:
            tsv = args.work / f"{db.stem}.tsv"
            tsv.unlink(missing_ok=True)
            run_tool([str(args.diamond_bin), "blastx", "--db", str(db),
                      "--query", str(fa), "--out", str(tsv), "--threads",
                      str(args.threads), "--outfmt", "6", "qseqid", "sseqid",
                      "pident", "bitscore", "--max-target-seqs", "1",
                      "--evalue", "10", "--sensitive", "--quiet"],
                     f"DIAMOND blastx against {db.name}")
            require_output(tsv, f"DIAMOND {db.name}", allow_empty=True)
            best: dict[str, tuple[float, float]] = {}
            with open(tsv) as h:
                for line in h:
                    q, _sid, pid, bits = line.rstrip("\n").split("\t")
                    nm = q.rsplit("/", 1)[0]
                    bits = float(bits)
                    if bits > best.get(nm, (0.0, 0.0))[0]:
                        best[nm] = (bits, float(pid))
            prot_best[db.stem] = best
            hit[f"prot_{db.stem}"] = {k for k, v in best.items()
                                      if v[0] >= args.diamond_bitscore}
            tsv.unlink(missing_ok=True)

        # tier 1 is every index the cascade SUBTRACTS; the two contig references it only
        # observes are reported separately, because subtracting them would make the
        # cross-sample and tier-2 fields zero by construction.
        tier1_all = set()
        for tag, sset in hit.items():
            if tag not in ("selfasm_unidentified", "cross_sample"):
                tier1_all |= sset

        rungs = []
        for div in args.divergence:
            rung = {"divergence": div, "genomes": {}}
            for gname in genomes:
                mine = [k for k, (d, g) in key_of.items() if d == div and g == gname]
                n = len(mine)
                mset = set(mine)
                survival = {
                    "pairs": n,
                    "removed_by_tier": {t: len(s & mset) / n for t, s in hit.items()},
                    "tier1_removed": len(tier1_all & mset) / n,
                    "survives_tier1": 1 - len(tier1_all & mset) / n,
                    "survivor_maps_unidentified_contig":
                        len((hit.get("selfasm_unidentified", set()) & mset) - tier1_all) / n,
                    "protein_bitscore_of_best_hit": _dist(
                        [v[0] for db in prot_best for k, v in prot_best[db].items()
                         if k in mset]),
                    "protein_pident_of_best_hit": _dist(
                        [v[1] for db in prot_best for k, v in prot_best[db].items()
                         if k in mset]),
                }
                surv = [k for k in mine if k not in tier1_all]
                fields = (_ladder([seqs[k] for k in surv[:500]], rng) if surv else None)
                rung["genomes"][gname] = {"survival": survival, "ladder_on_survivors": fields}
            rungs.append(rung)

    bg_fields = _ladder([(r[1], r[2]) for r in bg[:500]], rng) if bg else None

    # ---- abundance at the REAL rate, for EVERY spiked genome.
    #
    # This arm is measured at 1 in a million rather than at a convenient rate, because
    # recurrence is the one field whose value depends on how many copies are present.
    # It is also the arm that caught a mistake in the spike's design: two of the three
    # genomes used here are already IN this sample (pepper mild mottle virus assembles at
    # 285x coverage, human astrovirus as a 6,823 nt contig at 11.7x), so their reads
    # recur in the library for a real biological reason and not because the sketch is
    # colliding. Human rhinovirus 1 does not appear in this sample's assembly at all and
    # is therefore the clean holdout; reading the three side by side is what makes the
    # field interpretable, which is why all three are measured rather than the first.
    spike_in_library = max(int(round(args.real_rate * args.library_pairs)), 1)
    sk = LibraryAbundance(bits=28)
    for _nm, a, b, _q1, _q2 in read_pairs(args.r1, args.r2, args.library_pairs, 14, 0):
        sk.add_seq(a)
        sk.add_seq(b)
    bg_counts = []
    for _nm, a, b, _q1, _q2 in read_pairs(args.r1, args.r2, 2000, 14, 0):
        cs = sorted(sk.counts(a))
        if cs:
            bg_counts.append(cs[len(cs) // 2])
    per_genome = {}
    for gname, gseq in genomes.items():
        real_spike = simulate_pairs(gseq, spike_in_library, rng, args.read_len,
                                    insert_q, args.err1, args.err2)
        for _nm, a, b in real_spike:
            sk.add_seq(a)
            sk.add_seq(b)
        med = []
        for _nm, a, b in real_spike:
            cs = sorted(sk.counts(a) + sk.counts(b))
            if cs:
                med.append(cs[len(cs) // 2])
        per_genome[gname] = spread(med) if med else None

    out = {
        "seed": args.seed,
        "genomes": {k: len(v) for k, v in genomes.items()},
        "spike_pairs_per_rung": args.spike_pairs,
        "insert_quantiles": insert_q,
        "error_rates": [args.err1, args.err2],
        "rungs": rungs,
        "background_ladder": bg_fields,
        "abundance_at_real_rate": {
            "real_rate": args.real_rate,
            "library_pairs_in_sketch": args.library_pairs,
            "spike_pairs_in_sketch": spike_in_library,
            "spike_median_31mer_count_per_genome": per_genome,
            "background_median_31mer_count": spread(bg_counts) if bg_counts else None,
            "sketch_saturation": sk.saturation(),
            "note": "the sketch block is the SAME 400,000 pairs for every seed by "
                    "design - it stands for 'the library' - so seeds vary the spike "
                    "reads and not the background, and the background number here is "
                    "one measurement rather than five. "
                    "the count-min sketch only overstates, so a read it calls rare is "
                    "rare for certain; a genome already present in the library recurs "
                    "for a biological reason, which is why all three are reported",
        },
        "one_in_a_million_arithmetic": {
            "run_pairs": 8_731_224,
            "expected_spike_pairs_in_run": 8_731_224 * args.real_rate,
            "note": "per-read survival from the rungs above times this expectation gives "
                    "the expected number of surviving spike pairs in one run; nothing "
                    "here was measured at this rate directly",
        },
        "wall_seconds": time.perf_counter() - t_start,
    }
    (args.out / f"spike_s{args.seed}.json").write_text(json.dumps(out, indent=2))
    (args.out / f"manifest_s{args.seed}.json").write_text(json.dumps({
        "script": "experiments/spike_holdout.py",
        **provenance(__file__),
        # Every index in `parameters` was required to exist before the run started, so
        # this list is what tier 1 actually consisted of, not what was asked for.
        "tier1_indexes_consulted": sorted(hit) if not args.abundance_only else [],
        "parameters": {k: (str(v) if isinstance(v, Path) else
                           [str(x) for x in v] if isinstance(v, list) else v)
                       for k, v in vars(args).items()},
        "wall_seconds": time.perf_counter() - t_start,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "index_bytes": {str(p): dir_bytes(p) for p in
                        [args.silva, args.identified, args.unidentified] + list(args.dmnd)
                        if p.exists()},
    }, indent=2))

    for rung in rungs:
        for g, d in rung["genomes"].items():
            s = d["survival"]
            print(f"div {rung['divergence']:.2f}  {g:<14} tier1 removes "
                  f"{s['tier1_removed'] * 100:6.2f}%  survives {s['survives_tier1'] * 100:6.2f}%"
                  f"  of which on an unidentified contig "
                  f"{s['survivor_maps_unidentified_contig'] * 100:5.2f}%"
                  f"  best prot bits median "
                  f"{(s['protein_bitscore_of_best_hit'] or {}).get('median')}")


def _ladder(pairs, rng) -> dict | None:
    """The ladder fields, summarised, on a list of (r1, r2). None when there is nothing
    to summarise, which happens when a seed's block starts past the end of the file."""
    if not pairs:
        return None
    gc, ent, orf, ident, dup, cpg, upa, z = [], [], [], [], [], [], [], []
    for s1, s2 in pairs:
        art = artefact_fields(s1, s2)
        gc.append(art["gc"])
        ent.append(art["entropy3"])
        fr = frame_evidence(s1 + "NN" + s2)
        orf.append(fr["longest_orf_aa"])
        ident.append(1 if fr["frame_identifiable"] else 0)
        rep = repetitiveness(s1, s2)
        dup.append(rep["kmer16_duplicate_frac"])
        dn = dinucleotide_bias(s1 + s2)
        if dn["CpG_oe"] is not None:
            cpg.append(dn["CpG_oe"])
        if dn["UpA_oe"] is not None:
            upa.append(dn["UpA_oe"])
        f = folding(s1, rng, shuffles=8)
        if f["z"] is not None:
            z.append(f["z"])
    return {"n": len(pairs), "gc": spread(gc), "entropy3": spread(ent),
            "longest_orf_aa": spread(orf),
            "frame_identifiable_frac": sum(ident) / len(ident) if ident else None,
            "kmer16_duplicate_frac": spread(dup),
            "CpG_oe": spread(cpg) if cpg else None,
            "UpA_oe": spread(upa) if upa else None,
            "fold_z": spread(z) if z else None}


def _dist(xs):
    xs = sorted(xs)
    if not xs:
        return None
    return {"n": len(xs), "min": xs[0], "median": xs[len(xs) // 2], "max": xs[-1]}


if __name__ == "__main__":
    main()
