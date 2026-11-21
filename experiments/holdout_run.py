#!/usr/bin/env python3
"""The post-hoc holdout, run at read level against a date-keyed index.

WHAT THIS ASKS
--------------
`results/holdout_candidates/` found the target: a sequence whose first public
deposit is 2023-07-31, assembled from a wastewater RNA library collected
2022-05-10, sitting at 66.5% protein identity to its nearest relative over 80%
of its length and carrying 32,575 reads at 452x in that library. So for the
duration between collection and deposit the thing was in the water and in
nobody's database, and the archive supplies the blind for free.

The question here is the one SPEC s2.1(i) actually promises to answer:
**with the clock turned back, does the per-read readout report the true distance
to the nearest relative?** Not "is it flagged" - the spec's output is
identity-to-closest-relative with coverage beside it, deliberately
un-thresholded - but whether that number, computed from 145 nt reads against an
index that has never seen this organism, lands where the finished 4,147 nt
genome says it should.

That is a measurement with a right answer known in advance and not chosen by
us, which is what `results/spike_holdout/` could never be.

WHAT IS AND IS NOT RUN
----------------------
This is the identification-and-distance half of the cascade, at read level,
against the two indexes SPEC s2 keeps for it: viral RefSeq nucleotide and viral
RefSeq protein. It is NOT the full tier-1 cascade - SILVA, the Kraken2
databases and the sample's own assembly are not staged in this container, and
Kraken2 alone is 16 GB at ~1.3 MB/s. Those tiers change how much background is
*removed*; they do not change what distance is reported for a read that reaches
the protein tier, which is what is being measured. The removal fractions
therefore are not reported here, and the figure for running them is in the
report.

THE HOLDOUT, AND THE TRAP IT AVOIDS
-----------------------------------
Exclusion is keyed to the deposit date of the ORIGINAL GenBank accession behind
each RefSeq entry, recovered from the RefSeq record's own COMMENT ("the
reference sequence is identical to / was derived from X"), 99.9% of 19,641
entries. This matters because a RefSeq entry's own createdate is when the
RefSeq record was built, not when the sequence was first submitted: 3,926 viral
RefSeq entries carry a 2023 createdate while wrapping decades-older sequence
(`results/holdout_candidates/README.md` s3). Keying on that would delete
genuinely old relatives and inflate the reported distance - the dangerous
direction, since it would manufacture the result this experiment is looking
for.

The cutoff is the EARLIEST collection date among the three libraries of the
study, so the holdout is stricter than each library needs rather than looser.

FOUR ARMS, AND WHY EACH IS THERE
--------------------------------
Every number here is unreadable without the other three.

- **candidate**: reads of the held-out organism. The test.
- **leak**: reads of a sequence in the SAME library whose nearest relative is
  99.9% identical and was deposited long before the cutoff. After the holdout
  these MUST still be identified at ~99%. If they are not, the holdout has
  removed more than it should and every other arm reads low for the wrong
  reason. If the candidate arm instead reads ~99%, the holdout leaked.
- **background**: reads from the same library matching none of the study's
  deposited sequences. Says what the index reports for ordinary material, which
  is the "already seen" baseline the candidate must be read against.
- **shuffled**: dinucleotide-preserving shuffles of the candidate reads. Says
  what the identity readout reports on sequence carrying no information. A
  mononucleotide shuffle is not used: stop-codon density is a dinucleotide
  property, so it inflates ORFs and would make the null too easy
  (`results/residue_contig_annotation/`).

READING THE OUTPUT
------------------
Identity is reported only with its alignment length, query coverage and
E-value, and only for significant hits. At 145 nt a protein alignment is ~48
residues, where percent identity is quantised onto n/48 and a bare percentage
is meaningless; this repo has already had to retract an identity claim that was
chance matching over 40 bases.
"""
import argparse, csv, gzip, hashlib, json, os, random, subprocess, sys, time
from collections import Counter, defaultdict

FMT = ("6 qseqid qlen sseqid pident length mismatch gapopen qstart qend "
       "qcovhsp evalue bitscore stitle")


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------------- stage: holdout
def stage_holdout(out, nt_fa, prot_fa, src_map, src_dates, cutoff, prot_parent,
                  blastbin):
    """Rebuild both indexes as they stood on `cutoff`.

    Dropping an entry needs two things to be true, and the asymmetry is
    deliberate: an entry goes only if its source accession is dated AFTER the
    cutoff, so an entry whose date could not be recovered is KEPT. Keeping an
    undated entry can only leave a relative in the index that should have been
    hidden, which makes the candidate look less novel. Dropping it would do the
    opposite and fabricate novelty.
    """
    src = json.load(open(src_map))
    dates = json.load(open(src_dates))
    drop_date, keep, undated = set(), set(), 0
    for ver, acc in src.items():
        d = dates.get(acc)
        if not d:
            undated += 1
            continue
        if d.replace("/", "-") > cutoff:
            drop_date.add(ver.split(".")[0])
    log(f"  {len(drop_date)} RefSeq entries deposited after {cutoff}; "
        f"{undated} undated and therefore KEPT")

    # A protein accession is not its nucleotide accession and the defline does not
    # carry the parent, so filtering proteins by organism name would apply a
    # DIFFERENT rule to the two indexes and the comparison between them would stop
    # meaning anything. Each protein is instead mapped to the RefSeq nucleotide
    # record whose CDS declared it, and both indexes then obey one date rule.
    parent = json.load(open(prot_parent))
    stats = {}
    for label, fa, dbtype in [("nt", nt_fa, "nucl"), ("prot", prot_fa, "prot")]:
        kept = dropped = unmapped = 0
        o = os.path.join(out, f"holdout_{label}.fa")
        with open(o, "w") as w:
            emit = False
            for line in open(fa):
                if line.startswith(">"):
                    acc = line[1:].split()[0]
                    if label == "prot":
                        p = parent.get(acc) or parent.get(acc.split(".")[0])
                        if p is None:
                            unmapped += 1
                            bad = False          # unmapped is kept, as above
                        else:
                            bad = p.split(".")[0] in drop_date
                    else:
                        bad = acc.split(".")[0] in drop_date
                    emit = not bad
                    kept += emit
                    dropped += not emit
                if emit:
                    w.write(line)
        if unmapped:
            log(f"  {unmapped} {label} entries had no parent record and were kept")
        stats[label] = dict(kept=kept, dropped=dropped, unmapped_kept=unmapped)
        subprocess.run([os.path.join(blastbin, "makeblastdb"), "-in", o,
                        "-dbtype", dbtype, "-out", os.path.join(out, f"db_holdout_{label}")],
                       check=True, capture_output=True)
        log(f"  holdout {label}: kept {kept}, dropped {dropped}")
    json.dump(sorted(drop_date), open(os.path.join(out, "dropped_accessions.json"), "w"))
    json.dump(dict(cutoff=cutoff, dropped_by_date=len(drop_date), undated_kept=undated,
                   note=("the held-out organism has no entry of its own in viral RefSeq, "
                         "so the holdout is purely the date filter"),
                   per_index=stats),
              open(os.path.join(out, "holdout_manifest.json"), "w"), indent=1)
    return stats


# --------------------------------------------------------------- stage: reads
def dinuc_shuffle(seq, rng, tries=200):
    """Shuffle preserving dinucleotide composition (Altschul-Erikson).

    A mononucleotide shuffle destroys stop-codon density, which is a
    dinucleotide property, so it inflates open reading frames and makes the null
    easier to beat than it should be.

    The naive version of this walk - shuffle each vertex's out-edges, then walk -
    dead-ends on most real reads, because an Eulerian path only exists if, for
    every vertex but the terminal one, the LAST edge taken out of it leads
    (transitively) to the terminal vertex. Altschul and Erikson's fix is to
    choose those last edges first and reject the draw unless they form a tree
    rooted at the terminal vertex. The first version of this function omitted
    that and fell back to returning the input sequence unchanged on failure:
    70.4% of the shuffled arm was then unshuffled candidate reads, and the null
    read the same as the test. It now returns None on failure and the caller
    counts those, because a control that silently becomes a copy of the test is
    worse than no control.
    """
    if len(seq) < 4:
        return None
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    term = seq[-1]
    verts = [v for v in edges if v != term]
    chosen = None
    for _ in range(tries):
        pick = {v: rng.choice(edges[v]) for v in verts}
        ok = True
        for v in verts:                        # every vertex must reach term
            seen, cur = set(), v
            while cur != term:
                if cur in seen or cur not in pick:
                    ok = False
                    break
                seen.add(cur)
                cur = pick[cur]
            if not ok:
                break
        if ok:
            chosen = pick
            break
    if chosen is None:
        return None
    order = {}
    for v, lst in edges.items():
        rest = list(lst)
        if v in chosen:
            rest.remove(chosen[v])
            rng.shuffle(rest)
            rest.append(chosen[v])             # the reserved last edge
        else:
            rng.shuffle(rest)
        order[v] = rest
    at = {v: 0 for v in order}
    cur, outl = seq[0], [seq[0]]
    for _ in range(len(seq) - 1):
        lst = order.get(cur)
        if lst is None or at[cur] >= len(lst):
            return None
        nxt = lst[at[cur]]
        at[cur] += 1
        outl.append(nxt)
        cur = nxt
    return "".join(outl) if len(outl) == len(seq) else None


def shuffle_arm(take, rng):
    """Shuffle every read of the candidate subsample, and say what was lost.

    The shuffled arm is the SAME subsample as the candidate arm, not an
    independent draw from the same pool: the control's job is to say what the
    identity readout reports on these reads carrying no information, which a
    different set of reads cannot answer. Reads whose shuffle fails, or whose
    composition is degenerate enough that the shuffle returns the input, are
    DROPPED and counted rather than passed through.
    """
    outl, failed, identical = [], 0, 0
    for nm, sq in take:
        sh = dinuc_shuffle(sq, rng)
        if sh is None:
            failed += 1
        elif sh == sq:
            identical += 1
        else:
            outl.append((nm + "_shuf", sh))
    return outl, failed, identical


def stage_reads(out, candfa, fastq, candidate, leak, n_reads, seeds):
    """Pull each arm's reads out of the library by mapping, then subsample.

    Arms are defined by which deposited sequence a read's PRIMARY alignment goes
    to, so a read cannot be in two arms, and the background arm is reads that map
    to none of them - which is most of the library.
    """
    import mappy
    idx = mappy.Aligner(candfa, preset="sr", n_threads=os.cpu_count() or 4)
    pools = {"candidate": [], "leak": [], "background": []}
    nreads = 0
    for name, seq, qual in mappy.fastx_read(fastq):
        nreads += 1
        best = None
        for h in idx.map(seq):
            if h.is_primary and (best is None or h.mlen > best.mlen):
                best = h
        if best is None:
            # the background arm is huge; reservoir-free thinning keeps memory flat
            if len(pools["background"]) < n_reads * 20:
                pools["background"].append((name, seq))
        elif best.ctg == candidate:
            pools["candidate"].append((name, seq))
        elif best.ctg == leak:
            pools["leak"].append((name, seq))
    log(f"  {nreads} reads scanned; " +
        ", ".join(f"{k}={len(v)}" for k, v in pools.items()))

    counts, cand_take, shuffle_losses = {}, {}, {}
    for seed in seeds:
        rng = random.Random(seed)
        for arm, pool in pools.items():
            take = pool if len(pool) <= n_reads else rng.sample(pool, n_reads)
            if arm == "candidate":
                cand_take[seed] = take
            with open(os.path.join(out, f"reads_{arm}_s{seed}.fa"), "w") as f:
                for nm, sq in take:
                    f.write(f">{nm}\n{sq}\n")
            counts[f"{arm}_s{seed}"] = len(take)
        shuf, nfail, nsame = shuffle_arm(cand_take[seed], rng)
        with open(os.path.join(out, f"reads_shuffled_s{seed}.fa"), "w") as f:
            for nm, sq in shuf:
                f.write(f">{nm}\n{sq}\n")
        counts[f"shuffled_s{seed}"] = len(shuf)
        shuffle_losses[f"s{seed}"] = dict(failed=nfail, returned_input=nsame)
    json.dump(dict(reads_scanned=nreads, pool_sizes={k: len(v) for k, v in pools.items()},
                   subsampled=counts, shuffle_losses=shuffle_losses,
                   fastq=fastq, candidate=candidate, leak=leak),
              open(os.path.join(out, "reads_manifest.json"), "w"), indent=1)
    return counts


# -------------------------------------------------------------- stage: search
def stage_search(out, blastbin, dbs, arms, seeds, evalue="1e-5"):
    """Search every arm against every index, nucleotide and protein.

    The E-value floor is 1e-5 and it is doing real work: BLAST reports to E=10 by
    default, and this project has already published and retracted a distance that
    was chance matching at a median E-value of 0.24. Both indexes are searched
    with the full and the held-out version, because the comparison between those
    two IS the experiment - the same reads, the same code, one index with the
    organism's relatives of the future removed.
    """
    jobs = []
    for arm in arms:
        for seed in seeds:
            q = os.path.join(out, f"reads_{arm}_s{seed}.fa")
            if not os.path.exists(q):
                continue
            for dbname, (dbpath, prog) in dbs.items():
                tag = f"{arm}_s{seed}_{dbname}"
                o = os.path.join(out, f"hits_{tag}.tsv")
                if os.path.exists(o) and os.path.getsize(o) >= 0 and \
                        os.path.exists(o + ".done"):
                    continue
                jobs.append((tag, q, dbpath, prog, o))
    for i, (tag, q, dbpath, prog, o) in enumerate(jobs, 1):
        cmd = [os.path.join(blastbin, prog), "-query", q, "-db", dbpath,
               "-outfmt", FMT, "-max_target_seqs", "5", "-evalue", evalue,
               "-num_threads", str(os.cpu_count() or 4)]
        if prog == "blastn":
            cmd += ["-task", "blastn"]          # most sensitive nucleotide setting
        t0 = time.time()
        with open(o, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)
        open(o + ".done", "w").close()
        log(f"  [{i}/{len(jobs)}] {tag}: {sum(1 for _ in open(o))} hsps "
            f"in {time.time()-t0:.0f}s")


# -------------------------------------------------------------- stage: report
def stage_report(out, arms, seeds, dbs, truth, dropped=frozenset()):
    """Per arm, per index: what fraction of reads is identified at all, and what
    distance the identified ones report.

    Identification rate and reported distance are kept apart on purpose. They
    answer different questions and in this project they have moved in opposite
    directions before: a search can identify a read and report a distance that
    says "already seen", which is why `truth` - the distance measured from the
    finished genome - is carried in the output rather than left to the reader.
    """
    rows = []
    for arm in arms:
        for dbname in dbs:
            for seed in seeds:
                p = os.path.join(out, f"hits_{arm}_s{seed}_{dbname}.tsv")
                qp = os.path.join(out, f"reads_{arm}_s{seed}.fa")
                if not (os.path.exists(p) and os.path.exists(qp)):
                    continue
                nq = sum(1 for l in open(qp) if l.startswith(">"))
                best = {}
                for r in csv.reader(open(p), delimiter="\t"):
                    if len(r) < 13:
                        continue
                    q, bits = r[0], float(r[11])
                    if q not in best or bits > float(best[q][11]):
                        best[q] = r
                pid = sorted(float(r[3]) for r in best.values())
                aln = sorted(int(r[4]) for r in best.values())
                cov = sorted(float(r[9]) for r in best.values())
                med = lambda v: (v[len(v) // 2] if v else "")
                rows.append(dict(
                    arm=arm, index=dbname, seed=seed, reads=nq,
                    identified=len(best),
                    identified_pct=round(100 * len(best) / nq, 3) if nq else "",
                    median_identity=med(pid), median_aln_len=med(aln),
                    median_qcov=med(cov),
                    p10_identity=(pid[len(pid) // 10] if pid else ""),
                    p90_identity=(pid[9 * len(pid) // 10] if pid else ""),
                    future_subject_pct=round(100 * sum(
                        1 for r in best.values()
                        if r[2].split(".")[0] in dropped) / max(len(best), 1), 3),
                    truth_identity=truth.get((arm, dbname.split("_")[1]), "")))
    rows.sort(key=lambda r: (r["arm"], r["index"], r["seed"]))
    keys = list(rows[0].keys()) if rows else []
    with open(os.path.join(out, "readlevel.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # spread across seeds, reported as an observed min-max and never as a
    # dispersion estimate: five seeds give an honest range and a useless variance
    agg = defaultdict(list)
    for r in rows:
        agg[(r["arm"], r["index"])].append(r)
    out_rows = []
    for (arm, dbname), rs in sorted(agg.items()):
        ip = [x["identified_pct"] for x in rs if x["identified_pct"] != ""]
        mi = [x["median_identity"] for x in rs if x["median_identity"] != ""]
        out_rows.append(dict(arm=arm, index=dbname, seeds=len(rs),
                             reads_per_seed=rs[0]["reads"],
                             identified_pct_min=min(ip) if ip else "",
                             identified_pct_max=max(ip) if ip else "",
                             median_identity_min=min(mi) if mi else "",
                             median_identity_max=max(mi) if mi else "",
                             median_aln_len=rs[0]["median_aln_len"],
                             truth_identity=rs[0]["truth_identity"]))
    with open(os.path.join(out, "SUMMARY.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    for r in out_rows:
        log(f"  {r['arm']:11} {r['index']:14} id {r['identified_pct_min']}-"
            f"{r['identified_pct_max']}%  identity {r['median_identity_min']}-"
            f"{r['median_identity_max']}  truth {r['truth_identity']}")
    return out_rows


def manifest(out, args, timings):
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    json.dump(dict(script=os.path.basename(__file__),
                   script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                   git_sha=sha + ("-dirty" if dirty else ""),
                   argv=vars(args), timings_s=timings,
                   run_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
              open(os.path.join(out, f"manifest_s{min(args.seeds)}.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="results/holdout_run")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--stages", default="holdout,reads,search,report")
    ap.add_argument("--blastbin", required=True)
    ap.add_argument("--nt-fa", default="")
    ap.add_argument("--prot-fa", default="")
    ap.add_argument("--src-map", default="")
    ap.add_argument("--src-dates", default="")
    ap.add_argument("--full-nt-db", default="")
    ap.add_argument("--full-prot-db", default="")
    ap.add_argument("--cutoff", default="2022-05-10",
                    help="earliest collection date of the study's libraries")
    ap.add_argument("--prot-parent", default="",
                    help="protein accession -> parent RefSeq nucleotide record, from the "
                         "RefSeq release flatfile's CDS /protein_id entries")
    ap.add_argument("--candfa", default="results/holdout_candidates/candidates.fa")
    ap.add_argument("--fastq", default="")
    ap.add_argument("--candidate", default="OQ722337.1")
    ap.add_argument("--leak", default="OQ722330.1")
    ap.add_argument("--n-reads", type=int, default=2000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    stages = args.stages.split(",")
    timings, t0 = {}, time.time()

    if "holdout" in stages:
        t = time.time()
        stage_holdout(args.outdir, args.nt_fa, args.prot_fa, args.src_map,
                      args.src_dates, args.cutoff, args.prot_parent, args.blastbin)
        timings["holdout"] = round(time.time() - t, 1)
    if "reads" in stages:
        t = time.time()
        stage_reads(args.outdir, args.candfa, args.fastq, args.candidate, args.leak,
                    args.n_reads, args.seeds)
        timings["reads"] = round(time.time() - t, 1)

    dbs = {
        "full_nt": (args.full_nt_db, "blastn"),
        "holdout_nt": (os.path.join(args.outdir, "db_holdout_nt"), "blastn"),
        "full_prot": (args.full_prot_db, "blastx"),
        "holdout_prot": (os.path.join(args.outdir, "db_holdout_prot"), "blastx"),
    }
    arms = ["candidate", "leak", "background", "shuffled"]
    if "search" in stages:
        t = time.time()
        stage_search(args.outdir, args.blastbin, dbs, arms, args.seeds)
        timings["search"] = round(time.time() - t, 1)
    if "report" in stages:
        t = time.time()
        # distances measured from the finished genomes in
        # results/holdout_candidates/shortlist.tsv, which is what the read-level
        # readout is being asked to reproduce
        # truth is per LEVEL, not per index: the finished 4,147 nt genome reads
        # 69.72% nucleotide identity and 66.45% protein identity to its nearest
        # relative (results/holdout_candidates/shortlist.tsv), and the leak arm's
        # genome reads 99.92% / 99.81%. Keying these by index name rather than by
        # level was a bug that compared a protein readout against a nucleotide
        # truth.
        truth = {("candidate", "nt"): 69.72, ("candidate", "prot"): 66.45,
                 ("leak", "nt"): 99.92, ("leak", "prot"): 99.81}
        # which best hits in the FULL index are entries the holdout removed: the
        # direct measure of how much knowledge from after the collection date the
        # unrestricted index is contributing, rather than an inference from it
        dropped = set()
        hm = os.path.join(args.outdir, "dropped_accessions.json")
        if os.path.exists(hm):
            dropped = set(json.load(open(hm)))
        stage_report(args.outdir, arms, args.seeds, dbs, truth, dropped)
        timings["report"] = round(time.time() - t, 1)

    timings["total"] = round(time.time() - t0, 1)
    manifest(args.outdir, args, timings)


if __name__ == "__main__":
    main()
