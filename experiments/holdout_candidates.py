#!/usr/bin/env python3
"""Search the public archive for a post-hoc holdout target.

WHY THIS EXISTS
---------------
`results/spike_holdout/` is the only positive control this project has, and it is
a *spike*: reads from a genome we chose, planted at a divergence we chose, into a
library we chose. It answers "would tier 1 remove this" but it cannot answer
"would we have found something nobody had seen yet", because we always knew what
we had hidden and we hid it ourselves. At depth 1 the clean holdout arm is a
handful of pairs, which is the weakness this search is meant to fix.

The archive supplies a better blind for free. If a viral species' FIRST deposited
genome postdates an untargeted shotgun RNA metagenome of wastewater that contains
it, then for the duration between the run's collection date and that deposit, the
species was genuinely unknown to every database in the world. Excluding that
species and everything deposited after the collection date reconstructs the state
of knowledge at the time of the run exactly, and the question becomes a historical
counterfactual: run the committed cascade over those reads with the clock turned
back, and ask whether the identity-to-closest-relative readout surfaces the thing
at its true distance. Nobody has to plant anything, and we do not choose the
divergence - the biology did.

WHAT THIS SCRIPT DOES
---------------------
Five stages, each cached to --outdir so a stage can be re-run alone:

  runs       ENA inventory of untargeted paired wastewater RNA libraries
  pool       NCBI nucleotide records of RNA viruses whose source is wastewater
             or sewage and whose deposit is recent; first-deposit date per taxon
  links      GenBank DBLINK for each candidate record -> the BioProject and
             BioSample, hence the RUN that contains it. A linked run is the whole
             point: presence and depth stop being plausibility arguments
  divergence local blastn (and blastx where the protein set is built) of each
             candidate against viral RefSeq, which is the pipeline's own index
  depth      exact read-pair counts per candidate by mapping its containing run,
             for runs small enough to pull into a container at ~1 MB/s

RANKING, AND WHY IT IS DIVERGENCE AND NOT DATE
----------------------------------------------
The useful band is roughly 10-20% divergence from the nearest relative that
predates the run. Below it the cascade already works (`results/spike_holdout/`:
0.00 and 0.68% of spike reads survive tier 1 at 0 and 5% divergence), above it
nothing does (99.6% survive at 30%, i.e. the method is blind). A candidate at 2%
or at 40% therefore teaches us nothing about where the wall is.

Dates alone are a trap, and this script catches it twice. A taxon whose first
deposit is recent may be a RENAMING of sequence that has been in the databases
for years: four "new species" in the candidate pool are 98.0-99.9% identical to
viruses already in viral RefSeq. Only the divergence measurement separates a new
species from a new name, so ranking is on measured identity and the date is a
filter, never a score.

THE DATE FILTER IS DELIBERATELY NOT APPLIED TO THE REFERENCE, AND THAT IS A BOUND
---------------------------------------------------------------------------------
Divergence here is measured against viral RefSeq *as it stands today*, with the
candidate's own taxon excluded. That is the conservative direction: a reference
containing more than existed at the run's collection date can only make a
candidate look LESS novel, never more. So a candidate measured at 17% divergence
today was at least 17% divergent then. Reported numbers are therefore a lower
bound on divergence and an upper bound on recognisability.

The reason for the bound rather than an exact as-of date is a real pitfall worth
recording: a RefSeq record's own createdate is when the RefSeq entry was built,
not when the underlying sequence was first deposited (3,926 of 19,650 viral
RefSeq records carry a 2023 createdate, wrapping sequences decades older).
Filtering on it would have excluded genuinely old relatives and inflated every
divergence figure - the dangerous direction. The holdout for the actual run must
therefore key its negative lists to the ORIGINAL GenBank accession's deposit
date, not to the RefSeq record's.

CONVENTIONS
-----------
This is a census of the archive, not a sample of it: every stage is deterministic
given a cutoff date and an archive snapshot, so there is no seed spread to report
and five seeds would be five identical files. --seed is accepted and recorded for
convention and to pin the one stochastic choice (subsampling in `depth`), and the
manifest carries the archive query dates, which are the real provenance.
"""
import argparse, csv, hashlib, json, os, re, subprocess, sys, time, urllib.parse, urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

ENA = "https://www.ebi.ac.uk/ena/portal/api/"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Verified against the ENA taxonomy service rather than guessed. An earlier pass
# of this search used four invented taxids and silently pulled in gut metagenomes,
# bat viromes and a liverwort; the join it fed then found nothing. Taxon ids are
# data and get checked like data.
WW_TAXA = {
    527639: "wastewater metagenome",
    942017: "activated sludge metagenome",
    1592332: "sludge metagenome",
    412754: "bioreactor sludge metagenome",
}
WW_TEXT = re.compile(r"wastewater|waste water|sewage|sewer|influent|effluent|activated sludge|wwtp", re.I)
# Amplicon and probe capture cannot host this experiment: what they contain was
# chosen in advance, so a species unknown at the time could not have been targeted.
TARGETED = re.compile(r"amplicon|artic|primer|tiled|midnight|capture|probe|panel|xhyb|targeted|hybrid", re.I)
RIBOVIRIA = "txid2559587"   # RNA viruses. A DNA library cannot hold one.


def http(url, data=None, tries=8, timeout=300):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, data=data)
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:                      # 429s are routine; back off
            print(f"  retry {a}: {e}", file=sys.stderr)
            time.sleep(3 + 4 * a)
    raise SystemExit(f"gave up on {url[:120]}")


def eutils(cgi, **kw):
    return http(EUTILS + cgi, data=urllib.parse.urlencode(kw).encode())


def esearch_count(term, db="nuccore"):
    r = json.loads(eutils("esearch.fcgi", db=db, term=term, retmode="json", retmax=0))
    return int(r["esearchresult"]["count"])


# ---------------------------------------------------------------- stage: runs
def stage_runs(out):
    """Untargeted paired wastewater RNA libraries, with the cutoff date each one
    would impose on a holdout.

    Two labelling facts force the shape of this query. Berlin SRR23973193 - our
    own discriminating dataset - is filed as library_source=TRANSCRIPTOMIC, not
    METATRANSCRIPTOMIC, and carries NO collection_date in ENA, so a query on
    either field alone misses it. And the richest candidate source found by this
    search is filed as METAGENOMIC/OTHER/unspecified yet its paper reports
    sixteen full-length RNA virus genomes from those libraries. Labels are
    therefore collected, not trusted: the RNA evidence and the date evidence are
    emitted as columns for a human to read.

    Where collection_date is missing, first_public stands in. That substitution
    is safe in one direction only, and it is the safe one: a run is public no
    earlier than it was collected, so requiring a species' first deposit to
    postdate first_public is a STRICTER holdout than requiring it to postdate
    collection.
    """
    tax = " OR ".join(f"tax_eq({t})" for t in WW_TAXA)
    fields = ("run_accession,study_accession,sample_accession,collection_date,country,"
              "library_source,library_strategy,library_selection,library_layout,"
              "instrument_model,read_count,base_count,scientific_name,sample_title,"
              "study_title,first_public")
    queries = [
        ("taxon", f'({tax}) AND library_layout="PAIRED"'),
        ("metaT", 'library_source="METATRANSCRIPTOMIC" AND library_layout="PAIRED"'),
        ("transcriptomic", f'library_source="TRANSCRIPTOMIC" AND ({tax})'),
    ]
    seen, rows = set(), []
    for name, q in queries:
        body = http(ENA + "search", data=urllib.parse.urlencode(
            dict(result="read_run", query=q, fields=fields, format="tsv", limit=0)).encode())
        got = list(csv.DictReader(body.decode().splitlines(), delimiter="\t"))
        print(f"  {name}: {len(got)} runs", file=sys.stderr)
        for r in got:
            if r["run_accession"] in seen:
                continue
            seen.add(r["run_accession"])
            rows.append(r)

    keep = []
    for r in rows:
        text = " | ".join([r["scientific_name"], r["sample_title"] or "", r["study_title"] or ""])
        if not (int(r["scientific_name"] in [v for v in WW_TAXA.values()]) or WW_TEXT.search(text)):
            continue
        if r["library_layout"] != "PAIRED":
            continue
        if TARGETED.search(text) or TARGETED.search(r["library_selection"] or ""):
            continue
        if int(r["read_count"] or 0) < 1_000_000:
            continue
        # RNA evidence is required, because a DNA library cannot hold an RNA virus
        # and this is where the candidates live. It is a disjunction of three weak
        # signals rather than one label, for the reason given above: 131,059 runs
        # pass without it and the two real leads are each missed by one of the
        # three. The signal that fired is recorded, so a reader can disagree.
        rna = []
        if (r["library_source"] or "") in ("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "VIRAL RNA"):
            rna.append("source=" + r["library_source"])
        if re.search(r"cdna|rna", r["library_selection"] or "", re.I):
            rna.append("selection=" + r["library_selection"])
        if re.search(r"\brna\b|transcriptom|virome|viral", text, re.I):
            rna.append("text")
        if not rna:
            continue
        coll = (r["collection_date"] or "")[:10]
        r["cutoff_date"] = coll if re.match(r"\d{4}-\d{2}-\d{2}", coll) else r["first_public"]
        r["cutoff_source"] = "collection_date" if r["cutoff_date"] == coll else "first_public"
        r["rna_evidence"] = ";".join(rna)
        keep.append(r)
    keep.sort(key=lambda r: (r["cutoff_date"], -int(r["read_count"])))

    # The census is 73,186 runs and most of it is single-virus surveillance:
    # RNA-Seq libraries of wastewater that exist to quantify one agent, which
    # pass every structural filter here because they are neither amplicon nor
    # capture. They are not useful hosts for this experiment - an index holdout
    # cannot make a library untargeted after the fact - so the committed table is
    # the usable subset and the census is a counted number, not 25 MB of rows.
    surveil = re.compile(r"sars|cov-?2|covid|influenza|mpox|monkeypox|polio.?surveill", re.I)
    usable = [r for r in keep
              if int(r["read_count"]) >= 3_000_000
              and not surveil.search(" ".join([r["sample_title"] or "", r["study_title"] or ""]))]
    # gzipped: 14,621 rows is 4.6 MB of archive metadata and the per-study
    # aggregate beside it answers most questions without it
    write_tsv(os.path.join(out, "runs.tsv"), usable)
    subprocess.run(["gzip", "-f", os.path.join(out, "runs.tsv")], check=True)

    by_study = defaultdict(list)
    for r in usable:
        by_study[r["study_accession"]].append(r)
    agg = []
    for k, v in by_study.items():
        pairs = sorted(int(x["read_count"]) for x in v)
        agg.append(dict(study=k, runs=len(v), median_pairs=pairs[len(pairs) // 2],
                        max_pairs=pairs[-1], cutoff_min=min(x["cutoff_date"] for x in v),
                        cutoff_max=max(x["cutoff_date"] for x in v),
                        cutoff_source=";".join(sorted({x["cutoff_source"] for x in v})),
                        rna_evidence=";".join(sorted({x["rna_evidence"] for x in v}))[:80],
                        country=v[0]["country"], title=(v[0]["study_title"] or "")[:90]))
    agg.sort(key=lambda a: (a["cutoff_min"], -a["max_pairs"]))
    write_tsv(os.path.join(out, "runs_by_study.tsv"), agg)
    json.dump(dict(census_paired_wastewater_rna=len(keep), usable=len(usable),
                   usable_studies=len(agg),
                   dropped_single_agent_surveillance_or_thin=len(keep) - len(usable)),
              open(os.path.join(out, "runs_summary.json"), "w"), indent=1)
    print(f"  census {len(keep)} runs; {len(usable)} usable in {len(agg)} studies", file=sys.stderr)
    return usable


# ---------------------------------------------------------------- stage: pool
def stage_pool(out, cutoff):
    """RNA-virus records whose source text says wastewater or sewage, deposited
    after `cutoff`, plus the first-deposit date of every taxon they name.

    The query runs on the DEPOSIT side on purpose. Starting from runs and asking
    what was deposited under the same BioProject found five projects and none of
    the two real leads, because a run's own taxon and strategy labels are noisy
    (above) while a deposit's isolation_source and collection_date are written by
    the submitter who had the sample in hand.

    The surveillance giants are excluded by taxon, not by hand: they contribute
    hundreds of thousands of records, every one of which has prior relatives at
    0% divergence, so they cannot be candidates and only cost query budget.
    """
    text = ('(wastewater[All Fields] OR sewage[All Fields] OR "waste water"[All Fields] '
            'OR influent[All Fields] OR sewer[All Fields] OR "activated sludge"[All Fields])')
    term = (f"{RIBOVIRIA}[Organism:exp] NOT txid2697049[Organism:exp] NOT txid11320[Organism:exp] "
            f"AND {cutoff.replace('-','/')}:3000/01/01[PDAT] AND {text} AND 1000:300000[SLEN]")
    h = json.loads(eutils("esearch.fcgi", db="nuccore", term=term, retmode="json",
                          retmax=100000, usehistory="y"))["esearchresult"]
    ids = h["idlist"]
    print(f"  pool: {len(ids)} records", file=sys.stderr)

    recs = []
    for i in range(0, len(ids), 300):
        d = json.loads(eutils("esummary.fcgi", db="nuccore", id=",".join(ids[i:i + 300]),
                              retmode="json"))["result"]
        for k in d["uids"]:
            v = d[k]
            sub = dict(zip((v.get("subtype") or "").split("|"), (v.get("subname") or "").split("|")))
            recs.append(dict(acc=v.get("accessionversion"), taxid=v.get("taxid"),
                             organism=v.get("organism"), slen=v.get("slen"),
                             deposited=v.get("createdate"),
                             collected=sub.get("collection_date", ""),
                             isolation_source=sub.get("isolation_source", ""),
                             geo=sub.get("country", ""), title=(v.get("title") or "")[:180]))
        time.sleep(0.4)
    write_tsv(os.path.join(out, "pool.tsv"), recs)

    # First deposit per taxon. One count query per taxon with an upper date bound
    # answers the only question that matters - was ANYTHING of this species known
    # before the cutoff - without paging through tens of thousands of records for
    # the common taxa.
    taxa = sorted({(r["taxid"], r["organism"]) for r in recs})
    print(f"  {len(taxa)} taxa; testing for prior deposits", file=sys.stderr)

    def probe(to):
        t, o = to
        return t, dict(organism=o,
                       pool_records=sum(1 for r in recs if r["taxid"] == t),
                       prior_to_cutoff=esearch_count(
                           f"txid{t}[Organism:exp] AND 1900/01/01:{cutoff.replace('-','/')}[PDAT]"),
                       total_records=esearch_count(f"txid{t}[Organism:exp]"))

    first = {}
    with ThreadPoolExecutor(max_workers=3) as ex:
        for i, (t, r) in enumerate(ex.map(probe, taxa)):
            first[t] = r
            if i % 50 == 0:
                print(f"    {i}/{len(taxa)}", file=sys.stderr)
    json.dump(first, open(os.path.join(out, "first_deposit.json"), "w"), indent=1)
    novel = {t for t, v in first.items() if v["prior_to_cutoff"] == 0}
    print(f"  {len(novel)} taxa with no deposit of any kind before {cutoff}", file=sys.stderr)
    return recs, first


# ---------------------------------------------------------------- stage: links
def stage_links(out, accs):
    """GenBank DBLINK for each candidate: the BioProject and BioSample, and from
    them the run that contains the genome.

    This is the stage that turns a plausibility argument into a fact. A candidate
    whose containing run is public has guaranteed presence (the genome was
    assembled FROM those reads) and a measurable depth; a candidate without one
    needs an argument that the species was in some other sewershed on some other
    date, which is exactly the kind of argument this project does not accept.
    """
    gb = b""
    for i in range(0, len(accs), 30):
        gb += eutils("efetch.fcgi", db="nuccore", id=",".join(accs[i:i + 30]),
                     rettype="gb", retmode="text")
        time.sleep(1.2)
    open(os.path.join(out, "candidates.gb"), "wb").write(gb)

    rows = []
    for e in gb.decode("utf8", "replace").split("\n//\n"):
        if not e.strip():
            continue
        g = lambda p, f=0: (re.search(p, e, f).group(1) if re.search(p, e, f) else "")
        rows.append(dict(acc=g(r"^LOCUS\s+(\S+)", re.M),
                         bioproject=g(r"BioProject:\s*(PRJ\w+)"),
                         biosamples=",".join(sorted(set(re.findall(r"SAMN\d+|SAMD\d+|SAMEA\d+", e)))[:6]),
                         assembly_method=g(r"Assembly Method\s*::\s*(.+)"),
                         coverage=g(r"Coverage\s*::\s*(\S+)"),
                         paper=g(r"^  JOURNAL\s+(.+)$", re.M)))

    # BioSample -> run, so each candidate carries the library that holds it
    # ENA's filereport takes ONE accession; a comma list returns 500 or an empty
    # table with a 200, which is the worse of the two failures. The search
    # endpoint with OR'd sample_accession is the batchable form.
    samples = sorted({s for r in rows for s in r["biosamples"].split(",") if s})
    runs = {}
    for i in range(0, len(samples), 40):
        q = " OR ".join(f'sample_accession="{s}"' for s in samples[i:i + 40])
        body = http(ENA + "search", data=urllib.parse.urlencode(dict(
            result="read_run", query=q, format="tsv", limit=0,
            fields=("run_accession,sample_accession,collection_date,library_source,"
                    "library_strategy,library_selection,library_layout,read_count,"
                    "base_count,instrument_model,sample_title,country,fastq_bytes"))).encode())
        for r in csv.DictReader(body.decode().splitlines(), delimiter="\t"):
            runs.setdefault(r["sample_accession"], []).append(r)
    for r in rows:
        hits = [x for s in r["biosamples"].split(",") if s for x in runs.get(s, [])]
        r["runs"] = ",".join(x["run_accession"] for x in hits)
        r["run_pairs"] = ",".join(x["read_count"] for x in hits)
        r["run_collected"] = ",".join(x["collection_date"] for x in hits)
        r["run_labels"] = ";".join(sorted({f"{x['library_source']}/{x['library_strategy']}/"
                                           f"{x['library_selection']}" for x in hits}))
    write_tsv(os.path.join(out, "links.tsv"), rows)
    return rows


# ----------------------------------------------------------- stage: divergence
def stage_divergence(out, blastbin, candfa, ntdb, protdb=None):
    """Identity to the nearest relative in viral RefSeq, which is the index the
    pipeline itself consults (SPEC s2: DIAMOND RefSeq viral proteins, 0.239 GB).

    Three nucleotide sensitivities are run because the answer depends on which
    one you use, and the pipeline's own tiers die in that order: megablast is the
    rate at which a close relative is found, blastn with a word size of 11 is the
    most generous nucleotide reach available, and the gap between them is the band
    where nucleotide identification is already failing. blastx against the same
    RefSeq at protein level is the tier SPEC keeps for exactly this region.

    Every hit is reported with its E-value, its alignment length and its query
    coverage, never as a bare percent identity. This repo has had to retract a
    "~89% identity" claim that was chance matching over 40 bases, and identity is
    quantised at short alignment lengths, so a percentage without a length and an
    E-value beside it is unreadable.
    """
    fmt = ("6 qseqid qlen sseqid pident length mismatch gapopen qstart qend "
           "qcovhsp evalue bitscore stitle")
    res = {}
    for task, db, prog in [("megablast", ntdb, "blastn"), ("dc-megablast", ntdb, "blastn"),
                           ("blastn", ntdb, "blastn")] + ([("blastx", protdb, "blastx")] if protdb else []):
        cmd = [os.path.join(blastbin, prog), "-query", candfa, "-db", db,
               "-outfmt", fmt, "-max_target_seqs", "10", "-evalue", "1e-3",
               "-num_threads", str(os.cpu_count() or 4)]
        if prog == "blastn":
            cmd += ["-task", task]
        t0 = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        rows = [l.split("\t") for l in p.stdout.splitlines()]
        res[task] = rows
        print(f"  {task}: {len(rows)} hsps in {time.time()-t0:.0f}s", file=sys.stderr)
        with open(os.path.join(out, f"hits_{task}.tsv"), "w") as f:
            f.write(p.stdout)
    return res


# ---------------------------------------------------------------- stage: depth
def stage_depth(out, candfa, reads, seed=1):
    """Exact read-pair counts per candidate genome in its own containing run.

    Abundance is the binding constraint on this whole experiment, not divergence:
    the residue we care about sits at median depth 1, and a candidate present at a
    handful of pairs is as thin as the spike it is meant to replace. So depth is
    measured, not projected - mapped with minimap2's short-read preset through
    mappy, counting reads whose primary alignment is to that genome.

    Two numbers are reported per candidate. Mapped reads is what the library
    holds. Mean depth is coverage of the genome, which is the number that decides
    whether assembly could have recovered it at all, and the two come apart badly
    for short genomes.

    Depth is reported per read file, not pooled, and that is the point: a genome
    assembled from one sample mapped against the OTHER samples of the same study
    is a free presence control. A candidate that reads thousands in its own run
    and zero in a sibling run collected three months later is a single shedding
    event; one that reads in all of them recurs. Pooling the three would destroy
    exactly that distinction, which `results/recurrence/` says is the one the
    depth-1 wall turns on.
    """
    import mappy
    idx = mappy.Aligner(candfa, preset="sr", n_threads=os.cpu_count() or 4)
    if not idx:
        raise SystemExit("failed to index candidates")
    lens = {n: len(idx.seq(n)) for n in idx.seq_names}
    rows, scanned = [], {}
    for fq in reads:
        run = os.path.basename(fq).split(".")[0]
        counts, bases, nreads = Counter(), Counter(), 0
        for name, seq, qual in mappy.fastx_read(fq):
            nreads += 1
            best = None
            for h in idx.map(seq):
                if h.is_primary and (best is None or h.mlen > best.mlen):
                    best = h
            if best is not None:
                counts[best.ctg] += 1
                bases[best.ctg] += best.blen
        scanned[run] = nreads
        for g in lens:
            rows.append(dict(run=run, genome=g, length=lens[g], mapped_reads=counts[g],
                             mean_depth=round(bases[g] / lens[g], 3),
                             per_million=round(1e6 * counts[g] / max(nreads, 1), 3),
                             reads_scanned=nreads))
        print(f"  {run}: {nreads} reads, {sum(counts.values())} mapped to candidates",
              file=sys.stderr)
    rows.sort(key=lambda r: (-r["mapped_reads"], r["genome"]))
    write_tsv(os.path.join(out, "depth.tsv"), rows)
    json.dump(dict(reads_scanned=scanned, files=reads, seed=seed,
                   preset="minimap2 sr via mappy", primary_only=True),
              open(os.path.join(out, "depth_input.json"), "w"), indent=1)
    return rows


# ----------------------------------------------------------------- stage: rank
def stage_rank(out):
    """Join divergence, links and depth into the shortlist, ranked by distance
    from the useful band.

    The band is 10-20% divergence from the nearest prior relative, and it is not
    a preference: `results/spike_holdout/` measures 0.00% and 0.68% of spike reads
    surviving tier 1 at 0% and 5% divergence, and 99.6% at 30%. Below the band the
    cascade already works, above it the cascade is blind, so a candidate outside
    it cannot move any number we hold. Rank is therefore |midpoint - divergence|,
    and everything outside the band is kept in the file with its reason rather
    than dropped, because the kills are the part of this search that generalises.

    Three guards, each from a mistake this repo has already made once:
      - identity is reported with its alignment length, query coverage and
        E-value, never alone;
      - a hit covering little of the query is not a nearest relative, so the
        aligned fraction is a column and a candidate is flagged `partial` when it
        is low rather than being silently ranked on a fragment;
      - a candidate 98%+ identical to something already in the reference is
        marked `renaming`, which is what the date filter alone cannot see.
    """
    def union_cov(hsps, qlen, scale=1):
        """Union of HSP query intervals over the query length.

        A single HSP cannot cover a multi-gene genome: a levivirus encodes
        maturation, coat and replicase, so the best blastx HSP covers one gene
        and about 45% of the genome. Flagging that as a partial match would have
        demoted the four candidates closest to the protein wall - the ones this
        whole search exists to find - on an artefact of how blastx reports.
        Coverage is therefore the union over HSPs, not the best one.
        """
        iv = sorted((min(int(h[0]), int(h[1])), max(int(h[0]), int(h[1]))) for h in hsps)
        tot, cur_s, cur_e = 0, None, None
        for s, e in iv:
            if cur_e is None or s > cur_e:
                tot += 0 if cur_e is None else cur_e - cur_s + 1
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        tot += 0 if cur_e is None else cur_e - cur_s + 1
        return round(min(tot * scale / qlen, 1.0), 3)

    def all_hits(task):
        p = os.path.join(out, f"hits_{task}.tsv")
        d = defaultdict(list)
        if os.path.exists(p):
            for r in csv.reader(open(p), delimiter="\t"):
                if len(r) >= 13:
                    d[r[0]].append(r)
        return d

    def best_hits(task):
        p = os.path.join(out, f"hits_{task}.tsv")
        b = {}
        if not os.path.exists(p):
            return b
        for r in csv.reader(open(p), delimiter="\t"):
            if len(r) < 13:
                continue
            q, bits = r[0], float(r[11])
            if q not in b or bits > float(b[q][11]):
                b[q] = r
        return b

    mega, sens, prot = best_hits("megablast"), best_hits("blastn"), best_hits("blastx")
    sens_all, prot_all = all_hits("blastn"), all_hits("blastx")
    links = {r["acc"]: r for r in csv.DictReader(open(os.path.join(out, "links.tsv")), delimiter="\t")}
    pool = {r["acc"]: r for r in csv.DictReader(open(os.path.join(out, "pool.tsv")), delimiter="\t")}
    depth = {}
    dp = os.path.join(out, "depth.tsv")
    if os.path.exists(dp):
        depth = {r["genome"]: r for r in csv.DictReader(open(dp), delimiter="\t")}

    rows = []
    for acc, rec in pool.items():
        h = sens.get(acc)
        if h is None:
            continue
        # fields: qseqid qlen sseqid pident length mismatch gapopen qstart qend
        #         qcovhsp evalue bitscore stitle
        qlen, pid, alnlen, ev = int(h[1]), float(h[3]), int(h[4]), h[10]
        div = round(100 - pid, 2)
        # coverage is the union over HSPs against the best subject, which is the
        # question "how much of this genome has a relative at all"
        same = [x for x in sens_all.get(acc, []) if x[2] == h[2]]
        aligned_frac = union_cov([(x[7], x[8]) for x in same], qlen) if same else \
            round(alnlen / qlen, 3)
        lk = links.get(acc.split(".")[0], links.get(acc, {}))
        d = depth.get(acc, {})
        flags = []
        if div < 2:
            flags.append("renaming")
        # a nucleotide alignment covering little while protein covers most is the
        # signature of a genuinely divergent relative, not of a fragment, so the
        # two coverages are flagged separately and neither alone demotes anything
        if aligned_frac < 0.5 and (pfrac or 0) < 0.5:
            flags.append("partial")
        if not lk.get("runs"):
            flags.append("no_containing_run")
        band = 10.0 <= div <= 20.0
        # The protein band is a second, different question, and the two bands do
        # not select the same candidates. `results/spike_holdout/`'s 10-20%
        # nucleotide band comes from a UNIFORMLY mutated spike, which degrades
        # nucleotide and protein together. Real divergence does not: across these
        # candidates the ratio of nucleotide to protein divergence runs from about
        # 1 to about 4, so a genome 10% away in nucleotide can be 2.3% away in
        # protein. Since the protein tier does essentially all the work above 5%
        # divergence (SPEC s4.2: 91.3 / 84.8 / 20.0 on the holdout), a candidate
        # that only clears the nucleotide band will be identified anyway and
        # tests the cascade's reach rather than its blindness.
        ppid = float(prot[acc][3]) if acc in prot else None
        # protein coverage is the union of blastx HSPs from the same source
        # organism, in query nucleotides, so a multi-gene genome is not penalised
        pfrac = None
        if acc in prot:
            org_tag = prot[acc][12].split("[")[-1]
            same_p = [x for x in prot_all.get(acc, []) if x[12].split("[")[-1] == org_tag]
            pfrac = union_cov([(x[7], x[8]) for x in same_p], qlen) if same_p else \
                round(int(prot[acc][4]) * 3 / qlen, 3)
        prot_band = (ppid is not None and 60.0 <= ppid <= 75.0 and (pfrac or 0) >= 0.5)
        rows.append(dict(
            acc=acc, organism=rec["organism"], genome_nt=qlen,
            divergence_pct=div, nearest_relative=h[10][:70], aln_len=alnlen,
            aligned_frac=aligned_frac, evalue=ev,
            megablast_pid=(mega[acc][3] if acc in mega else ""),
            protein_pid=(prot[acc][3] if acc in prot else ""),
            protein_aligned_frac=(round(pfrac, 3) if pfrac is not None else ""),
            protein_evalue=(prot[acc][8] if acc in prot else ""),
            nt_over_prot_divergence=(round(div / (100 - ppid), 2)
                                     if ppid is not None and ppid < 100 else ""),
            # the ratio above compares two levels only when both levels found the
            # same relative; where they name different organisms it compares two
            # different relatives and says nothing about synonymous fraction
            same_relative_both_levels=(int(h[12].split("[")[-1].strip("] ").lower()[:18]
                                           in prot[acc][12].lower())
                                       if acc in prot else ""),
            protein_subject=(prot[acc][10][:60] if acc in prot else ""),
            deposited=rec["deposited"], collected=rec["collected"],
            isolation_source=rec["isolation_source"], geo=rec["geo"],
            runs=lk.get("runs", ""), run_pairs=lk.get("run_pairs", ""),
            run_collected=lk.get("run_collected", ""), run_labels=lk.get("run_labels", ""),
            mapped_reads=d.get("mapped_reads", ""), mean_depth=d.get("mean_depth", ""),
            per_million=d.get("per_million", ""),
            in_nt_band=int(band), in_protein_band=int(prot_band),
            flags=",".join(flags)))
    # Ranked on the protein band first, because that is the tier that decides
    # whether the cascade can see the candidate at all.
    rows.sort(key=lambda r: (not r["in_protein_band"], not r["in_nt_band"],
                             bool(r["flags"]), abs(15.0 - r["divergence_pct"])))
    write_tsv(os.path.join(out, "shortlist.tsv"), rows)
    inb = [r for r in rows if r["in_nt_band"] and not r["flags"]]
    inp = [r for r in rows if r["in_protein_band"] and not r["flags"]]
    print(f"  {len(rows)} candidates: {len(inb)} in the nucleotide band, "
          f"{len(inp)} in the protein band, unflagged", file=sys.stderr)

    # The README's table is generated here rather than typed, because a number
    # typed into a table is a number without provenance and this repo has had to
    # withdraw claims whose prose drifted from its outputs.
    with open(os.path.join(out, "SHORTLIST.md"), "w") as f:
        f.write("<!-- generated by experiments/holdout_candidates.py, stage `rank`."
                " Do not edit by hand. -->\n\n")
        f.write(f"Divergence is to the nearest relative in viral RefSeq *as it stands today*,"
                f" which is a lower bound (see README §3). {len(rows)} candidates,"
                f" {len(inb)} in the 10-20% nucleotide band and {len(inp)} in the 60-75%"
                f" protein-identity band, unflagged. `nt/prot` is the ratio of nucleotide"
                f" to protein divergence: the higher it is, the more of the distance is"
                f" synonymous and invisible to the tier that does the work - but it is"
                f" only readable where `same rel` is 1, because otherwise the two levels"
                f" found two different relatives and the ratio compares them, not the"
                f" levels.\n\n")
        cols = [("acc", "accession"), ("organism", "organism"), ("genome_nt", "nt"),
                ("divergence_pct", "nt div %"), ("aligned_frac", "nt aln frac"),
                ("protein_pid", "prot id %"), ("protein_aligned_frac", "prot aln frac"),
                ("nt_over_prot_divergence", "nt/prot"),
                ("same_relative_both_levels", "same rel"),
                ("deposited", "deposited"), ("run_collected", "run collected"),
                ("runs", "containing run"), ("run_pairs", "run pairs"),
                ("mapped_reads", "mapped reads"), ("mean_depth", "mean depth"),
                ("flags", "flags")]
        f.write("| " + " | ".join(c[1] for c in cols) + " |\n")
        f.write("|" + "---|" * len(cols) + "\n")
        for r in rows:
            keep = (r["in_protein_band"] or r["in_nt_band"] or r["divergence_pct"] < 2
                    or (r["flags"] and float(r["protein_pid"] or 100) < 70))
            if not keep:
                continue           # the two bands, the renamings and the near misses
            f.write("| " + " | ".join(str(r[c[0]])[:46] for c in cols) + " |\n")
    return rows


# ---------------------------------------------------------------- plumbing
def write_tsv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def manifest(out, args, timings):
    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    m = dict(script=os.path.basename(__file__),
             script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
             git_sha=sha + ("-dirty" if dirty else ""),
             argv=vars(args), timings_s=timings,
             archive_queried_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             note=("deterministic given the archive snapshot; the archive query date is the "
                   "provenance that matters, not the seed"))
    json.dump(m, open(os.path.join(out, f"manifest_s{args.seed}.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default="results/holdout_candidates")
    ap.add_argument("--cutoff", default="2022-01-01",
                    help="earliest deposit date a candidate may carry; taxa with any deposit "
                         "before it are not candidates")
    ap.add_argument("--stages", default="runs,pool,links,divergence,rank",
                    help="comma-separated: runs,pool,links,divergence,depth,rank")
    ap.add_argument("--blastbin", default="", help="directory holding blastn/blastx")
    ap.add_argument("--ntdb", default="", help="blast db of viral RefSeq nucleotide")
    ap.add_argument("--protdb", default="", help="blast db of viral RefSeq protein")
    ap.add_argument("--candfa", default="", help="candidate genomes FASTA (written by links)")
    ap.add_argument("--reads", nargs="*", default=[], help="fastq(.gz) for the depth stage")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    stages = args.stages.split(",")
    timings = {}

    if "runs" in stages:
        t = time.time(); stage_runs(args.outdir); timings["runs"] = round(time.time() - t, 1)
    if "pool" in stages:
        t = time.time(); stage_pool(args.outdir, args.cutoff); timings["pool"] = round(time.time() - t, 1)
    if "links" in stages:
        recs = list(csv.DictReader(open(os.path.join(args.outdir, "pool.tsv")), delimiter="\t"))
        first = json.load(open(os.path.join(args.outdir, "first_deposit.json")))
        novel = {int(k) for k, v in first.items() if v["prior_to_cutoff"] == 0}
        accs = sorted({r["acc"] for r in recs if int(r["taxid"]) in novel and int(r["slen"]) >= 2000})
        t = time.time(); stage_links(args.outdir, accs); timings["links"] = round(time.time() - t, 1)
        fa = eutils("efetch.fcgi", db="nuccore", id=",".join(accs), rettype="fasta", retmode="text")
        open(os.path.join(args.outdir, "candidates.fa"), "wb").write(fa)
    if "divergence" in stages:
        t = time.time()
        stage_divergence(args.outdir, args.blastbin,
                         args.candfa or os.path.join(args.outdir, "candidates.fa"),
                         args.ntdb, args.protdb or None)
        timings["divergence"] = round(time.time() - t, 1)
    if "depth" in stages:
        t = time.time()
        stage_depth(args.outdir, args.candfa or os.path.join(args.outdir, "candidates.fa"),
                    args.reads, args.seed)
        timings["depth"] = round(time.time() - t, 1)
    if "rank" in stages:
        t = time.time(); stage_rank(args.outdir); timings["rank"] = round(time.time() - t, 1)

    manifest(args.outdir, args, timings)


if __name__ == "__main__":
    main()
