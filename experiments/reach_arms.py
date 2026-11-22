"""Five ways of finding a distant relative, scored on one panel in one pass.

WHAT IS BEING MEASURED. For each query contig: the most distant relative each method
recovers, with the E-value, identity, aligned length and query coverage of that alignment;
whether the method found anything the baseline did not; and what it cost in wall clock. A
method that finds nothing new is the result, not a failure to report.

WHY THE SUBJECT UNIVERSE IS HELD FIXED WHERE IT CAN BE. Reach is a property of an
instrument *and* of the set it searches, and this project's own measured lesson is that the
set dominates: 1.4 MB of the sample's own contigs removed more of the library than 272 MB of
SILVA. So the arms that can share a subject set do:

  * Swiss-Prot, as sequence   -> DIAMOND (arm 1), BLASTP and exhaustive Smith-Waterman
                                 (arm 2), jackhmmer (arm 3's same-universe half)
  * Swiss-Prot, as structure  -> Foldseek against AlphaFold/Swiss-Prot (arm 4)

AlphaFold's Swiss-Prot set is Swiss-Prot folded, so arm 4 against it is the same subjects in
a different alphabet. **That is the one genuinely controlled reach comparison available
here**, and it is the comparison the headline table is built on. Pfam (arm 3's profile
library), PDB (arm 4's experimental half) and RefSeq viral genomes (arm 5) are different
universes and are reported as their own rows, never subtracted from the Swiss-Prot ones.

WHY THE SEVEN ARMS ARE NOT ALL "SEARCH". Arms 1 and 2 differ only in whether seeding is
heuristic, so the difference between them *is* the cost of the shortcut, measured rather
than assumed. Arm 3 replaces a sequence with a profile. Arm 4 replaces a sequence with a
structural alphabet. Arm 5 scores in codons rather than residues. Each one changes exactly
one thing about the comparison, which is why a difference can be attributed at all.

EVERY ARM IS RUN ON EVERY PILE, including the dinucleotide-shuffled null and the
near-reference negative arm, in the same pass with the same parameters. Reading the control
first is the rule here (CLAUDE.md): a method that reaches further on the near-reference arm
too is reading its own threshold. The null fixes the floor; the near-reference arm fixes the
direction.

NO THRESHOLDS ARE APPLIED. One row per alignment, with the numbers that let a floor be
drawn afterwards. `results/README.md` records six single-axis thresholds that died; the
deliverable this project keeps arriving at is a distance with a coverage beside it, so that
is what the rows carry.

ONE ARM'S FAILURE DOES NOT COST THE RUN. There is no SSH and no second chance, so each arm
is wrapped: an arm that dies records its traceback in its own timing file and the rest
continue. The script fails only if every arm failed, because that means the setup is wrong
rather than one tool.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import time
import traceback
from pathlib import Path

#: BLAST and DIAMOND tabular fields, in this order, for both programs so one parser reads
#: both. `qcovhsp` is here because it is the field that makes an identity readable: a 40%
#: alignment at 95% identity and a full-length one at 95% are the same `pident` and
#: completely different claims, and this project has had to retract a number for exactly
#: that reason.
TAB_FIELDS = ("qseqid", "sseqid", "pident", "length", "qlen", "slen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")

#: The genetic code the SUBJECT side is read with. Table 11 is the bacterial and phage
#: table, and the test arm is ssRNA bacteriophage, so 11 rather than 1.
#:
#: The QUERY side does not use this, and that is a correctness point rather than a detail.
#: Prodigal in metagenome mode selects a translation table per sequence from its own
#: training set, and on this panel it chose table 4 for some records - where TGA codes
#: tryptophan rather than stopping. Translating those genes with table 11 anyway cuts each
#: protein at its first internal TGA, which would have shortened real ORFs and handed every
#: downstream arm a truncated query. So each ORF is translated with the table its genes were
#: actually called with, and that table is recorded on the row.
SUBJECT_GENETIC_CODE = 11

#: Minimum ORF length in residues. 50 matches `results/divergent_known/` and
#: `results/divergent_contigs/`, so the ORF pile here is the same pile those reports
#: describe. Shorter ORFs are where identity goes quantised and optimistic.
MIN_ORF_AA = 50


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    print(f"+ {' '.join(str(c) for c in cmd)}", flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    """(name, description, sequence) in file order."""
    out, name, desc, chunks = [], None, "", []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, desc, "".join(chunks)))
                head = line[1:].rstrip("\n")
                parts = head.split(None, 1)
                name, desc, chunks = parts[0], (parts[1] if len(parts) > 1 else ""), []
            else:
                chunks.append(line.strip())
    if name is not None:
        out.append((name, desc, "".join(chunks)))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(
        f">{name}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n"
        for name, seq in records))


# ------------------------------------------------------------------ ORFs
def call_orfs(panel: Path, nulls: list[Path], out: Path) -> dict:
    """Prodigal in metagenome mode over every pile, via pyrodigal.

    pyrodigal rather than the Prodigal binary: it is the same algorithm and the same
    `meta` training, it installs from a wheel with no compiler, and it returns the CDS
    nucleotide beside the protein in one call. Arm 5 needs that nucleotide - a codon
    alignment cannot be rebuilt from a protein - and reading it back out of a GFF plus a
    contig is where an off-by-one would be invisible on a box nobody can log into.

    The shuffled pile goes through the identical call. A dinucleotide shuffle preserves
    stop-codon density, so it yields few and short ORFs: that is the honest behaviour of the
    control and the counts are recorded so the report can say how empty the null is rather
    than implying it was comparable.
    """
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    proteins: list[tuple[str, str]] = []
    rows: list[dict] = []
    for source, pile in [(panel, "real")] + [(n, "shuffled") for n in nulls]:
        for name, desc, seq in read_fasta(source):
            arm = desc.split()[0] if desc else pile
            genes = finder.find_genes(seq.encode())
            kept = 0
            tables = set()
            for gene in genes:
                aa = gene.translate().rstrip("*")
                if len(aa) < MIN_ORF_AA:
                    continue
                kept += 1
                orf = f"{name}__orf{kept}"
                tables.add(gene.translation_table)
                proteins.append((orf, aa))
                rows.append(dict(orf=orf, query=name, arm=arm, pile=pile,
                                 start=gene.begin, end=gene.end, strand=gene.strand,
                                 aa_len=len(aa), translation_table=gene.translation_table,
                                 partial=bool(gene.partial_begin or gene.partial_end),
                                 cds_nt=gene.sequence()))
            print(f"  {name:10s} {arm:18s} {len(seq):6d} nt -> {kept} ORFs >= "
                  f"{MIN_ORF_AA} aa (table {sorted(tables) or '-'}, "
                  f"bin {genes.metagenomic_bin.description if genes.metagenomic_bin else '-'})",
                  flush=True)
    write_fasta(out / "orfs.faa", proteins)
    (out / "orfs.json").write_text(json.dumps(rows, indent=1))
    by_pile: dict[str, int] = {}
    for row in rows:
        by_pile[row["pile"]] = by_pile.get(row["pile"], 0) + 1
    return dict(orfs=len(rows), by_pile=by_pile, residues=sum(r["aa_len"] for r in rows))


# ------------------------------------------------------------------ tabular arms
def parse_tab(path: Path, arm: str, program: str, db: str) -> list[dict]:
    rows = []
    if not path.exists():
        return rows
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(TAB_FIELDS):
                continue
            row = dict(zip(TAB_FIELDS, parts))
            rows.append(dict(
                arm=arm, program=program, db=db, orf=row["qseqid"],
                subject=row["sseqid"], pident=float(row["pident"]),
                aln=int(row["length"]), qlen=int(row["qlen"]), slen=int(row["slen"]),
                evalue=float(row["evalue"]), bitscore=float(row["bitscore"]),
                qcov=float(row["qcovhsp"]) / 100.0,
                qstart=int(row["qstart"]), qend=int(row["qend"]),
                sstart=int(row["sstart"]), send=int(row["send"])))
    return rows


def arm_diamond(bins: dict, dbs: dict, out: Path, threads: int, mode: str) -> list[dict]:
    """The production search, at one of its sensitive modes. The baseline to beat.

    `--max-target-seqs 500` and no E-value floor below 10: the question is the most distant
    relative recovered, so a floor drawn at the production operating point would answer a
    different question and flatter the baseline. The floor is applied in the report, to
    every arm at once.
    """
    tag = mode.replace("-", "_")
    tsv = out / f"diamond_{tag}.tsv"
    sh([bins["diamond"], "blastp", "--query", out / "orfs.faa", "--db", dbs["dmnd"],
        "--out", tsv, "--outfmt", "6", *TAB_FIELDS, "--threads", threads,
        f"--{mode}", "--max-target-seqs", 500, "--evalue", 10, "--quiet"])
    return parse_tab(tsv, f"diamond_{tag}", f"diamond blastp --{mode}", "swissprot+viral")


def arm_blastp(bins: dict, dbs: dict, out: Path, threads: int) -> list[dict]:
    """BLASTP on the same subjects. Still seeded, but with the exact Karlin-Altschul
    statistics and composition-based score adjustment that DIAMOND approximates.

    This is the middle term of arm 2, and it is the one that separates two different
    explanations for a reach limit: if BLASTP finds what DIAMOND missed, the limit was
    seeding; if it does not, the limit is the significance statistics, and no amount of
    exhaustive alignment will move it. `-comp_based_stats 1` is on because it is what makes
    an E-value against a compositionally biased query honest, and viral coat proteins are
    compositionally biased.
    """
    tsv = out / "blastp.tsv"
    sh([bins["blastp"], "-query", out / "orfs.faa", "-db", dbs["blast_prot"],
        "-out", tsv, "-outfmt", "6 " + " ".join(TAB_FIELDS), "-num_threads", threads,
        "-max_target_seqs", 500, "-evalue", 10, "-comp_based_stats", 1])
    return parse_tab(tsv, "blastp", "blastp -comp_based_stats 1", "swissprot+viral")


#: Shared state for the exhaustive arm's worker processes. A module-level pair rather than a
#: closure because `multiprocessing` pickles the callable it is given, and a closure over a
#: 230 MB subject dictionary is both unpicklable and the wrong thing to send: under `fork`
#: the children inherit this already parsed, at no copy cost.
_SW: dict = {}


def _sw_one_query(item: tuple[str, str]) -> list[dict]:
    """Every subject against one query, scored then traced back.

    The score pass is vectorised but single-threaded per call, so the parallelism is over
    queries - one process per core, each sweeping the whole subject set. That is the
    opposite of how a seeded searcher spends its threads, and it is the honest shape here:
    with no index there is nothing to share between queries, which is itself part of what
    the shortcut buys.
    """
    import math

    import parasail

    qname, qseq = item
    subjects = _SW["subjects"]
    matrix = parasail.blosum62
    lam, kay = _SW["lam"], _SW["kay"]
    started = time.time()
    scored: list[tuple[int, str]] = []
    for sname, sseq in subjects.items():
        score = parasail.sw_striped_16(qseq, sseq, 11, 1, matrix).score
        # 40 is below any floor this project has ever drawn, so nothing discarded here
        # could have been reported. It exists to keep the traceback pass small.
        if score >= 40:
            scored.append((score, sname))
    scored.sort(key=lambda t: (-t[0], t[1]))
    out_rows = []
    for score, sname in scored[:500]:
        # Traceback only for the survivors: it costs about an order of magnitude more than
        # the score alone, and only reported rows need identity and coverage.
        res = parasail.sw_trace_striped_16(qseq, subjects[sname], 11, 1, matrix)
        tb = res.get_traceback()
        q_aln, s_aln = tb.query, tb.ref
        aln = len(q_aln)
        ident = sum(1 for a, b in zip(q_aln, s_aln) if a == b and a != "-")
        covered = sum(1 for a in q_aln if a != "-")
        out_rows.append(dict(
            arm="sw_exhaustive", program="parasail sw_striped BLOSUM62 11/1",
            db="swissprot+viral", orf=qname, subject=sname,
            pident=100.0 * ident / max(aln, 1), aln=aln, qlen=len(qseq),
            slen=len(subjects[sname]),
            evalue=_SW["evalue_of"](score, len(qseq)),
            bitscore=(lam * score - math.log(kay)) / math.log(2),
            qcov=covered / max(len(qseq), 1),
            qstart=0, qend=0, sstart=0, send=0, raw_score=int(score)))
    print(f"    {qname} {time.time() - started:.1f}s {len(scored)} scores >= 40",
          flush=True)
    return out_rows


def arm_exhaustive(dbs: dict, out: Path, threads: int) -> list[dict]:
    """Full Smith-Waterman against every subject, no seeding at all.

    THE POINT OF THIS ARM IS TO BE OUTCLASSED, OR NOT. Both searches above find a candidate
    by matching short exact or near-exact words first; neither can report an alignment whose
    best local path contains no such word. Removing the seed entirely is the only way to
    measure what that costs, and the cost is widely assumed to be small without being
    measured on material at this distance.

    The score is compared on bitscore, computed from the raw score with BLOSUM62's
    published lambda and K for the same gap costs BLASTP uses, so a row here is on the same
    scale as a row from the two arms above. The E-value uses the same Karlin-Altschul form
    and the same effective search space, so an alignment that only this arm finds can be
    read as significant or not on the baseline's own terms rather than on a new one.

    Cost is the result as much as the hits are: the subject set is about 2.3e8 residues and
    the ORF pile a few thousand, so this is a few trillion cell updates, which is minutes
    with vectorised striped alignment and would be hours in any scalar implementation. That
    ratio is why nobody runs it, and it is cheap enough here to be worth knowing.
    """
    import math

    import parasail

    # BLOSUM62 with gap open 11 extend 1: the BLASTP defaults, and the parameter set
    # Altschul et al. publish lambda and K for. Hard-coded rather than fitted, because a
    # fitted lambda on a few thousand residues of query would be noise dressed as rigour.
    lam, kay = 0.267, 0.041
    matrix = parasail.blosum62
    subjects = {name: seq for name, _, seq in read_fasta(dbs["prot_fasta"])}
    total_subject_residues = sum(len(s) for s in subjects.values())
    n_subjects = len(subjects)
    queries = read_fasta(out / "orfs.faa")
    print(f"  exhaustive: {len(queries)} ORFs x {n_subjects} subjects "
          f"({total_subject_residues / 1e6:.1f} M residues)", flush=True)

    def evalue_of(score: int, qlen: int) -> float:
        """BLAST's own E-value for a database search, from the raw score.

        Two corrections matter and both are BLAST's, not inventions here. The edge
        correction shortens the effective query and subject, because an alignment cannot
        start arbitrarily close to either end and ignoring that overstates significance on
        short sequences - which is the regime a 130 aa coat protein sits in. The search
        space is then the whole database rather than one subject, so the E-value a row
        carries is directly comparable to the BLASTP row for the same pair.
        """
        avg = total_subject_residues / max(n_subjects, 1)
        corr = math.log(kay * qlen * avg) / lam
        m_prime = max(qlen - corr, 1 / kay)
        n_prime = max(total_subject_residues - n_subjects * corr, 1 / kay)
        return kay * m_prime * n_prime * math.exp(-lam * score)

    _SW.update(subjects=subjects, lam=lam, kay=kay, evalue_of=evalue_of)
    work = [(name, seq) for name, _, seq in queries]
    if threads > 1:
        import multiprocessing

        with multiprocessing.get_context("fork").Pool(threads) as pool:
            batches = pool.map(_sw_one_query, work, chunksize=1)
    else:
        batches = [_sw_one_query(item) for item in work]
    return [row for batch in batches for row in batch]


#: Zero-based column positions of `--domtblout`, named because two of them are the ones a
#: reader gets wrong. HMMER's "hmm coord" pair is on the QUERY and its "ali coord" pair is on
#: the TARGET - and which of those is our ORF depends on the program. hmmsearch searches a
#: profile library against our sequences, so the ORF is the target; jackhmmer searches our
#: sequence against a sequence database, so the ORF is the query. Using one mapping for both
#: silently puts the profile's length in the ORF's column, which the first run did: it
#: produced query coverages above 1.0, and 1.156 is the only reason it was noticed.
DOM_TARGET_NAME, DOM_TLEN = 0, 2
DOM_QUERY_NAME, DOM_QLEN = 3, 5
DOM_FULL_EVALUE, DOM_FULL_SCORE = 6, 7
DOM_C_EVALUE, DOM_I_EVALUE, DOM_DOM_SCORE = 11, 12, 13
DOM_HMM_FROM, DOM_HMM_TO = 15, 16          # on the query
DOM_ALI_FROM, DOM_ALI_TO = 17, 18          # on the target


def parse_domtbl(path: Path, *, arm: str, program: str, db: str,
                 orf_is_target: bool, accession_as_subject: bool = False) -> list[dict]:
    """One `--domtblout` file, with the ORF's own length and span in the ORF's columns.

    `orf_is_target` is the whole point of this function existing. Everything reported about
    the query here - `qlen`, `qstart`, `qend` and therefore `qcov` - is about our ORF, and
    everything reported about the subject is about the profile or the database sequence,
    whichever side that was. A coverage above 1.0 is the symptom of getting this wrong, so
    the caller is handed a value that cannot exceed 1.0 unless HMMER itself is inconsistent.
    """
    rows = []
    if not path.exists():
        return rows
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            if orf_is_target:
                orf, subject = f[DOM_TARGET_NAME], f[DOM_QUERY_NAME]
                orf_len, subject_len = int(f[DOM_TLEN]), int(f[DOM_QLEN])
                orf_from, orf_to = int(f[DOM_ALI_FROM]), int(f[DOM_ALI_TO])
                sub_from, sub_to = int(f[DOM_HMM_FROM]), int(f[DOM_HMM_TO])
                accession = f[4]
            else:
                orf, subject = f[DOM_QUERY_NAME], f[DOM_TARGET_NAME]
                orf_len, subject_len = int(f[DOM_QLEN]), int(f[DOM_TLEN])
                orf_from, orf_to = int(f[DOM_HMM_FROM]), int(f[DOM_HMM_TO])
                sub_from, sub_to = int(f[DOM_ALI_FROM]), int(f[DOM_ALI_TO])
                accession = f[1]
            span = orf_to - orf_from + 1
            rows.append(dict(
                arm=arm, program=program, db=db, orf=orf,
                subject=(accession if accession_as_subject and accession != "-"
                         else subject),
                profile_name=(subject if orf_is_target else None),
                pident=None, aln=span, qlen=orf_len, slen=subject_len,
                evalue=float(f[DOM_FULL_EVALUE]),
                domain_evalue=float(f[DOM_I_EVALUE]),
                conditional_evalue=float(f[DOM_C_EVALUE]),
                bitscore=float(f[DOM_FULL_SCORE]),
                domain_bitscore=float(f[DOM_DOM_SCORE]),
                qcov=span / max(orf_len, 1),
                qstart=orf_from, qend=orf_to, sstart=sub_from, send=sub_to))
    return rows


def arm_hmmsearch(bins: dict, dbs: dict, out: Path, threads: int) -> list[dict]:
    """Pfam-A at the curators' gathering cutoffs, and every domain rather than the best.

    WHY EVERY DOMAIN AND EVERY CLAN MEMBER. Callanan et al. (Sci Adv 2020) and Rumnieks et
    al. (Protein Sci 2026) both report that one profile per protein does not find the
    divergent members of this family, and that several profiles per role are needed. A
    best-hit-per-ORF readout would hide exactly that, so this keeps every domain above the
    cutoff and the report groups them by Pfam clan: whether a second profile in the same
    clan reaches further than the first is then a number rather than an expectation.

    `--cut_ga` is the operating point, not an E-value, and that is deliberate: measured in
    `results/profile_reach/`, the gathering thresholds gave identical figures across two
    passes while per-domain E-values moved by up to 0.42 points when the target set grew.
    GA is comparable between runs; an E-value cut is not. The E-value is recorded on every
    row anyway so a curve can be drawn.
    """
    rows = []
    for tag, flag in (("ga", "--cut_ga"), ("open", "-E")):
        dom = out / f"pfam_{tag}.domtbl"
        cmd = [bins["hmmsearch"], "--domtblout", dom, "--cpu", threads, "--noali"]
        cmd += [flag] if tag == "ga" else [flag, "10"]
        cmd += [dbs["pfam"], out / "orfs.faa"]
        sh(cmd, stdout=subprocess.DEVNULL)
        # The profile library is the query and our ORFs are the targets, which is the
        # opposite way round from jackhmmer below.
        rows += parse_domtbl(dom, arm=f"hmmsearch_{tag}", program=f"hmmsearch {flag}",
                             db="Pfam-A", orf_is_target=True, accession_as_subject=True)
    return rows


def arm_jackhmmer(bins: dict, dbs: dict, out: Path, threads: int) -> list[dict]:
    """An iterated profile built from the query, searched against the SAME subject set.

    This is the arm that separates "profiles help" from "Pfam helps". hmmsearch above
    searches a curated profile library, which is a different subject universe from
    Swiss-Prot and therefore cannot be differenced against the sequence arms. jackhmmer
    builds its profile out of the query's own hits and searches Swiss-Prot, so the subjects
    are identical to arms 1 and 2 and the only thing that changed is that a profile replaced
    a sequence. Any extra reach here is attributable.

    Five iterations. Three is the common default and the literature on this family reports
    reach still growing past it; five is where it is cheap enough to just run and see, and
    the per-iteration counts are kept so the curve is visible rather than a single endpoint.
    """
    tbl = out / "jackhmmer.tbl"
    # `-A` keeps the final iteration's alignment per query. That file is this run's only
    # source of a query profile built from real homologues, and the profile-to-profile arm
    # is built on it rather than on an hhblits pass against UniRef30 - which would be a 50
    # GB download and, more to the point, a different subject universe from the one every
    # other arm here searches.
    sh([bins["jackhmmer"], "-N", 5, "--cpu", threads, "--noali",
        "--domtblout", tbl, "-A", out / "jackhmmer.sto", "-E", 10, "--incE", 1e-3,
        out / "orfs.faa", dbs["prot_fasta"]], stdout=subprocess.DEVNULL)
    # Our ORF is the query here and the database sequences are the targets, so the
    # coordinate pairs mean the opposite of what they mean for hmmsearch above.
    return parse_domtbl(tbl, arm="jackhmmer", program="jackhmmer -N 5",
                        db="swissprot+viral", orf_is_target=False)


def arm_hhsearch(bins: dict, dbs: dict, out: Path, threads: int) -> list[dict]:
    """Profile against profile: the query's own HMM searched against Pfam's.

    WHY THIS IS A SEPARATE ARM FROM hmmsearch. hmmsearch compares a curated profile to the
    query's single sequence, so one side of the comparison carries no information about the
    query's family. HHsearch compares two profiles, which is the more sensitive instrument
    on paper and the one the literature on this family leans on. If profiles help here at
    all, this is where the help should be largest - and if it is not larger than hmmsearch,
    that is the measurement.

    WHERE THE QUERY PROFILE COMES FROM, AND WHY NOT hhblits. The standard recipe builds the
    query profile by iterating against UniRef30, a 50 GB download. That would also change
    the subject universe behind the query side, which is the one thing this comparison holds
    fixed: every other arm here derives whatever it knows about the query from Swiss-Prot
    plus RefSeq viral. So the profile is built from jackhmmer's own final alignment against
    exactly that set. The cost is a shallower profile than hhblits would give, which makes
    this arm's reach a lower bound rather than the method's ceiling, and that is said in the
    report rather than left for a reader to notice.

    The target side is Pfam's HHM database, so this arm's subjects are Pfam and no
    difference against the sequence baseline is computed for it.
    """
    sto = out / "jackhmmer.sto"
    if not sto.exists():
        raise FileNotFoundError(
            f"{sto} is missing: the profile arm needs jackhmmer's alignment, so the "
            f"jackhmmer arm has to run before it. Pass both, in that order.")

    # jackhmmer writes one Stockholm block per query, concatenated. Split on the terminator
    # rather than parsing: each block names its own query in the `#=GS` lines and hhmake
    # reads Stockholm directly, so nothing needs reformatting.
    blocks, current = [], []
    for line in sto.read_text().splitlines(keepends=True):
        current.append(line)
        if line.strip() == "//":
            blocks.append("".join(current))
            current = []
    print(f"  {len(blocks)} query alignments from jackhmmer", flush=True)

    hhm_dir = out / "hhm"
    hhm_dir.mkdir(exist_ok=True)
    queries = [name for name, _, _ in read_fasta(out / "orfs.faa")]
    pairs: list[tuple[str, Path]] = []
    for i, block in enumerate(blocks):
        # The query an alignment belongs to is the sequence jackhmmer searched with, which
        # appears in the block's `#=GS` or `#=GF` metadata depending on the version. Falling
        # back to file order is safe because jackhmmer emits one block per query in query
        # order, and the assertion below catches it if that ever stops being true.
        named = [q for q in queries if q in block]
        qname = named[0] if named else (queries[i] if i < len(queries) else f"block{i}")
        path = hhm_dir / f"{qname}.sto"
        path.write_text(block)
        pairs.append((qname, path))

    rows = []
    for qname, path in pairs:
        hhm = path.with_suffix(".hhm")
        try:
            sh([bins["hhmake"], "-i", path, "-o", hhm, "-M", "first", "-v", 0])
        except subprocess.CalledProcessError:
            # An alignment of one sequence is a profile with no information; hhmake refuses
            # it. That is a real outcome for an ORF jackhmmer found nothing for, and it is
            # recorded by the ORF's absence from these rows rather than by a crash.
            print(f"    {qname}: no profile (alignment too shallow)", flush=True)
            continue
        hhr = path.with_suffix(".hhr")
        sh([bins["hhsearch"], "-i", hhm, "-d", dbs["hh_pfam"], "-o", hhr,
            "-cpu", threads, "-e", 10, "-v", 0])
        for hit in _parse_hhr(hhr):
            rows.append(dict(
                arm="hhsearch", program="hhsearch, query HMM from jackhmmer", db="Pfam-A",
                orf=qname, subject=hit["subject"], pident=None, aln=hit["cols"],
                qlen=hit["qlen"], slen=hit["slen"], evalue=hit["evalue"],
                bitscore=hit["score"], probability=hit["prob"],
                qcov=hit["cols"] / max(hit["qlen"], 1),
                qstart=hit["qstart"], qend=hit["qend"],
                sstart=hit["sstart"], send=hit["send"]))
    return rows


def _parse_hhr(path: Path) -> list[dict]:
    """The hit table of one `.hhr`, which is fixed-width and not tabular.

    Parsed from the summary block rather than the alignments, and by the column positions
    hhsearch documents rather than by splitting on whitespace: a long Pfam description
    contains spaces and would shift every field after it.
    """
    text = path.read_text().splitlines()
    qlen = 0
    for line in text:
        if line.startswith("Match_columns"):
            qlen = int(line.split()[-1])
            break
    out, started = [], False
    for line in text:
        if line.startswith(" No Hit"):
            started = True
            continue
        if started:
            if not line.strip():
                break
            # ` No Hit  Prob E-value P-value  Score  SS Cols Query HMM Template HMM`
            tail = line[34:].split()
            if len(tail) < 7:
                continue
            subject = line[4:34].split()[0]
            try:
                prob, evalue, _pvalue, score = (float(tail[0]), float(tail[1]),
                                                float(tail[2]), float(tail[3]))
                cols = int(tail[5])
                qrange = tail[6].split("-")
                srange = tail[7].split("-") if len(tail) > 7 else ["0", "0"]
                slen_field = tail[8] if len(tail) > 8 else "(0)"
            except (ValueError, IndexError):
                continue
            out.append(dict(
                subject=subject, prob=prob, evalue=evalue, score=score, cols=cols,
                qlen=qlen, qstart=int(qrange[0]), qend=int(qrange[-1]),
                sstart=int(srange[0]), send=int(srange[-1]),
                slen=int(slen_field.strip("()") or 0)))
    return out


def arm_foldseek(bins: dict, dbs: dict, out: Path, threads: int, target: str) -> list[dict]:
    """Structure comparison, with the structures derived from sequence by ProstT5.

    WHY NOT FOLD THEM. Predicting ~35 structures with a folding model means a GPU instance
    and, for an MSA-based predictor, the MSA databases too - about $4 before anything is
    searched, and the search itself is then the same Foldseek call as here. ProstT5 writes
    the 3Di structural alphabet straight from sequence, which is what Foldseek actually
    aligns, so this arm asks the structure question at CPU cost. The model was already
    staged in this account's S3 (`corpora/models/prostt5/`), so it cost nothing to try.

    WHAT THAT COSTS IN HONESTY. A predicted 3Di string is not a structure and ProstT5 was
    trained on sequences with known folds, so on a protein with no structural relative in
    its training set the alphabet it writes is an informed guess. If this arm reaches
    further than the sequence arms, the folded version is worth the $4; if it does not, that
    is a weaker negative than a folded run would give, and the report says so rather than
    reading it as a kill.

    AlphaFold/Swiss-Prot is the target that makes this comparable: it is Swiss-Prot folded,
    so arm 4 against it is arms 1-2's subjects in a different alphabet. PDB is run as well
    because it is where a remote structural relative of a phage coat protein would actually
    be, and it is reported as its own universe.
    """
    qdb, resdb = out / f"fsq_{target}", out / f"fsr_{target}"
    tsv = out / f"foldseek_{target}.tsv"
    if not qdb.with_suffix(".dbtype").exists():
        # `--prostt5-model` wants the DIRECTORY foldseek's own `databases ProstT5` would
        # have written, not the weights file inside it. Passing the file is what the first
        # run did, and foldseek accepted the option, printed it back, and then exited 1 - so
        # the file form is tried as a fallback rather than assumed wrong, and whichever
        # worked goes in the row.
        model = dbs["prostt5"]
        try:
            sh([bins["foldseek"], "createdb", out / "orfs.faa", qdb,
                "--prostt5-model", model, "--threads", threads])
        except subprocess.CalledProcessError:
            inner = next(iter(sorted(Path(model).glob("*.gguf"))), None)
            if inner is None:
                raise
            sh([bins["foldseek"], "createdb", out / "orfs.faa", qdb,
                "--prostt5-model", inner, "--threads", threads])
    # `--exhaustive-search` skips the prefilter, which is the whole point of the arm: a
    # structural relative this far out is exactly what a k-mer prefilter would drop. It is
    # affordable only because the query set is 35 ORFs. If this build of foldseek will not
    # take the flag combination, fall back to its most sensitive prefiltered setting and say
    # which ran, rather than losing the arm over an option name.
    try:
        sh([bins["foldseek"], "search", qdb, dbs[f"fs_{target}"], resdb,
            out / f"fstmp_{target}", "-a", "--threads", threads, "-e", 10,
            "--max-seqs", 1000, "--exhaustive-search", 1])
        mode = "exhaustive-search"
    except subprocess.CalledProcessError:
        shutil.rmtree(out / f"fstmp_{target}", ignore_errors=True)
        sh([bins["foldseek"], "search", qdb, dbs[f"fs_{target}"], resdb,
            out / f"fstmp_{target}", "-a", "--threads", threads, "-e", 10,
            "--max-seqs", 1000, "-s", 9.5])
        mode = "prefiltered -s 9.5"
    print(f"  foldseek {target}: {mode}", flush=True)
    sh([bins["foldseek"], "convertalis", qdb, dbs[f"fs_{target}"], resdb, tsv,
        "--threads", threads, "--format-output",
        "query,target,fident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,"
        "qcov,alntmscore,lddt"])
    rows = []
    with open(tsv) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            rows.append(dict(
                arm=f"foldseek_{target}",
                program=f"foldseek search (ProstT5 3Di, {mode})",
                db=target, orf=f[0], subject=f[1], pident=100.0 * float(f[2]),
                aln=int(f[3]), qlen=int(f[4]), slen=int(f[5]),
                qstart=int(f[6]), qend=int(f[7]), sstart=int(f[8]), send=int(f[9]),
                evalue=float(f[10]), bitscore=float(f[11]), qcov=float(f[12]),
                tmscore=float(f[13]), lddt=float(f[14])))
    return rows


def arm_codon(bins: dict, dbs: dict, out: Path, threads: int) -> list[dict]:
    """Codon-space alignment in a clean reading frame, and the saturation test.

    WHY tblastn AND NOT A PROTEIN HIT. A codon alignment needs the subject's nucleotide, not
    its protein, and resolving a protein accession back to its coding nucleotide means a
    network round trip per subject and a different answer each time the taxonomy moves.
    tblastn against RefSeq viral *genomes* gives both at once: the alignment is scored in
    translated space, and its subject coordinates are nucleotide coordinates in a genome
    this run has on disk, so each aligned residue pair can be walked back to a codon pair
    exactly. The query's own codons come from Prodigal's CDS, which is why `call_orfs` keeps
    the nucleotide.

    WHAT IS EXPECTED TO HAPPEN, SO THAT IT CAN FAIL VISIBLY. At 30-45% amino acid identity
    synonymous sites are saturated: pS approaches its asymptote, dS stops being estimable,
    and any dN/dS computed there is a ratio with a meaningless denominator. The near-
    reference arm sits at 95-99% identity where the same estimator is well behaved. **So
    this arm is run on both bands on purpose, and the control is the measurement.** If
    saturation is real, the near-reference arm gives stable per-gene dN/dS across its seven
    contigs and the divergent arm does not, and the instability is the finding. If instead
    the divergent arm is also stable, this project's premise about where codon scoring helps
    is wrong and that is worth more than a confirmation.

    Nei-Gojobori counting with Jukes-Cantor correction, which is the estimator whose failure
    mode at saturation is known and visible: the correction's argument goes non-positive and
    dS becomes undefined rather than large. Undefined is recorded as `null` and counted; a
    silently clipped dS would turn saturation into a finite number and hide the whole point.
    """
    import math

    orfs = {r["orf"]: r for r in json.loads((out / "orfs.json").read_text())}
    tsv = out / "tblastn.tsv"
    fields = ("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore", "qcovhsp", "qseq", "sseq", "sframe")
    # Skipped when the table is already there. Two reasons, and the second is the one that
    # matters: a re-run of the analysis after a fix does not have to pay for the search
    # again, and the codon walk - which is the part of this file most able to be silently
    # wrong by one - becomes testable from a hand-written table with a known answer.
    # `tests/test_reach_codon.py` is that test.
    if not tsv.exists():
        sh([bins["tblastn"], "-query", out / "orfs.faa", "-db", dbs["viral_nt"],
            "-out", tsv, "-outfmt", "6 " + " ".join(fields), "-num_threads", threads,
            "-max_target_seqs", 50, "-evalue", 10, "-db_gencode", SUBJECT_GENETIC_CODE])
    else:
        print(f"  reusing {tsv.name}: {sum(1 for _ in open(tsv))} alignments", flush=True)

    #: Nei-Gojobori needs the number of synonymous and non-synonymous sites per codon, which
    #: depends only on the codon and the table. Built from the table rather than tabulated by
    #: hand, and per table rather than once: the query side can be table 4, where TGA codes
    #: tryptophan, and counting TGA as a stop there would drop real codon pairs from the
    #: densest part of the alignment.
    BASES = "TCAG"
    CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
    STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

    def aa_map(table: int) -> dict[str, str]:
        letters = list(STANDARD)
        if table == 4:                      # TGA is tryptophan, not a stop
            letters[CODONS.index("TGA")] = "W"
        return dict(zip(CODONS, letters))

    maps = {t: aa_map(t) for t in (1, 4, 11)}

    def site_counts(codon: str, aa_of: dict[str, str]) -> tuple[float, float]:
        """Synonymous and non-synonymous sites in one codon, by the standard count."""
        if codon not in aa_of or aa_of[codon] == "*":
            return 0.0, 0.0
        syn = 0.0
        for pos in range(3):
            for base in BASES:
                if base == codon[pos]:
                    continue
                alt = codon[:pos] + base + codon[pos + 1:]
                if aa_of.get(alt) == aa_of[codon]:
                    syn += 1 / 3
        return syn, 3.0 - syn

    rows = []
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < len(fields):
            continue
        rec = dict(zip(fields, f))
        orf = orfs.get(rec["qseqid"])
        if orf is None:
            continue
        q_aa = maps.get(orf.get("translation_table", 11), maps[11])
        s_aa = maps[SUBJECT_GENETIC_CODE]
        q_aln, s_aln = rec["qseq"], rec["sseq"]
        qcds = orf["cds_nt"].upper()
        sgenome = dbs["viral_seqs"].get(rec["sseqid"], "")
        qpos = (int(rec["qstart"]) - 1) * 3          # residue -> codon start in the CDS
        spos = int(rec["sstart"]) - 1                # nucleotide, 1-based in the genome
        step = 1 if int(rec["send"]) >= int(rec["sstart"]) else -1
        syn_sites = nonsyn_sites = 0.0
        syn_diff = nonsyn_diff = 0
        codon_pairs = stop_or_gap = 0
        for qa, sa in zip(q_aln, s_aln):
            if qa == "-":
                spos += 3 * step
                continue
            if sa == "-":
                qpos += 3
                continue
            qcod = qcds[qpos:qpos + 3]
            if step > 0:
                scod = sgenome[spos:spos + 3]
            else:
                chunk = sgenome[max(spos - 2, 0):spos + 1]
                scod = chunk[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
            qpos += 3
            spos += 3 * step
            if len(qcod) != 3 or len(scod) != 3 or qcod not in q_aa or scod not in s_aa:
                stop_or_gap += 1
                continue
            if q_aa[qcod] == "*" or s_aa[scod] == "*":
                stop_or_gap += 1
                continue
            codon_pairs += 1
            # Sites averaged over the pair, which is Nei-Gojobori's convention: the two
            # codons differ, so each has its own site count and neither alone is the right
            # denominator. Synonymy of a DIFFERENCE is judged on the query's table below.
            for codon in (qcod, scod):
                s, n = site_counts(codon, q_aa)
                syn_sites += s / 2
                nonsyn_sites += n / 2
            # Differences counted per position. With more than one difference in a codon the
            # path matters; averaging over the orderings is Nei-Gojobori's own convention and
            # the error it leaves is second order next to saturation, which is what this arm
            # is measuring.
            diffs = [i for i in range(3) if qcod[i] != scod[i]]
            if not diffs:
                continue
            if len(diffs) == 1:
                i = diffs[0]
                alt = qcod[:i] + scod[i] + qcod[i + 1:]
                if q_aa.get(alt) == q_aa[qcod]:
                    syn_diff += 1
                else:
                    nonsyn_diff += 1
            else:
                # Two or three differences in one codon: the order the substitutions
                # happened in changes whether each was synonymous, and the orderings are
                # unobservable. Each differing position is scored on the query background
                # and the count is apportioned, which is the first-order version of
                # Nei-Gojobori's path averaging. The error that leaves is second order
                # beside saturation, which is what this arm is measuring - and at the
                # distances where multi-hit codons are common, dS is undefined anyway and
                # the row reports that instead of a ratio.
                syn = sum(1 for i in diffs
                          if q_aa.get(qcod[:i] + scod[i] + qcod[i + 1:]) == q_aa[qcod])
                syn_diff += syn
                nonsyn_diff += len(diffs) - syn

        def jc(p: float) -> float | None:
            """Jukes-Cantor correction. None where the argument has gone non-positive.

            That is saturation, reported as undefined rather than clipped: a clipped value
            would be a finite dS where the data cannot support one, which is exactly the
            number this arm exists to refuse to produce.
            """
            arg = 1 - 4 * p / 3
            return None if arg <= 0 else -0.75 * math.log(arg)

        p_s = syn_diff / syn_sites if syn_sites else None
        p_n = nonsyn_diff / nonsyn_sites if nonsyn_sites else None
        d_s = jc(p_s) if p_s is not None else None
        d_n = jc(p_n) if p_n is not None else None
        rows.append(dict(
            arm="codon", program="tblastn + Nei-Gojobori", db="refseq viral genomes",
            orf=rec["qseqid"], subject=rec["sseqid"], pident=float(rec["pident"]),
            aln=int(rec["length"]), qlen=int(rec["qlen"]), slen=int(rec["slen"]),
            evalue=float(rec["evalue"]), bitscore=float(rec["bitscore"]),
            qcov=float(rec["qcovhsp"]) / 100.0, qstart=int(rec["qstart"]),
            qend=int(rec["qend"]), sstart=int(rec["sstart"]), send=int(rec["send"]),
            codon_pairs=codon_pairs, skipped_codons=stop_or_gap,
            syn_sites=round(syn_sites, 3), nonsyn_sites=round(nonsyn_sites, 3),
            p_syn=p_s, p_nonsyn=p_n, d_syn=d_s, d_nonsyn=d_n,
            dnds=(d_n / d_s if (d_n is not None and d_s not in (None, 0)) else None),
            synonymous_saturated=(d_s is None)))
    return rows


# ------------------------------------------------------------------ driver
#: THE ORDER IS A BUDGET DECISION, NOT A LIST. The driver runs arms in this order whatever
#: order they were requested in, because on a box with a dead-man timer the last arm is the
#: one that does not run, and which arm that is should be chosen here rather than by however
#: the staging script happened to accumulate its variable.
#:
#: Two facts set it. First, `jackhmmer` is the expensive arm by an order of magnitude: it
#: scans the whole subject set once per iteration per query, so five iterations over 200
#: query ORFs is hours where every other arm is minutes. Second, `hhsearch` builds its query
#: profile from the alignment `jackhmmer` writes, so it cannot precede it. That puts the two
#: of them at the end, and means a timeout costs the profile pair and nothing else - rather
#: than costing structure and codon space, which is what happened on the first run.
RUN_ORDER = ("diamond", "blastp", "hmmsearch", "exhaustive", "codon", "foldseek_afdb",
             "foldseek_pdb", "jackhmmer", "hhsearch")
ARMS = RUN_ORDER


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--null", type=Path, nargs="+", required=True,
                    help="one or more shuffled piles, one per null seed. Five is the "
                         "convention: the real arms are fixed sequence, so the shuffle is "
                         "the only place a seed means anything here")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--arms", nargs="+", default=list(ARMS), choices=ARMS)
    ap.add_argument("--bin-dir", type=Path, default=Path("/opt/hapax/bin"))
    ap.add_argument("--db-dir", type=Path, default=Path("/opt/hapax/db"))
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    bins = {name: (args.bin_dir / name) for name in
            ("diamond", "blastp", "tblastn", "hmmsearch", "jackhmmer", "foldseek",
             "hhsearch", "hhmake")}
    for name, path in bins.items():
        if not path.exists():
            found = shutil.which(name)
            if found:
                bins[name] = Path(found)
    dbs = dict(
        dmnd=args.db_dir / "prot.dmnd",
        blast_prot=args.db_dir / "prot",
        prot_fasta=args.db_dir / "prot.faa",
        pfam=args.db_dir / "Pfam-A.hmm",
        prostt5=args.db_dir / "prostt5",
        fs_afdb=args.db_dir / "fs_afdb_sprot" / "afdb_sprot",
        fs_pdb=args.db_dir / "fs_pdb" / "pdb",
        viral_nt=args.db_dir / "viral_nt",
        hh_pfam=args.db_dir / "hh_pfam" / "pfam",
    )

    orf_stats = call_orfs(args.panel, list(args.null), args.out)
    print(f"=== ORFs: {json.dumps(orf_stats)} ===", flush=True)

    # Loaded once and shared: arm 5 walks subject coordinates into genome sequence, and
    # re-reading a 5 GB FASTA per alignment would dominate the arm's own cost.
    if "codon" in args.arms:
        dbs["viral_seqs"] = {n: s for n, _, s in read_fasta(args.db_dir / "viral_nt.fna")}
        print(f"=== viral genomes in memory: {len(dbs['viral_seqs'])} ===", flush=True)

    plan = {
        "diamond": lambda: (arm_diamond(bins, dbs, args.out, args.threads, "sensitive")
                            + arm_diamond(bins, dbs, args.out, args.threads,
                                          "very-sensitive")
                            + arm_diamond(bins, dbs, args.out, args.threads,
                                          "ultra-sensitive")),
        "blastp": lambda: arm_blastp(bins, dbs, args.out, args.threads),
        "exhaustive": lambda: arm_exhaustive(dbs, args.out, args.threads),
        "hmmsearch": lambda: arm_hmmsearch(bins, dbs, args.out, args.threads),
        "jackhmmer": lambda: arm_jackhmmer(bins, dbs, args.out, args.threads),
        "hhsearch": lambda: arm_hhsearch(bins, dbs, args.out, args.threads),
        "foldseek_afdb": lambda: arm_foldseek(bins, dbs, args.out, args.threads, "afdb"),
        "foldseek_pdb": lambda: arm_foldseek(bins, dbs, args.out, args.threads, "pdb"),
        "codon": lambda: arm_codon(bins, dbs, args.out, args.threads),
    }

    failures, succeeded = {}, []
    ordered = [a for a in RUN_ORDER if a in args.arms]
    print(f"=== running in order: {' '.join(ordered)} ===", flush=True)
    for arm in ordered:
        started = time.time()
        timing = dict(arm=arm, seed=args.seed, threads=args.threads,
                      host=platform.node(), cpu_count=os.cpu_count())
        try:
            rows = plan[arm]()
            elapsed = time.time() - started
            with open(args.out / f"rows_{arm}.jsonl", "w") as handle:
                for row in rows:
                    handle.write(json.dumps(row) + "\n")
            timing |= dict(seconds=round(elapsed, 2), rows=len(rows), ok=True)
            succeeded.append(arm)
            print(f"=== {arm}: {len(rows)} rows in {elapsed:.1f}s ===", flush=True)
        except Exception:
            timing |= dict(seconds=round(time.time() - started, 2), ok=False,
                           traceback=traceback.format_exc())
            failures[arm] = traceback.format_exc()
            print(f"=== {arm} FAILED after {timing['seconds']}s ===\n"
                  f"{timing['traceback']}", flush=True)
        (args.out / f"timing_{arm}.json").write_text(json.dumps(timing, indent=1))

    (args.out / "arms_manifest.json").write_text(json.dumps(dict(
        seed=args.seed, arms_requested=args.arms, arms_run_order=ordered,
        arms_ok=succeeded,
        arms_failed=sorted(failures), orf_stats=orf_stats, threads=args.threads,
        min_orf_aa=MIN_ORF_AA, subject_genetic_code=SUBJECT_GENETIC_CODE,
        query_translation_table="per sequence, chosen by prodigal meta mode; recorded on "
                                 "each ORF row",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    ), indent=1))
    if not succeeded:
        print("=== every arm failed; the setup is wrong rather than one tool ===")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
