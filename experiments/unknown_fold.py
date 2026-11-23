#!/usr/bin/env python3
"""Can structural comparison supply an identity for the pile where sequence supplies none?

WHY THIS EXISTS, AND WHY IT IS NOT A RE-RUN OF THE REACH ARM. `results/reach/`'s structure
arm pointed ProstT5 plus Foldseek at seven divergent *Fiersviridae* genomes, which sequence
search can already name. It placed 8 of 24 genes where sequence placed 22, and the 8 were a
subset, so on that panel it added no identity. That is a fact about material sequence can
name. This experiment points the same instrument at the pile where sequence returns no name
at all - the cascade's 46,157 unidentifiable contigs - because "new" without "and it is a
polymerase from a family that infects X" is a sentence with its second half missing, and
structure is the only instrument we hold that reads something sequence does not.

Expect a low number. If structure placed a third of the genes sequence could name, it will
place fewer where sequence names none. A number near zero is the result, not a failure, and
it ships either way.

WHAT THE SUBJECT LIBRARY IS, AND WHY IT IS NOT THE REACH ARM'S. The reach arm's library was
every solved *virus* structure (27,439 entries). For a pile nobody has named, that library
can only answer "is this like a virus someone has solved", and a blank from it would read as
"structure does not work" when it meant "the relative was not a virus". Of the nameable
material in this sample 76% is bacterial and 8% viral, so the library here is the whole
databank, made non-redundant by RCSB's own 30%-sequence-identity clustering: one
representative per cluster, 40,066 distinct entries covering 272,663 clusters. Non-redundant
rather than complete on purpose - 263,157 experimental entries is 100 GB and mostly the same
folds over again, and a cluster representative set is the standard way to search fold space
once.

THE CONTROLS, READ FIRST.

*The comparator arm* is the identified contigs from the same assembly of the same library,
ribosomal ones removed, run through the identical path. Without it a placement rate is
unreadable: it says nothing about the pile unless you know what the instrument does on
material from the same water that sequence can name.

*The null* is a dinucleotide-preserving shuffle of each sampled contig, from
`shuffle_control.shuffle_arm`, ORF-called through the same caller and searched in the same
pass. It carries its census - offered, shuffled, dropped_identical - because on this project
a shuffle that quietly handed back its input supplied 100% of a null's hits.

*The matching* is the control the reach arm did not need and this one does. The two piles
are not comparable in ORF length: the unidentifiable pile's ORFs run 78 aa at the median
against the identified pile's 105, and only 23% of them reach 100 aa against 53%. Alignment
statistics are a function of query length, so an unmatched comparison would measure contig
length and report it as novelty. Every headline here is reported inside ORF-length bands,
and the overall figure is the band-weighted one with the raw one printed beside it.

WHAT "PLACED" MEANS. A query ORF is placed when it has an alignment to a subject chain that
clears the E-value floor. The rows are written un-thresholded and the readout applies the
floor, so a different floor can be read off the same run. No TM-score and no LDDT is
reportable from a ProstT5 query: both are superposition scores and the query side has no
C-alpha, which is the whole reason this runs on CPU in a container instead of on a GPU.

AND THEN WHAT IT IS. A placement rate on its own is only the first half of what Imo asked
for. For every placed ORF the subject entry's own description and source organism are
fetched from RCSB and carried into the rows, so the output says what the thing is like and
not merely that something matched. Those names are a database record's words about a
database entry, never a statement about the water.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import random
import statistics
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "experiments"))

from provenance import provenance  # noqa: E402
from shuffle_control import ShuffleFailed, shuffle_arm  # noqa: E402

#: Minimum ORF length in residues, the same floor `reach_arms.call_orfs` uses, so the two
#: structure arms are talking about the same kind of object.
MIN_ORF_AA = 50

#: ORF-length bands the comparison is reported inside. The edges are the quantiles where the
#: two piles' ORF-length distributions diverge, chosen once from the whole-pile ORF census
#: and then fixed, not tuned per run.
BANDS = [(50, 80), (80, 100), (100, 150), (150, 200), (200, 10**9)]

#: The E-value the search is run at. Rows are written un-thresholded down to this and the
#: readout applies its own floor, so the floor is a readout parameter and not a run one.
EVALUE_SEARCHED = 10.0

#: The floor the headline is read at, matching every other arm in this project.
EVALUE_FLOOR = 1e-5

RCSB_CLUSTERS = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt"
RCSB_FILE = "https://files.rcsb.org/download/{}.cif.gz"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
FETCH_PARALLEL = 8


# ----------------------------------------------------------------- small shared helpers
def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    """Run a command, keeping its own complaint.

    `check=True` on its own throws the tool's output away, and the exit code carries none of
    it. That is how a foldseek segfault spent four runs of the reach arm looking like a bare
    exit 1.
    """
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def read_fasta(path: Path):
    name = None
    buf: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def write_fasta(path: Path, items) -> None:
    with path.open("w") as fh:
        for name, seq in items:
            fh.write(f">{name}\n{seq}\n")


def band_of(aa_len: int) -> str:
    for lo, hi in BANDS:
        if lo <= aa_len < hi:
            return f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
    raise ValueError(f"ORF of {aa_len} aa falls in no band; BANDS must cover MIN_ORF_AA up")


# ----------------------------------------------------------------- the subject library
def cluster_representatives(work: Path) -> dict:
    """One PDB entry per 30%-identity cluster, from RCSB's own cluster file.

    The file is a line per cluster, members ordered by RCSB's own ranking, so the first
    token is the representative. Taking entries rather than entities means a representative
    drags its fellow chains in, which is wanted: a capsid's neighbours in the same
    crystal are legitimate subjects.
    """
    path = work / "clusters-30.txt"
    if not path.exists() or path.stat().st_size < 10_000_000:
        sh(["curl", "-sS", "--max-time", "600", "-o", path, RCSB_CLUSTERS])
    clusters = 0
    entries: set[str] = set()
    for line in path.open():
        tok = line.split()
        if not tok:
            continue
        clusters += 1
        entries.add(tok[0].split("_")[0])
    if clusters < 100_000:
        raise RuntimeError(
            f"the RCSB cluster file parsed to {clusters} clusters, which is an order of "
            f"magnitude short. A short subject library would read as structure finding "
            f"fewer relatives rather than as a short library, so this refuses instead.")
    ids = sorted(entries)
    (work / "pdb30_entry_ids.txt").write_text("\n".join(ids) + "\n")
    return dict(source=RCSB_CLUSTERS, clusters=clusters, representative_entries=len(ids),
                identity_cutoff=0.30,
                fetched_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))


def fetch_structures(ids_file: Path, cif_dir: Path) -> dict:
    """Download each representative's mmCIF, and decide on the files rather than the exit code.

    `foldseek databases` exits 0 when its download fails, untarring a 17-byte HTTP 429 body
    and reporting success, which produced three runs of the reach arm with no subjects at
    all. So nothing here trusts a return code: the resume test is a non-empty file, and the
    gzip magic is checked on a sample before the database is built.
    """
    cif_dir.mkdir(parents=True, exist_ok=True)
    script = (f'test -s {cif_dir}/$1.cif.gz || '
              f'curl -sS --max-time 90 --retry 3 -o {cif_dir}/$1.cif.gz '
              f'"{RCSB_FILE.format("$1")}"')
    started = time.time()
    subprocess.run(["xargs", "-P", str(FETCH_PARALLEL), "-I", "{}",
                    "sh", "-c", script.replace("$1", "{}"), "{}"],
                   stdin=ids_file.open(), check=False)
    wanted = [i for i in ids_file.read_text().split() if i]
    present = [i for i in wanted if (cif_dir / f"{i}.cif.gz").exists()
               and (cif_dir / f"{i}.cif.gz").stat().st_size > 0]
    missing = [i for i in wanted if i not in set(present)]
    bad = []
    for entry in present[:: max(1, len(present) // 200)]:
        with (cif_dir / f"{entry}.cif.gz").open("rb") as fh:
            if fh.read(2) != b"\x1f\x8b":
                bad.append(entry)
    if bad:
        raise RuntimeError(
            f"{len(bad)} of the sampled downloads are not gzip, so the server served "
            f"something other than a structure: {bad[:5]}. A rate-limit page untars to "
            f"nothing and reads as structure finding no relative.")
    return dict(wanted=len(wanted), present=len(present), missing=len(missing),
                missing_examples=missing[:10], magic_checked=len(present[:: max(1, len(present) // 200)]),
                seconds=round(time.time() - started, 1))


def build_subject_db(foldseek: Path, cif_dir: Path, db: Path, threads: int) -> dict:
    """Foldseek database from real coordinates, which is what makes it a structural target.

    `createdb` over mmCIF reads C-alpha and derives 3Di from geometry. The subject side
    therefore has coordinates and the ProstT5 query side does not, which is exactly why no
    superposition score is reportable from this arm.
    """
    db.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    done = sh([foldseek, "createdb", cif_dir, db, "--threads", threads],
              capture_output=True, text=True)
    log = (done.stdout or "") + (done.stderr or "")
    (db.parent / "createdb_subjects.log").write_text(log)
    index = db.with_suffix(".index")
    if not index.exists():
        raise RuntimeError("foldseek createdb wrote no subject index, so this arm has no "
                           "subjects. Its own output:\n" + log[-2500:])
    chains = sum(1 for _ in index.open())
    if chains < 50_000:
        raise RuntimeError(
            f"the subject database holds {chains} chains against 40,000 entries offered. "
            f"That is a staging failure, and an empty table would read as a result.")
    return dict(chains=chains, seconds=round(time.time() - started, 1))


# ----------------------------------------------------------------- the query side
def sample_blocks(pile: Path, seeds: list[int], per_seed: int) -> dict:
    """Disjoint random blocks of contigs, one per seed.

    DISJOINT, not five draws from the same pile. Overlapping samples would share contigs and
    the spread across seeds would understate itself, which on this project is the mistake
    that reversed a statistic's sign once already. The draw is one shuffle of the pile's
    contig names under a fixed seed, then consecutive slices, so a seed's membership is
    reproducible from this file plus the pile.
    """
    names = [n for n, _ in read_fasta(pile)]
    if len(names) < per_seed * len(seeds):
        raise RuntimeError(
            f"{pile.name} holds {len(names)} contigs, too few for {len(seeds)} disjoint "
            f"blocks of {per_seed}. Reduce --contigs-per-seed rather than letting the "
            f"blocks overlap silently.")
    rng = random.Random(20261123)
    order = list(names)
    rng.shuffle(order)
    return {seed: order[i * per_seed:(i + 1) * per_seed] for i, seed in enumerate(seeds)}


def call_orfs(records, pile: str, arm: str, rows: list[dict], proteins: list[tuple[str, str]],
              finder) -> None:
    """pyrodigal in metagenome mode, the same caller and the same floor as the reach arm.

    The null goes through this identical call. A dinucleotide shuffle preserves stop-codon
    density, so it yields fewer and shorter ORFs; that is the control behaving honestly and
    the counts are recorded so the readout can say how thin the null is rather than implying
    it was matched.
    """
    for name, seq in records:
        kept = 0
        for gene in finder.find_genes(seq.encode()):
            aa = gene.translate().rstrip("*")
            if len(aa) < MIN_ORF_AA:
                continue
            kept += 1
            orf = f"{name}__{arm}__orf{kept}"
            proteins.append((orf, aa))
            rows.append(dict(orf=orf, contig=name, pile=pile, arm=arm, aa_len=len(aa),
                             band=band_of(len(aa)), start=gene.begin, end=gene.end,
                             strand=gene.strand, partial=bool(gene.partial_begin or gene.partial_end)))


def build_queries(piles: dict[str, Path], seeds: list[int], per_seed: int,
                  out: Path) -> tuple[Path, dict]:
    """ORFs for every arm of every seed, in one FASTA, so one encoding pass covers them all.

    The shuffled null is drawn from the same sampled contigs as the real arm it controls, not
    from the whole pile, so the null is a shuffle of exactly the material searched.
    """
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    rows: list[dict] = []
    proteins: list[tuple[str, str]] = []
    blocks: dict[str, dict] = {}
    census: dict[str, dict] = {}
    failures: dict[str, int] = {}
    for pile, path in piles.items():
        blocks[pile] = {str(k): v for k, v in sample_blocks(path, seeds, per_seed).items()}
        wanted = {n for members in blocks[pile].values() for n in members}
        seqs = {n: s for n, s in read_fasta(path) if n in wanted}
        if len(seqs) != len(wanted):
            raise RuntimeError(f"{pile}: sampled {len(wanted)} contigs but read {len(seqs)}")
        for seed in seeds:
            members = blocks[pile][str(seed)]
            call_orfs(((n, seqs[n]) for n in members), pile, f"{pile}_s{seed}", rows,
                      proteins, finder)
            rng = random.Random(1_000_000 + seed)
            kept: list[tuple[str, str]] = []
            failed = 0
            for name in members:
                try:
                    arm = shuffle_arm([(name, seqs[name])], rng, kind="dinuc")
                except ShuffleFailed:
                    # Counted here, in the open, rather than inside the shared module: a
                    # walk that cannot shuffle a sequence is a repetitive sequence, which is
                    # also what throws chance alignments, so dropping them silently would
                    # bias the null in a direction nobody measured.
                    failed += 1
                    continue
                for shuffled_name, shuffled in arm["shuffled"].items():
                    kept.append((f"{shuffled_name}__shuf", shuffled))
                census.setdefault(f"{pile}_s{seed}", {"offered": 0, "shuffled": 0,
                                                      "dropped_identical": 0})
                c = census[f"{pile}_s{seed}"]
                c["offered"] += arm["census"]["offered"]
                c["shuffled"] += arm["census"]["shuffled"]
                c["dropped_identical"] += arm["census"]["dropped_identical"]
            failures[f"{pile}_s{seed}"] = failed
            call_orfs(kept, f"{pile}_null", f"{pile}_null_s{seed}", rows, proteins, finder)
    faa = out / "orfs.faa"
    write_fasta(faa, proteins)
    (out / "orfs.json").write_text(json.dumps(rows, indent=1))
    (out / "blocks.json").write_text(json.dumps(blocks, indent=1))
    by_pile = collections.Counter(r["pile"] for r in rows)
    stats = dict(orfs=len(rows), residues=sum(r["aa_len"] for r in rows),
                 by_pile=dict(by_pile), contigs_per_seed=per_seed, seeds=seeds,
                 shuffle_census=census, shuffle_failed=failures,
                 median_aa={p: statistics.median([r["aa_len"] for r in rows if r["pile"] == p])
                            for p in by_pile})
    return faa, stats


def shard_fasta(faa: Path, shards: int, work: Path) -> list[Path]:
    """Split the ORFs into independent FASTA shards, round-robin by sequence.

    ROUND-ROBIN, not consecutive blocks. Encoding cost rises faster than length - ProstT5 is
    a transformer and attention is quadratic in it - so consecutive blocks of a
    length-sorted-ish file would hand one shard every long ORF and the run would wait on it.
    Round-robin mixes the lengths, so the shards finish together.
    """
    paths = [work / f"shard_{i}.faa" for i in range(shards)]
    handles = [p.open("w") for p in paths]
    try:
        for i, (name, seq) in enumerate(read_fasta(faa)):
            handles[i % shards].write(f">{name}\n{seq}\n")
    finally:
        for h in handles:
            h.close()
    return [p for p in paths if p.stat().st_size > 0]


def build_query_dbs(foldseek: Path, faa: Path, model: Path, work: Path, threads: int,
                    shards: int, residues: int) -> tuple[list[Path], dict]:
    """ProstT5 writes Foldseek's 3Di alphabet straight from sequence: no folding, no GPU.

    `--prostt5-model` takes a directory and foldseek looks inside it for a weights file by
    name, and the name differs between releases: 9-427df8a reads `model.safetensors` or
    `model.gguf` and segfaults on the f16 gguf this account holds, while 10-941cd33 reads
    `prostt5-f16.gguf`. The pin is what makes this work.

    WHY SHARDED PROCESSES RATHER THAN MORE THREADS. The reach arm measured foldseek's own
    threading at 49.6 residues/s on four cores and 66.3 on eight - 1.34x for twice the
    cores, so the thread scaling is badly sublinear and simply asking for eight threads
    wastes most of them. Independent processes do not share that bottleneck. Each one loads
    its own 2.4 GB copy of the model, so the shard count is bounded by memory and not by
    cores, and the launcher picks an instance with room for it.
    """
    weights = [model / "prostt5-f16.gguf", model / "model.gguf", model / "model.safetensors"]
    if not any(w.exists() for w in weights):
        raise FileNotFoundError(
            f"{model} holds none of the names this foldseek build reads "
            f"({', '.join(w.name for w in weights)}); it holds "
            f"{sorted(q.name for q in model.glob('*'))[:6]}.")
    pieces = shard_fasta(faa, shards, work)
    per_shard = max(1, threads // max(1, len(pieces)))
    started = time.time()
    running = []
    for i, piece in enumerate(pieces):
        qdb = work / f"qdb_{i}"
        log = (work / f"createdb_queries_{i}.log").open("w")
        cmd = [str(foldseek), "createdb", str(piece), str(qdb), "--prostt5-model",
               str(model), "--threads", str(per_shard)]
        print("+ " + " ".join(cmd), flush=True)
        running.append((qdb, log, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)))
    qdbs, failed = [], []
    for qdb, log, proc in running:
        rc = proc.wait()
        log.close()
        if rc != 0 or not qdb.with_suffix(".index").exists():
            failed.append((qdb.name, rc))
            continue
        qdbs.append(qdb)
    if failed:
        raise RuntimeError(
            f"{len(failed)} of {len(running)} encoder shards failed: {failed}. A partial "
            f"query set would read as structure placing fewer genes rather than as a run "
            f"that encoded fewer, so this refuses. The shard logs are in {work}.")
    seconds = time.time() - started
    return qdbs, dict(seconds=round(seconds, 1), residues=residues, shards=len(qdbs),
                      threads_per_shard=per_shard,
                      residues_per_second=round(residues / seconds, 1) if seconds else None)


def search(foldseek: Path, qdbs: list[Path], sdb: Path, out: Path, threads: int,
           sensitivity: float) -> tuple[Path, dict]:
    """Prefiltered search, per shard, and the arm says so.

    The reach arm turned the prefilter off, which it could afford at 202 queries against
    178,660 chains. This arm has tens of thousands of queries against the whole databank's
    representatives, and exhaustive would be weeks. So the prefilter is on at the most
    sensitive setting the tool offers, and every placement rate here is a FLOOR: a structural
    relative distant enough to be dropped by the prefilter is exactly the one this arm wants
    and the one it may miss. That is recorded in the manifest rather than left to a reader.

    Per shard rather than over a merged database, because concatenating foldseek databases is
    a step that can half-succeed and a query silently missing from the merged set would
    lower the placement rate with nothing to show it. The shards' own row counts are summed
    and checked against the query count instead.
    """
    tsv = out / "foldseek_pdb30.tsv"
    started = time.time()
    per_shard = []
    with tsv.open("w") as merged:
        for i, qdb in enumerate(qdbs):
            resdb, tmp = out / f"fsr_{i}", out / f"fstmp_{i}"
            sh([foldseek, "search", qdb, sdb, resdb, tmp, "-a", "--threads", threads,
                "-e", EVALUE_SEARCHED, "--max-seqs", 300, "-s", sensitivity],
               capture_output=True, text=True)
            part = out / f"hits_{i}.tsv"
            sh([foldseek, "convertalis", qdb, sdb, resdb, part, "--threads", threads,
                "--format-output",
                "query,target,fident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,qcov"],
               capture_output=True, text=True)
            rows = 0
            with part.open() as fh:
                for line in fh:
                    merged.write(line)
                    rows += 1
            per_shard.append(dict(shard=i, queries=sum(1 for _ in qdb.with_suffix(".index").open()),
                                  rows=rows))
            part.unlink(missing_ok=True)
            shutil.rmtree(tmp, ignore_errors=True)
    return tsv, dict(mode=f"prefiltered -s {sensitivity}", exhaustive=False,
                     evalue_searched=EVALUE_SEARCHED, max_seqs=300, per_shard=per_shard,
                     rate_is_a_floor=True, seconds=round(time.time() - started, 1))


def annotate_subjects(targets: list[str], out: Path) -> dict:
    """What the placed thing is LIKE: the subject entry's own description and source.

    This is the half Imo asked for. A placement rate answers "is it new"; it takes the
    entry's own words to answer "and like what". Fetched in batches from RCSB's data API for
    the entries that were actually hit, so the cost is bounded by the result and not by the
    library.

    TWO AXES, NEVER ONE. `description` is what the record IS - a polymerase, a capsid - and
    `source_organism` is where that record came from, which is not what the query infects
    and not what is in the water. They are carried as separate fields for that reason.
    """
    entries = sorted({t.split("_")[0].split(".")[0].upper() for t in targets})
    rows: dict[str, dict] = {}
    for i in range(0, len(entries), 100):
        batch = entries[i:i + 100]
        ids = ", ".join(f'"{e}"' for e in batch)
        query = (f'{{entries(entry_ids:[{ids}]){{rcsb_id '
                 f'struct{{title}} '
                 f'polymer_entities{{rcsb_polymer_entity{{pdbx_description}} '
                 f'rcsb_entity_source_organism{{ncbi_scientific_name '
                 f'ncbi_parent_scientific_name}}}}}}}}')
        payload = out / f"gql_{i}.json"
        payload.write_text(json.dumps({"query": query}))
        done = subprocess.run(
            ["curl", "-sS", "--max-time", "180", "-X", "POST", RCSB_GRAPHQL,
             "-H", "Content-Type: application/json", "--data", f"@{payload}"],
            capture_output=True, text=True)
        if done.returncode != 0 or not done.stdout.strip():
            raise RuntimeError(
                f"the RCSB data API returned nothing for entries {i}+. Unannotated "
                f"placements would read as placements onto nothing, which is a different "
                f"claim. curl said: {done.stderr[-400:]}")
        body = json.loads(done.stdout)
        for entry in (body.get("data", {}).get("entries") or []):
            if not entry:
                continue
            descs, orgs = [], []
            for pe in (entry.get("polymer_entities") or []):
                d = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"))
                if d:
                    descs.append(d)
                for so in (pe.get("rcsb_entity_source_organism") or []):
                    n = so.get("ncbi_scientific_name")
                    if n:
                        orgs.append(n)
            rows[entry["rcsb_id"].upper()] = dict(
                title=(entry.get("struct") or {}).get("title"),
                descriptions=descs[:6], source_organisms=sorted(set(orgs))[:6])
        payload.unlink(missing_ok=True)
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--contigs-per-seed", type=int, required=True)
    ap.add_argument("--unidentified", type=Path, required=True)
    ap.add_argument("--identified", type=Path, required=True)
    ap.add_argument("--foldseek", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--sensitivity", type=float, default=9.5)
    ap.add_argument("--shards", type=int, default=2,
                    help="independent encoder processes; each loads its own 2.4 GB "
                         "copy of the model, so this is bounded by memory not cores")
    ap.add_argument("--skip-subjects", action="store_true",
                    help="reuse an already-built subject database in --work")
    args = ap.parse_args(argv)

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    manifest: dict = dict(provenance(Path(__file__)),
                          parameters=dict(seeds=seeds, contigs_per_seed=args.contigs_per_seed,
                                          threads=args.threads, sensitivity=args.sensitivity, shards=args.shards,
                                          min_orf_aa=MIN_ORF_AA, bands=[list(b) for b in BANDS],
                                          evalue_searched=EVALUE_SEARCHED,
                                          evalue_floor=EVALUE_FLOOR),
                          started_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))

    sdb = args.work / "subjects" / "pdb30"
    if not args.skip_subjects or not sdb.with_suffix(".index").exists():
        manifest["library"] = cluster_representatives(args.work / "subjects")
        manifest["structures"] = fetch_structures(args.work / "subjects" / "pdb30_entry_ids.txt",
                                                  args.work / "subjects" / "cif")
        manifest["subject_db"] = build_subject_db(args.foldseek, args.work / "subjects" / "cif",
                                                  sdb, args.threads)
    else:
        manifest["subject_db"] = dict(
            chains=sum(1 for _ in sdb.with_suffix(".index").open()), reused=True)

    piles = {"unident": args.unidentified, "ident": args.identified}
    faa, qstats = build_queries(piles, seeds, args.contigs_per_seed, args.out)
    manifest["queries"] = qstats
    qdbs, qmeta = build_query_dbs(args.foldseek, faa, args.model, args.work, args.threads,
                                  args.shards, qstats["residues"])
    manifest["query_db"] = qmeta
    tsv, smeta = search(args.foldseek, qdbs, sdb, args.work, args.threads, args.sensitivity)
    manifest["search"] = smeta
    searched = sum(s["queries"] for s in smeta["per_shard"])
    if searched != qstats["orfs"]:
        raise RuntimeError(
            f"{searched} ORFs reached the search but {qstats['orfs']} were called. A query "
            f"lost between the encoder and the search lowers the placement rate with "
            f"nothing to show it, so this refuses rather than writing the table.")

    hits = [line.rstrip("\n").split("\t") for line in tsv.open() if line.strip()]
    cols = ["query", "target", "fident", "alnlen", "qlen", "tlen", "qstart", "qend",
            "tstart", "tend", "evalue", "bits", "qcov"]
    rows = [dict(zip(cols, h)) for h in hits]
    manifest["rows"] = len(rows)
    with gzip.open(args.out / "hits.tsv.gz", "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for h in hits:
            fh.write("\t".join(h) + "\n")

    significant = [r for r in rows if float(r["evalue"]) <= EVALUE_FLOOR]
    annotations = annotate_subjects([r["target"] for r in significant], args.work) \
        if significant else {}
    (args.out / "subject_annotations.json").write_text(json.dumps(annotations, indent=1))
    manifest["annotated_entries"] = len(annotations)
    manifest["tmscore_available"] = False
    manifest["tmscore_reason"] = (
        "ProstT5 writes 3Di from sequence, so the query side has no C-alpha; TM-score and "
        "LDDT need coordinates on both sides")
    manifest["host"] = os.uname().nodename
    manifest["cpu_count"] = os.cpu_count()
    manifest["seconds"] = round(time.time() - started, 1)
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(json.dumps({k: manifest[k] for k in
                      ("subject_db", "queries", "query_db", "search", "rows")}, indent=1)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
