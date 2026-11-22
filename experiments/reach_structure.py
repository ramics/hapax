#!/usr/bin/env python3
"""The structure arm of the reach comparison, against solved virus structures.

WHY THIS IS A SEPARATE SCRIPT FROM `reach_arms.py`. The other eight arms need the subject
sets that only fit on an instance - Swiss-Prot plus RefSeq viral as BLAST volumes and a
DIAMOND index, Pfam-A pressed, RefSeq viral genomes. This arm needs about 10 GB of
coordinates and a 2.3 GB model, which fits in the container, so it runs here for nothing
instead of adding an hour of staging to a paid run. It calls `reach_arms.call_orfs`, so the
genes are the same genes, and the manifest records the hash of `orfs.faa` to prove it rather
than asserting it.

WHY NOT AlphaFold/Swiss-Prot, WHICH WAS THE PLAN. Two reasons, one practical and one that
matters.

The practical one: every prebuilt Foldseek database is served by one Cloudflare worker, that
worker answered HTTP 429 with a 17-byte body from three networks across several hours on
2026-11-22, and `foldseek databases` untars the error body and exits 0 anyway, so three runs
reported a structure arm that had never had subjects. The only bulk-servable alternative,
EBI's 28.6 GB `swissprot_pdb_v6.tar`, came down at 6.7 MB/s and then stalled outright at
33 GB of a 200 GB disk, taking a run's whole dead-man timer with it.

The one that matters: AlphaFold/Swiss-Prot was chosen because it is Swiss-Prot folded, which
would make this arm the sequence arms' own subjects in a different alphabet. But the shared
sequence database is Swiss-Prot PLUS the RefSeq viral proteins, and `report_reach.py` §5b
counts where the nearest relatives actually come from: 47 of the 51 ORFs with a significant
baseline hit find it in the RefSeq viral half, which is folded nowhere. A null from that
target would have been a statement about Swiss-Prot's coverage of phage proteins and not
about structure, and it would have read as the latter.

WHAT THIS TARGET IS, AND WHAT IT IS NOT. Every experimentally determined PDB entry whose
source organism lineage includes Viruses (NCBI taxid 10239), fetched by the RCSB search API
so the selection is a query and not a list somebody typed. That is the universe of solved
virus structures: chosen by what has been crystallised, which is independent of anything our
sequence arms found, and it is where a remote structural relative of a phage capsid,
maturation protein or replicase would actually sit.

It is NOT an unbiased sample of fold space, and it is smaller than the sequence arms'
subject set by two orders of magnitude. So this arm can show that structural comparison does
or does not place these genes at distances where sequence fails, which is a calibration of
the instrument. It cannot be read as "structure reaches further than sequence over a
comparable universe", because the universes are not comparable and the report says so.

THE QUERY SIDE HAS NO COORDINATES, which is the whole point of running it here. ProstT5
writes Foldseek's 3Di alphabet straight from sequence, so the arm costs CPU rather than a
GPU and a folding run. Foldseek then aligns 3Di plus amino acid. No TM-score and no LDDT are
reportable: both are superposition scores and need C-alpha on both sides, and foldseek says
so itself (`Cannot find <query> C-alpha or <target> C-alpha database`). Folding the 54 real
ORFs is the way to get those two numbers and is about $4 on a GPU instance; folding the
SUBJECT side - the only way to make this a controlled comparison against the sequence arms -
is 175,064,339 residues for the RefSeq viral half alone, which at this arm's own measured
66.3 residues/s on 8 vCPU is 30 days of CPU and therefore a GPU job of its own.

CONTROLS, in the same pass and through the same search, because this project's rule is that
the control is read first: the seven near-reference contigs from the same assemblies as the
negative arm, and five seeds of dinucleotide-preserving shuffles of every panel member. The
shuffle comes from `reach_queries.py`, which takes it from the one shared module, and the
null piles are the same files the paid arms searched.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: NCBI taxid for Viruses. The target set is defined by this query against the RCSB search
#: API rather than by an accession list, so the selection is reproducible and its date is
#: recorded: the PDB grows, and a run six months from now should get a bigger set and say so.
VIRUS_TAXID = "10239"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download/{}.cif.gz"

#: Entries per page of the search API, which caps a page at 10,000.
PAGE_ROWS = 10000

#: Concurrent fetches. Deliberately modest: this is a public archive and the measured rate at
#: eight is already about 13 entries a second, so the whole set is well under an hour.
FETCH_PARALLEL = 8

#: The E-value floor every arm in this comparison reports against, from `report_reach.py`.
#: Nothing here is thresholded on it - the rows are written un-thresholded and the readout
#: applies it - but it is recorded so the manifest says what the run was aimed at.
EVALUE_REPORTED = 10.0


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    """Run a command, failing loudly, with its own output kept.

    `check=True` alone discards the tool's complaint and the exit code carries none of it,
    which is how a foldseek segfault spent four runs looking like a bare exit 1.
    """
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def load_call_orfs():
    """Borrow the gene caller from `reach_arms.py` rather than reimplementing it.

    THE GENES MUST BE THE SAME GENES. A structure arm on its own ORF set would not be
    comparable with the sequence arms at all, and the failure would be invisible: both
    tables would be full and the contigs would have the same names. So this imports the one
    implementation by path, and the manifest carries the sha256 of the `orfs.faa` it
    produces, which is the observable version of the claim.
    """
    spec = importlib.util.spec_from_file_location(
        "reach_arms", ROOT / "experiments" / "reach_arms.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.call_orfs


def fetch_entry_ids(out: Path) -> dict:
    """Every experimental PDB entry whose source lineage includes Viruses.

    Paged, de-duplicated and sorted, so two runs on the same day produce the same file. The
    returned dict is provenance: the query, the date, and the count the server reported,
    against which the number of files actually fetched can be checked.
    """
    ids: list[str] = []
    total = None
    for start in range(0, 1_000_000, PAGE_ROWS):
        query = {
            "query": {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                "operator": "exact_match", "value": VIRUS_TAXID}},
            "return_type": "entry",
            "request_options": {"paginate": {"start": start, "rows": PAGE_ROWS},
                                "results_content_type": ["experimental"]}}
        page = out / f"rcsb_query_{start}.json"
        page.write_text(json.dumps(query))
        done = subprocess.run(
            ["curl", "-sS", "--max-time", "180", "-G", RCSB_SEARCH,
             "--data-urlencode", f"json@{page}"], capture_output=True, text=True)
        if done.returncode != 0 or not done.stdout.strip():
            raise RuntimeError(
                f"the RCSB search API returned nothing for rows {start}+. An empty target "
                f"set would read as structure finding no relative, which is a different "
                f"claim. curl said: {done.stderr[-400:]}")
        body = json.loads(done.stdout)
        total = body["total_count"]
        got = [r["identifier"] for r in body.get("result_set", [])]
        ids += got
        print(f"  rcsb rows {start}+: {len(got)} of {total}", flush=True)
        if not got or len(ids) >= total:
            break
    ids = sorted(set(ids))
    (out / "pdb_ids.txt").write_text("\n".join(ids) + "\n")
    return dict(taxid=VIRUS_TAXID, api=RCSB_SEARCH, reported_total=total,
                unique_entries=len(ids),
                fetched_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))


def fetch_structures(ids_file: Path, cif_dir: Path) -> dict:
    """Download each entry's mmCIF, skipping anything already on disk.

    Resumable on purpose: 27,000 files is long enough that a run should not start over, and
    the skip is on a non-empty file so a truncated download is retried rather than kept.
    """
    cif_dir.mkdir(parents=True, exist_ok=True)
    script = (f'test -s {cif_dir}/$1.cif.gz || '
              f'curl -sS --max-time 90 --retry 2 -o {cif_dir}/$1.cif.gz '
              f'"{RCSB_FILE.format("$1")}"')
    started = time.time()
    subprocess.run(["xargs", "-P", str(FETCH_PARALLEL), "-I", "{}",
                    "sh", "-c", script.replace("$1", "{}"), "{}"],
                   stdin=ids_file.open(), check=False)
    present = sorted(p for p in cif_dir.glob("*.cif.gz") if p.stat().st_size > 0)
    wanted = [i for i in ids_file.read_text().split() if i]
    missing = [i for i in wanted if not (cif_dir / f"{i}.cif.gz").exists()]
    return dict(wanted=len(wanted), present=len(present), missing=len(missing),
                missing_examples=missing[:10], seconds=round(time.time() - started, 1))


def build_subject_db(foldseek: Path, cif_dir: Path, db: Path, threads: int) -> dict:
    """Foldseek database from real coordinates, which is what makes it a structure target.

    `createdb` over mmCIF reads C-alpha and derives 3Di from geometry. The subject side
    therefore has `_ca` and the ProstT5 query side does not, which is exactly why no
    superposition score is reportable and the arm says so rather than omitting the columns.
    """
    db.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    done = sh([foldseek, "createdb", cif_dir, db, "--threads", threads],
              capture_output=True, text=True)
    log = (done.stdout or "") + (done.stderr or "")
    (db.parent / "createdb_subjects.log").write_text(log)
    index = db.with_suffix(".index")
    if not index.exists():
        raise RuntimeError(
            "foldseek createdb wrote no index over the structures, so this arm has no "
            "subjects. Its own output:\n" + log[-2500:])
    subjects = sum(1 for _ in index.open())
    # A floor, not an expectation. The query returns about 27,000 entries and each yields
    # one or more chains, so anything in the hundreds means the fetch or the parse failed
    # and the table would read as structure finding fewer relatives rather than as a short
    # target set.
    if subjects < 5000:
        raise RuntimeError(
            f"the subject database holds {subjects} chains, which is an order of magnitude "
            f"short of what this query returns. That is a staging failure, and an empty "
            f"table here would read as a result.")
    return dict(chains=subjects, seconds=round(time.time() - started, 1))


def build_query_db(foldseek: Path, orfs: Path, qdb: Path, model: Path,
                   threads: int) -> dict:
    """ProstT5 writes the 3Di alphabet straight from sequence.

    `--prostt5-model` takes a DIRECTORY and foldseek looks inside it for a weights file by
    name, and the name differs between releases: 9-427df8a reads `model.safetensors` or
    `model.gguf` and SEGFAULTS on the f16 gguf this account holds, while 10-941cd33 reads
    `prostt5-f16.gguf` and reads that same file. The names are all accepted here so the arm
    does not track the pin, but the pin is what makes it work.
    """
    weights = [model / "prostt5-f16.gguf", model / "model.gguf",
               model / "model.safetensors", model / "model" / "model.gguf"]
    if not any(w.exists() for w in weights):
        raise FileNotFoundError(
            f"{model} holds none of the names this foldseek build reads "
            f"({', '.join(w.name for w in weights)}). Its contents are "
            f"{sorted(q.name for q in model.glob('*'))[:6]}.")
    started = time.time()
    done = sh([foldseek, "createdb", orfs, qdb, "--prostt5-model", model,
               "--threads", threads], capture_output=True, text=True)
    (qdb.parent / "createdb_queries.log").write_text((done.stdout or "") + (done.stderr or ""))
    if not qdb.with_suffix(".index").exists():
        raise RuntimeError("foldseek createdb wrote no query index. Its own output:\n"
                           + ((done.stdout or "") + (done.stderr or ""))[-2500:])
    return dict(seconds=round(time.time() - started, 1))


def search(foldseek: Path, qdb: Path, sdb: Path, out: Path, threads: int) -> tuple[Path, str]:
    """Search with the prefilter off, which is the point of the arm.

    A structural relative this far out is exactly what a k-mer prefilter drops, and the query
    set is 202 ORFs, so the exhaustive mode is affordable. If a build will not take the flag
    the fallback says what it fell back FROM, because silently dropping to the prefilter
    changes what the arm measures and would leave no record of the swap.
    """
    resdb, tmp = out / "fsr", out / "fstmp"
    base = [foldseek, "search", qdb, sdb, resdb, tmp, "-a", "--threads", threads,
            "-e", EVALUE_REPORTED, "--max-seqs", 1000]
    try:
        sh(base + ["--exhaustive-search", 1], capture_output=True, text=True)
        mode = "exhaustive-search"
    except subprocess.CalledProcessError as exc:
        print("  foldseek --exhaustive-search refused; falling back to the most sensitive "
              "prefilter. Its own output:\n"
              + ((exc.stderr or "") + (exc.stdout or ""))[-1200:], flush=True)
        shutil.rmtree(tmp, ignore_errors=True)
        sh(base + ["-s", 9.5], capture_output=True, text=True)
        mode = "prefiltered -s 9.5"
    tsv = out / "foldseek_pdbvirus.tsv"
    sh([foldseek, "convertalis", qdb, sdb, resdb, tsv, "--threads", threads,
        "--format-output",
        "query,target,fident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,qcov"])
    return tsv, mode


def rows_from_tsv(tsv: Path, mode: str) -> list[dict]:
    """One row per alignment, un-thresholded, in the shape the other arms write.

    `tmscore` and `lddt` are carried as None rather than left out, so the readout can tell
    "this arm cannot report a superposition score" from "this arm was not run".
    """
    rows = []
    for line in tsv.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        rows.append(dict(
            arm="foldseek_pdbvirus",
            program=f"foldseek search (ProstT5 3Di, {mode})",
            db="pdbvirus", orf=f[0], subject=f[1], pident=100.0 * float(f[2]),
            aln=int(f[3]), qlen=int(f[4]), slen=int(f[5]),
            qstart=int(f[6]), qend=int(f[7]), sstart=int(f[8]), send=int(f[9]),
            evalue=float(f[10]), bitscore=float(f[11]), qcov=float(f[12]),
            tmscore=None, lddt=None))
    return rows


def git_sha() -> str:
    try:
        done = subprocess.run(["git", "-C", str(ROOT), "describe", "--always", "--dirty"],
                              capture_output=True, text=True, check=True)
        return done.stdout.strip()
    except Exception:
        return "unknown"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1,
                    help="the panel's seed. Carried into the manifest; the panel files "
                         "themselves are built by reach_queries.py.")
    ap.add_argument("--panel-dir", type=Path, required=True,
                    help="holds panel.fa and null_s*.fa from reach_queries.py")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cif-dir", type=Path, required=True,
                    help="where the mmCIF files live; reused between runs")
    ap.add_argument("--prostt5", type=Path, required=True,
                    help="directory holding the ProstT5 weights")
    ap.add_argument("--foldseek", type=Path, default=Path("foldseek"))
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--skip-fetch", action="store_true",
                    help="use the entry list and structures already on disk. The manifest "
                         "records that the set was not re-queried.")
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()

    call_orfs = load_call_orfs()
    nulls = sorted(args.panel_dir.glob("null_s*.fa"))
    orf_stats = call_orfs(args.panel_dir / "panel.fa", nulls, args.out)

    ids_file = args.out / "pdb_ids.txt"
    if args.skip_fetch and ids_file.exists():
        target = dict(skipped=True, unique_entries=len(ids_file.read_text().split()))
    else:
        target = fetch_entry_ids(args.out)
    fetched = fetch_structures(ids_file, args.cif_dir)

    subjects = build_subject_db(args.foldseek, args.cif_dir, args.out / "fs_pdbvirus",
                               args.threads)
    queries = build_query_db(args.foldseek, args.out / "orfs.faa", args.out / "fsq",
                             args.prostt5, args.threads)
    tsv, mode = search(args.foldseek, args.out / "fsq", args.out / "fs_pdbvirus",
                       args.out, args.threads)
    rows = rows_from_tsv(tsv, mode)
    with (args.out / "rows_foldseek_pdbvirus.jsonl").open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")

    version = subprocess.run([str(args.foldseek), "version"], capture_output=True,
                             text=True).stdout.strip()
    manifest = dict(
        seed=args.seed, git_sha=git_sha(), script_sha256=sha256(Path(__file__)),
        orfs_sha256=sha256(args.out / "orfs.faa"), orf_stats=orf_stats,
        foldseek_version=version, mode=mode, threads=args.threads,
        target=target, structures=fetched, subject_db=subjects, query_db=queries,
        rows=len(rows), evalue_searched=EVALUE_REPORTED,
        tmscore_available=False,
        tmscore_reason="ProstT5 writes 3Di from sequence, so the query database has no "
                       "C-alpha; TM-score and LDDT need coordinates on both sides",
        host=os.uname().nodename, cpu_count=os.cpu_count(),
        seconds=round(time.time() - started, 1))
    (args.out / "structure_manifest.json").write_text(json.dumps(manifest, indent=1))
    print(json.dumps(manifest, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
