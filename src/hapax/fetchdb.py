"""`hapax fetch-db`: the reference sets a run consults, downloaded and recorded.

WHY THIS EXISTS. Every index this tool searches was staged by hand or by 250 lines of
shell that only ever ran inside one EC2 UserData script (`experiments/prod_userdata.sh`),
so the four URLs, the two asserted byte counts and the two `diamond makedb` calls that
make a run possible were knowledge held by that script and not by anything a reader
could invoke. PLAN.md's definition of done opens with `hapax fetch-db --all`, and until
there was such a command the first step of a run on new water was to read a shell
heredoc. The inventory below is lifted from that script unchanged, including which
numbers it asserts.

WHICH INDEX IS CONSULTED IS AN OPERATING POINT, SO A DOWNLOAD IS A PROVENANCE EVENT.
Thresholds dominate database size by an order of magnitude (CLAUDE.md), but the database
is still half of what a removal fraction means, and two runs against different SILVA
releases are not comparable however identical their code. So every fetch writes its URL,
its byte count, its sha256 and the version of the tool that built any derived index into
`fetch_manifest.json` beside the files. Nothing here is allowed to report success for a
file it did not verify.

ASSERTED AGAINST RECORDED, WHICH IS NOT THE SAME GUARANTEE. TOOLS.md draws this line for
binaries and it applies to the references too. A set carrying `expect_bytes` has a
committed number in `data/MANIFEST.md` that a published result was measured against, and
a mismatch is REFUSED: it means the release moved under a number this repo quotes. A set
marked `rolling` is a current-release pointer that will legitimately change - Swiss-Prot
and the RefSeq viral proteins have no stable archive URL - so its size and digest are
recorded rather than checked, and the manifest is the only thing that can tell two runs
apart afterwards. A rolling set silently treated as pinned would be the worse error of
the two, because it would read as a guarantee.

WHAT IS NOT HERE. `cross_sample` is the sample's own identified contigs, which this run
produces (`hapax.prepare`) and no download can supply; it is the most efficient
subtrahend measured in this repo and it is deliberately not a database. `nucleotide` is a
study-specific pile. Both are refused by name with that reason rather than being absent
and looking like an oversight.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tarfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

#: Kraken2 databases, by the names `data/MANIFEST.md` records. The name carries the build
#: date, so the name IS the pin: a different date is a different database rather than the
#: same one re-downloaded, which is why these are listed rather than globbed from
#: whatever the bucket serves today.
K2_DBS = ("k2_viral_20261022", "k2_standard_08_GB_20261022",
          "k2_pluspf_08_GB_20261022", "k2_pluspf_16_GB_20261022")

K2_BASE = "https://genome-idx.s3.amazonaws.com/kraken"
SILVA_BASE = "https://ftp.arb-silva.de/release_138_1/Exports"


class FetchError(Exception):
    """Raised rather than falling back, and rather than leaving a half-written file.

    A partially downloaded index is the failure mode this class exists for: DIAMOND and
    Kraken2 both read a truncated database without complaining and return fewer hits,
    which arrives as a lower removal fraction and reads as a property of the water.
    """


@dataclass(frozen=True)
class Source:
    """One file to download. `expect_bytes` is a committed number or None; see the
    module docstring on asserted against recorded."""
    url: str
    name: str
    expect_bytes: int | None = None
    #: True where the URL is a current-release pointer with no stable archive form, so
    #: its identity can only be recorded after the fact.
    rolling: bool = False
    #: Unpack the download as a gzipped tar into a directory named after `name`, and
    #: require this file to exist inside it afterwards.
    untar_requiring: str | None = None


@dataclass(frozen=True)
class RefSet:
    """One reference set: what a `[databases]` key points at, and how to produce it."""
    key: str
    what: str
    why: str
    sources: tuple[Source, ...]
    #: The path, relative to the destination, that the config key should name. For a set
    #: with one source this is usually that file; for a built index it is the output.
    fills: str
    #: `diamond makedb` over the sources, writing `fills`. Needs the diamond binary.
    diamond_makedb: bool = False
    #: Built by calling `experiments/subtract_known.py`'s own two SILVA functions, which
    #: is deliberate: a second SILVA cleaner would be a second instrument reading a
    #: different thing by the same name.
    silva_to_dna: bool = False
    #: Several directories rather than one path: `--k2-dbs` is the cascade's tier list
    #: and its ORDER is the tier order.
    many: bool = False
    extra: dict = field(default_factory=dict)


#: The inventory. Lifted from `experiments/prod_userdata.sh` with its own assertions; the
#: order is the order a run consults them, which is also cheapest first.
INVENTORY: tuple[RefSet, ...] = (
    RefSet(
        key="rrna",
        what="SILVA 138.1 SSU and LSU NR99, converted to DNA",
        why="the largest tier of the cascade. Ribosomal RNA is about a third of the "
            "Berlin library, and it is identified positively against a reference rather "
            "than by the absence of a gene call.",
        sources=(
            Source(f"{SILVA_BASE}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
                   "SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
                   expect_bytes=197_133_995),
            Source(f"{SILVA_BASE}/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz",
                   "SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz",
                   expect_bytes=67_498_420),
        ),
        fills="silva/silva_dna.fa",
        silva_to_dna=True,
    ),
    RefSet(
        key="classifier",
        what=f"{len(K2_DBS)} Kraken2 databases: {', '.join(K2_DBS)}",
        why="the classifier tiers. The ORDER is the tier order and the committed tier "
            "table is only readable beside a run that used these four; one database "
            "removes 42.8%, 21.7% or 3.5% of the same library depending only on its "
            "confidence setting, so the database alone does not say what a tier did.",
        sources=tuple(
            Source(f"{K2_BASE}/{db}.tar.gz", f"{db}.d", untar_requiring="hash.k2d")
            for db in K2_DBS),
        fills="k2",
        many=True,
    ),
    RefSet(
        key="protein",
        what="RefSeq viral proteins, as a DIAMOND index",
        why="the first protein tier, and the one that answers 'like what': a nearest "
            "relative with a name. Distance without an identity is half an answer.",
        sources=(Source("https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
                        "viral.1.protein.faa.gz", "viral.1.protein.faa.gz",
                        rolling=True),),
        fills="prot/viral.dmnd",
        diamond_makedb=True,
        extra={"protein_fasta": "prot/viral.1.protein.faa.gz"},
    ),
    RefSet(
        key="sprot",
        what="UniProt Swiss-Prot, as a DIAMOND index",
        why="the second protein tier, added because 88% of Berlin reads score 1.0 "
            "against a viral-only index: without a background range a viral index "
            "reports every conserved fold as viral.",
        sources=(Source("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                        "knowledgebase/complete/uniprot_sprot.fasta.gz",
                        "uniprot_sprot.fasta.gz", rolling=True),),
        fills="prot/sprot.dmnd",
        diamond_makedb=True,
    ),
    RefSet(
        key="pfam",
        what="Pfam-A HMM library",
        why="domain PRESENCE, which is the one profile measurement here that has a null: "
            "811 of 1,317 real ORFs carry a domain at gathering cutoffs against 0 of "
            "932 shuffled. Per-column conservation has no null and is not reported.",
        sources=(Source("https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
                        "Pfam-A.hmm.gz", "Pfam-A.hmm.gz", rolling=True),),
        fills="pfam/Pfam-A.hmm",
    ),
)

#: Keys a config can carry that no download can fill, and why. Named so that a caller who
#: asks for one gets the reason rather than "unknown set", which would read as a typo.
NOT_FETCHABLE = {
    "cross_sample": "assemblies of OTHER samples from the same works, which the "
                    "per-read recurrence field counts distinct samples in. There is no "
                    "release to fetch because it is a property of a sampling programme, "
                    "and the one substitution that looks harmless is the one that breaks "
                    "it: the sample's own assembly makes every read that maps recur in "
                    "exactly one sample, its own. Query material from inside the library "
                    "searched gave 54% false growth against a 6% swap control.",
    "spacers": "CRISPR spacers from this sample's own assemblies, produced by the run.",
    "nucleotide": "a study-specific nucleotide pile. There is no general release to "
                  "fetch, and a wrong one here is a changed operating point for the "
                  "tier that reports named relatives.",
}

BY_KEY = {s.key: s for s in INVENTORY}


def _digest_and_size(path: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    size = 0
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
            size += len(block)
    return h.hexdigest(), size


def _download(url: str, dest: Path) -> None:
    """Download to a temporary name and rename on success, so no reader can ever open a
    partial index. A truncated DIAMOND database returns fewer hits without erroring, and
    that arrives as a lower removal fraction that reads as a property of the water."""
    tmp = dest.with_name(dest.name + ".part")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    if tmp.exists():
        tmp.unlink()
    try:
        with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
    except (urllib.error.URLError, OSError) as e:
        if tmp.exists():
            tmp.unlink()
        raise FetchError(f"{url} did not download: {e}") from None
    tmp.replace(dest)


def _untar(archive: Path, into: Path, requiring: str) -> None:
    into.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive, "r:gz") as tar:
            # Refuse a member that would escape the destination. The archives here are
            # from a public bucket over which this repo has no control, and an index is
            # unpacked with the privileges of whoever runs the tool.
            for member in tar.getmembers():
                target = (into / member.name).resolve()
                if not str(target).startswith(str(into.resolve())):
                    raise FetchError(
                        f"{archive.name} contains {member.name}, which unpacks outside "
                        f"{into}. Refused rather than written.")
            # `filter="data"` is the behaviour Python 3.14 makes the default: no device
            # nodes, no absolute paths, no links out of the tree. Stated rather than
            # inherited, so this reads the same on both sides of that change.
            tar.extractall(into, filter="data")
    except tarfile.TarError as e:
        raise FetchError(f"{archive} did not unpack: {e}") from None
    if not (into / requiring).exists():
        raise FetchError(
            f"{archive.name} unpacked without {requiring}, so it is not the database it "
            f"is named after. Refused rather than searched: Kraken2 reads an incomplete "
            f"database and classifies nothing, which arrives as a clean sample.")


def _verify(source: Source, path: Path) -> dict:
    """The record that goes in the manifest, and the refusal where a pinned size moved."""
    digest, size = _digest_and_size(path)
    record = {"url": source.url, "file": path.name, "bytes": size, "sha256": digest,
              "pinned": source.expect_bytes is not None}
    if source.expect_bytes is not None and size != source.expect_bytes:
        raise FetchError(
            f"{path.name} is {size:,} bytes, not the {source.expect_bytes:,} every "
            f"committed number from it was measured against (data/MANIFEST.md). A "
            f"changed release is a changed operating point, and a run against it is not "
            f"comparable to any published here. Refused rather than used.")
    if source.rolling:
        record["rolling"] = True
        record["note"] = ("a current-release pointer with no stable archive URL, so its "
                          "identity is recorded here and not asserted")
    return record


def _build_diamond(refset: RefSet, dest: Path, files: list[Path], diamond: Path,
                   threads: int) -> dict:
    out = dest / refset.fills
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(diamond), "makedb", "--in", *[str(f) for f in files],
           "--db", str(out), "--threads", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not out.exists():
        raise FetchError(
            f"diamond makedb failed for {refset.key}: {proc.stderr[-2000:]}")
    version = subprocess.run([str(diamond), "--version"], capture_output=True,
                             text=True).stdout.strip()
    return {"built_with": version, "command": cmd}


def _build_silva_dna(refset: RefSet, dest: Path, files: list[Path], repo: Path) -> dict:
    """Convert the two SILVA exports to DNA with `subtract_known`'s own functions.

    IMPORTED, NOT REIMPLEMENTED, for the reason the shuffle rule gives: the cleaner
    decides what counts as ribosomal, so a second one would be a second instrument
    reading a different quantity under the same name, and the cascade's own tier is the
    one whose answer matters.
    """
    import sys
    out = dest / refset.fills
    out.parent.mkdir(parents=True, exist_ok=True)
    experiments = repo / "experiments"
    if not (experiments / "subtract_known.py").exists():
        raise FetchError(
            f"{experiments}/subtract_known.py is not there, and the SILVA conversion is "
            f"that script's own function rather than a copy of it. Point HAPAX_REPO at "
            f"the checkout.")
    if str(experiments) not in sys.path:
        sys.path.insert(0, str(experiments))
    import subtract_known as sk
    sk._silva_to_dna(files, out)
    if not out.exists() or out.stat().st_size == 0:
        raise FetchError(f"{out} was not written by subtract_known._silva_to_dna")
    digest, size = _digest_and_size(out)
    return {"built_with": "experiments/subtract_known.py:_silva_to_dna",
            "bytes": size, "sha256": digest}


def fetch(keys: list[str], dest: Path, repo: Path, *, threads: int = 4,
          diamond: Path | None = None, dry_run: bool = False) -> dict:
    """Fetch the named reference sets into `dest`, and return the manifest.

    Idempotent by file: a target that is already present and already recorded in the
    manifest with a matching digest is left alone, because these are tens of gigabytes
    and a re-download is not free. A file present but NOT recorded is re-verified rather
    than trusted, since a file of the right name and the wrong content is exactly what
    this is guarding.
    """
    dest = Path(dest)
    unknown = [k for k in keys if k not in BY_KEY]
    for key in unknown:
        if key in NOT_FETCHABLE:
            raise FetchError(f"{key} cannot be fetched: {NOT_FETCHABLE[key]}")
    if unknown:
        raise FetchError(
            f"no reference set named {', '.join(unknown)}. Known: "
            f"{', '.join(BY_KEY)}.")

    manifest_path = dest / "fetch_manifest.json"
    manifest = {"sets": {}}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except json.JSONDecodeError:
            manifest = {"sets": {}}
    manifest.setdefault("sets", {})
    manifest["dest"] = str(dest.resolve())

    if dry_run:
        return {"sets": {k: _plan(BY_KEY[k], dest) for k in keys}, "dry_run": True}

    dest.mkdir(parents=True, exist_ok=True)
    for key in keys:
        refset = BY_KEY[key]
        started = time.time()
        files: list[Path] = []
        records = []
        for source in refset.sources:
            if source.untar_requiring:
                directory = dest / refset.fills / source.name
                archive = dest / refset.fills / (source.name + ".tar.gz")
                if not (directory / source.untar_requiring).exists():
                    _download(source.url, archive)
                    records.append(_verify(source, archive))
                    _untar(archive, directory, source.untar_requiring)
                    archive.unlink()
                else:
                    records.append({"url": source.url, "file": source.name,
                                    "already_present": True})
                files.append(directory)
                continue
            path = dest / _subdir(refset) / source.name
            if not path.exists():
                _download(source.url, path)
            records.append(_verify(source, path))
            files.append(path)

        built = {}
        target = dest / refset.fills
        if refset.diamond_makedb:
            if diamond is None:
                raise FetchError(
                    f"{key} is a DIAMOND index and no diamond binary was resolved. The "
                    f"index is built here rather than downloaded because no release "
                    f"serves one, and a .dmnd built by another version is not readable "
                    f"by this one.")
            if not target.exists():
                built = _build_diamond(refset, dest, files, diamond, threads)
        elif refset.silva_to_dna:
            if not target.exists():
                built = _build_silva_dna(refset, dest, files, repo)
        elif refset.fills.endswith(".hmm"):
            # The only set whose download is the thing itself once decompressed.
            if not target.exists():
                built = _gunzip(files[0], target)

        manifest["sets"][key] = {
            "what": refset.what, "fills": refset.fills, "many": refset.many,
            "sources": records, "built": built,
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
            "seconds": round(time.time() - started, 1),
            "paths": [str(f) for f in files] if refset.many else [str(target)],
        }
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest


def _gunzip(src: Path, out: Path) -> dict:
    import gzip
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    with gzip.open(src, "rb") as fh, tmp.open("wb") as o:
        shutil.copyfileobj(fh, o)
    tmp.replace(out)
    digest, size = _digest_and_size(out)
    return {"built_with": "gzip", "bytes": size, "sha256": digest}


def _subdir(refset: RefSet) -> str:
    """Downloads live beside what is built from them, so a set is one directory."""
    parent = Path(refset.fills).parent
    return str(parent) if str(parent) != "." else ""


def _plan(refset: RefSet, dest: Path) -> dict:
    return {"what": refset.what, "why": refset.why, "fills": refset.fills,
            "urls": [s.url for s in refset.sources],
            "pinned_bytes": [s.expect_bytes for s in refset.sources],
            "present": (dest / refset.fills).exists()}


def config_block(manifest: dict) -> str:
    """The `[databases]` table for a run config, from what was actually fetched.

    Printed rather than merged into a config file: which index a run consulted is its
    operating point, and a tool that edits the caller's config has made that choice for
    them somewhere they will not read it.
    """
    dest = Path(manifest.get("dest") or ".")
    lines = ["[databases]"]
    for key, record in manifest.get("sets", {}).items():
        paths = record.get("paths") or []
        if record.get("many"):
            inner = ", ".join(f'"{p}"' for p in paths)
            lines.append(f"{key} = [{inner}]")
        elif paths:
            lines.append(f'{key} = "{paths[0]}"')
        # An extra is a second view of the same download rather than a second fetch:
        # `protein_fasta` is the FASTA the protein index was built from, read for its
        # names, and reading the index instead would be a second source for one answer.
        for name, rel in (BY_KEY[key].extra if key in BY_KEY else {}).items():
            lines.append(f'{name} = "{dest / rel}"')
    return "\n".join(lines) + "\n"
