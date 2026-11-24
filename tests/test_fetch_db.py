"""`hapax fetch-db`, tested without the network.

WHY THESE TESTS AND NOT A DOWNLOAD. The thing that can go wrong here is not whether
`urlopen` works. It is whether a reference set that arrives wrong is USED: a release that
moved under a byte count this repo quotes, a Kraken2 tarball that unpacks without its
hash table, a half-written file left where a reader will open it. Each of those arrives
silently as a lower removal fraction, which reads as a property of the water rather than
as a broken index. So every test here drives the real fetcher over `file://` URLs served
out of a temp directory, and what it asserts is the refusals.

The inventory itself is also checked, because it is the part a future edit will touch:
every set must say what fills it, every pinned size must be a number this repo has
committed, and a key that cannot be fetched must be refused with its reason rather than
with "unknown", which would read as a typo.
"""
from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import pytest

from hapax import fetchdb
from hapax.config import ConfigError, load


def _as_url(path: Path) -> str:
    return path.resolve().as_uri()


def test_every_set_declares_what_fills_it():
    for refset in fetchdb.INVENTORY:
        assert refset.sources, f"{refset.key} has no source"
        assert refset.fills, f"{refset.key} names no path for a config to point at"
        assert refset.why.strip(), f"{refset.key} has no reason to exist"
        # Exactly one way of producing the target, or none. Two would mean the file
        # written depends on evaluation order of the branches in `fetch`.
        built = [refset.diamond_makedb, refset.silva_to_dna]
        assert sum(bool(b) for b in built) <= 1, refset.key


def test_a_pinned_size_is_one_this_repo_committed():
    """A number asserted here must be findable in data/MANIFEST.md.

    An assertion against a number nobody wrote down would refuse a correct download for
    a reason no reader could check, which is worse than recording the size.
    """
    manifest = (Path(__file__).resolve().parents[1] / "data" / "MANIFEST.md").read_text()
    for refset in fetchdb.INVENTORY:
        for source in refset.sources:
            if source.expect_bytes is None:
                continue
            assert f"{source.expect_bytes:,}" in manifest, (
                f"{refset.key} asserts {source.expect_bytes:,} bytes for {source.name} "
                f"and data/MANIFEST.md does not record that number")


def test_a_set_is_either_pinned_or_marked_rolling():
    """Neither asserted nor recorded-as-rolling is the one state that lies.

    TOOLS.md draws this line for binaries. A current-release URL whose size is checked
    would refuse every legitimate update; one whose size is neither checked nor flagged
    would read as a guarantee it cannot give.
    """
    for refset in fetchdb.INVENTORY:
        for source in refset.sources:
            if source.untar_requiring:
                continue    # the name carries the build date, so the name is the pin
            assert (source.expect_bytes is not None) != source.rolling, (
                f"{refset.key}/{source.name} is neither pinned nor marked rolling")


def test_a_key_that_cannot_be_fetched_is_refused_with_its_reason(tmp_path):
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["cross_sample"], tmp_path, Path("."))
    # The reason has to be the one that matters: not "we have no copy" but "the obvious
    # substitute fabricates the measurement".
    assert "OTHER samples" in str(e.value) and "false growth" in str(e.value)
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["nonsense"], tmp_path, Path("."))
    assert "no reference set named" in str(e.value)


def test_dry_run_downloads_nothing_and_names_what_is_asserted(tmp_path):
    plan = fetchdb.fetch(list(fetchdb.BY_KEY), tmp_path, Path("."), dry_run=True)
    assert plan["dry_run"] is True
    assert set(plan["sets"]) == set(fetchdb.BY_KEY)
    assert plan["sets"]["rrna"]["pinned_bytes"] == [197_133_995, 67_498_420]
    assert plan["sets"]["sprot"]["pinned_bytes"] == [None]
    assert not list(tmp_path.iterdir())


def _local_set(served: Path, key: str, payload: bytes, *, expect: int | None,
               rolling: bool, fills: str) -> fetchdb.RefSet:
    """One real set whose source is a file on disk, so `fetch` runs unchanged."""
    served.parent.mkdir(parents=True, exist_ok=True)
    served.write_bytes(payload)
    return fetchdb.RefSet(
        key=key, what="a local stand-in", why="tests the refusals, not the network",
        sources=(fetchdb.Source(_as_url(served), served.name, expect_bytes=expect,
                                rolling=rolling),),
        fills=fills)


def test_a_moved_release_is_refused_rather_than_used(tmp_path, monkeypatch):
    """The whole point of a pinned byte count: refuse, and say what it would cost."""
    payload = b"ACGT" * 10
    refset = _local_set(tmp_path / "served" / "ssu.fa.gz", "rrna", payload,
                        expect=len(payload) + 1, rolling=False, fills="silva/ssu.fa.gz")
    monkeypatch.setitem(fetchdb.BY_KEY, "rrna", refset)
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["rrna"], tmp_path / "db", Path("."))
    assert "changed operating point" in str(e.value)
    assert f"{len(payload):,}" in str(e.value)


def test_a_rolling_set_records_its_identity_instead_of_asserting_it(tmp_path,
                                                                   monkeypatch):
    payload = b">p\nMKV\n"
    refset = _local_set(tmp_path / "served" / "prot.faa.gz", "sprot", payload,
                        expect=None, rolling=True, fills="prot/prot.faa.gz")
    monkeypatch.setitem(fetchdb.BY_KEY, "sprot", refset)
    manifest = fetchdb.fetch(["sprot"], tmp_path / "db", Path("."))
    record = manifest["sets"]["sprot"]["sources"][0]
    assert record["bytes"] == len(payload)
    assert record["rolling"] is True and record["pinned"] is False
    assert len(record["sha256"]) == 64
    # The manifest is on disk, not only returned: a run directory that cannot say which
    # index it searched has no operating point.
    on_disk = json.loads((tmp_path / "db" / "fetch_manifest.json").read_text())
    assert on_disk["sets"]["sprot"]["sources"][0]["sha256"] == record["sha256"]


def test_a_failed_download_leaves_nothing_a_reader_could_open(tmp_path, monkeypatch):
    """No partial file, because DIAMOND and Kraken2 read a truncated database quietly."""
    refset = fetchdb.RefSet(
        key="sprot", what="a URL that is not there", why="tests the partial-file path",
        sources=(fetchdb.Source(_as_url(tmp_path / "absent.gz"), "absent.gz",
                                rolling=True),),
        fills="prot/absent.gz")
    monkeypatch.setitem(fetchdb.BY_KEY, "sprot", refset)
    with pytest.raises(fetchdb.FetchError):
        fetchdb.fetch(["sprot"], tmp_path / "db", Path("."))
    assert not list((tmp_path / "db" / "prot").glob("*")) or not any(
        p.name.endswith(".part") for p in (tmp_path / "db" / "prot").iterdir())


def test_a_database_that_unpacks_without_its_hash_table_is_refused(tmp_path,
                                                                   monkeypatch):
    """Kraken2 reads an incomplete database and classifies nothing, which arrives as a
    clean sample rather than as an error."""
    served = tmp_path / "served" / "k2_toy.tar.gz"
    served.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        info = tarfile.TarInfo("taxo.k2d")
        data = b"not the hash table"
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    served.write_bytes(buf.getvalue())
    refset = fetchdb.RefSet(
        key="classifier", what="a toy tarball", why="tests the completeness check",
        sources=(fetchdb.Source(_as_url(served), "k2_toy.d",
                                untar_requiring="hash.k2d"),),
        fills="k2", many=True)
    monkeypatch.setitem(fetchdb.BY_KEY, "classifier", refset)
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["classifier"], tmp_path / "db", Path("."))
    assert "hash.k2d" in str(e.value) and "clean sample" in str(e.value)


def test_a_member_that_would_unpack_outside_the_destination_is_refused(tmp_path,
                                                                      monkeypatch):
    served = tmp_path / "served" / "evil.tar.gz"
    served.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        info = tarfile.TarInfo("../../escaped")
        info.size = 0
        tar.addfile(info, io.BytesIO(b""))
    served.write_bytes(buf.getvalue())
    refset = fetchdb.RefSet(
        key="classifier", what="a tarball with a path traversal",
        why="an index is unpacked with the privileges of whoever runs the tool",
        sources=(fetchdb.Source(_as_url(served), "evil.d",
                                untar_requiring="hash.k2d"),),
        fills="k2", many=True)
    monkeypatch.setitem(fetchdb.BY_KEY, "classifier", refset)
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["classifier"], tmp_path / "db", Path("."))
    assert "unpacks outside" in str(e.value)
    assert not (tmp_path / "db" / "escaped").exists()


def test_a_diamond_set_with_no_diamond_refuses_naming_the_set(tmp_path, monkeypatch):
    served = tmp_path / "served" / "prot.faa.gz"
    served.parent.mkdir(parents=True, exist_ok=True)
    served.write_bytes(gzip.compress(b">p\nMKV\n"))
    refset = fetchdb.RefSet(
        key="sprot", what="a DIAMOND index", why="tests the missing-binary refusal",
        sources=(fetchdb.Source(_as_url(served), "prot.faa.gz", rolling=True),),
        fills="prot/sprot.dmnd", diamond_makedb=True)
    monkeypatch.setitem(fetchdb.BY_KEY, "sprot", refset)
    with pytest.raises(fetchdb.FetchError) as e:
        fetchdb.fetch(["sprot"], tmp_path / "db", Path("."), diamond=None)
    assert "sprot" in str(e.value) and "another version" in str(e.value)


def test_the_printed_config_block_is_a_config_this_package_accepts(tmp_path,
                                                                   monkeypatch):
    """What `fetch-db` prints must load, or the handover is prose.

    The block carries absolute paths out of the manifest, and `classifier` carries
    several in the order they were fetched, because that order is the tier order.
    """
    manifest = {
        "dest": str(tmp_path / "db"),
        "sets": {
            "classifier": {"many": True, "paths": [str(tmp_path / "db" / "k2" / "a.d"),
                                                   str(tmp_path / "db" / "k2" / "b.d")]},
            "protein": {"many": False,
                        "paths": [str(tmp_path / "db" / "prot" / "viral.dmnd")]},
        },
    }
    block = fetchdb.config_block(manifest)
    cfg_path = tmp_path / "run.toml"
    cfg_path.write_text("[inputs]\nlabel = \"x\"\n\n" + block)
    cfg = load(cfg_path)
    assert [p.name for p in cfg.database_list("classifier")] == ["a.d", "b.d"]
    assert cfg.database("protein").name == "viral.dmnd"
    # `protein_fasta` comes along with the protein index, because the name map reads the
    # FASTA the index was built from rather than the index.
    assert "protein_fasta" in cfg.databases
    # And `classifier_contigs` comes along with the classifier list, naming the one of the
    # four that identifies contigs. Without it in the printed block, a run config taken
    # from `fetch-db` unedited reaches the contig stage with a list where it needs one
    # path, which is how a five-seed launch died 57 minutes in.
    assert cfg.database("classifier_contigs").name == f"{fetchdb.K2_DB_CONTIGS}.d"


def test_the_contig_classifier_is_one_of_the_databases_actually_fetched():
    """A named operating point that nothing downloads is a path that will not exist.

    It is named rather than taken by position, because a position moves when the list
    does, and `experiments/ec2_prod.py` carries the same value for the launcher path.
    """
    import sys
    assert fetchdb.K2_DB_CONTIGS in fetchdb.K2_DBS
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
    import ec2_prod
    assert ec2_prod.K2_DB_CONTIGS == fetchdb.K2_DB_CONTIGS, (
        "the launcher and the package name different contig classifiers, so a run "
        "through one is not comparable to a run through the other")


def test_asking_for_one_path_where_the_config_names_several_is_refused(tmp_path):
    cfg_path = tmp_path / "run.toml"
    cfg_path.write_text('[databases]\nclassifier = ["a", "b"]\n')
    cfg = load(cfg_path)
    with pytest.raises(ConfigError) as e:
        cfg.database("classifier")
    assert "asked for one" in str(e.value)
    assert [str(p) for p in cfg.database_list("classifier")] == ["a", "b"]


def test_the_committed_example_config_loads_and_states_its_bases():
    """The README points a new user at `examples/berlin.toml`, so it has to load.

    A documented example that no longer parses is worse than none: it is the first thing
    a reader copies, and a stated library fact with no account of where it came from is
    the default this repo spent two days removing from twenty scripts.
    """
    from hapax.config import LIBRARY_FACTS
    path = Path(__file__).resolve().parents[1] / "examples" / "berlin.toml"
    cfg = load(path)
    assert cfg.label == "SRR23973193"
    for key in LIBRARY_FACTS:
        if key in cfg.library:
            assert cfg.library[f"{key}_basis"], key
    # The classifier list is the tier order, and the example names all four so that the
    # committed tier table is reproducible from it.
    assert len(cfg.database_list("classifier")) == 4
    # Every key the example names must be one `fetch-db` can produce or the run makes.
    for key in cfg.databases:
        assert (key in fetchdb.BY_KEY or key in fetchdb.NOT_FETCHABLE
                or any(key in s.extra for s in fetchdb.INVENTORY)), key


def test_a_run_names_every_absence_in_one_refusal_and_creates_nothing(tmp_path,
                                                                     monkeypatch,
                                                                     capsys):
    """One pass over the config, not one category per pass.

    The tool check and the path check used to be two passes with the driver's preflight
    in between, so a caller with a new machine and a new config fixed three installs,
    re-ran, and only then learned their FASTQ path had a typo. Nothing is created either:
    a run directory left behind by a refused run is a directory a later run will read.
    """
    from hapax.cli import main
    cfg = tmp_path / "run.toml"
    cfg.write_text(
        '[inputs]\nlabel = "x"\nr1 = "no/such_1.fq.gz"\nr2 = "no/such_2.fq.gz"\n\n'
        '[databases]\nrrna = "no/silva.fa"\nclassifier = ["no/a.d", "no/b.d"]\n')
    out = tmp_path / "run"
    monkeypatch.setenv("HAPAX_REPO", str(Path(__file__).resolve().parents[1]))
    code = main(["run", "--config", str(cfg), "--out", str(out), "--skip-prepare"])
    assert code == 2
    message = capsys.readouterr().err
    # Both categories in the one message, each naming the table it came from.
    for needle in ("[inputs] r1", "[inputs] r2", "[databases] rrna",
                   "[databases] classifier"):
        assert needle in message, message
    assert message.count("[databases] classifier") == 2, (
        "a key naming two indexes reports one absence; a reader fixes one path and "
        "finds the next on the following run")
    assert not out.exists(), (
        "the refused run created its output directory, which a later run would read as "
        "a run that happened")


def test_a_configs_sample_facts_are_keys_the_page_can_speak():
    """A fact the sentence has no clause for REFUSES the page, so it cannot ship in an
    example or be written by a launcher.

    `page_descriptor.check` raises on an unknown key, which is right: a fact somebody
    wrote down that the renderer drops without saying so is the asymmetry that once lost
    "South Korea" off the one library here not in Europe. This is that rule applied one
    layer out, to the places a `sample_facts` block is now AUTHORED: the committed
    example, and the TOML that the production launcher writes on the instance.
    """
    import sys
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "experiments"))
    import page_descriptor

    cfg = load(root / "examples" / "berlin.toml")
    page_descriptor.check("SRR23973193", cfg.inputs["sample_facts"])

    import tomllib

    import ec2_prod
    import library_facts
    for accession in ("SRR23973193", "SRR23702504", "SRR34217566"):
        block = ec2_prod._sample_facts_toml(
            library_facts.library(accession).get("sample") or {})
        facts = tomllib.loads(block)
        page_descriptor.check(accession, facts)
    # And the empty case, which is the production case: water nobody deposited from.
    facts = tomllib.loads(ec2_prod._sample_facts_toml({}))
    page_descriptor.check("nothing-known", facts)
    assert "does not know" in facts["basis"]
