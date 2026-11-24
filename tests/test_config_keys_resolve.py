"""Every database key the package asks for must resolve against a production config.

WHY THIS EXISTS. `prepare.split` asked for `classifier` as a single path. On the end-to-end
fixture that works, because the fixture has exactly one toy Kraken2 library and names it as a
bare string. On production water `classifier` is the cascade's ordered list of four, and
`Config.database` refuses a list rather than silently taking its first member. So the fault
could not appear in CI and appeared instead 57 minutes into a paid run, after the install,
39 GB of reference sets, the reads and the five-seed trim measurement had all succeeded.

The fixture being smaller than production is the point of a fixture. What was missing is
anything that reads the package's own asks against a production-shaped config. This does
that from the committed example, which is production-shaped by construction: the `classifier`
list has four members there because the committed tier table was produced with four.

It is a static read of the source, so it needs no instance, no network and no indexes.
"""

from __future__ import annotations

import re
from pathlib import Path

from hapax.config import ConfigError, load

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/berlin.toml"

SINGLE = re.compile(r"""\bcfg\.database\(\s*["'](?P<key>\w+)["']""")
SEVERAL = re.compile(r"""\bcfg\.database_list\(\s*["'](?P<key>\w+)["']""")


def asks() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    single: dict[str, set[str]] = {}
    several: dict[str, set[str]] = {}
    for path in sorted((ROOT / "src/hapax").glob("*.py")):
        text = path.read_text()
        for m in SINGLE.finditer(text):
            single.setdefault(m.group("key"), set()).add(path.name)
        for m in SEVERAL.finditer(text):
            several.setdefault(m.group("key"), set()).add(path.name)
    return single, several


def test_the_package_asks_for_keys_the_example_config_can_answer():
    cfg = load(EXAMPLE)
    single, several = asks()
    # `several` is read for shape only and is allowed to be empty: `cli.py` resolves its
    # index flags through a loop variable, so a literal-key search cannot see them. That
    # loop is covered by the flag assertion in the next test.
    assert single, "found no database asks, so this test reads nothing"

    problems = []
    for key, where in sorted(single.items()):
        if key not in cfg.databases:
            # Not every key is in every run: `cross_sample` and `nucleotide` are declared
            # absent with a reason on a sample with no panel. What must not happen is a
            # key the config DOES carry being asked for in the wrong shape.
            continue
        try:
            cfg.database(key)
        except ConfigError as exc:
            problems.append(f"{', '.join(sorted(where))} asks for one {key}: {exc}")
    for key, where in sorted(several.items()):
        if key not in cfg.databases:
            continue
        cfg.database_list(key)                   # a bare string is a list of one
    assert not problems, "\n".join(problems)


def test_the_two_classifier_keys_are_asked_for_in_the_shapes_they_hold():
    """The specific confusion, pinned: one list for the cascade, one path for contigs.

    Stated as its own test because the two keys differ by one word and the failure is
    silent on any config that names a single database.
    """
    single, several = asks()
    assert "classifier" not in single, (
        f"{sorted(single.get('classifier', ()))} asks for the cascade's tier list as a "
        f"single path. It has four members on production water.")
    assert "classifier_contigs" not in several, (
        f"{sorted(several.get('classifier_contigs', ()))} asks for the contig "
        f"classifier as a list. classify_contigs.py reads one database and writes its "
        f"name into every row.")
    # The cascade's own flag, which `cli.py` fills through a loop variable: the whole list
    # reaches `--k2-dbs`, because the driver's refusal on an empty one exists so that a
    # tier can never run against an index nobody chose.
    cli = (ROOT / "src/hapax/cli.py").read_text()
    assert '"--k2-dbs"' in cli and "database_list(name)" in cli, (
        "cli.py no longer fills --k2-dbs from database_list, so the cascade may be "
        "running against fewer indexes than the config names")
    cfg = load(EXAMPLE)
    assert len(cfg.database_list("classifier")) == 4
    assert cfg.database("classifier_contigs").name.endswith(".d")
