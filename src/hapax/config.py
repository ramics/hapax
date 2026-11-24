"""One TOML file per run: inputs, databases, operating points, and the library's facts.

WHY A CONFIG FILE. `hapax run` has about thirty flags, the per-library facts live in a
JSON registry keyed by accession, and the rest of what a run needs is typed into scripts.
A new input therefore needed a code or JSON edit before it would run at all, which is the
opposite of a tool that works on water nobody has looked at (PLAN.md 1.9).

WHY A STATED FACT STILL NEEDS A BASIS. The registry's rule is that a value which is a
property of one sequencing run never has a default, because a trim of 14 written as
`default=14` deletes fourteen real bases from every mate of a library that has no
non-template prefix, silently, in the direction of less recognition. Moving those values
into a config file would reintroduce exactly that if the file could carry a bare number:
a config copied from another run's config is a default with extra steps. So every library
fact stated here carries `*_basis`, saying how it was established, and this module
REFUSES one that does not. The basis is what a reader needs to know whether the number
belongs to their library.

WHAT OVERRIDES WHAT. A value stated in the config wins over a value measured off the
reads, and a measured value wins over nothing. There is no third source: a config that
names neither a value nor an accession to look one up in is refused rather than handed
the discriminating dataset's.
"""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

#: Library facts that are properties of a sequencing run, so each needs a basis when it
#: is stated by hand. `trim_r2` and `block` are the two this repo has been bitten by.
LIBRARY_FACTS = ("trim_r2", "block", "layout", "pairs")

#: The tables a config may contain, and the keys each may carry. Unknown keys are
#: REFUSED rather than ignored: a typo in a threshold is a silently different operating
#: point, and every number here ships with its operating point or is unreadable.
SCHEMA = {
    # `report_without` names per-read report fields to leave out ON PURPOSE, and the
    # reason it is a setting rather than a fallback is the distinction this whole file is
    # about: a field whose reference is absent stops the run, because a dropped field
    # reads as a zero. Saying so in the config makes the omission a stated operating
    # point that the page carries, instead of a hole.
    "run": {"seed", "seeds", "threads", "out", "pairs", "report_pairs", "stages",
            "allow_partial", "report_without"},
    # `sample_facts` is what the page says the water was: the matrix, the place, the
    # collection date, the molecule, each with how it was established. It is here rather
    # than in the library registry because production water has no registry entry, and
    # the page states its ignorance rather than guessing when it is absent.
    # `assembly` names an assembly produced elsewhere, which is used as it is: two hours
    # of MEGAHIT that only exists on a machine with a dead-man timer gets paid for twice.
    "inputs": {"r1", "r2", "label", "accession", "single_end", "sample_facts",
               "assembly"},
    "library": set(LIBRARY_FACTS) | {f"{k}_basis" for k in LIBRARY_FACTS},
    # `pfam` and `sprot` are the per-gene workup's; `protein_fasta` is the FASTA the
    # protein index was built from, which the name map reads rather than reading the
    # index - a second source for the same names could disagree with the database the
    # search actually used.
    # `classifier` is the read cascade's ORDERED TIER LIST and `classifier_contigs` is the
    # single database that identifies assembled contigs. They are two keys and not one
    # because they are two operating points: the cascade's number is only readable beside
    # the whole list it consulted, while a contig is classified against one database whose
    # name is written into every row of `contig_identification.json`. Reducing the list to
    # its first member for the contig stage, or running the contig stage four times, would
    # both report something other than what the config states.
    "databases": {"nucleotide", "protein", "rrna", "classifier", "classifier_contigs",
                  "cross_sample", "spacers", "pfam", "sprot", "protein_fasta"},
    "tools": {"diamond", "classify", "minimap2", "megahit"},
}


class ConfigError(Exception):
    """Raised rather than falling back. Every fallback in this file would be a default."""


@dataclass
class Config:
    path: Path
    run: dict = field(default_factory=dict)
    inputs: dict = field(default_factory=dict)
    library: dict = field(default_factory=dict)
    databases: dict = field(default_factory=dict)
    tools: dict = field(default_factory=dict)

    @property
    def seed(self) -> int:
        return int(self.run.get("seed", 1))

    @property
    def threads(self) -> int:
        return int(self.run.get("threads", 4))

    @property
    def label(self) -> str:
        """What this run calls itself. The accession where there is one, else the label.

        A fixture and a simulated library have no accession, and inventing one would put
        a string that looks like an archive record on a page.
        """
        return self.inputs.get("accession") or self.inputs.get("label") or self.path.stem

    def fact(self, key: str):
        """A library fact stated in the config, or None where it was not stated.

        None means "measure it", which is a stage's job (PLAN.md step 5). It never means
        another library's value.
        """
        return self.library.get(key)

    def database(self, name: str) -> Path:
        if name not in self.databases:
            raise ConfigError(
                f"[databases] has no {name!r} in {self.path}. A stage that consults an "
                f"index cannot be given one by this code: which index is consulted is "
                f"the operating point, and a run whose index was guessed is a number "
                f"without one.")
        value = self.databases[name]
        if isinstance(value, list):
            raise ConfigError(
                f"[databases] {name} names {len(value)} paths in {self.path}, and this "
                f"caller asked for one. Use database_list where several are meant: a "
                f"set of indexes silently reduced to its first member is a different "
                f"operating point from the one the config states.")
        return Path(value).expanduser()

    def database_list(self, name: str) -> list[Path]:
        """One or several paths for a key, in the order the config gives them.

        The order is load-bearing where a key names a cascade's tiers, so this preserves
        it and never sorts. A single string is a list of one rather than an error,
        because the common case is one index and a config should not have to say `[...]`
        to name it.
        """
        value = self.databases[name] if name in self.databases else None
        if value is None:
            self.database(name)          # raises with the reason
        values = value if isinstance(value, list) else [value]
        return [Path(v).expanduser() for v in values]


def load(path: Path) -> Config:
    """Read and check a run configuration. Refuses rather than filling anything in."""
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"{path} does not exist")
    try:
        raw = tomllib.loads(path.read_text())
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path} is not readable TOML: {e}") from None

    unknown_tables = sorted(set(raw) - set(SCHEMA))
    if unknown_tables:
        raise ConfigError(
            f"{path} has tables this version does not know: {', '.join(unknown_tables)}. "
            f"Known: {', '.join(sorted(SCHEMA))}. An ignored table is a setting the "
            f"caller believes is in force and is not.")
    for table, allowed in SCHEMA.items():
        unknown = sorted(set(raw.get(table) or {}) - allowed)
        if unknown:
            raise ConfigError(
                f"[{table}] in {path} has keys this version does not know: "
                f"{', '.join(unknown)}. Known: {', '.join(sorted(allowed))}.")

    cfg = Config(path=path, **{t: dict(raw.get(t) or {}) for t in SCHEMA})

    # Every stated library fact carries how it was established, or it is refused. This is
    # the registry's rule moved rather than dropped: a config copied from another run's
    # config is a default with extra steps, and the basis is what tells a reader whether
    # the number belongs to their library.
    for key in LIBRARY_FACTS:
        if key in cfg.library and not cfg.library.get(f"{key}_basis"):
            raise ConfigError(
                f"[library] states {key} = {cfg.library[key]!r} in {path} with no "
                f"{key}_basis. {key} is a property of a sequencing run, and a stated "
                f"value with no account of where it came from is the default this repo "
                f"spent two days removing from twenty scripts. Say how it was "
                f"established, or leave it out and let the run measure it.")
    strays = sorted(k for k in cfg.library
                    if k.endswith("_basis") and k[: -len("_basis")] not in cfg.library)
    if strays:
        raise ConfigError(
            f"[library] in {path} gives {', '.join(strays)} for values it does not "
            f"state. A basis for a value nobody set describes nothing.")
    return cfg
