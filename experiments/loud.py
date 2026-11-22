"""Fail loudly: the checks every stage that produces a committed number runs before it starts.

WHY THIS EXISTS. This project has been bitten three times by one failure mode: a stage that
did less than it was asked and reported success. `hapax.py` validated three indexes and
never passed them on, so a run used the defaults while its manifest named the indexes it
had been asked for. It never passed the protein index to the per-read report, so `strains`
wrote an empty file that read as "no divergent relatives" when nothing had been searched.
`residue_report.py` returned a motif rate of 0.0 when no profile database was on disk. In
every case the file was well-formed and the number was a lie.

THE RULES, which every helper below enforces:

  1. **Refuse to start.** A requested input, index or tool that is absent, empty or not
     executable stops the run in its first seconds, with every problem listed at once, not
     the first one an hour in. Asking for an index is a statement that the number depends on
     it; quietly skipping it changes the number and keeps the label.
  2. **Record what was consulted, not what was asked for.** `Consulted` writes one entry per
     index actually opened, with its size and modification time read at the moment it was
     opened, and how many queries went in and how many came back with a hit. A manifest's
     `parameters` say what was requested; `consulted` says what happened. A reader, or
     `hapax.py`, compares the two.
  3. **"Searched, zero hits" and "not searched" are different values.** The first is a
     count of 0 beside a query count; the second is a string saying why. Neither is ever a
     blank or a bare 0.0.
  4. **A tool's non-zero exit is fatal and its stderr is shown.** `subprocess.run(...,
     capture_output=True, check=True)` raises with the stderr swallowed, which leaves a log
     saying only that something failed.

A cached derivative of an input (a minimap2 index, a converted FASTA) is the same hazard in
another form: a work directory reused across runs serves the index built from LAST run's
reference under this run's name. `cache_is_current` keys the cache to its sources.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

#: Exit status for "refused to start". Distinct from 1 (a stage failed while running) so a
#: wrapper can tell a misconfigured run from a broken one without reading the log.
EXIT_PREFLIGHT = 2


class Preflight:
    """Collect every missing or unusable input, then refuse with all of them at once."""

    def __init__(self, script: str):
        self.script = script
        self.problems: list[str] = []

    def file(self, path: Path | None, what: str) -> None:
        if path is None:
            self.problems.append(f"{what}: no path given")
            return
        p = Path(path)
        if not p.exists():
            self.problems.append(f"{what}: {p} does not exist")
        elif p.is_dir():
            self.problems.append(f"{what}: {p} is a directory, expected a file")
        elif p.stat().st_size == 0:
            self.problems.append(f"{what}: {p} is empty (0 bytes)")

    def kraken2_db(self, db: Path) -> None:
        """A Kraken2 database is a directory of three files, and all three are read."""
        if not Path(db).is_dir():
            self.problems.append(f"Kraken2 database: {db} is not a directory")
            return
        for f in ("hash.k2d", "taxo.k2d", "opts.k2d"):
            self.file(Path(db) / f, f"Kraken2 database {Path(db).name}")

    def exe(self, path: Path | None, what: str) -> None:
        if path is None:
            self.problems.append(f"{what}: no path given")
            return
        p = Path(path)
        if not p.is_file():
            self.problems.append(f"{what}: {p} does not exist")
        elif not os.access(p, os.X_OK):
            self.problems.append(f"{what}: {p} is not executable")

    def module(self, name: str, what: str) -> None:
        import importlib.util
        if importlib.util.find_spec(name) is None:
            self.problems.append(f"{what}: python module `{name}` is not importable")

    def require(self, ok: bool, message: str) -> None:
        if not ok:
            self.problems.append(message)

    def done(self, preflight_only: bool = False) -> None:
        """Refuse if anything was wrong; with `--preflight`, stop here either way."""
        if self.problems:
            print(f"{self.script}: REFUSING TO START - {len(self.problems)} problem(s):",
                  file=sys.stderr, flush=True)
            for p in self.problems:
                print(f"  - {p}", file=sys.stderr, flush=True)
            raise SystemExit(EXIT_PREFLIGHT)
        if preflight_only:
            print(f"{self.script}: preflight ok", flush=True)
            raise SystemExit(0)


def describe(path: Path) -> dict:
    """Size and modification time of an input, read when it is opened.

    Not a content hash: hashing a 16 GB Kraken2 database per run costs minutes. Size plus
    mtime is enough to tell two different files apart in practice, which is the question.
    """
    p = Path(path)
    if p.is_dir():
        files = [f for f in p.rglob("*") if f.is_file()]
        return {"path": str(p), "bytes": sum(f.stat().st_size for f in files),
                "files": len(files),
                "mtime": max((int(f.stat().st_mtime) for f in files), default=None)}
    st = p.stat()
    return {"path": str(p), "bytes": st.st_size, "mtime": int(st.st_mtime)}


class Consulted:
    """What a stage actually opened, and what came back from each.

    `searched(...)` is called at the point an index is used, never at the point it is
    requested, so an index that was asked for and then skipped by any code path is absent
    from this record and the comparison against `parameters` exposes it.
    """

    def __init__(self):
        self.entries: list[dict] = []

    def searched(self, role: str, path: Path, queries: int, hits: int, **kw) -> None:
        if queries < 0 or hits < 0 or hits > queries:
            raise SystemExit(f"internal: impossible search record for {role} {path}: "
                             f"{hits} hits from {queries} queries")
        self.entries.append({"role": role, **describe(path), "state": "searched",
                             "queries": queries, "hits": hits, **kw})

    def not_searched(self, role: str, path: Path | None, why: str) -> None:
        self.entries.append({"role": role, "path": str(path) if path else None,
                             "state": "not searched", "why": why})

    def paths(self, state: str = "searched") -> set[str]:
        return {e["path"] for e in self.entries if e["state"] == state}

    def as_list(self) -> list[dict]:
        return list(self.entries)


def run_tool(cmd: list, what: str, **kw) -> subprocess.CompletedProcess:
    """Run an external tool; on non-zero exit, die with its stderr rather than without it."""
    kw.setdefault("stderr", subprocess.PIPE)
    proc = subprocess.run([str(c) for c in cmd], text=True, **kw)
    if proc.returncode != 0:
        err = (proc.stderr or "").strip().splitlines()[-20:]
        print(f"{what} FAILED with exit {proc.returncode}:\n  $ {' '.join(map(str, cmd))}",
              file=sys.stderr, flush=True)
        for line in err:
            print(f"  | {line}", file=sys.stderr, flush=True)
        raise SystemExit(f"{what} exited {proc.returncode}")
    return proc


def require_output(path: Path, what: str, allow_empty: bool = False) -> None:
    """A tool that exits 0 without writing its output has still failed."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"{what}: exited 0 but did not write {p}")
    if not allow_empty and p.stat().st_size == 0:
        raise SystemExit(f"{what}: wrote an empty {p}")


def cache_is_current(cache: Path, sources: list[Path]) -> bool:
    """Whether a cached derivative was built from exactly these sources.

    A sidecar `<cache>.source.json` records the sources' paths, sizes and mtimes when the
    cache is built. A cache with no sidecar predates this check and is treated as stale:
    rebuilding costs seconds to minutes, and serving a stale index under the current
    reference's name is the failure this module exists to stop.
    """
    side = Path(str(cache) + ".source.json")
    if not Path(cache).exists() or not side.exists():
        return False
    try:
        return json.loads(side.read_text()) == [describe(s) for s in sources]
    except Exception:
        return False


def mark_cache(cache: Path, sources: list[Path]) -> None:
    Path(str(cache) + ".source.json").write_text(json.dumps([describe(s) for s in sources]))
