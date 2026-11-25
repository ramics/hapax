"""Files the renderer reads from the REPO must be in the EC2 payload.

WHY THIS EXISTS. `make_web.py` reads two files by repo-relative path rather than out of the
run directory: the Methods table and the virus vocabulary. `experiments/` is shipped as a
directory, so the vocabulary travelled; `results/web/methods.json` is under `results/`,
which was not in the payload at all. The generator therefore threw a `FileNotFoundError`
after every stage of seed 1 had finished — 56 minutes of profile, cascade, reads, linkage
and the rest — and `hapax run` refused to write a page with no Methods section, which is
the correct refusal to the wrong problem.

The payload's own docstring says the lesson twice already: the thing that was missing was
the thing nobody listed, so the list is derived from the registry that names it. This does
that for the renderer, by reading its source for the paths it resolves against the repo.

It is a static read, so it needs no instance and no network.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

#: `ROOT / "<path>"` in the renderer, which is the only shape it uses for a repo file.
REPO_READ = re.compile(r'ROOT\s*/\s*"(?P<rel>[^"]+)"')
#: A plain `import x` of a module that is itself a file in `experiments/`.
LOCAL_IMPORT = re.compile(r"^(?:import|from)\s+(?P<mod>[a-z_][a-z0-9_]*)", re.M)


def renderer_modules() -> list[Path]:
    """The renderer AND the modules it imports from `experiments/`, transitively.

    Reading `make_web.py` alone was the narrow version of this test and it passed while
    the generator still refused on an instance: `web_controls.py` resolves
    `results/web/controls.json` against the repo in exactly the same shape, one import
    away, and `results/` is not shipped as a directory. A module that reads a repo file
    is not a property of which file holds the `main()`.
    """
    seen: dict[str, Path] = {}
    queue = ["make_web"]
    while queue:
        mod = queue.pop()
        path = ROOT / "experiments" / f"{mod}.py"
        if mod in seen or not path.exists():
            continue
        seen[mod] = path
        queue.extend(LOCAL_IMPORT.findall(path.read_text()))
    return [seen[m] for m in sorted(seen)]


def repo_reads() -> list[str]:
    """Committed repo files only.

    `data/` is gitignored and fetched onto the machine by `hapax fetch-db`, so a path
    under it is never a payload entry and the modules the walk reaches carry several as
    argparse defaults. The question this test asks is about files the COMMIT describes.
    A path with no suffix is an output directory an argparse default names, not a file
    anything reads.
    """
    found: set[str] = set()
    for path in renderer_modules():
        found.update(REPO_READ.findall(path.read_text()))
    return sorted(f for f in found
                  if not f.startswith("data/") and Path(f).suffix)


def test_the_renderer_and_its_modules_are_both_read():
    """The transitive walk is the point of the test, so a walk that collapsed fails."""
    names = [p.name for p in renderer_modules()]
    assert "make_web.py" in names and "web_controls.py" in names, names


def test_every_repo_file_the_renderer_reads_is_shipped():
    import ec2_prod

    wanted = repo_reads()
    assert wanted, "found no repo-relative reads in the renderer, so this test reads nothing"
    shipped = set(ec2_prod.payload_files(ROOT, "SRR34217566"))
    missing = [w for w in wanted if w not in shipped]
    assert not missing, (
        f"the renderer reads {missing} from the repo and the payload does not carry "
        f"them, so the generator will refuse on an instance after every stage has run")


def test_the_files_it_reads_are_actually_in_the_tree():
    """A shipped path that does not exist is a payload entry that cannot arrive."""
    wanted = repo_reads()
    absent = [w for w in wanted if not (ROOT / w).exists()]
    assert not absent, f"{absent} is read by the renderer and is not in the tree"
