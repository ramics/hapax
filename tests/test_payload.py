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


def test_every_repo_file_the_renderer_reads_is_shipped():
    import ec2_prod

    wanted = sorted(set(REPO_READ.findall((ROOT / "experiments/make_web.py").read_text())))
    assert wanted, "found no repo-relative reads in make_web.py, so this test reads nothing"
    shipped = set(ec2_prod.payload_files(ROOT, "SRR34217566"))
    missing = [w for w in wanted if w not in shipped]
    assert not missing, (
        f"the renderer reads {missing} from the repo and the payload does not carry "
        f"them, so the generator will refuse on an instance after every stage has run")


def test_the_files_it_reads_are_actually_in_the_tree():
    """A shipped path that does not exist is a payload entry that cannot arrive."""
    wanted = sorted(set(REPO_READ.findall((ROOT / "experiments/make_web.py").read_text())))
    absent = [w for w in wanted if not (ROOT / w).exists()]
    assert not absent, f"{absent} is read by the renderer and is not in the tree"
