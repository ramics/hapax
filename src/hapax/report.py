"""The report, built from a run directory, refusing rather than printing an absence.

WHAT CHANGES HERE. PLAN.md 1.4: the renderer treats a missing input as valid. `load()`
returns None for an absent file and the panels print what is missing, so a run that
produced almost nothing still built a page that rendered. On a fixture that is a test that
passes on a broken tool; on real water it tells a reader that something is missing from
THEIR SAMPLE when what is missing is a stage.

So "run complete" now means every declared section rendered from a file that meets its
schema, and anything less exits non-zero. `--allow-partial` is the way to look at an
incomplete run on purpose, and it marks the page rather than quietly allowing it.

WHY A RUN DIRECTORY IS NOW ENOUGH. The generator takes a page descriptor, and the
descriptor carried the sample's own facts, so a run on water nobody has deposited from
could not have a page: `experiments/libraries.json` raises for a library it does not know,
which is right, and a production sample is by definition one it does not know. The
descriptor is therefore assembled here from what the run itself recorded - its config, its
manifest - and the facts it carries are marked as supplied rather than derived. Nothing is
guessed: a fact the run does not carry is left out of the sentence, and the page says so.
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from . import scripts
from .registry import Stage
from .schemas import validate
from .stages import ORPHAN_REPORT_INPUTS, R, SECTIONS


class ReportError(Exception):
    """Raised instead of building a page that would state an absence as a measurement."""


@dataclass(frozen=True)
class InputState:
    """One declared output of one stage, as the report finds it."""
    stage: str
    filename: str
    section: str | None
    present: bool
    problems: tuple[str, ...] = ()
    optional: bool = False
    absence_means: str | None = None

    @property
    def usable(self) -> bool:
        return self.present and not self.problems

    @property
    def blocks_its_section(self) -> bool:
        """An absent enrichment does not stop the build; an absent measurement does.

        An output declared optional has a section that renders without it and says what
        that costs. A present-but-INVALID optional file is still a failure: the section
        has no way to say "this file is here and wrong", and silently treating a broken
        file as an absence would be the renderer's old behaviour with an extra step.
        """
        if self.usable:
            return False
        return not (self.optional and not self.present)

    def line(self) -> str:
        if not self.present:
            return (f"{self.filename} is absent. {self.stage} writes it"
                    + (f", and {self.section!r} renders it" if self.section else ""))
        return (f"{self.filename} does not meet its contract with the report: "
                + "; ".join(self.problems))


def _resolve(rundir: Path, filename: str, seed: int) -> Path | None:
    """The file a declared output names in this run directory, or None if absent.

    A declaration may name `x_s{seed}.json`; a run directory may also hold the file for
    another seed, and a glob would find it. It is NOT accepted: a page built from one
    seed's reads and another seed's profile is a mixture, and a mixture is not a result.
    """
    p = rundir / filename.format(seed=seed)
    return p if p.exists() else None


def inspect(rundir: Path, seed: int) -> list[InputState]:
    """Every declared stage output, present or not, checked against its schema."""
    rundir = Path(rundir)
    states: list[InputState] = []
    for name in R.order():
        stage: Stage = R.stages[name]
        if stage.validation or name == "report":
            # Validation output is not on any page (SPEC 13.4), and the report does not
            # render itself.
            continue
        for out in stage.outputs:
            path = _resolve(rundir, out.filename, seed)
            if path is None:
                states.append(InputState(name, out.filename.format(seed=seed),
                                         out.section, False, optional=out.optional,
                                         absence_means=out.absence_means))
                continue
            problems: tuple[str, ...] = ()
            if path.suffix == ".json":
                try:
                    data = json.loads(path.read_text())
                except json.JSONDecodeError as e:
                    problems = (f"is not readable JSON ({e})",)
                else:
                    problems = tuple(validate(out.schema, data))
            elif path.stat().st_size == 0:
                problems = ("is empty",)
            states.append(InputState(name, path.name, out.section, True, problems,
                                     optional=out.optional,
                                     absence_means=out.absence_means))
    return states


def missing_sections(states: list[InputState]) -> dict[str, list[InputState]]:
    """Report sections that cannot be drawn, each with the reason.

    A section is drawable when EVERY output declared as filling it is usable. Any, rather
    than every, would let a section with one of its two inputs render as half a
    measurement, and a half-drawn panel is the thing a reader cannot tell from a complete
    one.
    """
    out: dict[str, list[InputState]] = {}
    for section in SECTIONS:
        bad = [s for s in states if s.section == section and s.blocks_its_section]
        if bad:
            out[section] = bad
    return out


def descriptor(rundir: Path, out_path: Path) -> Path:
    """Assemble the page descriptor from what the run recorded about itself.

    Everything here comes out of the run directory. `sample_facts` is whatever the run's
    config stated about the water, passed to the generator as SUPPLIED rather than derived
    so the page marks it as an override; absent, the page says it knows nothing about the
    water, which is the honest reading for an input nobody has deposited from.
    """
    rundir = Path(rundir)
    manifest = json.loads((rundir / "manifest.json").read_text())
    invocation = {}
    if (rundir / "invocation.json").exists():
        invocation = json.loads((rundir / "invocation.json").read_text())
    accession = (manifest.get("parameters") or {}).get("library")
    if not accession:
        raise ReportError(
            f"{rundir}/manifest.json does not record which library it ran, so a page "
            f"built from it could not say whose water it is.")
    # None, not {}, when the run recorded nothing about the water. {} is a positive
    # claim that nothing is known, and it would silence the registry for a library the
    # registry does know: the first version of this did exactly that and titled the
    # Yeongdok page by its accession, dropping "Yeongdok plant effluent, South Korea"
    # from a run whose facts are committed. The generator consults the registry when
    # this is absent, and states its ignorance when the registry does not know either.
    facts = invocation.get("sample_facts")
    d = {
        "sample": {"accession": accession},
        "rundir": str(rundir),
        # The files outside the run directory that panels may read. Declared as not run,
        # with the reason, rather than left out: a declared absence and an undeclared one
        # look the same on a page and are not the same thing.
        "inputs": {
            "spectrum": {"not_run": "the divergent-contig workup is not yet a stage of "
                                    "hapax run (PLAN.md step 6)"},
            "spectrum_big_index": {"not_run": "no larger-index block was run against "
                                              "these contigs"},
        },
    }
    if facts is not None:
        d["sample_facts"] = facts
    out_path.write_text(json.dumps(d, indent=1) + "\n")
    return out_path


BANNER = """<div style="background:#7a1f1f;color:#fff;padding:14px 18px;font:14px/1.5
 system-ui,sans-serif" role="alert"><strong>Incomplete run, shown with --allow-partial.
 </strong> {n} of {total} sections could not be drawn from this run: {which}. What is
 missing is a stage of this tool, not a measurement of this water.</div>"""


def build(rundir: Path, out: Path | None = None, seed: int = 1,
          allow_partial: bool = False, repo: Path | None = None) -> int:
    """Check the run against the section contract, then build the page.

    The check comes FIRST and the page is not written at all when it fails. Writing the
    page and then reporting the problem leaves a file on disk that renders, and a file
    that renders is the thing somebody reads.
    """
    rundir = Path(rundir).resolve()
    if not (rundir / "manifest.json").exists():
        raise ReportError(f"{rundir} is not a hapax run directory: no manifest.json")
    states = inspect(rundir, seed)
    absent = missing_sections(states)
    if absent and not allow_partial:
        lines = [f"the report refuses: {len(absent)} of {len(SECTIONS)} sections cannot "
                 f"be drawn from {rundir}"]
        for section, bad in absent.items():
            lines.append(f"  {section}")
            lines += [f"    {s.line()}" for s in bad]
        lines.append("  Pass --allow-partial to build the page anyway; it is then marked "
                     "as incomplete, because a page that renders is a page somebody "
                     "reads.")
        if ORPHAN_REPORT_INPUTS:
            lines.append(f"  Separately, {len(ORPHAN_REPORT_INPUTS)} files the page reads "
                         f"have no stage at all (PLAN.md step 6): "
                         f"{', '.join(ORPHAN_REPORT_INPUTS)}")
        raise ReportError("\n".join(lines))

    repo = Path(repo) if repo else _repo()
    page = descriptor(rundir, rundir / "page_descriptor.json")
    out = Path(out) if out else rundir / "report.html"
    # Without -I, as before: the generator was always started this way, and isolation is
    # a change to what it can import that this refactor does not make.
    proc, generator = scripts.run(
        repo, "make_web.py", [page, "--rundir", rundir, "--out", out, "--seed", seed],
        isolated=False)
    sys.stderr.write(proc.stderr)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        raise ReportError(f"the generator refused to build {out}")

    if absent:
        _mark_partial(out, absent)
    # The render check, written where a test can read it. The generator runs one and
    # prints it; a printed check is not a check anybody else can assert on.
    check = _render_check(repo, out)
    (rundir / "render_check.json").write_text(json.dumps({
        "page": str(out), "seed": seed, "allow_partial": allow_partial,
        "generator": generator,
        "sections_declared": SECTIONS,
        "sections_not_drawable": sorted(absent),
        # Enrichments this run did not produce, each with what the page shows instead.
        # Recorded rather than left out, because an absence a reader cannot see is the
        # same as an absence nobody declared.
        "enrichments_absent": [{"stage": s.stage, "file": s.filename,
                                "section": s.section, "means": s.absence_means}
                               for s in states if s.optional and not s.present],
        "inputs": [{"stage": s.stage, "file": s.filename, "section": s.section,
                    "present": s.present, "problems": list(s.problems),
                    "optional": s.optional}
                   for s in states],
        **(check or {"top": None, "render_check": _why_no_render_check()}),
    }, indent=1) + "\n")
    if check is not None:
        undrawn = [s for s in SECTIONS if s not in (check.get("top") or [])]
        if undrawn and not allow_partial:
            raise ReportError(
                f"the page was written and {len(undrawn)} declared sections never "
                f"reached the DOM: {undrawn}. A section that resolved and did not render "
                f"is the one case file checks cannot see.")
    return 0


def _why_no_render_check() -> str:
    """Which of the two absences this is, because they are different machines to fix.

    The message said "node is not installed" for both, and on the machine where the job
    runs node IS installed and the Playwright module is not. One sentence naming the wrong
    cause sends the next reader to the wrong fix, which is the same failure as a census
    that cannot tell an empty table from a search that never ran.
    """
    if shutil.which("node") is None:
        return ("node is not installed, so the page was not run. The check loads the "
                "written page in a browser and reports which sections reached the DOM.")
    return ("node is installed and its playwright module is not, so the page was not "
            "run: `npm install playwright && npx playwright install chromium`.")


def _mark_partial(page: Path, absent: dict) -> None:
    """Put the incompleteness on the page itself, where its reader is.

    Injected after the fact rather than rendered by the generator, because the generator
    is still a script outside the package (PLAN.md step 7 moves it in, and the banner
    becomes native then). It is injected BEFORE the render check runs, so the check sees
    the page a reader will see.
    """
    html = page.read_text()
    banner = BANNER.format(n=len(absent), total=len(SECTIONS),
                           which=", ".join(sorted(absent)))
    # The page is a fragment with no <body>: it opens with its own <title> and <style>,
    # and the content hangs off one element with id="app" that the page's own script
    # fills. The banner goes immediately before that element, which is the top of what a
    # reader sees, and it is static markup so it survives the script that replaces the
    # element's children. Matched by pattern rather than by an exact string, because the
    # element is `<div class="wrap" id="app">` and an exact string would silently stop
    # matching the first time a class was added.
    m = re.search(r'<\w+[^>]*\bid="app"', html)
    if m is None:
        raise ReportError(f"{page} has no element with id=\"app\" to mark as incomplete")
    page.write_text(html[:m.start()] + banner + html[m.start():])


def _render_check(repo: Path, page: Path) -> dict | None:
    """Run the generator's own render check on the page as written."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "hapax_make_web", repo / "experiments/make_web.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    sys.path.insert(0, str(repo / "experiments"))
    spec.loader.exec_module(module)
    return module.render_check(page)


def _repo() -> Path:
    from .cli import repo_root

    return repo_root()
