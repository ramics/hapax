"""The stage registry: what every stage declares, and the order that follows from it.

WHY A REGISTRY RATHER THAN A LIST. The stage order was typed in three places and the
launcher's copy is the one that ran, so a validation stage ran in a production pass and
two stages left production in one copy and not the others. A list orders stages; it does
not say what a stage needs, what it writes, or which part of the report goes empty when
it does not run, so none of those could be checked. Every defect in PLAN.md 1.3 to 1.6 is
a thing nothing was in a position to check.

So a stage declares four things and the runner derives the rest:

  * what it NEEDS, as the names of other stages. The order is computed from that, so
    there is one order and it cannot disagree with itself.
  * what it WRITES, as filenames in the run directory with a schema name each. A stage
    that exits zero and writes none of them has failed, which is the distinction
    PLAN.md 1.5 is about.
  * which REPORT SECTION each output feeds, so that a section with no stage behind it is
    a build error rather than a page that renders a census of an absent file. Three
    sections had no producer at all and the page said so to a reader of their own water.
  * whether it is VALIDATION rather than output. A stage that scores the run against
    sequence whose distance was fixed before it has no answer on water nobody has
    deposited from, which is every sample production is for, so it is not in the default
    set (SPEC 13.4).

WHAT THIS MODULE DOES NOT DO. It does not run anything and it holds no paths. The
declarations here are about the shape of a run; what a stage costs, which database it
consults and where that database lives are in the run configuration, because they are
properties of a machine and an input rather than of the pipeline.
"""
from __future__ import annotations

from dataclasses import dataclass, field


class RegistryError(Exception):
    """Raised rather than falling back to declaration order, which is not an order."""


@dataclass(frozen=True)
class Output:
    """One file a stage writes, and what reads it.

    `filename` may contain `{seed}`, because most stages write one file per seed and the
    seed is a run parameter rather than a stage parameter.

    `section` names the report section this output fills, or is None for an output
    nothing on the page renders. None is a real answer and not an omission: the cascade's
    own tier table is read by other stages rather than drawn.
    """
    filename: str
    schema: str
    section: str | None = None
    #: An output whose absence the section HANDLES, with a stated consequence, rather
    #: than a hole where a measurement should be.
    #:
    #: The distinction is the whole of PLAN.md step 6 and it took writing the
    #: declarations down to see. Three report inputs had no stage and the sections that
    #: read them rendered as a census of an absent file: that is a hole in the tool, and
    #: the page told a reader of their own water that something was missing from it.
    #: Three more had no stage and their sections render WITHOUT them and say what that
    #: costs - the distance panel falls back to the whole-contig hit, the subject table
    #: shows an accession where it could show an organism. Those are enrichments.
    #:
    #: Treating the two alike in either direction is wrong. Calling an enrichment a hole
    #: would make every run that lacks a Pfam database refuse; calling a hole an
    #: enrichment is what the page was doing when a run with three empty sections built
    #: cleanly. So an optional output is recorded as absent WITH its consequence, and a
    #: required one stops the build.
    optional: bool = False
    #: What the page shows instead, where `optional`. Required with it: an absence with
    #: no stated consequence is the thing a reader cannot weigh.
    absence_means: str | None = None

    def __post_init__(self):
        if self.optional and not self.absence_means:
            raise ValueError(
                f"{self.filename} is declared optional with no absence_means. An output "
                f"whose absence is allowed has to say what the page shows instead, or "
                f"nobody reading the page can tell what they are missing.")


@dataclass(frozen=True)
class Stage:
    """One stage, as declared. `why` is the reason it exists, not a description of it."""
    name: str
    why: str
    outputs: tuple[Output, ...] = ()
    needs: tuple[str, ...] = ()
    validation: bool = False
    #: The experiments/ script this stage currently shells out to, or None where the
    #: driver does the work itself. PLAN.md step 7 replaces these with imports; until
    #: then the name is recorded here so that the dependency is declared rather than
    #: implicit in a subprocess call, which is PLAN.md 1.7.
    script: str | None = None
    #: Whether the legacy driver runs this stage, or the package does.
    #:
    #: Recorded rather than inferred, because the alternative is a stage that silently
    #: does not run: the stage list was typed in three places and the launcher's copy is
    #: the one that ran, so a stage left out of one copy ran nowhere and nothing said so.
    #: `tests/test_registry.py` checks this field against the driver's own list, so a
    #: stage marked as the driver's that the driver does not have fails the suite.
    #: Step 7 ends the split by moving everything into the package.
    in_driver: bool = True

    def files(self, seed: int) -> list[str]:
        """The distinct files this stage writes. Distinct, because one file can fill two
        report sections - the profile fills both the watchlist and the composition - and
        a stage that looked like it wrote the same file twice would read as a defect."""
        return list(dict.fromkeys(o.filename.format(seed=seed) for o in self.outputs))


@dataclass
class Registry:
    stages: dict[str, Stage] = field(default_factory=dict)

    def add(self, stage: Stage) -> Stage:
        if stage.name in self.stages:
            raise RegistryError(f"{stage.name} is declared twice")
        self.stages[stage.name] = stage
        return stage

    def order(self, names: list[str] | None = None) -> list[str]:
        """The order the declarations imply, deterministic and refusing a cycle.

        Kahn's algorithm with the ready set taken in DECLARATION order, so the result is
        one order rather than one of many: a run whose stage order depends on dictionary
        iteration is a run whose manifest cannot be compared with another one.

        A stage asked for whose need was not asked for is an error and not a silent
        inclusion. Pulling the need in would mean `--stages cascade` quietly running four
        other stages, and a subset run that is not the subset the caller named is how
        `--stages` came to hide the reason a stage was missing.
        """
        wanted = list(self.stages) if names is None else list(dict.fromkeys(names))
        unknown = [n for n in wanted if n not in self.stages]
        if unknown:
            raise RegistryError(f"no such stage: {', '.join(unknown)}. "
                                f"Declared: {', '.join(self.stages)}")
        chosen = set(wanted)
        missing = {n: [d for d in self.stages[n].needs if d not in chosen] for n in wanted}
        missing = {n: d for n, d in missing.items() if d}
        if missing:
            raise RegistryError(
                "these stages need stages that were not selected: " + "; ".join(
                    f"{n} needs {', '.join(d)}" for n, d in missing.items()))
        done: list[str] = []
        remaining = list(wanted)
        while remaining:
            ready = [n for n in remaining
                     if all(d in done for d in self.stages[n].needs)]
            if not ready:
                raise RegistryError(
                    f"the declarations have a cycle among: {', '.join(remaining)}")
            for n in ready:
                done.append(n)
                remaining.remove(n)
        return done

    def production(self) -> list[str]:
        """Every stage that is output of the tool rather than validation of it."""
        return self.order([n for n, s in self.stages.items() if not s.validation])

    def sections(self) -> dict[str, list[str]]:
        """Report section -> the stages that fill it, in run order."""
        out: dict[str, list[str]] = {}
        for name in self.order():
            for o in self.stages[name].outputs:
                if o.section:
                    out.setdefault(o.section, []).append(name)
        return out

    def output_of(self, filename: str) -> Stage | None:
        for stage in self.stages.values():
            for o in stage.outputs:
                if o.filename == filename:
                    return stage
        return None
