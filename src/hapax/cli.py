"""The `hapax` entry point: one command from a config, and the plan it would run.

WHAT THIS IS AND IS NOT, AT THIS STEP. PLAN.md step 3 is the package skeleton: the entry
point, the registry and the config, wrapping the existing scripts unchanged. So `hapax
run` here translates one config file into the legacy driver's thirty flags and calls it.
That is a real narrowing - the thirty flags become one file, and the order comes from the
registry rather than from a list typed in three places - and it is NOT yet the single
command of PLAN.md 4, because the driver still starts after assembly and demands
`--contigs`, `--unidentified-contigs` and `--cross-sample-ref` from its caller. Step 5
removes those, and `tests/test_end_to_end.py` is what says so out loud.

WHY THE DRIVER IS LOADED BY PATH. The legacy driver is `hapax.py` at the repo root and
this package is `hapax`, so importing it by name would resolve to whichever came first on
`sys.path` - a landmine that depends on test ordering. It is loaded from its own file
instead, which also keeps its content untouched: its `script_sha256` pins five committed
seed manifests, and editing it to make an import work would break the thing the hash is
for.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

from . import __version__
from .config import Config, ConfigError, load
from .registry import RegistryError
from .stages import ORPHAN_REPORT_INPUTS, R

#: The repo the legacy scripts live in. Found from the installed package when the repo is
#: a checkout, and from HAPAX_REPO otherwise, because a stage that shells out to
#: `experiments/subtract_known.py` needs to know where that is and this package does not
#: carry the scripts (PLAN.md step 7 is what ends that).
def repo_root() -> Path:
    env = os.environ.get("HAPAX_REPO")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "hapax.py").exists() and (parent / "experiments").is_dir():
            return parent
    raise ConfigError(
        "cannot find the hapax checkout that holds experiments/ and hapax.py. The "
        "package does not carry the stage scripts yet (PLAN.md step 7); set HAPAX_REPO "
        "to the checkout.")


def legacy_driver():
    """Load the root `hapax.py` from its own file, never by module name."""
    path = repo_root() / "hapax.py"
    spec = importlib.util.spec_from_file_location("hapax_legacy_driver", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def legacy_argv(cfg: Config, r1: Path, r2: Path | None, out: Path,
                stages: list[str]) -> list[str]:
    """The config, as the flags the legacy driver takes.

    Every value comes from the config or is left out. Nothing here supplies a default for
    a per-library fact: the driver and the scripts under it refuse what they were not
    given, and this function passing a plausible number would turn that refusal into a
    wrong answer.
    """
    argv = ["run", str(r1)]
    if r2 is not None:
        argv.append(str(r2))
    argv += ["--seed", str(cfg.seed), "--out", str(out), "--threads", str(cfg.threads)]
    # The legacy driver requires --library and resolves anything it was not given from
    # the accession. A config that states its facts never reaches that lookup, which is
    # how a fixture with no archive accession runs at all.
    argv += ["--library", cfg.label]
    for key, flag in (("trim_r2", "--trim-r2"), ("block", "--block")):
        value = cfg.fact(key)
        if value is not None:
            argv += [flag, str(value)]
    if cfg.inputs.get("single_end"):
        argv.append("--single-end")
    for name, flag in (("nucleotide", "--cross-sample-ref"), ("rrna", "--silva"),
                       ("protein", "--dmnd"), ("classifier", "--k2-dbs"),
                       ("spacers", "--spacers")):
        if name in cfg.databases:
            argv += [flag, str(cfg.database(name))]
    for name, flag in (("classify", "--classify-bin"), ("diamond", "--diamond-bin"),
                       ("minimap2", "--minimap2")):
        if name in cfg.tools:
            argv += [flag, str(Path(cfg.tools[name]).expanduser())]
    if "pairs" in cfg.run:
        argv += ["--pairs", str(cfg.run["pairs"])]
    if "report_pairs" in cfg.run:
        argv += ["--report-pairs", str(cfg.run["report_pairs"])]
    argv += ["--stages", *[s for s in stages if s != "report"]]
    return argv


def selected_stages(cfg: Config) -> list[str]:
    """The stages this run performs, ordered by the registry.

    `[run] stages` names a subset. Its order is ignored, because the order is a property
    of the declarations and not of the caller: a caller who lists them in the wrong order
    is asking for a run that cannot work, and silently reordering is better than running
    it.
    """
    asked = cfg.run.get("stages")
    if asked:
        return R.order(list(asked))
    return R.production()


def cmd_stages(args) -> int:
    """Print the plan: the order, what each stage writes, and which section reads it."""
    cfg = load(args.config) if args.config else None
    names = selected_stages(cfg) if cfg else R.production()
    print(f"{len(names)} stages, ordered by the registry's declarations:\n")
    for i, name in enumerate(names, 1):
        stage = R.stages[name]
        needs = ", ".join(stage.needs) or "nothing"
        print(f"{i:3}. {name}")
        print(f"     needs    {needs}")
        print(f"     writes   {', '.join(stage.files(cfg.seed if cfg else 1)) or 'nothing'}")
        for o in stage.outputs:
            if o.section:
                print(f"     fills    {o.section}")
        if stage.script:
            print(f"     script   experiments/{stage.script}")
    held = [n for n in R.stages if n not in names]
    if held:
        print(f"\nnot in this run: {', '.join(held)}")
    validation = [n for n, s in R.stages.items() if s.validation]
    print(f"validation rather than output, so never in the default set: "
          f"{', '.join(validation)}")
    if ORPHAN_REPORT_INPUTS:
        print(f"\n{len(ORPHAN_REPORT_INPUTS)} files the report reads that no stage "
              f"writes (PLAN.md step 6):")
        for f, who in ORPHAN_REPORT_INPUTS.items():
            print(f"  {f}\n    would come from {who}")
    return 0


def cmd_run(args) -> int:
    cfg = load(args.config)
    out = Path(args.out or cfg.run.get("out") or ".").resolve()
    stages = selected_stages(cfg)
    r1 = Path(args.reads[0]) if args.reads else Path(cfg.inputs["r1"])
    r2 = (Path(args.reads[1]) if len(args.reads) > 1
          else Path(cfg.inputs["r2"]) if cfg.inputs.get("r2") else None)
    argv = legacy_argv(cfg, r1, r2, out, stages)
    out.mkdir(parents=True, exist_ok=True)
    # The translation is written into the run directory before the run, so that a run
    # directory says what command produced it even when the config has since changed.
    # The config file itself is copied for the same reason: a path to a config is not a
    # record of one.
    (out / "invocation.json").write_text(json.dumps({
        "hapax_version": __version__,
        "config": str(cfg.path),
        "config_text": cfg.path.read_text(),
        "stages": stages,
        "legacy_argv": argv,
        "note": "the package translates one config into the legacy driver's flags; "
                "PLAN.md step 7 replaces the subprocess with imports",
    }, indent=1) + "\n")
    root = repo_root()
    proc = subprocess.run([sys.executable, str(root / "hapax.py"), *argv], cwd=root)
    return proc.returncode


def cmd_report(args) -> int:
    """Rebuild the page from a run directory alone.

    The strict version of this - a missing section is a build error, a schema per output,
    `--allow-partial` for a banner instead - is PLAN.md step 4. Today it calls the
    existing generator, which needs a page descriptor, and says so rather than pretending
    a run directory is enough.
    """
    root = repo_root()
    rundir = Path(args.rundir).resolve()
    if args.page is None:
        print("hapax report needs --page until PLAN.md step 4: the generator takes a page "
              "descriptor (results/web/pages/<sample>.json) and a run directory, and a "
              "run directory alone does not yet say what the sample is.", file=sys.stderr)
        return 2
    cmd = [sys.executable, str(root / "experiments/make_web.py"), str(args.page),
           "--rundir", str(rundir)]
    if args.out:
        cmd += ["--out", str(args.out)]
    return subprocess.run(cmd, cwd=root).returncode


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="hapax",
        description="Subtract what is identifiable in a metagenome, then characterise "
                    "what is left.")
    ap.add_argument("--version", action="version", version=f"hapax {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="reads to report, from one config file")
    r.add_argument("reads", nargs="*", type=Path,
                   help="R1 and R2. Omitted, they are read from [inputs] in the config")
    r.add_argument("--config", type=Path, required=True)
    r.add_argument("--out", type=Path, default=None)

    p = sub.add_parser("report", help="rebuild the page from a run directory")
    p.add_argument("rundir", type=Path)
    p.add_argument("--page", type=Path, default=None)
    p.add_argument("--out", type=Path, default=None)

    s = sub.add_parser("stages", help="print the stages, their order and what they fill")
    s.add_argument("--config", type=Path, default=None)

    args = ap.parse_args(argv)
    try:
        return {"run": cmd_run, "report": cmd_report, "stages": cmd_stages}[
            args.command](args)
    except (ConfigError, RegistryError) as e:
        print(f"hapax {args.command}: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
