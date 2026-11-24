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
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__, tools
from .config import Config, ConfigError, load
from .prepare import PrepareError, prepare
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
                stages: list[str], prepared=None,
                binaries: dict[str, Path] | None = None) -> list[str]:
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
    # THE TRIM AND THE CONTIG PILES COME FROM THE RUN, NOT FROM THE CALLER. This is the
    # whole of PLAN.md step 5: the driver required --trim-r2 and four paths, and the only
    # thing that knew how to produce them was 250 lines of shell, so the tool could not
    # be run on a new input by reading the repo. `prepare` performs those steps and this
    # passes what they produced.
    if prepared is not None:
        argv += ["--trim-r2", str(prepared.trim_r2),
                 "--assembly", str(prepared.assembly),
                 "--contigs", str(prepared.rrna), str(prepared.identified_norrna),
                 "--identified-contigs", str(prepared.identified_norrna),
                 "--rrna-contigs", str(prepared.rrna),
                 "--unidentified-contigs", str(prepared.unidentified),
                 # The sample's own identified contigs are the most efficient subtrahend
                 # measured here, +30.3 points, and they are the sample's own: nothing
                 # from another library reaches this flag.
                 "--cross-sample-ref", str(prepared.identified_norrna)]
        if prepared.removal_evidence:
            argv += ["--removal-evidence", str(prepared.removal_evidence)]
        block = cfg.fact("block") or (prepared.measured or {}).get("block")
        if block:
            argv += ["--block", str(block)]
    else:
        for key, flag in (("trim_r2", "--trim-r2"), ("block", "--block")):
            value = cfg.fact(key)
            if value is not None:
                argv += [flag, str(value)]
    if cfg.inputs.get("single_end"):
        argv.append("--single-end")
    for name, flag in (("rrna", "--silva"), ("classifier", "--k2-dbs"),
                       ("spacers", "--spacers"), ("pfam", "--pfam")):
        if name in cfg.databases:
            # `--k2-dbs` is the cascade's TIER LIST and its order is the tier order, so
            # `classifier` may name several databases. One database removes 42.8%, 21.7%
            # or 3.5% of the same library depending only on its confidence setting, so a
            # config able to name only one could not reproduce a committed tier table.
            argv += [flag, *[str(p) for p in cfg.database_list(name)]]
    # `--dmnd` is an ORDERED pair, viral then Swiss-Prot, and the distance stage refuses a
    # list of any other length: it reports which of the two produced the nearest relative,
    # and a single index cannot answer that. The order is the config's two keys rather
    # than a list under one, so neither can be given without saying which it is.
    protein = [cfg.database(n) for n in ("protein", "sprot") if n in cfg.databases]
    if protein:
        argv += ["--dmnd", *[str(p) for p in protein]]
    if cfg.run.get("report_without"):
        argv += ["--report-without", *cfg.run["report_without"]]
    # Resolved by the caller and always passed on, rather than left to each script's
    # default. Nine of them default `--classify-bin` to a path inside one container, which
    # is a machine fact in the place a parameter goes; see hapax.tools for why that is the
    # same fault as a trim reaching a library it was not measured on. It is resolved in
    # `cmd_run` rather than here so that this function stays a pure translation of the
    # config and the unit suite keeps running with no external tool installed (TOOLS.md).
    for name, flag in (("classify", "--classify-bin"), ("diamond", "--diamond-bin"),
                       ("minimap2", "--minimap2")):
        if (binaries or {}).get(name):
            argv += [flag, str(binaries[name])]
    if "pairs" in cfg.run:
        argv += ["--pairs", str(cfg.run["pairs"])]
    if "report_pairs" in cfg.run:
        argv += ["--report-pairs", str(cfg.run["report_pairs"])]
    # Only the stages the DRIVER has. `report` is the package's, and so are the three
    # enrichments adopted in step 6, which run after it; passing one of those on
    # `--stages` is an argparse choice error rather than a skipped stage, so the registry's
    # `in_driver` is what decides rather than a name-by-name exclusion that has to be kept
    # in step with the adoptions.
    driver_stages = [s for s in stages if R.stages[s].in_driver]
    argv += ["--stages", *driver_stages]
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
    out.mkdir(parents=True, exist_ok=True)
    work = Path(args.work).resolve() if args.work else out / "work"
    prepared = None
    if not args.skip_prepare:
        prepared = prepare(cfg, r1, r2, out, work, repo_root())
    # The external binaries, found once for the whole run: two tiers of one cascade run by
    # two builds of Kraken2 is an operating point nobody could write down.
    binaries = tools.resolve_all(cfg, ("classify", "diamond", "minimap2"))
    argv = legacy_argv(cfg, r1, r2, out, stages, prepared, binaries)
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
        "prepared": prepared.as_manifest() if prepared else None,
        "sample_facts": cfg.inputs.get("sample_facts"),
        # Which binary each name resolved to, recorded because the same tool name on two
        # machines is two builds, and a version is part of an operating point.
        "binaries": {k: str(v) for k, v in binaries.items()},
        "note": "the package translates one config into the legacy driver's flags; "
                "PLAN.md step 7 replaces the subprocess with imports",
    }, indent=1) + "\n")
    root = repo_root()
    proc = subprocess.run([sys.executable, str(root / "hapax.py"), *argv], cwd=root)
    if proc.returncode != 0:
        return proc.returncode
    # The run's manifest describes the whole run, and until step 7 the run is two
    # programs: the driver writes manifest.json and knows nothing of the steps that
    # produced its inputs. So the assembler's record is merged in here. Without it a
    # reader of the manifest cannot tell a run that assembled its own library from one
    # handed contigs by a caller, which is the difference PLAN.md step 5 is about.
    manifest = out / "manifest.json"
    if prepared is not None and manifest.exists():
        doc = json.loads(manifest.read_text())
        doc["assembly"] = prepared.assembly_record
        doc["prepared"] = prepared.as_manifest()
        manifest.write_text(json.dumps(doc, indent=1) + "\n")
    # The enrichment stages, which the driver does not have (PLAN.md step 6). They run
    # after it because each reads a file the driver's stages wrote, and their absences are
    # recorded rather than raised: the sections that read them render without them and say
    # what that costs.
    from .enrich import enrich

    for record in enrich(cfg, out, root, cfg.seed):
        if record["status"] != "ok":
            print(f"  {record['stage']}: {record['status']} - {record['reason']}",
                  file=sys.stderr)
    if not args.report:
        return 0
    # One command from reads to report (PLAN.md 4). The report refuses an incomplete run,
    # so a zero here means every declared section was drawn from a file this run wrote.
    from .report import ReportError, build

    try:
        return build(out, seed=cfg.seed, repo=root)
    except ReportError as e:
        print(f"hapax run: {e}", file=sys.stderr)
        return 1


def cmd_report(args) -> int:
    """Rebuild the page from a run directory alone, refusing an incomplete one.

    PLAN.md step 4. A run directory is now the whole input: the page descriptor is
    assembled from what the run recorded about itself, every declared stage output is
    checked against its schema before anything is written, and a section that cannot be
    drawn stops the build unless the caller asked for a partial page on purpose.
    """
    from .report import ReportError, build

    try:
        return build(Path(args.rundir), out=args.out, seed=args.seed,
                     allow_partial=args.allow_partial, repo=repo_root())
    except ReportError as e:
        print(f"hapax report: {e}", file=sys.stderr)
        return 1


def cmd_fetch_db(args) -> int:
    """Download the reference sets a run consults, and record what was downloaded.

    PLAN.md's definition of done opens with this command. Until it existed the four
    URLs and the two asserted byte counts lived only in an EC2 UserData script, so the
    first step of a run on new water was to read a shell heredoc.
    """
    from . import fetchdb

    keys = list(fetchdb.BY_KEY) if args.all else list(args.only or [])
    if not keys:
        print("hapax fetch-db: name a set with --only, or --all. Known sets:",
              file=sys.stderr)
        for refset in fetchdb.INVENTORY:
            print(f"  {refset.key:<11} {refset.what}", file=sys.stderr)
        for key, why in fetchdb.NOT_FETCHABLE.items():
            print(f"  {key:<11} NOT FETCHABLE: {why.split('.')[0]}.", file=sys.stderr)
        return 2
    # The binary that builds the protein indexes. A stated path is the override and is
    # used as it is; otherwise it is found on PATH, and left as None when it is not
    # there, so that fetching the sets that need no DIAMOND still works on a machine
    # that has none. The refusal then names only the set that actually needed it.
    diamond = args.diamond or shutil.which("diamond")
    try:
        manifest = fetchdb.fetch(keys, Path(args.dest), repo_root(),
                                 threads=args.threads,
                                 diamond=Path(diamond) if diamond else None,
                                 dry_run=args.dry_run)
    except fetchdb.FetchError as e:
        print(f"hapax fetch-db: {e}", file=sys.stderr)
        return 1
    if args.dry_run:
        for key, plan in manifest["sets"].items():
            print(f"{key}: {plan['what']}")
            print(f"  why      {plan['why']}")
            print(f"  fills    {plan['fills']}"
                  f"{'  (already present)' if plan['present'] else ''}")
            for url, pinned in zip(plan["urls"], plan["pinned_bytes"]):
                size = f"{pinned:,} bytes, asserted" if pinned else "size recorded"
                print(f"  {url}  [{size}]")
        return 0
    print(fetchdb.config_block(manifest))
    print(f"recorded in {Path(args.dest) / 'fetch_manifest.json'}", file=sys.stderr)
    return 0


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
    r.add_argument("--work", type=Path, default=None,
                   help="where the assembly and the indexes are built. Default: work/ "
                        "inside the run directory")
    r.add_argument("--skip-prepare", action="store_true",
                   help="do not measure, assemble or split; the config must then state "
                        "every value the stages need. For re-running the cascade over an "
                        "assembly this run directory already has")
    r.add_argument("--report", action="store_true",
                   help="build the page when the stages finish, which is what makes this "
                        "one command from reads to report")

    p = sub.add_parser("report", help="rebuild the page from a run directory")
    p.add_argument("rundir", type=Path)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--allow-partial", action="store_true",
                   help="build the page from an incomplete run. It is marked as "
                        "incomplete on the page itself, because a page that renders is "
                        "a page somebody reads")

    s = sub.add_parser("stages", help="print the stages, their order and what they fill")
    s.add_argument("--config", type=Path, default=None)

    f = sub.add_parser("fetch-db", help="download the reference sets a run consults")
    f.add_argument("--all", action="store_true", help="every fetchable set")
    f.add_argument("--only", nargs="+", default=None, metavar="SET",
                   help="named sets only, e.g. --only rrna protein")
    f.add_argument("--dest", type=Path, default=Path("data/db"),
                   help="where the sets are written. Default: data/db")
    f.add_argument("--threads", type=int, default=4)
    f.add_argument("--diamond", type=Path, default=None,
                   help="the diamond binary that builds the protein indexes. A .dmnd "
                        "built by one version is not readable by another, so this is "
                        "the same binary the run will search with")
    f.add_argument("--dry-run", action="store_true",
                   help="print the inventory, the URLs and which sizes are asserted, "
                        "and download nothing")

    args = ap.parse_args(argv)
    try:
        return {"run": cmd_run, "report": cmd_report, "stages": cmd_stages,
                "fetch-db": cmd_fetch_db}[args.command](args)
    except (ConfigError, RegistryError, PrepareError) as e:
        print(f"hapax {args.command}: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
