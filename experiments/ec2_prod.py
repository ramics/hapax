"""Launch a production end-to-end pass over one library on EC2 (`prod_userdata.sh`).

A copy of `ec2_crossplant.py` with the payload and template swapped; the spend gate, the
dirty-tree refusal, the bootstrap split and terminate-in-finally are unchanged, because
every one of them encodes a refusal that cost somebody a launch.

WHY A PRODUCTION RUN NEEDS AN INSTANCE AT ALL. The development container holds none of the
cascade's indexes and none of its binaries - no Kraken2, no MEGAHIT, no DIAMOND, no SILVA,
no reads - and cannot reach NCBI. Nothing about that is new; what is new is that this is
the first launcher whose instance stages the WHOLE index set rather than one search's, and
the first whose per-library facts are measured on the instance instead of read from a
committed entry, because for a library nobody has touched there is no entry to read.

WHY MEMORY RATHER THAN COMPUTE. The account ceiling is eight vCPU (`ec2_common.py`), so
the choice is only how much RAM sits behind them. A whole-library metatranscriptome
assembly is the memory-heavy step and PlusPF-16 is a 16 GB hash the cascade memory-maps
five times over, so `r6i.2xlarge` (64 GB) is the right buy at $0.50/h against the
compute-optimised type's $0.34.

THE SIZE DECISION, STATED. A full `hapax run` over a Berlin-sized library is tens of hours
and this library is 3.7x Berlin's size, so the sizing that matters is what the cascade
reads rather than what the assembler reads: five seeds of 200,000 pairs each, the
committed cascade sample size, at a stride sized to this library. The assembly is the whole
library, because the sample's own contigs are the largest subtrahend the cascade has
(+30.3 points) and an assembly of a fifth of the library is a different, weaker instrument.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ec2_common as ec2lib
from ec2_bigindex import AUTO_APPROVED_USD, projected_cost, render

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# THE STAGE LIST HAS ONE HOME. It was typed in three places - hapax.py's `--stages`
# default, prod_run.sh's `STAGES`, and here - and this launcher's copy is the one that
# actually ran, so the first production pass ran `panel` after that stage had been taken
# out of production everywhere else. Nothing broke, because the stage refuses on a library
# with no truth panel and said so; but a validation stage ran in a production pass because
# a list was typed a third time, which is the same family as a default where a refusal
# belongs. Imported now, so there is one list and the launcher cannot disagree with the
# tool it launches.
# Loaded from its path rather than by name, because `hapax` is the installed package now
# and `hapax.py` is the driver beside it: which one `import hapax` reaches would depend on
# sys.path order. The driver's content stays untouched, since its script_sha256 pins five
# committed seed manifests.
def _driver():
    import importlib.util
    path = Path(__file__).resolve().parent.parent / "hapax.py"
    spec = importlib.util.spec_from_file_location("hapax_legacy_driver", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PRODUCTION_STAGES = _driver().PRODUCTION_STAGES

#: Production stages this launcher's instance has no reference data for, each with the
#: reason, subtracted from the list above rather than a second list typed beside it.
#:
#: WHY THIS IS NOT A SHORTER STAGE LIST. Importing hapax's production list, correctly, put
#: `divergent` back into the pass - and it refused in 0.1 s on a missing `Pfam-A.hmm`,
#: which is the designed refusal and not a defect. What IS a defect is a launcher that
#: ships a run it knows will refuse: the instance stages no Pfam, so the stage can only
#: fail, and a run directory with a red stage in it every time teaches a reader to ignore
#: red stages. So the omission is stated here, by name, with what it would need - and
#: `--stages` still accepts it explicitly for a machine that has the index.
NO_REFERENCE_ON_THIS_INSTANCE = {
    "divergent": "experiments/divergent_known.py needs data/pfam/Pfam-A.hmm for its "
                 "domain-presence calibration, which this instance does not stage, and "
                 "its gene-level workup is hours of BLAST over the removed pile",
    "calibration": "the divergence ladder draws its spike genome from the RefSeq viral "
                   "nucleotide release, which this template does not stage and "
                   "`hapax fetch-db --all` does. It is a positive control and belongs in "
                   "every pass, so --oneshot runs it; here it could only refuse",
}
LAUNCHER_STAGES = [s for s in PRODUCTION_STAGES if s not in NO_REFERENCE_ON_THIS_INSTANCE]

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = Path(__file__).resolve().parent / "prod_userdata.sh"
#: The same launch, running the documented install and the documented command instead of
#: 450 lines of hand-rolled staging. Selected with `--oneshot`; see that template's header
#: for why it exists beside the other rather than replacing it.
ONESHOT_TEMPLATE = Path(__file__).resolve().parent / "prod_oneshot.sh"
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

#: micromamba, which creates the environment from the committed `environment.yml`. Pinned
#: for the reason every version here is pinned: the solve decides which build of MEGAHIT,
#: DIAMOND, Kraken2 and minimap2 a number was produced with, and a run's own log records
#: what it resolved to.
MICROMAMBA_VERSION = "2.0.8"

#: Pinned for the reason every version in this repo is pinned: a number produced by
#: "latest" cannot be reproduced by anyone reading it later.
SRA_VERSION = "3.1.1"
VDB_GUID = "00000000-0000-0000-0000-00000hapax01"
DIAMOND_VERSION = "2.1.10"            # TOOLS.md: the production search
KRAKEN2_VERSION = "2.17.2"            # TOOLS.md: the cascade calls its `classify` directly
MEGAHIT_VERSION = "1.2.9"             # TOOLS.md records it for every committed assembly
MINIMAP2_VERSION = "2.31"             # TOOLS.md: the same release as mappy 2.31

#: The three binaries no committed script has ever downloaded. Kraken2, MEGAHIT and
#: minimap2 were hand-placed in a staging directory, which is why `$S/k2pkg/root/share/...`
#: appears as a bare path in nine scripts with nothing in the repo that builds it. The
#: bioconda packages are the only pinned, verifiable source for the exact versions above -
#: upstream MEGAHIT ships no 1.2.9 linux tarball at a stable URL, and Kraken2's `classify`
#: only exists inside a package.
KRAKEN2_URL = ("https://conda.anaconda.org/bioconda/linux-64/"
               "kraken2-2.17.2-pl5321h3be2455_0.conda")
MEGAHIT_URL = ("https://conda.anaconda.org/bioconda/linux-64/"
               "megahit-1.2.9-haf24da9_8.tar.bz2")
MINIMAP2_URL = ("https://conda.anaconda.org/bioconda/linux-64/"
                "minimap2-2.31-h118bc1c_0.conda")

#: SILVA 138.1 Exports. The byte counts in `data/MANIFEST.md` are asserted on the instance:
#: a changed release is a changed operating point for the largest tier in the cascade, and
#: 138.2 is not present at that path.
SILVA_SSU_URL = ("https://ftp.arb-silva.de/release_138_1/Exports/"
                 "SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz")
SILVA_LSU_URL = ("https://ftp.arb-silva.de/release_138_1/Exports/"
                 "SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz")

#: The cascade's Kraken2 tiers, in the order and at the sizes `data/MANIFEST.md` records,
#: so this run's tier table is readable beside the committed one.
K2_DBS = ("k2_viral_20261022", "k2_standard_08_GB_20261022",
          "k2_pluspf_08_GB_20261022", "k2_pluspf_16_GB_20261022")
#: The contig split consults one database, and it is the largest: on contigs the thresholds
#: and not the size decide, and PlusPF-16 is the one the committed split used.
K2_DB_CONTIGS = "k2_pluspf_16_GB_20261022"

#: 64 GB behind eight vCPU; see the module docstring.
DEFAULT_INSTANCE_TYPE = "r6i.2xlarge"

#: GB, gp3. 33 GB of Kraken2 databases, ~11 GB of reads as gzip, ~20 GB of trimmed copies
#: for the assembler, MEGAHIT's own intermediates, and the two protein indexes. 300 GB is
#: several times the peak and costs $0.03/h.
DEFAULT_VOLUME_GB = 300
DEFAULT_VOLUME_THROUGHPUT = 250
DEFAULT_VOLUME_IOPS = 3000

#: Minutes. The dead-man timer, which is also the cost projection, so it is the number to
#: argue about: 11 hours is an assembly of a 9.7 Gbase library plus five cascade seeds with
#: room for one of them to be slower than measured, and nothing beyond it.
DEFAULT_HARD_STOP_MINUTES = 660
HEARTBEAT_SECONDS = 300

PAYLOAD_MODULES = ("hapax.py", "requirements.txt", "SPEC.md", "blind_spots.json",
                   "experiments/prod_run.sh",
                   # The two files that make `--oneshot`'s install the DOCUMENTED install
                   # rather than something this launcher assembled. The instance asserts
                   # both by name before it spends anything.
                   "pyproject.toml", "environment.yml", "README.md",
                   # The renderer reads this one from the repo, not from the run
                   # directory: the six fixed Methods rows, stored per accession, printed
                   # "Not measured on this library" for any other. Nothing under
                   # `results/` was in the payload, so the generator threw a
                   # FileNotFoundError after all 14 stages of seed 1 had finished, and
                   # `hapax run` correctly refused to write a page without its Methods
                   # section. `tests/test_payload.py` derives this from the renderer's
                   # own source rather than trusting the list.
                   "results/web/methods.json",
                   # Read the same way, one import away, and missed by the first version
                   # of that test because it read `make_web.py` alone: the control
                   # declaration every figure's control row is drawn from. Its absence
                   # cost the sixth launch of this library its page after every stage of
                   # seed 1 had run.
                   "results/web/controls.json")


def payload_files(root: Path, library: str | None = None) -> list[str]:
    """Everything the instance runs: the listed files, all of `experiments/`, the one
    module the current code still imports from the ARCHIVE, and the committed reference
    tables the selected stages read.

    The reference project twice shipped a payload missing a module and the instance died
    having read nothing, which is why `ec2_crossplant.py` lists its four files by hand.
    That list cannot be written for this run: `hapax run` dispatches to eleven scripts
    which import each other and read three JSON registries, and a missing one is a stage
    that refuses mid-run after the expensive half. So the unit is the directory, and the
    named files above are asserted on the instance one by one.

    THE DIRECTORY IS NOT THE WHOLE UNIT, MEASURED. The first production pass assembled
    for two hours and then refused every seed, on two files this rule does not catch:

    - `experiments/codon_assembler/experiments/common.py`. `readlevel_features.py` imports
      `common`, which lives only in the archived tree, and `residue_report.py` puts that
      tree on `sys.path` at import time to read one constant from it. The archive is
      excluded here to keep its results and its crate out of the payload - so the exclusion
      removed a module the live code imports. It is shipped by name, and the coupling is
      recorded in `experiments/SHARED.md` rather than hidden by a wider glob.
    - `results/crossplant/panel.tsv` and the rest of the panel the library resolves to.
      Those are committed RESULTS, read as reference data by a stage, and nothing under
      `results/` was in the payload at all.

    Both are the same shape as the rest of this file's lessons: the thing that was missing
    was the thing nobody listed, so the list is derived from the registry that names it.
    """
    files = list(PAYLOAD_MODULES)
    # THE PACKAGE, for the --oneshot template, whose whole install is `pip install .`.
    # Shipped as a directory for the reason the docstring above gives for `experiments/`:
    # a list of modules typed by hand is a list that loses one. `pyproject.toml` and
    # `environment.yml` are in PAYLOAD_MODULES and asserted on the instance by name,
    # because without either of them the install is not the documented install.
    for p in sorted((root / "src").rglob("*")):
        # `.egg-info` is build metadata left behind by an editable install in whatever
        # tree this launcher runs in. It is gitignored, so the dirty-tree refusal cannot
        # see it, and shipping it would put files in the payload that the commit does not
        # describe - the one thing the commit is there to guarantee.
        if p.is_dir() or "__pycache__" in p.parts or ".egg-info" in str(p):
            continue
        files.append(str(p.relative_to(root)))
    for p in sorted((root / "experiments").rglob("*")):
        if p.is_dir() or "codon_assembler" in p.parts or "__pycache__" in p.parts:
            continue
        if p.suffix in (".py", ".json", ".sh", ".md", ".toml", ".lock", ".rs"):
            files.append(str(p.relative_to(root)))
    files.extend(ARCHIVE_MODULES)
    files.extend(reference_files(root, library))
    return sorted(set(files))


#: Modules the live code imports from the archived assembler tree. Shipped by name so that
#: a payload cannot quietly lose one, and short enough to read: if it grows, the live code
#: has a dependency on an archive that should be lifted out of it instead.
ARCHIVE_MODULES = ("experiments/codon_assembler/experiments/common.py",)


def reference_files(root: Path, library: str | None) -> list[str]:
    """Committed files under `results/` that this library's stages read as reference data.

    Resolved from `experiments/stage_references.json`, so a panel that gains a file ships
    it without anybody editing this launcher. Only the production references are taken:
    `truth_panel` is validation and a production run does not select that stage.
    """
    if library is None:
        return []
    sys.path.insert(0, str(root / "experiments"))
    import stage_refs
    out = []
    r = stage_refs.resolve_optional(library, "crossplant_panel")
    for path in (r.get("files") or {}).values():
        rel = Path(path).resolve().relative_to(root)
        out.append(str(rel))
    return out


def build_payload(root: Path, files: list[str]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for rel in files:
            path = root / rel
            if not path.exists():
                raise SystemExit(f"{rel} is in the payload list and not on disk")
            tar.add(path, arcname=rel)
    return buf.getvalue()


def git_head(root: Path, shipped: list[str]) -> tuple[str, bool, bool]:
    """HEAD, whether the shipped files are dirty, and whether the wider tree is.

    Two flags rather than one, for the reason `ec2_bigindex.py` gives: what has to be clean
    is what the instance will run, and that is checkable exactly, while a session has
    result files landing in the tree the whole time.
    """
    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True,
                              text=True, check=True).stdout.strip()

    watched = [*shipped, "experiments/ec2_prod.py", "experiments/prod_userdata.sh",
               "experiments/prod_oneshot.sh", "experiments/ec2_common.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *watched)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


def _sample_facts_toml(sample: dict) -> str:
    """The library's `sample` block as TOML lines for `[inputs.sample_facts]`.

    Each fact carries the same `basis` from the entry, because the entry records one basis
    for the block rather than one per clause and inventing a per-clause basis here would
    be prose this code wrote. A library with no sample block produces a `basis` line
    saying so, which is what the page then prints: it states its ignorance rather than
    describing water nobody recorded.
    """
    if not sample:
        return ('basis = "the archive record carries no matrix, place, date or molecule '
                'for this library, so the page says what it does not know"')
    # ONE `basis` FOR THE BLOCK, and only keys the sentence has a clause for.
    # `page_descriptor.check` RAISES on a key it cannot speak, which is right: a fact
    # somebody wrote down that the renderer silently drops is the asymmetry that lost
    # "South Korea" off the one library not in Europe. So a per-clause `*_basis` here
    # would refuse the page rather than decorate it.
    import page_descriptor

    lines = []
    for key, _, _ in page_descriptor.CLAUSES:
        if sample.get(key):
            lines.append(f'{key} = "{str(sample[key]).replace(chr(34), chr(39))}"')
    basis = (sample.get("basis") or "").replace('"', "'").replace("\n", " ")
    lines.append(f'basis = "{basis}"')
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", required=True,
                    help="the run accession. It must have an entry in "
                         "experiments/libraries.json and one in "
                         "experiments/stage_references.json, because the stages refuse "
                         "rather than guess which panel is ground truth for it")
    ap.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE,
                    choices=ec2lib.ALLOWED_TYPES)
    ap.add_argument("--pairs", type=int, default=200_000,
                    help="pairs per cascade seed; the committed cascade sample size")
    ap.add_argument("--seeds", default="1 2 3 4 5")
    ap.add_argument("--report-pairs", type=int, default=2_000,
                    help="pairs through the per-read report. hapax defaults to 20,000; the "
                         "fold field alone costs 290 ms/read, so 20,000 is ~97 min per "
                         "seed. 2,000 is what results/holdout_cascade/ used")
    ap.add_argument("--prefix-seeds", default="1 2 3 4 5",
                    help="seeds for the mate-2 prefix measurement. More than one is not "
                         "optional: a mode that differs between seeds is not a length to "
                         "cut, and that is the condition trim_from_prefix.py refuses on")
    ap.add_argument("--prefix-pairs", type=int, default=100_000)
    ap.add_argument("--trim-r2", type=int, default=None,
                    help="state the mate-2 trim instead of measuring it. Stating a value "
                         "always wins and is never checked, because stating it IS the "
                         "override (SPEC 4.1); the measurement still runs and its JSON "
                         "still ships, so the two can be read against each other. Use it "
                         "only with a basis written into experiments/libraries.json - on a "
                         "library with little ribosomal content the SILVA instrument "
                         "cannot anchor enough reads to measure anything")
    ap.add_argument("--trim-r2-basis", default=None,
                    help="why that value, in a sentence, carried into the run config and "
                         "so onto the page. Required with --trim-r2, for the same reason "
                         "the config refuses a library fact with no basis: a stated trim "
                         "is a person's judgement, and a judgement nobody wrote down is "
                         "indistinguishable from another library's default")
    ap.add_argument("--stages", default=None,
                    help="hapax stages. The default template drops the stages it stages "
                         "no reference data for (NO_REFERENCE_ON_THIS_INSTANCE says which "
                         "and why); --oneshot keeps them, because `hapax fetch-db --all` "
                         "downloads the two indexes that were missing. `panel` is "
                         "validation and is not in the production list at all: a run "
                         "cannot create a depositor (SPEC 13.4)")
    ap.add_argument("--megahit-memory", type=float, default=0.85,
                    help="MEGAHIT --memory as a fraction of RAM")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--volume-size", type=int, default=DEFAULT_VOLUME_GB)
    ap.add_argument("--volume-throughput", type=int, default=DEFAULT_VOLUME_THROUGHPUT)
    ap.add_argument("--volume-iops", type=int, default=DEFAULT_VOLUME_IOPS)
    ap.add_argument("--hard-stop-minutes", type=int, default=DEFAULT_HARD_STOP_MINUTES)
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--no-spot", action="store_true", default=True,
                    help="on demand, and the default here: spot was reclaimed twice in "
                         "four hours, and this run has one unresumable step in the middle "
                         "of it, so a reclaim costs the assembly rather than one library")
    ap.add_argument("--spot", dest="no_spot", action="store_false",
                    help="pass only if losing the assembly mid-run is acceptable")
    ap.add_argument("--allow-price-change", action="store_true",
                    help="pass only with Imo's explicit yes on the figure printed below")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--timeout-minutes", type=float, default=700.0)
    ap.add_argument("--poll-seconds", type=float, default=300.0)
    ap.add_argument("--assembly-from", default=None,
                    help="s3:// URI of a gzipped, already-relabelled assembly of THIS "
                         "library to use instead of assembling. Two hours and about a "
                         "dollar of a 32 M-pair pass is MEGAHIT, and the first production "
                         "run spent them and then refused at the next step: the assembly "
                         "survived in the store, so the next pass reads it. It is used "
                         "as-is and no assembly manifest is written, so the run's "
                         "provenance for that step is the pass that produced it")
    ap.add_argument("--oneshot", action="store_true",
                    help="run `conda env create -f environment.yml && pip install .`, "
                         "`hapax fetch-db --all` and `hapax run --config --report` and "
                         "nothing else, which is what PLAN.md's definition of done is "
                         "judged on. The default template stages the indexes and the "
                         "binaries by hand and calls a shell driver, all of which is now "
                         "a command in the package")
    ap.add_argument("--report-without", nargs="*", default=None,
                    help="per-read report fields to leave out ON PURPOSE (--oneshot "
                         "only). Defaults to the fields whose reference is a property of "
                         "a sampling programme or a profile library this pass does not "
                         "stage. A field whose reference is absent and NOT named here "
                         "stops the run, because a dropped field reads as a zero")
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

    # WHICH STAGES, RESOLVED AFTER PARSING rather than as a default, because the answer
    # depends on the template. `divergent` was dropped from production launches for one
    # reason: the instance staged no Pfam-A and the stage refuses without it. `hapax
    # fetch-db --all` downloads it, so under --oneshot the exclusion no longer describes
    # anything and keeping it would drop a stage the report needs a section from.
    if args.stages is None:
        args.stages = " ".join(PRODUCTION_STAGES if args.oneshot else LAUNCHER_STAGES)
    if args.report_without is not None and not args.oneshot:
        raise SystemExit("--report-without is only read by the --oneshot template, which "
                         "writes it into the run configuration; the default template "
                         "takes it from prod_run.sh's own environment")
    if args.report_without is None:
        # Each of these is a field whose reference this pass does not have, named so that
        # the omission is a stated operating point rather than a zero. cross_sample needs
        # assemblies of OTHER samples from the same works, which is a property of a
        # sampling programme and which the sample's own assembly cannot stand in for;
        # rfam and pfam each need a profile library and, for rfam, Infernal besides;
        # crispr needs spacers called from this sample's assemblies.
        args.report_without = ["cross_sample", "rfam", "crispr", "pfam"]

    # A LIBRARY WITH NO ENTRY IS THE CASE THIS IS FOR, so a missing one is not a refusal.
    # The registries are still read here rather than on the instance, because where an
    # entry DOES exist its arithmetic can be checked before paying for anything; where it
    # does not, the instance measures the pair count, the read lengths, the layout and the
    # stride off the reads themselves (`measure_library.py`) and records them beside the
    # run. Only `crossplant_panel` is resolved: `truth_panel` is a validation reference and
    # a production run does not ask for one (SPEC 13.4).
    import library_facts
    import stage_refs
    try:
        facts = library_facts.library(args.library)
    except library_facts.UnknownLibrary:
        facts = {}
        print(f"{args.library} has no entry in experiments/libraries.json. Nothing here "
              f"needs one: the instance measures what the file can answer and the only "
              f"value a person supplies is --trim-r2.")
    block = facts.get("block")
    r = stage_refs.resolve_optional(args.library, "crossplant_panel")
    print("crossplant_panel: " + r["state"]
          + (f" ({r['panel']})" if r["state"] == "resolved" else f" ({r.get('reason', '')})"))
    if facts and facts.get("layout") not in (None, "paired"):
        raise SystemExit(f"{args.library} is {facts.get('layout')} in libraries.json; this "
                         f"template stages two mates and measures a mate-2 prefix")
    last = max(int(s) for s in args.seeds.split())
    if block is None:
        print(f"{args.library}: block stride measured on the instance from the counted "
              f"pairs; {last} seeds x {args.pairs:,} pairs")
    else:
        need = (last - 1) * block + args.pairs
        if facts.get("pairs") and need > facts["pairs"]:
            raise SystemExit(
                f"refusing to launch: seed {last} would start at {(last - 1) * block:,} in "
                f"a library of {facts['pairs']:,} pairs. The instance checks this again "
                f"against the counted reads, but there is no reason to pay for it.")
        print(f"{args.library}: {facts.get('pairs', 'unknown'):,} pairs, block {block:,}, "
              f"{last} seeds x {args.pairs:,} needs {need:,}")
    if args.trim_r2 is not None:
        print(f"trim_r2 STATED as {args.trim_r2}; libraries.json records "
              f"{facts.get('trim_r2')!r} with "
              f"{'an unverified flag' if facts.get('trim_r2_unverified') else 'its basis'}. "
              f"The measurement still runs and ships, and is not used.")
    elif facts.get("trim_r2") is None and not facts:
        print("this library has no entry at all, so the instance measures the trim over "
              f"{len(args.prefix_seeds.split())} seeds and refuses to run if no single "
              "mode holds, which is the intended path for a new library")
    elif facts.get("trim_r2") is None:
        print("trim_r2 is unmeasured in libraries.json: the instance measures it over "
              f"{len(args.prefix_seeds.split())} seeds and refuses to run if no single "
              "mode holds, which is the intended path for a new library")

    cost, cost_note = projected_cost(args.instance_type, args.hard_stop_minutes,
                                     not args.no_spot, args.volume_size)
    if args.attach:
        cost, cost_note = 0.0, f"watching {args.attach}; this launches nothing"
    print(f"projected worst case: {cost_note}")
    if not (cost == cost):
        raise SystemExit(cost_note)
    if cost > AUTO_APPROVED_USD and not args.allow_price_change:
        raise SystemExit(
            f"refusing to launch: {cost_note}, over the ${AUTO_APPROVED_USD:.0f} Imo "
            f"pre-approved for a single run. Put that figure in front of him and pass "
            f"--allow-price-change once he has answered, or shorten --hard-stop-minutes: "
            f"the projection is the timer's worst case, not the expected bill.")

    files = payload_files(ROOT, args.library)
    commit, shipped_dirty, tree_dirty = git_head(ROOT, files)
    if shipped_dirty and not args.allow_dirty:
        raise SystemExit(
            "refusing to launch: the code that would run on the instance is not committed. "
            "The instance has no checkout, so the commit is the only thing that makes its "
            "output traceable. Commit the launcher, the template and the payload first.")

    bucket, prefix = ec2lib.BUCKET, ec2lib.PREFIX
    code_key = f"{prefix}/code/prod-{args.library}-{commit}.tar.gz"
    tokens = dict(
        COMMIT=commit, BUCKET=bucket, PREFIX=prefix, REGION=ec2lib.REGION,
        HARD_STOP_MINUTES=args.hard_stop_minutes, HEARTBEAT_SECONDS=HEARTBEAT_SECONDS,
        CODE_KEY_NAME=Path(code_key).name, ACC=args.library,
        PAIRS=args.pairs, SEEDS=args.seeds,
        # Empty, not "None": prod_run.sh reads an empty BLOCK as "measure it",
        # and a string "None" would reach the arithmetic as a stated stride.
        BLOCK=("" if block is None else block),
        REPORT_PAIRS=args.report_pairs, STAGES=args.stages,
        PREFIX_SEEDS=args.prefix_seeds, PREFIX_PAIRS=args.prefix_pairs,
        TRIM_R2_STATED="" if args.trim_r2 is None else args.trim_r2,
        THREADS=args.threads, MEGAHIT_MEMORY=args.megahit_memory,
        SRA_VERSION=SRA_VERSION, VDB_GUID=VDB_GUID,
        DIAMOND_VERSION=DIAMOND_VERSION, KRAKEN2_VERSION=KRAKEN2_VERSION,
        KRAKEN2_URL=KRAKEN2_URL, MEGAHIT_URL=MEGAHIT_URL, MINIMAP2_URL=MINIMAP2_URL,
        SILVA_SSU_URL=SILVA_SSU_URL, SILVA_LSU_URL=SILVA_LSU_URL,
        K2_DBS=" ".join(K2_DBS), K2_DB_CONTIGS=K2_DB_CONTIGS,
        ASSEMBLY_FROM=(args.assembly_from or ""),
    )
    if args.oneshot:
        tokens.update(
            MICROMAMBA_VERSION=MICROMAMBA_VERSION,
            SEEDS_COUNT=len(args.seeds.split()),
            STAGES_TOML=", ".join(f'"{s}"' for s in args.stages.split()),
            REPORT_WITHOUT_TOML=", ".join(f'"{f}"' for f in args.report_without),
            # The page's one-line description of the water, from this library's own entry
            # rather than from anything this code knows. A clause whose fact is absent is
            # left out of the sentence instead of guessed, which is why this is built from
            # the keys that are there rather than from a fixed list.
            SAMPLE_FACTS=_sample_facts_toml(facts.get("sample") or {}),
            TRIM_R2_BASIS_STATED=(args.trim_r2_basis or ""),
        )
        if args.trim_r2 is not None and not args.trim_r2_basis:
            raise SystemExit(
                "refusing to launch: --trim-r2 states a library fact, and a stated fact "
                "with no basis is what experiments/libraries.json and the run config both "
                "refuse. Pass --trim-r2-basis saying why that value, in a sentence; it is "
                "carried into the run config and onto the page.")
    run_script = render(ONESHOT_TEMPLATE if args.oneshot else TEMPLATE, None, **tokens)
    # The two templates write to two different prefixes on purpose: the bucket has no
    # delete, and one of them is the committed producer of results/prod_SRR23702504/.
    # The launch manifest named only the one, so for a --oneshot run it recorded a
    # location nothing would ever be written to, and the watcher built from it reported a
    # live run as silent. Read against the rendered script rather than kept in step by
    # hand, because the pointer and the thing it points at were already out of step.
    seg = "oneshot" if args.oneshot else "prod"
    for want in (f"/logs/{seg}-", f"/results/{seg}/"):
        if want not in run_script:
            raise SystemExit(
                f"the rendered run script does not write under {want}, so the launch "
                f"manifest would record a prefix nothing uploads to")
    run_key = (f"{prefix}/code/prod-run-{commit[:12]}-"
               f"{hashlib.sha256(run_script.encode()).hexdigest()[:12]}.sh")
    script = render(BOOTSTRAP, 16_384, COMMIT=commit, BUCKET=bucket, REGION=ec2lib.REGION,
                    HARD_STOP_MINUTES=args.hard_stop_minutes, RUN_SCRIPT_KEY=run_key)

    if args.print_user_data:
        print(script)
        print(f"\n### ---- run script, fetched from s3://{bucket}/{run_key} ----\n")
        print(run_script)
        print(f"\n### {len(script.encode())} bytes of UserData, "
              f"{len(run_script.encode()):,} of run script, "
              f"{len(files)} payload files", file=sys.stderr)
        return 0

    import boto3

    s3 = ec2lib.s3_client()
    ec2 = boto3.client("ec2", region_name=ec2lib.REGION)
    print(f"commit {commit[:12]}, {len(files)} payload files", flush=True)

    if args.attach:
        return watch(ec2, args.attach, args)

    if not args.dry_run:
        s3.put_object(Bucket=bucket, Key=code_key, Body=build_payload(ROOT, files))
        print(f"payload at s3://{bucket}/{code_key}", flush=True)
        s3.put_object(Bucket=bucket, Key=run_key, Body=run_script.encode())
        print(f"run script at s3://{bucket}/{run_key}", flush=True)

    ssm = boto3.client("ssm", region_name=ec2lib.REGION)
    ami = ssm.get_parameter(Name=ec2lib.AMI_PARAMETER)["Parameter"]["Value"]
    market = ec2lib.market_options(not args.no_spot)

    kwargs: dict[str, Any] = dict(
        ImageId=ami, InstanceType=args.instance_type, MinCount=1, MaxCount=1,
        IamInstanceProfile={"Name": args.profile},
        InstanceInitiatedShutdownBehavior="terminate",
        BlockDeviceMappings=[{"DeviceName": "/dev/xvda",
                              "Ebs": {"VolumeSize": args.volume_size,
                                      "VolumeType": "gp3",
                                      "Iops": args.volume_iops,
                                      "Throughput": args.volume_throughput,
                                      "DeleteOnTermination": True}}],
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [{"Key": "Project", "Value": "ClaudeSpike"},
                                     {"Key": "Owner", "Value": args.owner},
                                     {"Key": "Task", "Value": f"hapax-prod-{args.library}"}]}],
        UserData=script,
    )
    if market is not None:
        kwargs["InstanceMarketOptions"] = market

    if args.dry_run:
        try:
            ec2.run_instances(DryRun=True, **kwargs)
        except Exception as exc:                      # DryRunOperation is the success case
            print(f"dry run: {exc}")
        return 0

    try:
        reply = ec2.run_instances(**kwargs)
    except Exception as exc:
        message = ec2lib.spot_refusal(str(exc)) if not args.no_spot else None
        if message:
            raise SystemExit(message) from exc
        raise
    iid = reply["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    record = json.dumps({
        "instance_id": iid, "commit": commit, "shipped_dirty": shipped_dirty,
        "wider_tree_dirty": tree_dirty, "instance_type": args.instance_type,
        "library": args.library, "pairs_per_seed": args.pairs, "seeds": args.seeds,
        "block": block, "report_pairs": args.report_pairs, "stages": args.stages,
        "prefix_seeds": args.prefix_seeds, "prefix_pairs": args.prefix_pairs,
        "trim_r2": (f"stated as {args.trim_r2} on the command line; the measurement runs "
                    f"and ships beside it but is not used"
                    if args.trim_r2 is not None else
                    "measured on the instance, and the run refuses if no single mode holds"),
        "volume_gb": args.volume_size, "spot": not args.no_spot,
        "hard_stop_minutes": args.hard_stop_minutes, "projected_usd": round(cost, 2),
        "sra_toolkit": SRA_VERSION, "diamond": DIAMOND_VERSION,
        "kraken2": KRAKEN2_VERSION, "megahit": MEGAHIT_VERSION,
        "minimap2": MINIMAP2_VERSION, "kraken2_dbs": list(K2_DBS),
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/{seg}/{args.library}/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/{seg}-{args.library}-{iid}/",
    }, indent=2)
    out = ROOT / f"results/prod_{args.library}"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"launch_{iid}.json").write_text(record + "\n")

    if args.detach:
        print(f"detached: {iid} is held by its own {args.hard_stop_minutes}-minute timer "
              f"and by terminate-on-shutdown.\n"
              f"  watch with --attach {iid}; terminate explicitly to stop it early.",
              flush=True)
        return 0

    try:
        return watch(ec2, iid, args)
    finally:
        try:
            ec2.terminate_instances(InstanceIds=[iid])
            print(f"terminate requested for {iid}", flush=True)
        except Exception as exc:
            print(f"terminate failed for {iid}: {exc}", file=sys.stderr)


def watch(ec2: Any, iid: str, args: Any) -> int:
    deadline = time.time() + args.timeout_minutes * 60
    while time.time() < deadline:
        state = ec2lib.state_of(ec2, iid)
        print(f"  {iid} {state} at {time.strftime('%H:%M:%SZ', time.gmtime())}", flush=True)
        if state in {"terminated", "stopped", "shutting-down"}:
            print(f"  reason: {ec2lib.transition_reason(ec2, iid)}", flush=True)
            return 0
        time.sleep(args.poll_seconds)
    print("timed out watching; the instance terminates itself on its own timer",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
