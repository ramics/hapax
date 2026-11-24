"""Measure the library, assemble it, and split the assembly by identification.

WHAT THIS REMOVES. `hapax run` started after the hardest three steps. A caller had to
supply `--contigs`, `--unidentified-contigs`, `--identified-contigs`, `--rrna-contigs`,
`--assembly`, `--cross-sample-ref` and a measured `--trim-r2`, and the only thing that knew
how to produce them was 250 lines of shell in `experiments/prod_run.sh` plus EC2 UserData.
So the tool could not be run on a new input by reading the repo, and the one number that
decides how much of every read is thrown away arrived as a flag somebody typed (PLAN.md
1.2, 1.9). This module performs those steps as declared stages of the run.

WHAT IT DOES NOT CHANGE. It drives the same scripts with the same arguments and the same
operating points as the shell did: the mate-2 trim and the MEGAHIT relabelling (lifted into
`hapax.assembly_input`, byte-equal to `trim_mate2.py` and `megahit_rename.py` by test),
`path_settings.py --stages assemble --arms baseline`, the SILVA index built from `subtract_known`'s own
two functions, and `classify_contigs.py`. Rewriting any of them here would be a second
instrument reading the same quantity, and the committed numbers were produced by these.
PLAN.md step 7 lifts their code into the package; this step moves the ORCHESTRATION, which
is the part that only existed as shell.

THE SPLIT IS BY IDENTIFICATION, NEVER BY ASSEMBLY. A contig that assembled and could not
be identified stays in the residue: subtracting unidentified contigs inflated removal to
85.9% against the true 77.0%. And the sample's own identified contigs are both the most
efficient subtrahend measured here, +30.3 points, and the trap: a divergent genome's own
contig that counts as identified is subtracted as background in tier 1 and never reaches
the part of the report anyone reads. That is what happened to the holdout positive control.
`removal_evidence.tsv` is what keeps such a removal from being silent about how far away
the thing it removed was.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import scripts, tools
from .assembly_input import rename_megahit, trim_mate2
from .measure import measure_to
from .config import Config, ConfigError


class PrepareError(Exception):
    """Raised rather than carrying on with an input one of these steps did not produce."""


@dataclass
class Prepared:
    """What the three steps produce, as the paths the cascade needs."""
    trim_r2: int
    trim_r2_source: str
    assembly: Path
    identified_norrna: Path
    rrna: Path
    unidentified: Path
    removal_evidence: Path | None
    measured: dict = field(default_factory=dict)
    #: Empty where the caller staged an assembly, and the assembler's own record of the
    #: run otherwise. `assembled_at_run_time` is the existing manifest convention, and
    #: what it is worth is visible in `results/holdout_cascade/`, whose assembly recorded
    #: it as false and pointed at a file that does not exist: only the contig names
    #: survive, so the command line behind a committed result cannot be recovered.
    assembly_record: dict = field(default_factory=dict)

    def as_manifest(self) -> dict:
        return {
            "assembly_record": self.assembly_record,
            "trim_r2": self.trim_r2,
            "trim_r2_source": self.trim_r2_source,
            "assembly": str(self.assembly),
            "identified_contigs_norrna": str(self.identified_norrna),
            "rrna_contigs": str(self.rrna),
            "unidentified_contigs": str(self.unidentified),
            "removal_evidence": (str(self.removal_evidence)
                                 if self.removal_evidence else None),
            "measured_off_the_reads": self.measured,
            "note": "produced by hapax run, not supplied by a caller. The split is by "
                    "identification and never by assembly: a contig that assembled and "
                    "could not be identified stays in the residue.",
        }


def _run(repo: Path, name: str, args: list, what: str, out: Path) -> str:
    """One declared script through `hapax.scripts`, its code record kept in the run.

    The record goes to `script_calls.jsonl` in the run directory, so a run that assembled
    its own library names the assembler driver's hash and every module it imports.
    """
    started = time.time()
    proc, code = scripts.run(repo, name, args)
    with open(out / "script_calls.jsonl", "a") as fh:
        fh.write(json.dumps({"step": what, **code}) + "\n")
    if proc.returncode != 0:
        raise PrepareError(
            f"{what} failed after {time.time() - started:.0f}s:\n"
            f"  {' '.join(code['argv'])}\n{proc.stdout[-2000:]}\n{proc.stderr[-3000:]}")
    return proc.stdout


def measure(cfg: Config, r1: Path, r2: Path | None, out: Path, repo: Path) -> dict:
    """Read the library's own facts off its reads.

    Imo, on being shown the ledger of what a production run had to be told by hand: "why
    are we not discovering the facts from the sample itself?" The pair count, the modal
    read length per mate, the layout and the seed stride are all in the file, and
    `measure_library.py` takes one pass over it for all four, refusing a stride that would
    make the five seeds overlap.

    The trim is NOT among them, and that is the honest part. Measuring a non-template 5'
    prefix needs reads anchored to a reference, and the instrument that does it
    (`nontemplate_prefix_map.py`, floor 2,000 anchored mates) refuses below that floor -
    correctly: at lower depth the histogram describes the reference's reach into the water
    rather than the library's chemistry. So the trim is resolved separately, below, and the
    run records which of the three possible sources it came from.
    """
    manifest = out / "library_measured.json"
    pairs_per_seed = int(cfg.run.get("pairs", 200_000))
    # HOW MANY SEEDS THE LIBRARY HAS TO HOLD, stated rather than assumed at five. The
    # stride is what makes the seeds read DISJOINT blocks, so it depends on the count, and
    # `measure_library.py` refuses a stride that would make them overlap - which is right:
    # five blocks of 1.5 M do not fit in a 4.29 M-read library, and that produced a
    # directory that looked like a five-seed result while seeds 4 and 5 read nothing. A
    # headline result still needs five (CLAUDE.md); a fixture is not a result, and stating
    # one here is how a 9,987-pair library runs at all instead of refusing.
    seeds = int(cfg.run.get("seeds", 5))
    try:
        return measure_to(r1, r2, pairs_per_seed, seeds, manifest,
                          accession=cfg.inputs.get("accession"))
    except ValueError as e:
        raise PrepareError(f"measuring the library off its reads failed: {e}") from e


def resolve_trim(cfg: Config, measured: dict) -> tuple[int, str]:
    """The mate-2 trim, and where it came from. Refuses rather than assuming one.

    Three sources, in order, and the run records which one answered:

      1. the config, which must carry `trim_r2_basis`. Stating it is the override.
      2. the library registry, for a library it knows.
      3. a refusal naming both.

    There is deliberately no fourth. A trim of 14 as a default deletes fourteen real bases
    from every mate of a library with no non-template prefix, silently, in the direction of
    less recognition, and it was written into roughly twenty scripts that way. A library
    whose layout is single-end is the one case with an answer that needs no measurement:
    there is no mate 2 to trim.
    """
    stated = cfg.fact("trim_r2")
    if stated is not None:
        return int(stated), f"stated in {cfg.path.name}: {cfg.library['trim_r2_basis']}"
    if (measured.get("layout") or cfg.inputs.get("single_end")) == "single" or \
            cfg.inputs.get("single_end"):
        return 0, "single-end: there is no mate 2 to trim"
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "experiments"))
    try:
        import library_facts

        value = library_facts.fact(cfg.label, "trim_r2")
        basis = library_facts.library(cfg.label).get("trim_r2_basis", "")
        return int(value), f"experiments/libraries.json for {cfg.label}: {basis}"
    except Exception:
        pass
    raise ConfigError(
        f"the mate-2 trim is not stated in {cfg.path.name} and {cfg.label} is not a "
        f"library experiments/libraries.json knows, so there is nothing to read it from. "
        f"It is a property of a sequencing run and this code has no right one to assume: "
        f"state it in [library] with its trim_r2_basis, or measure it with "
        f"experiments/nontemplate_prefix_map.py, which refuses below 2,000 anchored "
        f"mates rather than reading the reference's reach into the water as chemistry.")


def assemble(cfg: Config, r1: Path, r2: Path | None, trim_r2: int, work: Path,
             out: Path, repo: Path) -> Path:
    """Assemble the whole library, at the committed baseline setting.

    Through `path_settings.py --arms baseline` rather than a bare megahit line, because
    that script is the one place in this repo that records the full argv, the MEGAHIT
    version and the timing of an assembly, and `results/path_settings/` measured the whole
    ladder of alternatives against this arm. Baseline is MEGAHIT's defaults plus
    `--min-contig-len 200`.

    An assembly already staged by the caller is used as it is, and that is not an
    optimisation: a 32 M-pair library is two hours of the pass, and the first production
    run spent them, refused at the next step, and lost the assembly with the instance.
    """
    staged = cfg.inputs.get("assembly")
    asm_dir = work / "asm_baseline"
    contigs = asm_dir / "sample_contigs.fa"
    if contigs.exists() and contigs.stat().st_size:
        return contigs
    if staged:
        asm_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(Path(staged).expanduser(), contigs)
        return contigs

    # Read-level stages trim as they stream; MEGAHIT reads files, so it needs a trimmed
    # copy. The copy is an intermediate and nothing reads it afterwards.
    asm_in = work / "asm_in"
    try:
        trim_mate2(r1, r2 or r1, trim_r2, asm_in, cfg.label,
                   manifest=out / "trim_manifest.json")
    except ValueError as e:
        raise PrepareError(f"trimming mate 2 for the assembler failed: {e}") from e

    found = tools.resolve_all(cfg, ("megahit", "minimap2"))
    megahit, minimap2 = found["megahit"], found["minimap2"]
    _run(repo, "path_settings.py", ["--seed", 1,
          "--stages", "assemble", "--arms", "baseline",
          "--reads", asm_in / f"{cfg.label}_1.fq", "--r2", asm_in / f"{cfg.label}_2.fq",
          "--megahit", megahit, "--minimap2", minimap2,
          "--threads", cfg.threads, "--no-truth-panel",
          "--work", work, "--out", out / "assembly"],
         "assembling the library with MEGAHIT", out)

    final = asm_dir / "final.contigs.fa"
    if not (final.exists() and final.stat().st_size):
        raise PrepareError(f"the assembler wrote no contigs at {final}")
    # MEGAHIT's header into the SPAdes-like form the committed scripts parse. NAMED
    # sample_contigs.fa deliberately: crossplant_stage.py's queries arm requires
    # sample_contigs_unidentified.fa and sample_contigs_identified.fa in its contig
    # directory BY NAME, so calling the relabelled assembly anything else makes a stage
    # that would otherwise run refuse for a reason that has nothing to do with the sample.
    try:
        rename_megahit(final, contigs, cfg.label, manifest=out / "rename_manifest.json")
    except ValueError as e:
        raise PrepareError(f"relabelling the assembly failed: {e}") from e
    return contigs


def split(cfg: Config, contigs: Path, work: Path, out: Path, repo: Path) -> dict:
    """Identify the contigs, and split them by that and never by assembly.

    `classify_contigs.py` consults the classifier, the protein index and SILVA and writes
    three piles plus `removal_evidence.tsv`. The SILVA index it reads is built from
    `subtract_known`'s own two functions rather than a copy of them: a second SILVA cleaner
    or a different minimap2 preset would be a second instrument reading one quantity.
    """
    silva_mmi = work / "silva.mmi"
    if not silva_mmi.exists():
        rrna = cfg.database("rrna")
        if rrna.suffix == ".mmi":
            silva_mmi = rrna
        else:
            _build_silva_index(rrna, work, cfg.threads, repo)

    result = out / "contig_identification.json"
    if not result.exists():
        _run(repo, "classify_contigs.py", ["--seed", 1,
              "--contigs", contigs,
              # The single contig-identification database, not the cascade's tier list:
              # `classify_contigs.py` reads one `hash.k2d` and writes that database's name
              # into every row it produces. See CONFIG_KEYS in config.py.
              "--k2-db", cfg.database("classifier_contigs"),
              "--dmnd", cfg.database("protein"),
              "--silva-mmi", silva_mmi,
              "--threads", cfg.threads, "--work", work / "contigs", "--out", out,
              # Always passed, never left to the script's default, which is the layout of
              # one container written where a parameter goes. See hapax.tools.
              "--classify-bin", tools.resolve(cfg, "classify"),
              "--diamond-bin", tools.resolve(cfg, "diamond")],
             "identifying the contigs", out)

    asm_dir = contigs.parent
    piles = {"identified_norrna": asm_dir / "sample_contigs_identified_norrna.fa",
             "rrna": asm_dir / "sample_contigs_rrna.fa",
             "unidentified": asm_dir / "sample_contigs_unidentified.fa"}
    absent = [str(p) for p in piles.values() if not p.exists()]
    if absent:
        raise PrepareError(
            f"the identification step wrote no {', '.join(absent)}. The three piles are "
            f"the cascade's whole input, and an absent pile would be subtracted as an "
            f"empty one, which reads as a library with nothing identifiable in it.")
    return piles


def _build_silva_index(rrna: Path, work: Path, threads: int, repo: Path) -> Path:
    """The rRNA index, built by the cascade's own functions.

    Imported rather than reimplemented, for the reason the shell gave when it did the
    same: two cleaners or two minimap2 presets are two instruments reading one quantity.
    """
    sys.path.insert(0, str(repo / "experiments"))
    import subtract_known as sk

    work.mkdir(parents=True, exist_ok=True)
    fa = work / "silva_dna.fa"
    if not fa.exists():
        sk._silva_to_dna([rrna], fa)
    mmi = work / "silva.mmi"
    sk._mm_index(fa, mmi, threads)
    return mmi


def prepare(cfg: Config, r1: Path, r2: Path | None, out: Path, work: Path,
            repo: Path) -> Prepared:
    """The three steps, in order, writing what they produced into the run directory."""
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    measured = measure(cfg, r1, r2, out, repo)
    trim_r2, source = resolve_trim(cfg, measured)
    contigs = assemble(cfg, r1, r2, trim_r2, work, out, repo)
    # The assembler's own record, carried into the run's manifest rather than left in a
    # work directory: whether the run built its own contigs or was handed them is the
    # difference between one command from reads and a pipeline with a missing first half,
    # and a reader cannot tell them apart from the contigs.
    record = out / "assembly" / "assembly_baseline.json"
    assembly_record = {"assembled_at_run_time": False,
                       "why": "the configuration named an assembly under [inputs], which "
                              "is used as it is"} if cfg.inputs.get("assembly") else {}
    if record.exists():
        assembly_record = {"assembled_at_run_time": True,
                           **json.loads(record.read_text())}
    piles = split(cfg, contigs, work, out, repo)
    evidence = out / "removal_evidence.tsv"
    prepared = Prepared(
        trim_r2=trim_r2, trim_r2_source=source, assembly=contigs,
        identified_norrna=piles["identified_norrna"], rrna=piles["rrna"],
        unidentified=piles["unidentified"],
        removal_evidence=evidence if evidence.exists() else None,
        measured=measured, assembly_record=assembly_record)
    (out / "prepared.json").write_text(json.dumps(prepared.as_manifest(), indent=1) + "\n")
    return prepared
