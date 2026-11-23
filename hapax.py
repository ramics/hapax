"""`hapax run` - the production tool's command line, per SPEC.md 8.1.

WHAT THIS IS. An orchestrator, deliberately thin. Every measurement it reports already
exists as a script in `experiments/`, each with its own seed, manifest and calibration,
because that is how this project establishes that a number is real. The production tool is
those scripts run in a fixed order into one directory with one manifest - not a
reimplementation of them. A parallel implementation would be a second instrument, and this
project's own rule is never to compare numbers from two instruments.

WHAT IT ADDS, which is the part worth having:

  1. **One run directory with one provenance record.** Each stage's own manifest is read
     back and its `script_sha256` recorded, so the run says which code produced every
     file in it. A `-dirty` tree is recorded as dirty rather than hidden: it names HEAD,
     not what ran, and two runs can carry the same SHA and different code.
  2. **Enforcement of the three interface rules** from SPEC.md 8.3, at write time rather
     than in review. No tier row leaves this script without its operating point; no pile
     fraction leaves it without an interval; no field leaves it without its calibration or
     an explicit statement that it has none. `cannot discriminate here` and `not measured`
     are values, written into the file, never blanks.
  3. **The declared blind spots travel with the run**, copied from the committed
     `blind_spots.json` with its hash, so a run directory cannot be read without them.

WHY THE RUN REFUSES RATHER THAN DEGRADES. Before any stage runs, every selected stage is
preflighted: hapax checks its own inputs, and each experiment script is called with
`--preflight` and the exact argv it will get, so the script's own list of required indexes
and tools is the one checked. Anything absent, empty or not executable stops the run in its
first seconds with every problem listed, exit status 2, and a manifest saying it refused.
A stage is `skipped` only when the command line asked for that (`--stages`, no
`--calibrate`, no `--unidentified-contigs`); an input that vanishes mid-run is `failed`.
Neither falls back to a cheaper index - the project has already had to withdraw a figure
whose index silently differed from the one named.

WHAT WAS CONSULTED, NOT WHAT WAS ASKED FOR. hapax validated three indexes and never passed
them on, so runs used the defaults while the manifest named the indexes requested; the
existence check made it harder to see, because it proved the files were there. So every
path a stage is asked to use is asserted (a) to appear in the argv actually built and (b)
where the script writes a `consulted` record, to appear there as opened. A mismatch fails
the stage.

THE RESIDUE IS LAYERED, NOT SCORED. `layers.json` orders the residue from most to least
worth a human's attention using only fields whose direction has been measured on known
material: whether a relative exists at all and how far, whether anything recurs, whether
the read's own mate places, whether a profile matches, and whether the pair is artefact-
shaped. A layer is a sort order with its evidence attached, and emphatically not a novelty
score: an index can say what is already known and cannot say what is new.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "experiments"))

ROOT = Path(__file__).resolve().parent
PY = sys.executable


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return ""


def wilson(k: int, n: int, z: float = 1.96) -> dict | None:
    """Wilson score interval for a count out of a total.

    Every pile in this project ships with an interval, and a single-seed run has no
    across-seed spread to quote, so the interval has to come from the counts. Wilson rather
    than normal-approximation because the interesting fractions here are small and a few of
    them are zero, where the normal interval is a point.
    """
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {"point": round(p, 6), "lo": round(max(0.0, centre - half), 6),
            "hi": round(min(1.0, centre + half), 6), "k": k, "n": n,
            "method": "wilson_95"}


class Run:
    def __init__(self, args):
        self.args = args
        self.dir: Path = args.out
        self.stages: list[dict] = []
        self.t0 = time.perf_counter()
        (self.dir / "stages").mkdir(parents=True, exist_ok=True)
        self.work: Path = args.work or (self.dir / "work")
        self.work.mkdir(parents=True, exist_ok=True)
        #: True while every stage is being checked before any runs: `call` then runs each
        #: script with --preflight, and `missing` ignores what an earlier stage will write.
        self.preflighting = False
        self.problems: list[str] = []

    # -- stage plumbing -------------------------------------------------------
    def record(self, name: str, status: str, **kw) -> dict:
        if self.preflighting:
            if status == "failed":
                self.problems.append(f"{name}: {kw.get('reason')}")
            return {"stage": name, "status": status, **kw}
        entry = {"stage": name, "status": status, **kw}
        self.stages.append(entry)
        print(f"[{name}] {status}" + (f": {kw.get('reason')}" if kw.get("reason") else ""),
              flush=True)
        return entry

    def missing(self, name: str, paths: list[Path]) -> list[str]:
        """Inputs are checked by existence before a stage runs, never assumed.

        A stage that starts without its reference data fails deep inside a subprocess,
        often hours in, and on a machine that cannot be debugged. `test -f` first is the
        cheapest rule this project keeps. An EMPTY file is as absent as a missing one: a
        zero-byte index is a failed download, not an index with nothing in it.

        During preflight, a path under the work directory is an earlier stage's output and
        is not expected to exist yet, provided that stage is selected.
        """
        out = []
        for p in map(Path, paths):
            if self.preflighting and "cascade" in self.args.stages and \
                    self.work in p.parents:
                continue
            if not p.exists():
                out.append(str(p))
            elif p.is_file() and p.stat().st_size == 0:
                out.append(f"{p} (0 bytes)")
        return out

    def reference(self, library: str, kind: str) -> dict:
        """Resolve a stage's per-library reference, keeping the three answers apart.

        `resolved` runs the stage; `absent` is a MEASURED absence and the stage is skipped
        with the committed reason; an unknown library is a `failed` stage and, during
        preflight, a problem that stops the run. Those three were one condition before, and
        collapsing them is how a run comes to report ground truth as checked when no panel
        for this water exists - or scores it against another sample's genomes.
        """
        import stage_refs
        try:
            r = stage_refs.resolve(library, kind)
        except (stage_refs.UnknownLibrary, stage_refs.MissingReference) as e:
            return {"state": "unknown", "status": "failed", "reason": str(e)}
        if r["state"] == "absent":
            return {"state": "absent", "status": "skipped", "reason": r["reason"]}
        gone = stage_refs.missing_files(r)
        if gone:
            return {"state": "incomplete", "status": "failed",
                    "reason": f"{kind} {r['panel']}: {gone}"}
        return {"state": "resolved", "status": "ok", "panel": r["panel"],
                "record": r["record"], "files": r["files"]}

    def call(self, name: str, script: str, argv: list[str],
             expect: dict[str, Path], manifest: Path | None,
             requested: list[Path] | None = None) -> dict:
        """Run one experiment script and adopt its outputs into the run directory.

        `requested` is every index or reference this stage was asked to use. Each must
        appear in `argv` - the check that would have caught the cascade being validated
        against three indexes and then run without them - and, where the script records
        what it opened, in that record.
        """
        # A script's built-in default for a sample-specific input is one particular
        # sample's file (Berlin's assembly, Berlin's CRISPR arrays, Berlin's trim). On a
        # machine where that file exists the stage runs happily against the wrong sample
        # and nothing in the output says so, which is worse than a missing file. So hapax
        # must pass every such flag itself; one absent from the argv is a refusal.
        # AND IT MUST PASS THE LIBRARY ITSELF, TO EVERY SCRIPT THAT TAKES ONE. `hapax run`
        # requires `--library` and then handed it to two of the thirteen stages. The other
        # scripts resolve their per-library facts through `library_facts`, whose `--library`
        # defaults to the discriminating dataset so that every committed command reproduces
        # unchanged - so a stage that is not told which library it is looking at asks for
        # Berlin's facts. Measured, on the first production library: the `linkage` stage
        # needed a block stride, was not given one, and refused because the read path named
        # SRR23702504 while `--library` said SRR23973193. That refusal is the path check
        # doing its job one layer in; on a library whose reads sit under a neutral path it
        # would have used Berlin's stride silently, which is the original bug.
        if script in LIBRARY_AWARE and "--library" not in argv:
            argv = [*argv, "--library", self.args.library]
        without = set(getattr(self.args, "report_without", []) or [])
        defaulted = [flag for flag, field in SAMPLE_SPECIFIC.get(script, [])
                     if flag not in argv and not (field and field in without)]
        if defaulted:
            return self.record(name, "failed", script=script,
                               reason=f"{defaulted} not passed, so {script} would use its "
                                      f"built-in default, which is one specific sample's "
                                      f"file or setting rather than this run's")
        dropped = [str(p) for p in (requested or []) if str(p) not in argv]
        if dropped:
            return self.record(name, "failed", script=script,
                               reason=f"asked to use {dropped} but they are not in the "
                                      f"argv this stage would be run with, so the script "
                                      f"would fall back to its own defaults")
        if self.preflighting:
            return self._preflight_script(name, script, argv)
        t = time.perf_counter()
        log = self.dir / "stages" / f"{name}.log"
        cmd = [PY, "-I", str(ROOT / "experiments" / script), *argv]
        with open(log, "w") as lg:
            lg.write("$ " + " ".join(cmd) + "\n\n")
            lg.flush()
            proc = subprocess.run(cmd, cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT)
        seconds = round(time.perf_counter() - t, 1)
        if proc.returncode != 0:
            tail = log.read_text().splitlines()[-12:]
            return self.record(name, "failed", script=script, seconds=seconds,
                               argv=argv, log=str(log.relative_to(self.dir)),
                               reason=f"exit {proc.returncode}", tail=tail)
        # A stage's manifest is its provenance, and a stage without one has produced
        # numbers nobody can attribute to code. That used to be recorded as an empty dict.
        if manifest is None or not manifest.exists():
            return self.record(name, "failed", script=script, seconds=seconds, argv=argv,
                               reason=f"exit 0 but its manifest {manifest} was not written")
        try:
            m = json.loads(manifest.read_text())
        except ValueError as e:
            return self.record(name, "failed", script=script, seconds=seconds, argv=argv,
                               reason=f"its manifest {manifest} is not JSON: {e}")
        prov = {k: m[k] for k in ("git_sha", "tree_dirty", "script_sha256") if k in m}
        if "script_sha256" not in prov:
            return self.record(name, "failed", script=script, seconds=seconds, argv=argv,
                               reason=f"its manifest carries no script_sha256, so the code "
                                      f"that produced this stage cannot be pinned")
        consulted = m.get("consulted")
        if consulted is not None and requested:
            seen = {e.get("path") for e in consulted}
            never = [str(p) for p in requested if str(p) not in seen]
            if never:
                return self.record(name, "failed", script=script, seconds=seconds,
                                   argv=argv, consulted=consulted,
                                   reason=f"asked to use {never}; the stage's own record "
                                          f"of what it opened does not include them")
        outputs = {}
        for label, src in expect.items():
            if not src.exists():
                return self.record(name, "failed", script=script, seconds=seconds,
                                   argv=argv, reason=f"{src} was not written")
            if src.stat().st_size == 0:
                return self.record(name, "failed", script=script, seconds=seconds,
                                   argv=argv, reason=f"{src} was written empty")
            dst = self.dir / label
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            outputs[label] = {"bytes": dst.stat().st_size, "sha256_16": sha256(dst)[:16],
                              "from": str(src)}
        extra = {"consulted": consulted} if consulted is not None else {
            "consulted": "this script does not record what it opened; every path in argv "
                         "was checked to exist before the run"}
        return self.record(name, "ok", script=script, seconds=seconds, argv=argv,
                           provenance=prov, outputs=outputs, **extra,
                           log=str(log.relative_to(self.dir)))

    def _preflight_script(self, name: str, script: str, argv: list[str]) -> dict:
        """Ask the script itself whether it could start with this argv.

        The script owns the list of what it needs, so hapax does not keep a second copy
        that drifts. Scripts that have no --preflight are checked by hapax alone.
        """
        if script not in PREFLIGHT_SCRIPTS:
            return {"stage": name, "status": "ok"}
        cmd = [PY, "-I", str(ROOT / "experiments" / script), *argv, "--preflight"]
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        if proc.returncode != 0:
            lines = [ln.strip() for ln in (proc.stderr or proc.stdout).splitlines()
                     if ln.strip()][-15:]
            self.problems.append(f"{name} ({script} --preflight, exit {proc.returncode}):"
                                 + "".join(f"\n      {ln}" for ln in lines))
        return {"stage": name, "status": "ok"}

    def preflight(self) -> None:
        """Check every selected stage before any of them runs, and refuse if one cannot.

        A run used to discover a missing index when it reached the stage that needed it,
        which for the per-read report is after the cascade - tens of minutes in, on a
        machine that may have been launched for that run alone.
        """
        a = self.args
        for p, what in ((a.r1, "r1"), (a.r2, "r2")):
            if self.missing("inputs", [p]):
                self.problems.append(f"inputs: {what} {p} is absent or empty")
        self.preflighting = True
        try:
            for name in STAGE_ORDER:
                if name in a.stages and name in SCRIPT_STAGES:
                    getattr(self, name)()
        finally:
            self.preflighting = False
        if a.unidentified_contigs and "contigs" in a.stages:
            for p in self.missing("contigs", a.unidentified_contigs):
                self.problems.append(f"contigs: {p} is absent or empty")

    # -- stages ---------------------------------------------------------------
    def cascade(self) -> dict:
        """Tier 1: remove everything identifiable at all, and record what each index cost.

        This is the subtrahend and it is IDENTIFICATION, not assembly. Its own operating
        points matter more than its database sizes: the same Kraken2 database removes 42.8%
        or 3.5% of this library depending only on confidence.
        """
        a = self.args
        if a.contigs is None:
            return self.record("cascade", "failed",
                               reason="--contigs not given. It has no default because any "
                                      "default is one sample's own assembly; give this "
                                      "sample's identified contigs, or `--contigs` with "
                                      "no paths to run the cascade without them")
        need = [a.r1, a.r2, *a.k2_dbs, *a.silva, *a.dmnd, *(a.contigs or [])]
        gone = self.missing("cascade", need)
        if gone:
            return self.record("cascade", "failed", reason=f"absent inputs: {gone}")
        out = self.work / "cascade"
        argv = ["--seed", str(a.seed), "--pairs", str(a.pairs),
                "--r1", str(a.r1), "--r2", str(a.r2), "--trim-r2", str(a.trim_r2),
                "--threads", str(a.threads), "--rows-out",
                "--work", str(self.work / "cascade_work"), "--out", str(out)]
        # The block stride WAS NOT EXPOSED HERE AT ALL, so `hapax run` could only ever use
        # subtract_known.py's default of 1,500,000 - the stride for SRR23973193's 8,731,224
        # pairs. On a smaller library the later seeds skip past the end of the FASTQ and the
        # cascade exits with "no pairs read". Measured: on SRR23421544's 4,290,429 reads,
        # seeds 4 and 5 start at 4.5 M and 6.0 M and read nothing, so a five-seed run
        # produced three seeds and two failures while the run directory still looked like
        # five. A parameter that decides WHICH READS ARE MEASURED does not belong behind a
        # default tuned to one library.
        if a.block:
            argv += ["--block", str(a.block)]
        # These three were VALIDATED above and then never passed, so the cascade silently
        # fell back to subtract_known.py's own default paths whatever the run was given.
        # On a machine where the defaults are absent the stage dies; on one where they
        # exist it would quietly measure against the wrong indexes and report the run's
        # arguments in the manifest - the index named in the provenance not being the one
        # consulted, which is the failure this repo already had to withdraw a figure for.
        if a.k2_dbs:
            argv += ["--k2-dbs", *[str(d) for d in a.k2_dbs]]
        if a.silva:
            argv += ["--silva", *[str(d) for d in a.silva]]
        if a.dmnd:
            argv += ["--dmnd", *[str(d) for d in a.dmnd]]
        if a.classify_bin:
            argv += ["--classify-bin", str(a.classify_bin)]
        if a.diamond_bin:
            argv += ["--diamond-bin", str(a.diamond_bin)]
        if a.k2_no_memory_map:
            argv += ["--no-k2-memory-map"]
        # Always passed, even empty: an absent --selfasm makes subtract_known.py fall back
        # to Berlin's assembly.
        argv += ["--selfasm", *[str(c) for c in a.contigs]]
        if a.unidentified_contigs:
            argv += ["--observe", *[str(c) for c in a.unidentified_contigs]]
        return self.call("cascade", "subtract_known.py", argv,
                         {"tiers.json": out / f"subtract_s{a.seed}.json",
                          "stages/cascade_manifest.json": out / f"manifest_s{a.seed}.json",
                          "stages/cascade_rows.jsonl":
                              out / f"residue_rows_s{a.seed}.jsonl"},
                         out / f"manifest_s{a.seed}.json",
                         requested=[*a.k2_dbs, *a.silva, *a.dmnd, *(a.contigs or []),
                                    *(a.unidentified_contigs or [])])

    def reads(self) -> dict:
        """Tier 3: the per-read report, un-thresholded, every field beside its calibration."""
        a = self.args
        res = self.work / "cascade_work" / f"residue_s{a.seed}.fa"
        sub = self.work / "cascade_work" / f"subtracted_s{a.seed}.fa"
        gone = self.missing("reads", [res, sub, a.r1, a.r2])
        if gone:
            return self.record("reads", "failed",
                               reason=f"the cascade's residue and its subtracted pile are "
                                      f"the inputs, and {gone} is absent or empty. An "
                                      f"empty subtracted pile leaves every field without "
                                      f"its known-positive control")
        out = self.work / "reads"
        argv = ["--seed", str(a.seed), "--residue", str(res), "--reads", str(res),
                "--max-pairs", str(a.report_pairs), "--threads", str(a.threads),
                "--library-r1", str(a.r1), "--library-r2", str(a.r2),
                "--work", str(self.work / "reads_work"), "--out", str(out)]
        # The subtracted pile is the known-positive control for every field, so it is
        # always passed, and its absence fails the stage rather than silently leaving the
        # report with no control.
        argv += ["--subtracted-reads", str(sub), "--trim-r2", str(a.trim_r2)]
        # The references that belong to ONE sample. Each is this run's or the field is
        # left out by name; the script's own defaults are Berlin's files.
        without = set(a.report_without or [])
        for field, flag, value, what in (
                ("chimaera", "--chimaera-ref", a.assembly,
                 "--assembly (this sample's own full assembly)"),
                ("cross_sample", "--cross-sample-ref", a.cross_sample_ref,
                 "--cross-sample-ref (assemblies of other samples from the same works)"),
                ("crispr", "--spacers", a.spacers,
                 "--spacers (CRISPR spacers from this sample's assemblies)")):
            if field in without:
                continue
            if value is None:
                return self.record("reads", "failed",
                                   reason=f"the {field} field needs {what}, which has no "
                                          f"default because any default is another "
                                          f"sample's file; give it, or pass "
                                          f"--report-without {field}")
            argv += [flag, str(value)]
        # The protein indexes were NOT passed through before, so the report fell back to
        # its own defaults. Where those are absent every `closest_*` field is missing, and
        # the `strains` stage - which reads exactly those fields - writes an empty file
        # that looks like "no divergent relatives" rather than "nothing was searched".
        if a.dmnd:
            argv += ["--dmnd", *[str(d) for d in a.dmnd]]
        if a.diamond_bin:
            argv += ["--diamond-bin", str(a.diamond_bin)]
        if a.single_end:
            argv += ["--single-end"]
        if a.report_without:
            # Fields left out on purpose are named, and the per-read report records each
            # as "not searched: --without" - never as a zero.
            argv += ["--without", *a.report_without]
        return self.call("reads", "residue_report.py", argv,
                         {"reads.jsonl": out / f"rows_s{a.seed}.jsonl",
                          "fields.json": out / f"calibration_s{a.seed}.json"},
                         out / f"manifest_s{a.seed}.json",
                         requested=list(a.dmnd) + [
                             v for f, v in (("chimaera", a.assembly),
                                            ("cross_sample", a.cross_sample_ref))
                             if v is not None and f not in (a.report_without or [])])

    #: Why these two stages refuse on a single-end library rather than running: both are
    #: mate-linkage instruments whose ground truth is "these two sequences came from one
    #: insert". With no second mate the only way to run them is to make mate 2 a copy of
    #: mate 1, and then `linkage` reports near-perfect friendship and `census` labels every
    #: residue read's neighbourhood with its own copy. Those are FABRICATED positives, not
    #: degraded measurements, and a degraded measurement is the only kind this tool is
    #: allowed to report. So they are skipped with the reason in the manifest, and the
    #: merge, mate disagreement and mate linkage of SPEC 4.3 are absent from the run.
    SINGLE_END_SKIP = ("the library has no second mate. This is a mate-linkage instrument, "
                       "and duplicating mate 1 would report a fabricated positive rather "
                       "than a degraded measurement, so it is not run. SPEC 4.3's merge, "
                       "mate disagreement and mate linkage are therefore not measured in "
                       "this run")

    def linkage(self) -> dict:
        """The grouping calibration: can a read find a read it is genuinely related to?

        Mate linkage is the only label-free ground truth for "same molecule" available in a
        real library, so it is what the grouping step is calibrated on. It is an UPPER bound
        on friendship resolution - two mates share an insert - and the planted-genome purity
        in results/residue_structure/ is the lower bound. Both belong in the run directory
        or the grouping in the web view has no operating point.
        """
        a = self.args
        if a.single_end:
            return self.record("linkage", "skipped", reason=self.SINGLE_END_SKIP)
        tool = ROOT / "experiments/kmer_coverage/target/release/kmer_coverage"
        gone = self.missing("linkage", [tool, a.r1, a.r2])
        if gone:
            return self.record("linkage", "failed",
                               reason=f"absent inputs: {gone}. cargo build --release in "
                                      f"experiments/kmer_coverage builds the counter")
        out = self.work / "linkage"
        argv = ["--seed", str(a.seed), "--pairs", str(a.linkage_pairs),
                "--r1", str(a.r1), "--r2", str(a.r2), "--trim-r2", str(a.trim_r2),
                "--work", str(self.work / "linkage_work"), "--out", str(out)]
        return self.call("linkage", "mate_friends.py", argv,
                         {"linkage.json": out / f"mate_friends_s{a.seed}.json"},
                         out / f"mate_friends_s{a.seed}.json")

    def census(self) -> dict:
        """What the residue sits NEXT TO, from the adjacency a paired library gives free.

        Coverage says how often a read's sequence recurs and cannot say what it is attached
        to. This finds the library pairs covering the same sequence as a residue read and
        labels THEIR other mates, which splits the residue into reads that are unidentified
        because the reference is a gene catalogue and reads whose neighbourhood is as
        unknown as they are. The control pile goes through the same pass, and a query with
        no witnesses is reported apart rather than folded into the unknown.
        """
        a = self.args
        if a.single_end:
            return self.record("census", "skipped", reason=self.SINGLE_END_SKIP)
        # All three label classes are this sample's own contigs. The script's defaults
        # are Berlin's, so on any other library every neighbour would have been labelled
        # against the wrong assembly.
        absent = [f for f, v in (("--rrna-contigs", a.rrna_contigs),
                                 ("--identified-contigs", a.identified_contigs),
                                 ("--unidentified-contigs", a.unidentified_contigs))
                  if not v]
        if absent:
            return self.record("census", "failed",
                               reason=f"the census labels neighbours against this "
                                      f"sample's own contigs and needs {absent}")
        res = self.work / "cascade_work" / f"residue_s{a.seed}.fa"
        binary = ROOT / "experiments/kmer_coverage/target/release/census"
        gone = self.missing("census", [binary, res, a.r1, a.r2])
        if gone:
            return self.record("census", "failed", reason=f"absent inputs: {gone}")
        out = self.work / "census"
        argv = ["--seed", str(a.seed), "--queries", str(a.census_queries),
                "--cascade-work", str(self.work / "cascade_work"),
                "--r1", str(a.r1), "--r2", str(a.r2), "--trim-r2", str(a.trim_r2),
                "--threads", str(a.threads),
                "--work", str(self.work / "census_work"), "--out", str(out)]
        # The census labels a neighbour against the contigs that assembled without being
        # identified, so it takes the first of these sets only.
        argv += ["--unidentified-contigs", str(a.unidentified_contigs[0]),
                 "--rrna-contigs", str(a.rrna_contigs),
                 "--identified-contigs", str(a.identified_contigs)]
        return self.call("census", "linkage_census.py", argv,
                         {"census.json": out / f"linkage_census_s{a.seed}.json"},
                         out / f"linkage_census_s{a.seed}.json")

    def calibration(self) -> dict:
        """Per-sample divergence ladder: a holdout genome spiked at known divergence.

        Imo's requirement, and the only positive control this tool has. Without it every
        sensitivity statement in the run is an assertion: the ladder is what converts
        "survived subtraction" into "would have survived at 20% divergence".
        """
        a = self.args
        if not a.calibrate:
            return self.record("calibration", "skipped",
                               reason="--calibrate not given. The run's divergence ladder "
                                      "is therefore absent rather than assumed, and no "
                                      "sensitivity-at-divergence figure may be read off it")
        # The ladder's tier 1 has to be THIS run's tier 1: the same Kraken2 databases and
        # protein indexes as the cascade, and this sample's own contigs. None of them was
        # passed before, so the ladder measured Berlin's subtrahend for any library. The
        # spike genome is sample-specific too: it is a clean holdout only if it is absent
        # from this sample's assembly.
        absent = [f for f, v in (("--spike-genome", a.spike_genome),
                                 ("--contigs", a.contigs),
                                 ("--unidentified-contigs", a.unidentified_contigs),
                                 ("--cross-sample-ref", a.cross_sample_ref)) if not v]
        if absent:
            return self.record("calibration", "failed",
                               reason=f"the ladder needs this sample's {absent}")
        silva = self.work / "cascade_work" / "silva.fa"
        gone = self.missing("calibration", [a.r1, a.r2, a.spike_genome, silva, *a.contigs])
        if gone:
            return self.record("calibration", "failed", reason=f"absent inputs: {gone}")
        # The ladder takes one identified-contig file and the cascade subtracts several,
        # so the cascade's set is concatenated: the ladder's tier 1 is then the cascade's.
        ident = self.work / "calibration_identified.fa"
        if not self.preflighting:
            with open(ident, "w") as h:
                for c in a.contigs:
                    h.write(Path(c).read_text())
        out = self.work / "calibration"
        argv = ["--seed", str(a.seed), "--out", str(out),
                "--genome", str(a.spike_genome),
                "--r1", str(a.r1), "--r2", str(a.r2), "--threads", str(a.threads),
                "--trim-r2", str(a.trim_r2),
                "--silva", str(silva),
                "--identified", str(ident),
                "--unidentified", str(a.unidentified_contigs[0]),
                "--cross-sample-ref", str(a.cross_sample_ref),
                "--k2-dbs", *[str(d) for d in a.k2_dbs],
                "--dmnd", *[str(d) for d in a.dmnd],
                "--work", str(self.work / "calibration_work")]
        if a.classify_bin:
            argv += ["--classify-bin", str(a.classify_bin)]
        if a.diamond_bin:
            argv += ["--diamond-bin", str(a.diamond_bin)]
        return self.call("calibration", "spike_holdout.py", argv,
                         {"calibration_ladder.json": out / f"spike_s{a.seed}.json",
                          "stages/calibration_manifest.json":
                              out / f"manifest_s{a.seed}.json"},
                         out / f"manifest_s{a.seed}.json")

    # -- assembled-but-unidentifiable ----------------------------------------
    def contigs(self) -> dict:
        """Tier 2: what assembled but could not be identified.

        The assembly is an input rather than a step here, and its provenance is recorded as
        such. A production run assembles at run time; this version refuses to pretend it
        did, because a contig set of unknown origin in a run directory is exactly the
        unresolved provenance this project has a rule against.
        """
        a = self.args
        if not a.unidentified_contigs:
            return self.record("contigs", "skipped",
                               reason="--unidentified-contigs not given, so tier 2 has no "
                                      "population in this run")
        p = Path(a.unidentified_contigs[0])
        if self.missing("contigs", [p]):
            return self.record("contigs", "failed", reason=f"{p} is absent or empty")
        rows = self.dir / "contigs.jsonl"
        n = 0
        # ONE HASH, NOT ONE PER ROW. Every row carries the source file's digest, and the
        # digest was taken inside the per-contig writer, so this re-read the whole file
        # once per contig: on the first production library's 104 MB of unidentifiable
        # contigs, 194,703 contigs x 104 MB is about 20 TB of SHA-256 for a stage whose
        # own arithmetic is a length and a GC fraction. Measured: the stage had run 1h45
        # and was roughly two thirds through, against seconds for the parse, and it is a
        # PER-SEED cost, so five seeds would have spent most of a day on it. The file does
        # not change while it is being read, so the digest is taken once here.
        source_sha = sha256(p)[:16]
        with open(p) as fh, open(rows, "w") as out:
            name, buf = None, []

            def flush():
                nonlocal n
                if name is None:
                    return
                seq = "".join(buf).upper()
                if not seq:
                    return
                acgt = sum(1 for c in seq if c in "ACGT") or 1
                out.write(json.dumps({
                    "contig": name, "length": len(seq),
                    "gc": round(sum(1 for c in seq if c in "GC") / acgt, 4),
                    "n_fraction": round(1 - acgt / len(seq), 4),
                    "class": "assembled_unidentifiable",
                    "operating_point": {
                        "source": str(p), "source_sha256_16": source_sha,
                        "note": "class is the cascade's own label at its 40-base minimum "
                                "match length; no classifier was run on the contig here"},
                }) + "\n")
                n += 1

            for line in fh:
                if line.startswith(">"):
                    flush()
                    name, buf = line[1:].strip().split()[0], []
                else:
                    buf.append(line.strip())
            flush()
        if n == 0:
            return self.record("contigs", "failed",
                               reason=f"{p} holds no sequence, so tier 2 would read as an "
                                      f"empty population rather than an unreadable file")
        return self.record("contigs", "ok", rows=n,
                           outputs={"contigs.jsonl": {
                               "bytes": rows.stat().st_size, "rows": n,
                               "sha256_16": sha256(rows)[:16]}},
                           assembly={"source": str(p), "sha256_16": source_sha,
                                     "assembled_at_run_time": False,
                                     "note": "supplied contig set; see assembly.json"})

    # -- the distance a known thing was removed at ---------------------------
    def divergent(self) -> dict:
        """Every identified contig, with the distance at which it was identified.

        WHY THIS STAGE EXISTS AND WHY IT IS NOT `strains`. `strains` groups the RESIDUE by
        the subject it is nearest to, which finds a divergent relative only when the
        cascade failed to remove it. The case it cannot see is the opposite one: a contig
        that IS removed, correctly, by a tier that recognises it, from a long way away.
        The cascade's subtrahend is one bit, so a contig matching its reference at 99.9%
        and one matching it at 66.5% record the same word and neither is reported again.

        That is not hypothetical. A withheld genome assembles out of its library as a
        single contig covering 96% of itself, hits the held-out reference at 66.5%
        protein identity over 620 aa with a bitscore of 840 against the protein tier's
        threshold of 50, and is removed in tier 1. Not a marginal removal, so **raising
        the tier threshold cannot be the fix**; the fix is to keep the distance.

        The stage is skipped rather than approximated when the identified contig set is
        not given, like every other stage here.
        """
        a = self.args
        src = getattr(a, "identified_contigs", None)
        if not src:
            return self.record("divergent", "skipped",
                               reason="--identified-contigs not given, so there is no "
                                      "identified pile to measure the distance over")
        src = Path(src)
        gone = self.missing("divergent", [src, a.r1, a.r2])
        if gone:
            return self.record("divergent", "failed", reason=f"absent inputs: {gone}")
        # The script takes its two protein indexes as separate flags, and was given
        # NEITHER, so it used its own defaults whatever --dmnd said - the same dropped
        # argument as the cascade's. hapax's --dmnd is (viral, Swiss-Prot) in that order.
        if len(a.dmnd) != 2:
            return self.record("divergent", "failed",
                               reason=f"needs exactly two --dmnd indexes (viral, then "
                                      f"Swiss-Prot); got {[str(d) for d in a.dmnd]}")
        out = self.work / "divergent"
        if not self.preflighting:
            out.mkdir(parents=True, exist_ok=True)
        argv = ["--seed", str(a.seed), "--contigs", str(src),
                "--r1", str(a.r1), "--r2", str(a.r2),
                "--trim-r2", str(a.trim_r2),
                "--dmnd", str(a.dmnd[0]), "--dmnd-sprot", str(a.dmnd[1]),
                "--work", str(out / "work"),
                "--out", str(out / f"divergent_known_s{a.seed}.json")]
        if a.diamond_bin:
            argv += ["--diamond-bin", str(a.diamond_bin)]
        return self.call("divergent", "divergent_known.py", argv,
                         {f"divergent_known_s{a.seed}.json":
                              out / f"divergent_known_s{a.seed}.json"},
                         out / f"divergent_known_s{a.seed}.json",
                         requested=[src, *a.dmnd])

    # -- what the removal pile hides -----------------------------------------
    def shortlist(self) -> dict:
        """Rank the pile tier 1 removed by how much of each contig its best hit explains.

        WHY THIS IS A STAGE. The one measured mechanism that catches a genuinely novel
        genome being subtracted as known. On the only positive this project has - a virus
        deposited after the library that contains it - identity ordering puts it at the
        93.46th percentile of the removed pile, so `flag the low-identity removals` reaches
        7,104 other contigs first; adding the fraction of the contig the best hit explains
        turns the same pile into a cell of 69 contigs, 0.91% of the removals, which no
        recognisable contig from the same water enters. Until now that lived in a report
        over one directory and a new sample got none of it.

        IT READS `removal_evidence.tsv`, which the cascade already writes and then nobody
        reads. So the stage is skipped, not approximated, when the contig-identification
        step that writes that file was not run in this run.
        """
        a = self.args
        src = getattr(a, "removal_evidence", None)
        if not src:
            return self.record("shortlist", "skipped",
                               reason="--removal-evidence not given, so there is no "
                                      "recorded removal evidence to rank. It is written by "
                                      "experiments/classify_contigs.py beside the contig "
                                      "split this run's --contigs came from")
        src = Path(src)
        gone = self.missing("shortlist", [src])
        if gone:
            return self.record("shortlist", "failed", reason=f"absent inputs: {gone}")
        out = self.work / "shortlist"
        argv = ["--seed", str(a.seed), "--evidence", str(src), "--out", str(out)]
        if a.shuffled_removal_evidence:
            argv += ["--shuffled-evidence", str(a.shuffled_removal_evidence)]
        else:
            # The null is dropped BY NAME or not at all: a band quoted with a silently
            # absent null is the shape this project has withdrawn figures over.
            argv += ["--without-null"]
        return self.call("shortlist", "removal_shortlist.py", argv,
                         {"removal_shortlist.json": out / "removal_shortlist.json",
                          "removal_shortlist.tsv": out / "removal_shortlist.tsv",
                          "stages/removal_shortlist_manifest.json":
                              out / "removal_shortlist_manifest.json"},
                         out / "removal_shortlist_manifest.json", requested=[src])

    # -- ground truth from the same water ------------------------------------
    def panel(self) -> dict:
        """Score this run's assembly against the depositors' own genomes from this water.

        WHY EVERY RUN GETS IT. The residue fraction - the number this tool reports about
        itself - moves 8.79 points on assembler flags that do not improve recovery of a
        single real genome, so a residue figure without a recovery figure from the same pass
        is unreadable: a setting that shrinks the residue by fusing novel sequence into an
        abundant known contig reads as an improvement on every summary number here. The
        recovery half cannot be computed from the residue, which is why it is a stage and
        not a sweep.

        A LIBRARY WITH NO PANEL IS SKIPPED WITH ITS REASON, which is Berlin's case: no
        depositor genome set exists for that water, so its positive control is the spiked
        ladder (`--calibrate`). That is a different fact from a panel that was scored and
        found nothing, and the stage's output says which.
        """
        a = self.args
        # Resolved HERE as well as inside the script, so a library with no panel is a
        # skip recorded without starting a subprocess, and an unknown library is a problem
        # the preflight lists beside the others rather than an exit status to interpret.
        ref = self.reference(a.library, "truth_panel")
        if ref.get("state") != "resolved":
            return self.record("panel", ref["status"], reason=ref["reason"])
        argv = ["--seed", str(a.seed), "--library", a.library,
                "--work", str(self.work / "panel_work"),
                "--out", str(self.work / "panel")]
        if a.assembly:
            argv += ["--contigs", str(a.assembly)]
            if not a.minimap2:
                return self.record("panel", "failed",
                                   reason="--assembly given but --minimap2 is not, so the "
                                          "panel could be composed and nothing scored "
                                          "against it")
            argv += ["--minimap2", str(a.minimap2), "--threads", str(a.threads)]
        out = self.work / "panel"
        return self.call("panel", "truth_panel.py", argv,
                         {"truth_panel.json": out / "truth_panel.json",
                          "stages/truth_panel_manifest.json":
                              out / "truth_panel_manifest.json"},
                         out / "truth_panel_manifest.json",
                         requested=[Path(a.assembly)] if a.assembly else None)

    # -- other plants' water -------------------------------------------------
    def crossplant(self) -> dict:
        """Subtract a panel of other plants' libraries, with that panel's rules as checks.

        WHAT IT ENFORCES, because these are what get lost when a measurement is reused:
        query reads from OUTSIDE the library searched (a library inside its own panel is a
        refusal, not a caveat), paired per plant rather than pooled, contig queries in the
        panel's own 90 nt window so the two arms stay one instrument, and both operating
        points travelling together because they disagree by a factor of four.

        THE PANEL IS REUSED, NOT REBUILT. It resolves from the committed registry with its
        digest; what is staged in S3 and what each pass still re-downloads are both recorded,
        so a fresh download cannot read as reuse.
        """
        a = self.args
        ref = self.reference(a.library, "crossplant_panel")
        if ref.get("state") != "resolved":
            return self.record("crossplant", ref["status"], reason=ref["reason"])
        argv = ["--seed", str(a.seed), "--library", a.library,
                "--work", str(self.work / "crossplant_work"),
                "--out", str(self.work / "crossplant")]
        stages = []
        if a.unidentified_contigs and a.identified_contigs:
            stages.append("queries")
            argv += ["--contig-dir", str(Path(a.unidentified_contigs[0]).parent)]
        if a.crossplant_run_dir:
            stages.append("report")
            argv += ["--run-dir", str(a.crossplant_run_dir), "--arm", a.crossplant_arm]
        if not stages:
            return self.record("crossplant", "skipped",
                               reason="neither a contig set to cut into panel-width query "
                                      "windows (--unidentified-contigs with "
                                      "--identified-contigs) nor a completed pass to "
                                      "report (--crossplant-run-dir) was given")
        argv += ["--stages", *stages]
        out = self.work / "crossplant"
        return self.call("crossplant", "crossplant_stage.py", argv,
                         {"crossplant_stage.json": out / "crossplant_stage.json",
                          "stages/crossplant_stage_manifest.json":
                              out / "crossplant_stage_manifest.json"},
                         out / "crossplant_stage_manifest.json")

    # -- new strains of known viruses ----------------------------------------
    def strains(self) -> dict:
        """Group the residue by the relative it is nearest to, which is where a new strain
        of a known virus shows up.

        THE SHAPE OF THE ANSWER. A new strain is not a read with a low score; it is a
        SUBJECT that several residue reads point at, at an identity too low for the
        nucleotide tiers to have removed them, covering enough of the read to mean
        something. So the unit here is the subject, not the read, and the three numbers that
        make it readable travel together: how many residue reads point at it, how far they
        are from it, and how much of each read the alignment covers.

        WHY THIS IS NOT A NOVELTY SCORE. Ordering subjects by identity would order them by
        database proximity, which this project has measured repeatedly and refuses to do.
        What is reported is a distance with its coverage and its E-value, per subject, and
        the already-identified pile's count on the same subject beside it. A subject that
        the identified pile also points at is a virus that is present and partly removed;
        one that only the residue points at is a relative of something in the database and
        nothing more until it assembles.

        EVERY SUBJECT CARRIES ITS E-VALUE DISTRIBUTION. A 107 nt query against a large
        index finds about 40 bases at about 88% by chance, and the chance band reads
        identically on residue and identified reads, so a subject built out of chance hits
        would look exactly like a divergent strain without one.
        """
        p = self.dir / "reads.jsonl"
        if not p.exists():
            return self.record("strains", "skipped" if "reads" not in self.args.stages
                               else "failed", reason="reads.jsonl is absent")
        # WHETHER ANYTHING WAS SEARCHED decides what an empty table means, and it is read
        # from the per-read report's own record of what it opened - never inferred from
        # the rows, where "searched, no hit" and "never searched" look identical (a
        # missing `closest_*` key). This stage once wrote an empty strains.json that read
        # as "no divergent relatives found" when the protein index had never been passed.
        fp = self.dir / "fields.json"
        calib = json.loads(fp.read_text()) if fp.exists() else {}
        cons = [e for e in calib.get("consulted", []) or []
                if str(e.get("role", "")).startswith("closest_")]
        res_searched = [e for e in cons if e.get("state") == "searched"
                        and e["role"].endswith(" residue")]
        known_cons = [e for e in cons if e.get("state") == "searched"
                      and e["role"].endswith(" known")]
        if not res_searched:
            searched = {"searched": False, "indexes": [],
                        "why": ("no --dmnd given, so no read was searched for a relative"
                                if not self.args.dmnd else
                                "the per-read report records no protein search, although "
                                f"{[str(d) for d in self.args.dmnd]} was requested")}
            (self.dir / "strains.json").write_text(json.dumps({
                "state": "not searched", **searched, "subjects": None,
                "subjects_total": None,
                "reading": "NOT a negative result: nothing was searched, so this run says "
                           "nothing about divergent relatives"}, indent=2))
            if self.args.dmnd:
                return self.record("strains", "failed", reason=searched["why"])
            return self.record("strains", "skipped", reason=searched["why"])
        by: dict[str, dict] = {}
        for line in open(p):
            r = json.loads(line)
            pile = "residue" if r.get("set", "residue") == "residue" else "known"
            best = None
            for k, v in r.items():
                if k.startswith("closest_") and isinstance(v, dict) \
                        and v.get("evalue") is not None:
                    if best is None or v["evalue"] < best["evalue"]:
                        best = {**v, "index": k[len("closest_"):]}
            if best is None or best["evalue"] >= SIGNIFICANCE_EVALUE:
                continue
            s = by.setdefault(best["subject"], {
                "subject": best["subject"], "index": best["index"],
                "residue_reads": 0, "known_reads": 0,
                "pident": [], "aln_len": [], "evalue": [], "bitscore": []})
            s[f"{pile}_reads"] += 1
            if pile == "residue":
                for f in ("pident", "aln_len", "evalue", "bitscore"):
                    if best.get(f) is not None:
                        s[f].append(best[f])
        # Whether the control pile was searched at all decides what a zero in its column
        # means, so it is counted rather than assumed. The same check saved the layer
        # table from reporting a separation nobody had measured. It used to count known
        # rows WITH a hit, so a control searched with no hits read as "not searched".
        known_searched = max((e["queries"] for e in known_cons), default=0)
        rows = []
        for s in by.values():
            if not s["pident"]:
                continue
            def med(xs, sig=False):
                v = sorted(xs)[len(xs) // 2]
                return float(f"{v:.3g}") if sig else round(v, 4)
            rows.append({
                "subject": s["subject"], "index": s["index"],
                "residue_reads": s["residue_reads"], "known_reads": s["known_reads"],
                "median_identity_percent": med(s["pident"]),
                "identity_range": [min(s["pident"]), max(s["pident"])],
                "median_aligned_aa": med(s["aln_len"]),
                "median_evalue": med(s["evalue"], sig=True),
                "worst_evalue": float(f"{max(s['evalue']):.3g}"),
                "reading": ("present and partly identified: the already-subtracted pile "
                            "points at this subject too" if s["known_reads"] else
                            "only the residue points here, which is a relative in the "
                            "database and not an identification" if known_searched else
                            "not comparable in this run: the already-subtracted pile was "
                            "not searched, so its zero is absence and not a measurement"),
            })
        rows.sort(key=lambda d: (-d["residue_reads"], d["median_identity_percent"]))
        out = {
            "state": "searched",
            "searched": True,
            "indexes": [{"path": e["path"], "queries": e["queries"], "hits": e["hits"]}
                        for e in res_searched],
            "significance_evalue": SIGNIFICANCE_EVALUE,
            "subjects": rows[:500],
            "subjects_total": len(rows),
            "known_pile_searched_reads": known_searched,
            "unit": "one row per subject, not per read",
            "not_a_score": "subjects are not ordered by novelty. Identity is a distance "
                           "from a database and score ordering tracks database proximity, "
                           "so the ordering here is by how many residue reads point at the "
                           "subject",
            "caveat": "a 107 nt query against a large index finds about 40 bases at about "
                      "88% by chance and the chance band reads identically on residue and "
                      "identified reads, so every row carries its worst E-value as well as "
                      "its median",
        }
        (self.dir / "strains.json").write_text(json.dumps(out, indent=2))
        multi = [r for r in rows if r["residue_reads"] >= 3]
        return self.record("strains", "ok", subjects=len(rows),
                           subjects_with_3_or_more_reads=len(multi))

    # -- derived views --------------------------------------------------------
    def intervals(self) -> dict:
        """Attach a Wilson interval and an operating point to every tier row.

        The rule is SPEC.md 8.3's: no pile without its interval. A single-seed run cannot
        quote across-seed spread, so the interval comes from the counts, and the file says
        which it is.
        """
        p = self.dir / "tiers.json"
        if not p.exists():
            return self.record("intervals", "skipped" if "cascade" not in self.args.stages
                               else "failed", reason="tiers.json is absent")
        d = json.loads(p.read_text())
        n = d.get("pairs") or d.get("pairs_scored") or 0
        if not n:
            return self.record("intervals", "failed",
                               reason="tiers.json carries no pair count, so no interval "
                                      "can be computed; refusing rather than writing nulls")
        if d.get("union_removed") is None:
            return self.record("intervals", "failed",
                               reason="tiers.json carries no union_removed")
        # The operating point lives in the cascade's own manifest rather than in its result
        # file, so it is read from there. Reading it from the result file would have
        # produced a row of nulls and called it an operating point.
        mp = self.dir / "stages" / "cascade_manifest.json"
        if not mp.exists():
            return self.record("intervals", "failed",
                               reason="the cascade's manifest is absent, so the tier rows "
                                      "would ship without their operating point")
        par = json.loads(mp.read_text()).get("parameters", {})
        op = {"k2_confidence": par.get("k2_operating_confidence"),
              "min_match_length_nt": par.get("min_mlen"),
              "diamond_bitscore": par.get("diamond_bitscore"),
              "diamond_sensitivity": par.get("diamond_sensitivity"),
              "entropy_min": par.get("entropy_min"),
              "homopolymer_max": par.get("homopolymer_max"),
              "trim_r2": par.get("trim_r2")}
        tiers = d.get("tiers", {})
        rows = []
        for c in d.get("curve", []):
            row = dict(c)
            for key in ("tier_alone_removed", "marginal_removed", "cumulative_removed",
                        "residue"):
                if isinstance(c.get(key), (int, float)) and n:
                    row[f"{key}_interval"] = wilson(int(round(c[key] * n)), int(n))
            row["operating_point"] = op
            # The tier's own threshold sweep travels with the row, because the sweep IS the
            # finding: one database removes 42.8% or 3.5% of the same library on its
            # confidence setting alone, and a reader given one of those numbers without the
            # others has been told something unreadable.
            t = tiers.get(c.get("tier"), {})
            sweep = {k: v for k, v in t.items() if k.endswith("_sweep")}
            if sweep:
                row["threshold_sweep"] = sweep
                row["sweep_span_points"] = {
                    name: round((max(s.values()) - min(s.values())) * 100, 3)
                    for name, s in sweep.items() if s}
            missing = [k for k, v in op.items() if v is None]
            if missing:
                row["operating_point_incomplete"] = missing
            rows.append(row)
        # The three-way breakdown, which is the bar at the top of the web view. Tier 2 is
        # the OBSERVED population - contigs that assembled without being identified - and
        # it is deliberately computed from the observe tiers rather than from the curve,
        # because anything in the curve has been subtracted and tier 2 has not.
        comp = d.get("residue_composition") or {}
        tier3 = d.get("tier3_pairs")
        resid_pairs = d.get("residue_pairs")
        # Tier 2 is the UNION of the observe tiers, taken as residue minus tier 3 rather
        # than as the largest single observe tier. The observe sets overlap, so picking one
        # of them would leave the three segments not summing to the library and quietly
        # move the difference into tier 3 - the segment whose size is the whole point.
        tier2_pairs = (resid_pairs - tier3) if (resid_pairs is not None
                                                and tier3 is not None) else None
        breakdown = {
            "identified": wilson(int(round((d.get("union_removed") or 0) * n)), int(n))
            if n else None,
            "assembled_unidentified": wilson(int(tier2_pairs), int(n))
            if (tier2_pairs is not None and n) else None,
            "residual": wilson(int(tier3), int(n)) if (tier3 is not None and n) else None,
            "tier2_from": ("the union of " + ", ".join(sorted(comp)) if comp else
                           "no observe tier in this run, so tier 2 has no population and "
                           "the bar has two segments rather than three"),
            "tier2_per_observe_set": comp,
            "note": "tier 2 is observed and never subtracted. The three fractions share "
                    "one denominator and sum to the library; the per-set fractions below "
                    "do not, because the observe sets overlap",
        }
        # The reads scored here are a contiguous block of the FASTQ, so the Wilson
        # interval above is too narrow. Both are carried: Wilson because it is what the
        # counts support, and the widened one because it is what the sampling supports.
        # Neither replaces the other silently.
        for key in ("identified", "assembled_unidentified", "residual"):
            iv = breakdown.get(key)
            if not iv:
                continue
            de = DESIGN_EFFECT.get(key, DESIGN_EFFECT["_default"])
            half = 0.5 * (iv["hi"] - iv["lo"]) * de
            iv["design_effect"] = de
            iv["lo_block_corrected"] = max(0.0, iv["point"] - half)
            iv["hi_block_corrected"] = min(1.0, iv["point"] + half)
        breakdown["design_effect"] = DESIGN_EFFECT
        # Advisory, not fatal: the factor is the only estimate there is, and it is a
        # property of flow-cell layout rather than of the classifier. But on any library
        # other than the one it was measured on it is BORROWED, and the file says so.
        breakdown["design_effect_borrowed"] = "SRR23973193" not in str(self.args.r1)
        if breakdown["design_effect_borrowed"]:
            breakdown["design_effect_note"] = (
                "measured on SRR23973193 and applied here unmeasured; the block-corrected "
                "intervals are an estimate borrowed from another library")
        out = {"pairs": n, "operating_point": op, "breakdown": breakdown, "curve": rows,
               "interval_method": "wilson_95_on_this_run's_counts; a multi-seed run "
                                  "reports across-seed spread instead",
               "source": "tiers.json", "source_sha256_16": sha256(p)[:16]}
        (self.dir / "tier_table.json").write_text(json.dumps(out, indent=2))
        return self.record("intervals", "ok", rows=len(rows),
                           operating_point_complete=not any(v is None for v in op.values()))

    def layers(self) -> dict:
        """Order the residue most to least worth a human's attention.

        A LAYER IS A SORT ORDER WITH ITS EVIDENCE, NOT A SCORE. Each field used here has a
        measured direction on known material, which is the only reason it may be used at
        all: a field whose polarity is unmeasured would be ordering reads by database
        proximity and calling it interest. The layers, most interesting first:

          L1 recurs_and_unplaceable  no relative at significance, and the sequence recurs
                                     in the library. Rare AND unsequenced is the narrow
                                     honest claim; a singleton is more likely an artefact.
          L2 distant_relative        a relative exists but far away. This is where a new
                                     strain or a divergent homolog of a known pathogen
                                     lives, and it is the layer with an actual number
                                     attached: identity to closest relative, with coverage.
          L3 profile_only            no sequence relative but a family profile matches - a
                                     model reaches further than a sequence.
          L4 placed                  a close relative exists. Already known; kept so the
                                     reader can see the pile it was separated from.
          L5 artefact_shaped         low entropy, long homopolymer, adapter-like, or N-rich.
                                     Last on purpose and never silently dropped.

        A read with no field determinable lands in `undetermined`, which is a layer and not
        an error.
        """
        p = self.dir / "reads.jsonl"
        if not p.exists():
            return self.record("layers", "skipped" if "reads" not in self.args.stages
                               else "failed", reason="reads.jsonl is absent")
        order = ["recurs_and_unplaceable", "distant_relative", "profile_only", "placed",
                 "artefact_shaped", "undetermined"]
        counts: dict[str, dict[str, int]] = {}
        totals: dict[str, int] = {}
        coverage: dict[str, dict[str, int]] = {}
        out_rows = self.dir / "layers.jsonl"
        with open(p) as fh, open(out_rows, "w") as out:
            for line in fh:
                r = json.loads(line)
                pile = "residue" if r.get("set", "residue") == "residue" else "known"
                cov = coverage.setdefault(pile, {})
                for f in {f for fs in LAYER_FIELDS.values() for f in fs}:
                    if r.get(f) is not None or (f == "closest_*" and
                                                any(k.startswith("closest_") and
                                                    r[k] is not None for k in r)):
                        cov[f] = cov.get(f, 0) + 1
                layer, why = self._layer_of(r)
                counts.setdefault(pile, {})
                counts[pile][layer] = counts[pile].get(layer, 0) + 1
                totals[pile] = totals.get(pile, 0) + 1
                out.write(json.dumps({"name": r.get("name"), "pile": pile,
                                      "layer": layer, "evidence": why}) + "\n")
        # THE LAYERING IS LAID BESIDE ITS CONTROL, which is the whole reason it may be
        # trusted at all: the same rules are applied to reads the cascade DID identify. A
        # layer that fires as often on known reads as on the residue is ordering the reader
        # by something other than interest, and in this project the control moving with the
        # test is the tell that has caught every bad separator so far. Read it first.
        res_n, kn_n = totals.get("residue", 0), totals.get("known", 0)
        comparison = {}
        for k in order:
            a = wilson(counts.get("residue", {}).get(k, 0), res_n) if res_n else None
            b = wilson(counts.get("known", {}).get(k, 0), kn_n) if kn_n else None
            entry = {"residue": a, "known_positive": b,
                     "depends_on_fields": LAYER_FIELDS[k],
                     "field_coverage": {pile: {f: coverage.get(pile, {}).get(f, 0)
                                               for f in LAYER_FIELDS[k]}
                                        for pile in ("residue", "known")}}
            # A ZERO IN THE CONTROL IS ONLY A MEASUREMENT IF THE FIELD WAS MEASURED. The
            # per-read report currently runs its search fields on the residue alone, so for
            # any layer that depends on them the control pile reads zero by construction.
            # Calling that a separation would have been this project's exact failure mode -
            # a separator that reads the same on both piles describes the database - only
            # worse, because here the control was never taken.
            unmeasured = [f for f in LAYER_FIELDS[k]
                          if kn_n and coverage.get("known", {}).get(f, 0) == 0]
            if unmeasured and k != "artefact_shaped":
                entry["separates"] = None
                entry["reading"] = (
                    f"not comparable in this run: {', '.join(unmeasured)} was not computed "
                    f"on the known pile, so its zero is absence rather than a measurement")
            elif a and b:
                entry["separates"] = bool(a["lo"] > b["hi"] or b["lo"] > a["hi"])
                direction = ("enriched in the residue" if a["point"] > b["point"]
                             else "fires MORE on already-identified reads, which is a sign "
                                  "flip: whatever it reads, it is not interest")
                entry["reading"] = (f"separates at 95%, {direction}" if entry["separates"]
                                    else "cannot discriminate here: the intervals overlap")
            else:
                entry["reading"] = "known-positive not measured in this run"
            comparison[k] = entry
        summary = {
            "layer_order_most_to_least_interesting": order,
            "counts": counts, "totals": totals,
            "against_known_positive": comparison,
            "rule": dict(LAYER_RULES),
            "not_a_score": "layers are a sort order with evidence attached. Searching only "
                           "ever removes background; it never scores novelty, because "
                           "score ordering tracks database proximity",
        }
        (self.dir / "layers.json").write_text(json.dumps(summary, indent=2))
        return self.record("layers", "ok", totals=totals,
                           counts=counts.get("residue", {}),
                           separating_layers=[k for k, v in comparison.items()
                                              if v.get("separates")])

    @staticmethod
    def _layer_of(r: dict) -> tuple[str, dict]:
        """Assign one read to a layer, returning the evidence that put it there."""
        art = (r.get("entropy3") is not None
               and r.get("entropy3") < LAYER_RULES["artefact_entropy3_below"]) or \
              (r.get("homopolymer") or 0) > LAYER_RULES["artefact_homopolymer_above"] or \
              bool(r.get("adapter")) or bool(r.get("primer_dimer")) or \
              (r.get("n_frac") or 0) > 0.05
        if art:
            return "artefact_shaped", {
                "entropy3": r.get("entropy3"), "homopolymer": r.get("homopolymer"),
                "adapter": r.get("adapter"), "primer_dimer": r.get("primer_dimer"),
                "n_frac": r.get("n_frac")}

        best = None
        for key in [k for k in r if k.startswith("closest_")]:
            h = r.get(key)
            if isinstance(h, dict) and h.get("evalue") is not None:
                if best is None or h["evalue"] < best["evalue"]:
                    best = {**h, "field": key}
        sig = best is not None and best["evalue"] < SIGNIFICANCE_EVALUE
        close = bool(sig and (best.get("pident") or 0) >= CLOSE_IDENTITY_PERCENT)
        ab = r.get("library_abundance") or {}
        recurs = ab.get("median_31mer_count")
        profile = r.get("motif") or r.get("rfam")

        if sig and close:
            return "placed", {"closest_relative": best}
        if sig:
            return "distant_relative", {"closest_relative": best,
                                        "identity_percent": best.get("pident"),
                                        "alignment_length": best.get("aln_len")}
        if profile:
            return "profile_only", {"profile": profile, "no_sequence_relative_at": {
                "evalue": SIGNIFICANCE_EVALUE}}
        if recurs is not None:
            if recurs >= RECURRENCE_MIN_OCCURRENCES:
                return "recurs_and_unplaceable", {"occurrences": recurs,
                                                  "threshold": RECURRENCE_MIN_OCCURRENCES}
            return "undetermined", {"occurrences": recurs,
                                    "reason": "no relative at significance and the "
                                              "sequence does not recur, so an artefact "
                                              "and a singleton molecule cannot be "
                                              "separated here"}
        return "undetermined", {"reason": "neither a relative nor an abundance field was "
                                          "determinable for this read"}

    # -- the manifest ---------------------------------------------------------
    def finish(self) -> int:
        a = self.args
        bs = ROOT / "blind_spots.json"
        if bs.exists():
            shutil.copyfile(bs, self.dir / "blind_spots.json")
        # Field calibration is checked at write time rather than in review, because a field
        # printed without its calibration is the failure this project is built around. A
        # null says what the field reads on meaningless sequence; the known-positive says
        # whether it reads differently on sequence the cascade DID identify. The second
        # comes from the rows file, where subtracted reads travel beside residue reads, so
        # both files are consulted and neither absence is inferred from the other.
        fields = self.dir / "fields.json"
        rows = self.dir / "reads.jsonl"
        known_positive_rows = 0
        if rows.exists():
            # Parsed, not string-matched: the old test counted a row as known whenever the
            # text "residue" appeared nowhere in it, so a subject name or any field
            # containing that word moved a row between piles.
            with open(rows) as fh:
                for line in fh:
                    if json.loads(line).get("set", "residue") != "residue":
                        known_positive_rows += 1
        field_check: dict | str = {}
        if fields.exists():
            d = json.loads(fields.read_text())
            for name, v in (d.get("fields", d) or {}).items():
                if not isinstance(v, dict):
                    continue
                null = next((v[k] for k in ("null", "shuffled", "shuffled_null")
                             if k in v), None)
                real = v.get("residue")
                state = "null only" if null is not None else \
                    "cannot discriminate here: no null in this run"
                if null is not None and known_positive_rows:
                    state = "calibrated"
                note = None
                if isinstance(null, dict) and isinstance(real, dict) and \
                        null.get("mean") == real.get("mean") and \
                        null.get("sd") == real.get("sd"):
                    # Not a bug here but it must be stated: the shuffle preserves
                    # mononucleotide composition, so a composition field's null is the field
                    # itself. A reader shown the two side by side without this would read an
                    # agreement as a result.
                    note = ("null identical to the measurement by construction: the "
                            "shuffle preserves this field, so it has no null here and "
                            "cannot discriminate")
                    state = "cannot discriminate here: the null is the field itself"
                field_check[name] = {"state": state, "note": note} if note else state
        else:
            field_check = "fields.json absent in this run"
        bad = [s for s in self.stages if s["status"] == "failed"]
        skipped = [s for s in self.stages if s["status"] == "skipped"]
        manifest = {
            "tool": "hapax run", "spec": "SPEC.md 8.1",
            "script": "hapax.py", "script_sha256": sha256(ROOT / "hapax.py"),
            "git_sha": git("rev-parse", "HEAD"),
            "tree_dirty": bool(git("status", "--porcelain")),
            "seed": a.seed,
            "inputs": {"r1": str(a.r1), "r2": str(a.r2),
                       "r1_bytes": a.r1.stat().st_size if a.r1.exists() else None,
                       "r2_bytes": a.r2.stat().st_size if a.r2.exists() else None,
                       "trim_r2": a.trim_r2,
                       "trim_note": (
                           "measured, not detected: the non-template prefix is a "
                           "distribution with mode 14 (results/nontemplate_prefix/). The "
                           "trim is TWO operating points and this run uses one. 14 is "
                           "right for anything that reads sequence, because those bases "
                           "are real template; 16 is the better Kraken2 setting, because "
                           "confidence is a fraction of a read's minimizers, and it "
                           "identifies 1.03-1.11 points more on this library and "
                           "1.60-1.84 more on the 2022 one (5 seeds each). Passing 14 "
                           "everywhere therefore costs about a point of identification, "
                           "recorded here rather than left implicit.")},
            "parameters": {k: (str(v) if isinstance(v, Path) else
                               [str(x) for x in v] if isinstance(v, list) else v)
                           for k, v in vars(a).items()},
            "stages": self.stages,
            "complete": not bad and not skipped,
            "incomplete_because": ([f"{s['stage']} failed" for s in bad] +
                                   [f"{s['stage']} skipped: {s.get('reason')}"
                                    for s in skipped]) or None,
            "field_calibration": field_check,
            "known_positive_rows": known_positive_rows,
            "known_positive_note": ("subtracted reads travel in reads.jsonl as the "
                                    "known-positive comparison; a field that reads the "
                                    "same on both carries no information about the "
                                    "residue" if known_positive_rows else
                                    "no subtracted reads in this run, so every field's "
                                    "known-positive comparison is not measured"),
            "blind_spots": {"file": "blind_spots.json",
                            "sha256_16": sha256(bs)[:16] if bs.exists() else None},
            "interface_rules": [
                "no figure without its operating point",
                "no pile without its interval",
                "no score without the calibration run behind it"],
            "wall_seconds": round(time.perf_counter() - self.t0, 1),
            "python": platform.python_version(), "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
        }
        (self.dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"\n{self.dir}")
        for s in self.stages:
            print(f"  {s['stage']:<12} {s['status']:<8} "
                  f"{s.get('seconds', '')}{'s' if s.get('seconds') else ''}"
                  f"{'  ' + str(s.get('reason')) if s.get('reason') else ''}")
        print(f"  complete: {manifest['complete']}")
        return 1 if bad else 0


#: Written into layers.json so the rule is in the run directory rather than only in code.
#: How much wider than Wilson an interval has to be, because the reads scored here are a
#: CONTIGUOUS block of the FASTQ and not a random sample of the library. Measured over all
#: 8,731,224 pairs of SRR23973193 in `results/block_position_bias/`: the observed
#: between-block standard deviation divided by the binomial one a Wilson interval assumes,
#: with uniform random samples of the same size as the control reading 0.92-1.10.
#:
#: The structure sits in the abundant classes - adjacent blocks correlate at +0.25 and
#: +0.29 on ribosomal and identified - and is near absent from the residue classes, whose
#: autocorrelation is within noise of zero. So the residue, which is what this tool is
#: about, is very nearly a simple random sample and the correction to it is small.
#:
#: These are PAIR-LEVEL factors: a pair counts for a class when either mate carries it,
#: which is how the cascade removes a pair, and the residual factor is measured on pairs
#: whose BOTH mates are unplaced, which is tier 3 as this tool defines it. The mate-1 and
#: mate-2 figures differ and must not be substituted - on the tier-2 class they are 1.13
#: and 1.39, and a disagreement with another thread over exactly that class turned out to
#: be this rule and not a discrepancy in the quantity.
#:
#: Carried across from a label-based instrument (SILVA plus the sample's own contigs at 40
#: bases, no Kraken2 and no DIAMOND) to cascade numbers as the best available estimate,
#: not as the same measurement. It is a property of how the reads are arranged on the flow
#: cell rather than of the classifier, which is what makes the transfer defensible, and
#: every number that uses it carries the citation.
DESIGN_EFFECT = {
    "identified": 1.58,
    "assembled_unidentified": 1.27,
    "residual": 1.09,
    "_default": 1.27,
    "_source": "results/block_position_bias/REPORT.md, 5 seeds, measured on SRR23973193",
    "_why": "the cascade reads a contiguous block of the FASTQ, so a Wilson interval, "
            "which assumes a simple random sample of the library, is too narrow",
    "_control": "uniform random samples of the same size read 0.92-1.10 on the same "
                "statistic, so the instrument is sound and the blocks are the finding",
    "_never_from_seeds": "these come from all 87 windows of the library and must never be "
                         "re-estimated from a run's own five seeds: computed that way the "
                         "same statistic reads 0.65 and 0.44 and reverses its sign, "
                         "because five numbers cannot measure a variance. Reproduced "
                         "independently at 1.44 and 1.61 by experiments/block_dispersion.py",
}

SIGNIFICANCE_EVALUE = 1e-5
CLOSE_IDENTITY_PERCENT = 90.0
RECURRENCE_MIN_OCCURRENCES = 2
LAYER_RULES = {
    "significance_evalue": SIGNIFICANCE_EVALUE,
    "significance_note": "1e-5, because the chance band reads identically on residue and "
                         "already-identified reads: at E<1e-5 14.6% of residue reads have "
                         "a relative covering ~96% of the read at 90.1% identity, and a "
                         "weaker hit describes the database rather than the read "
                         "(results/partial_match_subjects/)",
    "close_identity_percent": CLOSE_IDENTITY_PERCENT,
    "recurrence_min_occurrences": RECURRENCE_MIN_OCCURRENCES,
    "recurrence_note": "a 1-in-a-million virus is a singleton against ~68 occurrences for "
                       "background, so recurrence is weak evidence per read and is used "
                       "only to separate 'rare and unsequenced' from 'could be artefact'",
    "artefact_entropy3_below": 0.5,
    "artefact_homopolymer_above": 25,
}

#: Which row fields each layer's assignment actually reads. Recorded so the control can be
#: checked for having been measured at all: a layer whose fields were computed on the
#: residue alone cannot be compared against the known pile, and its zero there is absence.
LAYER_FIELDS = {
    "recurs_and_unplaceable": ["closest_*", "library_abundance"],
    "distant_relative": ["closest_*"],
    "profile_only": ["closest_*", "motif"],
    "placed": ["closest_*"],
    "artefact_shaped": ["entropy3", "homopolymer", "adapter"],
    "undetermined": ["closest_*", "library_abundance"],
}

#: Stages that run an experiment script, and so are preflighted by calling it.
SCRIPT_STAGES = ["cascade", "reads", "linkage", "census", "calibration", "divergent",
                 "shortlist", "panel", "crossplant"]
#: Flags whose script default is one specific sample's file or setting, which hapax must
#: therefore always pass. The second element names the per-read field that may be left out
#: with --report-without instead.
SAMPLE_SPECIFIC = {
    "subtract_known.py": [("--r1", None), ("--r2", None), ("--trim-r2", None),
                          ("--selfasm", None)],
    "residue_report.py": [("--library-r1", "abundance"), ("--library-r2", "abundance"),
                          ("--trim-r2", None), ("--chimaera-ref", "chimaera"),
                          ("--cross-sample-ref", "cross_sample"), ("--spacers", "crispr")],
    "linkage_census.py": [("--r1", None), ("--r2", None), ("--trim-r2", None),
                          ("--rrna-contigs", None), ("--identified-contigs", None),
                          ("--unidentified-contigs", None)],
    "spike_holdout.py": [("--r1", None), ("--r2", None), ("--trim-r2", None),
                         ("--genome", None), ("--identified", None),
                         ("--unidentified", None), ("--cross-sample-ref", None),
                         ("--silva", None)],
    "mate_friends.py": [("--r1", None), ("--r2", None), ("--trim-r2", None)],
    "divergent_known.py": [("--r1", None), ("--r2", None), ("--trim-r2", None),
                           ("--contigs", None)],
    # The three promoted stages. Each resolves its reference data from the library, and
    # `--library` is the flag that decides WHICH sample's panel is ground truth: left out,
    # the stage would score this run against whatever panel happened to be committed.
    "removal_shortlist.py": [("--evidence", None)],
    "truth_panel.py": [("--library", None)],
    "crossplant_stage.py": [("--library", None)],
}

#: Scripts that accept `--library`, and are therefore told which one every time they run.
#:
#: `hapax run` requires `--library` of its caller and passed it on to two of them. The rest
#: take it through `library_facts.add_argument`, whose default is the discriminating
#: dataset, so a stage that was not told asked for Berlin's facts. `Run.call` now appends
#: it for every script in this set. The set is asserted against the scripts' own source by
#: `tests/test_hapax_stages.py`, so a script that becomes library-aware and is not added
#: here fails the suite rather than quietly inheriting Berlin.
LIBRARY_AWARE = {"subtract_known.py", "residue_report.py", "linkage_census.py",
                 "spike_holdout.py", "mate_friends.py", "divergent_known.py",
                 "truth_panel.py", "crossplant_stage.py"}

#: Scripts that implement --preflight: they check their own required inputs and tools.
PREFLIGHT_SCRIPTS = {"subtract_known.py", "residue_report.py", "linkage_census.py",
                     "spike_holdout.py", "removal_shortlist.py", "truth_panel.py",
                     "crossplant_stage.py"}

STAGE_ORDER = ["cascade", "contigs", "reads", "linkage", "census", "calibration",
               "intervals", "layers", "strains", "divergent", "shortlist", "panel",
               "crossplant"]

#: Stages that score the run against sequence whose distance was fixed BEFORE it, and so
#: are validation of the tool rather than output of it.
#:
#: Imo, 2026-11-23: *"why are we asking for depositor assembled genomes IN PROD THERE IS NO
#: DEPOSITOR"*, and earlier the same day: *"why would we, with unknown wastewater, have a
#: genome with a distance fixed before the run? by default we know nothing. we always know
#: nothing."* `panel` scores recovery against the depositors' own assemblies from that
#: water; `calibration` scores a spiked divergence ladder. Both are real and both are how
#: this tool is known to work at all - and neither has an answer on water nobody has
#: deposited from, which is every sample production is for. So they are not in the default
#: set: a production run does not ask the question, rather than asking it and recording an
#: absence. They stay selectable by name, which is how a fixture library runs them.
#:
#: This is the same conclusion the report reached from the other end (SPEC 13.4): no
#: ground-truth section on any page, because a section that renders empty on real water
#: tells its reader something is missing from their sample.
VALIDATION_STAGES = ["calibration", "panel"]
PRODUCTION_STAGES = [s for s in STAGE_ORDER if s not in VALIDATION_STAGES]


def main() -> int:
    ap = argparse.ArgumentParser(prog="hapax", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="subtract the known, then characterise what is left")
    r.add_argument("r1", type=Path)
    r.add_argument("r2", type=Path)
    r.add_argument("--seed", type=int, required=True)
    r.add_argument("--out", type=Path, required=True, help="the run directory")
    r.add_argument("--pairs", type=int, default=200_000,
                   help="pairs through the cascade. The block is seed-determined, so five "
                        "seeds are five disjoint blocks of the same library")
    r.add_argument("--block", type=int, default=None,
                   help="pairs skipped per seed, so seed N reads the block starting at "
                        "(N-1)*block. MUST BE SET FOR A LIBRARY SMALLER THAN BERLIN: the "
                        "stage default of 1,500,000 sends seed 4 past the end of a 4.3 M "
                        "read library, which fails the seed loudly and leaves the run "
                        "directory looking complete. Left unset, the stage's default is "
                        "used and recorded in its manifest")
    r.add_argument("--report-pairs", type=int, default=20_000,
                   help="residue pairs through the per-read report, whose per-read fields "
                        "are far costlier than the cascade's")
    r.add_argument("--report-without", nargs="+", default=[],
                   choices=["abundance", "cross_sample", "rfam", "crispr", "pfam",
                            "chimaera"],
                   help="per-read report fields to leave out ON PURPOSE. Without this, a "
                        "field whose reference or tool is absent stops the run before it "
                        "starts, rather than being dropped and read as a zero")
    r.add_argument("--linkage-pairs", type=int, default=20_000)
    r.add_argument("--census-queries", type=int, default=3_000,
                   help="pairs per pile through the linkage census. Both piles "
                        "always run: the control is what makes the residue "
                        "number readable")
    r.add_argument("--trim-r2", type=int, required=True,
                   help="non-template bases to cut from the start of mate 2. REQUIRED: it "
                        "is a property of the library kit (14 on SRR23973193, 0 on most "
                        "others), and the old default of 14 silently cut 14 real bases off "
                        "every other library's mate 2")
    r.add_argument("--library", required=True,
                   help="the run accession, as experiments/stage_references.json spells "
                        "it. REQUIRED for the same reason --trim-r2 is: the panels a run "
                        "is checked against are properties of the library, and a default "
                        "would score this sample against another one's ground truth. An "
                        "unknown library refuses to start; a known library with no panel "
                        "is skipped with the committed reason, which is a different fact")
    r.add_argument("--removal-evidence", type=Path, default=None,
                   help="removal_evidence.tsv from experiments/classify_contigs.py: the "
                        "identity, coverage, E-value and bitscore that justified removing "
                        "each contig, which the cascade computes and used to discard. The "
                        "`shortlist` stage ranks that pile by how much of each contig its "
                        "best hit EXPLAINS - the one mechanism measured to catch a novel "
                        "genome being subtracted as known")
    r.add_argument("--shuffled-removal-evidence", type=Path, default=None,
                   help="the same table over shuffled contigs: the shortlist's null. "
                        "Absent, the stage records the absence by name and carries the "
                        "nulls measured elsewhere as citations rather than implying one")
    r.add_argument("--minimap2", type=Path, default=None,
                   help="required when --assembly is given, for the `panel` stage: the "
                        "assembly is aligned to the truth panel with SECONDARY ALIGNMENTS "
                        "ON, because suppressing them credits one member of each "
                        "near-identical pair and read recovery as 9/25 instead of 20/25")
    r.add_argument("--crossplant-run-dir", type=Path, default=None,
                   help="a completed cross-plant pass downloaded from S3 (packed/, "
                        "sample/). Without it the `crossplant` stage cuts this run's "
                        "contigs into the panel's own query windows and stops there")
    r.add_argument("--crossplant-arm", default="reads", choices=["reads", "contigs"],
                   help="which arm of the pass in --crossplant-run-dir to report")
    r.add_argument("--single-end", action="store_true",
                   help="the library has no second mate, so r2 is r1 again and --trim-r2 "
                        "must be 0. Every cascade tier is a min, max or union over the "
                        "two mates, so a duplicated mate gives exactly the read-level "
                        "answer there; the per-read report is told separately to drop "
                        "mate 2 rather than duplicate it, and the two mate-linkage stages "
                        "refuse. Costs twice the search time for the same answer, which "
                        "is cheaper than a second code path through the cascade")
    r.add_argument("--diamond-bin", type=Path, default=None)
    r.add_argument("--classify-bin", type=Path, default=None,
                   help="Kraken2's raw `classify` binary, not the kraken2 wrapper")
    r.add_argument("--threads", type=int, default=4)
    r.add_argument("--k2-no-memory-map", action="store_true",
                   help="load each Kraken2 hash into RAM rather than memory-mapping it. "
                        "Identical calls and confidences; wall time and peak RSS only, "
                        "so those two fields are comparable only across runs that agreed "
                        "on this")
    r.add_argument("--contigs", type=Path, nargs="*", default=None,
                   help="REQUIRED for the cascade, no default: the default was Berlin's "
                        "assembly, which would have been subtracted from any library. "
                        "Pass it with no paths to run without. "
                        "The sample's own IDENTIFIED contigs, the most efficient "
                        "subtrahend this project has measured by orders of magnitude - "
                        "1.4 MB of its own ribosomal contigs beats all 272 MB of SILVA. "
                        "It must be the identified ones only: the subtrahend is "
                        "IDENTIFICATION, not assembly, and subtracting contigs that "
                        "nothing identified would remove the tier 2 population into tier "
                        "1 and inflate every removal figure. Measured here: the whole "
                        "contig set removes 85.9%% against 77.0%% for the identified "
                        "half")
    r.add_argument("--unidentified-contigs", type=Path, nargs="+", default=None,
                   help="contigs that assembled but nothing identified: tier 2's "
                        "population, OBSERVED and never subtracted. Subtracting these "
                        "would move tier 2 into tier 1 and inflate every removal figure")
    r.add_argument("--identified-contigs", type=Path, default=None,
                   help="contigs the cascade DID identify, excluding ribosomal ones. The "
                        "`divergent` stage measures how far each one was from the "
                        "reference that identified it, which the cascade computes and "
                        "then discards. A ribosomal match says what the sequence is and "
                        "not whose, so that pile is a different question")
    r.add_argument("--rrna-contigs", type=Path, default=None,
                   help="this sample's own ribosomal contigs; the census's rRNA class")
    r.add_argument("--assembly", type=Path, default=None,
                   help="this sample's own full assembly: the per-read chimaera test's "
                        "reference. No default: the old one was Berlin's assembly, which "
                        "on any other library is the wrong sample")
    r.add_argument("--cross-sample-ref", type=Path, default=None,
                   help="assemblies of OTHER samples from the same works, for cross-sample "
                        "recurrence (per-read report) and the ladder")
    r.add_argument("--spacers", type=Path, default=None,
                   help="CRISPR spacers extracted from this sample's assemblies")
    r.add_argument("--k2-dbs", type=Path, nargs="+",
                   default=[ROOT / "data/k2_viral", ROOT / "data/k2_standard08",
                            ROOT / "data/k2_pluspf08", ROOT / "data/k2_pluspf16"])
    r.add_argument("--silva", type=Path, nargs="*",
                   default=[ROOT / "data/silva/SSU.fa.gz", ROOT / "data/silva/LSU.fa.gz"])
    r.add_argument("--dmnd", type=Path, nargs="*",
                   default=[ROOT / "data/prot/viral.dmnd", ROOT / "data/prot/sprot.dmnd"])
    r.add_argument("--spike-genome", type=Path, default=None,
                   help="REQUIRED with --calibrate, no default: it was chosen as absent "
                        "from Berlin's assembly, which says nothing about another sample's. "
                        "Genomes for the divergence ladder. Only a genome ABSENT from "
                        "this sample's own assembly is a clean holdout - a resident one "
                        "is removed by the sample's own contigs at every rung and reads "
                        "as the tool working when it is the sample working "
                        "(results/spike_holdout/)")
    r.add_argument("--calibrate", action="store_true",
                   help="run the per-sample divergence ladder. Without it the run carries "
                        "no positive control and says so")
    r.add_argument("--stages", nargs="+", default=PRODUCTION_STAGES, choices=STAGE_ORDER,
                   help="default: every production stage. The validation stages "
                        + ", ".join(VALIDATION_STAGES) + " are selectable by name and are "
                        "not run by default, because they score the run against sequence "
                        "whose distance was fixed before it and unknown water has none")
    r.add_argument("--work", type=Path, default=None)
    r.add_argument("--plan", action="store_true", help="print the stages and exit")
    args = ap.parse_args()
    # Every stage runs with cwd at the repo root, so a relative path that hapax checked
    # from the caller's directory would name a different file (or none) inside the stage.
    # Made absolute once here, so what is checked, passed and recorded is one path.
    for k, v in vars(args).items():
        if isinstance(v, Path):
            setattr(args, k, v.absolute())
        elif isinstance(v, list) and v and all(isinstance(x, Path) for x in v):
            setattr(args, k, [x.absolute() for x in v])

    if args.plan:
        print("hapax run would execute, in order:")
        for s in args.stages:
            print(f"  {s:<12} {getattr(Run, s).__doc__.splitlines()[0]}")
        return 0

    run = Run(args)
    run.preflight()
    if run.problems:
        print("hapax run: REFUSING TO START - the run would not have done what was asked:",
              file=sys.stderr)
        for p in run.problems:
            print(f"  - {p}", file=sys.stderr)
        (run.dir / "manifest.json").write_text(json.dumps({
            "tool": "hapax run", "script_sha256": sha256(ROOT / "hapax.py"),
            "git_sha": git("rev-parse", "HEAD"), "seed": args.seed,
            "status": "refused to start", "complete": False,
            "problems": run.problems,
            "parameters": {k: (str(v) if isinstance(v, Path) else
                               [str(x) for x in v] if isinstance(v, list) else v)
                           for k, v in vars(args).items()}}, indent=2))
        return 2
    for name in STAGE_ORDER:
        if name not in args.stages:
            run.record(name, "skipped", reason="not selected with --stages")
            continue
        getattr(run, name)()
    return run.finish()


if __name__ == "__main__":
    raise SystemExit(main())
