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

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = Path(__file__).resolve().parent / "prod_userdata.sh"
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

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
                   "experiments/prod_run.sh")


def payload_files(root: Path) -> list[str]:
    """Everything the instance runs: the listed files plus all of `experiments/`.

    The reference project twice shipped a payload missing a module and the instance died
    having read nothing, which is why `ec2_crossplant.py` lists its four files by hand.
    That list cannot be written for this run: `hapax run` dispatches to eleven scripts
    which import each other and read three JSON registries, and a missing one is a stage
    that refuses mid-run after the expensive half. So the unit is the directory, and the
    named files above are asserted on the instance one by one.
    """
    files = list(PAYLOAD_MODULES)
    for p in sorted((root / "experiments").rglob("*")):
        if p.is_dir() or "codon_assembler" in p.parts or "__pycache__" in p.parts:
            continue
        if p.suffix in (".py", ".json", ".sh", ".md", ".toml", ".lock", ".rs"):
            files.append(str(p.relative_to(root)))
    return sorted(set(files))


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
               "experiments/ec2_common.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *watched)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


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
    ap.add_argument("--stages", default="cascade contigs reads linkage census intervals "
                                        "layers strains shortlist panel crossplant",
                    help="hapax stages. `calibration` is left out because the spiked "
                         "ladder needs a cross-sample reference assembly this run does "
                         "not stage, and `divergent` because its gene-level workup is "
                         "hours of BLAST over the removed pile")
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
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

    # The two committed registries decide whether this library can be run at all, and they
    # are read HERE rather than on the instance so that a missing entry costs nothing.
    import library_facts
    import stage_refs
    facts = library_facts.library(args.library)          # raises UnknownLibrary
    block = library_facts.fact(args.library, "block")    # raises if absent
    for kind in ("truth_panel", "crossplant_panel"):
        r = stage_refs.resolve(args.library, kind)
        print(f"{kind}: {r['state']}"
              + (f" ({r['panel']})" if r["state"] == "resolved" else ""))
    if facts.get("layout") != "paired":
        raise SystemExit(f"{args.library} is {facts.get('layout')} in libraries.json; this "
                         f"template stages two mates and measures a mate-2 prefix")
    last = max(int(s) for s in args.seeds.split())
    need = (last - 1) * block + args.pairs
    if facts.get("pairs") and need > facts["pairs"]:
        raise SystemExit(
            f"refusing to launch: seed {last} would start at {(last - 1) * block:,} in a "
            f"library of {facts['pairs']:,} pairs. The instance checks this again against "
            f"the counted reads, but there is no reason to pay for the arithmetic.")
    print(f"{args.library}: {facts.get('pairs', 'unknown'):,} pairs, block {block:,}, "
          f"{last} seeds x {args.pairs:,} needs {need:,}")
    if facts.get("trim_r2") is None:
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

    files = payload_files(ROOT)
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
        PAIRS=args.pairs, SEEDS=args.seeds, BLOCK=block,
        REPORT_PAIRS=args.report_pairs, STAGES=args.stages,
        PREFIX_SEEDS=args.prefix_seeds, PREFIX_PAIRS=args.prefix_pairs,
        THREADS=args.threads, MEGAHIT_MEMORY=args.megahit_memory,
        SRA_VERSION=SRA_VERSION, VDB_GUID=VDB_GUID,
        DIAMOND_VERSION=DIAMOND_VERSION, KRAKEN2_VERSION=KRAKEN2_VERSION,
        KRAKEN2_URL=KRAKEN2_URL, MEGAHIT_URL=MEGAHIT_URL, MINIMAP2_URL=MINIMAP2_URL,
        SILVA_SSU_URL=SILVA_SSU_URL, SILVA_LSU_URL=SILVA_LSU_URL,
        K2_DBS=" ".join(K2_DBS), K2_DB_CONTIGS=K2_DB_CONTIGS,
    )
    run_script = render(TEMPLATE, None, **tokens)
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
        "trim_r2": "measured on the instance; libraries.json has none for this library",
        "volume_gb": args.volume_size, "spot": not args.no_spot,
        "hard_stop_minutes": args.hard_stop_minutes, "projected_usd": round(cost, 2),
        "sra_toolkit": SRA_VERSION, "diamond": DIAMOND_VERSION,
        "kraken2": KRAKEN2_VERSION, "megahit": MEGAHIT_VERSION,
        "minimap2": MINIMAP2_VERSION, "kraken2_dbs": list(K2_DBS),
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/prod/{args.library}/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/prod-{args.library}-{iid}/",
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
