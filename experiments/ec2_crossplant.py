"""Launch the cross-plant pass on EC2 (`crossplant_sample.py` then `crossplant_pass.py`).

A copy of `ec2_recurrence.py` with the panel, payload and template swapped; the spend gate,
dirty-tree refusal and terminate-in-finally are unchanged.
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
TEMPLATE = Path(__file__).resolve().parent / "crossplant_userdata.sh"
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

#: Pinned, for the reason every version in this repo is pinned: a number produced by
#: "latest" cannot be reproduced by anyone reading it later.
SRA_VERSION = "3.1.1"

#: The toolkit refuses to run without a settings file and generates a GUID when a human
#: answers its prompt. A fixed one is written instead; it identifies nothing and exists so
#: an unattended instance never waits for an answer that cannot come.
VDB_GUID = "00000000-0000-0000-0000-00000hapax01"

#: Eight vCPU is the account ceiling and this pass is scan-bound rather than memory-bound -
#: the k-mer table holds only the query set's distinct k-mers, a few hundred MB - so the
#: compute-optimised type is the right buy, unlike the big-index run's memory-heavy one.
DEFAULT_INSTANCE_TYPE = "c6i.2xlarge"

#: GB, gp3. Sized by the largest single subject: 3.6 Gbase as FASTQ is about 9 GB, plus its
#: SRA object and the toolkit's temp space. 150 GB is several times that, which is the right
#: margin for something that deletes each subject after scanning it; it is also cheap.
DEFAULT_VOLUME_GB = 150
DEFAULT_VOLUME_THROUGHPUT = 250
DEFAULT_VOLUME_IOPS = 3000

#: Minutes. Short, because the pilot is a rate measurement and a dead-man timer that is also
#: the cost projection should not be generous.
DEFAULT_HARD_STOP_MINUTES = 180
HEARTBEAT_SECONDS = 300

#: What the instance runs, listed by hand and checked by the UserData's own `test -f` lines,
#: because the reference project twice shipped a payload missing a module and the instance
#: died having read nothing.
PAYLOAD_MODULES = ("experiments/recurrence_pass.py",
                   "experiments/crossplant_sample.py",
                   "experiments/crossplant_pass.py",
                   "experiments/crossplant_reduce.py",
                   "experiments/kmer_coverage/Cargo.toml",
                   "experiments/kmer_coverage/Cargo.lock",
                   "experiments/kmer_coverage/src/main.rs")


def build_payload(root: Path) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for rel in PAYLOAD_MODULES:
            path = root / rel
            if not path.exists():
                raise SystemExit(f"{rel} is in PAYLOAD_MODULES and not on disk")
            tar.add(path, arcname=rel)
    return buf.getvalue()


def git_head(root: Path) -> tuple[str, bool, bool]:
    """HEAD, whether the shipped files are dirty, and whether the wider tree is.

    Two flags rather than one, for the reason `ec2_bigindex.py` gives: what has to be clean
    is what the instance will run, and that is checkable exactly, while this session has
    result files landing in the tree the whole time.
    """
    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True,
                              text=True, check=True).stdout.strip()

    shipped = [*PAYLOAD_MODULES, "experiments/ec2_crossplant.py",
               "experiments/crossplant_userdata.sh", "experiments/ec2_common.py",
               "experiments/crossplant_panel.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *shipped)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE,
                    choices=ec2lib.ALLOWED_TYPES)
    ap.add_argument("--queries", type=Path,
                    default=ROOT / "results/recurrence/queries_s1.tsv",
                    help="unused; the Berlin query set is read from scratch/recurrence-s1-5/")
    ap.add_argument("--panel", type=Path, default=ROOT / "results/crossplant/panel.tsv")
    ap.add_argument("--query-set", default="cp1",
                    help="names the uploaded set. The bucket has no delete, so a set is "
                         "written once under a fixed key and a changed set gets a new name")
    ap.add_argument("--query-run", default="SRR23973193",
                    help="the sample the queries come from; excluded from its own panel, "
                         "which is the one exclusion that would silently invent recurrence")
    ap.add_argument("--subjects", type=int, default=8,
                    help="how many panel libraries to scan, spread evenly over the panel's "
                         "dates. 0 means all of them")
    ap.add_argument("--budget-hours", type=float, default=2.0,
                    help="the pass stops starting subjects after this long, inside the "
                         "dead-man timer, so an underestimated rate costs a short table "
                         "rather than a dead run")
    ap.add_argument("--depth", type=int, default=2_000_000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--shards", type=int, default=1,
                    help="how many subjects to process at once. The scan is one "
                         "single-threaded sweep per library, so one subject at a time "
                         "leaves most of eight vCPU idle; each shard also holds its own "
                         "copy of the k-mer table and its own FASTQ on disk, which is what "
                         "bounds this rather than the cores")
    ap.add_argument("--volume-size", type=int, default=DEFAULT_VOLUME_GB)
    ap.add_argument("--volume-throughput", type=int, default=DEFAULT_VOLUME_THROUGHPUT)
    ap.add_argument("--volume-iops", type=int, default=DEFAULT_VOLUME_IOPS)
    ap.add_argument("--hard-stop-minutes", type=int, default=DEFAULT_HARD_STOP_MINUTES)
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--no-spot", action="store_true",
                    help="on demand. Spot was reclaimed twice in four hours mid-search, and "
                         "this pass is a long sequence of short independent scans, so a "
                         "reclaim costs one library - but on demand is also only ~$0.34/h")
    ap.add_argument("--allow-price-change", action="store_true",
                    help="pass only with Imo's explicit yes on the figure printed below")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--timeout-minutes", type=float, default=200.0)
    ap.add_argument("--poll-seconds", type=float, default=300.0)
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

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
            f"pre-approved for a single run. Put that figure on a decision card and pass "
            f"--allow-price-change once he has answered, or shorten --hard-stop-minutes: "
            f"the projection is the timer's worst case, not the expected bill.")

    commit, shipped_dirty, tree_dirty = git_head(ROOT)
    if shipped_dirty and not args.allow_dirty:
        raise SystemExit(
            "refusing to launch: the code that would run on the instance is not committed. "
            "The instance has no checkout, so the commit is the only thing that makes its "
            "output traceable. Commit the launcher, the template and the payload first.")

    for path in (args.panel,):
        if not path.exists():
            raise SystemExit(f"{path} is missing; the instance would start on nothing")

    bucket, prefix = ec2lib.BUCKET, ec2lib.PREFIX
    code_key = f"{prefix}/code/crossplant-{commit}.tar.gz"
    tokens = dict(
        COMMIT=commit, BUCKET=bucket, PREFIX=prefix, REGION=ec2lib.REGION,
        HARD_STOP_MINUTES=args.hard_stop_minutes, HEARTBEAT_SECONDS=HEARTBEAT_SECONDS,
        CODE_KEY_NAME=Path(code_key).name, QUERY_SET=args.query_set,
        QUERY_RUN=args.query_run, SUBJECT_COUNT=args.subjects,
        BUDGET_HOURS=args.budget_hours, DEPTH=args.depth, THREADS=args.threads,
        SHARDS=args.shards,
        # fasterq-dump is the only threaded step, so the threads are divided between the
        # shards rather than each shard asking for all of them. The reference project lost a
        # fourteen-hour ceiling to eight workers each asking for three threads on eight
        # cores, with nothing in the launcher able to notice (ec2_common.py).
        SHARD_THREADS=max(1, args.threads // max(1, args.shards)),
        SRA_VERSION=SRA_VERSION, VDB_GUID=VDB_GUID,
    )
    run_script = render(TEMPLATE, None, **tokens)
    run_key = (f"{prefix}/code/crossplant-run-{commit[:12]}-"
               f"{hashlib.sha256(run_script.encode()).hexdigest()[:12]}.sh")
    script = render(BOOTSTRAP, 16_384, COMMIT=commit, BUCKET=bucket, REGION=ec2lib.REGION,
                    HARD_STOP_MINUTES=args.hard_stop_minutes, RUN_SCRIPT_KEY=run_key)

    if args.print_user_data:
        print(script)
        print(f"\n### ---- run script, fetched from s3://{bucket}/{run_key} ----\n")
        print(run_script)
        print(f"\n### {len(script.encode())} bytes of UserData, "
              f"{len(run_script.encode()):,} of run script", file=sys.stderr)
        return 0

    import boto3

    s3 = ec2lib.s3_client()
    ec2 = boto3.client("ec2", region_name=ec2lib.REGION)
    print(f"commit {commit[:12]}, {args.subjects or 'all'} subjects, "
          f"panel {args.panel}", flush=True)

    if args.attach:
        return watch(ec2, args.attach, args)

    if not args.dry_run:
        s3.put_object(Bucket=bucket, Key=code_key, Body=build_payload(ROOT))
        print(f"payload at s3://{bucket}/{code_key}", flush=True)
        s3.put_object(Bucket=bucket, Key=run_key, Body=run_script.encode())
        print(f"run script at s3://{bucket}/{run_key}", flush=True)
        for path, name in ((args.panel, "panel.tsv"),):
            key = f"{prefix}/scratch/crossplant-{args.query_set}/{name}"
            s3.put_object(Bucket=bucket, Key=key, Body=path.read_bytes())
            print(f"  {key}  {path.stat().st_size / 1e6:.1f} MB", flush=True)

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
                                     {"Key": "Task", "Value": "hapax-crossplant"}]}],
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
        "query_set": args.query_set, "query_run": args.query_run,
        "subjects": args.subjects, "depth": args.depth, "budget_hours": args.budget_hours,
        "volume_gb": args.volume_size, "spot": not args.no_spot,
        "hard_stop_minutes": args.hard_stop_minutes, "projected_usd": round(cost, 2),
        "sra_toolkit": SRA_VERSION,
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/crossplant/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/crossplant-{iid}/",
    }, indent=2)
    (ROOT / f"results/crossplant/launch_{iid}.json").write_text(record + "\n")

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
