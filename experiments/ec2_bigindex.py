"""Launch the one EC2 run Imo approved: the big-index pass. Plan in
`results/subtract_known/BIGINDEX_RUN_PLAN.md`.

WHAT THIS BUYS, and it is three things that cannot be bought in this container. Requirement
1.iii - "the sequences with no match are really out at <80% identity" - is currently
answered against indexes that fit on a 15 GB disk plus 60 reads per class sent to the public
NCBI service, which is rate-limited and runs against a database that moves under us with no
version to cite. The far end of the removal curve is extrapolated from two database sizes
whose linear and superlinear fits disagree by a factor of fifteen. And the residue name
lists a sibling thread needs, to find out whether its wall is real or clutter from known
material, need a residue defined against something bigger than Swiss-Prot.

WHY THE RUN'S FIRST ACT IS A MEASUREMENT OF ITSELF. There is no honest figure for a blastn
pass against a 318 GB core_nt from inside here, so the script does not pretend to have one:
it searches a small pilot set, projects the full pass from the rate it measures, and refuses
the full pass if the projection does not fit the budget. A pass that cannot finish produces
no table and spends the money anyway, which is the failure this ordering exists to prevent.

COPIED FROM THE GUIDE'S WORKING EXAMPLE rather than written fresh, as the guide instructs:
`ramics/evo-training-check`'s `tools/ec2_detect_census.py`, with its EC2 facts in
`experiments/ec2_common.py`. The shape that matters is the one every lesson in that file is
about - dry run first, dead man's switch before anything that can fail, terminate in a
`finally`, and the instance pulling pinned code from S3 rather than being handed a working
tree.

THE ONE PLACE THIS DIVERGES FROM THE GUIDE is the volume: the guide says `<= 200 GB`, and
Imo raised the ceiling to 2 TB on 2026-11-20 for this run specifically, which is what makes
staging core_nt and nr possible at all. His message is the newer instruction; the guide's
number is recorded here so the difference is visible rather than silent.
"""

from __future__ import annotations

import argparse
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

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = Path(__file__).resolve().parent / "bigindex_userdata.sh"

#: The pinned NCBI copy. `latest-dir` in that bucket read `2026-10-22-01-05-01` on
#: 2026-11-20; pinning it rather than following `latest-dir` is what makes a number from
#: this run citable six months from now, and it is the entire reason for staging our own
#: copy instead of continuing to query the public service.
NCBI_BUCKET = "ncbi-blast-databases"
NCBI_PREFIX = "2026-10-22-01-05-01"

#: BLAST+ reads the staged volumes and nothing else does. DIAMOND is pinned to the version
#: every timing in `results/diamond_scaling/` was measured with, because a third scaling
#: point measured on a different binary is not a third point on the same curve.
BLAST_VERSION = "2.17.0"
DIAMOND_VERSION = "2.1.10"

#: Memory is why this type: BLAST memory-maps its volumes, so 64 GB is comfort rather than
#: requirement for the search, but the UniRef100 DIAMOND build that may follow wants it.
#: Eight vCPU is the account ceiling either way, so the cheaper types buy nothing but risk.
DEFAULT_INSTANCE_TYPE = "r6i.2xlarge"

#: GB, gp3. core_nt 319 + nr 928 + a UniRef100 build with headroom, under Imo's 2 TB.
DEFAULT_VOLUME_GB = 1700

#: Staging 1.25 TB at gp3's 125 MB/s baseline is about three hours of a 48 hour budget, and
#: throughput is the cheapest thing on this machine to buy: 600 MB/s costs about $0.04/h
#: and turns that into forty minutes.
DEFAULT_VOLUME_THROUGHPUT = 600
DEFAULT_VOLUME_IOPS = 6000

#: Minutes. The approved ceiling is 48 hours and this is the switch that enforces it without
#: depending on anything in the run working, including this launcher staying alive.
DEFAULT_HARD_STOP_MINUTES = 2880

HEARTBEAT_SECONDS = 300

#: Everything the instance runs. Derived by hand and checked by the PRECHECK block in the
#: UserData, because the reference project twice added a module to a tool's imports and not
#: to a list like this, and the instance died on an ImportError having read nothing.
PAYLOAD_MODULES = ("experiments/bigindex_search.py",)

STAGES = ("core_nt", "pilot", "core_nt_full", "nr")


def render_user_data(**tokens: Any) -> str:
    """The UserData, with @@TOKENS@@ filled in. Every token must be consumed.

    An unreplaced token is a bash line that silently does the wrong thing - `--region
    @@REGION@@` would be a literal - and on an instance nobody can log into, that is a
    wasted launch. So this refuses rather than renders a template with a token left in it.
    """
    text = TEMPLATE.read_text()
    for key, value in tokens.items():
        text = text.replace(f"@@{key}@@", str(value))
    leftover = sorted({t.split("@@")[0] for t in text.split("@@")[1::2]})
    if leftover:
        raise SystemExit(f"UserData still has unreplaced tokens: {leftover}")
    if len(text.encode()) > 16_384:
        raise SystemExit(f"UserData is {len(text.encode())} bytes; EC2's limit is 16,384")
    return text


def build_payload(root: Path) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for rel in PAYLOAD_MODULES:
            tar.add(root / rel, arcname=rel)
    return buf.getvalue()


def git_head(root: Path) -> tuple[str, bool, bool]:
    """HEAD, whether the SHIPPED files are dirty, and whether the wider tree is.

    Two flags rather than one, and the distinction is the point. The run plan makes a clean
    tree a property of the run: because the code is pinned, every number the instance
    produces is reproducible from a resolved SHA instead of a `-dirty` one. But this session
    has five-seed jobs writing into `results/` the whole time, so the wider tree is almost
    never clean, and refusing on that would mean either waiting hours or committing
    half-finished result files to buy a cosmetic check.

    What has to be clean is what the instance will run: the payload, this launcher and the
    UserData template. That is checkable exactly - `git status` over those paths - and it is
    the claim actually being made. The wider tree's state is recorded beside it rather than
    enforced, so a reader can see it was dirty and in which sense.
    """
    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True,
                              text=True, check=True).stdout.strip()

    shipped = [*PAYLOAD_MODULES, "experiments/ec2_bigindex.py",
               "experiments/bigindex_userdata.sh", "experiments/ec2_common.py",
               "experiments/bigindex_queries.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *shipped)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


def upload_queries(s3: Any, bucket: str, prefix: str, query_dir: Path,
                   query_set: str) -> int:
    """Put the query FASTAs where the instance will look for them.

    Under `scratch/`, not `results/`: the bucket has no delete, so a query set uploaded
    under a results key would be a permanent object that is not a result. Keys are fixed
    before launch for the same reason.
    """
    sent = 0
    for path in sorted(query_dir.glob("*.fa")):
        key = f"{prefix}/scratch/bigindex-queries-{query_set}/{path.name}"
        s3.put_object(Bucket=bucket, Key=key, Body=path.read_bytes())
        print(f"  {key}  {path.stat().st_size / 1e6:.1f} MB", flush=True)
        sent += 1
    if not sent:
        raise SystemExit(f"no *.fa in {query_dir}; the instance would start on nothing")
    return sent


def count_queries(query_dir: Path, pattern: str) -> int:
    n = 0
    for path in sorted(query_dir.glob(pattern)):
        with open(path) as h:
            n += sum(1 for line in h if line.startswith(">"))
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE,
                    choices=ec2lib.ALLOWED_TYPES)
    ap.add_argument("--stages", nargs="+", default=list(STAGES), choices=STAGES,
                    help="which blocks of the UserData run. The default is the whole "
                         "approved run; the pilot stage is what gates the rest")
    ap.add_argument("--query-dir", type=Path,
                    default=Path("/tmp/claude-0/wk/bigindex_queries"))
    ap.add_argument("--query-set", default="v1",
                    help="names the scratch prefix the queries go to. Change it rather "
                         "than overwriting: the bucket has no delete, so a second upload "
                         "under the same name leaves two objects and the instance would "
                         "read whichever won")
    ap.add_argument("--skip-query-upload", action="store_true",
                    help="the set is already up under this name")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--volume-size", type=int, default=DEFAULT_VOLUME_GB)
    ap.add_argument("--volume-throughput", type=int, default=DEFAULT_VOLUME_THROUGHPUT)
    ap.add_argument("--volume-iops", type=int, default=DEFAULT_VOLUME_IOPS)
    ap.add_argument("--hard-stop-minutes", type=int, default=DEFAULT_HARD_STOP_MINUTES)
    ap.add_argument("--full-budget-hours", type=float, default=20.0,
                    help="the wall clock the full core_nt pass is allowed to project to "
                         "before the run refuses it and stops at the pilot")
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--no-spot", action="store_true",
                    help="on demand is 3-4x the spot price Imo approved, so this needs "
                         "--allow-price-change and a new approval")
    ap.add_argument("--allow-price-change", action="store_true")
    ap.add_argument("--timeout-minutes", type=float, default=100.0)
    ap.add_argument("--poll-seconds", type=float, default=300.0)
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true",
                    help="launch, print the id, stop watching. The instance is then held "
                         "only by its own timer; this session's background jobs cannot "
                         "outlive a run measured in hours, and a watcher killed mid-wait "
                         "would fire its finally and terminate a healthy instance")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

    if args.no_spot and not args.allow_price_change:
        raise SystemExit(
            "refusing to launch on demand: the $25 Imo approved is a spot figure "
            "(r6i.2xlarge at $0.29/h). On demand is 3-4x that, which is a new approval, "
            "not a flag. Pass --allow-price-change once he has given one.")

    commit, shipped_dirty, tree_dirty = git_head(ROOT)
    if shipped_dirty and not args.allow_dirty:
        raise SystemExit(
            "refusing to launch: the code that would run on the instance is not committed. "
            "The instance has no checkout, so the commit is the only thing that makes its "
            "tables traceable - and the run plan makes that a property of the run rather "
            "than a detail. Commit the launcher, the template and the payload first.")

    bucket, prefix = ec2lib.BUCKET, ec2lib.PREFIX
    code_key = f"{prefix}/code/bigindex-{commit}.tar.gz"
    full_queries = count_queries(args.query_dir, "block_*.fa") if args.query_dir.exists() else 0

    tokens = dict(
        COMMIT=commit, BUCKET=bucket, PREFIX=prefix, REGION=ec2lib.REGION,
        HARD_STOP_MINUTES=args.hard_stop_minutes, HEARTBEAT_SECONDS=HEARTBEAT_SECONDS,
        STAGES=" ".join(args.stages), BLAST_VERSION=BLAST_VERSION,
        DIAMOND_VERSION=DIAMOND_VERSION, CODE_KEY_NAME=Path(code_key).name,
        QUERY_SET=args.query_set, THREADS=args.threads,
        NCBI_BUCKET=NCBI_BUCKET, NCBI_PREFIX=NCBI_PREFIX,
        FULL_QUERIES=full_queries, FULL_BUDGET_HOURS=args.full_budget_hours,
    )
    script = render_user_data(**tokens)

    if args.print_user_data:
        print(script)
        print(f"\n### {len(script.encode())} bytes of UserData; "
              f"{full_queries:,} queries in the full pass", file=sys.stderr)
        return 0

    import boto3

    s3 = ec2lib.s3_client()
    ec2 = boto3.client("ec2", region_name=ec2lib.REGION)

    print(f"commit {commit[:12]}, stages: {' '.join(args.stages)}, "
          f"{full_queries:,} queries in the full pass", flush=True)

    if args.attach:
        return watch(ec2, args.attach, args)

    # Not on a dry run: a dry run authorizes and executes nothing, and this bucket has no
    # delete, so an upload here would leave permanent objects behind for a launch that
    # never happened.
    if not args.dry_run:
        s3.put_object(Bucket=bucket, Key=code_key, Body=build_payload(ROOT))
        print(f"payload at s3://{bucket}/{code_key}", flush=True)
        if not args.skip_query_upload:
            n = upload_queries(s3, bucket, prefix, args.query_dir, args.query_set)
            print(f"{n} query file(s) uploaded as set {args.query_set}", flush=True)

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
                                     {"Key": "Task", "Value": "hapax-bigindex"}]}],
        UserData=script,
    )
    if market is not None:
        kwargs["InstanceMarketOptions"] = market

    if args.dry_run:
        try:
            ec2.run_instances(DryRun=True, **kwargs)
        except Exception as exc:  # DryRunOperation is the success case
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
    (ROOT / "results/subtract_known/bigindex_launch.json").write_text(json.dumps({
        "instance_id": iid, "commit": commit, "shipped_dirty": shipped_dirty,
        "wider_tree_dirty": tree_dirty,
        "instance_type": args.instance_type, "stages": args.stages,
        "query_set": args.query_set, "full_queries": full_queries,
        "volume_gb": args.volume_size, "spot": not args.no_spot,
        "hard_stop_minutes": args.hard_stop_minutes,
        "ncbi_prefix": NCBI_PREFIX, "blast": BLAST_VERSION, "diamond": DIAMOND_VERSION,
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/bigindex/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/bigindex-{iid}/",
    }, indent=2))

    if args.detach:
        print(f"detached: {iid} is held by its own {args.hard_stop_minutes}-minute timer "
              f"and by terminate-on-shutdown.\n"
              f"  watch with --attach {iid}; terminate explicitly to stop it early.",
              flush=True)
        return 0

    try:
        return watch(ec2, iid, args)
    finally:
        # Both, not either: the boot timer covers a launcher that dies, this covers a
        # UserData that fails before reaching its own shutdown.
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
