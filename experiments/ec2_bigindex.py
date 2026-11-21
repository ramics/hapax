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

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = Path(__file__).resolve().parent / "bigindex_userdata.sh"

#: What actually goes in UserData. EC2 caps that at 16,384 bytes and the run script passed
#: 15,100 with four stages, so the run script is uploaded to our own prefix and this fetches
#: it. See the comment at the top of the bootstrap for what the split costs.
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

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

#: US-east-1 on-demand rates, for the pre-launch cost projection only. Rough and rounded up;
#: the projection exists to catch an order-of-magnitude mistake before it bills, not to
#: predict the invoice. Spot is taken as 0.35x on-demand, which is conservative against the
#: ~0.3x actually seen.
ON_DEMAND_USD_PER_HOUR = {"c6i.2xlarge": 0.34, "c6i.4xlarge": 0.68, "c6i.8xlarge": 1.36,
                          "r6i.2xlarge": 0.50, "r6i.4xlarge": 1.01, "m6i.2xlarge": 0.38,
                          "m6i.4xlarge": 0.77}
SPOT_FRACTION = 0.35
#: Imo's standing project rule, 2026-11-21: "a single EC2 spend under $10 is autoapproved".
#: Above it, a launch still needs his explicit yes on a figure.
AUTO_APPROVED_USD = 10.0


def projected_cost(instance_type: str, hard_stop_minutes: int, spot: bool,
                   volume_gb: int) -> tuple[float, str]:
    """Worst-case dollars for one run: the dead-man timer is the ceiling on its life.

    The timer is the only thing that bounds the bill without depending on anything in the
    run working, so it is what the projection uses - not the budget hours the stages are
    gated on, which a crash loop could sit inside without ever finishing.
    """
    hours = hard_stop_minutes / 60.0
    rate = ON_DEMAND_USD_PER_HOUR.get(instance_type)
    if rate is None:
        return float("nan"), (f"{instance_type} is not in ON_DEMAND_USD_PER_HOUR, so no "
                              f"projection is possible; add it rather than guessing")
    compute = rate * (SPOT_FRACTION if spot else 1.0) * hours
    # gp3 at $0.08/GB-month, charged for the hours the volume exists.
    storage = volume_gb * 0.08 * hours / 730.0
    total = compute + storage
    return total, (f"{instance_type} {'spot' if spot else 'on demand'} for up to "
                   f"{hours:.1f}h (the dead-man timer) is about ${compute:.2f} of compute "
                   f"plus ${storage:.2f} for {volume_gb} GB = **${total:.2f}**")

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

#: Seconds. CHUNK TARGET, AND WHY IT FOLLOWS THE MARKET RATHER THAN THE WORK.
#:
#: Chunking exists for one reason: spot reclaimed two instances in four hours, rows upload
#: only when a search finishes, so a reclaim costs whatever chunk it lands in. 1500 s was
#: chosen to make that loss about twenty-five minutes.
#:
#: But a chunk pays the database scan again, and the scan is the dominant term. Measured on
#: this box, core_nt at 318.8 GB over 8 threads: scan 1082 s fixed, marginal 1.973 s/read.
#: At a 1500 s target that sizes a chunk to 211 reads, so 72% of every chunk is re-scanning
#: a database whose contents have not changed, and the effective rate is 7.11 s/read against
#: a marginal 1.97. At a three-hour target the chunk is ~4,900 reads, re-scan is 10%, and the
#: effective rate is 2.19 s/read - 3.2x the reads for the same money.
#:
#: On demand is not reclaimed, so that 3.2x buys nothing but exposure to our own hard stop,
#: which the budget gate already handles. Paying the spot insurance premium on an on-demand
#: box is simply a 3.2x bill, and that is what the first on-demand launch did: it was given
#: 1500 s because the flag defaulted on the work rather than on the market. Hence the default
#: is resolved from --no-spot below, and the overriding flag stays for a measurement that
#: deliberately wants small chunks.
ON_DEMAND_CHUNK_TARGET_SECONDS = 10_800.0
SPOT_CHUNK_TARGET_SECONDS = 1_500.0

#: Everything the instance runs. Derived by hand and checked by the PRECHECK block in the
#: UserData, because the reference project twice added a module to a tool's imports and not
#: to a list like this, and the instance died on an ImportError having read nothing.
PAYLOAD_MODULES = ("experiments/bigindex_search.py",)

STAGES = ("core_nt", "pilot", "core_nt_full", "core_nt_fast", "nt_viruses",
          "alignment", "nr")


def render(template: Path, limit: int | None, **tokens: Any) -> str:
    """One template with @@TOKENS@@ filled in. Every token must be consumed.

    An unreplaced token is a bash line that silently does the wrong thing - `--region
    @@REGION@@` would be a literal - and on an instance nobody can log into, that is a
    wasted launch. So this refuses rather than renders a template with a token left in it.

    `limit` is EC2's UserData cap and applies only to the bootstrap now. The run script is
    fetched from S3, so it has no ceiling - which is the whole reason for the split, since
    the alternative was deleting the reasoning from the artefact that gets reviewed.
    """
    text = template.read_text()
    for key, value in tokens.items():
        text = text.replace(f"@@{key}@@", str(value))
    leftover = sorted({t.split("@@")[0] for t in text.split("@@")[1::2]})
    if leftover:
        raise SystemExit(f"{template.name} still has unreplaced tokens: {leftover}")
    if limit is not None and len(text.encode()) > limit:
        raise SystemExit(f"{template.name} is {len(text.encode())} bytes; the limit "
                         f"is {limit:,}")
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
    ap.add_argument("--chunk", type=int, default=400,
                    help="queries per core_nt chunk. Small because spot "
                         "reclaimed two instances in four hours and a reclaim "
                         "loses the chunk it is in, not the run")
    ap.add_argument("--chunk-target-seconds", type=float, default=None,
                    help="how long one chunk should take; default 1500 on spot and "
                         f"{ON_DEMAND_CHUNK_TARGET_SECONDS:.0f} on demand. The instance "
                         "sizes the chunk itself from its own measured scan and marginal "
                         "cost. See CHUNK_TARGET for why the default follows the market")
    ap.add_argument("--nr-chunk", type=int, default=500,
                    help="queries per blastx-against-nr chunk. Small, because nothing in "
                         "this project has ever measured that rate and the first chunk is "
                         "what measures it")
    ap.add_argument("--nr-budget-hours", type=float, default=8.0)
    ap.add_argument("--fast-budget-hours", type=float, default=4.0,
                    help="dc-megablast on the same reads. Cheap by design: if it needs a large budget it is not the cheap task and the comparison says so")
    ap.add_argument("--full-budget-hours", type=float, default=20.0,
                    help="the wall clock the full core_nt pass is allowed to project to "
                         "before the run refuses it and stops at the pilot")
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--no-spot", action="store_true",
                    help="launch on demand rather than spot. No longer gated on the "
                         "market: the gate is the projected cost against the $10 Imo "
                         "pre-approved for a single run. Spot was reclaimed twice in four "
                         "hours, and a reclaim costs whatever chunk it lands in")
    ap.add_argument("--allow-price-change", action="store_true",
                    help="proceed with a run projected over the pre-approved ceiling. Only "
                         "after Imo has said yes to that figure")
    ap.add_argument("--timeout-minutes", type=float, default=100.0)
    ap.add_argument("--poll-seconds", type=float, default=300.0)
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true",
                    help="launch, print the id, stop watching. The instance is then held "
                         "only by its own timer; this session's background jobs cannot "
                         "outlive a run measured in hours, and a watcher killed mid-wait "
                         "would fire its finally and terminate a healthy instance")
    ap.add_argument("--viral-budget-hours", type=float, default=3.0,
                    help="budget for the nt_viruses pass. Smaller than core_nt's because "
                         "the database is 69.8 GB against 318.8, so the scan floor is "
                         "about four minutes rather than eighteen")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

    # THE GATE IS THE COST, NOT THE MARKET. This used to refuse `--no-spot` outright, on the
    # grounds that the $25 Imo approved was a spot figure. That was the wrong variable: the
    # market is not what he was protecting against, the bill is. His standing rule of
    # 2026-11-21 is that a single run under $10 is pre-approved, so the gate is now a
    # projection against that ceiling and on demand needs no special pleading when it is
    # cheap. The projection uses the dead-man timer rather than the stage budgets, because
    # the timer is the only bound that does not depend on the run working.
    cost, cost_note = projected_cost(args.instance_type, args.hard_stop_minutes,
                                     not args.no_spot, args.volume_size)
    print(f"projected worst case: {cost_note}")
    if not (cost == cost):                                   # NaN: unknown instance type
        raise SystemExit(cost_note)
    if cost > AUTO_APPROVED_USD and not args.allow_price_change:
        raise SystemExit(
            f"refusing to launch: {cost_note}, over the ${AUTO_APPROVED_USD:.0f} Imo "
            f"pre-approved for a single run. Above that he wants an explicit yes on a "
            f"figure, so put the number above on a decision card and pass "
            f"--allow-price-change once he has given one. Shortening "
            f"--hard-stop-minutes is the other way down: the projection is the timer's "
            f"worst case, not the expected bill.")

    if args.chunk_target_seconds is None:
        args.chunk_target_seconds = (ON_DEMAND_CHUNK_TARGET_SECONDS if args.no_spot
                                     else SPOT_CHUNK_TARGET_SECONDS)
        print(f"chunk target {args.chunk_target_seconds:.0f}s "
              f"({'on demand, not reclaimed' if args.no_spot else 'spot, reclaim-sized'})")

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
        NR_CHUNK=args.nr_chunk, NR_BUDGET_HOURS=args.nr_budget_hours,
        CHUNK=args.chunk, CHUNK_TARGET_SECONDS=args.chunk_target_seconds,
        FAST_BUDGET_HOURS=args.fast_budget_hours,
        VIRAL_BUDGET_HOURS=args.viral_budget_hours,
    )
    # The run script is content-addressed by the commit AND by a hash of its own rendered
    # text, because two launches off one commit can differ in their stages and budgets, and
    # a key that collided would silently hand the second launch the first one's script. The
    # bucket has no delete, so a collision could not be repaired either.
    run_script = render(TEMPLATE, None, **tokens)
    run_key = (f"{prefix}/code/bigindex-run-{commit[:12]}-"
               f"{hashlib.sha256(run_script.encode()).hexdigest()[:12]}.sh")
    script = render(BOOTSTRAP, 16_384, COMMIT=commit, BUCKET=bucket, REGION=ec2lib.REGION,
                    HARD_STOP_MINUTES=args.hard_stop_minutes, RUN_SCRIPT_KEY=run_key)

    if args.print_user_data:
        print(script)
        print(f"\n### ---- run script, fetched from s3://{bucket}/{run_key} ----\n")
        print(run_script)
        print(f"\n### {len(script.encode())} bytes of UserData, "
              f"{len(run_script.encode()):,} of run script; "
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
        s3.put_object(Bucket=bucket, Key=run_key, Body=run_script.encode())
        print(f"run script at s3://{bucket}/{run_key}", flush=True)
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
    # Per instance, not one file: the first launch's record was overwritten by the second
    # within the hour, and a launch record that only remembers the newest launch cannot
    # answer "what was running when this number was produced".
    record = json.dumps({
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
    }, indent=2)
    (ROOT / f"results/subtract_known/bigindex_launch_{iid}.json").write_text(record)

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
