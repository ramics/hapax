"""Launch the reach comparison: seven search methods, one panel, one instance.

WHAT THIS BUYS THAT THE CONTAINER CANNOT. Nothing here is compute-bound by much - the panel
is 32 contigs and about 50 real ORFs. What the container cannot do is hold the subject sets:
Swiss-Prot plus RefSeq viral as both BLAST volumes and a DIAMOND index, Pfam-A pressed,
RefSeq viral genomes as a nucleotide database, and two Foldseek structure databases, which
is tens of gigabytes against a 15 GB disk whose egress to NCBI and EBI is about 1 MB/s. In
us-east-1 the NCBI databases are an in-region S3 copy and the rest is a CDN download, so the
staging is minutes rather than a day. **The money buys the databases, not the cycles.**

COPIED FROM `ec2_bigindex.py` rather than written fresh, which is what Imo's launch guide
says to do and what that module says about its own ancestor: dry run first, dead man's
switch before anything that can fail, terminate in a `finally`, and the instance pulling
pinned code from S3 rather than being handed a working tree. The differences from that
launcher are three, and all three are because this run is many small things instead of one
big thing:

  * **No chunking and no budget gates.** Those exist because a spot reclaim loses whatever
    block it lands in and a pass that cannot finish spends the money anyway. Here every arm
    is minutes to an hour and uploads on completion, so the unit of loss is already small.
  * **A smaller volume.** 150 GB against 1,700: there is no core_nt and no nr here.
  * **Graded staging.** The run script drops an arm whose database did not stage rather than
    dying, so the table arrives with an honest gap instead of not arriving.

WHY ON DEMAND. Spot was reclaimed twice in four hours on 2026-11-21. At this size the
insurance is not worth the exposure: the whole run is projected at about two dollars on
demand, which is inside the pre-approved ceiling, so there is nothing to save.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ec2_common as ec2lib

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = Path(__file__).resolve().parent / "reach_userdata.sh"
#: The same bootstrap the big-index run uses, token for token. EC2 caps UserData at 16,384
#: bytes and this run script is well over that, so the bootstrap arms the timer and fetches
#: the script from S3. Reused rather than copied: it is the one file whose failure mode is
#: "the machine never terminates", and a second copy of it is a second thing to get wrong.
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

#: Pinned NCBI snapshot, the same one `results/bigindex/` used, so the Swiss-Prot volumes
#: here are the volumes that run searched. `latest-dir` would make a number from this run
#: uncitable in six months.
NCBI_BUCKET = "ncbi-blast-databases"
NCBI_PREFIX = "2026-10-22-01-05-01"

#: ProstT5, already staged in this account's S3. Its presence is the reason the structure
#: arm costs CPU minutes rather than a GPU instance, so the key is pinned here rather than
#: discovered on the instance.
PROSTT5_KEY = "corpora/models/prostt5/b072042c7d5b/prostt5-f16.gguf"

BLAST_VERSION = "2.17.0"
#: The version every timing in `results/diamond_scaling/` was measured with. The baseline
#: arm of this run is the production search, so it has to be the production binary.
DIAMOND_VERSION = "2.1.10"
HMMER_VERSION = "3.4"
HHSUITE_VERSION = "3.3.0"
#: The Steinegger lab's mirror of Pfam as HHM profiles. Pinned to a release rather than a
#: `current` path: the profile-to-profile arm's reach depends on which Pfam it searched, and
#: a floating URL would make that number uncitable.
HH_PFAM_URL = ("https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
               "pfamA_35.0.tar.gz")
#: WHY THIS VERSION AND NOT THE EARLIER PIN, which cost four runs. Release 9-427df8a reads
#: ProstT5 weights only as `model.safetensors` or `model.gguf`, and the only gguf its own
#: downloader ships is a q4_0 quantisation. Handed the f16 gguf this account holds - the
#: right model, the right format, 2.4 GB - it SEGFAULTS inside the model load, which on EC2
#: surfaced as a bare exit 1 with no message. Release 10-941cd33 looks for
#: `prostt5-f16.gguf`, which is what that file is, and writes 3Di from it on CPU. Verified
#: in the container against the real panel before this pin was changed, which is also why
#: the fifth run was the first cheap one: foldseek plus the weights is a 20-second download
#: here, so the arm is debuggable locally and never needed an instance to fail on.
#:
#: Foldseek's own weights distribution is not an alternative: both
#: `prostt5-f16-safetensors.tar.gz` and `prostt5-q4_0-gguf.tar.gz` answered Cloudflare 1027
#: (daily quota exhausted) from two networks on 2026-11-22, and a 429 body is 17 bytes that
#: `foldseek databases ProstT5` untars, fails on, and still exits 0 for.
FOLDSEEK_VERSION = "10-941cd33"

#: Eight vCPU is the account ceiling - every `4xlarge` and larger came back
#: `UnauthorizedOperation`, which reads exactly like having no EC2 permission at all. c6i
#: rather than r6i because nothing here is memory-bound: the largest thing held in memory is
#: the RefSeq viral genome set for the codon arm, a couple of gigabytes.
DEFAULT_INSTANCE_TYPE = "c6i.2xlarge"

#: GB, gp3. Pfam pressed is about 9, the two Foldseek databases about 25, ProstT5 2.4,
#: Swiss-Prot and the viral sets about 5. 150 leaves room for a Foldseek database that has
#: grown since this was measured, and at $0.08/GB-month the headroom costs under a cent an
#: hour.
DEFAULT_VOLUME_GB = 200
DEFAULT_VOLUME_THROUGHPUT = 250
DEFAULT_VOLUME_IOPS = 3000

#: Minutes. The only bound on the bill that does not depend on anything in the run working,
#: including this launcher staying alive. Five hours: staging is projected at under an hour
#: and the arms at two, so this is roughly double and still inside the pre-approved ceiling.
DEFAULT_HARD_STOP_MINUTES = 300
HEARTBEAT_SECONDS = 300

ON_DEMAND_USD_PER_HOUR = {"c6i.xlarge": 0.17, "c6i.2xlarge": 0.34, "m6i.xlarge": 0.19,
                          "m6i.2xlarge": 0.38, "r6i.2xlarge": 0.50}
SPOT_FRACTION = 0.35
#: Imo's standing project rule: a single EC2 run under $10 is pre-approved. Above it, a
#: launch needs his explicit yes on a figure.
AUTO_APPROVED_USD = 10.0

#: Everything the instance runs. Checked by the PRECHECK block in the run script, because
#: the reference project twice added a module to a tool's imports and not to a list like
#: this, and the instance died on an ImportError having read nothing.
PAYLOAD_MODULES = ("experiments/reach_arms.py",)


def projected_cost(instance_type: str, hard_stop_minutes: int, spot: bool,
                   volume_gb: int) -> tuple[float, str]:
    """Worst-case dollars for one run: the dead man's timer is the ceiling on its life.

    The timer rather than the expected run time, because the timer is the only thing that
    bounds the bill without depending on the run working - a crash loop could sit inside a
    stage budget forever.
    """
    hours = hard_stop_minutes / 60.0
    rate = ON_DEMAND_USD_PER_HOUR.get(instance_type)
    if rate is None:
        return float("nan"), (f"{instance_type} is not in ON_DEMAND_USD_PER_HOUR, so no "
                              f"projection is possible; add it rather than guessing")
    compute = rate * (SPOT_FRACTION if spot else 1.0) * hours
    storage = volume_gb * 0.08 * hours / 730.0           # gp3 $0.08/GB-month, prorated
    total = compute + storage
    return total, (f"{instance_type} {'spot' if spot else 'on demand'} for up to "
                   f"{hours:.1f}h (the dead-man timer) is about ${compute:.2f} of compute "
                   f"plus ${storage:.2f} for {volume_gb} GB = **${total:.2f}**")


def render(template: Path, limit: int | None, **tokens: Any) -> str:
    """One template with @@TOKENS@@ filled in. Every token must be consumed.

    An unreplaced token is a bash line that silently does the wrong thing - `--region
    @@REGION@@` would be a literal - and on an instance nobody can log into, that is a
    wasted launch. So this refuses rather than renders a template with a token left in it.
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

    Two flags rather than one. What has to be clean is what the instance will run - the
    payload, this launcher and the run script - because the commit is the only thing that
    makes the instance's tables traceable. The wider tree is almost never clean in a session
    with result files landing in it, and refusing on that would mean committing
    half-finished output to buy a cosmetic check, so it is recorded rather than enforced.
    """
    import subprocess

    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True,
                              text=True, check=True).stdout.strip()

    shipped = [*PAYLOAD_MODULES, "experiments/ec2_reach.py", "experiments/reach_userdata.sh",
               "experiments/bigindex_bootstrap.sh", "experiments/ec2_common.py",
               "experiments/reach_queries.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *shipped)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


def upload_panel(s3: Any, bucket: str, prefix: str, panel_dir: Path, name: str) -> int:
    """Put the panel where the instance will look for it.

    Under `scratch/`, not `results/`: the bucket has no delete, so a query set uploaded
    under a results key would be a permanent object that is not a result. `panel.json` goes
    with the FASTAs so the instance's output prefix carries the panel that produced it.
    """
    sent = 0
    for path in sorted(panel_dir.iterdir()):
        if path.suffix not in {".fa", ".json"}:
            continue
        key = f"{prefix}/scratch/reach-panel-{name}/{path.name}"
        s3.put_object(Bucket=bucket, Key=key, Body=path.read_bytes())
        print(f"  {key}  {path.stat().st_size / 1e3:.1f} kB", flush=True)
        sent += 1
    # Checked AFTER the uploads only because the uploads are idempotent under a fixed key;
    # the point of the check is that the instance cannot be launched against a panel with no
    # test arm or no control, and it names which of the two is missing because those are
    # different mistakes.
    if not (panel_dir / "panel.fa").exists():
        raise SystemExit(f"no panel.fa in {panel_dir}; the instance would start on nothing")
    nulls = sorted(panel_dir.glob("null_s*.fa"))
    if not nulls:
        raise SystemExit(
            f"no null_s*.fa in {panel_dir}; the instance would run every arm with no "
            f"shuffled control, which is the one thing no table here is readable without. "
            f"`reach_queries.py` writes one pile per --null-seeds value")
    print(f"  {len(nulls)} null pile(s): {', '.join(p.name for p in nulls)}", flush=True)
    return sent


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE,
                    choices=ec2lib.ALLOWED_TYPES)
    ap.add_argument("--panel-dir", type=Path,
                    default=Path("/tmp/claude-0/wk/reach_queries"))
    ap.add_argument("--panel-set", default="r1",
                    help="names the scratch prefix the panel goes to. Change it rather "
                         "than overwriting: the bucket has no delete, so a second upload "
                         "under the same name leaves two objects and the instance would "
                         "read whichever won")
    ap.add_argument("--skip-panel-upload", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--volume-size", type=int, default=DEFAULT_VOLUME_GB)
    ap.add_argument("--volume-throughput", type=int, default=DEFAULT_VOLUME_THROUGHPUT)
    ap.add_argument("--volume-iops", type=int, default=DEFAULT_VOLUME_IOPS)
    ap.add_argument("--hard-stop-minutes", type=int, default=DEFAULT_HARD_STOP_MINUTES)
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--spot", action="store_true",
                    help="launch on spot rather than on demand. Off by default: spot was "
                         "reclaimed twice in four hours and this run is cheap enough that "
                         "the insurance is not worth the exposure")
    ap.add_argument("--allow-price-change", action="store_true",
                    help="proceed with a run projected over the pre-approved ceiling. Only "
                         "after Imo has said yes to that figure")
    ap.add_argument("--timeout-minutes", type=float, default=100.0)
    ap.add_argument("--poll-seconds", type=float, default=300.0)
    ap.add_argument("--attach", help="watch an already-running instance")
    ap.add_argument("--detach", action="store_true",
                    help="launch, print the id, stop watching. The instance is then held "
                         "only by its own timer; a watcher killed mid-wait would fire its "
                         "finally and terminate a healthy instance")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="authorize only, launch nothing")
    ap.add_argument("--print-user-data", action="store_true")
    args = ap.parse_args(argv)

    cost, cost_note = projected_cost(args.instance_type, args.hard_stop_minutes,
                                     args.spot, args.volume_size)
    if args.attach:
        # --attach launches nothing, so it must not be gated: the instance it watches was
        # approved when it was launched, and refusing to WATCH a running box would leave it
        # billing unobserved, which is the opposite of what a spend gate is for.
        cost, cost_note = 0.0, f"watching {args.attach}; this launches nothing"
    print(f"projected worst case: {cost_note}")
    if not (cost == cost):                                   # NaN: unknown instance type
        raise SystemExit(cost_note)
    if cost > AUTO_APPROVED_USD and not args.allow_price_change:
        raise SystemExit(
            f"refusing to launch: {cost_note}, over the ${AUTO_APPROVED_USD:.0f} Imo "
            f"pre-approved for a single run. Put the figure above in front of him and "
            f"pass --allow-price-change once he has said yes. Shortening "
            f"--hard-stop-minutes is the other way down: the projection is the timer's "
            f"worst case, not the expected bill.")

    commit, shipped_dirty, tree_dirty = git_head(ROOT)
    if shipped_dirty and not args.allow_dirty:
        raise SystemExit(
            "refusing to launch: the code that would run on the instance is not committed. "
            "The instance has no checkout, so the commit is the only thing that makes its "
            "tables traceable. Commit the launcher, the run script and the payload first.")

    bucket, prefix = ec2lib.BUCKET, ec2lib.PREFIX
    code_key = f"{prefix}/code/reach-{commit}.tar.gz"
    tokens = dict(
        COMMIT=commit, BUCKET=bucket, PREFIX=prefix, REGION=ec2lib.REGION,
        HARD_STOP_MINUTES=args.hard_stop_minutes, HEARTBEAT_SECONDS=HEARTBEAT_SECONDS,
        THREADS=args.threads, SEED=args.seed, PANEL_SET=args.panel_set,
        CODE_KEY=Path(code_key).name, INSTANCE_TYPE=args.instance_type,
        BLAST_VERSION=BLAST_VERSION, DIAMOND_VERSION=DIAMOND_VERSION,
        HMMER_VERSION=HMMER_VERSION, FOLDSEEK_VERSION=FOLDSEEK_VERSION,
        HHSUITE_VERSION=HHSUITE_VERSION, HH_PFAM_URL=HH_PFAM_URL,
        NCBI_BUCKET=NCBI_BUCKET, NCBI_PREFIX=NCBI_PREFIX, PROSTT5_KEY=PROSTT5_KEY,
    )
    # The run script is content-addressed by the commit AND by a hash of its own rendered
    # text, because two launches off one commit can differ in their parameters, and a key
    # that collided would silently hand the second launch the first one's script. The bucket
    # has no delete, so a collision could not be repaired either.
    run_script = render(TEMPLATE, None, **tokens)
    run_key = (f"{prefix}/code/reach-run-{commit[:12]}-"
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
    print(f"commit {commit[:12]}, panel set {args.panel_set}, seed {args.seed}", flush=True)

    if args.attach:
        return watch(ec2, args.attach, args)

    # Not on a dry run: a dry run authorizes and executes nothing, and this bucket has no
    # delete, so an upload here would leave permanent objects behind for a launch that never
    # happened.
    if not args.dry_run:
        s3.put_object(Bucket=bucket, Key=code_key, Body=build_payload(ROOT))
        print(f"payload at s3://{bucket}/{code_key}", flush=True)
        s3.put_object(Bucket=bucket, Key=run_key, Body=run_script.encode())
        print(f"run script at s3://{bucket}/{run_key}", flush=True)
        if not args.skip_panel_upload:
            n = upload_panel(s3, bucket, prefix, args.panel_dir, args.panel_set)
            print(f"{n} panel file(s) uploaded as set {args.panel_set}", flush=True)

    ssm = boto3.client("ssm", region_name=ec2lib.REGION)
    ami = ssm.get_parameter(Name=ec2lib.AMI_PARAMETER)["Parameter"]["Value"]
    market = ec2lib.market_options(args.spot)

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
                                     {"Key": "Task", "Value": "hapax-reach"}]}],
        UserData=script,
    )
    if market is not None:
        kwargs["InstanceMarketOptions"] = market

    if args.dry_run:
        try:
            ec2.run_instances(DryRun=True, **kwargs)
        except Exception as exc:                 # DryRunOperation is the success case
            print(f"dry run: {exc}")
        return 0

    try:
        reply = ec2.run_instances(**kwargs)
    except Exception as exc:
        message = ec2lib.spot_refusal(str(exc)) if args.spot else None
        if message:
            raise SystemExit(message) from exc
        raise
    iid = reply["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)

    # Per instance, not one file: a launch record that only remembers the newest launch
    # cannot answer "what was running when this number was produced".
    out = ROOT / "results/reach"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"launch_{iid}.json").write_text(json.dumps({
        "instance_id": iid, "commit": commit, "shipped_dirty": shipped_dirty,
        "wider_tree_dirty": tree_dirty, "instance_type": args.instance_type,
        "panel_set": args.panel_set, "seed": args.seed, "threads": args.threads,
        "volume_gb": args.volume_size, "spot": args.spot,
        "hard_stop_minutes": args.hard_stop_minutes,
        "projected_worst_case_usd": round(cost, 2),
        "ncbi_prefix": NCBI_PREFIX, "blast": BLAST_VERSION, "diamond": DIAMOND_VERSION,
        "hmmer": HMMER_VERSION, "foldseek": FOLDSEEK_VERSION, "prostt5_key": PROSTT5_KEY,
        "hhsuite": HHSUITE_VERSION, "hh_pfam": HH_PFAM_URL,
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/reach/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/reach-{iid}/",
        "run_script": f"s3://{bucket}/{run_key}",
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
