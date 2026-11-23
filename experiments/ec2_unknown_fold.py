"""Launch the structure arm on the unnameable pile: one instance, one arm, cycles not databases.

WHAT THIS RUN BUYS, AND WHY IT IS THE FIRST OF ITS KIND HERE. Every earlier launcher in this
repo bought *databases*: core_nt at 319 GB, Swiss-Prot plus RefSeq viral as BLAST volumes,
Pfam pressed - tens of gigabytes that a 15 GB container with 1 MB/s egress cannot hold. This
run stages one 2.4 GB model and downloads a public CDN archive, both of which the container
can do. What the container cannot do is the arithmetic: ProstT5 encodes about 34 residues a
second on its four contended cores, measured here on this run's own ORFs, and this run has
of the order of half a million residues to encode. **The money buys the cycles.**

WHY r6i RATHER THAN c6i, WHICH EVERY OTHER LAUNCHER USES. Foldseek's own threading goes
badly sublinear - the reach arm measured 49.6 residues/s on four cores and 66.3 on eight,
1.34x for twice the cores - so asking for eight threads wastes most of them. The encoding is
split across independent processes instead, and each one loads its own 2.4 GB copy of the
model. Eight shards is 20 GB of resident model, which does not fit c6i.2xlarge's 16 GB and
fits r6i.2xlarge's 64 GB with room to spare. The 16 cents an hour buys roughly twice the
throughput, so it is cheaper per residue and not only faster.

COPIED FROM `ec2_reach.py` rather than written fresh, which is what Imo's launch guide says
to do: dry run first, dead man's switch before anything that can fail, terminate in a
`finally`, the instance pulling pinned code from S3 rather than being handed a working tree,
and a launch record per instance. Two things are deliberately simpler than that ancestor.
There is no graded staging - that run dropped an arm whose database failed and still
produced a six-arm table, where this run is one arm and a staging failure has nothing to
fall back to, so every stage is fatal. And there is no chunking or budget gate, because the
whole run is hours rather than days.

WHY ON DEMAND. Spot was reclaimed twice in four hours on 2026-11-21. The whole run is
projected at a couple of dollars on demand, inside the pre-approved ceiling, so there is
nothing to save and a reclaim would cost the encoding done so far.
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
TEMPLATE = Path(__file__).resolve().parent / "unknown_fold_userdata.sh"
#: The same bootstrap every run here uses, token for token. EC2 caps UserData at 16,384
#: bytes; the run script is over it, so the bootstrap arms the timer and fetches the script
#: from S3. Reused rather than copied: it is the one file whose failure mode is "the machine
#: never terminates", and a second copy of it is a second thing to get wrong.
BOOTSTRAP = Path(__file__).resolve().parent / "bigindex_bootstrap.sh"

#: ProstT5, already staged in this account's S3. Its presence is why this arm costs CPU
#: hours rather than a GPU instance and a folding run, so the key is pinned here rather than
#: discovered on the instance.
PROSTT5_KEY = "corpora/models/prostt5/b072042c7d5b/prostt5-f16.gguf"

#: WHY THIS VERSION AND NOT THE EARLIER PIN, which cost the reach run four launches. Release
#: 9-427df8a reads ProstT5 weights only as `model.safetensors` or `model.gguf`, and handed
#: the f16 gguf this account holds it SEGFAULTS inside the model load, which on EC2 surfaced
#: as a bare exit 1 with no message. 10-941cd33 reads `prostt5-f16.gguf`, which is what that
#: file is. The run script checks the binary reports this commit before it encodes anything.
FOLDSEEK_VERSION = "10-941cd33"
FOLDSEEK_COMMIT = "941cd33ff0771cd2e3f144e3293e22a2b87e9fda"

#: Eight vCPU is the account ceiling: every `4xlarge` and larger came back
#: `UnauthorizedOperation`, which reads exactly like having no EC2 permission at all.
DEFAULT_INSTANCE_TYPE = "r6i.2xlarge"
#: One per vCPU. Each is a separate process with its own copy of the model, so this is
#: bounded by the instance's memory, which is why the default type is the 64 GB one.
DEFAULT_SHARDS = 8

#: GB, gp3. The subject library is 40,066 mmCIF files at about 20 GB plus a Foldseek
#: database of the same order, ProstT5 2.4, and the piles are megabytes. 150 leaves room for
#: a databank that has grown since this was measured.
DEFAULT_VOLUME_GB = 150
DEFAULT_VOLUME_THROUGHPUT = 250
DEFAULT_VOLUME_IOPS = 3000

#: Minutes. The only bound on the bill that does not depend on anything in the run working,
#: including this launcher staying alive. Five hours: the library is about 40 minutes of
#: downloading plus a quarter of an hour to format, the encoding is projected at under two
#: hours at eight shards, and the search is minutes. Better than double the projection,
#: which is the margin the earlier runs found they needed.
DEFAULT_HARD_STOP_MINUTES = 300
HEARTBEAT_SECONDS = 300

ON_DEMAND_USD_PER_HOUR = {"c6i.xlarge": 0.17, "c6i.2xlarge": 0.34, "m6i.xlarge": 0.19,
                          "m6i.2xlarge": 0.38, "r6i.2xlarge": 0.50}
SPOT_FRACTION = 0.35
#: Imo's standing project rule: a single EC2 run under $10 is pre-approved. Above it, a
#: launch needs his explicit yes on a figure.
AUTO_APPROVED_USD = 10.0

#: Everything the instance runs. Checked by the PRECHECK block in the run script, because
#: this project twice added a module to a script's imports and not to a list like this, and
#: the instance died on an ImportError having read nothing.
PAYLOAD_MODULES = ("experiments/unknown_fold.py", "experiments/shuffle_control.py",
                   "experiments/provenance.py")

#: The contig sets the cascade itself produced, already in this account's S3 under the
#: shared prefix. Copied to a scratch prefix per pile set rather than read from there
#: directly, so the instance's output prefix records exactly which files it read.
PILE_SOURCE_PREFIX = "shared/contig_sets"
PILE_FILES = ("sample_contigs_unidentified.fa", "sample_contigs_identified_norrna.fa")


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
    with result files landing in it, and refusing on that would mean committing half-finished
    output to buy a cosmetic check, so it is recorded rather than enforced.
    """
    import subprocess

    def git(*a: str) -> str:
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True,
                              text=True, check=True).stdout.strip()

    shipped = [*PAYLOAD_MODULES, "experiments/ec2_unknown_fold.py",
               "experiments/unknown_fold_userdata.sh", "experiments/bigindex_bootstrap.sh",
               "experiments/ec2_common.py"]
    try:
        return (git("rev-parse", "HEAD"),
                bool(git("status", "--porcelain", "--", *shipped)),
                bool(git("status", "--porcelain")))
    except Exception:
        return "unknown", True, True


def copy_piles(s3: Any, bucket: str, prefix: str, pile_set: str) -> list[str]:
    """Put the two contig piles where the instance will look for them.

    SERVER-SIDE COPY, not a download and re-upload: these are 32 MB that are already in this
    bucket, and the container's egress is slow enough that round-tripping them would be the
    longest step in the launch. Under `scratch/`, not `results/`: the bucket has no delete,
    so a query set uploaded under a results key would be a permanent object that is not a
    result.
    """
    keys = []
    for name in PILE_FILES:
        src = f"{prefix}/{PILE_SOURCE_PREFIX}/{name}"
        dst = f"{prefix}/scratch/unknown-fold-piles-{pile_set}/{name}"
        head = s3.head_object(Bucket=bucket, Key=src)
        s3.copy_object(Bucket=bucket, Key=dst, CopySource={"Bucket": bucket, "Key": src})
        print(f"  {dst}  {head['ContentLength'] / 1e6:.1f} MB", flush=True)
        keys.append(dst)
    return keys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--contigs-per-seed", type=int, default=400,
                    help="contigs drawn per seed per pile, as five DISJOINT blocks. 400 is "
                         "about 330 ORFs a seed from the unnameable pile, which measures a "
                         "placement rate to about a point and a half")
    ap.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE,
                    choices=ec2lib.ALLOWED_TYPES)
    ap.add_argument("--shards", type=int, default=DEFAULT_SHARDS)
    ap.add_argument("--pile-set", default="p1",
                    help="names the scratch prefix the piles go to. Change it rather than "
                         "overwriting: the bucket has no delete, so a second copy under the "
                         "same name leaves two objects and the instance would read "
                         "whichever won")
    ap.add_argument("--skip-pile-copy", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--volume-size", type=int, default=DEFAULT_VOLUME_GB)
    ap.add_argument("--volume-throughput", type=int, default=DEFAULT_VOLUME_THROUGHPUT)
    ap.add_argument("--volume-iops", type=int, default=DEFAULT_VOLUME_IOPS)
    ap.add_argument("--hard-stop-minutes", type=int, default=DEFAULT_HARD_STOP_MINUTES)
    ap.add_argument("--profile", default="AmazebioStager")
    ap.add_argument("--owner", default="hapax")
    ap.add_argument("--spot", action="store_true",
                    help="launch on spot rather than on demand. Off by default: spot was "
                         "reclaimed twice in four hours and a reclaim here loses the "
                         "encoding done so far")
    ap.add_argument("--allow-price-change", action="store_true",
                    help="proceed with a run projected over the pre-approved ceiling. Only "
                         "after Imo has said yes to that figure")
    ap.add_argument("--timeout-minutes", type=float, default=330.0)
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

    if args.shards > 8 and args.instance_type != "r6i.2xlarge":
        raise SystemExit(
            f"{args.shards} shards is {args.shards * 2.4:.0f} GB of resident model, which "
            f"only r6i.2xlarge has room for. An OOM-killed shard reads as structure "
            f"placing fewer genes, so this refuses rather than letting the kernel decide.")

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
            f"pre-approved for a single run. Put the figure above in front of him and pass "
            f"--allow-price-change once he has said yes. Shortening --hard-stop-minutes is "
            f"the other way down: the projection is the timer's worst case, not the "
            f"expected bill.")

    commit, shipped_dirty, tree_dirty = git_head(ROOT)
    if shipped_dirty and not args.allow_dirty:
        raise SystemExit(
            "refusing to launch: the code that would run on the instance is not committed. "
            "The instance has no checkout, so the commit is the only thing that makes its "
            "tables traceable. Commit the launcher, the run script and the payload first.")

    bucket, prefix = ec2lib.BUCKET, ec2lib.PREFIX
    code_key = f"{prefix}/code/unknown-fold-{commit}.tar.gz"
    tokens = dict(
        COMMIT=commit, BUCKET=bucket, PREFIX=prefix, REGION=ec2lib.REGION,
        HARD_STOP_MINUTES=args.hard_stop_minutes, HEARTBEAT_SECONDS=HEARTBEAT_SECONDS,
        THREADS=args.threads, SHARDS=args.shards, SEEDS=args.seeds,
        CONTIGS_PER_SEED=args.contigs_per_seed, PILE_SET=args.pile_set,
        CODE_KEY=Path(code_key).name, FOLDSEEK_VERSION=FOLDSEEK_VERSION,
        FOLDSEEK_COMMIT=FOLDSEEK_COMMIT, PROSTT5_KEY=PROSTT5_KEY,
    )
    # The run script is content-addressed by the commit AND by a hash of its own rendered
    # text, because two launches off one commit can differ in their parameters, and a key
    # that collided would silently hand the second launch the first one's script. The bucket
    # has no delete, so a collision could not be repaired either.
    run_script = render(TEMPLATE, None, **tokens)
    run_key = (f"{prefix}/code/unknown-fold-run-{commit[:12]}-"
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
    print(f"commit {commit[:12]}, pile set {args.pile_set}, seeds {args.seeds}, "
          f"{args.contigs_per_seed} contigs a seed, {args.shards} shards", flush=True)

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
        if not args.skip_pile_copy:
            copy_piles(s3, bucket, prefix, args.pile_set)

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
                                     {"Key": "Task", "Value": "hapax-unknown-fold"}]}],
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
    out = ROOT / "results/unknown_fold"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"launch_{iid}.json").write_text(json.dumps({
        "instance_id": iid, "commit": commit, "shipped_dirty": shipped_dirty,
        "wider_tree_dirty": tree_dirty, "instance_type": args.instance_type,
        "pile_set": args.pile_set, "seeds": args.seeds,
        "contigs_per_seed": args.contigs_per_seed, "threads": args.threads,
        "shards": args.shards, "volume_gb": args.volume_size, "spot": args.spot,
        "hard_stop_minutes": args.hard_stop_minutes,
        "projected_worst_case_usd": round(cost, 2),
        "foldseek": FOLDSEEK_VERSION, "prostt5_key": PROSTT5_KEY,
        "launched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results_prefix": f"s3://{bucket}/{prefix}/results/unknown_fold/{iid}/",
        "logs_prefix": f"s3://{bucket}/{prefix}/logs/unknown_fold-{iid}/",
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
