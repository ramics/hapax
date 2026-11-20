"""Launch-time facts about this account's EC2, and the one way to build an S3 client.

COPIED, DELIBERATELY, from `ramics/evo-training-check` (`src/amazebio/ec2.py`, commit
c6cb64a) because Imo's launch guide says to copy that launcher rather than write one fresh:
every function here encodes a refusal that cost somebody a launch, and a rewrite from the
guide's prose would lose the half of each lesson that is in the wording. What is NOT copied
is the parts specific to that project's corpora.

The one addition is `s3_client`. The guide calls it out as the mistake other projects will
make: a bare `boto3.client("s3", region_name="us-east-1")` resolves to the legacy global
hostname `s3.amazonaws.com`, which the egress proxy does not sign, and S3 answers an
unsigned caller with a bare `AccessDenied` that looks exactly like a missing IAM grant. The
same trap caught this thread's first attempt to list the NCBI public bucket from the
container, from the other side: there, the signed regional path fails and the UNSIGNED
global one works, because the bucket is public and our principal has no grant on it.
"""

from __future__ import annotations

from typing import Any

REGION = "us-east-1"
AMI_PARAMETER = "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
S3_ENDPOINT = "https://s3.us-east-1.amazonaws.com"

#: The shared bucket. No `s3:DeleteObject` exists on it, so every write is permanent and the
#: keys a run will write are fixed before launch (guide, "The shared bucket").
BUCKET = "amazebio-corpora-897641797669"

#: Our prefix inside it. `corpora/` belongs to amazebio and is left alone.
PREFIX = "projects/hapax"

#: What the account's policy permits, measured by dry run in the reference project on
#: 2026-10-20: the ceiling is eight vCPU, and every `4xlarge` and larger came back
#: `UnauthorizedOperation`, which reads exactly like having no EC2 permission at all.
ALLOWED_TYPES = ("c6i.xlarge", "c6i.2xlarge", "m6i.xlarge", "m6i.2xlarge", "r6i.2xlarge")

#: vCPUs per permitted type, so a launcher can check its own thread arithmetic. The
#: reference project lost a fourteen-hour ceiling to eight workers asking for three threads
#: each on eight cores, with nothing in the launcher able to notice.
VCPUS = {"c6i.xlarge": 4, "c6i.2xlarge": 8, "m6i.xlarge": 4, "m6i.2xlarge": 8,
         "r6i.2xlarge": 8}


def s3_client(*, config: Any = None) -> Any:
    """An S3 client pinned to the regional endpoint. Build them nowhere else.

    `region_name` alone is not enough; see the module docstring for why the failure is
    misleading rather than loud.
    """
    import boto3

    return boto3.client("s3", region_name=REGION, endpoint_url=S3_ENDPOINT, config=config)


def spot_refusal(text: str) -> str | None:
    """Why a spot launch was refused, in terms the next person can act on. `None` if unclear.

    All three refusals this distinguishes mean nothing launched and nothing was billed, and
    they need different responses - which is why they are not one message. The two that will
    not clear by retrying are the expensive ones to misread: a retry loop against a quota
    cannot win, and a capacity refusal sent to an IAM page wastes somebody's afternoon.
    """
    if "ServiceLinkedRoleCreationNotPermitted" in text:
        return (
            f"spot is unavailable in this account: {text[:160]}\n"
            "  the EC2 Spot service-linked role does not exist and iam:* is denied, so it "
            "cannot be created from here.\n"
            "  nothing launched, nothing billed. Either someone with IAM permissions "
            "creates AWSServiceRoleForEC2Spot once, or the run goes on demand at roughly "
            "3-4x the price - which needs Imo's approval again, because the figure he "
            "approved is a spot figure."
        )
    if "MaxSpotInstanceCountExceeded" in text or "Max spot instance count" in text:
        return (
            f"spot refused by an account limit: {text[:160]}\n"
            "  nothing launched, nothing billed, and nothing in the account is holding it - "
            "the reference project saw this with zero spot instances and zero open "
            "requests.\n"
            "  it is a quota, so retrying will not clear it. Either someone raises the spot "
            "limit for this instance family, or the run goes on demand - which is a price "
            "change and therefore a new approval."
        )
    if "Capacity" in text or "capacity" in text:
        return (
            f"spot capacity refused: {text[:200]}\n"
            "  nothing launched and nothing billed. Retry, try another permitted instance "
            "type, or go on demand at the higher price."
        )
    return None


def market_options(spot: bool) -> dict[str, Any] | None:
    """What to ask EC2 for, market-wise. `None` means say nothing and get on demand.

    Three choices, each with a cheaper-looking alternative that is wrong here. ONE-TIME
    rather than persistent: a reclaim on a persistent request re-launches the instance and
    restarts a multi-hour pass with nobody watching, and bills for it. TERMINATE rather than
    stop or hibernate on interruption: stop keeps billing the volume - 1.6 TB of it here,
    which is most of this run's cost - and leaves a resource to forget. NO MaxPrice:
    omitting it caps the bid at on demand, and a bid under the going rate converts a run
    that would have cost money into a launch that silently never happens.
    """
    if not spot:
        return None
    return {"MarketType": "spot",
            "SpotOptions": {"SpotInstanceType": "one-time",
                            "InstanceInterruptionBehavior": "terminate"}}


def transition_reason(ec2: Any, iid: str) -> str:
    """Why EC2 says the instance changed state. Empty when it does not say.

    Worth reading rather than inferring, and worth knowing what it cannot tell you: a
    boot-armed shutdown and a manual termination both report `User initiated`, while a spot
    reclaim reports `Service initiated`. So this separates a reclaim from a shutdown and does
    not separate a shutdown from a human.
    """
    try:
        inst = ec2.describe_instances(InstanceIds=[iid])["Reservations"][0]["Instances"][0]
        return str(inst.get("StateTransitionReason", ""))
    except Exception:
        return ""


def state_of(ec2: Any, iid: str) -> str:
    """The instance's state name, or `"unknown"` when the call fails.

    Never raises: this is read on the way out of a launch, including from a `finally`, and a
    transient API failure there must not mask the outcome it is reporting on.
    """
    try:
        return ec2.describe_instances(InstanceIds=[iid])["Reservations"][0][
            "Instances"][0]["State"]["Name"]
    except Exception:
        return "unknown"
