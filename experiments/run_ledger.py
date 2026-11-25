"""The cost and wall clock of every launch behind one production run, generated.

WHY THIS IS A SCRIPT. The previous library's report carries this table typed by hand, one
row per launch with a wall clock and a dollar figure read off a console. Two of its rows
are estimates marked with a tilde and none of them can be re-derived, which is the shape
`CLAUDE.md` rules out for figures and which a cost table has no more right to than a
figure: the spend is the number Imo approves against, so it needs provenance like any
other.

WHERE EACH NUMBER COMES FROM. The launch time, instance type, commit and projected spend
are in the launch manifest the launcher commits. The end of a launch is NOT recoverable
from EC2: a terminated instance leaves `describe_instances` within hours, and every one of
these had already gone when this ran. What survives is the instance's own writes to S3,
under its results prefix and its log prefix, so the last write is the floor on when it
stopped and the figure here is a LOWER BOUND on wall clock and therefore on cost. The
bound is stated rather than padded, because a guess at the gap between the last write and
the terminate call would be the only unprovenanced number in the table.

THE RATE IS A CONSTANT WITH A DATE. On-demand pricing changes, so the rate used is written
down with where it came from rather than queried: a table regenerated next year must not
silently re-price last year's run.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "experiments"))

import ec2_common  # noqa: E402
from provenance import provenance  # noqa: E402

#: USD per hour, us-east-1 on-demand, as of 2026-11-25. Written down rather than queried
#: so that regenerating this table later cannot re-price an older run.
RATES_USD_PER_HOUR = {"r6i.2xlarge": 0.504, "r6i.4xlarge": 1.008,
                      "r6i.8xlarge": 2.016, "m6i.2xlarge": 0.384}
RATES_AS_OF = "2026-11-25"

#: gp3, USD per GB-month, same source and date. Charged for the life of the instance
#: because these volumes are delete-on-termination.
EBS_USD_PER_GB_MONTH = 0.08


class Unreachable(SystemExit):
    """The bucket could not be listed, which is not the same as a launch writing nothing.

    CI runs this suite with no AWS credentials, so the first version of the refusal below
    never fired: botocore raised `NoCredentialsError` out of the paginator and the script
    exited on a traceback. An environment that cannot reach the bucket and an instance
    that wrote nothing to it produce the same empty listing and are different facts, so
    they get different exits.
    """


def _window(s3, bucket: str, prefixes: list[str]):
    """First and last write under any of `prefixes`, and how many keys there are."""
    first = last = None
    n = 0
    for pre in prefixes:
        try:
            pages = list(s3.get_paginator("list_objects_v2").paginate(Bucket=bucket,
                                                                      Prefix=pre))
        except Exception as e:  # botocore raises several unrelated types here
            raise Unreachable(
                f"cannot list s3://{bucket}/{pre}, so no wall clock can be measured "
                f"and no cost can be stated: {type(e).__name__}: {e}") from e
        for page in pages:
            for o in page.get("Contents", []):
                n += 1
                t = o["LastModified"]
                first = t if first is None or t < first else first
                last = t if last is None or t > last else last
    return first, last, n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", required=True)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--kind", default="oneshot")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    res = a.results or ROOT / f"results/prod_{a.library}"
    launches = sorted(res.glob("launch_i-*.json"))
    if not launches:
        raise SystemExit(f"{res} holds no launch_i-*.json, so there is no ledger to write")

    try:
        s3 = ec2_common.s3_client()
    except Exception as e:
        # Building the client fails before any listing when there is no profile or
        # credential to build it from, so this exit has to sit here as well.
        raise Unreachable(
            f"cannot open an S3 client, so no wall clock can be measured and no cost "
            f"can be stated: {type(e).__name__}: {e}") from e
    rows = []
    for f in launches:
        d = json.loads(f.read_text())
        iid = f.name.removeprefix("launch_").removesuffix(".json")
        # PREFIX has no trailing slash, and the first version of this joined it without
        # one: every listing came back empty and the table rendered a column of dashes
        # and a $0.00 total, which is the shape this repo already withdrew a claim over.
        # An empty listing is now a refusal below rather than a row of its own.
        base = ec2_common.PREFIX.rstrip("/")
        first, last, keys = _window(s3, ec2_common.BUCKET, [
            f"{base}/results/{a.kind}/{a.library}/{iid}/",
            f"{base}/logs/{a.kind}-{a.library}-{iid}/"])
        itype = d.get("instance_type")
        rate = RATES_USD_PER_HOUR.get(itype)
        hours = None
        if last is not None and d.get("launched"):
            from datetime import datetime
            t0 = datetime.fromisoformat(d["launched"].replace("Z", "+00:00"))
            hours = (last - t0).total_seconds() / 3600
        cost = None
        if hours is not None and rate is not None:
            cost = hours * (rate + (d.get("volume_gb") or 0)
                            * EBS_USD_PER_GB_MONTH / 730)
        rows.append({"instance": iid, "commit": (d.get("commit") or "")[:7],
                     "instance_type": itype, "launched": d.get("launched"),
                     "last_write": last.isoformat() if last else None,
                     "keys_written": keys, "hours_to_last_write": hours,
                     "projected_usd": d.get("projected_usd"),
                     "usd_to_last_write": cost,
                     "rate_usd_per_hour": rate, "volume_gb": d.get("volume_gb")})
    rows.sort(key=lambda r: r["launched"] or "")
    # A LEDGER WITH NO TIMINGS IS NOT A CHEAP LEDGER. If nothing was found for any
    # launch the prefix is wrong or the bucket is unreachable, and a table of dashes
    # totalling $0.00 reads like a free run rather than a failed measurement.
    if not any(r["keys_written"] for r in rows):
        raise SystemExit(
            f"no S3 key was found for any of the {len(rows)} launches under "
            f"{ec2_common.BUCKET}/{ec2_common.PREFIX.rstrip('/')}/{{results,logs}}/"
            f"{a.kind}*{a.library}*, so no wall clock can be measured. A ledger of "
            "dashes would read as a free run.")
    silent = [r["instance"] for r in rows if not r["keys_written"]]

    L = [f"# Launch ledger — {a.library}", "",
         f"Generated by `experiments/run_ledger.py`. Rates are {RATES_AS_OF} us-east-1 "
         f"on-demand, written into the script rather than queried, so regenerating this "
         f"table cannot re-price an older run. **Wall clock is measured to each "
         f"instance's last write to S3**, which is a lower bound: a terminated instance "
         f"leaves `describe_instances` within hours and every one of these had already "
         f"gone, so the gap between the last write and the terminate call is not "
         f"recoverable and is not estimated here.", "",
         "| instance | commit | launched | to last write | keys | cost to last write | "
         "projected |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        h = f"{r['hours_to_last_write']:.2f} h" if r["hours_to_last_write"] else "-"
        c = f"${r['usd_to_last_write']:.2f}" if r["usd_to_last_write"] else "-"
        pj = f"${r['projected_usd']:.2f}" if r.get("projected_usd") else "-"
        L += [f"| `{r['instance']}` | `{r['commit']}` | {r['launched']} | {h} | "
              f"{r['keys_written']} | {c} | {pj} |"]
    total = sum(r["usd_to_last_write"] or 0 for r in rows)
    L += ["", f"**At least ${total:.2f} across {len(rows)} launches.** Every figure is a "
              f"lower bound for the reason above, and the last row is the pass that "
              f"produced the committed result."]
    if silent:
        L += ["", f"{len(silent)} launch{'' if len(silent) == 1 else 'es'} wrote nothing "
                  f"to S3 and so contribute no measured time or cost: "
                  + ", ".join(f"`{i}`" for i in silent)
                  + ". Each died before its first upload, which bounds it by the "
                    "interval to the next launch rather than by anything it recorded."]
    L += [""]

    out = a.out or res / "LEDGER.md"
    out.write_text("\n".join(L) + "\n")
    (out.parent / "ledger.json").write_text(json.dumps(
        {"library": a.library, "rows": rows,
         "rates_usd_per_hour": RATES_USD_PER_HOUR, "rates_as_of": RATES_AS_OF,
         "ebs_usd_per_gb_month": EBS_USD_PER_GB_MONTH,
         "wall_clock_basis": "launch time to the instance's last write to S3; a lower "
                             "bound, because a terminated instance is purged from "
                             "describe_instances before this can read it",
         **provenance(__file__)}, indent=2) + "\n")
    print(f"written {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
