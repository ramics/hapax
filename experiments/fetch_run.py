"""Copy one production run's outputs out of S3 into a results directory.

WHY THIS IS A SCRIPT. Every committed production directory in this repo arrived by
somebody typing boto3 calls at a prompt once, which makes the provenance of the files
themselves unrecorded: nothing says which instance, which prefix, or what was left behind.
A run's outputs are the evidence under a REPORT.md, so how they got here is part of the
record, and `fetch_manifest.json` beside them says it.

WHAT IT REFUSES. A local file that already exists with different content is not
overwritten: a results directory holding two instances' idea of `tiers.json` cannot be
read, and the one that would win is whichever ran last. A file over the size cap is
skipped and NAMED, because a quiet omission and a run that never wrote it look the same
afterwards.

    python3 -I experiments/fetch_run.py --instance i-0123 --library SRR34217566
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ec2_common as ec2lib                                      # noqa: E402
import provenance                                                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--instance", required=True, help="the instance id the run ran on")
    ap.add_argument("--library", required=True)
    ap.add_argument("--kind", default="oneshot", choices=("oneshot", "prod"),
                    help="which template wrote it: the two use separate prefixes, and "
                         "the bucket has no delete, so they are never merged")
    ap.add_argument("--out", type=Path, default=None,
                    help="default results/prod_<library>/")
    ap.add_argument("--max-mb", type=float, default=50.0,
                    help="skip a file larger than this and name it in the manifest")
    ap.add_argument("--include-logs", action="store_true",
                    help="also copy the instance's own logs, which are not results")
    # Present so the convention holds and so a reader can see it changes nothing: this
    # copies bytes and draws no sample.
    ap.add_argument("--seed", type=int, default=1,
                    help="recorded only; nothing here samples anything")
    a = ap.parse_args()

    out = a.out or ROOT / f"results/prod_{a.library}"
    out.mkdir(parents=True, exist_ok=True)
    s3 = ec2lib.s3_client()
    prefixes = [f"{ec2lib.PREFIX}/results/{a.kind}/{a.library}/{a.instance}/"]
    if a.include_logs:
        prefixes.append(f"{ec2lib.PREFIX}/logs/{a.kind}-{a.library}-{a.instance}/")

    copied, skipped, differed = [], [], []
    for prefix in prefixes:
        token = None
        while True:
            kw = {"Bucket": ec2lib.BUCKET, "Prefix": prefix}
            if token:
                kw["ContinuationToken"] = token
            page = s3.list_objects_v2(**kw)
            for obj in page.get("Contents", []):
                rel = obj["Key"][len(prefix):]
                if not rel:
                    continue
                if obj["Size"] > a.max_mb * 1e6:
                    skipped.append({"key": rel, "bytes": obj["Size"],
                                    "why": f"over --max-mb {a.max_mb}"})
                    continue
                body = s3.get_object(Bucket=ec2lib.BUCKET,
                                     Key=obj["Key"])["Body"].read()
                dest = out / rel
                if dest.exists() and _sha(dest.read_bytes()) != _sha(body):
                    differed.append({"key": rel, "local": str(dest)})
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(body)
                copied.append({"key": rel, "bytes": len(body),
                               "sha256_16": _sha(body)[:16]})
            if not page.get("IsTruncated"):
                break
            token = page.get("NextContinuationToken")

    manifest_path = out / "fetch_manifest_local.json"
    manifest_path.write_text(json.dumps({
        "script": "experiments/fetch_run.py",
        **provenance.provenance(__file__),
        "parameters": {k: (str(v) if isinstance(v, Path) else v)
                       for k, v in vars(a).items()},
        "instance": a.instance, "library": a.library, "prefixes": prefixes,
        "copied": copied, "skipped": skipped, "differed": differed,
        "note": "skipped and differed are censuses, not omissions: a file left behind "
                "and a file the run never wrote are indistinguishable afterwards",
    }, indent=2) + "\n")
    print(f"{len(copied)} files into {out}")
    for row in skipped:
        print(f"  skipped {row['key']} ({row['bytes'] / 1e6:.1f} MB): {row['why']}")
    for row in differed:
        print(f"  DIFFERS from the local copy, left alone: {row['key']}")
    print(f"  manifest {manifest_path}")
    if differed:
        raise SystemExit(f"{len(differed)} files differ from what is already here. A "
                         f"results directory holding two runs' answers to one question "
                         f"cannot be read; move or delete the local copies first.")


if __name__ == "__main__":
    main()
