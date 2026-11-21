#!/bin/bash
# The UserData, which is now only a bootstrap: it arms the dead man's switch, fetches the
# real run script from S3 and hands over.
#
# WHY THE SPLIT EXISTS. EC2 caps UserData at 16,384 bytes. The run script passed 15,100 with
# four stages and the fifth took it over, so the choice was to delete the reasoning comments
# from the one artefact that is reviewed before a launch nobody can log into, or to stop
# shipping the script through UserData at all. The instance already copies its code and its
# query sets from our own prefix with the instance profile's credentials, so the script is
# one more object on a path that is already load-bearing - and the ceiling is gone for good
# rather than pushed a few hundred bytes away.
#
# WHAT IS KEPT HERE AND WHY. The dead man's switch is armed BEFORE the fetch, because a
# bootstrap that cannot reach S3 must still terminate the machine; that is the one failure
# this file has to survive on its own. Everything after the handover is the run script's
# problem, and it re-arms the timer itself.
#
# THE ONE THING WE LOSE: if the fetch fails there is no script to upload a log with, so the
# failure leaves nothing in the store and has to be read from the console. That is why the
# echo below names the exact key - a missing object is then one listing away from obvious.
set -xeuo pipefail
exec > >(tee /dev/console /var/log/hapax-bootstrap.log) 2>&1
echo "=== hapax bootstrap $(date -Is) commit @@COMMIT@@ ==="

shutdown -h +@@HARD_STOP_MINUTES@@ &

REGION="@@REGION@@"
KEY="s3://@@BUCKET@@/@@RUN_SCRIPT_KEY@@"
echo "=== fetching run script from $KEY ==="

mkdir -p /opt/hapax
if ! aws s3 cp "$KEY" /opt/hapax/run.sh --region "$REGION"; then
  echo "=== the run script did not copy from $KEY; shutting down ==="
  shutdown -c || true
  shutdown -h now || systemctl --no-block poweroff || halt -p || true
  exit 1
fi
chmod +x /opt/hapax/run.sh
# `exec` so that the run script is this process: the timer stays armed against it and
# nothing is left waiting on a parent that has gone.
exec bash /opt/hapax/run.sh
