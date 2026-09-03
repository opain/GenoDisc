#!/bin/bash
# Called by Snakemake as `slurm-status.sh <jobid>`. Must print exactly one of
# running/success/failed to stdout - nothing else.
set -euo pipefail

jobid="$1"
state=$(sacct -j "$jobid" --format=State --noheader --parsable2 | head -n1 | tr -d ' ')

case "$state" in
  COMPLETED)
    echo success
    ;;
  RUNNING|PENDING|COMPLETING|CONFIGURING)
    echo running
    ;;
  FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|DEADLINE)
    echo failed
    ;;
  *)
    # Unknown/transient sacct state (e.g. not yet visible in accounting) -
    # let Snakemake keep polling rather than false-failing the job.
    echo running
    ;;
esac
