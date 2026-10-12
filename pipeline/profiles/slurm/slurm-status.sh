#!/bin/bash
# Called by Snakemake as `slurm-status.sh <jobid>`. Must print exactly one of
# running/success/failed to stdout - nothing else.
#
# sacct queries the SLURM accounting DB (slurmdbd). That service can briefly
# refuse connections ("Resource temporarily unavailable"), and Snakemake treats
# ANY non-zero exit from this hook as "Failed to obtain job status", which aborts
# the ENTIRE workflow (not just the one job). So a single transient DB blip must
# not be fatal: retry a few times, and if sacct is still unreachable report
# "running" so Snakemake keeps polling instead of tearing the run down. Genuinely
# failed jobs are still caught once sacct responds; restart-times re-runs them.
#
# NOTE: intentionally NO `set -e` - a failing sacct must fall through to the
# retry/"running" logic below, not abort the script.
set -uo pipefail

jobid="$1"

state=""
for attempt in 1 2 3; do
  # On a sacct connection failure the substitution is empty / the pipeline exits
  # non-zero; either way `state` stays empty and we retry.
  state=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null | head -n1 | tr -d ' ') || true
  [ -n "$state" ] && break
  [ "$attempt" -lt 3 ] && sleep 5
done

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
    # Empty (sacct unreachable after retries, or job not yet visible in
    # accounting) or any unexpected state - keep polling rather than
    # false-failing the whole workflow.
    echo running
    ;;
esac
