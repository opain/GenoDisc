#!/bin/bash
# Wrapper around `sbatch` for Snakemake's generic-cluster `cluster:` submit
# command. Mirrors slurm-status.sh's philosophy for the SUBMIT path.
#
# The SLURM controller (slurmctld) can briefly time out / refuse RPCs when it is
# overloaded (many concurrent orchestrators under one shared account). sbatch
# then exits non-zero with e.g.
#   "Batch job submission failed: Socket timed out on send/recv operation"
#   "Unable to contact slurm controller (connect failure)"
# Snakemake treats ANY non-zero submit as a failed job, and restart-times is far
# too small to ride out a sustained overload window - so a transient controller
# blip would abort the whole multi-hour run even though nothing is actually
# wrong (zero rules failed). A busy controller should mean "wait", not "fail".
#
# So: retry TRANSIENT controller/network failures with backoff; surface a
# GENUINE rejection (bad partition, quota, malformed jobscript) immediately so a
# real error isn't masked by endless retrying.
#
# On success `sbatch --parsable` prints "<jobid>[;<cluster>]" to stdout; this
# script passes that through UNCHANGED - Snakemake parses stdout for the job id,
# so stdout must carry nothing else. All diagnostics go to stderr.
#
# Invoked via `bash <script>` from the cluster: command because the CIFS mount
# doesn't reliably expose the execute bit (same reason as slurm-status.sh and
# run_genodisc.sh). No `set -e`: sbatch failures are handled explicitly below.
set -uo pipefail

# Transient conditions worth waiting out (matched case-insensitively against
# sbatch's stderr). Anything NOT matching is treated as a real rejection.
transient_re='Socket timed out|Unable to contact slurm controller|Connection timed out|Connection refused|Resource temporarily unavailable|slurm_persist_conn|Transport endpoint is not connected|try again|timed out on send/recv'

max_attempts="${GENODISC_SBATCH_RETRIES:-15}"
attempt=1
while :; do
  err_file=$(mktemp)
  out=$(sbatch "$@" 2>"$err_file")
  rc=$?
  err=$(cat "$err_file"); rm -f "$err_file"

  if [ "$rc" -eq 0 ] && [ -n "$out" ]; then
    printf '%s\n' "$out"        # the external job id, for Snakemake
    exit 0
  fi

  if printf '%s' "$err" | grep -qiE "$transient_re"; then
    if [ "$attempt" -ge "$max_attempts" ]; then
      echo "sbatch-submit: controller still unreachable after ${attempt} attempts; giving up: ${err}" >&2
      exit 1
    fi
    sleep_s=$(( attempt * 5 )); [ "$sleep_s" -gt 30 ] && sleep_s=30
    echo "sbatch-submit: transient controller error (attempt ${attempt}/${max_attempts}), retrying in ${sleep_s}s: ${err}" >&2
    attempt=$(( attempt + 1 ))
    sleep "$sleep_s"
    continue
  fi

  # Non-transient failure (genuine rejection) or a zero exit with empty output:
  # surface it rather than loop forever on a real error. Never exit 0 here, or
  # Snakemake would read an empty job id.
  echo "sbatch-submit: submission failed (rc=${rc}), not a recognised transient error: ${err:-<no stderr>}" >&2
  [ "$rc" -eq 0 ] && rc=1
  exit "$rc"
done
