#!/bin/bash
# ============================================================================
# run_genodisc.sh  (WORKER, versioned)  -- STAGED FOR REVIEW
#
# On activation this REPLACES repo/current/pipeline/misc/run_genodisc.sh (commit
# it into dev BEFORE cutting a release, so the tag captures it). Each cut release
# then carries a frozen copy at repo/versions/<version>/pipeline/misc/run_genodisc.sh.
#
# This is the real orchestrator: it activates the env and runs Snakemake for one
# submission. Unlike the old script it does NOT hardcode repo/current - it
# resolves the Snakefile RELATIVE TO ITS OWN LOCATION, so whichever version's
# worktree it lives in is the version that runs. The stable dispatcher at
# $GD/bin/run_genodisc.sh (bootstrap) picks which worker to exec.
#
# No #SBATCH directives here: this script is exec'd by the bootstrap (which is
# the submitted batch script and carries the scheduling directives), not
# sbatch'd directly.
#
# Contract with the web app (unchanged):
#   INPUT:  <job_dir>/input/<filename>.gz, <job_dir>/gwas_list.txt, <job_dir>/config.yaml
#   OUTPUT: <job_dir>/results/..., results/package/manifest.json, results/bundle.tar.gz
#   Success: results/bundle.tar.gz exists AND exit 0.
#   Exit codes: 0 ok; 1 pipeline failed; 2 setup error.
#
# Usage: run_genodisc.sh <absolute_path_to_job_dir>
# ============================================================================

set -uo pipefail

# ---- Config -----------------------------------------------------------------
GD="/scratch/prj/neurohackpain/GenoDisc"
ACTIVATE_SCRIPT="${GD}/bin/activate_genodisc.sh"
CONDA_PREFIX_DIR="${GD}/conda-envs"

# Self-locate: this file lives at <pipeline>/misc/run_genodisc.sh, so the
# Snakefile is one level up from misc/. This is what makes the worker run its
# OWN version's code instead of a hardcoded repo/current.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
SNAKEFILE="${PIPELINE_DIR}/Snakefile"

# ---- Plain logging (no colour; goes to SLURM stdout) ------------------------
ts()  { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { echo "[$(ts)] [INFO]  $*"; }
err() { echo "[$(ts)] [ERROR] $*" >&2; }

# ---- Arg parsing ------------------------------------------------------------
if [[ $# -ne 1 ]]; then
    err "Usage: $(basename "$0") <absolute_path_to_job_dir>"
    exit 2
fi

JOB_DIR="$1"
JOB_UUID=$(basename "$JOB_DIR")
GWAS_LIST="${JOB_DIR}/gwas_list.txt"
CONFIG_FILE="${JOB_DIR}/config.yaml"
TARGET="${JOB_DIR}/results/bundle.tar.gz"

# ---- Pre-flight checks ------------------------------------------------------
log "GenoDisc pipeline run starting"
log "  job_uuid:    ${JOB_UUID}"
log "  job_dir:     ${JOB_DIR}"
log "  pipeline:    ${PIPELINE_DIR}"
log "  slurm_jobid: ${SLURM_JOB_ID:-<not in slurm>}"
log "  host:        $(hostname)"

[[ -d "$JOB_DIR"          ]] || { err "Job directory does not exist: $JOB_DIR";        exit 2; }
[[ -f "$GWAS_LIST"        ]] || { err "Missing input: $GWAS_LIST";                    exit 1; }
[[ -f "$CONFIG_FILE"      ]] || { err "Missing input: $CONFIG_FILE";                  exit 1; }
[[ -f "$SNAKEFILE"        ]] || { err "Missing pipeline Snakefile: $SNAKEFILE";       exit 1; }
[[ -f "$ACTIVATE_SCRIPT"  ]] || { err "Missing activation script: $ACTIVATE_SCRIPT";  exit 1; }

# Check at least one input file has been staged.
INPUT_COUNT=$(find "${JOB_DIR}/input" -maxdepth 1 \( -type l -o -type f \) 2>/dev/null | wc -l)
if [[ "$INPUT_COUNT" -eq 0 ]]; then
    err "No input file staged under ${JOB_DIR}/input/"
    exit 1
fi

log "Pre-flight checks passed."

# ---- Activate the env -------------------------------------------------------
log "Sourcing $ACTIVATE_SCRIPT"
# shellcheck source=/dev/null
source "$ACTIVATE_SCRIPT" || { err "activate_genodisc.sh failed"; exit 2; }
log "Snakemake version: $(snakemake --version 2>&1)"

# ---- Invoke Snakemake -------------------------------------------------------
DEFAULT_CONFIG="${PIPELINE_DIR}/config.yaml"

log "Invoking Snakemake"
log "  target:          $TARGET"
log "  --snakefile:     $SNAKEFILE"
log "  --directory:     $JOB_DIR"
log "  default config:  $DEFAULT_CONFIG"
log "  override config: $CONFIG_FILE"
log "  --profile:       ${PIPELINE_DIR}/profiles/slurm"

# IMPORTANT: --configfile takes BOTH paths as values of ONE flag (Snakemake
# 7.32 silently overrides on repeated --configfile flags; see 05-status.md).
# IMPORTANT: --directory <job_dir> requires the pipeline patches from
# commits 7a916db, 3201b37, f1a2b09, d22e016, 1883a84, 2b54dbb. Without
# those, the pipeline assumes CWD is the pipeline folder and breaks.

# A killed run (SLURM TIMEOUT, kill signal, power loss) leaves the
# .snakemake lock in place, which makes the next run refuse to start with
# LockException. We unlock automatically so Modify & Rerun "just works" -
# but only after confirming via squeue that no OTHER active SLURM job is
# still targeting this same job_dir. Every submission/rerun for a given job
# shares job_dir but gets its own --job-name (original: JOB_UUID; reruns:
# JOB_UUID_r<timestamp> - see job_service.py), so matching on that prefix
# and excluding ourselves catches genuine concurrent runs. If one is found,
# the lock is left alone and Snakemake will correctly refuse to start.
OTHER_ACTIVE=$(squeue -u "$(whoami)" -h -o "%i %j" | awk -v uuid="$JOB_UUID" -v self="${SLURM_JOB_ID:-}" \
    '$2 == uuid || $2 ~ ("^" uuid "_r") { if ($1 != self) print $1 }')

if [[ -n "$OTHER_ACTIVE" ]]; then
    log "Another active SLURM job (${OTHER_ACTIVE}) is already targeting this job directory - leaving any lock in place."
else
    log "Unlocking any stale Snakemake lock (safe no-op if none present)"
    snakemake \
        --snakefile  "$SNAKEFILE" \
        --directory  "$JOB_DIR" \
        --configfile "$DEFAULT_CONFIG" "$CONFIG_FILE" \
        --unlock
fi

START_TIME=$(date +%s)

snakemake \
    --snakefile      "$SNAKEFILE" \
    --directory      "$JOB_DIR" \
    --configfile     "$DEFAULT_CONFIG" "$CONFIG_FILE" \
    --conda-prefix   "$CONDA_PREFIX_DIR" \
    --profile        "${PIPELINE_DIR}/profiles/slurm" \
    "$TARGET"
SNAKEMAKE_EXIT=$?

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

if [[ $SNAKEMAKE_EXIT -ne 0 ]]; then
    err "Snakemake exited with code ${SNAKEMAKE_EXIT} after ${DURATION}s"
    exit 1
fi

log "Snakemake completed in ${DURATION}s"

# ---- Verify the expected output ---------------------------------------------
if [[ ! -f "$TARGET" ]]; then
    err "Snakemake reported success but ${TARGET} is missing"
    exit 1
fi
log "Output present: $TARGET"
log "Done."
exit 0
