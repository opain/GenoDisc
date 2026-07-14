#!/bin/bash
#SBATCH --partition=neurohack_cpu
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --time=04:00:00
# ============================================================================
# run_genodisc.sh
#
# THIS IS A TRACKED COPY. The web app actually sbatches the live copy at
# GenoDisc/bin/run_genodisc.sh (sibling of activate_genodisc.sh), not this
# file directly. They're separate copies, NOT a symlink (symlinks across this
# CIFS mount don't resolve correctly when accessed from the HPC side - sbatch
# reads the raw mfsymlinks reparse-point data instead of the script and
# rejects it). Edit one, then manually copy the change to the other.
# Whichever you edit first, do the copy in the same sitting - don't let them
# drift.
#
# SLURM job script that runs the GenoDisc pipeline for one submission.
# Invoked by the web app via `sbatch run_genodisc.sh <job_dir>`.
#
# Contract with the web app:
#
#   INPUT (web app must have written these before invoking):
#     - <job_dir>/input/<filename>.gz   uploaded sumstats file
#     - <job_dir>/gwas_list.txt         space-delimited, one row of metadata
#     - <job_dir>/config.yaml           override config with the 5 MVP keys
#
#   OUTPUT (this script produces):
#     - <job_dir>/results/...                       raw Snakemake outputs
#     - <job_dir>/results/results_package.rds       the .rds the web app reads
#     - SLURM stdout/stderr                         logs for human debugging
#
# Success signal: presence of <job_dir>/results/results_package.rds AND
# exit code 0. The web app determines this by polling sacct and stat-ing
# the .rds file.
#
# Exit codes:
#   0   Pipeline completed and results_package.rds exists.
#   1   Pipeline failed at some stage. Check SLURM stdout/stderr.
#   2   Setup error (job dir missing, activation failed, etc.).
#
# Usage:
#   run_genodisc.sh <absolute_path_to_job_dir>
#
# Typical SLURM invocation by the web app (paramiko ssh into HPC):
#   sbatch \
#     --job-name=genodisc-<uuid> \
#     --chdir=${GD}/jobs/<uuid> \
#     --output=${GD}/jobs/<uuid>/slurm.out \
#     ${GD}/bin/run_genodisc.sh ${GD}/jobs/<uuid>
#
# The web app can override partition, time, cpus-per-task, etc. by passing
# the corresponding sbatch flags — they win over the #SBATCH directives above.
# ============================================================================

set -uo pipefail

# ---- Config -----------------------------------------------------------------
GD="/scratch/prj/neurohackpain/GenoDisc"
SNAKEFILE="${GD}/repo/current/pipeline/Snakefile"
ACTIVATE_SCRIPT="${GD}/bin/activate_genodisc.sh"
CONDA_PREFIX_DIR="${GD}/conda-envs"

# Cores: SLURM tells us via SLURM_CPUS_PER_TASK; fall back to 5 for manual runs.
CORES="${SLURM_CPUS_PER_TASK:-5}"

# Memory (MB): SLURM tells us via SLURM_MEM_PER_NODE when --mem is set; fall back to
# 100000 (100G, matching the #SBATCH --mem default above) for manual runs. Passed to
# snakemake as --resources mem_mb so it actually throttles parallel rule scheduling to
# fit the allocation - rules declare their own mem_mb (e.g. run_twas wants 20000 per
# panel), but without this budget Snakemake schedules purely by CPU-slot availability
# and ignores those declarations entirely, which is what caused OOM kills once several
# 20G-per-panel rules landed in the same parallel wave (see misc/05-status.md).
#
# Two non-obvious gotchas hit while adding this (see misc/05-status.md):
# - --resources must come before -c/the target positional arg, not after - it takes
#   nargs='+' and will otherwise greedily swallow the next bare argument (the target
#   path) and crash trying to parse it as a key=value pair.
# - --default-resources mem_mb=... disk_mb=... is required alongside --resources, or
#   rules that don't declare their own mem_mb/disk_mb get an internal "TBD" placeholder
#   that can't be compared against the --resources budget and crashes the scheduler.
MEM_MB="${SLURM_MEM_PER_NODE:-100000}"

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
TARGET="${JOB_DIR}/results/results_package.rds"

# ---- Pre-flight checks ------------------------------------------------------
log "GenoDisc pipeline run starting"
log "  job_uuid:    ${JOB_UUID}"
log "  job_dir:     ${JOB_DIR}"
log "  cores:       ${CORES}"
log "  mem_mb:      ${MEM_MB}"
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
PIPELINE_DIR="$(dirname "$SNAKEFILE")"
DEFAULT_CONFIG="${PIPELINE_DIR}/config.yaml"

log "Invoking Snakemake"
log "  target:          $TARGET"
log "  --directory:     $JOB_DIR"
log "  default config:  $DEFAULT_CONFIG"
log "  override config: $CONFIG_FILE"
log "  --resources:     mem_mb=$MEM_MB"

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
    --use-conda \
    --conda-frontend mamba \
    --conda-prefix   "$CONDA_PREFIX_DIR" \
    --default-resources mem_mb=4000 disk_mb=4000 \
    --resources      mem_mb="$MEM_MB" \
    --rerun-incomplete \
    -c "$CORES" \
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
