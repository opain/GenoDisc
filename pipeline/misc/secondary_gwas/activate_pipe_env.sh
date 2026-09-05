#!/usr/bin/env bash
# Sources the snakemake-built GenoDiscPipe conda env so standalone scripts
# (build_targets.R, munge_one.sh) can call the same Rscript / GenoUtils the
# pipeline uses.
#
# Discovers the env dynamically because snakemake bumps the hashed path
# whenever pipeline/envs/main.yaml changes. Picks the most-recently-modified
# env under conda-envs/ that has both GenoUtils and data.table installed.
#
# Usage: source /scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas/activate_pipe_env.sh
#
# Do NOT `set -u` here — the conda `activate.d` hooks reference variables
# (ADDR2LINE, host_alias, etc.) without defaulting them, and would fail
# fatally under nounset.
set -eo pipefail

# Bring in the base conda toolchain (mostly for PATH sanitisation).
source /scratch/prj/neurohackpain/GenoDisc/bin/activate_genodisc.sh

ENVS_DIR="/scratch/prj/neurohackpain/GenoDisc/conda-envs"
PIPE_ENV=""
for f in $(ls -1td "$ENVS_DIR"/*/); do
  # Strip trailing slash for `conda activate --prefix`
  f_noslash="${f%/}"
  if [[ -x "$f_noslash/bin/Rscript" ]] && \
     [[ -d "$f_noslash/lib/R/library/GenoUtils" ]] && \
     [[ -d "$f_noslash/lib/R/library/data.table" ]]; then
    PIPE_ENV="$f_noslash"
    break
  fi
done

if [[ -z "$PIPE_ENV" ]]; then
  echo "ERROR: no conda env under $ENVS_DIR has both GenoUtils and data.table." \
       "Run the pipeline once to have snakemake build it." >&2
  return 1 2>/dev/null || exit 1
fi

conda activate "$PIPE_ENV"
export PATH="$CONDA_PREFIX/bin:$PATH"
export GENODISC_PIPE_ENV="$PIPE_ENV"
