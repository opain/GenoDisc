#!/usr/bin/env bash
# SLURM array wrapper. Submit from a login node (compute nodes don't have
# sbatch on KCL ERC-HPC):
#
#   cd /scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas
#   N=$(($(wc -l < targets.csv) - 1))
#   sbatch --array=1-${N}%10 submit_munge.sh
#
# One array task per row in targets.csv (after the header). Concurrency
# capped at 10 by default via %10 — bump/lower to taste.
#
# Partition follows the same convention as the pipeline's snakemake slurm
# profile (profiles/slurm/config.yaml): $GENODISC_PARTITIONS if set, else
# neurohack_cpu.
#SBATCH --job-name=gd-munge-secondary
#SBATCH --partition=neurohack_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/munge-%A_%a.out
#SBATCH --error=logs/munge-%A_%a.err
set -euo pipefail

MISC_DIR="/scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas"
# Staging lives under the panel dir (see munge_one.sh) to keep the 10+GB
# of intermediate .cleaned.* outputs out of the repo.
STAGING="/scratch/prj/neurohackpain/GenoDisc/pipeline_resources/data/secondary_gwas/staging"
cd "$MISC_DIR"
mkdir -p logs "$STAGING/work"

IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID unset — run under sbatch --array=...}"
# Skip the header (NR=1); array index i maps to file row i+1.
CODE=$(awk -F',' -v i="$IDX" 'NR==i+1 {print $1; exit}' targets.csv)
if [[ -z "$CODE" ]]; then
  echo "ERROR: no code at row $IDX of targets.csv" >&2
  exit 2
fi

exec ./munge_one.sh "$CODE"
