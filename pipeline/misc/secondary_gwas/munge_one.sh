#!/usr/bin/env bash
# Per-code worker: munge one secondary GWAS through the pipeline's own
# sumstat_cleaner.R (GenoUtils) with --merge_alleles, then move the
# merge-alleles output into the staging dir as <code>.sumstats.gz.
#
# Invoked one of two ways:
#   ./munge_one.sh <code>                                  # local / xargs
#   ./munge_one.sh $(awk NR==<i> targets.csv)              # SLURM array
#
# Reads: targets.csv (in the same dir) columns code,source_path,n
# Writes: <STAGING>/<code>.sumstats.gz, plus per-code work products under
#         <STAGING>/work/<code>.*
#
# Exits non-zero on failure so SLURM records FAILED. Success = mergealleles
# file present, non-empty, first data row has 5 tab-separated fields.
#
# `set -u` is deliberately off: activate_pipe_env.sh's conda hooks reference
# variables like ADDR2LINE / host_alias without defaulting them, and would
# fail fatally under nounset. Everything below still runs under -e -o pipefail.
set -eo pipefail

CODE="${1:?usage: munge_one.sh <code>}"

MISC_DIR="/scratch/prj/neurohackpain/GenoDisc/repo/current/pipeline/misc/secondary_gwas"
RESDIR="/scratch/prj/neurohackpain/GenoDisc/pipeline_resources"
STAGING="${MISC_DIR}/staging"
WORK="${STAGING}/work/${CODE}"
TARGETS="${MISC_DIR}/targets.csv"
MERGE_ALLELES="${RESDIR}/data/ldsc/w_hm3.snplist"
REF_CHR="${RESDIR}/data/1kg/1KG.Phase3.MAF_001.chr"

source "${MISC_DIR}/activate_pipe_env.sh"

# Look up source path + fallback N from targets.csv (CSV: code,source_path,n,...)
row=$(awk -F',' -v c="$CODE" 'NR>1 && $1==c {print; exit}' "$TARGETS")
if [[ -z "$row" ]]; then
  echo "ERROR: $CODE not found in $TARGETS" >&2
  exit 2
fi
SRC=$(awk -F',' -v c="$CODE" 'NR>1 && $1==c {print $2; exit}' "$TARGETS")
N_FALLBACK=$(awk -F',' -v c="$CODE" 'NR>1 && $1==c {print $3; exit}' "$TARGETS")

if [[ ! -r "$SRC" ]]; then
  echo "ERROR: source not readable: $SRC" >&2
  exit 3
fi

mkdir -p "$WORK"

# Resolve the pipeline's installed cleaner (same script sumstat_qc.smk uses).
CLEANER=$(Rscript --vanilla -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
if [[ -z "$CLEANER" || ! -r "$CLEANER" ]]; then
  echo "ERROR: could not resolve GenoUtils::sumstat_cleaner.R" >&2
  exit 4
fi

# Only pass --n if we have a fallback (the R script accepts NA and prefers
# the per-SNP N column when present in the input).
N_ARG=()
if [[ -n "$N_FALLBACK" && "$N_FALLBACK" != "NA" ]]; then
  N_ARG=(--n "$N_FALLBACK")
fi

Rscript --vanilla "$CLEANER" \
  --sumstats "$SRC" \
  --ref_chr "$REF_CHR" \
  --population EUR \
  --munged T \
  --merge_alleles "$MERGE_ALLELES" \
  "${N_ARG[@]}" \
  --output "$WORK/${CODE}.cleaned"

OUT="$WORK/${CODE}.cleaned.munged.mergealleles.sumstats.gz"
if [[ ! -s "$OUT" ]]; then
  echo "ERROR: expected output missing/empty: $OUT" >&2
  exit 5
fi

# Sanity: header must be the pipeline's mergealleles schema.
# Use `read -r … < <(…)` — a plain `zcat … | head -n 1` under
# `set -o pipefail` returns 141 (SIGPIPE) because `head` closes the pipe
# after one line while `zcat` still has more to emit. Same header the
# primary GWAS side of ldsc.py --rg uses.
IFS= read -r hdr < <(zcat "$OUT") || true
EXPECTED_HDR=$'CHR\tSNP\tBP\tA1\tA2\tZ\tN'
if [[ "$hdr" != "$EXPECTED_HDR" ]]; then
  echo "ERROR: $OUT unexpected header: '$hdr'" >&2
  exit 6
fi

install -D -m 0644 "$OUT" "$STAGING/${CODE}.sumstats.gz"
echo "OK: ${CODE} -> $STAGING/${CODE}.sumstats.gz"
