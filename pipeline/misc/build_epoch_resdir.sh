#!/usr/bin/env bash
#
# build_epoch_resdir.sh — construct a per-reference-epoch GenoDisc resource
# directory as a symlink farm over a shared "common" store.
#
# WHY. GenoDisc historically used ONE shared resdir (pipeline_resources) for
# every pipeline version. The 1kg+HGDP reference migration (v1.3.6) rebuilt
# reference-DERIVED resources in place, which (a) left old pinned versions
# reading resources inconsistent with their code and (b) let a re-triggered
# NA-source rule delete a pre-provisioned output (the 2026-11-03 prod
# k2260822 rosmap_smr .besd.epi incident). This tool splits the store into:
#
#     <root>/common/        reference-INDEPENDENT resources, built ONCE, shared
#     <root>/v-<epoch>/     per-epoch resdir: real reference-DERIVED dirs +
#                           symlinks into common/ for everything shareable
#
# Epochs map from pipeline version at the v1.3.6 boundary:
#     legacy  = versions < v1.3.6  (old 1kg reference, data/1kg)
#     dense   = versions >= v1.3.6 (dense 1kg+HGDP reference, data/1kg_hgdp)
#
# SAFETY. Jobs only ever READ resources, so a v-<epoch> symlink farm is safe:
# Snakemake follows symlinks to the target's mtime, and the rosmap-smr rule now
# wraps its reference inputs in ancient() (committed to dev) so nothing
# re-triggers. Provisioning writes only into the real stores (common/ and the
# per-epoch real dirs), never through a symlink. This script DOES NOT run the
# Snakemake DAG and DOES NOT touch the Django resdir wiring — it reorganises
# existing resource data. Default MODE is dry-run: it prints the plan and the
# classification audit and changes nothing.
#
# USAGE
#   build_epoch_resdir.sh --epoch dense|legacy [options]
#
#   --epoch   dense|legacy                 (required)
#   --gd      <GenoDisc base dir>          (default: /scratch/prj/neurohackpain/GenoDisc)
#   --src     <existing resdir to migrate> (default: <gd>/pipeline_resources)
#   --root    <new resources root>         (default: <gd>/resources)
#   --mode    dry-run|link|move|copy       (default: dry-run)
#               dry-run  print the plan + audit, change nothing
#               link     hardlink files from --src into common/ & the epoch dir
#                        (instant, no extra bytes; needs hardlink support on the FS)
#               move     rename from --src (instant; DESTRUCTIVE to --src, so run
#                        ONLY when no jobs are reading --src)
#               copy     copy from --src (safe, leaves --src intact; slow on CIFS)
#   -h|--help
#
# Run dry-run first and read the audit. The fail-on-unclassified guard means any
# resource subdir not listed below aborts the run, so new resources get reviewed
# rather than silently mis-placed.

set -euo pipefail

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
GD="/scratch/prj/neurohackpain/GenoDisc"
SRC=""
ROOT=""
EPOCH=""
MODE="dry-run"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --epoch) EPOCH="$2"; shift 2;;
    --gd)    GD="$2"; shift 2;;
    --src)   SRC="$2"; shift 2;;
    --root)  ROOT="$2"; shift 2;;
    --mode)  MODE="$2"; shift 2;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "ERROR: unknown argument: $1" >&2; exit 2;;
  esac
done

SRC="${SRC:-$GD/pipeline_resources}"
ROOT="${ROOT:-$GD/resources}"
COMMON="$ROOT/common"
VIEW="$ROOT/v-$EPOCH"

case "$EPOCH" in dense|legacy) ;; *) echo "ERROR: --epoch must be 'dense' or 'legacy'" >&2; exit 2;; esac
case "$MODE" in dry-run|link|move|copy) ;; *) echo "ERROR: --mode must be dry-run|link|move|copy" >&2; exit 2;; esac
[[ -d "$SRC" ]] || { echo "ERROR: --src not a directory: $SRC" >&2; exit 2; }

# ---------------------------------------------------------------------------
# Classification (data/ entries). Derived from the resource-rule audit.
#   SHARED = reference-INDEPENDENT (external download / reference-independent
#            formatting) -> symlinked to common/, identical across epochs.
#   REF    = reference-DERIVED -> a real directory inside each v-<epoch>.
# software/ is ALL shared (downloads/installs; nothing reference-derived).
# Any data/ or software/ entry NOT matched here aborts the run (see guard).
# ---------------------------------------------------------------------------

# Reference-DERIVED data/ entries (real per-epoch). Both epochs' reference dirs
# are listed; only the one matching the epoch is expected to exist in --src.
REF_DATA=(
  1kg_hgdp            # dense reference panel (prep_1kg; epoch=dense)
  1kg                 # legacy reference panel (prep_1kg; epoch=legacy)
  magma_ref           # MAGMA LD ref: 1kg_hgdp_* is reference-derived (orphan g1000_eur download also lives here)
  magma               # NCBI37.3.genes.annot is reference-derived (NCBI37.3.gene.loc download also lives here)
  predicted_expression # feature_pred / format_pred / build_twas_gsea_cormat: predicted into the reference sample
  rosmap_smr          # besd .esi MAF is reference-derived (.besd/.epi restored 2026-11-03)
  ref                 # transient staging (prep_1kg deletes it; usually empty)
)

# Reference-DERIVED data/ FLAG FILES (bare files at data/ top level).
REF_DATA_FLAGS=(
  prep_1kg.done          # legacy marker (dense marker lives INSIDE 1kg_hgdp/)
  format_psychencode.done
)

# Reference-epoch SOURCE artifacts (dense-only raw tarballs).
REF_DATA_FILES=(
  genopred_dense_1kg_hgdp.tar.gz
  genopred_dense_incl_ambig_1kg_hgdp.tar.gz
)

# Shared data/ subdirectories (symlink to common/).
SHARED_DATA=(
  GCSC MetaBrain alphamissense atc banner_twas biomart csq drug_targetor
  eqtlgen glist-hg19 gtex ld_scores ldsc liftover msigdb pops psychencode_smr
  recomb_maps rosmap_twas sbayesrc_annot sbayesrc_r_ref sbayesrc_ref
  secondary_gwas sumher_tagging
)

# A shared subdir that contains ONE reference-derived child: fusion_snp_weights
# holds 56 downloaded FUSION weight panels (shared) + 'psychencode', which
# format_psychencode re-matches against the reference .bim (reference-derived).
# It is built as a REAL dir whose children are symlinked to common/ EXCEPT the
# real-per-epoch children listed here.
#   REVIEW: confirm format_psychencode's script actually reads the reference
#   .bim before trusting 'psychencode' as epoch-specific; if it does not, move
#   fusion_snp_weights into SHARED_DATA instead.
declare -A SHARED_WITH_REAL_CHILDREN=(
  [fusion_snp_weights]="psychencode"
)

# Shared "marker file" GLOBS at data/ top level (Snakemake .done flags for
# reference-independent rules). Matched by prefix; symlinked to common/.
SHARED_DATA_FLAG_GLOBS=(
  'download_fusion_gtex_*.done'
  'download_non_gtex_*.done'
  'download_psychENCODE_weights.done'
  'insert_n_nongtex_*.done'
  'update_gtex_coord_*.done'
  'eqtlgen.done'
  'MetaBrain_download.out'
)

# Per-epoch real top-level entries (NOT shared, NOT symlinked).
REAL_TOPLEVEL=( logs benchmarks last_version.txt )

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
say() { printf '%s\n' "$*"; }
run() {
  # Execute (or, in dry-run, just print) a mutating command.
  if [[ "$MODE" == "dry-run" ]]; then
    printf '  [dry-run] %s\n' "$*"
  else
    printf '  [run] %s\n' "$*"
    "$@"
  fi
}

# Relocate a path from --src into a destination store, honouring MODE.
relocate() {
  local src="$1" dst="$2"
  [[ -e "$dst" || -L "$dst" ]] && { say "  [skip] already present: $dst"; return 0; }
  case "$MODE" in
    dry-run) printf '  [dry-run] relocate(%s): %s -> %s\n' "$MODE" "$src" "$dst";;
    link)    run cp -al "$src" "$dst";;
    move)    run mv "$src" "$dst";;
    copy)    run cp -a "$src" "$dst";;
  esac
}

# Create a relative symlink dst -> target (target expressed relative to dst's dir).
symlink_rel() {
  local target_abs="$1" dst="$2"
  local rel; rel="$(realpath -m --relative-to="$(dirname "$dst")" "$target_abs")"
  [[ -L "$dst" || -e "$dst" ]] && { say "  [skip] link exists: $dst"; return 0; }
  run ln -s "$rel" "$dst"
}

is_in() { local n="$1"; shift; local x; for x in "$@"; do [[ "$x" == "$n" ]] && return 0; done; return 1; }

matches_shared_flag() {
  local n="$1" g
  for g in "${SHARED_DATA_FLAG_GLOBS[@]}"; do
    # shellcheck disable=SC2053
    [[ "$n" == $g ]] && return 0
  done
  return 1
}

# ---------------------------------------------------------------------------
# 1. Classification audit + fail-on-unclassified guard
# ---------------------------------------------------------------------------
say "=============================================================="
say "build_epoch_resdir.sh  epoch=$EPOCH  mode=$MODE"
say "  src    = $SRC"
say "  common = $COMMON"
say "  view   = $VIEW"
say "=============================================================="

UNCLASSIFIED=()
declare -A PLAN_DATA   # name -> SHARED|REF|SHARED_MIXED|SHARED_FLAG|REF_FLAG|REF_FILE

for entry in "$SRC"/data/*; do
  [[ -e "$entry" || -L "$entry" ]] || continue
  n="$(basename "$entry")"
  if   [[ -n "${SHARED_WITH_REAL_CHILDREN[$n]+x}" ]]; then PLAN_DATA[$n]=SHARED_MIXED
  elif is_in "$n" "${SHARED_DATA[@]}";      then PLAN_DATA[$n]=SHARED
  elif is_in "$n" "${REF_DATA[@]}";         then PLAN_DATA[$n]=REF
  elif is_in "$n" "${REF_DATA_FLAGS[@]}";   then PLAN_DATA[$n]=REF_FLAG
  elif is_in "$n" "${REF_DATA_FILES[@]}";   then PLAN_DATA[$n]=REF_FILE
  elif matches_shared_flag "$n";            then PLAN_DATA[$n]=SHARED_FLAG
  else UNCLASSIFIED+=("data/$n")
  fi
done

# software/ : everything shared (dirs, binaries, and install_*.done markers).
SOFT_ENTRIES=()
for entry in "$SRC"/software/*; do
  [[ -e "$entry" || -L "$entry" ]] || continue
  SOFT_ENTRIES+=("$(basename "$entry")")
done

if (( ${#UNCLASSIFIED[@]} )); then
  say ""
  say "ERROR: ${#UNCLASSIFIED[@]} unclassified resource entr(ies) — refusing to run."
  say "Add each to the SHARED_* or REF_* lists (or a SHARED_WITH_REAL_CHILDREN entry):"
  printf '  - %s\n' "${UNCLASSIFIED[@]}"
  exit 3
fi

say ""
say "--- Classification audit (data/) ---"
for n in $(printf '%s\n' "${!PLAN_DATA[@]}" | sort); do
  printf '  %-14s %s\n' "${PLAN_DATA[$n]}" "data/$n"
done
say "  SHARED  software/* : ${#SOFT_ENTRIES[@]} entries (all symlinked to common)"
say ""

# ---------------------------------------------------------------------------
# 2. Populate common/ (shared, reference-independent)  [from --src]
# ---------------------------------------------------------------------------
say "--- Step 1: populate common store ($COMMON) ---"
run mkdir -p "$COMMON/data" "$COMMON/software"

for n in "${!PLAN_DATA[@]}"; do
  case "${PLAN_DATA[$n]}" in
    SHARED|SHARED_FLAG) relocate "$SRC/data/$n" "$COMMON/data/$n";;
    SHARED_MIXED)
      # Real parent in common with its shareable children; the reference-derived
      # child(ren) are NOT placed in common (they are per-epoch, see step 3).
      run mkdir -p "$COMMON/data/$n"
      read -r -a real_children <<< "${SHARED_WITH_REAL_CHILDREN[$n]}"
      for child in "$SRC/data/$n"/*; do
        [[ -e "$child" || -L "$child" ]] || continue
        cn="$(basename "$child")"
        is_in "$cn" "${real_children[@]}" && continue
        relocate "$child" "$COMMON/data/$n/$cn"
      done
      ;;
  esac
done
for s in "${SOFT_ENTRIES[@]}"; do relocate "$SRC/software/$s" "$COMMON/software/$s"; done

# ---------------------------------------------------------------------------
# 3. Build the per-epoch view (real reference-derived + symlinks to common)
# ---------------------------------------------------------------------------
say ""
say "--- Step 2: build epoch view ($VIEW) ---"
run mkdir -p "$VIEW/data" "$VIEW/software"

# 3a. Shared data/ + software/ -> relative symlinks into common/.
for n in "${!PLAN_DATA[@]}"; do
  case "${PLAN_DATA[$n]}" in
    SHARED|SHARED_FLAG) symlink_rel "$COMMON/data/$n" "$VIEW/data/$n";;
  esac
done
for s in "${SOFT_ENTRIES[@]}"; do symlink_rel "$COMMON/software/$s" "$VIEW/software/$s"; done

# 3b. Reference-derived data/ -> real per-epoch (relocated from --src).
for n in "${!PLAN_DATA[@]}"; do
  case "${PLAN_DATA[$n]}" in
    REF|REF_FLAG|REF_FILE) relocate "$SRC/data/$n" "$VIEW/data/$n";;
  esac
done

# 3c. Mixed dirs: real parent, shareable children symlinked to common, the
#     reference-derived child(ren) real per-epoch.
for n in "${!SHARED_WITH_REAL_CHILDREN[@]}"; do
  [[ "${PLAN_DATA[$n]:-}" == SHARED_MIXED ]] || continue
  run mkdir -p "$VIEW/data/$n"
  read -r -a real_children <<< "${SHARED_WITH_REAL_CHILDREN[$n]}"
  for child in "$SRC/data/$n"/*; do
    [[ -e "$child" || -L "$child" ]] || continue
    cn="$(basename "$child")"
    if is_in "$cn" "${real_children[@]}"; then
      relocate "$child" "$VIEW/data/$n/$cn"                         # reference-derived
    else
      symlink_rel "$COMMON/data/$n/$cn" "$VIEW/data/$n/$cn"         # shareable
    fi
  done
done

# 3d. Per-epoch real top-level dirs/files (Snakemake writes logs/benchmarks here).
for t in "${REAL_TOPLEVEL[@]}"; do
  if [[ -d "$SRC/$t" ]]; then relocate "$SRC/$t" "$VIEW/$t"
  elif [[ -f "$SRC/$t" ]]; then relocate "$SRC/$t" "$VIEW/$t"
  else run mkdir -p "$VIEW/$t"; fi
done

# ---------------------------------------------------------------------------
# 4. Verify (skipped in dry-run — nothing was created)
# ---------------------------------------------------------------------------
if [[ "$MODE" != "dry-run" ]]; then
  say ""
  say "--- Step 3: verify ---"
  broken=0
  while IFS= read -r -d '' l; do
    if [[ ! -e "$l" ]]; then say "  BROKEN SYMLINK: $l"; broken=$((broken+1)); fi
  done < <(find "$VIEW" -type l -print0)
  if (( broken )); then say "  FAIL: $broken broken symlink(s) in $VIEW"; exit 4; fi
  say "  OK: all symlinks in $VIEW resolve."
  say "  Reference-derived (real) data dirs present:"
  for n in "${!PLAN_DATA[@]}"; do
    case "${PLAN_DATA[$n]}" in REF|REF_FLAG|REF_FILE)
      [[ -e "$VIEW/data/$n" && ! -L "$VIEW/data/$n" ]] && printf '    %s\n' "data/$n";;
    esac
  done
fi

say ""
say "Done (mode=$MODE). Next: set the job's resdir to $VIEW for epoch '$EPOCH'"
say "(Django: resource_epoch_for_version() -> resdir=<root>/v-<epoch>; see the plan)."
[[ "$MODE" == "dry-run" ]] && say "This was a DRY RUN — nothing changed. Re-run with --mode link|move|copy to apply."
