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
SHARED_STORE="copy"   # copy = build a real common/ store; src = symlink shared straight to --src
RUN_BUILD=0
BUILD_CORES=4
# Pipeline dir used for the Step-4 build command. Defaults to this script's own
# pipeline (it lives in <pipeline>/misc/); override per epoch with --pipeline.
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd || echo "$GD/repo/current/pipeline")"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --epoch) EPOCH="$2"; shift 2;;
    --gd)    GD="$2"; shift 2;;
    --src)   SRC="$2"; shift 2;;
    --root)  ROOT="$2"; shift 2;;
    --mode)  MODE="$2"; shift 2;;
    --shared-store) SHARED_STORE="$2"; shift 2;;
    --pipeline) PIPE="$2"; shift 2;;
    --cores) BUILD_CORES="$2"; shift 2;;
    --run-build) RUN_BUILD=1; shift;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "ERROR: unknown argument: $1" >&2; exit 2;;
  esac
done

SRC="${SRC:-$GD/pipeline_resources}"
ROOT="${ROOT:-$GD/resources}"
VIEW="$ROOT/v-$EPOCH"
# Where shared (reference-independent) resources live that the view symlinks to:
#   copy -> a real <root>/common store populated from --src (Step 1).
#   src  -> --src itself (pipeline_resources); shared entries symlink straight
#           there, Step 1 is skipped. Avoids copying the huge small-file trees
#           (GCSC/twas_weights, fusion_snp_weights) that make cp pathologically
#           slow on cephfs; pipeline_resources stays as the shared backing store.
if [[ "$SHARED_STORE" == "src" ]]; then COMMON="$SRC"; else COMMON="$ROOT/common"; fi

case "$EPOCH" in dense|legacy) ;; *) echo "ERROR: --epoch must be 'dense' or 'legacy'" >&2; exit 2;; esac
case "$MODE" in dry-run|link|move|copy) ;; *) echo "ERROR: --mode must be dry-run|link|move|copy" >&2; exit 2;; esac
case "$SHARED_STORE" in copy|src) ;; *) echo "ERROR: --shared-store must be copy|src" >&2; exit 2;; esac
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

# Reference-derived resources to BUILD FRESH rather than relocate: these are
# produced by the resource DAG and can be STALE/WRONG in --src after a reference
# change (e.g. rosmap_smr's besd .esi was left at the 2024 EUR-reference values
# when the v1.3.6 dense migration's rebuild failed - the k2260822 incident). For
# each, the script creates an empty real dir and emits a snakemake command
# (Step 4) that rebuilds it into the view against this epoch's reference + the raw
# sources in common/. Other REF resources (the 36G 1kg_hgdp reference itself,
# magma, predicted_expression, ...) were built correctly by the migration and are
# relocated as-is - rebuilding them would needlessly re-download/recompute. Extend
# this list for a fuller clean-room rebuild (add e.g. magma magma_ref
# predicted_expression; note psychencode lives under the fusion_snp_weights mixed
# dir, handled separately).
REBUILD_REF=(
  rosmap_smr
)

# Shared data/ subdirectories (symlink to common/).
SHARED_DATA=(
  GCSC MetaBrain alphamissense atc banner_twas biomart csq drug_targetor
  eqtlgen glist-hg19 gtex ld_scores ldsc liftover msigdb pops psychencode_smr
  recomb_maps rosmap_twas sbayesrc_annot sbayesrc_r_ref sbayesrc_ref
  secondary_gwas sumher_tagging
  # Raw external source datasets (INPUTS to the reference-derived build rules):
  # the ROSMAP/Banner panels and the raw dense-reference tarballs. Reference-
  # INDEPENDENT - they are the inputs FROM which resources are built, not derived
  # from the reference - so they live once in common/ and the build rules read
  # them via resdir; the reference-derived PRODUCTS (rosmap_smr besd) are per-epoch.
  ROSMAP_pQTL ROSMAP-PWAS Banner-PWAS
  genopred_dense_1kg_hgdp.tar.gz genopred_dense_incl_ambig_1kg_hgdp.tar.gz
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
declare -A PLAN_DATA   # name -> SHARED|REF|SHARED_MIXED|SHARED_FLAG|REF_FLAG

for entry in "$SRC"/data/*; do
  [[ -e "$entry" || -L "$entry" ]] || continue
  n="$(basename "$entry")"
  if   [[ -n "${SHARED_WITH_REAL_CHILDREN[$n]+x}" ]]; then PLAN_DATA[$n]=SHARED_MIXED
  elif is_in "$n" "${SHARED_DATA[@]}";      then PLAN_DATA[$n]=SHARED
  elif is_in "$n" "${REF_DATA[@]}";         then PLAN_DATA[$n]=REF
  elif is_in "$n" "${REF_DATA_FLAGS[@]}";   then PLAN_DATA[$n]=REF_FLAG
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
if [[ "$SHARED_STORE" == "src" ]]; then
  say "--- Step 1: SKIPPED (shared-store=src: shared resources stay in $SRC; the view symlinks straight to them) ---"
else
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
fi

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

# 3b. Reference-derived data/ -> real per-epoch. Entries listed in REBUILD_REF are
#     created EMPTY to be rebuilt fresh by Step 4 (they can be stale in --src); the
#     rest are relocated as-is (already correct; cheaper than recomputing).
# The two epochs use DISTINCT reference-panel dirs (dense=1kg_hgdp, legacy=1kg) -
# include only this epoch's, so the dense view doesn't carry the 24G legacy panel
# (and vice-versa). Other reference-derived dirs (magma, magma_ref,
# predicted_expression) share a path across epochs: dense copies the current ones;
# legacy lists them in REBUILD_REF to rebuild from the old reference.
case "$EPOCH" in dense) _other_panel="1kg" ;; legacy) _other_panel="1kg_hgdp" ;; esac
for n in "${!PLAN_DATA[@]}"; do
  [[ "$n" == "$_other_panel" ]] && continue   # skip the other epoch's reference panel
  case "${PLAN_DATA[$n]}" in
    REF|REF_FLAG)
      if is_in "$n" "${REBUILD_REF[@]}"; then
        say "  [rebuild] data/$n -> empty real dir (built by Step 4, not relocated)"
        run mkdir -p "$VIEW/data/$n"
      else
        relocate "$SRC/data/$n" "$VIEW/data/$n"
      fi
      ;;
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
    case "${PLAN_DATA[$n]}" in REF|REF_FLAG)
      [[ -e "$VIEW/data/$n" && ! -L "$VIEW/data/$n" ]] && printf '    %s\n' "data/$n";;
    esac
  done
fi

# ---------------------------------------------------------------------------
# 5. Build step: (re)build the REBUILD_REF resources into the view from source.
#    Only the empty REBUILD_REF dirs (and anything else missing) are produced -
#    relocated REF resources and common/ symlinks are already present so the DAG
#    skips them. The raw sources (ROSMAP/Banner) resolve via resdir -> common/.
# ---------------------------------------------------------------------------
BUILD_CMD=(snakemake prepare_resources
  --snakefile "$PIPE/Snakefile"
  --configfile "$PIPE/prepare_resources.yaml"
  --config "resdir=$VIEW"
  --use-conda -c"$BUILD_CORES")
say ""
say "--- Step 4: build reference-derived resources  (REBUILD_REF: ${REBUILD_REF[*]}) ---"
say "Run on the cluster (dry-run FIRST to confirm only the intended rebuilds fire):"
say "  source $GD/bin/activate_genodisc.sh"
say "  ${BUILD_CMD[*]} -n      # review the plan"
say "  ${BUILD_CMD[*]}         # build"
if [[ "$RUN_BUILD" == "1" && "$MODE" != "dry-run" ]]; then
  say "  [--run-build] activating env and building now..."
  # shellcheck disable=SC1090
  source "$GD/bin/activate_genodisc.sh"
  "${BUILD_CMD[@]}"
else
  say "  (pass --run-build to invoke it from here; omitted by default so the"
  say "   heavy build runs in a deliberate, supervised window.)"
fi

say ""
say "Done (mode=$MODE)."
say "ARM (go live) only AFTER the Step-4 build succeeds and you have verified the view."
say "The web app switches a job to v-<epoch> once its marker exists:"
say "  touch $VIEW/.ready          # arms epoch '$EPOCH'"
if [[ "$SHARED_STORE" == "copy" ]]; then
  say "  touch $COMMON/.ready        # arms the shared common store (copy mode only)"
fi
say "ROLLBACK at any time: rm $VIEW/.ready  (jobs fall back to pipeline_resources)."
[[ "$MODE" == "dry-run" ]] && say "This was a DRY RUN — nothing changed. Re-run with --mode link|move|copy to apply."
