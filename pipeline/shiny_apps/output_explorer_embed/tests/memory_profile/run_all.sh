#!/bin/bash
# Scale the app's example bundle to N ∈ {1,2,4,8,16} GWAS, run one profile
# per size, aggregate results.
#
# Env vars (all optional — see run_one.sh for its own set):
#   RSCRIPT       - Rscript binary                (default: Rscript on PATH)
#   SOURCE_BUNDLE - single-GWAS .tar.gz to scale  (default: <app>/data/als_bundle.tar.gz)
#   BUNDLES_DIR   - where scaled bundles are cached (default: <script dir>/bundles)
#   N_LIST        - space-separated list of Ns    (default: "1 2 4 8 16")
set -uo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RSCRIPT="${RSCRIPT:-Rscript}"
APP_DIR="${APP_DIR:-$( cd "$SCRIPT_DIR/../.." && pwd )}"
SOURCE_BUNDLE="${SOURCE_BUNDLE:-$APP_DIR/data/als_bundle.tar.gz}"
BUNDLES_DIR="${BUNDLES_DIR:-$SCRIPT_DIR/bundles}"
N_LIST="${N_LIST:-1 2 4 8 16}"

mkdir -p "$BUNDLES_DIR"

for N in $N_LIST; do
  OUT="$BUNDLES_DIR/als_N$(printf '%02d' $N).tar.gz"
  if [ ! -f "$OUT" ]; then
    echo "--- scaling to N=$N ---"
    "$RSCRIPT" "$SCRIPT_DIR/scale_bundle.R" --src "$SOURCE_BUNDLE" --n "$N" --out "$OUT"
  fi
done

for N in $N_LIST; do
  BUNDLE="$BUNDLES_DIR/als_N$(printf '%02d' $N).tar.gz"
  echo "########## N=$N ##########"
  bash "$SCRIPT_DIR/run_one.sh" "$N" "$BUNDLE"
  sleep 3
done

echo "--- aggregating ---"
APP_DIR="$APP_DIR" BUNDLES_DIR="$BUNDLES_DIR" "$RSCRIPT" "$SCRIPT_DIR/summary.R"
echo "all runs complete"
