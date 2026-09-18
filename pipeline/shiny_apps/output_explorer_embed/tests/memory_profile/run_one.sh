#!/bin/bash
# Run one profile: start app, sampler, driver; collect into $RESULTS_DIR/N<N>/.
# Usage: bash run_one.sh <N> <bundle_tarball>
#
# Env vars (all optional):
#   RSCRIPT     - Rscript binary (default: Rscript on PATH)
#   PYTHON      - python binary  (default: python on PATH)
#   PORT        - app port       (default: 4321)
#   APP_DIR     - shiny app dir  (default: two levels up from this script)
#   RESULTS_DIR - output dir     (default: <script dir>/results)
#   FIREFOX_BIN, GECKODRIVER_BIN - passed through to profile_run.py
#   PROFILE_LD_LIBRARY_PATH - prepended to LD_LIBRARY_PATH before running python
#                             (needed on conda-based Firefox: point at env/lib)
set -euo pipefail

N="${1:?N required}"
BUNDLE="${2:?bundle path required}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON="${PYTHON:-python}"
PORT="${PORT:-4321}"
APP_DIR="${APP_DIR:-$( cd "$SCRIPT_DIR/../.." && pwd )}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"

OUTDIR="$RESULTS_DIR/N$(printf '%02d' "$N")"
mkdir -p "$OUTDIR"

APP_LOG="$OUTDIR/shiny.log"
SAMPLER_CSV="$OUTDIR/samples.csv"

# Per-run TMPDIR so we can scope /tmp footprint measurement to just this
# shiny R process. R uses TMPDIR as the parent of its RtmpXXXX session dir.
RUN_TMP="/tmp/shiny_profile_run_$$_N${N}"
mkdir -p "$RUN_TMP"

# Kill any prior app on the port
pkill -f "shiny::runApp.*port = ${PORT}" 2>/dev/null || true
sleep 1

# Start app in background with scoped TMPDIR
env TMPDIR="$RUN_TMP" "$RSCRIPT" \
  -e "shiny::runApp('$APP_DIR', port = ${PORT}, host = '127.0.0.1', launch.browser = FALSE)" \
  > "$APP_LOG" 2>&1 &
APP_PID=$!
echo "started app pid=$APP_PID  log=$APP_LOG"

for i in $(seq 1 60); do
  if curl -s "http://127.0.0.1:${PORT}" > /dev/null 2>&1; then break; fi
  if ! kill -0 "$APP_PID" 2>/dev/null; then
    echo "app died during startup; log tail:"; tail -30 "$APP_LOG"; exit 1
  fi
  sleep 1
done
echo "app up"

bash "$SCRIPT_DIR/profile_sampler.sh" "$APP_PID" "$SAMPLER_CSV" 0.5 "$RUN_TMP" &
SAMPLER_PID=$!
echo "started sampler pid=$SAMPLER_PID"

sleep 3  # baseline samples

# Prepend PROFILE_LD_LIBRARY_PATH if set (used on conda-managed Firefox builds)
if [ -n "${PROFILE_LD_LIBRARY_PATH:-}" ]; then
  export LD_LIBRARY_PATH="$PROFILE_LD_LIBRARY_PATH${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

"$PYTHON" "$SCRIPT_DIR/profile_run.py" \
    --url "http://127.0.0.1:${PORT}" \
    --bundle "$BUNDLE" \
    --outdir "$OUTDIR" 2>&1 | tee "$OUTDIR/driver.log"

sleep 3  # let sampler catch a couple final ticks

kill "$SAMPLER_PID" 2>/dev/null || true
wait "$SAMPLER_PID" 2>/dev/null || true
kill "$APP_PID" 2>/dev/null || true
wait "$APP_PID" 2>/dev/null || true
rm -rf "$RUN_TMP" 2>/dev/null || true

echo "run N=$N done  ->  $OUTDIR"
