#!/bin/bash
# Poll /proc/<pid>/status VmRSS + VmData and a scoped TMPDIR footprint.
# Usage: bash profile_sampler.sh <PID> <OUT_CSV> [interval_sec] [rtmp]
# Writes header + rows every interval_sec. Stops when PID exits or on SIGTERM.
set -eu
PID="${1:?pid required}"
OUT="${2:?out csv required}"
INTERVAL="${3:-0.5}"
RTMP="${4:-}"

echo "epoch_ms,vmrss_kb,vmdata_kb,vmsize_kb,vmpeak_kb,tmp_bundle_bytes,tmp_upload_bytes,tmp_total_bytes" > "$OUT"

# If no explicit tmpdir was passed, try /proc/PID/environ as a fallback.
if [ -z "$RTMP" ] && [ -r "/proc/$PID/environ" ]; then
  RTMP=$(tr '\0' '\n' < "/proc/$PID/environ" 2>/dev/null | awk -F= '$1=="R_SESSION_TMPDIR"{print $2; exit}')
fi
echo "sampler: PID=$PID RTMP=${RTMP:-<not-found>}" >&2

cleanup() { exit 0; }
trap cleanup TERM INT

while kill -0 "$PID" 2>/dev/null; do
  status="/proc/$PID/status"
  if [ ! -r "$status" ]; then break; fi
  vmrss=$(awk '/^VmRSS:/  {print $2}' "$status" 2>/dev/null || echo 0)
  vmdata=$(awk '/^VmData:/ {print $2}' "$status" 2>/dev/null || echo 0)
  vmsize=$(awk '/^VmSize:/ {print $2}' "$status" 2>/dev/null || echo 0)
  vmpeak=$(awk '/^VmPeak:/ {print $2}' "$status" 2>/dev/null || echo 0)
  bundle_b=0; upload_b=0; total_b=0
  if [ -n "$RTMP" ] && [ -d "$RTMP" ]; then
    bundle_b=$(du -sb "$RTMP"/Rtmp*/gd_bundle_* "$RTMP"/gd_bundle_* 2>/dev/null | awk '{s+=$1} END{print s+0}')
    upload_b=$(du -sb "$RTMP"/Rtmp*/file* "$RTMP"/file* 2>/dev/null | awk '{s+=$1} END{print s+0}')
    total_b=$(du -sb "$RTMP" 2>/dev/null | awk '{print $1+0}')
  fi
  ts=$(date +%s%3N)
  echo "$ts,${vmrss:-0},${vmdata:-0},${vmsize:-0},${vmpeak:-0},$bundle_b,$upload_b,$total_b" >> "$OUT"
  sleep "$INTERVAL"
done
