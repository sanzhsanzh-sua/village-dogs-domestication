#!/bin/bash
# Auto-restart wrapper for 01_download_and_align.sh. The pipeline has
# been repeatedly killed mid-run by this environment under disk
# pressure (confirmed 3x: correlates with fqd_tmp/aligned BAM
# accumulation eating free disk, not a bug in the pipeline's own logic
# -- see aligned/align.master.log history around 2026-08-17/18). When
# that happens the killed process can't log its own failure or resume
# itself; something has to notice and relaunch it.
#
# 01_download_and_align.sh exits 0 only when it completes a full pass
# over the sample sheet without being killed (every sample either gets
# a verified BAM or permanently gives up after its own 2 in-script
# retry attempts). Any other exit (killed by signal, disk-full crash,
# etc.) is treated as abnormal and triggers a restart, after clearing
# fqd_tmp (always-safe scratch) to maximize the next attempt's disk
# headroom. Capped at MAX_RESTARTS to avoid an infinite loop if
# something is fundamentally broken rather than transiently
# disk-starved.
set -uo pipefail
cd "$(dirname "$0")/.."

SAMPLE_SHEET="${1:-metadata/pilot6_samples.tsv}"
LOG=aligned/supervisor.log
MAX_RESTARTS=50

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

mkdir -p aligned

RESTART_COUNT=0
while true; do
  RESTART_COUNT=$((RESTART_COUNT + 1))
  if [ "$RESTART_COUNT" -gt "$MAX_RESTARTS" ]; then
    log "*** supervisor: reached $MAX_RESTARTS restarts without a clean finish -- stopping, needs manual attention ***"
    exit 1
  fi

  AVAIL_GB=$(df -g . | tail -1 | awk '{print $4}')
  log "=== supervisor: launch attempt $RESTART_COUNT (disk avail: ${AVAIL_GB}GB) ==="
  rm -rf fqd_tmp/* 2>/dev/null

  bash scripts/01_download_and_align.sh "$SAMPLE_SHEET"
  EXIT_CODE=$?

  if [ "$EXIT_CODE" -eq 0 ]; then
    log "=== supervisor: 01_download_and_align.sh completed a full clean pass (exit 0) -- done ==="
    break
  fi

  log "=== supervisor: 01_download_and_align.sh exited abnormally (code $EXIT_CODE) -- restarting in 15s ==="
  sleep 15
done

log "=== supervisor: finished ==="
