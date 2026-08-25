#!/bin/bash
# External watchdog, meant to run from cron (independent of any Claude
# Code session process tree). The in-process supervisor
# (02_supervisor.sh) auto-restarts 01_download_and_align.sh when THAT
# gets killed, but this environment has at least once (2026-08-19/20)
# killed the supervisor's own bash loop too -- something under disk
# pressure appears to terminate the whole process group, not just the
# heaviest child. A supervisor can't restart itself once it's dead; only
# something outside its process tree can notice and relaunch it. Cron
# fits that role exactly.
#
# Idempotent and safe to run every few minutes: no-ops if the pipeline
# (supervisor or the inner download/align script) is already running,
# and no-ops once every sample in the sheet has a BAM.
cd /Users/sanzaruathanov/Desktop/village_dogs_data || exit 1

LOG=aligned/watchdog.log
SAMPLE_SHEET=metadata/pilot6_samples.tsv
mkdir -p aligned

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "$LOG"; }

if /usr/bin/pgrep -f "scripts/02_supervisor.sh" > /dev/null 2>&1; then
  exit 0
fi
if /usr/bin/pgrep -f "scripts/01_download_and_align.sh" > /dev/null 2>&1; then
  exit 0
fi

ALL_DONE=1
while IFS=$'\t' read -r RUN _REST; do
  [ "$RUN" = "Run" ] && continue
  [ -z "$RUN" ] && continue
  if [ ! -f "aligned/${RUN}.sorted.bam" ]; then
    ALL_DONE=0
    break
  fi
done < "$SAMPLE_SHEET"

if [ "$ALL_DONE" -eq 1 ]; then
  exit 0
fi

log "=== watchdog: no pipeline process alive and sample sheet incomplete -- relaunching supervisor ==="
rm -rf fqd_tmp/* 2>/dev/null

nohup /bin/bash scripts/02_supervisor.sh "$SAMPLE_SHEET" >> aligned/supervisor_stdout.log 2>&1 &
disown
NEWPID=$!
log "=== watchdog: relaunched supervisor as pid $NEWPID ==="
