#!/bin/bash
# Third request (Task H): H1 j=0.10, 0.125 -> H2 j=0.10, 0.15, 0.20 -> H1 j=0.175. Resumable (cached).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export N27_WORK=${N27_WORK:-/tmp/n27_work}
LOG="$N27_WORK/run_third.log"
cd "$HERE"
echo "=== H start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) threads $(nproc)" >> "$LOG"
python3 -u tasks.py H >> "$LOG" 2>&1
echo "=== H exit $? $(date -u +%FT%TZ)" >> "$LOG"
