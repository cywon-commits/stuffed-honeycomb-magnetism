#!/bin/bash
# Detached, resumable: unit tests then production tasks. Rerun to resume (solver results are cached).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export N27_WORK=${N27_WORK:-/tmp/n27_work}
mkdir -p "$N27_WORK"
LOG="$N27_WORK/run_all.log"
cd "$HERE"
echo "=== start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) threads $(nproc)" >> "$LOG"
python3 -u tests.py ${TESTS:-U3 U4 U5} >> "$LOG" 2>&1
python3 -u tasks.py ${TASKS:-A B D C} >> "$LOG" 2>&1
echo "=== exit $? $(date -u +%FT%TZ)" >> "$LOG"
