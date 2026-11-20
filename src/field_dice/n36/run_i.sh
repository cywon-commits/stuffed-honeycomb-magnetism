#!/bin/bash
# Task I (36 spins): unit tests, I1 (r=0.8), U14, I2 (r=0.6), then optional I3 (r=0.9) and S^z=4 at r=0.8.
# Resumable: rerun the same command (solver results are cached in $N36_WORK).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export N36_WORK=${N36_WORK:-/tmp/n36_work}
mkdir -p "$N36_WORK"
LOG="$N36_WORK/run_i.log"
cd "$HERE"
echo "=== start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) threads $(nproc)" >> "$LOG"
python3 -u task_i.py ${STEPS:-U11 U12 I1 U14 I2 U14 I3 U14 I1b U14} >> "$LOG" 2>&1
echo "=== exit $? $(date -u +%FT%TZ)" >> "$LOG"
