#!/bin/bash
# Second request (E-G).  PART=E runs U8, task E, U9 immediately (2 threads, alongside task C).
# PART=GF waits for a running run_all.sh (task C) to finish, then regenerates the task C summary
# notes (cached), runs task G and task F (+F2 unless EXTRA is set to empty).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export N27_WORK=${N27_WORK:-/tmp/n27_work}
LOG="$N27_WORK/run_second.log"
cd "$HERE"
if [ "${PART:-GF}" = "E" ]; then
  export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
  echo "=== E start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) threads $OMP_NUM_THREADS" >> "$LOG"
  python3 -u tests.py U8 >> "$LOG" 2>&1
  python3 -u tasks.py E >> "$LOG" 2>&1
  python3 -u tests.py U9 >> "$LOG" 2>&1
  echo "=== E exit $? $(date -u +%FT%TZ)" >> "$LOG"
else
  while pgrep -f "run_all.sh" > /dev/null || pgrep -f "PART=E" > /dev/null; do sleep 30; done
  echo "=== GF start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) threads $(nproc)" >> "$LOG"
  python3 -u tasks.py C >> "$LOG" 2>&1
  python3 -u tasks.py G F ${EXTRA-F2} >> "$LOG" 2>&1
  echo "=== GF exit $? $(date -u +%FT%TZ)" >> "$LOG"
fi
