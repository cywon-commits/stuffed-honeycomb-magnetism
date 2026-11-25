#!/bin/bash
# Usage: run_n36.sh SECTOR_SPEC R TAG WORKDIR
# Restartable: rerunning the same command resumes from WORKDIR/TAG checkpoints.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SPEC=$1; R=$2; TAG=$3; WORK=$4
mkdir -p "$WORK" "$HERE/results/N36"
LOG="$HERE/results/N36/${TAG}.log"
echo "=== start $(date -u +%FT%TZ) host_threads=$(nproc)" >> "$LOG"
OMP_NUM_THREADS=$(nproc) stdbuf -oL "$HERE/ed_sym" "$SPEC" --r "$R" --tol 1e-6 --maxit 800 \
  --ckpt "$WORK/$TAG" --ckevery 10 --out "$WORK/$TAG.vec" \
  --pairs "$(dirname "$SPEC")/pairs.txt" >> "$LOG" 2>&1
echo "=== exit $? $(date -u +%FT%TZ)" >> "$LOG"
