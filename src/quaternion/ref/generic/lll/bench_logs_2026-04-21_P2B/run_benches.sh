#!/usr/bin/env bash
# Phase 2 candidate B: mpz_realloc2 prealloc scaffold — bench matrix.
#
# For each (level, mode) we run twice: baseline (prealloc_mode=0) and prealloc
# (prealloc_mode=1). Output log file names encode the variant so the diff can
# be read directly.
#
# Alg 2 runs are capped at 1k trials because most wall time is the MLLL Cohen
# reference path (not prealloc-affected). For Alg 3 we run 10k trials since
# GRAM dominates.
set -u
BUILD=/mnt/c/Users/htelr/the-sqisign/build_wsl
LOG=/mnt/c/Users/htelr/the-sqisign/src/quaternion/ref/generic/lll/bench_logs_2026-04-21_P2B
cd "$BUILD"

run_one() {
  local L=$1 M=$2 IT=$3 P=$4
  local tag="L${L}_${M}_${IT}_${P}"
  local F="$LOG/${tag}.log"
  echo "=== ${tag} ===" > "$F"
  local args="--level=$L --iterations=$IT --mode=$M"
  [ "$P" = "prealloc" ] && args="$args --prealloc"
  /usr/bin/time -v ./src/quaternion/ref/generic/test/sqisign_bm_mlll $args >> "$F" 2>&1
  echo "done ${tag}"
}

for L in 1 3 5; do
  # Alg 2 at 1k (GRAM contribution small vs MLLL Cohen reference).
  run_one "$L" alg2 1000 baseline
  run_one "$L" alg2 1000 prealloc
  # Alg 3 at 10k (GRAM dominates).
  run_one "$L" alg3 10000 baseline
  run_one "$L" alg3 10000 prealloc
done
echo FINISHED
