#!/usr/bin/env bash
set -u
BUILD=/mnt/c/Users/htelr/the-sqisign/build_wsl
LOG=/mnt/c/Users/htelr/the-sqisign/src/quaternion/ref/generic/lll/bench_logs_2026-04-19
cd "$BUILD"
for L in 3 5; do
  for M in alg2 alg3; do
    F="$LOG/L${L}_${M}_10k.log"
    echo "=== L${L} ${M} 10k ===" > "$F"
    /usr/bin/time -v ./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=$L --iterations=10000 --mode=$M >> "$F" 2>&1
    echo "done L${L} ${M}"
  done
done
echo FINISHED
