#!/bin/bash
# Phase 2 primary-decision sweep: ibz / prealloc / fp across L1/L3/L5 × Alg2/Alg3.
# Driven from WSL Ubuntu because the build cache is Linux-toolchain.
#
# 2026-04-22 P2.1-time-rerun: iter 분할 (alg2=500, alg3=5000). 이전 30 iter
# 에서는 alg3 timing variance 크고 fp/ibz 비율 해석 불안정. dispatcher
# threshold 정정(§6.8) 직후라 L3/L5 fp 는 이번 sweep 이 최초 실측.
set -e

BUILD=/mnt/c/Users/htelr/the-sqisign/build
BENCH=$BUILD/src/quaternion/ref/generic/test/sqisign_bm_mlll
LOGDIR=/mnt/c/Users/htelr/the-sqisign/src/quaternion/ref/generic/lll/bench_logs/p2_phase2_sweep

ITER_ALG2=${ITER_ALG2:-500}
ITER_ALG3=${ITER_ALG3:-5000}

if [ ! -x "$BENCH" ]; then
    echo "ERROR: benchmark binary not found at $BENCH"
    echo "Rebuild with: cmake --build $BUILD --target sqisign_bm_mlll"
    exit 1
fi

printf '%-4s %-5s %-10s %-6s %-9s %-9s\n' "Lvl" "Mode" "Variant" "Iter" "GRAMms" "MLLLms"
printf '%-4s %-5s %-10s %-6s %-9s %-9s\n' "---" "----" "-------" "----" "------" "------"

for L in 1 3 5; do
    for M in alg2 alg3; do
        if [ "$M" = "alg2" ]; then
            ITER=$ITER_ALG2
        else
            ITER=$ITER_ALG3
        fi
        for V in baseline prealloc fp; do
            flags=""
            case "$V" in
                prealloc) flags="--prealloc" ;;
                fp)       flags="--fp" ;;
            esac
            OUT=$LOGDIR/L${L}_${M}_${V}.log
            $BENCH --level=$L --iterations=$ITER --mode=$M $flags > $OUT 2>&1

            # Parse "Total time: HNF=Xms  MLLL=Yms  GRAM=Zms"
            line=$(grep "^Total time:" $OUT || true)
            gram_ms=$(echo "$line" | sed -n 's/.*GRAM=\([0-9.]*\)ms.*/\1/p')
            mlll_ms=$(echo "$line" | sed -n 's/.*MLLL=\([0-9.]*\)ms.*/\1/p')
            printf '%-4s %-5s %-10s %-6s %-9s %-9s\n' "L$L" "$M" "$V" "$ITER" "${gram_ms:-?}" "${mlll_ms:-?}"
        done
    done
done

echo
echo "Logs: $LOGDIR"
