#!/bin/bash
# Phase 2 primary-decision sweep: ibz / prealloc / fp across L1/L3/L5 × Alg2/Alg3.
# Driven from WSL Ubuntu because the build cache is Linux-toolchain.
set -e

BUILD=/mnt/c/Users/htelr/the-sqisign/build
BENCH=$BUILD/src/quaternion/ref/generic/test/sqisign_bm_mlll
LOGDIR=/mnt/c/Users/htelr/the-sqisign/src/quaternion/ref/generic/lll/bench_logs/p2_phase2_sweep

ITER=${ITER:-30}

if [ ! -x "$BENCH" ]; then
    echo "ERROR: benchmark binary not found at $BENCH"
    echo "Rebuild with: cmake --build $BUILD --target sqisign_bm_mlll"
    exit 1
fi

printf '%-4s %-5s %-10s %-9s %-9s\n' "Lvl" "Mode" "Variant" "GRAMms" "MLLLms"
printf '%-4s %-5s %-10s %-9s %-9s\n' "---" "----" "-------" "------" "------"

for L in 1 3 5; do
    for M in alg2 alg3; do
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
            printf '%-4s %-5s %-10s %-9s %-9s\n' "L$L" "$M" "$V" "${gram_ms:-?}" "${mlll_ms:-?}"
        done
    done
done

echo
echo "Logs: $LOGDIR"
