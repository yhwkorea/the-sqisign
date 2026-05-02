#!/bin/bash
# Phase 2 fp capacity verification: 공간 딱 정해놓고 1k iter 돌려서
# 실제 peak bitsize가 선언 budget을 초과하는지 확인.
#
# Declared budget (overflow trap, 2026-04-22 정정): nwords*64 - 2 bits
# (부호 비트 1 + add/sub carry 비트 1).
#   L1 vec: 5*64-2  = 318   Gram: 9*64-2  = 574
#   L3 vec: 7*64-2  = 446   Gram: 13*64-2 = 830
#   L5 vec: 9*64-2  = 574   Gram: 17*64-2 = 1086
#
# 초안 공식 (nwords-1)*64 는 "top limb 전체를 sign-extension margin으로
# 예약"이라는 잘못된 narrative. Phase 1 peak(L1 vec 259 등)가 이미 top
# limb의 일부 비트를 쓰고 있었으므로 margin 전제가 성립하지 않는다.
# 상세: FIXED_PRECISION_DECISION_KO.md §6.4-c.
#
# Pass criteria:
#   - trap abort 0건 (실패 시 stderr + abort로 프로세스 종료)
#   - tracker vec/gram max ≤ 선언 budget
set -e

BUILD=/mnt/c/Users/htelr/the-sqisign/build
BENCH=$BUILD/src/quaternion/ref/generic/test/sqisign_bm_mlll
LOGDIR=/mnt/c/Users/htelr/the-sqisign/src/quaternion/ref/generic/lll/bench_logs/p2_fp_capacity_verify

ITER_ALG2=${ITER_ALG2:-500}
ITER_ALG3=${ITER_ALG3:-5000}

if [ ! -x "$BENCH" ]; then
    echo "ERROR: benchmark binary not found at $BENCH"
    exit 1
fi

mkdir -p "$LOGDIR"

declare -A VEC_BUDGET=(
    [1]=318 [3]=446 [5]=574
)
declare -A GRAM_BUDGET=(
    [1]=574 [3]=830 [5]=1086
)

printf '%-4s %-5s %-6s %-10s %-10s %-10s %-10s %-8s\n' \
    "Lvl" "Mode" "Iter" "vec_max" "vec_bdgt" "gram_max" "gram_bdgt" "Verdict"
printf '%-4s %-5s %-6s %-10s %-10s %-10s %-10s %-8s\n' \
    "---" "----" "----" "-------" "--------" "--------" "---------" "-------"

OVERALL_PASS=1

for L in 1 3 5; do
    for M in alg2 alg3; do
        if [ "$M" = "alg2" ]; then
            ITER=$ITER_ALG2
        else
            ITER=$ITER_ALG3
        fi
        OUT=$LOGDIR/L${L}_${M}.log
        # fp ON, trap ON (default)
        $BENCH --level=$L --iterations=$ITER --mode=$M --fp > "$OUT" 2>&1 || {
            echo "  L$L $M: ABORTED (trap likely fired)"
            grep "quat_fp_.*overflow" "$OUT" | head -3
            OVERALL_PASS=0
            continue
        }

        # Parse "GRAM vec coords: ... max=X" and "GRAM Gram entries: ... max=Y"
        vec_max=$(grep "GRAM vec coords:" "$OUT" | sed -n 's/.*max=\([0-9]*\).*/\1/p')
        gram_max=$(grep "GRAM Gram entries:" "$OUT" | sed -n 's/.*max=\([0-9]*\).*/\1/p')

        vb=${VEC_BUDGET[$L]}
        gb=${GRAM_BUDGET[$L]}

        verdict="PASS"
        if [ -z "$vec_max" ] || [ -z "$gram_max" ]; then
            verdict="PARSE?"
            OVERALL_PASS=0
        elif [ "$vec_max" -gt "$vb" ] || [ "$gram_max" -gt "$gb" ]; then
            verdict="FAIL"
            OVERALL_PASS=0
        fi

        printf '%-4s %-5s %-6s %-10s %-10s %-10s %-10s %-8s\n' \
            "L$L" "$M" "$ITER" "${vec_max:-?}" "$vb" "${gram_max:-?}" "$gb" "$verdict"
    done
done

echo
if [ "$OVERALL_PASS" = "1" ]; then
    echo "=== VERIFICATION PASS ==="
    echo "Logs: $LOGDIR"
    exit 0
else
    echo "=== VERIFICATION FAIL ==="
    echo "Logs: $LOGDIR"
    exit 1
fi
