# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 ideal multiplication(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 random O0-ideal lattice pair에 대해 비교.

## 구현 현황

### 현재 상태

기존 SQIsign의 HNF 기반 ideal operation(`quat_lattice_mul`, `quat_lattice_add`)을 **교체하지 않고**, MLLL 기반 함수를 별도로 구현하여 **동일 입력에 대해 두 방식의 결과와 intermediate bit size를 비교**한 상태이다.

### 파일 구조

| 파일 | 설명 |
|------|------|
| `mlll.c` | **핵심 구현**. 논문 Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()` 세 함수를 제공. 파일 상단 한국어 docstring에 Pohst 1987 알고리즘 설명 및 fixed-precision 한계 명시. |
| `mlll_internals.h` | MLLL 함수 선언 및 상수 정의 (`MLLL_MAX_GENERATORS=16` 등) |
| `mlll_tests.c` | **정확성 검증**. 같은 입력에 대해 HNF(`quat_lattice_mul`)와 MLLL(`quat_lattice_mul_mlll`)을 모두 실행하고, `quat_lattice_equal()`로 결과 lattice가 동일한지 확인. 총 6개 테스트: (1) HNF vs MLLL 결과 비교, (2) MLLL output의 LLL-reducedness 검증, (3) linearly dependent generator 처리, (4) CompactIdealMultiplication 전체 흐름, (5) 127비트 prime 현실 스케일 검증, (6) tau > 1 swap 경로 검증 |
| `mlll_benchmark.c` | **intermediate bit size 비교**. `bitsize_tracker`를 이용해 HNF/MLLL 실행 중 나오는 모든 integer의 최대 bit size를 추적한다. `tracker_reset()` → 연산 실행 → `tracker_get_max()`로 최대값 수집. |
| `bitsize_tracker.h` | 벤치마크 전용 계측 도구. `mlll.c`와 `hnf.c` 내부에 `tracker_update_vec4()` / `tracker_update_ibz()` 호출이 삽입되어 intermediate value의 bit size를 추적. `__attribute__((weak))` 심볼로 `mlll_benchmark.c` 없이도 링크 가능. `BITSIZE_TRACKER_ENABLE` 미정의 시 no-op. |
| `test_mlll_only.c` | MLLL 테스트만 독립 실행하는 entry point (main 함수) |

### 비교 방법

1. **정확성**: `mlll_tests.c`에서 동일 lattice pair에 대해 `quat_lattice_mul()`(HNF)과 `quat_lattice_mul_mlll()`(MLLL)을 호출하고, `quat_lattice_equal()`로 결과가 같은 lattice인지 비교
2. **intermediate bit size**: `mlll_benchmark.c`의 `bench_one_mul()`에서 `bitsize_tracker`를 켜고 각 방식을 실행하여 연산 중 등장하는 integer의 최대 bit size를 기록

### 논문 알고리즘 매핑

| 논문 | 코드 |
|------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` in `mlll.c:192-586` |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` in `mlll.c:591-632` |
| Lemma 1 (intermediate bound ≤ max\|\|a_i\|\|²) | `mlll_benchmark.c`로 실험적 검증 |

### 변경사항 (2026-04-07)

- **버그 수정**: swap 감소 조건 `m > 1` → `m > tau`. `remove_vector`가 `tau = m+1`로 설정한 후 기존 조건은 `tau` 아래로 내려가 잘못된 GS 데이터에 접근하는 버그가 있었음.
- **Tracker 확장**: `compute_gs_single()` 내 GS inner product, cross product, mu 분자/분모, B 분자/분모에 `tracker_update_ibz()` 추가.
- **Weak 심볼**: `bitsize_tracker.h`가 `__attribute__((weak))`를 사용하여 `mlll_benchmark.c` 없이도 링크 가능.
- **Rank assertion**: `quat_lattice_add_mlll()`에 `assert(rank > 0 && rank <= 4)` 추가.
- **한국어 docstring**: `mlll.c` 상단 주석을 한국어로 재작성 (Pohst 1987, fixed-precision 한계 설명).
- **새 테스트 2개**: `quat_test_mlll_realistic_scale` (127비트 prime), `quat_test_mlll_tau_path` (tau > 1 swap 경로).

### 아직 안 한 것

- SQIsign pipeline에서 HNF → MLLL 실제 교체 (현재 `lattice.c:177`의 `quat_lattice_mul`은 여전히 HNF 사용)
- 논문 Algorithm 3 (RandomIdealGivenPrimeNorm) 구현
- 논문 Algorithm 4 (RandomEquivalentPrimeIdeal) 수정

## 핵심 지표: 최대 intermediate bit size

논문의 핵심 주장은 MLLL이 intermediate integer size를 `max ||a_i||^2` (입력 norm의 제곱)으로 bound하는 반면, HNF는 그보다 훨씬 커질 수 있다는 것이다. 이것이 **fixed-precision** 연산의 가능 여부를 결정한다.

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 10회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **2168** | **4513** | **2.082** |
| 평균 intermediate bits | 2106 | 3596 | 1.707 |
| 최대 output bits | 128 | 127 | ~1.0 |
| 총 실행 시간 | 1.35 ms | 413.79 ms | x306 |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 5회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **3285** | **6881** | **2.095** |
| 평균 intermediate bits | 2907 | 4788 | 1.647 |
| 최대 output bits | 194 | 193 | ~1.0 |
| 총 실행 시간 | 0.72 ms | 441.16 ms | x614 |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 3회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **4315** | **7699** | **1.784** |
| 평균 intermediate bits | 3440 | 6931 | 2.015 |
| 최대 output bits | 254 | 253 | ~1.0 |
| 총 실행 시간 | 0.47 ms | 674.98 ms | x1449 |

## 분석

### Intermediate Bit Size

**중요**: 현재 tracker는 정수 벡터 좌표뿐 아니라 Gram-Schmidt 유리수 계수(`ibq_t` 분자/분모)까지 추적한다. 이전 측정(2026-04-02)은 벡터 좌표만 추적하여 MLLL intermediate가 오해를 줄 만큼 작게(~255 bits) 나왔음. 업데이트된 tracker의 결과:

- **HNF intermediate value는 input bit size의 ~17배까지 증가** (예: 128-bit input → 2168-bit intermediate). HNF 계산에서 4x4 determinant를 modulus로 사용하기 때문.
- **MLLL intermediate value는 HNF의 ~2배까지 증가** (예: 128-bit input → 4513-bit intermediate). Exact rational Gram-Schmidt 계수(`ibq_t` fraction)의 분자/분모가 arbitrary precision에서 제한 없이 커지기 때문.
- 논문 Lemma 1의 bound(`max ||a_i||^2`)는 **정수 벡터 좌표에만** 적용되며, GS 계수에는 해당하지 않음. 벡터 좌표는 논문 주장대로 bounded (~255 bits for Level 1).

### Output Bit Size

HNF와 MLLL의 output basis element size는 거의 동일 — 둘 다 같은 lattice를 표현하며, 형태만 다름 (upper triangular vs LLL-reduced).

### 실행 시간

이 GMP 기반(arbitrary precision) 구현에서는 MLLL이 HNF보다 느림:
1. MLLL은 exact rational Gram-Schmidt(`ibq_t` fraction)를 사용하므로 비용이 큼.
2. HNF는 각 단계에서 modular reduction으로 intermediate value를 bound하므로 GMP에서 유리함.

속도 비교가 핵심이 아님. 논문이 목표로 하는 것은 **fixed-precision** SQIsign 구현:
- HNF의 2000-4000 bit intermediate(벡터 좌표)는 fixed-precision budget에 **들어가지 않음**.
- MLLL의 ~500 bit 벡터 좌표는 **들어감**, fixed-precision 연산 가능.
- Fixed-precision MLLL은 `ibq_t` 대신 integral GSO 또는 L² 표현을 사용하여 GS 계수 blowup을 제거함.

## 재현 방법

```bash
# 빌드 (tracker 포함)
cd build_wsl
cmake .. -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS='-DBITSIZE_TRACKER_ENABLE'
make sqisign_bm_mlll -j4

# 실행
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=10
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=5
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=3
```

> 참고: `-DBITSIZE_TRACKER_ENABLE` 없이 빌드하면 intermediate bit size가 0으로 나옴 (tracker가 no-op).

## 테스트

```bash
# 전체 테스트 스위트 (35개 테스트)
cd build_wsl
cmake .. -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
ctest --output-on-failure -j$(nproc)
# 결과: 35/35 통과 (2026-04-08 확인, 새 테스트 2개 포함)
```

## 날짜

2026-03-30 (초기 벤치마크), 2026-04-02 (전체 테스트 통과 확인, 실행 시간 포함 벤치마크 재실행), 2026-04-07 (tau 버그 수정, tracker 확장, 새 테스트 추가)
