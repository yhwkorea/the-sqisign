# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 ideal multiplication(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 random O0-ideal lattice pair에 대해 비교.

## 구현 현황

### 현재 상태

기존 SQIsign의 HNF 기반 ideal operation(`quat_lattice_mul`, `quat_lattice_add`)을 **교체하지 않고**, MLLL 기반 함수를 별도로 구현하여 **동일 입력에 대해 두 방식의 결과와 intermediate bit size를 비교**한 상태이다.

### 파일 구조

| 파일 | 설명 |
|------|------|
| `mlll.c` | **핵심 구현**. 논문 Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()` 세 함수를 제공한다. |
| `mlll_internals.h` | MLLL 함수 선언 및 상수 정의 (`MLLL_MAX_GENERATORS=16` 등) |
| `mlll_tests.c` | **정확성 검증**. 같은 입력에 대해 HNF(`quat_lattice_mul`)와 MLLL(`quat_lattice_mul_mlll`)을 모두 실행하고, `quat_lattice_equal()`로 결과 lattice가 동일한지 확인한다. 총 4개 테스트: (1) HNF vs MLLL 결과 비교, (2) MLLL output의 LLL-reducedness 검증, (3) linearly dependent generator 처리, (4) CompactIdealMultiplication 전체 흐름 |
| `mlll_benchmark.c` | **intermediate bit size 비교**. `bitsize_tracker`를 이용해 HNF/MLLL 실행 중 나오는 모든 integer의 최대 bit size를 추적한다. `tracker_reset()` → 연산 실행 → `tracker_get_max()`로 최대값 수집. |
| `bitsize_tracker.h` | 벤치마크 전용 계측 도구. `mlll.c`와 `hnf.c` 내부에 `tracker_update_vec4()` 호출이 삽입되어 있어, 연산 중 intermediate value의 bit size를 global 변수로 추적한다. 일반 빌드에서는 `BITSIZE_TRACKER_ENABLE` 미정의 시 no-op. |
| `test_mlll_only.c` | MLLL 테스트만 독립 실행하는 entry point (main 함수) |

### 비교 방법

1. **정확성**: `mlll_tests.c`에서 동일 lattice pair에 대해 `quat_lattice_mul()`(HNF)과 `quat_lattice_mul_mlll()`(MLLL)을 호출하고, `quat_lattice_equal()`로 결과가 같은 lattice인지 비교
2. **intermediate bit size**: `mlll_benchmark.c`의 `bench_one_mul()`에서 `bitsize_tracker`를 켜고 각 방식을 실행하여 연산 중 등장하는 integer의 최대 bit size를 기록

### 논문 알고리즘 매핑

| 논문 | 코드 |
|------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` in `mlll.c:177-571` |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` in `mlll.c:576-617` |
| Lemma 1 (intermediate bound ≤ max\|\|a_i\|\|²) | `mlll_benchmark.c`로 실험적 검증 |

### 아직 안 한 것

- SQIsign pipeline에서 HNF → MLLL 실제 교체 (현재 `lattice.c:177`의 `quat_lattice_mul`은 여전히 HNF 사용)
- 논문 Algorithm 3 (RandomIdealGivenPrimeNorm) 구현
- 논문 Algorithm 4 (RandomEquivalentPrimeIdeal) 수정

## 핵심 지표: 최대 intermediate bit size

논문의 핵심 주장은 MLLL이 intermediate integer size를 `max ||a_i||^2` (입력 norm의 제곱)으로 bound하는 반면, HNF는 그보다 훨씬 커질 수 있다는 것이다. 이것이 **fixed-precision** 연산의 가능 여부를 결정한다.

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 10회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **2162** | **255** | **0.118** |
| 평균 intermediate bits | 2073 | 253 | 0.122 |
| 최대 output bits | 128 | 127 | ~1.0 |
| 총 실행 시간 | 1.16 ms | 4102 ms | x3527 |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 5회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **3290** | **388** | **0.118** |
| 평균 intermediate bits | 3265 | 385 | 0.118 |
| 최대 output bits | 194 | 194 | ~1.0 |
| 총 실행 시간 | 0.74 ms | 10269 ms | x13914 |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 3회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 intermediate bits** | **4325** | **510** | **0.118** |
| 평균 intermediate bits | 3975 | 508 | 0.128 |
| 최대 output bits | 255 | 254 | ~1.0 |
| 총 실행 시간 | 0.51 ms | 6511 ms | x12743 |

## 분석

### Intermediate Bit Size

- **HNF intermediate value는 input bit size의 ~17배까지 증가** (예: 128-bit input → 2162-bit intermediate). HNF 계산에서 4x4 determinant를 modulus로 사용하기 때문.
- **MLLL intermediate value는 input bit size의 ~2배로 유지** (예: 128-bit input → 255-bit intermediate). 논문 Lemma 1의 `max ||a_i||^2` bound와 일치.
- MLLL/HNF 비율은 모든 security level에서 일관되게 **~0.118**.

### Output Bit Size

HNF와 MLLL의 output basis element size는 거의 동일 — 둘 다 같은 lattice를 표현하며, 형태만 다름 (upper triangular vs LLL-reduced).

### 실행 시간

이 GMP 기반(arbitrary precision) 구현에서는 MLLL이 HNF보다 상당히 느림:
1. MLLL은 exact rational Gram-Schmidt(`ibq_t` fraction)를 사용하므로 비용이 큼.
2. HNF는 각 단계에서 modular reduction으로 intermediate value를 bound하므로 GMP에서 유리함.

하지만 속도 비교가 핵심이 아님. 논문이 목표로 하는 것은 **fixed-precision** SQIsign 구현:
- HNF의 2000-4000 bit intermediate는 fixed-precision budget에 **들어가지 않음**.
- MLLL의 ~500 bit intermediate는 **들어감**, fixed-precision 연산 가능.

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
# 결과: 35/35 통과 (2026-04-02 확인)
```

## 날짜

2026-03-30 (초기 벤치마크), 2026-04-02 (전체 테스트 통과 확인, 실행 시간 포함 벤치마크 재실행)
