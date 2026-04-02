# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 이데알 곱셈(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 랜덤 O0-이데알 격자 쌍에 대해 비교.

## 구현 현황

### 현재 상태

기존 SQIsign의 HNF 기반 이데알 연산(`quat_lattice_mul`, `quat_lattice_add`)을 **교체하지 않고**, MLLL 기반 함수를 별도로 구현하여 **동일 입력에 대해 두 방식의 결과와 중간값 비트 크기를 비교**한 상태이다.

### 파일 구조

| 파일 | 설명 |
|------|------|
| `mlll.c` | **핵심 구현**. 논문 Algorithm 1 (MLLL 알고리즘) + Algorithm 2 (CompactIdealMultiplication). `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()` 세 함수를 제공한다. |
| `mlll_internals.h` | MLLL 함수 선언 및 상수 정의 (`MLLL_MAX_GENERATORS=16` 등) |
| `mlll_tests.c` | **정확성 검증**. 같은 입력에 대해 HNF(`quat_lattice_mul`)와 MLLL(`quat_lattice_mul_mlll`)을 모두 실행하고, `quat_lattice_equal()`로 결과 격자가 동일한지 확인한다. 총 4개 테스트: (1) HNF vs MLLL 결과 비교, (2) MLLL 출력의 LLL-reducedness 검증, (3) 선형 종속 벡터 처리, (4) CompactIdealMultiplication 전체 흐름 |
| `mlll_benchmark.c` | **중간값 비트 크기 비교**. `bitsize_tracker`를 이용해 HNF/MLLL 실행 중 나오는 모든 정수의 최대 비트 크기를 추적한다. `tracker_reset()` → 연산 실행 → `tracker_get_max()`로 최대값 수집. |
| `bitsize_tracker.h` | 벤치마크 전용 계측 도구. `mlll.c`와 `hnf.c` 내부에 `tracker_update_vec4()` 호출이 삽입되어 있어, 연산 중 중간값의 비트 크기를 글로벌 변수로 추적한다. 일반 빌드에서는 `BITSIZE_TRACKER_ENABLE` 미정의 시 no-op. |
| `test_mlll_only.c` | MLLL 테스트만 독립 실행하는 진입점 (main 함수) |

### 비교 방법

1. **정확성**: `mlll_tests.c`에서 동일 격자 쌍에 대해 `quat_lattice_mul()`(HNF)과 `quat_lattice_mul_mlll()`(MLLL)을 호출하고, `quat_lattice_equal()`로 결과가 같은 격자인지 비교
2. **중간값 비트 크기**: `mlll_benchmark.c`의 `bench_one_mul()`에서 `bitsize_tracker`를 켜고 각 방식을 실행하여 연산 중 등장하는 정수의 최대 비트 크기를 기록

### 논문 알고리즘 매핑

| 논문 | 코드 |
|------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` in `mlll.c:177-571` |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` in `mlll.c:576-617` |
| Lemma 1 (중간값 바운드 ≤ max\|\|a_i\|\|²) | `mlll_benchmark.c`로 실험적 검증 |

### 아직 안 한 것

- SQIsign 파이프라인에서 HNF → MLLL 실제 교체 (현재 `lattice.c:177`의 `quat_lattice_mul`은 여전히 HNF 사용)
- 논문 Algorithm 3 (RandomIdealGivenPrimeNorm) 구현
- 논문 Algorithm 4 (RandomEquivalentPrimeIdeal) 수정

## 핵심 지표: 최대 중간값 비트 크기

논문의 핵심 주장은 MLLL이 중간 정수 크기를 `max ||a_i||^2` (입력 노름의 제곱)으로 바운드하는 반면, HNF는 그보다 훨씬 커질 수 있다는 것이다. 이것이 **고정 정밀도(fixed-precision)** 연산의 가능 여부를 결정한다.

### NIST Level 1 (p ~ 2^253, 노름 비트 크기 = 127, 10회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **2162** | **255** | **0.118** |
| 평균 중간값 비트 | 2073 | 253 | 0.122 |
| 최대 출력 비트 | 128 | 127 | ~1.0 |
| 총 실행 시간 | 1.16 ms | 4102 ms | x3527 |

### NIST Level 3 (p ~ 2^381, 노름 비트 크기 = 193, 5회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **3290** | **388** | **0.118** |
| 평균 중간값 비트 | 3265 | 385 | 0.118 |
| 최대 출력 비트 | 194 | 194 | ~1.0 |
| 총 실행 시간 | 0.74 ms | 10269 ms | x13914 |

### NIST Level 5 (p ~ 2^509, 노름 비트 크기 = 254, 3회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **4325** | **510** | **0.118** |
| 평균 중간값 비트 | 3975 | 508 | 0.128 |
| 최대 출력 비트 | 255 | 254 | ~1.0 |
| 총 실행 시간 | 0.51 ms | 6511 ms | x12743 |

## 분석

### 중간값 비트 크기

- **HNF 중간값은 입력 비트 크기의 ~17배까지 증가** (예: 128비트 입력 -> 2163비트 중간값). HNF 계산에서 4x4 행렬식을 모듈러스로 사용하기 때문.
- **MLLL 중간값은 입력 비트 크기의 ~2배로 유지** (예: 128비트 입력 -> 255비트 중간값). 논문 Lemma 1의 `max ||a_i||^2` 바운드와 일치.
- MLLL/HNF 비율은 모든 보안 레벨에서 일관되게 **~0.118**.

### 출력 비트 크기

HNF와 MLLL의 출력 basis 원소 크기는 거의 동일 -- 둘 다 같은 격자를 표현하며, 형태만 다름 (상삼각 행렬 vs LLL-reduced).

### 실행 시간

이 GMP 기반(임의 정밀도) 구현에서는 MLLL이 HNF보다 상당히 느림:
1. MLLL은 정확한 유리수 Gram-Schmidt(`ibq_t` 분수)를 사용하므로 비용이 큼.
2. HNF는 각 단계에서 모듈러 축약으로 중간값을 바운드하므로 GMP에서 유리함.

하지만 속도 비교가 핵심이 아님. 논문이 목표로 하는 것은 **고정 정밀도** SQIsign 구현:
- HNF의 2000-4000비트 중간값은 고정 정밀도 예산에 **들어가지 않음**.
- MLLL의 ~500비트 중간값은 **들어감**, 고정 정밀도 연산 가능.

## 재현 방법

```bash
# 빌드
cd build_wsl
cmake .. -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release
make sqisign_bm_mlll -j4

# 실행
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=10
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=5
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=3
```

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
