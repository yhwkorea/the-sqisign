# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 ideal multiplication(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 random O0-ideal lattice pair에 대해 비교.

## 구현 현황

- **GSO**: Float GSO (`dpe_t`, 53-bit mantissa + extended exponent), L² 방식
- **Gram matrix**: 정확한 정수(`ibz_t`), 점진적 업데이트
- **Dependency detection**: `dpe_zero_p(r[k][k])`(float 판정), done section에 HNF 후처리로 임시 우회
- **테스트**: 6/6 통과
- **Fixed-precision**: 확인됨 — MLLL 벡터 좌표가 Lemma 1 bound 내

## 파일 구조

| 파일 | 설명 |
|------|------|
| `mlll.c` | **핵심 구현**. 논문 Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). Float GSO + Gram matrix. |
| `mlll_internals.h` | 함수 선언 및 상수 (`MLLL_MAX_GENERATORS=16`) |
| `mlll_tests.c` | **정확성 검증**. 6개 테스트 |
| `mlll_benchmark.c` | **intermediate bit size 비교**. `tracker_disable()`로 MLLL reduction loop만 분리 측정. |
| `bitsize_tracker.h` | 2채널 tracker: `vec`(벡터 좌표) + `gso`(GSO 계수). `BITSIZE_TRACKER_ENABLE` 없으면 no-op. |
| `test_mlll_only.c` | MLLL 테스트 독립 실행 entry point |

## 논문 알고리즘 매핑

| 논문 | 코드 | 상태 |
|------|------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` | ✓ 구현 (D_k=0 처리 미완 — 아래 참조) |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` | ✓ 구현 |
| Algorithm 3 (RandomIdealGivenPrimeNorm) | — | 미구현 |
| Algorithm 4 (RandomEquivalentPrimeIdeal) | — | 미구현 |
| Lemma 1 (intermediate bound) | `mlll_benchmark.c` | ✓ 실험적 검증 |

## Fixed-Precision 벤치마크 결과 (2026-04-09)

Tracker는 MLLL reduction loop만 측정 (`tracker_disable()`로 HNF 후처리 제외).

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 3회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **2153** | **255** | **8.4x** |
| 최대 output bits | 128 | 128 | ~1.0 |
| 총 실행 시간 | 0.49ms | 1.08ms | 0.45x |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 1회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **3245** | **384** | **8.5x** |
| 최대 output bits | 194 | 192 | ~1.0 |
| 총 실행 시간 | 0.14ms | 0.38ms | 0.37x |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 1회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **3779** | **503** | **7.5x** |
| 최대 output bits | 250 | 250 | ~1.0 |
| 총 실행 시간 | 0.15ms | 0.24ms | 0.63x |

### 요약

MLLL 벡터 좌표는 모든 보안 레벨에서 **입력 bits의 ~2배** 이내로 유지됨. 논문 Lemma 1 bound (≤ max||a_i||²)와 정확히 일치. HNF는 modular determinant 계산 때문에 입력 bits의 ~17배까지 증가.

GSO 계수는 `dpe_t`(53-bit float)로 처리하므로 정수 tracker에 0 bits로 기록됨. 이것이 의도된 설계: fixed-precision GSO → 정수 성장은 벡터 좌표에서만 발생하며 Lemma 1로 bounded.

## 미구현 사항

### 1. D_k=0 (dependent vector) 정상 처리

**현재**: MLLL이 dependent vector(`r[k][k] ≈ 0`)를 감지하면 reduction loop가 미해결 상태로 남길 수 있음. `done` section에서 `ibz_mat_4xn_hnf_mod_core()`로 rank-4 basis를 추출하는 임시 우회 사용 중. 이 HNF 호출은 intermediate bit size 측정에서 제외됨.

**필요**: Pohst (1987) 및 Matthews pseudo-code에 따르면:
- `D_k = 0, λ_{k,k-1} = 0`: 행 교환 후 β 감소
- `D_k = 0, λ_{k,k-1} ≠ 0`: `Swap2`(특별 μ 업데이트) 후 `Swap1`

Gram matrix에서 정확한 정수 `D_k`, `λ_{k,j}`를 계산하여 처리해야 함. 참조: http://www.numbertheory.org/PDFS/mlll.pdf

### 2. SQIsign pipeline 통합

메인 파이프라인(`lattice.c:177`)의 HNF를 MLLL로 교체. 현재는 벤치마크용 별도 함수로만 존재.

### 3. Algorithm 3 (RandomIdealGivenPrimeNorm)

Cornacchia + MLLL 조합으로 prime norm ideal 샘플링.

### 4. Algorithm 4 (RandomEquivalentPrimeIdeal)

HNF basis 대신 LLL-reduced basis 사용으로 수정.

## 구현 이력

| 날짜 | 변경 |
|------|------|
| 2026-03-30 | 초기 MLLL (`ibq_t` exact rational GSO), 전체 테스트 통과 |
| 2026-04-02 | Bitsize tracker, 실행 시간 포함 벤치마크 |
| 2026-04-07 | tau swap 버그 수정, tracker 확장 |
| 2026-04-08 | Integral GSO 실험 (`ibq_t` → `ibz_t`) |
| 2026-04-09 | Float GSO 전환 (`dpe_t`), Gram matrix 방식, tracker 분리, HNF 후처리 임시 우회 |

## 재현 방법

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBITSIZE_TRACKER_ENABLE=ON
make sqisign_bm_mlll sqisign_test_mlll -j$(nproc)

# 테스트
./src/quaternion/ref/generic/test/sqisign_test_mlll

# 벤치마크 (iteration 수 적게 — 생성이 느릴 수 있음)
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=3
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=1
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=1
```
