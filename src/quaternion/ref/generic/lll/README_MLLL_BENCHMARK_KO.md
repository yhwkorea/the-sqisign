# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 ideal multiplication(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 random O0-ideal lattice pair에 대해 비교.

## 구현 현황

- **GSO**: 정확한 정수 d[]/lam[][] (Cohen 2.6.7 방식). d[k+1] = det(Gram(b[0..k])), lam[i][j] = d[j+1]*μ[i][j].
- **Gram matrix**: 정확한 정수(`ibz_t`), size-reduce 시 점진적 업데이트, 행 이동 후 재계산
- **Dependency detection**: `d[m+1] == 0` (정확한 정수 테스트). dependent vector는 forced-swap cascade + REDI로 제거. HNF 후처리 불필요.
- **테스트**: 6/6 통과
- **Fixed-precision**: 확인됨 — MLLL 벡터 좌표가 입력의 ~3배 이내

## 파일 구조

| 파일 | 설명 |
|------|------|
| `mlll.c` | **핵심 구현**. 논문 Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). 정확한 정수 d[]/lam[][], swap마다 전체 재계산. |
| `mlll_internals.h` | 함수 선언 및 상수 (`MLLL_MAX_GENERATORS=16`) |
| `mlll_tests.c` | **정확성 검증**. 6개 테스트 |
| `mlll_benchmark.c` | **intermediate bit size 비교**. `bitsize_tracker`로 계산 중 최대 bit size 측정. |
| `bitsize_tracker.h` | 2채널 tracker: `vec`(벡터 좌표) + `gso`(GSO 계수). `BITSIZE_TRACKER_ENABLE` 없으면 no-op. |
| `test_mlll_only.c` | MLLL 테스트 독립 실행 entry point |

## 논문 알고리즘 매핑

| 논문 | 코드 | 상태 |
|------|------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` | 완료 (정확한 정수 d/lam, D_k=0 정상 처리) |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` | 완료 |
| Algorithm 3 (RandomIdealGivenPrimeNorm) | — | 미구현 |
| Algorithm 4 (RandomEquivalentPrimeIdeal) | — | 미구현 |
| Lemma 1 (intermediate bound) | `mlll_benchmark.c` | 실험적 검증 완료 |

## D_k=0 처리 (Matthews/Pohst)

`d[m+1] = 0` (dependent vector 감지, 정확한 정수 테스트):

1. **Forced swap**: b[m]과 b[m-1] 교환, tau 장벽 무시 (`m = max(m-1, 1)`)
2. **Cascade**: dependent vector가 위치 1까지 내려가며, REDI로 더 짧은 벡터와 비교하여 제거
3. **Recompute**: 매 swap 후 Gram matrix에서 d[]/lam[][] 전체 재계산 (O(beta^2), beta ≤ 16이므로 무시 가능)
4. **최종 reduction**: dependent vector 제거 후 m=1부터 전체 LLL reduction 재실행

이전 HNF 후처리 우회를 수학적으로 정확한 방식으로 대체.

## 정확한 정수 벤치마크 결과 (2026-04-09)

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 3회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **2160** | **389** | **5.5x** |
| 최대 output bits | 128 | 127 | ~1.0 |
| 총 실행 시간 | 0.46ms | 8.17ms | 0.06x |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 1회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **3238** | **571** | **5.7x** |
| 최대 output bits | 189 | 189 | ~1.0 |
| 총 실행 시간 | 0.14ms | 4.74ms | 0.03x |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 1회)

| 지표 | HNF | MLLL vec | 개선 |
|------|-----|----------|------|
| **최대 intermediate bits** | **4307** | **758** | **5.7x** |
| 최대 output bits | 252 | 252 | ~1.0 |
| 총 실행 시간 | 0.15ms | 9.77ms | 0.02x |

### 요약

MLLL 벡터 좌표는 모든 보안 레벨에서 **입력 bits의 ~3배** 이내로 유지됨. 이전 float GSO 버전의 ~2배보다 높지만, HNF 후처리 없이 정확성을 보장하는 대가.

HNF는 modular determinant 계산으로 입력의 ~17배까지 증가하므로, MLLL의 intermediate bit 개선은 여전히 **5.5-5.7배**.

### 성능 참고

recompute-from-scratch 방식 (매 swap마다 d[]/lam[][] 전체 재계산)은 이전 float GSO 버전보다 느림:
- Float GSO: dpe_t 연산 O(1), swap 횟수 적음
- 정확한 정수: 재계산 O(beta^2), dependent vector cascade로 forced swap 다수

향후 최적화: Matthews의 증분 Swap1/Swap2 공식으로 swap당 비용을 O(beta^2) → O(beta)로 감소 가능.

## 미구현 사항

### 1. SQIsign pipeline 통합
메인 파이프라인(`lattice.c:177`)의 HNF를 MLLL로 교체.

### 2. 증분 Swap1/Swap2 (성능 최적화)
recompute-from-scratch를 Matthews의 증분 업데이트 공식으로 대체.

### 3. Algorithm 3 (RandomIdealGivenPrimeNorm)
Cornacchia + MLLL 조합으로 prime norm ideal 샘플링.

### 4. Algorithm 4 (RandomEquivalentPrimeIdeal)
HNF basis 대신 LLL-reduced basis 사용.

## 구현 이력

| 날짜 | 변경 |
|------|------|
| 2026-03-30 | 초기 MLLL (`ibq_t` exact rational GSO), 전체 테스트 통과 |
| 2026-04-02 | Bitsize tracker, 실행 시간 포함 벤치마크 |
| 2026-04-07 | tau swap 버그 수정, tracker 확장 |
| 2026-04-08 | Integral GSO 실험 (`ibq_t` → `ibz_t`) |
| 2026-04-09 | Float GSO 전환 (`dpe_t`), Gram matrix 방식, HNF 후처리 임시 우회 |
| 2026-04-09 | **정확한 정수 d[]/lam[][] 재작성**. D_k=0 정상 처리 (forced-swap cascade). HNF 후처리 제거. 전체 재계산 방식. 6/6 테스트 통과. |

## 재현 방법

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBITSIZE_TRACKER_ENABLE=ON
make sqisign_bm_mlll sqisign_test_mlll -j$(nproc)

# 테스트
./src/quaternion/ref/generic/test/sqisign_test_mlll

# 벤치마크
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=3
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=1
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=1
```
