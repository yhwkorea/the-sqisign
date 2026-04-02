# MLLL vs HNF 벤치마크 결과

"Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - 고려대학교) 논문 구현.

HNF 기반 vs MLLL 기반 이데알 곱셈(`quat_lattice_mul` vs `quat_lattice_mul_mlll`)을 랜덤 O0-이데알 격자 쌍에 대해 비교.

## 핵심 지표: 최대 중간값 비트 크기

논문의 핵심 주장은 MLLL이 중간 정수 크기를 `max ||a_i||^2` (입력 노름의 제곱)으로 바운드하는 반면, HNF는 그보다 훨씬 커질 수 있다는 것이다. 이것이 **고정 정밀도(fixed-precision)** 연산의 가능 여부를 결정한다.

### NIST Level 1 (p ~ 2^253, 노름 비트 크기 = 127, 10회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **2163** | **255** | **0.118** |
| 평균 중간값 비트 | 1759 | 253 | 0.144 |
| 최대 출력 비트 | 128 | 127 | ~1.0 |

### NIST Level 3 (p ~ 2^381, 노름 비트 크기 = 193, 5회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **3279** | **387** | **0.118** |
| 평균 중간값 비트 | 3267 | 385 | 0.118 |
| 최대 출력 비트 | 194 | 193 | ~1.0 |

### NIST Level 5 (p ~ 2^509, 노름 비트 크기 = 254, 3회 시행)

| 지표 | HNF | MLLL | 비율 (MLLL/HNF) |
|------|-----|------|-----------------|
| **최대 중간값 비트** | **4312** | **510** | **0.118** |
| 평균 중간값 비트 | 3971 | 508 | 0.128 |
| 최대 출력 비트 | 255 | 254 | ~1.0 |

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

2026-03-30 (초기 벤치마크), 2026-04-02 (전체 테스트 통과 확인)
