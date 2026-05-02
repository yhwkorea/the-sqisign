# SQIsign Fixed-Precision을 위한 MLLL 구현

참조 논문:
- **주 논문**: "SQIsign with Fixed-Precision Integer Arithmetic" (김원, 이정환, 김현학, 이창민 - 고려대학교, 2025) — `klkl25.pdf`
- **MLLL 알고리즘**: Pohst (1987), Matthews pseudo-code, Cohen Algorithm 2.6.7

## 동기: 왜 HNF를 MLLL로 교체하는가?

논문 Section 1 "Challenges"에서 **IdealMultiplication**을 fixed-precision의 병목으로 지목:

> "HNF 중간값의 크기는 4×4 minor들의 gcd로 제어할 수 있으나, (16 choose 4)개의 minor 계산은 비효율적. 원래 Round-2 알고리즘은 4개의 minor만 계산하므로 gcd가 커지고, 과대한 modulus로 HNF를 수행해 numerator가 급격히 팽창"

논문은 **Algorithm 9** (Modified IdealMultiplication)으로 modulus를 `m = nrd(r1*I1)^2 * nrd(r2*I2)^2`로 줄임 (Lemma 11). 우리 접근: **HNF를 완전히 MLLL로 대체**, modular arithmetic과 m^2 중간값 폭발을 원천 차단.

### 이론적 Bound 비교

| 접근 | IdealMul 중간값 bound | 출처 |
|------|---------------------|------|
| 원본 HNF (Alg. 2) | `2^10 * p^10 * r1^32 * r2^8 * nrd(I1)^16 * nrd(I2)^4` | Lemma 2 |
| 수정 HNF (Alg. 9) | `r1^8 * r2^8 * nrd(I1)^4 * nrd(I2)^4` | Lemma 12 |
| **MLLL (본 구현)** | `max ||a_i||^2` (a_i = generator) | Lemma 1 + Pohst |

MLLL은 modulus를 사용하지 않으므로 m^2 blowup이 없음. 중간값은 generator의 norm에 비례하며, 이는 수정 HNF bound보다 **이차적으로 작음**.

### nrd(I)의 역할

논문의 bound는 근본적으로 입력 ideal의 **reduced norm** `nrd(I)`로 표현됩니다.

**Lemma 12** (Modified IdealMul bound): 중간값 ≤ `r1^8 * r2^8 * nrd(I1)^4 * nrd(I2)^4`

- `nrd(I)` = ideal I의 reduced norm = gcd{nrd(α) | α ∈ I}
- `r_i` = basis 원소 분모의 lcm (Lemma 5: `r_α ≤ 2*nrd(connecting_ideal)`)

**nrd(I)가 중간값 크기를 결정하는 방식**: modulus m = nrd(r1·I1)^2·nrd(r2·I2)^2. HNF over Z_m의 중간값 ≤ m^2. 따라서 중간값은 **nrd(I)의 4제곱**에 비례.

**파이프라인 합성** (nrd(I)가 전파되는 과정):

| 서브루틴 | 관련 Ideal | nrd bound | 유도 |
|---------|-----------|-----------|------|
| KeyGen Line 2 | I_sk (비밀) | nrd(I_sk) ≤ D_mix ≈ 2^(8λ) | RandomIdealGivenNorm (Alg. 6) |
| KeyGen Line 4 | I_sk (재샘플링) | nrd(I_sk) ≤ p (가설) | RandomEquivalentPrimeIdeal |
| KeyGen Line 5 | J_t · I_sk | nrd(J_t) < (log p)^2 | IdealToIsogeny → SuitableIdeals |
| Sign Line 13 | [I_sk]* I'_chl | nrd ≤ nrd(I_sk)·nrd(I'_chl) ≤ p·2^f | Lemma 9, 10 |
| Sign Line 14 | I_sk · I_chl | nrd(I_chl) = nrd(I'_chl) ≤ 2^f | Lemma 12 적용 |
| Sign Line 24 | I_com,rsp ∩ I_aux | nrd ≤ q_rsp·D_mix · √p | Lemma 9 |

**Theorem 1** (Modified KeyGen): 모든 nrd(I) ≤ p이므로, Lemma 12에 의해 ≤ r^8·p^4. 전체 KeyGen bound 2^68·p^9은 RandomIdealGivenNorm의 `32·p·D_mix^2 ≈ 2^68·p^9`이 지배항.

**Theorem 2** (Modified Sign): 최악의 경우는 Sign Line 13에서 nrd(I_sk)^20·nrd(I'_chl)^4 ≤ p^20·2^(4f) ≤ p^28. 상수 포함 = **2^16·p^28**.

**MLLL의 경우**: m^2를 사용하지 않으므로, nrd(I)는 generator norm을 통해서만 중간값에 영향:

```
MLLL 중간값 ≤ max_s,t { nrd(r1·α_s · r2·β_t) }
             = max { r1^2·r2^2 · nrd(α_s) · nrd(β_t) }
```

LLL-reduced basis에서 nrd(α_s) ≤ 8·nrd(I) (Minkowski bound)이므로:

```
MLLL 중간값 ≤ 64 · r1^2 · r2^2 · nrd(I1) · nrd(I2)
```

이는 Modified HNF bound의 **제곱근** (Lemma 12: r^8·nrd^4 vs MLLL: r^4·nrd^2). L1 실측 GRAM vec 259 bits가 HNF 실측 2160 bits보다 훨씬 작은 이유.

### 실험적 검증 (2026-04-19 C1 감사 반영, 10k trials × 3 levels, `--mode=alg2`)

| Level | 입력 bits | GRAM vec max (Lemma 1) | GRAM Gram max (Lemma 3) | HNF int. max | Cohen MLLL max |
|-------|----------|-----------------------|--------------------------|--------------|----------------|
| L1 (p~2^253) | 127 | **259** | **518** | 2160 | 3543 |
| L3 (p~2^381) | 193 | **391** | **782** | 3293 | 5389 |
| L5 (p~2^509) | 254 | **513** | **1026** | 4330 | 7102 |

**읽는 법**:
- **GRAM vec**은 Lemma 1의 벡터 좌표 bound (~`nrd(α)·nrd(β)`의 제곱근 + p 보정).
- **GRAM Gram entry**는 Lemma 3의 Gram 원소 bound (~`2·max‖a‖²`)에 정확히 포화. Fixed-precision typedef 결정의 근거.
- **Cohen MLLL max**는 비교용(d/Λ integral GSO 경로). Lemma 3 밖이라 ~6× 폭주.

**C1 픽스 영향** (감사 `mlll.c` `track_intermediate` 5곳 누락):
- Cohen 경로: L1 3028 → **3543** bits (+515). 추적 구멍이 그 크기만큼 실재.
- GRAM 경로: L1 516 → **518** (+2), L3 780 → 782 (+2), L5 1025 → 1026 (+1) — 거의 다 추적돼 있어 변동 미미.
- 결론: **B\*_{1,3,5} typedef 결정엔 영향 없음**.

## 논문 알고리즘 → 코드 매핑

### Algorithm 9 (Modified IdealMultiplication) → `quat_lattice_mul_mlll()`

```
논문 (Algorithm 9)                     코드 (mlll.c:437-479)
──────────────────                     ─────────────────────
Line 1: r1 = lcm(r_alpha_i)           (caller에서 처리 — lattice denom)
Line 2: alpha_i *= r1, beta_j *= r2   (lattice basis 표현에 암묵적)
Line 3: m = nrd(r1*I1)^2*nrd(r2*I2)^2 ← 대체됨: modulus 불필요
Line 4: HNF(alpha_i*beta_j) over Z_m  ← 대체됨: quat_mlll(generators, 16)
Line 5: return M / (r1*r2)            res->denom = lat1->denom * lat2->denom
```

**핵심 차이**: HNF modulo m (m^2 중간값) 대신 16개 generator를 직접 MLLL에 투입, Lemma 1에 의해 벡터가 bounded.

```c
// mlll.c:454-468 — 16개 quaternion 곱 생성 (논문의 행렬 M과 동일)
for (int k = 0; k < 4; k++) {          // alpha_k (I1의 basis)
    for (int i = 0; i < 4; i++) {      // beta_i (I2의 basis)
        quat_alg_coord_mul(&elem_res, &elem1, &elem2, alg);
        // generators[4*k+i] = alpha_k * beta_i   (논문: M의 l=4k+i 열)
    }
}
quat_mlll(&(res->basis), &rank, generators, 16, alg);  // HNF → MLLL 대체
```

### Lemma 1 (Generator Entry Bound) → 벤치마크로 검증

논문 Lemma 1:
```
|M_{1,l}|, |M_{2,l}| ≤ sqrt(nrd(alpha) * nrd(beta))                 ← 좌표 0,1
|M_{3,l}|, |M_{4,l}| ≤ (1/sqrt(p)) * sqrt(nrd(alpha) * nrd(beta))   ← 좌표 2,3
```

nrd(α·β) = M₁² + M₂² + p(M₃² + M₄²) ≤ nrd(α)·nrd(β) 로부터 각 항의 norm-form 기여분을 분리한 bound.

#### 이론 vs 실측 (10k trials, 2026-04-19 C1 픽스 후)

| Level | nrd ≈ | p ≈ | 좌표 0,1 이론 bits | 좌표 2,3 이론 bits | GRAM vec 실측 max |
|-------|-------|-----|---------------------|---------------------|-------------------|
| L1 | 2^127 | 2^253 | 127 | ≈ 0 | **259** |
| L3 | 2^193 | 2^381 | 193 | ≈ 2 | **391** |
| L5 | 2^254 | 2^509 | 254 | ≈ 0 | **513** |

> 좌표 2,3 이론 bits = `(nrd(α) + nrd(β) − log₂ p) / 2` = L1: `(127+127−253)/2 ≈ 0.5`, L3: `(193+193−381)/2 ≈ 2.5`, L5: `(254+254−509)/2 ≈ −0.5` → 사실상 0 (거의 0벡터).

#### 실측이 이론 bound의 ~2× 커지는 이유

Lemma 1은 **단일 기약 원소** α·β의 좌표 bound. 실측은 MLLL intermediate 벡터 `b[]`의 좌표 최대치로, 다음 요인이 겹쳐짐:

1. **Denom scale**: `lat1.denom = r_1`이 basis 좌표에 곱해진 상태로 `quat_alg_coord_mul`에 들어감. r·α·β의 좌표는 r²만큼 커짐 (~`r_1·r_2` 추가 factor).
2. **16개 generator 조합**: 4×4 basis 쌍 중 non-min원소(LLL-reduced basis 내부 b_3, b_4 등)가 b_1보다 큼. Minkowski 2nd/4th successive minima까지 들어가면 `√(nrd) → √(γ_4·nrd) ≈ 2√(nrd)` 수준 팽창.
3. **Size-reduce 전 일시 팽창**: MLLL Pohst 루프에서 dependent 벡터 제거 전 `b[m] -= q·b[l]` 직후 최대치가 순간적으로 2× 부근까지 오를 수 있음 (tracker는 `ibz_mul` 출력을 찍으므로 이 순간 포착).

결과적으로 실측 **vec bits ≈ 2 · log₂(nrd) + O(log p^(1/2))** 관계. 차수는 Lemma 1 bound와 동일(= log nrd + constant).

#### 핵심 관찰: Lemma 3이 fixed-precision 결정 요인

Lemma 3: `|Gram entry <b_i, b_j>| ≤ 2·max‖a‖²`. 여기서 `‖·‖²`는 norm form `a₀² + a₁² + p(a₂² + a₃²)`. GRAM 경로에서 LLL-reduced basis는 `a_2, a_3 ≈ 0` 성분이 지배하므로 사실상 `‖·‖² ≈ a₀² + a₁² ≈ max|a_k|² · 4`. 결과적으로:

```
Gram_bits  ≤  2·vec_bits + 3    (Cauchy-Schwarz: <b_i,b_j>² ≤ ‖b_i‖²·‖b_j‖²)
```

| 대상 | 이론 | 실측 (L1) | 실측 (L3) | 실측 (L5) |
|------|------|-----------|-----------|-----------|
| vec max | `≤ 2·log₂(nrd)+O(1)` | **259** | **391** | **513** |
| Gram max | `≤ 2·vec + 3` | **518** = 2·259+0 | **782** = 2·391+0 | **1026** = 2·513+0 |

실측 Gram bits가 정확히 `2·vec bits` — Cauchy-Schwarz bound의 **상한 타이트**. p 인자가 들어간 항(p·a₂·b₂)이 기여하지 않는 이유는 reduced basis에서 a₂, a₃ 좌표가 1~2 bits 수준이라 a₀·b₀ + a₁·b₁ 항이 완전 지배하기 때문. Fixed-precision typedef 폭(Gram 9/13/17 × u64)은 이 2× 관계를 기준으로 잡혀 있으며, 10k trials 전체에서 이 bound 위로 벗어난 샘플이 **0건**이었음.

**Why this matters**: Cohen integral GSO 경로는 `d[k]/Λ[i][j]` 값이 Lemma 3 밖이라 L1에서 3543 bits, L5에서 7102 bits까지 폭주. GRAM 경로(L²-style Gram + dpe)만 Lemma 3에 의해 bounded — fixed-precision 전환의 **유일한 안전 경로**.

## MLLL 내부 알고리즘

### `quat_mlll()` 핵심 함수 매핑 (mlll.c:212-434)

| 구성요소 | 논문/참조 | 함수 | 줄 |
|---------|----------|------|-----|
| **내적** | `<a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3` (quaternion norm form) | `ibz_vec_4_dot_quat()` | 32-45 |
| **Gram 행렬** | `G[i][j] = <b[i], b[j]>` — 정확한 정수 | `gram_compute_row()`, `gram_size_reduce()`, `gram_swap()` | 68-111 |
| **d[k] / lam[i][j]** | Cohen Alg. 2.6.7: d[0]=1, d[k+1]=det(Gram(b[0..k])) | `compute_lam_d_row()`, `recompute_all_lam_d()` | 160-208 |
| **최근접 정수 나눗셈** | REDI 몫: q = round(lam[m][l] / d[l+1]) | `rounded_div()` | 118-139 |
| **정확 나눗셈** | Cohen 점화식: num = (d[s+1]*num - lam*lam) / d[s] | `exact_div()` | 144-152 |

### MLLL 단계 (Pohst 1987 / Matthews)

| 단계 | Matthews | 코드 | 설명 |
|------|----------|------|------|
| **LOAD** | "generator를 b[beta]에 적재" | 245-268 | generator 복사 → Gram 행 계산 → d/lam 재계산 → d[beta]=0 확인 |
| **REDI** | "b[m]을 b[l]에 대해 size-reduce" (Cohen 2.6.7) | 283-309 | \|2*lam\| > d이면: q=round(lam/d), b[m]-=q*b[l], Gram+lam 업데이트 |
| **Lovász 검사** | "4*(d[m-1]*d[m+1]+lam^2) < 3*d[m]^2 → swap" | 327-339 | 표준 LLL 조건 (delta=3/4) |
| **D_k=0 감지** | Matthews: "if D_m=0" (Pohst 1987 §3) | 320-325 | 정확한 정수 테스트: `ibz_is_zero(&d[m+1])` |
| **강제 swap (dep.)** | Matthews flagg=1 / Swap2 | 352-360 | `do_swap_dep`: tau 무시 (`m=max(m-1,1)`), 위치 1까지 cascade |
| **Lovász swap** | Matthews Swap1 | 363-371 | `do_swap_lov`: tau 존중 (`m=max(m-1,tau)`) |
| **영벡터 제거** | "b[m]=0 → shift, beta--" | 374-401 | alpha>=g 후: m=1부터 reduction 재실행 |

### D_k=0 처리 상세

```
Matthews (1987)                        우리 구현
──────────────                         ──────────
Case 1 (λ=0, flagg=1):               ┐
  d[m-1]=0, lambda들 zero화            │ 두 경우 통합:
  m = max(m-1, tau)                   │   물리적 swap b[m] ↔ b[m-1]
                                      │   Gram에서 d/lam 전체 재계산
Case 2 (λ≠0):                        │   m = max(m-1, 1)  ← tau 무시!
  Swap2 (특수 lam 업데이트)             │
  Swap1 (물리적 swap)                  │ 전체 재계산으로 복잡한 Swap2
  m = max(m-1, tau)                   ┘ 증분 공식 불필요
```

**tau 무시 이유**: dependent vector가 먼 위치의 벡터에 의존할 수 있음 (예: b[5]=(3,0,0,0)이 b[1]=(6,0,0,0)에 의존). `max(m-1,tau)`로는 cascade가 tau에서 멈춰 무한 루프. `max(m-1,1)`로 위치 1까지 내려가야 REDI로 제거 가능.

## 논문 Table 2와의 비교

논문의 전체 파이프라인 bound:

| | 원본 (Prop. 1,2) | 수정 (Thm. 1,2) |
|---|---|---|
| **KeyGen** | 2^50 * p^14 * (log p)^96 | **2^68 * p^9** |
| **Sign** | 2^50 * p^95 | **2^16 * p^28** |

NIST-I (p ≈ 2^253):
- 수정 KeyGen bound: 2^68 * 2^(253×9) = 2^2345 bits
- 수정 Sign bound: 2^16 * 2^(253×28) = 2^7100 bits
- 논문 실험 최대 (Table 2): 7,026 bits → 7,040 necessary

우리 MLLL 중간값 (389 bits, L1)은 **단일 IdealMultiplication 호출**의 측정값. 논문의 7,026 bits는 Sign 전체 체인 (다수의 IdealMultiplication 합성)을 포함.

## 테스트 → 논문 검증

| 테스트 | 검증 내용 | 논문 관련성 |
|--------|---------|------------|
| `dependent_generators` | Z-dependent 입력 처리 | D_k=0 case (Pohst §3) |
| `is_lll_reduced` | delta-LLL 출력 (delta=3/4) | 논문은 LLL-reduced basis 가정 (Lemma 2) |
| `vs_hnf_lattice_mul` | MLLL lattice = HNF lattice | Alg. 9 정확성 (같은 ideal, 다른 basis) |
| `compact_ideal_multiplication` | IdealMul 전체 파이프라인 | Algorithm 9 end-to-end |
| `realistic_scale` | 127-bit prime norm | NIST Level 1 parameter regime |
| `tau_path` | generator 적재 순서 | MLLL load/reduce/load cycle (Pohst §2) |

## 현재 상태

- **테스트**: MLLL 단위 6/6 통과, SQIsign 전체 35/35 통과
- **SQIsign 파이프라인**: lvl1/3/5 각 10회 반복 KeyGen+Sign+Verify 전부 통과
- **HNF 의존성**: MLLL 경로에서 완전 제거

### SQIsign 전체 파이프라인 타이밍 (10회 반복, KeyGen+Sign+Verify)

| 레벨 | 소수 | 10회 총 시간 | 1회당 | HNF 1회 (ctest) |
|------|------|-------------|-------|-----------------|
| L1 | p ≈ 2^253 | 0.79 s | **79 ms** | 81 ms |
| L3 | p ≈ 2^381 | 1.98 s | **198 ms** | 206 ms |
| L5 | p ≈ 2^509 | 3.37 s | **337 ms** | 344 ms |

파이프라인 수준에서는 HNF와 비슷 — isogeny 연산이 지배적이므로 MLLL 오버헤드는 미미함.

### Alg 2 CompactIdealMultiplication 벤치 (10k trials, 2026-04-19 C1 픽스 후)

| 레벨 | 입력 bits | HNF avg/max | GRAM avg/max | Cohen avg/max | 출력 bits |
|------|----------|-------------|--------------|---------------|-----------|
| L1 | 127 | — / 2160 | — / 518 | — / 3543 | — |
| L3 | 193 | 3096 / 3293 | 769 / 782 | 5244 / 5389 | 192-194 |
| L5 | 254 | 4075 / 4330 | 1014 / 1026 | 6883 / 7102 | 253-255 |

### Alg 3 RandomIdealGivenPrimeNorm 벤치 (10k trials)

| 레벨 | 입력 bits | HNF avg/max | GRAM vec avg/max (Lemma 1) | GRAM Gram avg/max (Lemma 3) | HNF ms / GRAM ms (GRAM/HNF) |
|------|----------|-------------|-----------------------------|-------------------------------|------------------------------|
| L3 | 193 | 772 / 778 | 194 / 195 | 767 / 771 | 191.05 / 119.13 (0.62) |
| L5 | 254 | 1016 / 1022 | 254 / 256 | 1012 / 1015 | 210.23 / 128.14 (0.61) |

**Fixed-precision 요구사항** (GRAM 경로):
- L1: vec ≤ 259 → **5×u64 (320 bits)**, Gram ≤ 518 → **9×u64 (576 bits)**
- L3: vec ≤ 391 → **7×u64 (448 bits)**, Gram ≤ 782 → **13×u64 (832 bits)**
- L5: vec ≤ 513 → **9×u64 (576 bits)**, Gram ≤ 1026 → **17×u64 (1088 bits)**

Cohen 경로는 d/Λ integral GSO 추적 폭이 Lemma 3 밖이라 ~6× 더 큼(L5 7102 bits). 프로덕션은 **GRAM 기본**, Cohen은 비교/교차검증 전용.

### 시간 성능 (Alg 2, 10k trials 총합 ms)

| 레벨 | HNF | Cohen MLLL | GRAM | GRAM/Cohen | GRAM/HNF |
|------|-----|------------|------|------------|----------|
| L3 | 2098 | 127270 | 7043 | 5.5% | 3.36x |
| L5 | 2332 | 238530 | 9454 | 4.0% | 4.05x |

- Cohen 경로는 recompute-from-scratch O(β²) 때문에 HNF 대비 ~60-100× 느림. GRAM 경로는 ~3-4× 느림 수준으로 실용적.
- 레벨↑ GRAM 상대 비중↑ (L3 5.5% → L5 4.0%) — bitsize 누적 효과가 Cohen 쪽에서 더 가파르게 커짐.

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

## 구현 이력

| 날짜 | 변경 |
|------|------|
| 2026-03-30 | 초기 MLLL (`ibq_t` exact rational GSO) |
| 2026-04-02 | Bitsize tracker, 벤치마크 |
| 2026-04-07 | tau swap 버그 수정, tracker 확장 |
| 2026-04-08 | Integral GSO 실험 (`ibq_t` → `ibz_t`) |
| 2026-04-09 | Float GSO (`dpe_t`), Gram matrix, HNF 후처리 우회 |
| 2026-04-09 | **정확한 정수 d[]/lam[][] 재작성. D_k=0 정상 처리. HNF 제거.** |
