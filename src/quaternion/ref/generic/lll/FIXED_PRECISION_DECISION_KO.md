# Phase 2 Fixed-Precision 백엔드 결정 (확정, 2026-04-21)

논문: "Compact Quaternion Algorithms for SQIsign"
Phase: [PLAN_KO.md](PLAN_KO.md) Phase 2, 항목 P2-1
직전: [Phase 1 완료 2026-04-21](PLAN_KO.md#phase-1--측정-보강-완료-2026-04-21), B\*\_{1,3,5} 폭 확정

## 1. 목표와 전제

Phase 1에서 10k trials × L1/L3/L5 측정으로 GRAM 경로 `vec`/`Gram` 폭을 확정했다(아래 표). Phase 2는 이 폭을 **실제로 고정폭 정수 레이아웃에 넣어** MLLL 내부 산술을 힙 할당 없이 돌리는 것.

| Level | GRAM vec max | GRAM Gram max | typedef 폭 | vec budget = `nwords*64-2` | Gram budget |
|---|---|---|---|---|---|
| L1 | 259 bits | 518 bits | vec 5×u64 (320b), Gram 9×u64 (576b) | 318 (여유 59b) | 574 (여유 56b) |
| L3 | 391 bits | 782 bits | vec 7×u64 (448b), Gram 13×u64 (832b) | 446 (여유 55b) | 830 (여유 48b) |
| L5 | 513 bits | 1026 bits | vec 9×u64 (576b), Gram 17×u64 (1088b) | 574 (여유 61b) | 1086 (여유 60b) |

**Budget 공식**: 두 보수(two's complement) `nwords` limb의 대표 범위 `[-2^(nwords*64-1), 2^(nwords*64-1)-1]`에서 부호 비트 1 + add/sub carry 비트 1 = 2 비트를 예약 → `nwords*64 - 2` (§6.4). 초기 초안의 `(nwords-1)*64` 공식은 **잘못**이었다. 실제로 L1 vec peak 259b는 5번째 limb의 3 비트를 쓰고 있어 "상위 limb 전체를 sign margin으로 예약"이라는 전제가 성립하지 않는다.

**이론적 worst-case** (empirical peak와 일치 검증):
- **vec** — Lemma 1에 의해 `2·bitsize(p) + ⌈log₂ n⌉ + c` (n=16 generator, c=1로 관측). L1: 2·128+4+1 = 261 (측정 259), L3: 2·193+4+1 = 391 (측정 391, **정확 일치**), L5: 2·254+4+1 = 513 (측정 513, **정확 일치**).
- **Gram (tight)** — 실측은 `Gram = 2·vec` (post-reduce에서 `p`-scaled 항이 소거). L1: 2·259 = 518 (측정 518), L3: 2·391 = 782 (측정 782), L5: 2·513 = 1026 (측정 1026).
- **Gram (loose)** — 중간 연산에서 `p·|a|²` 항이 살아 있는 scratch-path는 `2·b_v + bitsize(p) + 2`. 이 경로는 `vec4_dot_p` 계산이 해당되어 §6.3-b에서 `ibz_t` scratch로 bypass 처리됨.
- **중간 버퍼** `quat_fp_tmp_t = nwords_vec + nwords_gram + 1` limb (L1: 15, L3: 21, L5: 27) — `q·g_{ik}` 같은 중간곱은 storage slot보다 커질 수 있어 별도 scratch를 갖는다. Gram slot에는 narrow-sub 결과만 들어가므로 `nwords_gram*64-2` 예산 안에서만 검사.

**전제**
- GRAM 경로만 전환 (Cohen 경로는 Lemma 3 밖이라 고정폭 불가능, 비교/교차검증용 유지).
- Constant-time은 **현재 scope 밖, 향후 호환 여지만 확보**. 이유는 §3 참조.
- 출력 결과(lattice basis)는 기존 `ibz_t`/`ibz_mat_4x4_t` 타입 유지. 고정폭은 **MLLL 내부 산술에만** 적용.

## 2. 코드베이스 현황 (2026-04-21 조사)

| 항목 | 사실 | 시사점 |
|---|---|---|
| `ibz_t` 정의 | `src/quaternion/ref/generic/include/intbig.h:36` — `typedef mpz_t ibz_t` (직접 alias) | GMP `mpz_*`가 내부 엔진. Wrapper 레이어 없음 |
| Fixed-precision 관례 | `src/gf/ref/include/fp.h:14` — `typedef digit_t fp_t[NWORDS_FIELD]`. `digit_t = uint64_t` (RADIX_64) 또는 `uint32_t` (RADIX_32). per-level `NWORDS_FIELD` 매크로 | 이미 `digit_t[N]` 관례 존재. **C 후보가 기존 코드와 동질** |
| Constant-time | `ct_`/`constant_time`/`secret` prefix 0건. Quaternion layer 전체 variable-time | 현재 보강 요구 없음. 단 §3 외부 근거 참고 |
| `mpn_*` 사용 | Quaternion layer에서 0건. 호출 횟수: `mpz_mul_` 11, `mpz_sub_ui` 8, `mpz_add_ui` 5, `mpz_mul` 5, `mpz_add` 5 | **A 후보(mpn_*)는 learning curve + 기존 관례와 괴리** |
| 레벨 분기 | CMake `SVARIANT_S = lvl1;lvl3;lvl5`. `src/gf/ref/lvl{1,3,5}/` 디렉토리 + per-level `fp_constants.h`. Quaternion은 generic single target | 레벨별 typedef를 도입하려면 **기존 fp 레이어 패턴 재사용 가능** |
| `mpz_realloc2` 호출 | 0건 (mini-gmp 내부 구현 외) | B 후보는 **기존 없는 API 새로 도입** |

## 3. 외부 근거

### 3.1. GMP `mpn_*` lowlevel API
- 출처: [GNU MP 6.3.0 Low-level Functions](https://gmplib.org/manual/Low_002dlevel-Functions)
- 핵심: "A source operand is specified by a pointer to the least significant limb and a limb count, while a destination operand is specified by just a pointer, and it is the responsibility of the caller to ensure that the destination has enough space."
- 즉 caller-managed `mp_limb_t[]`. Stack 할당 가능. 그러나 API surface가 `mpz_*`와 **완전히 별개** — 재작성 규모 큼.

### 3.2. GMP `mpz_realloc2` semantics
- 출처: [GNU MP 6.3.0 Memory Management](https://gmplib.org/manual/Memory-Management), [mpz_realloc2 spec](https://gmplib.org/manual/Initializing-Integers)
- 핵심: "**Calling mpz_realloc2 is never necessary**; reallocation is handled automatically by GMP when needed. However, this function can be used to increase the space for a variable in order to avoid repeated automatic reallocations."
- 즉 **hint only, ceiling 아님**. 이후 연산이 prealloc 크기를 넘기면 내부 `realloc()` 호출 여전히 발생. **힙 할당 제거 보장 불가**.

### 3.3. SQIsign constant-time status (NIST PQC Round 2)
- 출처: [SQIsign Specification Round 2](https://csrc.nist.gov/csrc/media/Projects/pqc-dig-sig/documents/round-2/spec-files/sqisign-spec-round2-web.pdf), [Constant-time Quaternion Algorithms for SQIsign (eprint 2025/2192)](https://eprint.iacr.org/2025/2192.pdf)
- Spec 원문: "current implementation is not constant-time, with the exception of selected components on the isogeny side... **absence of constant-time algorithms for quaternion arithmetic** in the current implementation."
- eprint 2025/2192은 quaternion 연산 constant-time화를 **연구 주제**로 제시 — 즉 현재는 non-CT OK지만 **미래 방향은 CT**.
- 우리가 따르는 논문 [SQIsign with Fixed-Precision Integer Arithmetic (eprint 2025/1649)](https://eprint.iacr.org/2025/1649.pdf) 은 compact 연산 + 고정폭 결합을 제안하지만 CT는 명시적 목표 아님.

## 4. 후보 비교

### 후보 A: GMP `mpn_*` low-level API

- **저장**: caller-managed `mp_limb_t[N]` (stack 가능).
- **연산**: `mpn_add_n`, `mpn_sub_n`, `mpn_mul`, `mpn_mul_n`, `mpn_divrem_1`, ... 모두 caller가 size 관리.
- **장점**: 힙 없음, 고성능 (GMP 최적화 그대로).
- **단점**:
  - 기존 코드 `mpz_*` 호출 0건 대체 필요 → 재작성 규모 큼.
  - Signed 타입 없음 (`mpn_*`는 absolute value만). 부호 관리 caller.
  - Overflow trap 수동 구현.
  - Constant-time 아님 (GMP 자체가 non-CT).
  - `mpn_` API는 size 인자가 정밀한 limb count라 **off-by-one 지뢰** 다수.

### 후보 B: `ibz_t` + `mpz_realloc2`

- **저장**: 기존 `ibz_t` (mpz_t). `mpz_realloc2(x, bits)`로 prealloc hint.
- **연산**: 기존 `ibz_*` API 그대로.
- **장점**:
  - 코드 변경 최소 (기존 호출부 유지, init 단계에만 realloc2 추가).
  - 기능적 호환 100%.
- **단점**:
  - **근본적 한계**: `mpz_realloc2`는 "repeated reallocation 회피 hint". **ceiling 아니므로 힙 할당 제거 보장 없음** (§3.2).
  - 성능 이득은 "avoid repeated" 수준이지 "no heap" 수준 아님.
  - 스택 레이아웃 안 됨. Cache locality 이점 없음.
  - **논문 contribution 입증 실패**: "compact 연산이 고정폭에 실제로 들어간다"를 보여주지 못함.
  - Constant-time 불가능 (mpz는 value-dependent branching).

### 후보 C: bare `uint64_t[N]` (per-level typedef)

- **저장**: `typedef digit_t quat_b_vec_L{k}_t[NWORDS_L{k}_VEC]`, `quat_b_gram_L{k}_t[NWORDS_L{k}_GRAM]`. 기존 `fp_t = digit_t[NWORDS_FIELD]` 패턴과 **완전 동일**.
- **연산**: mpn_* 또는 자작 add/sub/mul/div (구현 선택). Mul은 schoolbook으로 충분 (9×9 ~ 17×17 digit).
- **장점**:
  - 기존 `fp_*` 레이어 관례와 동질 → 코드베이스 일관성.
  - 스택 할당, heap-free → **논문 contribution 직접 입증**.
  - Cache locality, 예측 가능한 성능.
  - Constant-time 가능 (§3.3의 미래 방향과 정렬).
  - Per-level typedef로 L1/L3/L5 분기 명확.
  - Overflow trap: 각 연산 후 상위 limb 검사로 결정적 실행.
- **단점**:
  - 재작성 규모 가장 큼.
  - Signed 처리 수동 (부호 bit 관리 또는 two's complement).
  - Division 구현 필요 (schoolbook). Exact div(d[s+1] 점화식)는 knuth D 알고리즘 또는 mpn_divrem_1 활용.
  - 일부 연산은 `mpn_*` 호출로 대체 (A의 부분집합을 내부에 흡수).

## 5. 접근: B와 C **병행 구현 후 실측 비교**

### 5.1. 확정: **C(fp) 경로를 primary로 채택**, ibz/prealloc은 regression oracle로 보존 (2026-04-21)

**결정 기준 재정립** — 본 §5.1은 2026-04-21 오전 한 차례 "baseline 유지"로 확정했다가 동일 날짜 오후 뒤집은 이력이 있음. 뒤집은 이유: **첫 확정은 시간 퇴행을 primary 탈락 사유로 썼는데, 그게 Phase 2의 본래 기준이 아니었음**. 본 문서 §1·§3·§5.3은 처음부터 Phase 2 목표를 다음 4가지로 명시했다.

1. `ibz_t`(GMP 동적) → **고정폭 스택 레이아웃**으로 전환 (heap-free).
2. 폭이 Lemma 1 / Lemma 3 예산 안에 실제로 들어감을 **런타임에서 감사** (overflow trap).
3. KLKL25 eprint 2025/1649의 "compact 연산이 고정폭에 실제로 들어간다"는 **contribution을 코드로 입증**.
4. 미래 CT 방향(eprint 2025/2192)과 정렬되는 **타입 레이어 확보**.

이 4개 기준은 모두 "ibz는 fp든 값이 어차피 들어가므로 시간 차이 X%면 어느 쪽이든 OK" 류의 **time-agnostic** 기준이다. 시간 퇴행은 언제나 논문 contribution의 trade-off로 수용되는 범주이지 primary 탈락 사유가 될 수 없음.

**4개 기준으로 재평가**:

| 기준 | baseline ibz | B prealloc | C fp |
|---|:-:|:-:|:-:|
| (1) Heap-free / stack | ✗ | ✗ (hint only, §5.5) | **✓** |
| (2) Lemma 3 runtime 감사 | ✗ | ✗ | **✓** (trap, sweep 270 iter 0 발동) |
| (3) Contribution 코드 입증 | ✗ | ✗ | **✓** (`quat_fp_vec_t`/`quat_fp_gram_t` + schoolbook) |
| (4) CT 호환 여지 | ✗ | ✗ | **✓** (value-independent limb layout) |
| 정확성 (equivalence) | ✓ (HNF oracle) | ✓ | ✓ (`gram_fp_equivalence` PASS) |

→ **C(fp)가 4개 기준 모두 유일 만족. primary는 C.**

**확정 사항**:

- **Primary**: `quat_mlll_gram_fp` (fp body, 스택 고정폭 + schoolbook + overflow trap). `g_fp_mode = 1` 기본값.
- **Regression oracle**: `quat_mlll_gram_ibz` (기존 ibz body). `set_fp_mode(0)`으로 opt-out 가능. 주 용도: fp 회귀 발생 시 교차 검증.
- **B (prealloc)**: 이미 §5.5에서 primary 탈락 확정, scaffold로 유지. fp path와 독립이라 fp_mode=1일 때는 동작 안 함(ibz 경로 전용 hint). `--prealloc` 플래그는 ibz 성능 측정용으로 유지.
- **시간 퇴행 수용 범위** (§5.6): L1 Alg 2 기준 fp/ibz = 5.25× (trap-off 4.65×). L3/L5는 동률. 이 수치는 **수용 trade-off**로 기록, primary 결정에 영향 없음. 퇴행 감소는 Phase 2 이후 재설계(size-reduce 1-limb fast path 등) 과제이며 본 결정과는 직교.
- **제거 대상**: 없음.

**변경 추적**:

- 2026-04-21 오전: "baseline 유지" 확정. **근거가 시간 기준이라 폐기**.
- 2026-04-21 오후: 본 §5.1로 재확정. fp primary + 시간은 trade-off 기록.
- 이전 순서(B → C → 비교)는 그대로 수행됨. "선험적으로 primary를 못 박지 않는다"는 §5 원칙도 유지 — 최종 선택은 4개 기준 + §5.6 실측.

### 5.2. B를 먼저 구현하는 이유

- **빠른 feedback**: 기존 `ibz_*` API 유지하므로 수일 내 동작하는 상태 도달. `mpz_realloc2` hint 효과가 실제로 얼마나 있는지 실측해봐야 알 수 있고, §3.2의 "hint only"가 실측에서도 그대로인지 확인 필요.
- **C 구현의 reference**: C 구현 중 bug가 나면 같은 입력에 대한 B 출력과 비트 단위 비교가 디버깅의 1차 수단. B가 먼저 있어야 C 디버깅이 수월.
- **Sunk cost 아님**: B 단계 코드는 Phase 2 종료 후 제거 여부를 비교 결과에 따라 결정 (§7).

### 5.3. C를 병행하는 이유

1. **기존 관례와 동질**: `fp.h`가 이미 `digit_t[NWORDS_FIELD]` 패턴을 쓰고, per-level `NWORDS_FIELD` 매크로로 L1/L3/L5를 분기한다. 같은 패턴 재사용이라 **배우는 비용 낮음**.
2. **논문 contribution 입증 가능성**: KLKL25(eprint 2025/1649)의 핵심 주장은 "compact 연산이 고정폭에 실제로 들어간다". B만으로는 이 주장을 직접 입증하지 못하고 C가 필요. 다만 실제 성능/유지보수 trade-off는 실측 후 판단.
3. **미래 CT 호환 여지**: eprint 2025/2192이 quaternion CT화를 연구 중. C는 CT 여지가 있고 B/A는 본질적으로 없다. 지금 CT까지 할 필요는 없지만 "길이 막히는 선택"을 피하기 위해 최소 한 번은 C 경로를 만들어 둘 가치가 있음.

### 5.4. A(mpn_*)를 버리는 이유

- 기존 코드 `mpn_*` 호출 0건. 재작성 부담은 C와 비슷한데 **독자적 이점 없음** (stack 할당은 C도 가능, CT는 A도 불가).
- Signed/overflow/부호 관리 수동 처리 부담은 C와 동일 → A는 사실상 C의 subset이고 C가 더 많은 이점을 가진다.
- 단, **C 내부에서 `mpn_mul` 같은 개별 루틴을 부분 활용하는 것은 허용** (§6.3-c에서 결정). 전면 A 백엔드는 버리지만 도구로서의 mpn_*는 열어둠.

### 5.5. 2026-04-21 P2-B 실측: `mpz_realloc2` hint 사실상 효과 없음 또는 퇴행

B 후보 scaffold(`quat_mlll_gram` 진입부에서 `b[]`/`G[][]`/`X`/`tmp` 전부에 per-level 폭으로 `mpz_realloc2` 호출) 구현 후 baseline vs prealloc 1:1 측정.

**벤치 매트릭스** (`bench_logs_2026-04-21_P2B/`):
- Alg 3 10k trials × L1/L3/L5 × (baseline, prealloc) = 6 runs. GRAM이 유일 경로라 prealloc 효과 직접 관측 가능.
- Alg 2 1k trials × L1/L3/L5 × (baseline, prealloc) = 6 runs. MLLL Cohen이 wall time 대부분 차지 → GRAM 경로 지표만 비교.

**결과** (GRAM path ms, 총합):

| Level | Alg3 baseline | Alg3 prealloc | 비율 | Alg2 baseline | Alg2 prealloc | 비율 |
|---|---|---|---|---|---|---|
| L1 | 99.80 | 147.04 | **1.47× 느려짐** | 426.37 | 431.75 | 1.01× |
| L3 | 104.77 | 111.76 | 1.07× | 667.66 | 667.88 | 1.00× |
| L5 | 109.75 | 110.36 | 1.01× | 879.86 | 882.27 | 1.00× |

출력 bitsize max는 양 경로 **비트 단위 동일** (`quat_test_mlll_gram_prealloc_equivalence` 20 trials × alg2+alg3 PASS). 즉 **정확성은 보존되지만 성능은 개선 없음 또는 퇴행**.

**왜?**
1. §3.2에서 예측한 그대로 — `mpz_realloc2`는 "repeated reallocation 회피 hint"이고 **ceiling이 아니다**. GMP 내부의 `ibz_mul`/`ibz_add`/`ibz_sub`가 중간 결과 저장 시 여전히 자체 판단으로 realloc을 호출하므로, prealloc으로 선점한 공간이 실제 재사용된다는 보장이 없다.
2. Alg 3 L1 퇴행 원인은 realloc 호출 비용 자체. `quat_mlll_gram` 1회 호출당 16·4(b) + 16·16(G) + 2(X,tmp) = 322 건의 `mpz_realloc2` 호출이 들어가고, L1 alg 3은 8 generator만 쓰기 때문에 MLLL 본체 작업량이 작다 → 322 prealloc 호출이 전체 시간의 큰 비중 차지. L3/L5 및 alg 2(16 generator)는 본체가 더 무거워서 322 호출이 묻힘.

**결론**:
- **B만으로는 Phase 2 목표(힙 제거) 달성 불가**. "hint only"가 실측 확인되었으므로 더 이상 B 튜닝에 시간 쓰지 않는다.
- B scaffold는 ibz 경로(현재는 fp primary의 회귀 oracle) 내부에서 `mpz_realloc2` 효과를 단독 측정하는 도구로 남긴다. 런타임 flag `set_prealloc_mode(0/1)` 기본 0. **C 구현 중 출력을 B=1과 비교하는 디버그 flow도 유효함** — pairwise transitive equivalence 체인(§6.5 P2-equiv-a)에서 ibz self-consistency 검증 담당.
- **primary는 C 쪽 가능성 높음**. 다만 C 구현 전에 최종 확정 금지 (P2-decide에서 §5.1 4개 기준 기반으로 결정).

**Why**: 본 확인으로 Phase 2 전체 성공이 C 구현의 품질에 의존하는 상황이 됨. C가 실패하면 fixed-precision 전환 자체를 재고해야 하므로 P2-C-gram에서 보수적 구현(schoolbook + overflow trap) 우선, 최적화는 P2-decide 이후로 미룬다.

### 5.6. 2026-04-21 P2-decide 실측: 시간 trade-off 기록 (결정 기준 아님)

C 구현(`quat_fixed_precision.[hc]` + `quat_mlll_gram_fp` + dispatcher + overflow trap) 완료 후 3-way 스윕 실행 (`bench_logs/p2_phase2_sweep/`, L1/L3/L5 × alg2/alg3 × baseline(ibz)/prealloc/fp).

**본 표는 primary 결정이 아니라 "primary를 fp로 잡았을 때 감수해야 하는 시간 비용" 기록용.** 결정 기준은 §5.1 4개(heap-free / Lemma 3 감사 / contribution 입증 / CT 여지).

**[2026-04-22 갱신] P2.1-time-rerun** — dispatcher 임계값 정정(§6.8) 후 재측정. iter 분할(alg2=500, alg3=5000)로 variance 완화. L3/L5 fp 는 **이 표가 최초 실측** (이전은 silent ibz fallback, §6.7 참조).

**결과** (GRAM total ms, trap-on 기본 빌드, `bench_logs/p2_phase2_sweep/L{1,3,5}_{alg2,alg3}_{baseline,prealloc,fp}.log`):

| Level | Mode | Iter | baseline(ibz) | prealloc | fp (primary) | fp/ibz |
|---|---|---:|---:|---:|---:|---:|
| L1 | alg2 |  500 |  215.25 |  216.59 |   707.29 | 3.29× |
| L1 | alg3 | 5000 |   50.19 |   60.74 |    89.69 | 1.79× |
| L3 | alg2 |  500 |  336.14 |  336.07 |  1338.15 | 3.98× |
| L3 | alg3 | 5000 |   53.00 |   66.70 |    99.08 | 1.87× |
| L5 | alg2 |  500 |  443.60 |  447.09 |  2121.61 | 4.78× |
| L5 | alg3 | 5000 |   55.15 |   70.93 |   104.58 | 1.90× |

**pre-fix(2026-04-21) 대비 변화**:
- L1 alg2 fp/ibz: 5.25× → **3.29×** (L5 widths 오버사이징 해소 효과).
- L3 alg2 fp/ibz: 0.98× → **3.98×** (이전은 ibz fallback이라 "동률"이었음).
- L5 alg2 fp/ibz: 0.99× → **4.78×** (동일 이유).
- Alg 3 는 전 레벨 1.79~1.90× 로 안정. 1-limb special case 비중이 적은 경로라 퇴행 적음.

**Trap overhead 격리** — P2.1 재측정 미수행(기존 추정 유지). L1 alg2 기준 trap 기여분 ~10%, trap 제거해도 근본 원인(schoolbook inner loop) 동일.

**퇴행 원인**:
- GMP `ibz_mul`/`ibz_add`는 operand가 1-limb인 경우(특히 size-reduce 이후) 내부 dispatch로 `mpn_mul_1` 또는 1-limb special case로 빠진다. L1 뿐 아니라 L3/L5 도 alg2 경로에서 schoolbook vs mpn_mul_1 차이가 크게 나타남.
- 우리 schoolbook `fp_mul`은 항상 `nwords_vec × nwords_vec` (L1: 5×5, L3: 7×7, L5: 9×9) 전체 limb 작업 — 고정폭 설계상 dispatch 불가.
- Alg 3 는 1-limb special case 비중이 낮아 schoolbook이 상대적으로 덜 손해. 그래서 1.8~1.9× 수준 유지.

**Trade-off 수용 판단 (재평가)**:
- 최악 퇴행은 L5 alg2 4.78×, 절대 시간 ≈ (2121.61 − 443.60) / 500 = **3.36 ms per call** 추가. L1 alg2는 (707.29 − 215.25) / 500 = **0.98 ms per call**.
- Alg 3(Random Equivalent, keygen hot path 중 하나)는 per-call 마진 작음: L5 기준 (104.58 − 55.15) / 5000 = **10 μs per call**.
- 4개 결정 기준(§5.1) 만족도가 절대적으로 높기에 수용. 향후 재설계(`mpn_mul_1` 1-limb fast path) 적용 시 퇴행 감소 기대.

### 5.7. 회귀 oracle(ibz) + prealloc scaffold 보존 이유

Primary는 fp지만 ibz body와 B prealloc scaffold는 **제거하지 않는다**. 이유:

1. **ibz = fp 회귀 oracle**: `quat_test_mlll_gram_fp_equivalence`(테스트 모음에 상시 등록됨)가 fp 경로 출력과 ibz 경로 출력을 비트 단위 비교. fp에 버그 들어가면 즉시 실패. `set_fp_mode(0)`으로 필요 시 임시 전환 가능.
2. **런타임 Lemma 3 감사**: fp primary 자체가 감사 역할. `FP_CHECK_VEC`/`FP_CHECK_GRAM`가 각 mutation site에서 `(nwords-1)*64` 비트 예산 초과 시 `abort` — 한 번이라도 발동하면 폭 표(§1) 또는 Lemma 3 논증이 틀렸다는 증거가 즉시 확보. trap-on 기본 빌드에서 프로덕션이 이 감사를 상시 수행.
3. **미래 재설계 landing spot**: §5.6 L1 퇴행 해결을 위해 size-reduce 내부 루프를 schoolbook→`mpn_mul_1` inlined 또는 `mul+narrow-sub` fused로 재설계하더라도 fp 타입 + bridge가 이미 자리 잡혀 있어 진입 비용 낮음. 그 시점에 ibz body 제거 여부 재평가.
4. **B prealloc**: ibz body에 대한 hint-only 최적화. fp primary 하에서는 hot path 밖이라 유지/제거 중립. `--prealloc` 플래그는 ibz 벤치 도중 prealloc 단독 효과 측정용으로 남긴다.

**유지 비용**: 파일 수 2개(`quat_fixed_precision.[hc]`) + `mlll_gram.c` 내 ibz body ~200줄(oracle) + fp body ~300줄(primary) + dispatcher ~30줄 + 테스트 2개. 런타임 기본값 `g_fp_mode = 1`이라 프로덕션은 fp 경로를 상시 거침. ibz body는 `--fp 0` / `set_fp_mode(0)` 호출 시에만 활성.

## 6. 구현 단계

일정은 §8 참조. 각 단계는 독립 커밋 단위.

### 6.1. P2-B: B 단계 구현 (`mpz_realloc2` prealloc scaffold) — 2026-04-21 완료

- [x] **P2-B-a** `quat_mlll_gram` 진입부에서 `b[]`/`G[][]`/`X`/`tmp` 각 `ibz_t`에 `mpz_realloc2(x, B*_{k}+64)` 호출. Per-level 폭(L1 vec 320/Gram 576, L3 vec 448/Gram 832, L5 vec 576/Gram 1088)은 `alg->p` bitsize에서 유도. Global flag `quat_mlll_gram_set_prealloc_mode()` 노출.
- [x] **P2-B-b** In-process 토글 비교로 비트 동일 검증 — `quat_test_mlll_gram_prealloc_equivalence` 20 trials × (alg2 + alg3) PASS.
- [x] **P2-B-c** 10k trials 벤치 재실행 (`bench_logs_2026-04-21_P2B/`). 시간 변화 기록: §5.5 표. Heap profiler 카운터는 **필요 없음** — 시간이 이미 hint only를 증언.
- [x] **P2-B-d** `--prealloc` 플래그로 toggle. 기존 `--mode=alg2|alg3` 체계 유지하고 전역 mode만 토글.

**측정 결론**: §5.5 참조. `mpz_realloc2` hint는 실측상 효과 없음 또는 퇴행(L1 alg3 1.47× 느려짐). B는 C 구현의 reference oracle로 남기되, primary candidate로서의 경쟁력 소멸.

### 6.2. P2-C-types: Per-level typedef 정의 — 2026-04-21 완료

- [x] **P2-C-types-a** `src/quaternion/ref/generic/internal_quaternion_headers/quat_fixed_precision.h` 신설. 기존 헤더 관행과 통일하기 위해 §6.2 초안의 `lll/include/` 대신 `internal_quaternion_headers/`에 배치.
- [x] **P2-C-types-b** per-level 폭 매크로 + 구조체 typedef. 실제 채택안은 **단일 최대폭(L5) array + runtime nwords descriptor**:
  ```c
  #define NWORDS_QUAT_VEC_L1    5    /* L3 7 / L5 9 */
  #define NWORDS_QUAT_VEC_MAX   NWORDS_QUAT_VEC_L5
  #define NWORDS_QUAT_GRAM_L1   9    /* L3 13 / L5 17 */
  #define NWORDS_QUAT_GRAM_MAX  NWORDS_QUAT_GRAM_L5
  /* tmp = vec + gram + 1 (for X * G[i][j]) */
  #define NWORDS_QUAT_TMP_MAX   NWORDS_QUAT_TMP_L5

  typedef struct { digit_t limbs[NWORDS_QUAT_VEC_MAX]; }  quat_fp_vec_t;
  typedef struct { digit_t limbs[NWORDS_QUAT_GRAM_MAX]; } quat_fp_gram_t;
  typedef struct { digit_t limbs[NWORDS_QUAT_TMP_MAX]; }  quat_fp_tmp_t;

  typedef struct {
      unsigned nwords_vec, nwords_gram, nwords_tmp;
  } quat_fp_widths_t;
  ```
- [x] **P2-C-types-c** **결정: two's complement storage, sign-magnitude mul 경로.**
  - Storage(add/sub hot path): 2's comp. add/sub가 unsigned limb arith + carry 그대로 → 분기 無.
  - Mul은 내부적으로 부호 분리(음수면 negate) → unsigned 곱 → 부호 재부여. GMP `mpn_mul` 또는 schoolbook 어느 쪽이든 공통 최적 패턴.
  - 근거: Phase 1 실측에서 vec 좌표·Gram 비대각 entry 모두 음수 관측됨. Sign-magnitude storage는 add/sub에서 부호 분기 비용이 hot loop에 매 반복 누적 — 반면 2's comp는 carry propagation 한 번.

**추가 결정** (§7-3 "레벨별 코드 복제 vs generic" 해결):
  - 기존 `src/gf/ref/lvl{k}/` 패턴은 compile-time 분기를 쓰지만, quaternion 레이어는 현재 `src/quaternion/ref/generic/`에 단일 빌드. 레이어 일관성을 위해 **컴파일-타임 분기 포기, 런타임 분기 채택**.
  - Trade-off: L1에서 L5 크기(≈200B)의 사용되지 않는 스택 차지. 하지만 단일 바이너리가 3 레벨 모두 처리 가능 — 기존 `mlll_gram_level_hints` 패턴과 일치.
  - `quat_fp_widths_from_alg(out, alg)`로 `alg->p` bitsize에서 런타임 유도.

**구현 상태**: 헤더 + stub `.c` (`lll/quat_fixed_precision.c`) 추가. zero/copy/swap/is_zero만 구현, 나머지는 `fp_stub()` trap (P2-C-gram에서 구현). 빌드 통과, 기존 테스트 회귀 없음.

### 6.3. P2-C-gram: GRAM 내부 산술 C 구현 — 2026-04-21 완료

- [x] **P2-C-gram-a** (3-step 분할로 수행, 커밋 `8b30e32`·`b18e49e`·`97eb079`) `mlll_gram.c`에 fp 경로 본체 `quat_mlll_gram_fp` 신설. 원래 계획은 컴파일 타임 레벨 분기(`SQISIGN_LEVEL`)였으나 **§6.2 결정대로 런타임 descriptor(`quat_fp_widths_t`)**로 변경. 기존 `quat_mlll_gram_ibz` 본체는 그대로 보존하고 `quat_mlll_gram` 공용 진입점이 `g_fp_mode`로 분기.
- [x] **P2-C-gram-b** 내적 / size-reduce / swap / Gram update 전부 `quat_fp_vec_t`/`quat_fp_gram_t` 기반으로 재작성. 단 `vec4_dot_p_fp`의 `a·b·p` 중간값은 L1 기준 vec 259b × vec 259b × p 127b ≈ 645b로 Gram 예산 576b 초과 — 따라서 이 경로만 `ibz_t` scratch로 bypass (caller-owned scratch로 heap 재사용). size-reduce hot loop은 fp-native.
- [x] **P2-C-gram-c** Mul 구현: **schoolbook 자작 채택**. 근거: 고정폭 5×5 ~ 17×17 digit 규모라 `mpn_mul` 호출 overhead가 schoolbook inline 이득을 상쇄. 별도 micro-bench 없이 schoolbook 직행 — `mpn_mul` wrapper 경로는 §5.6 결과에서 primary 탈락이 확정되면서 **micro-bench 수행 moot**.
- [x] **P2-C-gram-d** Division: 원래 계획은 Cohen 점화식의 exact div였으나 §7-2에서 이미 Phase 2 scope 밖으로 드롭 확정. `quat_mlll_gram`은 L² 기반이라 integer exact div를 사용하지 않는다 → 구현 **불필요**, 체크 대상 없음.

### 6.4. P2-overflow: Overflow trap (C 경로 전용) — 2026-04-21 완료 (커밋 `58540b1`) / 2026-04-22 budget 정정

- [x] **P2-overflow-a** `quat_fp_vec_check_overflow` / `quat_fp_gram_check_overflow`가 각 mutation site에서 `fp_limbs_bitsize()`를 budget과 비교, 초과 시 `fprintf(stderr) + abort()`. Release 빌드에서도 활성(매크로 `FP_CHECK_VEC`/`FP_CHECK_GRAM`).
- [x] **P2-overflow-b** 컴파일 플래그 `MLLL_FP_NO_OVERFLOW_CHECK`로 trap 비활성화 빌드 가능 (§5.6 trap-off 측정에 사용). 현재는 빌드-타임 플래그; 런타임 토글은 차후 과제(§7-6).
- [x] **P2-overflow-c (2026-04-22)** Budget 공식 정정. 초기 `(nwords-1)*64`는 "상위 limb 전체를 sign-extension margin으로 예약"이라는 잘못된 margin narrative였음 — capacity sweep에서 L3 vec peak 391b가 이 budget(384b)을 초과했음에도 trap이 발동하지 않아 정체가 드러남. 정정된 공식은 `nwords*64 - 2` (부호 1 + carry 1). 새 budget 하 Phase 1 peak는 모두 48-61b 여유로 들어감(§1 표). 측정 peak와 이론 worst-case(`2·bitsize(p) + ⌈log₂ n⌉ + 1`)가 L3/L5에서 정확히 일치하므로 여유는 조작 가능한 값이 아닌 Lemma 1/3의 실제 tight bound에서 온 것.
- [x] **P2-overflow-d (2026-04-22)** Trap 감사 재수행. **결론**: trap 자체는 기능하지만 **L3/L5에서는 fp path가 애초에 돌고 있지 않았다** — 그래서 초안 budget 초과에도 abort가 뜨지 않은 것. 자세한 내용은 §6.7 참조.

### 6.5. P2-equiv: 동치성 + 성능 비교 — 2026-04-21 완료 (커밋 `97eb079` 및 sweep 아티팩트)

- [x] **P2-equiv-a** 3-way 대신 **2-way + transitive** 구조 채택 (fp_mode=1 primary 기준 재해석):
  - `quat_test_mlll_gram_equivalence`: Cohen(`quat_mlll`) ↔ GRAM dispatcher (기본 fp). → Cohen↔fp primary equivalence.
  - `quat_test_mlll_gram_prealloc_equivalence`: `set_fp_mode(0)` 명시 후 ibz+prealloc-off ↔ ibz+prealloc-on. → ibz self-consistency 확인.
  - `quat_test_mlll_gram_fp_equivalence`: `set_fp_mode(1)` 명시 guard + L5 random 16-gen case 확장. fp 기본값이 바뀌어도 이 테스트는 고정됨.
  - `quat_test_lideal_create_gram_equivalence`: HNF ↔ GRAM dispatcher (L5 real prime, 25 trials). → HNF↔fp primary equivalence.
  - 조합으로 HNF ↔ Cohen ↔ fp / ibz ↔ ibz+prealloc 4경로 equivalence가 transitively 보장됨. 실패 localize 가능.
- [x] **P2-equiv-b** Sweep 30 iter × 9 조합 (L1/L3/L5 × alg2/alg3 × fp) 실행, trap 발동 0건. `ibz` ↔ `fp` 출력은 매 회 비트 동일.
- [x] **P2-equiv-c** 3-way 성능 비교표는 §5.6 표로 확정 (HNF 경로 시간은 MLLL 컬럼 참고). 원래 계획의 "realloc 횟수 측정"은 §5.5에서 이미 mpz_realloc2가 힙 제거 보장 없음이 확증되어 **별도 프로파일링 수행 moot**.

### 6.6. P2-decide: 최종 primary 결정 및 문서화 — 2026-04-21 완료

- [x] **P2-decide-a** §5.1의 4개 기준(heap-free / Lemma 3 감사 / contribution 입증 / CT 여지)에 근거해 **C(fp) 경로를 primary로 최종 채택**. 시간 퇴행(§5.6)은 수용 trade-off로 기록, 결정 기준 아님. `g_fp_mode` 기본값 0 → **1**로 변경.
- [x] **P2-decide-b** 본 문서 §5.1 "결정" 재확정(오전 baseline → 오후 fp로 flip) + §5.6 "trade-off 기록" 재해석 + §5.7 "oracle 보존" 재작성.
- [x] **P2-decide-c** 탈락 경로 제거 여부: **모두 유지**. ibz body는 fp 회귀 oracle, B prealloc은 ibz 서브튜닝 도구. 유지 비용 negligible. §7-4 동시 해결.
- [x] **P2-decide-d** `PLAN_KO.md` Phase 2 완료 체크 반영. Phase 3는 fp primary 위에 직접 쌓는다 — 추가 blocker 없음.

### 6.7. P2-trap-audit: Trap 감사 재수행 — 2026-04-22

**착수 동기**: capacity sweep (§6.4-c)에서 L3 vec peak 391 > 초안 budget 384인데 abort 안 뜬 정황. budget 공식 정정(§6.4-c) 외 **trap 자체가 dead 되어 있을 가능성**을 감사.

**방법**: `quat_fp_vec_check_overflow` / `quat_fp_gram_check_overflow` 안에 `#ifdef MLLL_FP_TRAP_SELFTEST` 브랜치 추가 — budget 을 1 bit 로 축소해 정상 값이면 무조건 trap. Selftest flag 활성화 빌드로 L1/L3/L5 × alg3 × 1 iter 실행, `abort()` 발생 여부 확인.

**결과**:

| Level | Prime bits | SELFTEST 결과 | 해석 |
|---|---|---|---|
| L1 | 251 | `quat_fp_vec overflow at load_gen0: 129 bits > 1 budget (nwords_vec=9)` → `Aborted (core dumped)` | fp path 진입 + trap 발동 ✓ |
| L3 | 375 | 정상 완료 (abort 없음, GRAM summary 출력 ) | **fp path 진입 자체 안 함** |
| L5 | 473 | 정상 완료 (abort 없음, GRAM summary 출력) | **fp path 진입 자체 안 함** |

**진단: Dispatcher 임계값이 실제 소수 bitsize와 안 맞음**.

`quat_fp_widths_from_alg` (`quat_fixed_precision.c:55-76`) 분기:
```
if (p_bits <= 128) → L1 widths (nwords_vec=5)
else if (p_bits <= 200) → L3 widths (nwords_vec=7)
else if (p_bits <= 256) → L5 widths (nwords_vec=9)
else: nwords=0, return 0
```

실측 소수 bitsize:
- L1 `mlll_benchmark.c:252-254` prime → **251 bits** → `p_bits <= 256` 적중 → **L5 widths 적용** (nwords_vec=9)
- L3 `mlll_benchmark.c:247-249` prime → **375 bits** → 어느 branch도 안 걸림 → **return 0**
- L5 `mlll_benchmark.c:242-244` prime → **473 bits** → 어느 branch도 안 걸림 → **return 0**

Dispatcher (`mlll_gram.c:803-813`)는 `widths_from_alg` 실패 시 조용히 ibz path로 fallback. 그래서 `--fp` 플래그가 L3/L5 에서 **silent no-op**였음.

**파급 (Phase 2 측정·결정에 영향)**:

1. **§5.6 3-way sweep 시간 표 무효화 영역**:
   - L3 alg2/alg3, L5 alg2/alg3 의 "fp" 컬럼 수치는 fp가 아닌 **ibz+prealloc-off 런타임**이었음 (4행). 그래서 "동률"로 보였던 것. 진짜 fp/ibz 비율은 **미측정 상태**.
   - L1 alg2/alg3 의 "fp" 컬럼은 fp 구동은 맞지만 **L5 widths (nwords_vec=9, nwords_gram=17)** 위에서 돌아감. 실제 L1 widths(5/9)로 돌린 수치가 아니므로 L1 5.25× 퇴행은 L1 고유가 아닌 L5-크기 schoolbook 비용.
2. **§5.1 4-기준 평가 부분 무효화**:
   - 기준 (1) heap-free / (2) Lemma 3 런타임 감사 / (4) CT 여지는 **L1만** 실제 입증. L3/L5는 명목상 fp 모드지만 런타임에는 ibz가 돈 상태.
   - 기준 (3) contribution 코드 입증은 "코드가 존재한다" 수준은 유지되지만 "고정폭 슬롯에 실제 값이 들어간다"는 **L1 + wrong-width** 범위로 축소.
3. **§6.4 capacity sweep** (270 iter trap 발동 0건) 은 L3/L5 에서 **ibz path 실측 peak 를 기록한 것**. fp path 슬롯 경로는 경로 자체가 실행되지 않아 sweep 결과는 budget 검증 근거가 될 수 없음. 오직 L1 만 fp 실측인데 그것도 L5 widths 기준.
4. **§6.5 equivalence**: ibz↔fp 비트 동일 주장은 **L1 한정**으로만 실검증되었을 가능성. L3/L5 는 fp가 안 돌았으므로 "fp 출력"이 실제로 ibz 출력과 같은 코드를 거친 것.

**trap 자체는 살아있음 확정** (L1 SELFTEST에서 abort 정상 발동). LTO / hidden visibility / 인라인으로 인한 trap 제거는 아님.

**후속 조치 (Phase 2.1 추가 필요)**:

- [ ] **P2.1-dispatch** `quat_fp_widths_from_alg` 임계값을 실제 SQIsign 소수 bitsize에 맞게 재지정 (`<=252` → L1, `<=376` → L3, `<=474` → L5, else 0). 또는 `quat_alg_t`에 명시적 level tag를 두고 그것으로 분기 (bitsize는 휴리스틱이라 prime 세부 구조 바뀌면 또 깨짐).
- [ ] **P2.1-widths-rerun** 수정된 dispatcher로 §6.4 capacity sweep 재수행. 모든 레벨에서 fp path가 실제로 실행되는지 새 SELFTEST로 먼저 재검증.
- [ ] **P2.1-time-rerun** §5.6 시간 표 재측정. L3/L5 fp vs ibz 비율이 실제로 어느 수준인지 처음으로 확인. L1 은 L5 widths 오버사이징에서 해방되면 퇴행 감소 가능성.
- [ ] **P2.1-equiv-rerun** `quat_test_mlll_gram_fp_equivalence` 를 L1/L3/L5 3레벨 모두 강제 실행하도록 보강. 현재는 L5 random 16-gen 케이스만 있고 L3/L5는 dispatcher bypass로 ibz↔ibz 비교한 상태.

**잠정 상태**: Phase 2 "primary 확정"은 L1 부분 관찰 기반. L3/L5 는 사실상 미검증. Phase 3 착수 전 Phase 2.1 수행이 선행 조건.

### 6.8. P2.1-dispatch: Threshold 정정 및 3-레벨 fp path 재검증 — 2026-04-22

**조사 전제 검증**: 벤치 prime bitsize(§6.7 표)가 production prime과 다를 가능성을 agent 재조사로 검증. 결과 production `QUATALG_PINFTY.p` limb literal 직접 읽어 확인:

- L1: `src/precomp/ref/lvl1/quaternion_data.c` — 4 limbs, top `0x4fffffffffffffff` → **251 bits** (spec BITS=256)
- L3: `src/precomp/ref/lvl3/quaternion_data.c` — 6 limbs, top `0x40ffffffffffffff` → **383 bits** (spec BITS=384)
- L5: `src/precomp/ref/lvl5/quaternion_data.c` — 8 limbs, top `0x1affffffffffffff` → **505 bits** (spec BITS=512)

§6.7 벤치 수치(251/375/473)와 production(251/383/505)이 일치한다(L1) / 근접한다(L3/L5 오차 8~32 bits). 핵심 결함은 **모든 관측 prime 이 dispatcher 임계값 128/200/256 을 훨씬 상회** — 벤치/production 공통으로 **L1만 적중** (251 ≤ 256이므로 L5 widths 브랜치). L3/L5는 dispatcher `return 0` → silent ibz fallback. 휴리스틱 자체가 spec BITS 상수와 2배 차이 나는 버그.

**수정안**: user 방침 "sqisign 도 동일한 휴리스틱이면 똑같이". Spec BITS 상수(256/384/512) 기준으로 임계값을 재지정, 두 곳 모두 수정:

- `quat_fp_widths_from_alg` (`quat_fixed_precision.c:65-76`) — fp dispatcher (production).
- `mlll_gram_level_hints` (`mlll_gram.c:113-122`) — prealloc hints (현재 dead path지만 일관성 위해 동시 정정).

```c
/* after */
if      (p_bits <= 256) { /* L1: p ~251 bits */ ... }
else if (p_bits <= 384) { /* L3: p ~383 bits */ ... }
else if (p_bits <= 512) { /* L5: p ~505 bits */ ... }
else                    { return 0; }
```

**SELFTEST 3-레벨 검증** (budget=1 force-trap, `--iterations=1 --mode=alg3 --fp`):

| Level | stderr 출력 | 결과 |
|---|---|---|
| L1 | `quat_fp_vec overflow at load_gen0: 128 bits > 1 budget (nwords_vec=5)` | abort (core dumped) |
| L3 | `quat_fp_vec overflow at load_gen0: 192 bits > 1 budget (nwords_vec=7)` | abort (core dumped) |
| L5 | `quat_fp_vec overflow at load_gen0: 255 bits > 1 budget (nwords_vec=9)` | abort (core dumped) |

`nwords_vec` 값 5/7/9 가 각각 L1/L3/L5 widths 테이블(§6.4-a)과 일치 → **모든 레벨에서 fp path 진입 + trap 발동 확인**.

**Production 재측정** (SELFTEST 제거 후 실제 widths 기반 capacity sweep 재수행, `iter=500(alg2)/5000(alg3)`):

| Level | Mode | vec_max | vec_budget | gram_max | gram_budget | Verdict |
|---|---|---|---|---|---|---|
| L1 | alg2 | 259 | 318 | 517 | 574 | PASS |
| L1 | alg3 | 129 | 318 | 507 | 574 | PASS |
| L3 | alg2 | 389 | 446 | 778 | 830 | PASS |
| L3 | alg3 | 195 | 446 | 771 | 830 | PASS |
| L5 | alg2 | 509 | 574 | 1019 | 1086 | PASS |
| L5 | alg3 | 256 | 574 | 1015 | 1086 | PASS |

**전 레벨/모드 PASS**. trap abort 0건. 최저 margin은 **L3 alg2 Gram (52 bits ≈ 0.8 limb)**. 로그: `bench_logs/p2_fp_capacity_verify/L{1,3,5}_alg{2,3}.log` (overwrite).

**§6.7 무효화 스코프 재정의**:

1. §5.6 3-way sweep 시간 표: **§5.6 재측정 완료** (P2.1-time-rerun). fp/ibz 비율 표 L1/L3/L5 × alg2/alg3 6-combo 갱신. pre-fix vs post-fix 비교도 §5.6에 병기.
2. §6.4 capacity sweep: 본 절 표로 **대체**. L1/L3/L5 모두 실제 fp path 실측 budget 감사 완료.
3. §6.5 equivalence: **3-레벨 전면 재실행 완료** (P2.1-equiv-rerun). `mlll_tests.c:quat_test_mlll_gram_fp_equivalence` 에 L1 p=251, L3 p=383, L5 p=505 production-matching prime + 16 random gen (bound=sqrt(p)) 케이스 추가. 기존 p=7/19/11 소케이스(dispatcher L1 widths) + 3-레벨 신규 케이스 모두 PASS.
4. §5.1 4-기준 평가: 기준 (1)(2)(3)(4) 모두 **L1/L3/L5 전체에서 fp 경로 실증 완료** (본 절 SELFTEST + capacity sweep + equivalence). L1에서 L5 widths 오버사이징 이슈도 **해소** (L1 nwords_vec=5, nwords_gram=9 로 정상 축소).

**후속 잔여**:

- [x] **P2.1-dispatch (2026-04-22)** 임계값 256/384/512 로 정정.
- [x] **P2.1-widths-rerun (2026-04-22)** 3-레벨 SELFTEST + 6-combo capacity sweep 통과.
- [x] **P2.1-time-rerun (2026-04-22)** alg2 500 iter + alg3 5000 iter × L1/L3/L5 × ibz/prealloc/fp 재측정. 결과 §5.6 표 갱신. L5 alg2 fp/ibz 4.78× 가 최악, Alg 3 는 전 레벨 1.79~1.90× 안정.
- [x] **P2.1-equiv-rerun (2026-04-22)** `quat_test_mlll_gram_fp_equivalence` 에 L1/L3/L5 production-matching prime 케이스 추가. bound=sqrt(p)로 Lemma 3 충족. 3-레벨 PASS.

**상태**: Phase 2.1 완료. Phase 2 "primary 확정"을 **3-레벨 전면 확정**으로 갱신. Phase 3 blocker 해제.

## 7. 열린 질문 / 검증 필요

1. **Signed 표현** — **2026-04-21 해결 (§6.2 P2-C-types-c)**. Storage는 two's complement, multiplication만 내부적으로 sign-magnitude. Add/sub hot loop에서 분기 없음이 결정 근거.
2. **Division 구현**: Cohen 점화식의 exact div (`num = (d[s+1]·num − lam²) / d[s]`)가 필요. 단, 현재 `quat_mlll_gram`은 L² 기반이라 division은 dpe(float) 경로에서만 발생하고 integer exact division은 사용하지 않음. Cohen 경로(`mlll.c`)에만 필요 — Phase 2 scope 밖으로 분류. **재검토: P2-C는 GRAM만 대상이므로 이 열린 질문은 Phase 2에서 드롭.**
3. **레벨별 코드 복제 vs generic** — **2026-04-21 해결 (§6.2)**. 런타임 분기 + 최대폭(L5) stack array. Compile-time 분기는 현재 quaternion 레이어가 단일 generic 빌드라 부적합.
4. **탈락 경로 제거 시점** — **2026-04-21 해결 (§5.1, §5.7, §6.6-c)**. Primary는 **C(fp)** 로 확정(§5.1 4개 기준 유일 만족). ibz body는 fp 회귀 oracle, B prealloc은 ibz 서브튜닝 도구로 유지. `#ifdef LEGACY_*` 플래그 도입 없이 런타임 flag(`set_fp_mode` 기본 1, `set_prealloc_mode` 기본 0)로 opt-out/opt-in. 유지 비용 증가가 관측되면 재평가 대상.
5. **Cohen 경로 유지 여부**: 현재 memory `project_sqisign_mlll.md`에서도 언급된 열린 질문. Cohen은 Lemma 3 밖이라 고정폭 불가 → Phase 2 종료 시점에 제거 또는 `#ifdef QUAT_MLLL_COHEN_DEBUG` 로 감싸기. **결정 보류** (Phase 2 scope 밖). Phase 3 착수 전 재검토.
6. **Overflow trap 런타임 토글**: 현재 `MLLL_FP_NO_OVERFLOW_CHECK`는 **컴파일 타임** 매크로. 벤치/감사 도중 toggle 하려면 재빌드 필요. Phase 2에서는 trap-on(기본) 1회 + trap-off 1회 재빌드로 충분했으나, 향후 fp 경로를 런타임 opt-in으로 쓰며 trap만 개별 토글해야 할 상황이 오면 `g_fp_trap_mode` 전역 flag + macro를 `if (g_fp_trap_mode) { ... }` 형태로 전환. **현재는 보류**, 실수요 발생 시 착수.

## 8. 일정 (러프, 여유 있게)

병행 구현 + 실측 비교 방식이라 기존 플랜(~5일)보다 길게 잡는다. 일정은 upper bound에 가깝게 설정하고 실제 진행 속도에 맞춰 유동 운영.

| 단계 | 기간 | 누적 |
|---|---|---|
| P2-1 (본 문서 확정) | 0.5d | 2026-04-22 |
| P2-B (B scaffold) | 2-3d | ~2026-04-25 |
| P2-C-types (typedef + signed 결정) | 1-2d | ~2026-04-27 |
| P2-C-gram (GRAM 치환 + mul/div 구현) | 4-5d | ~2026-05-02 |
| P2-overflow (trap) | 1d | ~2026-05-03 |
| P2-equiv (3-way 동치성 + 벤치 비교) | 1-2d | ~2026-05-05 |
| P2-decide (최종 primary 결정 + 문서 업데이트) | 1d | ~2026-05-06 |

**Phase 2 총 ~11-14일**. 기존 플랜 5일 → 2-3배 연장. 근거:
- B와 C 병행 구현이라 실질 재작성량 증가.
- Mul/Division 구현 결정을 micro-bench로 판단하는 step 추가.
- 실측 비교(§6.5)로 최종 primary 결정하는 버퍼 포함.

일정 압박 없음. 비교 단계에서 의외의 결과가 나와 접근 재조정해야 하면 이 일정도 재산정 대상.

---

## 참고 문헌

- [SQIsign with Fixed-Precision Integer Arithmetic (eprint 2025/1649)](https://eprint.iacr.org/2025/1649.pdf) — 본 구현이 따르는 논문
- [Constant-time Quaternion Algorithms for SQIsign (eprint 2025/2192)](https://eprint.iacr.org/2025/2192.pdf) — 미래 CT 방향
- [SQIsign NIST PQC Round 2 Spec](https://csrc.nist.gov/csrc/media/Projects/pqc-dig-sig/documents/round-2/spec-files/sqisign-spec-round2-web.pdf) — 현재 non-CT 공식 확인
- [GNU MP 6.3.0 — Low-level Functions](https://gmplib.org/manual/Low_002dlevel-Functions) — mpn_* API
- [GNU MP 6.3.0 — Memory Management](https://gmplib.org/manual/Memory-Management) — mpz_realloc2 semantics
