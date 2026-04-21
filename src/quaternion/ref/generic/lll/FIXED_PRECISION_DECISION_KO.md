# Phase 2 Fixed-Precision 백엔드 결정 (초안)

논문: "Compact Quaternion Algorithms for SQIsign"
Phase: [PLAN_KO.md](PLAN_KO.md) Phase 2, 항목 P2-1
직전: [Phase 1 완료 2026-04-21](PLAN_KO.md#phase-1--측정-보강-완료-2026-04-21), B\*\_{1,3,5} 폭 확정

## 1. 목표와 전제

Phase 1에서 10k trials × L1/L3/L5 측정으로 GRAM 경로 `vec`/`Gram` 폭을 확정했다(아래 표). Phase 2는 이 폭을 **실제로 고정폭 정수 레이아웃에 넣어** MLLL 내부 산술을 힙 할당 없이 돌리는 것.

| Level | GRAM vec max | GRAM Gram max | 제안 typedef |
|---|---|---|---|
| L1 | 259 bits | 518 bits | vec 5×u64 (320b), Gram 9×u64 (576b) |
| L3 | 391 bits | 782 bits | vec 7×u64 (448b), Gram 13×u64 (832b) |
| L5 | 513 bits | 1026 bits | vec 9×u64 (576b), Gram 17×u64 (1088b) |

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

### 5.1. 채택: B 먼저 구현 → C 구현 → 벤치/동치성 비교 → 최종 primary 결정

이번 단계에서는 **primary를 선험적으로 못 박지 않는다**. 대신 B와 C를 모두 동작하는 상태까지 만들고, 3-way 벤치(기존 ibz_t / B prealloc / C fixed-width)로 (성능, 힙 할당 횟수, overflow trap 빈도, 코드 변경 규모)를 실측한 뒤 **Phase 2 종료 시점**에 primary 결정.

**순서**:

1. **B 단계 (먼저)** — `ibz_t` 유지한 채 `mpz_realloc2` prealloc을 `quat_mlll_gram` 진입부에 삽입. 힙 할당 빈도 감소가 우선 목표. 결과 동치성 검증 scaffold 역할.
2. **C 단계** — Per-level `quat_b_vec_L{k}_t`/`quat_b_gram_L{k}_t` 도입, GRAM 내부 산술을 고정폭 `digit_t[N]`으로 치환. B 단계 결과와 비트 단위 동치성 검증.
3. **비교 단계** — 같은 입력 10k trials로 B와 C의 벤치 수치 나란히 뽑고 장단점 실측 근거로 정리. 최종 primary 채택.

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
- **B만으로는 Phase 2 목표(힙 할당 제거, 고정폭 성능 이득) 달성 불가**. 실측으로 "hint only"가 확인되었으므로 더 이상 B 튜닝에 시간 쓰지 않는다.
- B scaffold는 C 구현의 reference oracle로 남긴다 (`quat_mlll_gram_set_prealloc_mode(0/1)`로 토글; 기본 0). C 구현 중 산출물이 B=1 출력과 비트 단위 동일해야 정상.
- **primary는 C 쪽 가능성 높음**. 다만 C 구현 전에 최종 확정 금지 (P2-decide에서 3-way 비교 후 결정).

**Why**: 본 확인으로 Phase 2 전체 성공이 C 구현의 품질에 의존하는 상황이 됨. C가 실패하면 fixed-precision 전환 자체를 재고해야 하므로 P2-C-gram에서 보수적 구현(schoolbook + overflow trap) 우선, 최적화는 P2-decide 이후로 미룬다.

## 6. 구현 단계

일정은 §8 참조. 각 단계는 독립 커밋 단위.

### 6.1. P2-B: B 단계 구현 (`mpz_realloc2` prealloc scaffold) — 2026-04-21 완료

- [x] **P2-B-a** `quat_mlll_gram` 진입부에서 `b[]`/`G[][]`/`X`/`tmp` 각 `ibz_t`에 `mpz_realloc2(x, B*_{k}+64)` 호출. Per-level 폭(L1 vec 320/Gram 576, L3 vec 448/Gram 832, L5 vec 576/Gram 1088)은 `alg->p` bitsize에서 유도. Global flag `quat_mlll_gram_set_prealloc_mode()` 노출.
- [x] **P2-B-b** In-process 토글 비교로 비트 동일 검증 — `quat_test_mlll_gram_prealloc_equivalence` 20 trials × (alg2 + alg3) PASS.
- [x] **P2-B-c** 10k trials 벤치 재실행 (`bench_logs_2026-04-21_P2B/`). 시간 변화 기록: §5.5 표. Heap profiler 카운터는 **필요 없음** — 시간이 이미 hint only를 증언.
- [x] **P2-B-d** `--prealloc` 플래그로 toggle. 기존 `--mode=alg2|alg3` 체계 유지하고 전역 mode만 토글.

**측정 결론**: §5.5 참조. `mpz_realloc2` hint는 실측상 효과 없음 또는 퇴행(L1 alg3 1.47× 느려짐). B는 C 구현의 reference oracle로 남기되, primary candidate로서의 경쟁력 소멸.

### 6.2. P2-C-types: Per-level typedef 정의

- [ ] **P2-C-types-a** `src/quaternion/ref/generic/lll/include/quat_fixed_precision.h` 신설.
- [ ] **P2-C-types-b** per-level 매크로:
  ```c
  // Level 1: vec 5×u64 (320 bits ≥ 259), Gram 9×u64 (576 bits ≥ 518)
  #define NWORDS_QUAT_VEC_L1  5
  #define NWORDS_QUAT_GRAM_L1 9
  #define NWORDS_QUAT_VEC_L3  7
  #define NWORDS_QUAT_GRAM_L3 13
  #define NWORDS_QUAT_VEC_L5  9
  #define NWORDS_QUAT_GRAM_L5 17

  typedef digit_t quat_b_vec_L1_t[NWORDS_QUAT_VEC_L1];
  typedef digit_t quat_b_gram_L1_t[NWORDS_QUAT_GRAM_L1];
  // ... L3, L5 동일 패턴
  ```
- [ ] **P2-C-types-c** Signed 표현 결정: MSB = 부호 bit 2's complement, 또는 별도 `sign` field. Phase 1 실측에서 vec 좌표가 음수도 나오므로 signed 필수. **결정은 P2-C-types 작업 시점에 실측으로 판단** (초기 후보: two's complement).

### 6.3. P2-C-gram: GRAM 내부 산술 C 구현

- [ ] **P2-C-gram-a** `mlll_gram.c`의 `b[]`, `G[][]` 타입을 `quat_b_vec_L{k}_t`/`quat_b_gram_L{k}_t`로 교체. Level 분기는 컴파일 타임 매크로 (`SQISIGN_LEVEL`).
- [ ] **P2-C-gram-b** 내적 `vec4_dot_p`, size-reduce, swap, Gram update 모두 고정폭 mul/add/sub으로 재작성.
- [ ] **P2-C-gram-c** Mul 구현 결정: schoolbook 자작 vs `mpn_mul` wrapper. **P2-C-types 완료 후 micro-bench로 판단** (mpn_mul 호출 overhead vs schoolbook inline 비교).
- [ ] **P2-C-gram-d** Division: Cohen 점화식의 exact div. schoolbook(Knuth D) vs `mpn_divrem_1` wrapper. P2-C-gram-c와 동일 판단 프로세스.

### 6.4. P2-overflow: Overflow trap (C 경로 전용)

- [ ] **P2-overflow-a** 각 연산 후 MSB가 sign bit + reserved margin을 초과하면 `fprintf(stderr) + abort()`. release 빌드에서도 활성.
- [ ] **P2-overflow-b** 별도 `--enable-fatal-overflow=off` 빌드 플래그로 벤치 중 비활성 가능 (성능 측정용).

### 6.5. P2-equiv: 3-way 동치성 + 성능 비교

- [ ] **P2-equiv-a** `quat_test_mlll_gram_equivalence`를 **3-way**로 확장: HNF (reference) vs ibz_t GRAM (baseline) vs B prealloc vs C fixed-width. 100 trials.
- [ ] **P2-equiv-b** 10k trials 벤치에서 C 경로 overflow trap 0회 발생 확인.
- [ ] **P2-equiv-c** 같은 입력 10k trials로 4경로(HNF / ibz_t GRAM / B prealloc / C fixed-width) 시간/메모리/realloc 횟수 비교표 생성.

### 6.6. P2-decide: 최종 primary 결정 및 문서화

- [ ] **P2-decide-a** P2-equiv-c 표를 근거로 **B vs C primary 최종 채택**.
- [ ] **P2-decide-b** 본 문서(`FIXED_PRECISION_DECISION_KO.md`) §5.1 "결정" 업데이트 (초안 → 확정).
- [ ] **P2-decide-c** 채택되지 않은 경로의 제거 여부 판단 (B 채택 시 C 유지할지, 또는 반대). 열린 질문 §7-4 동시 해결.
- [ ] **P2-decide-d** `PLAN_KO.md` Phase 2 완료 체크 + Phase 3 착수 조건 정리.

## 7. 열린 질문 / 검증 필요

1. **Signed 표현**: two's complement vs `{digit_t[N]; int sign}` 분리. Two's complement는 sub 구현 단순, 그러나 sign extension 로직 필요. Phase 2 초반에 결정.
2. **Division 구현**: Cohen 점화식의 exact div (`num = (d[s+1]·num − lam²) / d[s]`)가 필요. schoolbook division(Knuth D) vs `mpn_divrem_1` 호출. 후자가 간단하지만 mpn_* 의존성 생김.
3. **레벨별 코드 복제 vs generic**: per-level typedef를 매크로로 하면 한 소스에서 컴파일 타임 분기 가능. 하지만 현재 `src/gf/ref/lvl{k}/` 구조는 소스 복제 방식. 통일성 위해 우리도 복제할지, 아니면 lll 레이어는 generic + 매크로로 갈지. **결정: 일단 generic + 매크로** (복제는 유지보수 부담).
4. **B 또는 C 중 탈락 경로 제거 시점**: §6.6 최종 결정 후 탈락 경로는 `-DLEGACY_*=1` 플래그로 유지할지 제거할지 판단. 장기적으로는 단일 경로만 유지(P5 PR 분리 전).
5. **Cohen 경로 유지 여부**: 현재 memoroy `project_sqisign_mlll.md`에서도 언급된 열린 질문. Cohen은 Lemma 3 밖이라 고정폭 불가 → Phase 2 종료 시점에 제거 또는 `#ifdef QUAT_MLLL_COHEN_DEBUG` 로 감싸기. **결정 보류**.

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
