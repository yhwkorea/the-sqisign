# MLLL 작업 계획 (2026-04-19 기준)

논문: "Compact Quaternion Algorithms for SQIsign"
브랜치: `feat/mlll-ideal-operations`
직전 마일스톤: 2026-04-19 MLLL 모듈 감사 완료 ([AUDIT_2026-04-19_KO.md](AUDIT_2026-04-19_KO.md))

## 2026-04-30 정정

이 문서의 이전 판본이 두 가지 사실 오류를 포함하고 있었음:

1. **"Alg 1 IdealFiltration"은 논문에 존재하지 않음.** 논문의 algorithm 목록은 Alg 1 MLLL(커널) / Alg 2 CompactIdealMultiplication / Alg 3 RandomIdealGivenPrimeNorm / Alg 4 RandomEquivalentPrimeIdeal 4개이며 "Filtration"은 03Ideal.tex, 04Sampling.tex, append.tex 어디에도 등장하지 않음. P3-1, P3-2 항목 및 의사결정 로그 #2는 유령 항목으로 폐기.
2. **Alg 4 RandomEquivalentPrimeIdeal MLLL 버전은 이미 구현됨.** `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` 이 `lll_applications.c:213`에 존재하며 단위 테스트 PASS. Phase 3-1 커밋(2026-04-29)에 포함됨. P3-3 항목은 ✅ 완료로 정정.

진짜 미완성은 **hot-path 라우팅(P3-2 재정의)** 1건과 **end-to-end 동치성 테스트(P3-5)** 1건. 아래 로드맵의 Phase 3 절은 이 정정에 따라 다시 읽어야 함 — 원문은 history 보존 목적으로 유지.

## 현재 상태 스냅샷

### 구현 완료

| 논문 Alg | 구현 함수 | 상태 |
|---|---|---|
| — (커널) | `quat_mlll` (Cohen 정수 GSO) | ✅ 감사 완료 |
| — (커널) | `quat_mlll_gram` (L²-style Gram+dpe) | ✅ 감사 완료, **프로덕션 경로** |
| Alg 2 CompactIdealMultiplication | `quat_lattice_mul_mlll_gram` | ✅ 감사 완료 |
| (Alg 2 보조) | `quat_lattice_add_mlll_gram` | ✅ 감사 완료 |
| Alg 3 RandomIdealGivenPrimeNorm 부분 | `quat_lideal_create_mlll_gram` | ✅ 감사 완료 |
| (HNF 비교용) | `quat_lattice_mul_mlll`, `quat_lattice_add_mlll` | ✅ 유지 |

### 테스트/벤치

- **테스트 12/12 PASS** (감사 후 M5 커버리지 보강 포함)
- **L1/L3/L5 10k trials 실측 완료** (C1 픽스 반영 후) — vec 259/391/513, Gram 518/782/1026 bits
- `sqisign_bm_mlll --mode=alg2|alg3 --iterations=N` CLI 동작

### HNF 기반으로 남아있는 것 (2026-04-30 정정)

| 논문 Alg | 기존 HNF 함수 | MLLL 버전 |
|---|---|---|
| Alg 4 RandomEquivalentPrimeIdeal | `quat_lideal_prime_norm_reduced_equivalent` (`lll_applications.c:120`) | ✅ `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` (`lll_applications.c:213`) |
| (보조) `quat_lideal_reduce_basis` | `lll_applications.c:7` | ✅ `quat_lideal_reduce_basis_mlll_gram` (`lll_applications.c:41`) |

알고리즘 본체 단계에서 4개 paper algorithm 모두 MLLL_GRAM 구현 존재. 남은 작업은 hot-path 라우팅 + e2e 동치성 검증 — 아래 Phase 3 재정의 참고.

## 로드맵

### Phase 1 — 측정 보강 (완료: 2026-04-21)

**목표**: C1 픽스 반영된 정확한 bitsize 표를 L1/L3/L5 전체에 대해 확보.

- [x] **P1-1** L3 10k trial 벤치 (`--mode=alg2`, `--mode=alg3`) — `bench_logs_2026-04-19/L3_alg{2,3}_10k.log`
- [x] **P1-2** L5 10k trial 벤치 (동일) — `bench_logs_2026-04-19/L5_alg{2,3}_10k.log`
- [x] **P1-3** `project_sqisign_mlll.md` 메모리 B*_{1,3,5} 표 갱신 (2026-04-21)
- [x] **P1-4** `README_MLLL_BENCHMARK_KO.md` 숫자 업데이트 + C1 영향 명시 (2026-04-21)

**측정 결과 요약** (GRAM 경로, max bits, C1 픽스 후):

| Level | vec (Lemma 1) | Gram (Lemma 3) | Cohen 비교 | HNF 비교 |
|---|---|---|---|---|
| L1 | 259 | 518 | 3543 | 2160 |
| L3 | 391 | 782 | 5389 | 3293 |
| L5 | 513 | 1026 | 7102 | 4330 |

- C1 픽스는 Cohen 경로에 +515 bits (L1) 드러남, GRAM은 +1~2 bits — **typedef 폭 결정 유효**.
- 시간비 GRAM/Cohen: L3 5.5%, L5 4.0% (레벨↑ GRAM 우세↑).
- Alg 3 GRAM/HNF 시간비: L3 0.62, L5 0.61 (이미 HNF 대비 빠름).

**Why**: 감사 C1 픽스로 Cohen path max가 3028→3543 bits로 이동. L3/L5도 같은 폭 누락 가능성 있어 재측정. GRAM 경로는 거의 영향 없음을 실측 확인 — Phase 2 typedef 폭 그대로 사용 가능.

### Phase 2 — Fixed-precision 전환 (완료: 2026-04-21)

**목표**: `ibz_t`(GMP 동적) → 고정폭 산술로 내부 전환. B와 C 둘 다 구현하고 실측으로 primary 결정. 상세: [FIXED_PRECISION_DECISION_KO.md](FIXED_PRECISION_DECISION_KO.md)

- [x] **P2-1** 백엔드 결정 문서 초안 — `FIXED_PRECISION_DECISION_KO.md` (2026-04-22). A 탈락, B와 C 병행 구현 방침 확정.
- [x] **P2-B** B 구현: `ibz_t` + `mpz_realloc2` prealloc scaffold (2026-04-21). 동치성 검증 PASS, 벤치 결과 hint only 확인 — B 단독으로는 목표 미달성, reference oracle로 유지. 상세: `FIXED_PRECISION_DECISION_KO.md` §5.5.
- [x] **P2-C-types** Per-level 폭 매크로 + 런타임 descriptor(`quat_fp_widths_t`) + 단일 최대폭 array. Two's complement storage + sign-magnitude mul 확정. 상세: `FIXED_PRECISION_DECISION_KO.md` §6.2.
- [x] **P2-C-gram** `quat_mlll_gram_fp` 본체 + dispatcher + fp-native size-reduce/swap/Gram update. `vec4_dot_p`만 ibz scratch bypass (L1 중간값 예산 초과). 상세: §6.3.
- [x] **P2-overflow** `FP_CHECK_VEC`/`FP_CHECK_GRAM` 매크로로 mutation site 마다 budget 초과 시 `abort()`. `MLLL_FP_NO_OVERFLOW_CHECK` 빌드 플래그로 비활성 가능. **Budget 정정 (2026-04-22)**: 초기 공식 `(nwords-1)*64`는 margin narrative가 틀렸음이 확인됨(L3 vec peak 391b > 384b인데 trap 미발동). 정정된 공식 `nwords*64 - 2` (부호 1 + carry 1)로 L1/L3/L5 모두 48-61b 여유. 상세: `FIXED_PRECISION_DECISION_KO.md` §1, §6.4-c. Trap 감사 재수행은 §6.4-d에서 진행.
- [x] **P2-equiv** 2-way pairwise 테스트 3개(ibz↔HNF, ibz↔B, ibz↔fp)로 HNF↔ibz↔B↔fp transitive equivalence. `quat_test_mlll_gram_fp_equivalence` PASS.
- [x] **P2-decide** Primary = **C fp 경로** (`quat_mlll_gram_fp`) 최종 채택 (2026-04-21). `g_fp_mode = 1` 기본값. 결정 기준은 논문 contribution 4축(heap-free / Lemma 3 runtime 감사 / compact 연산 고정폭 입증 / CT 호환 여지) — fp가 유일 만족. ibz body는 fp 회귀 oracle로 보존. 시간 퇴행(L1 Alg 2 5.25×)은 수용 trade-off.

**결정 요약** (FIXED_PRECISION_DECISION_KO.md §5.1, §5.6, §5.7):
- Primary: `quat_mlll_gram_fp` (stack 고정폭 + schoolbook + overflow trap). `g_fp_mode = 1`.
- 회귀 oracle: `quat_mlll_gram_ibz` (기존 ibz body). `set_fp_mode(0)` 또는 `--fp 0`로 명시 opt-out 가능.
- B prealloc(`--prealloc`): ibz body 한정 서브튜닝 도구로 보존, 기본 OFF.
- **시간 trade-off**: L1 Alg 2 fp/ibz = 5.25× (trap-off 4.65×), L3/L5 동률. GMP의 1-limb special case 상실이 주요인. 향후 `mpn_mul_1` inline 재설계 시 fp body에 직접 삽입.

**Why**: KLKL25(eprint 2025/1649)의 핵심 주장 "compact 연산이 고정폭에 실제로 들어간다"를 코드로 입증하는 것이 Phase 2의 본래 기준. 시간은 이 기준의 일부가 아님 — primary는 4축 만족 여부로 결정. Phase 3는 fp primary 위에 직접 쌓는다.

**전제**: Phase 1 완료 (typedef 폭이 L3/L5 재측정으로 확정됨). ✅

### Phase 2.1 — Dispatcher 임계값 수정 + fp 경로 실재성 재검증 ✅ (완료, 2026-04-22)

**Why**: 2026-04-22 trap 감사(`FIXED_PRECISION_DECISION_KO.md` §6.7)에서 `quat_fp_widths_from_alg`의 임계값 128/200/256이 실제 SQIsign 소수 bitsize 251/383/505와 맞지 않아 **L3/L5에서 fp path가 조용히 미동작**(ibz fallback), **L1은 L5 widths로 오버사이징 동작**임이 SELFTEST 로 확정. Phase 2 측정 및 primary 결정 중 L3/L5 부분은 사실상 **fp 미검증 상태**였음.

- [x] **P2.1-dispatch (2026-04-22)** `quat_fp_widths_from_alg` 임계값을 spec BITS 상수 기준 256/384/512로 정정. `mlll_gram_level_hints` 도 일관성 위해 동시 정정. §6.8 참조.
- [x] **P2.1-widths-rerun (2026-04-22)** SELFTEST 3-레벨 모두 abort 확인 (L1:128b/nw5, L3:192b/nw7, L5:255b/nw9). 후속 capacity sweep 6-combo PASS (`bench_logs/p2_fp_capacity_verify/`). 최저 margin L3 alg2 Gram 52b.
- [x] **P2.1-time-rerun (2026-04-22)** §5.6 3-way 시간 표 재측정 (alg2=500 iter, alg3=5000 iter). L1 alg2 fp/ibz 5.25×→3.29× (L5 오버사이징 해소 효과), L3/L5 alg2 는 최초 실측치 (3.98×, 4.78×), alg3 는 전 레벨 1.79~1.90× 안정.
- [x] **P2.1-equiv-rerun (2026-04-22)** `quat_test_mlll_gram_fp_equivalence` 에 L1/L3/L5 production-matching prime (p=251/383/505) + 16 random gen (bound=sqrt(p)) 케이스 추가. 3-레벨 PASS.

**결과**: Phase 3 blocker 해제. Phase 2 primary 확정이 3-레벨 전면 실증 기반으로 갱신됨.

### Phase 3 — MLLL_GRAM 핫 패스 라우팅 + 동치성 (재정의, 2026-04-30)

**원문 (history)**: 원래 P3-1/P3-2는 "Alg 1 IdealFiltration" 설계/구현이었으나 해당 algorithm은 논문에 존재하지 않음(상단 정정 참조)이라 폐기. P3-3 (Alg 4 MLLL 변종)은 2026-04-29 Phase 3-1 커밋에 이미 구현됨.

**재정의된 목표**: paper 4개 algorithm의 MLLL_GRAM 구현이 sign/keygen/id2iso 핫 패스에서 실제로 호출되도록 라우팅하고, HNF 빌드와 MLLL_GRAM 빌드가 동치임을 e2e로 검증.

- [x] ~~**P3-1** Alg 1 IdealFiltration 설계 검토~~ — 폐기 (유령 algorithm)
- [x] ~~**P3-2** `quat_lideal_filtration_mlll_gram` 구현~~ — 폐기
- [x] **P3-3** Alg 4 RandomEquivalentPrimeIdeal MLLL 버전 — `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` (`lll_applications.c:213`) 구현 + 단위 테스트 5/5 PASS
- [x] **P3-2′** Hot-path 라우팅 — **완료 (2026-04-30)**. `SQISIGN_USE_MLLL_GRAM=ON` 시 `quaternion.h` 매크로 alias로 sign.c / keygen.c / encode_signature.c / id2iso.c / dim2id2iso.c 의 11개 호출부(직접 ideal API 8건 + `quat_lideal_lideal_mul_reduced` 3건)가 무수정 라우팅. 4개 함수 alias 대상: `quat_lideal_create`, `quat_lideal_reduce_basis`, `quat_lideal_prime_norm_reduced_equivalent`, `quat_lideal_lideal_mul_reduced`. 본체 정의 .c 파일 3개(`ideal.c`, `lll_applications.c`, `mlll_gram.c`) + 명시적 HNF↔MLLL 비교 파일 2개(`mlll_benchmark.c`, `mlll_tests.c`) + HNF normal-form 단언 테스트(`test/ideal.c`)에 `SQISIGN_MLLL_GRAM_IMPL` 가드. `sqisign_namespace.h`에 `_mlll_gram` 변종 4개 등록(per-level mangling 일관성). CI workflow `mlll-routing-check.yml` 신설(ubuntu-latest, matrix={OFF,ON}, nm symbol-level routing 검증). 두 빌드 PASS는 commit `7e8f0eb` 기준.
  - **Acceptance (macro alias 단계)**: 단순 source-grep은 매크로 치환 전이라 통과 안 됨. 대신 다음 셋 중 하나로 검증:
    1. **Preprocessed output**: `cmake -DSQISIGN_USE_MLLL_GRAM=ON` 빌드 후 `gcc -E src/signature/ref/lvlx/sign.c | grep -c quat_lideal_create_mlll_gram` 1+건 (가장 직접적)
    2. **Symbol reference**: `nm build/src/signature/.../sign.c.o | grep -c '_mlll_gram'` 1+건
    3. **Build-test**: `SQISIGN_USE_MLLL_GRAM=ON` + `OFF` 두 빌드 모두 `make test` PASS, **그리고** ON 빌드의 `sqisign_test_mlll_gram_e2e_equiv` (P3-5에서 신설) PASS
  - **빌드 검증 명령** (사람이 즉석 확인용):
    ```bash
    cmake -B build_mlll -DSQISIGN_BUILD_TYPE=ref -DSQISIGN_USE_MLLL_GRAM=ON
    make -C build_mlll && ctest --test-dir build_mlll
    ```
- [ ] **P3-4** 벤치 모드: `--mode=alg2|alg3|alg4` 통일 (alg1 = MLLL 커널은 unit-level만, e2e bench 대상 아님)
- [x] **P3-5** 동치성 테스트 — **완료 (2026-04-30, CI에서 검증)**:
  - 단위 레벨: `quat_lideal_reduce_basis` ↔ `_mlll_gram` lattice equality 25 trials PASS (`sqisign_test_mlll`).
  - e2e self-consistency: ON 빌드의 ctest 통과 = SQIsign protocol 정상 동작:
    - `sqisign_test_signature_lvl{1,3,5}` (keygen+sign+verify roundtrip)
    - `sqisign_test_nistapi_lvl{1,3,5}` (NIST API 계약)
    - `sqisign_test_threadsafety_lvl{1,3,5}` (동시 sign N회)
    - `sqisign_lvl{1,3,5}_SELFTEST` (random keygen+sign+verify cycles)
    - `sqisign_test_id2iso_lvl{1,3,5}`
  - byte-identical 비교는 의도적으로 안 함 (alg 4 reduce가 다른 representative를 뽑으면 sk/pk/sig 갈림). e2e가 통과한다는 사실 자체가 충분.
  - 잔여 위험: ON 빌드는 현재 `g_fp_mode=0` (ibz_t 백엔드)로 fallback 중. fp 백엔드는 P3-2′ 라우팅 후 production lattice 크기에서 width 부족 발견(749/1137/1507 bit gram vs 518/782/1026 budget) — fp 재산정은 별도 task로 분리.

**Why**: paper completeness 기준은 이미 알고리즘 본체 단계에서 충족. 핫 패스 미통합 상태에서는 Phase 1/2의 모든 측정/감사가 "보조 호출 측정"에 그침 (README "현재 상태" 표 참조). P3-2′ 통과 후에야 SQIsign 서명 성능에 대한 paper 주장이 빌드 가능 코드로 증명됨.

### Phase 4 — SQIsign 상위 통합 (2-3주)

**목표**: sign/verify 경로에서 MLLL_GRAM 사용, 실제 서명 성능 측정.

- [ ] **P4-1** `quat_lideal_*` 호출부에서 HNF → MLLL_GRAM 전환 (feature flag 우선)
- [ ] **P4-2** 전체 `test_sqisign` 통과 확인 (L1/L3/L5)
- [ ] **P4-3** KAT 벡터 재생성 및 검증
- [ ] **P4-4** sign/verify 성능 벤치 (cycle count)
- [ ] **P4-5** 기존 HNF 기반 PQC reference 대비 메모리/시간 회귀 없음 확인

### Phase 5 — PR 분리 & upstream (1주)

**목표**: merge-ready 단위로 쪼개서 SQIsign upstream 또는 fork 유지.

- [ ] **P5-1** PR 1: MLLL 커널 (`mlll.c`, `mlll_gram.c`, `l2.c` 변경, 테스트/벤치 포함)
- [ ] **P5-2** PR 2: Alg 2/3 래퍼 (`lattice_*_mlll_gram`, `lideal_create_mlll_gram`)
- [ ] **P5-3** PR 3: Alg 1/4 래퍼
- [ ] **P5-4** PR 4: Fixed-precision 전환 (별도, 큰 표면적)
- [ ] **P5-5** PR 5: SQIsign 상위 호출부 전환

## 의사결정 로그에 남겨야 할 것

아직 미결정 항목. 진행 전에 별도 문서/메모리로 기록해야 함.

1. ~~**Phase 2 백엔드** (A/B/C 중 어느 것)~~ — **2026-04-21 해결: baseline ibz_t 유지** (DECISION §5.1).
2. ~~**Alg 1 기존 HNF 구현 유무** (확인 필요, P3-1)~~ — **2026-04-30 폐기**: 논문에 IdealFiltration algorithm 자체가 없음. 상단 정정 참조.
3. **Constant-time 요구 수준** — Phase 2 결과 무관하게 fp scaffold가 CT 여지를 확보해 둠(`g_fp_mode` opt-in). 실수요 생길 때 재검토.
4. **PR 1회 큰 덩어리 vs 5개 분할** — upstream 리뷰어 피드백 받아본 후 결정
5. **Gram 외 Cohen 경로 유지 여부** — 비교/교차검증용으로 남길지, 제거할지. Phase 3 착수 전 결정 필요.

## 논문 매핑 업데이트 (참고)

감사 때 확인된 논문 대응 (변경 없음):
- `mlll.c` ↔ `03Ideal.tex:13-55` Cohen Alg 2.6.3 변형
- `mlll_gram.c` ↔ L²(NS09) + Pohst 생성자 loading
- `quat_lattice_mul_mlll_gram` ↔ Alg 2 CompactIdealMultiplication
- `quat_lideal_create_mlll_gram` ↔ Alg 3 내부 sublogic

## 일정 (러프)

| Phase | 기간 | 누적 | 상태 |
|---|---|---|---|
| P1 측정 | ~2d | 2026-04-21 | ✅ 완료 |
| P2 fixed-precision | ~11-14d → 1d | 2026-04-21 | ⚠️ **부분완료** (L1 한정 검증, 2026-04-22 trap 감사로 L3/L5 fp 미동작 발견) |
| P2.1 dispatcher 수정 | ~1-2d → 1d | 2026-04-22 | ✅ **완료** (dispatch/widths/time/equiv 전 축 재검증, 3-레벨 PASS) |
| P3 핫패스 라우팅 + 동치성 (재정의) | ~3-5d | ~2026-05-05 | P3-3 ✅, P3-2′/4/5 착수 가능 |
| P4 SQIsign 통합 | ~3w | ~2026-05-26 | |
| P5 PR 분리 | ~1w | ~2026-06-02 | |

Phase 2가 예상 11-14일 대신 1일에 종결된 이유: B(prealloc)가 §5.5에서 즉시 탈락했고 C(fp) 전체 구현 후 sweep에서도 승격 조건 미달이라 추가 튜닝 사이클 없이 decision 확정. Phase 3 일정은 2026-04-30 정정으로 P3-1/P3-2(IdealFiltration)가 폐기되어 단축됨 — 남은 작업은 매크로 alias + e2e 동치성 스크립트.
