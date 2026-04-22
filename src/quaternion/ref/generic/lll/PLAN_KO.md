# MLLL 작업 계획 (2026-04-19 기준)

논문: "Compact Quaternion Algorithms for SQIsign"
브랜치: `feat/mlll-ideal-operations`
직전 마일스톤: 2026-04-19 MLLL 모듈 감사 완료 ([AUDIT_2026-04-19_KO.md](AUDIT_2026-04-19_KO.md))

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

### HNF 기반으로 남아있는 것

| 논문 Alg | 기존 HNF 함수 | MLLL 버전 |
|---|---|---|
| Alg 1 IdealFiltration | (확인 필요) | **미구현** |
| Alg 4 RandomEquivalentPrimeIdeal | `quat_lideal_prime_norm_reduced_equivalent` (`lll_applications.c:48`) | **미구현** |

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

### Phase 2.1 — Dispatcher 임계값 수정 + fp 경로 실재성 재검증 (진행 중, 2026-04-22)

**Why**: 2026-04-22 trap 감사(`FIXED_PRECISION_DECISION_KO.md` §6.7)에서 `quat_fp_widths_from_alg`의 임계값 128/200/256이 실제 SQIsign 소수 bitsize 251/383/505와 맞지 않아 **L3/L5에서 fp path가 조용히 미동작**(ibz fallback), **L1은 L5 widths로 오버사이징 동작**임이 SELFTEST 로 확정. Phase 2 측정 및 primary 결정 중 L3/L5 부분은 사실상 **fp 미검증 상태**.

- [x] **P2.1-dispatch (2026-04-22)** `quat_fp_widths_from_alg` 임계값을 spec BITS 상수 기준 256/384/512로 정정. `mlll_gram_level_hints` 도 일관성 위해 동시 정정. §6.8 참조.
- [x] **P2.1-widths-rerun (2026-04-22)** SELFTEST 3-레벨 모두 abort 확인 (L1:128b/nw5, L3:192b/nw7, L5:255b/nw9). 후속 capacity sweep 6-combo PASS (`bench_logs/p2_fp_capacity_verify/`). 최저 margin L3 alg2 Gram 52b.
- [ ] **P2.1-time-rerun** §5.6 3-way 시간 표 재측정. L3/L5 fp vs ibz 실비율 최초 확인
- [ ] **P2.1-equiv-rerun** `quat_test_mlll_gram_fp_equivalence` 3레벨 강제 실행 보강

**blocker**: Phase 3 착수 선행 조건. Phase 2 결론 중 L3/L5 관련 수치 중 time/equiv 2축만 남은 상태 (capacity/기준평가는 §6.8에서 해소).

### Phase 3 — Alg 1/4 MLLL 버전 구현 (1-2주)

**목표**: 논문 Alg 1/4를 MLLL 경로로 구현. 기존 HNF 경로와 동치성 테스트.

- [ ] **P3-1** Alg 1 IdealFiltration 설계 검토 (논문 03Ideal.tex + 04Sampling.tex 재정독). HNF 대조 기준이 필요하면 그때 `grep quat_lideal_filtration` 한 번이면 충분 — 별도 조사 step 아님.
- [ ] **P3-2** `quat_lideal_filtration_mlll_gram` 구현
- [ ] **P3-3** Alg 4 RandomEquivalentPrimeIdeal MLLL 버전
  - 기존 `quat_lideal_prime_norm_reduced_equivalent` (`lll_applications.c`) → `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` 분기
- [ ] **P3-4** 각각에 대해 벤치 모드 추가 (`--mode=alg1`, `--mode=alg4`)
- [ ] **P3-5** 동치성 테스트 (HNF ↔ MLLL_GRAM)

**Why**: 논문의 completeness. Alg 2/3만 MLLL로 바꾸고 Alg 1/4 HNF에 남겨두면 signature 성능 이득이 부분적.

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
2. **Alg 1 기존 HNF 구현 유무** (확인 필요, P3-1)
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
| P2.1 dispatcher 수정 | ~1-2d | 2026-04-23~24 | **착수 필요** (Phase 3 blocker) |
| P3 Alg 1/4 | ~2w | ~2026-05-05 | P2.1 완료 후 착수 |
| P4 SQIsign 통합 | ~3w | ~2026-05-26 | |
| P5 PR 분리 | ~1w | ~2026-06-02 | |

Phase 2가 예상 11-14일 대신 1일에 종결된 이유: B(prealloc)가 §5.5에서 즉시 탈락했고 C(fp) 전체 구현 후 sweep에서도 승격 조건 미달이라 추가 튜닝 사이클 없이 decision 확정. 실제 일정은 Alg 1 기존 구현 조사 결과에 따라 달라짐.
