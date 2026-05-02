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

### Phase 2 — Fixed-precision 전환 (3-5일)

**목표**: `ibz_t`(GMP 동적) → 레벨별 고정폭 정수로 내부 산술 전환. Constant-time 가능성 확보 + 힙 할당 제거.

- [ ] **P2-1** 백엔드 선택 결정 문서 (별도 `FIXED_PRECISION_DECISION_KO.md`)
  - 후보 A: GMP `mpn_*` lowlevel API (익숙, heap 유지)
  - 후보 B: `ibz_t` + `mpz_realloc2` (최소 변경, heap 유지)
  - 후보 C: bare `uint64_t[N]` (constant-time 가능, 재작성 많음)
- [ ] **P2-2** 레벨별 typedef 정의 (`quat_b_vec_Lk_t`, `quat_b_gram_Lk_t`)
  - L1: vec 5×u64(320b), Gram 9×u64(576b)
  - L3: vec 7×u64(448b), Gram 13×u64(832b)
  - L5: vec 9×u64(576b), Gram 17×u64(1088b)
- [ ] **P2-3** `quat_mlll_gram` 내부 산술을 typedef로 치환 (GRAM 경로만 우선)
- [ ] **P2-4** Overflow trap: B* 초과 시 `abort()` (릴리즈 빌드에서도 활성)
- [ ] **P2-5** 동치성 테스트: 기존 `ibz_t` 경로와 결과 일치 (`quat_test_mlll_gram_equivalence` 확장, trials 100+)

**Why**: 논문 contribution의 핵심이 "compact 연산이 고정폭에 실제로 들어간다"는 실증. 현 구현은 여전히 `ibz_t`라 bitsize만 작을 뿐 메모리 레이아웃은 동일. 고정폭 전환 후에야 성능/메모리 이득이 드러남.

**전제**: Phase 1 완료 (typedef 폭이 L3/L5 재측정으로 확정되어야 함).

### Phase 3 — Alg 1/4 MLLL 버전 구현 (1-2주)

**목표**: 논문 Alg 1/4를 MLLL 경로로 구현. 기존 HNF 경로와 동치성 테스트.

- [ ] **P3-1** Alg 1 IdealFiltration 설계 검토 (논문 03Ideal.tex + 04Sampling.tex 재정독)
  - 기존 HNF 기반 filtration 함수 존재 여부 먼저 조사
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

1. **Phase 2 백엔드** (A/B/C 중 어느 것)
2. **Alg 1 기존 HNF 구현 유무** (확인 필요, P3-1)
3. **Constant-time 요구 수준** — Phase 2에서 C 백엔드로 갈지 여부의 전제
4. **PR 1회 큰 덩어리 vs 5개 분할** — upstream 리뷰어 피드백 받아본 후 결정
5. **Gram 외 Cohen 경로 유지 여부** — 비교/교차검증용으로 남길지, 제거할지

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
| P2 fixed-precision | ~5d | 2026-04-26 | 진행 예정 (백엔드 결정 문서 착수) |
| P3 Alg 1/4 | ~2w | 2026-05-10 | |
| P4 SQIsign 통합 | ~3w | 2026-05-31 | |
| P5 PR 분리 | ~1w | 2026-06-07 | |

실제 일정은 Phase 2 백엔드 결정과 Alg 1 기존 구현 조사 결과에 따라 크게 달라짐.
