# CLAUDE.md — 이 프로젝트 전용 규칙

이 파일은 `the-sqisign` 저장소에서 작업하는 모든 Claude 세션에 적용되는 **프로젝트 전용 규칙**입니다.
글로벌 `~/.claude/CLAUDE.md` 의 규칙을 모두 그대로 따르되, 아래 항목이 **추가 + 우선**합니다.

---

## 0. 사건 기록 — 왜 이 규칙이 생겼나

2026-04-29, 외부 감사로 다음 사실이 발각:
- README, commit msg에 "MLLL을 단독 기본 경로로 사용 / paper-strict 100%" 라고 적혀 있었음
- 실제로는 `grep -rn '_mlll_gram\b' src/signature/ src/id2iso/ src/sqisign.c` 결과 = **빈 출력**
- 즉 sign.c, keygen.c 어디에서도 MLLL 변종을 호출하지 않음
- README 와 빌드 현실이 **모순** → GitHub에 거짓 사실 push 된 상태로 발견됨

→ **자가검증 없이 "구현 완료" 단정한 결과**.

---

## 1. README는 항상 "현재 계획 + 실제 상태" 를 반영해야 한다 (HARD RULE)

### 1.1 README 최상단에 의무 섹션
모든 README 작성/수정 시 **최상단에 다음 4개 섹션** 을 둔다:

```markdown
# <프로젝트명>

## 현재 작업 계획 (yyyy-mm-dd 갱신)
[3–10 줄. 무엇을 / 언제까지 / 어느 phase / 다음 행동.]

## 현재 상태 — 완료 vs 미완료
| 항목 | 상태 | 검증 명령 |
|---|---|---|
| ... | ✅/🟡/❌ | grep / make 등 재현 가능한 한 줄 |

## 마지막 검증 시점
[yyyy-mm-dd. 무엇을 어떻게 검증했는지. 어떤 빌드 구성에서 측정했는지.]

## 알려진 미통합/거짓 위험 항목
[정직하게 적기. "이렇게 들리지만 실제로는 …" 형식.]
```

### 1.2 README 갱신 트리거
다음 중 하나라도 발생하면 **README 즉시 갱신** (별도 commit):
- 새 commit이 push 됨
- 빌드 default 옵션이 바뀜
- 호출 경로가 바뀜 (HNF↔MLLL 등)
- 새 알고리즘 본체가 추가됨
- 알고리즘과 핫 패스 통합 상태가 바뀜
- 검증 결과 (KAT, scheme test 등) 가 갱신됨
- 외부 감사로 거짓 주장이 발견됨

### 1.3 README에 쓸 때 절대 쓰지 말 것
다음 표현은 **알고리즘 한정 의미**가 빌드 전체로 확대될 위험이 크므로 **금지**:
- ❌ "단독 기본 경로로 사용" (MLLL 한정 표현 금지)
- ❌ "paper-strict 100%"
- ❌ "다 구현 완료"
- ❌ "NIST 통과" (그건 NIST가 평가)

대신 **범위를 명시**:
- ✅ "Algorithm 본체 완성, 핫 패스 통합 Phase 3 미수행"
- ✅ "단위 테스트 12/12 PASS (HNF/MLLL 4-way 동치 포함)"
- ✅ "KAT round-trip 통과 — HNF 빌드 기준"
- ✅ "Lemma `mlll-bound` 페이퍼 상한 만족 (L1=506b ≤ 514b)"

---

## 2. "구현 완료" 단정 금지 — 자가 grep 검증 의무 (HARD RULE)

### 2.1 "구현했다 / 끝났다 / 통과했다" 라고 보고하기 전 의무 검증
다음 명령들을 **반드시 실행하고 결과를 보고서에 첨부**:

```bash
# A. 새 함수가 호출되는지 확인 (가장 중요)
grep -rn '<새함수이름>\b' src/<핵심경로>/ --include='*.c'
# 결과가 비어있으면 = 호출 0건 = "구현 완료" 라고 말하면 안 됨

# B. 새 옵션 매크로가 분기에 쓰이는지 확인
grep -rn '<NEW_MACRO>' src/ --include='*.c'
# 결과가 비어있으면 = 매크로 dead

# C. 함수 본체 .c 파일 존재 확인
grep -rn '<함수이름>\b' src/ --include='*.c' | grep -v 'test\|benchmark'
# 본체 0건 = 헤더 선언만 있는 상태

# D. KAT/벤치 파일 출처 확인
git log --oneline apps/PQCgenKAT_sign.c
# 마지막 변경이 새 빌드 분기를 추가했는지
```

### 2.2 보고 시 다음 두 단어 중 하나로 범위 명시
- **"알고리즘 한정"** — 본체 + 단위 테스트만 검증됨, 핫 패스 미통합
- **"빌드 통합"** — 핫 패스에서 실제 호출되며 round-trip 검증 통과

### 2.3 검증 결과를 README ↔ commit msg ↔ 보고 셋이 일치
- README 의 "완료/미완료" 표
- commit msg 의 범위 표현
- Claude 응답의 단정 표현

세 곳이 **동일한 단어** 를 써야 함. 한 곳만 다르면 거짓 위험.

---

## 3. Phase 분리 명시 (HARD RULE)

이 프로젝트는 다음 phase 구분을 따른다. **Phase 끝났다고 단정하기 전 자가 검증 의무**:

| Phase | 의미 | "끝났다" 의 검증 |
|---|---|---|
| Phase 1 | 알고리즘 본체 작성 | `mlll_tests.c` 단위 테스트 PASS |
| Phase 2 | 백엔드/정밀도/dispatcher | 4-way 동치 + Lemma bound 측정 |
| Phase 3-1 | 누락 함수 본체 작성 | `grep -rn '함수\\b' src/ --include='*.c'` 본체 1+건 |
| Phase 3-2 | 핫 패스 라우팅 | `grep -rn '_mlll_gram\\b' src/signature/` 1+건 |
| Phase 3-3 | KAT 생성 binary | `apps/PQCgenKAT_sign.c` 에 USE_MLLL_GRAM 분기 |
| Phase 3-4 | 문서 동기화 | README ↔ PLAN ↔ 빌드 일치 |
| Phase 3-5 | Algorithm 1 변종 | (선택) |

다른 phase는 끝났는데 다음 phase 미시작 상태면 **README 의 phase 명시 + 정직 표기 의무**.

---

## 4. 거짓 발견 시 대응 (HARD RULE)

거짓 사실이 외부 감사 / 리뷰로 발견된 경우:

1. **즉시 README 정정 commit** — 거짓 표현 삭제, 정확한 범위로 교체
2. **공개 정정 commit msg** — "정정: <commit hash> 의 <표현> 은 <범위 한정> 의미였음" 형식
3. **`PRESENTATION_KO.md` 또는 `AUDIT_<날짜>_KO.md` 에 사건 기록** — 무엇을, 언제, 어떻게 발견했는지
4. **CLAUDE.md (이 파일) 의 규칙 강화 검토** — 같은 사고가 다시 안 일어나도록

거짓 commit 을 **rebase / amend 로 숨기지 말 것**. 정정 commit으로 명시 추적.

---

## 5. 코드 품질 (글로벌 CLAUDE.md 보완)

### 5.1 헤더에 함수 선언만 추가하지 말 것
- 새 함수 헤더 선언 (`*.h`) 은 **반드시 같은 PR/commit 안에 .c 본체와 함께** 작성
- 헤더만 있고 본체 없는 상태로 commit/push 금지
- (예외: 인터페이스 설계 단계라면 PLAN_KO.md 에 명시)

### 5.2 "사용처 0건" 옵션 금지
- CMake `option(...)` 또는 `#define` 매크로 추가 시 **반드시 같은 PR에 분기 사용처 ≥ 1건**
- dead option 발견 시 즉시 제거 또는 활성 분기 추가

### 5.3 KAT/측정 데이터 출처 commit msg 의무
- KAT 파일 (`*.rsp*`) 생성/갱신 시 **어느 binary 가, 어느 빌드 구성에서, 언제 만들었는지** commit msg 에 명시
- 측정값 (`L1=506b` 등) 보고 시 **빌드 구성 명시** (예: "build_default = USE_MLLL_GRAM=OFF")

---

## 6. 검증 명령 보존 (재현 가능성)

다음 명령들은 **항상 재현 가능** 해야 한다. 변경 시 README 부록에 동기화:

```bash
# Lemma bound 측정
cd build_default/test && ./sqisign_test_scheme_lvl{1,3,5}

# KAT round-trip
./sqisign_test_kat_lvl{1,3,5}

# 핫 패스 라우팅 검증 (Phase 3-2 후 1+건 나와야 함)
grep -rn '_mlll_gram\b' src/signature/ src/id2iso/ src/sqisign.c

# 함수 본체 존재 (Phase 3-1 후 본체 1+건 나와야 함)
grep -rn 'quat_lideal_reduce_basis_mlll_gram' src/ --include='*.c'
grep -rn 'quat_lideal_prime_norm_reduced_equivalent_mlll_gram' src/ --include='*.c'
```

---

## 7. 세션 시작 시 의무 점검

새 Claude 세션이 이 저장소에서 시작되면, **첫 응답 전에**:

1. `git log --oneline -5` — 최근 commit 확인
2. `grep -rn '_mlll_gram\b' src/signature/ src/id2iso/ src/sqisign.c` — 핫 패스 라우팅 상태
3. README 최상단 "현재 작업 계획" 섹션 갱신일자 확인
4. 위 셋 중 하나라도 README 와 모순되면 **사용자에게 보고 + 정정 제안**

---

**이 규칙은 글로벌 CLAUDE.md 보다 우선합니다.** (충돌 시 이쪽이 이김)

마지막 갱신: 2026-04-29 (Phase 3 진입 + 거짓 발견 + 규칙 신설)
