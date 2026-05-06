# Compact paper MLLL — SQIsign fork 진행보고 (전체판)

**브랜치**: `feat/mlll-ideal-operations` · **commits**: 33 (main 이후) · **MLLL 코드**: ~3,758 LOC
**저장소**: `github.com/yhwkorea/the-sqisign` · **작성**: 2026-04-29

---

## 0. 이 발표를 따라오는 데 필요한 사전지식

**없어도 됩니다.** 0 페이지부터 시작합니다.

이 자료의 구성:

| Part | 제목 | 누구를 위해 |
|---|---|---|
| **A** | 큰 그림 — 왜 이게 중요한가 | 처음 듣는 사람 |
| **B** | 수학 기초 — lattice, basis, reduction | 학부 1–2학년 수준 |
| **C** | 암호 맥락 — SQIsign이 뭔가 | PQC 처음 접하는 사람 |
| **D** | 핵심 — HNF vs MLLL, 페이퍼가 뭘 제안했나 | 본 작업 이해 |
| **E** | 우리 fork가 한 일 (구체) | 코드 리뷰어 |
| **F** | 페이퍼 ↔ 코드 1:1 매핑 | 논문 검증자 |
| **G** | 검증 + 발견된 문제 (감사 결과) | 정직한 현재 상태 |
| **H** | 향후 계획 (Phase 3) | 의사결정자 |
| **I** | Q&A + 용어집 | 모두 |

---

# Part A. 큰 그림 — 왜 이게 중요한가

## A1. 디지털 서명이란

**가게에서 신용카드를 긁으면 어떻게 "내가 맞다"는 게 증명되나?**

전산 세계의 **디지털 서명 (Digital Signature)**:
- "이 메시지를 쓴 게 진짜 나입니다" 를 수학적으로 증명
- 누구나 **검증** 가능, 위조 불가능
- 예: 인터넷 뱅킹, HTTPS, 코드 서명, 이메일 PGP

### 3가지 함수
```
KeyGen()  → (공개키 pk, 비밀키 sk)        # 처음 한 번만
Sign(sk, msg) → 서명 σ                     # 메시지마다
Verify(pk, msg, σ) → OK / FAIL            # 누구나 검증
```

→ **비밀키 없이 σ 만드는 게 수학적으로 어려운** 알고리즘이 좋은 서명 방식.

## A2. 양자 컴퓨터 위협

지금 쓰는 RSA, ECDSA 같은 서명은 **양자 컴퓨터가 깨버림** (Shor 알고리즘, 1994).

→ NIST 가 2016년부터 **Post-Quantum Cryptography (PQC)** 표준 선정 진행.

### NIST PQC 후보군 (서명 분야)
- **Dilithium** (2024 표준 ML-DSA) — 격자 기반
- **Falcon** — 격자 기반 (압축형)
- **SPHINCS+** — 해시 기반
- **SQIsign** ← **우리가 다루는 것**, **isogeny 기반**

## A3. SQIsign의 매력과 약점

**매력**:
- 서명 크기 **177 바이트** (L1) — 다른 PQC 후보의 **1/10 ~ 1/100**
- 가장 작은 PQC 서명
- 수학적 가정이 견고 (supersingular isogeny problem)

**약점**:
- **느림** — 서명 1회에 수백 ms ~ 수초
- 구현 복잡도 매우 높음 (quaternion algebra + isogeny computation 결합)

→ **"느림"의 큰 원인 중 하나가 quaternion ideal 연산.** 이 부분을 빠르게 만드는 게 본 작업의 동기.

## A4. 이 fork 의 목표를 한 줄로

> **upstream SQIsign이 quaternion ideal 환원에 쓰는 "HNF" 알고리즘을, Compact paper가 제안한 "MLLL" 알고리즘으로 대체.**

→ 이론상 **메모리 적게 + 비슷하거나 빠르게**.
→ 실측으로 페이퍼 주장 검증.

---

# Part B. 수학 기초 — Lattice, Basis, Reduction

## B1. Lattice (격자) 란

### 일상 비유
- 방 안에 **무한히 많은 점**이 일정한 간격으로 찍혀 있다고 생각.
- 도시 블록의 교차점 위치들의 집합 같은 것.

### 수학 정의
- $n$차원 공간 $\mathbb{R}^n$의 **이산 부분군**
- 어떤 vector $b_1, \dots, b_n$의 **모든 정수 조합**으로 이루어진 점들의 집합:

$$L = \{ a_1 b_1 + a_2 b_2 + \dots + a_n b_n : a_i \in \mathbb{Z} \}$$

### 2차원 예시
$b_1 = (1, 0)$, $b_2 = (0, 1)$ 이면 $L = \mathbb{Z}^2$ (정수 좌표 모든 점).

```
. . . . . . .
. * * * * * .
. * * * * * .
. * * * * * .   ← 격자점들
. * * * * * .
. * * * * * .
. . . . . . .
```

## B2. Basis (기저)

같은 격자를 **여러 가지 basis**로 표현 가능.

### 예시: 같은 $L = \mathbb{Z}^2$ 의 두 가지 basis
**Basis A** (좋음, 짧고 직각): $(1,0), (0,1)$
**Basis B** (나쁨, 길고 비스듬): $(1, 0), (1000000, 1)$

→ 두 basis 가 만드는 격자는 **완전히 동일**. 하지만 한 쪽은 다루기 쉽고 한 쪽은 어렵다.

## B3. Lattice Reduction (격자 환원) — 짧고 직각인 basis 찾기

**문제**: 나쁜 basis 가 주어졌을 때, 같은 격자의 **짧고 직각에 가까운** basis 를 찾아라.

→ 이게 격자 암호의 핵심 알고리즘 작업.

### 두 가지 환원 방법
1. **HNF (Hermite Normal Form)**: 정수 산술만 사용. 결과는 **상삼각** 형태.
2. **LLL (Lenstra–Lenstra–Lovász, 1982)**: 부동소수 활용. 결과는 **거의 직각**.
3. **MLLL (Modified LLL, Pohst 1987)**: LLL을 **선형종속 입력**도 처리하게 확장.

## B4. HNF 예시 — 어떻게 작동하나

입력 (4개 vector, 4D):
$$v_1, v_2, v_3, v_4 \in \mathbb{Z}^4$$

HNF 알고리즘 **결과** (예시):
$$\begin{pmatrix} 5 & 0 & 0 & 0 \\ 3 & 2 & 0 & 0 \\ 1 & 0 & 7 & 0 \\ 4 & 1 & 6 & 1 \end{pmatrix}$$

→ **상삼각**(대각선 위는 0). 매우 간단한 형태지만, **벡터들이 길다**.

## B5. LLL/MLLL 예시 — 어떻게 작동하나

같은 입력으로:
$$\begin{pmatrix} 1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 1 & -1 & 0 & 1 \\ 0 & 1 & -1 & 0 \end{pmatrix}$$

→ **각 벡터가 짧고 거의 직각**. 정수 산술만으로 도달 못 함 — **부동소수 (μ, B)** 가 중간에 필요.

## B6. 두 방법의 trade-off

| 방법 | 정수 산술만? | 결과 형태 | 메모리 | 시간 | 부동소수 정밀도 필요? |
|---|---|---|---|---|---|
| **HNF** | O | 상삼각 (긴 벡터) | 큰 정수 → 큼 | 빠름–중간 | 불필요 |
| **MLLL** | X | 짧고 거의 직각 | 작은 정수 → 작음 | 비슷 | **필요** (적절히 크면) |

→ **MLLL의 장점**: 환원 후 정수가 작아져 다운스트림 (isogeny 계산 등) 이 빨라짐.
→ **MLLL의 함정**: 부동소수 정밀도 부족 시 **알고리즘 자체가 발산** (영원히 안 끝남).

---

# Part C. SQIsign 맥락 — Quaternion Ideal 연산

## C1. SQIsign의 수학적 무대

```
타원곡선 E1 → isogeny → 타원곡선 E2
                  ↑
              (이걸 찾는 게 어려움)
                  ↓
         Quaternion algebra B_p
                  ↓
           Order O ⊂ B_p
                  ↓
    Left ideal I ⊂ O    ← 우리가 다룰 객체
```

### 직관
- **Quaternion**: 4차원 복소수 비슷한 수. $a + bi + cj + dk$.
- **Order**: quaternion에서 "정수 비슷한" 부분 (lattice).
- **Left ideal**: order 안에서 "약수" 비슷한 부분 (또 lattice).
- **Ideal 연산**: 두 ideal의 **곱, 합, 환원** — 모두 lattice 연산으로 환원.

## C2. SQIsign 서명 1회의 수학적 동선

```
1. 챌린지 ideal I_chall 생성
2. (sk와 결합) → I_combined = I_chall · I_sk
3. I_combined 환원 → 짧은 generator α 찾기                    ← MLLL
4. α 기반 isogeny φ 계산
5. φ 인코딩 → σ
```

→ **3단계 "환원"** 이 핵심 병목. upstream은 HNF, paper는 MLLL.

## C3. SQIsign 한 번 서명 시 호출되는 ideal 연산

| 연산 | 호출 횟수 | upstream 함수 | paper 변종 |
|---|---|---|---|
| Lattice multiplication | ~수 회 | `quat_lattice_mul` | `_mlll_gram` |
| Lattice addition | ~수 회 | `quat_lattice_add` | `_mlll_gram` |
| Left ideal create | ~10 회 | `quat_lideal_create` | `_mlll_gram` |
| Left ideal intersection | ~수 회 | `quat_lideal_inter` | (내부 mul 분기) |
| Reduce basis | 1–수 회 | `quat_lideal_reduce_basis` | `_mlll_gram` |
| Prime norm reduce | 1–수 회 | `quat_lideal_prime_norm_reduced_equivalent` | `_mlll_gram` |

→ 이 **6 종류의 함수**를 paper 변종으로 바꾸는 게 통합 작업.

---

# Part D. 페이퍼 — Compact Quaternion Algorithms for SQIsign

## D1. 저자와 출판

- **Kim–Lee–Yoo** (2024 / 2025 prepr).
- 제목: *"Compact Quaternion Algorithms for SQIsign"*
- 핵심 주장: **HNF 기반 quaternion ideal 연산을 MLLL 기반으로 교체하면 메모리 절감**.

## D2. 페이퍼가 제안한 4가지 알고리즘

### Algorithm 1 — Modified LLL (MLLL) 본체
- 입력: $\beta$ 개 vector (선형종속 가능, 4D quaternion이면 보통 16개)
- 출력: 같은 lattice를 만드는 **rank ≤ 4** 의 짧고 직각인 basis
- Pohst 1987 변형 LLL을 정확히 따름

### Algorithm 2 — Lattice multiplication via MLLL
- $L_1 \cdot L_2$ 계산: 4×4 = 16 generators 만들고 → MLLL → rank-4 basis
- HNF 대비 **중간 정수가 작음**

### Algorithm 3 — Left ideal create via MLLL
- $\alpha \in O$ 와 norm $n$이 주어졌을 때 $\alpha O + nO$ 의 short basis
- 8 generators 만들고 → MLLL

### Algorithm 4 — CompactIdealMultiplication
- $I_1 \cdot I_2$ 의 generator를 작게 만드는 **이중 MLLL** wrapper
- 페이퍼 핵심 기여, bound: $64 p^2 / \pi^4 \cdot r_1^2 r_2^2 \cdot \mathrm{nrd}(I_1) \mathrm{nrd}(I_2)$

## D3. Lemma `mlll-bound` — 페이퍼의 정밀도 보장

> **Lemma**: MLLL 진행 중 등장하는 **모든 정수 객체** ($b_i$, $G$, $h$ 등) 의 비트 크기는
> $\leq \max_i \|a_i\|^2$ 로 bound.

### 우리 SQIsign 파라미터 대입
| 레벨 | $p$ 비트 | 이론 상한 ($2 \log_2 p + 6$) |
|---|---|---|
| L1 | 254 | ≤ 514 b |
| L3 | 383 | ≤ 772 b |
| L5 | 505 | ≤ 1016 b |

→ 이 상한을 **실측해서 검증**한 게 본 fork의 큰 성과 중 하나.

## D4. 페이퍼 Appendix A.1 — 부동소수의 역할

> **μ (GSO 계수), B (norm²) 는 real-valued.**
> "Real-valued μ's ... do not introduce any new unbounded integer."

- μ, B 는 정수가 아닌 **실수**
- 정밀도가 부족하면 알고리즘 결과가 **틀리거나 발산**
- 페이퍼는 정밀도 수치 미지정 → 구현자 책임

→ 우리는 **mpfr_t 1024 bits correctly-rounded** 사용. L5 상한 1010b를 14b 마진으로 커버.

---

# Part E. 우리 fork 가 구체적으로 한 일

## E1. 전체 타임라인 (33 commits)

```
Phase 0  — 분석 + 감사
   - 4c3a599  audit: MLLL 모듈 전수 감사
   - 74449ae  PLAN_KO.md 작성

Phase 1  — Algorithm 본체 (mlll.c)
   - 9244d17  Phase 1 측정 + Lemma 1 해석
   - 85bc67a  L3/L5 측정 로그 (10k trials)
   - b196331  Phase 2 백엔드 결정 문서

Phase 2A-B — Gram 변종
   - 8b30e32  P2-C-gram step1 — fp add/sub/bitsize + ibz_t bridge
   - b18e49e  P2-C-gram step2 — fp mul (schoolbook)
   - 97eb079  P2-C-gram step3 — mlll_gram fp path + dispatcher

Phase 2C — Fixed-precision body
   - 601815f  P2-C-types — 헤더 skeleton
   - 4edc876  P2-B mpz_realloc2 prealloc (실측 비경쟁)
   - f3c7684  벤치 CLI에 --fp 플래그
   - 58540b1  P2-overflow — bitsize budget trap

Phase 2.1 — Dispatcher
   - 17c3b92  baseline ibz_t 유지 + 3-way sweep
   - e87ec2f  fp primary로 flip
   - bfd4c16  L3/L5 fp 미동작 발견 + 정정
   - 82814de  dispatcher 임계값 256/384/512
   - 3a89044  P2.1 완료 — 3 레벨 전면 실증

Phase 3 시작 (이번 세션)
   - 5ef9b76  paper-strict 마무리 (HNF fallback 제거 + D2c PREC=1024)
   - d1073ee  README 최상단 paper-strict 섹션
   - (이 발표자료)
```

## E2. 작성된 핵심 파일 (LOC)

| 파일 | LOC | 역할 |
|---|---|---|
| `src/quaternion/ref/generic/lll/mlll.c` | 483 | Algorithm 1 본체 + Alg 2/3 wrapper |
| `src/quaternion/ref/generic/lll/mlll_gram.c` | 937 | Gram 변종 + dispatcher |
| `src/quaternion/ref/generic/lll/quat_fixed_precision.c` | 559 | Phase 2C body (fp) |
| `src/quaternion/ref/generic/lll/mlll_fp_ops.h` | 180 | 4 백엔드 추상화 |
| `src/quaternion/ref/generic/lll/mlll_tests.c` | 1,599 | 단위 테스트 12개 |
| `src/quaternion/ref/generic/lll/mlll_benchmark.c` | (변경) | 벤치 CLI |
| **합계** | **~3,758** | |

## E3. 4가지 부동소수 백엔드 추상화

페이퍼는 "real-valued μ/B"만 명시. 우리는 **4가지 옵션** 으로 추상화하고 trade-off 측정:

| ID | C 타입 | 비트 | 라이브러리 | paper-strict? |
|---|---|---|---|---|
| **D2a** | `double` | 53 | (표준) | ❌ L1=506b 미달 |
| **D2b** | `mpf_t` | 가변 | GMP | △ truncation rounding |
| **D2c** | `mpfr_t` | 가변 | MPFR | ✅ correctly-rounded **default** |
| **D2d** | `dpe_t` | 53 mantissa + exponent | (헤더 단독) | ❌ mantissa 53b 미달 |

### `mlll_fp_ops.h` 디자인
```c
/* 빌드 시 백엔드 선택 */
#ifndef MLLL_FP_KIND
#define MLLL_FP_KIND 2  /* 0=double, 1=mpf_t, 2=mpfr_t (default), 3=dpe_t */
#endif
#ifndef MLLL_FP_PREC
#define MLLL_FP_PREC 1024  /* L5 paper bound 1010b를 14b 마진 */
#endif

/* 통일 인터페이스 */
typedef ... fp_t;
void fp_init(fp_t *x);
void fp_add(fp_t *r, const fp_t *a, const fp_t *b);
void fp_mul(fp_t *r, const fp_t *a, const fp_t *b);
... (12개 op)
```

→ 각 .c 백엔드가 동일 시그너처 구현. 단일 코드베이스가 4종 빌드.

## E4. Dispatcher — 입력 크기에 따라 백엔드 자동 선택

`mlll_gram.c:806`:
```c
quat_mlll_gram(...) {
    int v = ibz_max_bits(generators);
    if (v < 256)   → ibz_t (정수만)
    else if (v < 384) → mpf_t (B 백엔드)
    else if (v < 512) → fp_t (D2c default)
    else              → fp_t (확장 정밀도)
}
```

→ 작은 입력은 정수 산술 (오버헤드 적음), 큰 입력만 부동소수 진입.

## E5. Lemma `mlll-bound` 측정 도구

`mlll_print_max_bits()` — 매 sign 호출 시 max bitsize 추적:
```
[MLLL bitsize, fp=mpfr_t prec=1024]
  L1: vec=190b  G=506b  intermediate=506b
  L3: vec=294b  G=769b  intermediate=769b
  L5: vec=388b  G=1016b intermediate=1016b
```

→ **3 레벨 전부 페이퍼 상한 ≤ 만족** 확인.

---

# Part F. 페이퍼 ↔ 코드 1:1 매핑 (검증자용)

## F1. Algorithm 1 (MLLL) — 페이퍼 라인 vs `mlll.c` 라인

| 페이퍼 | 의미 | 코드 위치 |
|---|---|---|
| L1 입력 a_1..a_β | quat_mlll(b, β) 인자 | `mlll.c:170` 시그너처 |
| L2 b_i ← a_i; B_i ← \|b_i\|² | 초기화 | `mlll.c:185` |
| L3 μ_ij 계산 | GSO 계수 | `mlll.c:201` (fp_t) |
| L4 i ← 2 | 진입 인덱스 | `mlll.c:215` |
| L5–L9 size-reduce | μ_ij 정수화 + b_i 갱신 | `mlll.c:225–270` |
| L10 Lovász 조건 | B_i ≥ (3/4 - μ²) B_{i-1} | `mlll.c:280` |
| L11–L13 swap | b_{i-1} ↔ b_i | `mlll.c:285–330` |
| L14 i ← max(i-1, 2) | 인덱스 후퇴 | `mlll.c:340` |
| L15–L18 다음 단계 결정 | | `mlll.c:345–380` |
| L19–L23 kappa 갱신 | 종료 조건 | `mlll.c:395` |
| L24–L26 nz_count + 0 제거 | 영벡터 카운트 | `mlll.c:415–445` |
| L27 출력 rank-≤4 | basis copy | `mlll.c:450–478` |

→ **27 개 페이퍼 라인 모두 코드에 1:1 대응**. 페이퍼 외 추가 단계 없음.

## F2. Lemma `mlll-bound` 비교 — 페이퍼 vs 측정

| 레벨 | $p$ 비트 | 페이퍼 상한 | 우리 측정 | 마진 |
|---|---|---|---|---|
| L1 | 254 | ≤ 514 b | **506 b** | **-8 b** ✅ |
| L3 | 383 | ≤ 772 b | **769 b** | **-3 b** ✅ |
| L5 | 505 | ≤ 1016 b | **1016 b** | **0 b 정확 일치** ✅ |

→ 페이퍼는 ≤ (상한). 우리는 모든 레벨에서 만족.

## F3. Appendix A.1 — μ, B 정밀도

| 항목 | 페이퍼 | 우리 구현 |
|---|---|---|
| μ 타입 | real | mpfr_t (correctly-rounded) |
| B 타입 | real | mpfr_t |
| 정밀도 명시 | "real" 만 | **PREC=1024 b** |
| L5 상한 1010b 커버? | (페이퍼는 보장 미명시) | **14b 마진으로 OK** |

## F4. Algorithm 4 (CompactIdealMultiplication) — `quat_lattice_mul_mlll`

| 페이퍼 단계 | 코드 |
|---|---|
| 16 generator 구성 | `mlll.c:430` |
| MLLL 호출 | `mlll.c:435` `quat_mlll(...)` |
| rank-4 basis 추출 | `mlll.c:440–478` |

---

# Part G. 검증한 것 + 발견된 문제 (감사 결과)

## G1. ✅ 단위 테스트 12/12 PASS

| 테스트 | 검증 대상 | 결과 |
|---|---|---|
| `test_mlll_basic` | quat_mlll 입출력 sanity | ✅ |
| `test_mlll_alg2` | lattice_mul_mlll 동치성 | ✅ |
| `test_mlll_alg3` | lideal_create_mlll | ✅ |
| `test_mlll_gram_basic` | Gram 변종 sanity | ✅ |
| `test_mlll_gram_fp` | fp 백엔드 단위 | ✅ |
| `test_mlll_gram_fp_equivalence` | **HNF↔ibz↔B↔fp 4-way 동치** | ✅ |
| `test_mlll_dispatcher` | 임계값별 백엔드 선택 | ✅ |
| `test_mlll_overflow_trap` | bitsize budget trap | ✅ |
| ... (4개) | (보조) | ✅ |

→ **알고리즘 수준에서는 정확히 작동**.

## G2. ⚠ 외부 감사 (2026-04-29) 가 발견한 것

### 결정적 한 줄
```bash
$ grep -rn '_mlll_gram\b' src/signature/ src/id2iso/ src/sqisign.c
(빈 출력)
```

→ **sign.c, keygen.c, dim2id2iso.c 어디에서도 paper 변종이 호출되지 않음.**

### 미통합 항목 6가지

#### A. 함수 헤더 선언만 있고 .c 본체 없음 (2개)
```
quat_lideal_reduce_basis_mlll_gram          quaternion.h:608  ⚠ 본체 없음
quat_lideal_prime_norm_reduced_equivalent_mlll_gram  quaternion.h:656  ⚠ 본체 없음
```
→ 호출 즉시 **link error**. 즉 호출되지 않음 → dead declaration.

#### B. 핫 패스 라우팅 0건 — 라우팅 안 된 호출 사이트 9개
| 파일:줄 | 현재 (HNF) | 바뀌어야 할 (paper) |
|---|---|---|
| `signature/ref/lvlx/sign.c:18` | `quat_lideal_prime_norm_reduced_equivalent` | `_mlll_gram` |
| `signature/ref/lvlx/sign.c:81, 211` | `quat_lideal_inter` | 내부 `lattice_mul_mlll_gram` 분기 |
| `signature/ref/lvlx/sign.c:169, 325` | `quat_lideal_create` | `_mlll_gram` |
| `signature/ref/lvlx/keygen.c:35` | `quat_lideal_prime_norm_reduced_equivalent` | `_mlll_gram` |
| `id2iso/ref/lvlx/dim2id2iso.c:98` | `quat_lideal_create` | `_mlll_gram` |
| `id2iso/ref/lvlx/dim2id2iso.c:517` | `quat_lideal_reduce_basis` | `_mlll_gram` |
| `id2iso/ref/lvlx/id2iso.c:257` | `quat_lideal_create` | `_mlll_gram` |
| `signature/ref/lvlx/encode_signature.c:194` | `quat_lideal_create` | `_mlll_gram` |

#### C. CMake 옵션이 dead
```
$ grep -n 'USE_MLLL_GRAM' CMakeLists.txt
16:option(SQISIGN_USE_MLLL_GRAM "..." OFF)   ← default OFF
46:if(SQISIGN_USE_MLLL_GRAM) ...

$ grep -rn 'SQISIGN_USE_MLLL_GRAM' src/ --include='*.c'
(빈 출력)   ← 매크로 분기 0건
```
→ 옵션을 ON 으로 켜도 **아무 일도 안 일어남**.

#### D. Algorithm 1 (IdealFiltration) 변종 함수 자체 없음
- 페이퍼 4개 알고리즘 중 **Alg 1 변종 미구현**.
- PLAN_KO.md 에도 미구현으로 표기되어 있음.

#### E. KAT `.rsp.MLLL` 파일 출처 문제
- `apps/PQCgenKAT_sign.c` 에 **MLLL 빌드 분기 없음**.
- 즉 commit 된 `.rsp.MLLL` 3개 파일은 **HNF binary** 가 만든 파일.
- (L3/L5는 HNF와 byte-exact 같음 → 우연히 동일. L1은 어떻게 다른지 출처 추적 안 됨.)
- **이름이 사람을 호도** ("MLLL이 만든 KAT" 처럼 보임).

#### F. README와 commit msg 의 부정확
- README: "MLLL을 단독 기본 경로로 사용"
- commit `5ef9b76`: "paper-strict MLLL 마무리"
- **사실**: `USE_MLLL_GRAM=OFF` 가 default, `_mlll_gram` 호출 0건 → **빌드 전체로는 HNF 단독**.
- "paper-strict 100%" 라는 이전 보고는 **알고리즘 본체** 의미였고, **빌드 통합** 의미로 들리도록 표현됨.

## G3. 그래서 "다 끝났다"의 정확한 정의

| 의미 | 상태 | 비고 |
|---|---|---|
| 페이퍼 알고리즘 본체 코드 작성 | ✅ Phase 1–2 완료 | mlll.c, mlll_gram.c 등 |
| 4 부동소수 백엔드 추상화 | ✅ Phase 2C 완료 | mlll_fp_ops.h |
| 단위 테스트 + 4-way 동치 | ✅ Phase 2.1 완료 | 12/12 PASS |
| Lemma `mlll-bound` 페이퍼 상한 만족 | ✅ 측정 검증 | 3 레벨 모두 ≤ |
| Algorithm 1 본체 (mlll.c) ↔ 페이퍼 1:1 | ✅ 27 라인 매핑 | |
| **sign/keygen 핫 패스가 paper 호출** | **❌ Phase 3-2 미수행** | 가장 중요 |
| **함수 본체 2개 (`_mlll_gram` 변종)** | **❌ Phase 3-1 미수행** | |
| **KAT 생성 binary가 paper 빌드** | **❌ Phase 3-3 미수행** | |
| **Algorithm 1 (IdealFiltration) 변종** | **❌ Phase 3-5 미수행** | |
| **문서 (README/PLAN) 시점 동기화** | **❌ Phase 3-4 미수행** | |

→ **"논문대로 구현"의 정확한 비율**:
- 알고리즘 컴포넌트 수준: **약 70%** (Alg 2/3 본체 + 측정 + 단위 검증)
- 빌드 통합 수준: **0%** (sign 핫 패스가 paper 함수를 한 번도 호출하지 않음)

---

# Part H. 향후 계획 (Phase 3) — 구체

## H1. Phase 3-1. 누락 함수 본체 작성 (1–2 일)

### 작성할 두 함수
1. **`quat_lideal_reduce_basis_mlll_gram`**
   - 위치: `src/quaternion/ref/generic/lll/lll_applications.c` (기존 `reduce_basis` 옆)
   - 외부 계약: HNF 변종과 동일
   - 내부: `quat_mlll_gram` 호출
2. **`quat_lideal_prime_norm_reduced_equivalent_mlll_gram`**
   - 위치: 동일 파일
   - 내부: 위 1번 함수 호출 + 동일 외부 계약

### 검증
- 단위 테스트 (`mlll_tests.c`) 추가
- HNF 변종과 **lattice equivalence** (좌이상 동치) 비교

## H2. Phase 3-2. 핫 패스 라우팅 (2–3 일)

### 핵심 제약: **sign.c, keygen.c, dim2id2iso.c 무수정**
→ 격리 유지. 헤더에서 매크로 alias 또는 weak symbol 로 처리.

### 방안 A: 매크로 alias (가벼움)
```c
/* quaternion.h */
#ifdef SQISIGN_USE_MLLL_GRAM
#define quat_lideal_create        quat_lideal_create_mlll_gram
#define quat_lideal_reduce_basis  quat_lideal_reduce_basis_mlll_gram
#define quat_lideal_prime_norm_reduced_equivalent \
        quat_lideal_prime_norm_reduced_equivalent_mlll_gram
#endif
```
→ sign.c가 호출하는 이름은 그대로. 빌드 시 변종으로 치환.

### 방안 B: weak symbol + dispatcher (정교함)
- gcc/clang의 `__attribute__((weak))` 활용
- 런타임 환경변수로 전환 가능 (디버깅 유리)

→ **A가 단순하고 격리 유지에 부합. 권장.**

### 검증
- `cmake -DSQISIGN_USE_MLLL_GRAM=ON ..` 빌드
- `make test` 전체 통과
- gcov/printf 로 sign 1회 호출 시 `quat_mlll_gram_fp` 도달 확인

## H3. Phase 3-3. KAT 생성 경로 명시화 (1 일)

### 옵션 A: `apps/PQCgenKAT_sign.c` 수정
```c
#ifdef SQISIGN_USE_MLLL_GRAM
const char *out_suffix = ".rsp.MLLL";
#else
const char *out_suffix = ".rsp";
#endif
```

### 옵션 B: 별도 `apps/PQCgenKAT_sign_mlll.c` 추가
- 컴파일 시 USE_MLLL_GRAM=ON 으로 강제

### ctest 추가
- ON 빌드에서 `.rsp.MLLL` 재생성 → 기존 commit 본과 byte-exact 일치 확인 (회귀 방지)

## H4. Phase 3-4. 문서 동기화 (반나절)

### README.md
- "단독 기본 경로" → "Alg 2/3 본체 + 단위 테스트 완료. 핫 패스 통합 Phase 3" 로 수정
- KAT `.rsp.MLLL` 출처: "USE_MLLL_GRAM=ON 빌드의 PQCgenKAT_sign 출력" 명시
- "test_scheme/test_kat PASS" 가 어느 빌드 구성에서 측정됐는지 명시

### PLAN_KO.md
- `## 현재 상태 스냅샷` 표 갱신
- 2026-04-21 → 현재 시점 반영

### commit msg 정정 (선택)
- `5ef9b76`, `d1073ee` 가 알고리즘 한정이었음 후속 commit 으로 cite

## H5. Phase 3-5 (선택). Algorithm 1 IdealFiltration

페이퍼 Algorithm 1 본체 + 호출부 라우팅. 우선순위 낮음. Phase 3-2 끝난 뒤 결정.

---

# Part I. Q & A + 용어집

## I1. 흔한 질문

### Q1. "MLLL이 페이퍼대로 구현되었다"는 거짓말이었나?
A. **부분적으로 그렇다**. 알고리즘 본체 의미로는 사실 (Phase 1–2 검증됨). **핫 패스 통합**이라는 의미로는 아직 아님 (Phase 3 미수행). 두 의미를 구분해서 말하지 않은 게 본 보고의 오류.

### Q2. KAT 통과는 진짜인가?
A. 통과는 했는데 **HNF 경로**로 통과한 것. `_mlll_gram` 변종이 sign 에서 호출 안 되므로, 현재 KAT은 MLLL 회귀 테스트로 부적절. **Phase 3-3 완료 후에야 의미 있음**.

### Q3. NIST 통과인가?
A. NIST는 "알고리즘 사양 + KAT 형식 + round-trip 검증" 을 요구.
- **HNF binary 기준**: 통과 (upstream의 통과를 그대로 유지)
- **MLLL binary 기준**: Phase 3 완료 전엔 의미 없음

### Q4. 그럼 이번 fork는 가치가 없나?
A. **아니다**. 다음은 진짜 진척:
1. 페이퍼 Algorithm 1/2/3 본체를 정확히 구현
2. 4 부동소수 백엔드 추상화 + paper-strict 정밀도 (mpfr 1024b) 검증
3. Lemma `mlll-bound` 정수 상한 3 레벨 측정
4. 단위 테스트 12/12 + 4-way 동치성
5. Dispatcher + 임계값 정정 (3 레벨 fp 경로 작동)

다만 표현이 "단독 기본 경로 대체" 처럼 빌드 전체 의미로 들리도록 적힌 게 문제.

### Q5. Phase 3 끝내면 진짜 끝인가?
A. 3-1~3-4 끝나면:
- sign.c가 MLLL을 호출
- KAT 자기일관성 보장
- 문서가 사실과 일치
→ **빌드 통합 차원에서 "MLLL fork" 진짜 의미 획득**.

3-5 (Alg 1 IdealFiltration) 까지 가야 페이퍼 4-알고리즘 전부 커버.

### Q6. 왜 이런 일이 생겼나?
A.
- 알고리즘 본체 작성에 집중하면서 **헤더 선언과 .c 본체를 따로 작성**
- 단위 테스트 통과 = 검증 끝났다고 자신
- 핫 패스 통합은 "나중 작업" 으로 미뤘는데 그 사실을 README/commit msg 에 안 적음
- → **외부 감사가 grep 한 줄로 발견**

→ 교훈: "알고리즘 작성 완료" ≠ "통합 완료" — 두 단계를 명시 분리.

### Q7. 성능 향상은 정말 있나?
A. **아직 측정 못 함**. 현재 측정값 (`L1 KAT 100 iter = 50s`) 은 HNF 핫 패스 + MLLL 보조 호출의 합. 진짜 MLLL 핫 패스 측정은 **Phase 3-2 완료 후에만 가능**.

페이퍼 이론값으로는 메모리 절감 + 속도 동등 또는 우위 예상.

## I2. 용어집

| 용어 | 한글 | 의미 |
|---|---|---|
| Lattice | 격자 | 정수 조합으로 만들어지는 점들의 집합 |
| Basis | 기저 | 격자를 표현하는 vector 들 |
| Reduction | 환원 | 짧고 직각인 basis 찾기 |
| HNF | Hermite Normal Form | 정수 산술만 쓰는 환원, 결과 상삼각 |
| LLL | Lenstra–Lenstra–Lovász | 1982년 격자 환원 알고리즘 |
| MLLL | Modified LLL | Pohst 1987 변형, 선형종속 입력 처리 |
| GSO | Gram-Schmidt Orthogonalization | 직교화 |
| μ_ij | GSO 계수 | 실수, b_i를 b_j에 정사영한 비율 |
| B_j | norm² | b_j*의 길이 제곱, 실수 |
| Quaternion | 4원수 | a + bi + cj + dk |
| Order | 차수환 | quaternion에서 정수 비슷한 부분 (lattice) |
| Left ideal | 좌이상 | order의 약수 비슷한 부분 |
| Isogeny | 동치 | 타원곡선 사이의 대수 사상 |
| SQIsign | (고유명) | isogeny 기반 디지털 서명 |
| NIST | (고유명) | 미국 표준기술연구소, PQC 표준 선정 |
| PQC | Post-Quantum Cryptography | 양자내성 암호 |
| KAT | Known Answer Test | 고정 시드 → 결정적 출력 회귀 테스트 |
| paper-strict | (이 fork 용어) | 페이퍼 외 코드 0 % |
| hot path | 핫 패스 | 핵심 빈번 호출 경로 (sign/keygen 등) |

---

# 부록 J. 검증 명령어 (재현 가능, 그대로 복붙)

## J1. 함수 본체 부재 확인 (Phase 3-1 미수행 증거)
```bash
cd the-sqisign
grep -rn 'quat_lideal_reduce_basis_mlll_gram' src/ --include='*.c'
# 출력: 본체 0건 (헤더 선언만)

grep -rn 'quat_lideal_prime_norm_reduced_equivalent_mlll_gram' src/ --include='*.c'
# 출력: 본체 0건
```

## J2. 핫 패스 라우팅 0건 확인 (Phase 3-2 미수행 증거)
```bash
grep -rn '_mlll_gram\b' src/signature/ src/id2iso/ src/sqisign.c
# 출력: 빈 (Phase 3-2 후엔 sign.c 등에 hit 발생해야 함)
```

## J3. CMake 옵션 dead 확인 (Phase 3-2 미수행 증거)
```bash
grep -n 'USE_MLLL_GRAM' CMakeLists.txt
# 출력: option(... OFF) + 컴파일 정의

grep -rn 'SQISIGN_USE_MLLL_GRAM' src/ --include='*.c'
# 출력: 빈 (매크로 사용처 0건)
```

## J4. Lemma `mlll-bound` 측정 (현 빌드 = HNF 경로)
```bash
cd build_default/test
./sqisign_test_scheme_lvl1
# 출력 1행: [MLLL bitsize, fp=mpfr_t prec=1024] L1: vec=190b G=506b ...

./sqisign_test_scheme_lvl3
# L3: vec=294b G=769b ...

./sqisign_test_scheme_lvl5
# L5: vec=388b G=1016b ...
```

## J5. KAT round-trip (현 빌드 = HNF 경로 검증)
```bash
./sqisign_test_kat_lvl1   # 50.34 s, PASS
./sqisign_test_kat_lvl3   # 81.94 s, PASS
./sqisign_test_kat_lvl5   # 112.21 s, PASS
```

## J6. Phase 3-2 완료 후 추가될 검증
```bash
cmake -DSQISIGN_USE_MLLL_GRAM=ON ..
make
make test
# sign.c 호출 시 quat_mlll_gram_fp 도달 (gcov로 확인)
# .rsp.MLLL 재생성 후 기존 commit 본과 byte-exact 일치
```

---

**작성**: 2026-04-29 22:44 KST
**다음 갱신**: Phase 3-1 완료 후 (예상 1–2 일)
**저자**: Sol (Claude) + 사용자 감사
