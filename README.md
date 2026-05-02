# SQIsign — Compact paper MLLL fork

> ⚠ **거짓 표현 정정 (2026-04-29)**: 이전 README는 "MLLL을 단독 기본 경로로 사용" 이라고
> 적혀 있었으나 외부 감사로 **거짓**임이 발각. 실제로는 sign/keygen 핫 패스가 여전히 HNF만
> 호출. 자세한 사고 기록은 `PRESENTATION_KO.md` 와 `CLAUDE.md` 참고.

## 현재 작업 계획 (2026-05-02 갱신)

**현 위치**: Phase 3 (3-1, 3-2′, 3-3, 3-5) **완료**. P3-4 (벤치 모드 `--mode=alg4`) + Phase 4 (fp budget protocol-level analysis) 대기.

**다음 행동 (우선순위 순)**:
1. ~~`quat_lideal_reduce_basis_mlll_gram` + `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` **본체 작성**~~ → ✅ Phase 3-1 완료 (`lll_applications.c` +94/+98 라인, 단위 테스트 17/17 PASS)
2. ~~`SQISIGN_USE_MLLL_GRAM=ON` 시 매크로 alias 로 sign.c/keygen.c/id2iso/dim2id2iso/encode_signature **무수정 라우팅**~~ → ✅ Phase 3-2′ 완료 (`quaternion.h:777-791` 4개 함수 alias, 11개 호출부, CI `mlll-routing-check.yml` 양쪽 빌드 PASS)
3. ~~e2e 동치성 (HNF↔MLLL_GRAM 양쪽 빌드 self-consistency)~~ → ✅ Phase 3-5 완료 (signature / nistapi / threadsafety / SELFTEST / id2iso lvl1/3/5 PASS, commit `7e8f0eb` 기준)
4. **P3-4** 벤치 모드 `--mode=alg4` 추가 (alg2/alg3 와 통일) ← 다음 단계 (작은 작업)
5. **Phase 4** fp 백엔드 폭을 SQIsign protocol bound 기반으로 산정 (현재 ON 빌드는 `g_fp_mode=0` ibz fallback. 자세한 finding: `src/quaternion/ref/generic/lll/PLAN_KO.md` "새 research finding (2026-04-30)" 절)
6. **Phase 5** PR 분리 / upstream merge

> **2026-04-30 정정**: 이전 판본의 5번 항목 "Algorithm 1 IdealFiltration 변종"은 **폐기**. 논문 *Compact Quaternion Algorithms for SQIsign*은 Alg 1=MLLL 커널 / Alg 2=CompactIdealMultiplication / Alg 3=RandomIdealGivenPrimeNorm / Alg 4=RandomEquivalentPrimeIdeal 4개로 구성되며 "IdealFiltration"이라는 algorithm은 paper LaTeX 어디에도 등장하지 않음. 자세한 내용은 `src/quaternion/ref/generic/lll/PLAN_KO.md` 상단 "2026-04-30 정정" 절 참조.

**Phase 3-2′ 검증 (alias 방식이라 source-grep 부적절. 셋 중 하나)**:

```bash
# (1) preprocessed output (가장 직접적)
gcc -E -I... build_mlll/.../sign.c | grep -c 'quat_lideal_create_mlll_gram'  # 1+건

# (2) symbol-level routing (오브젝트 파일)
nm build_mlll/src/signature/.../sign.c.o | grep -cE '_mlll_gram$'  # 1+건

# (3) 양쪽 빌드 ctest PASS — CI mlll-routing-check.yml 가 자동 수행
```

## 현재 상태 — 완료 vs 미완료

| 항목 | 상태 | 검증 명령 |
|---|---|---|
| Algorithm 1/2/3 본체 작성 | ✅ | `wc -l src/quaternion/ref/generic/lll/mlll*.c` |
| 4 부동소수 백엔드 (D2a–D2d) | ✅ | `cat src/quaternion/ref/generic/lll/mlll_fp_ops.h` |
| 단위 테스트 17/17 PASS (Phase 3-1 추가 2개 포함) | ✅ | `./build_default/src/quaternion/ref/generic/test/sqisign_test_mlll` |
| 4-way 동치성 (HNF/ibz/B/fp) | ✅ | `test_mlll_gram_fp_equivalence` |
| Lemma `mlll-bound` 페이퍼 상한 만족 | ✅ | L1=506b ≤ 514b, L3=769b ≤ 772b, L5=1016b = 1016b |
| 함수 본체 4개 (`reduce_basis_mlll_gram`, `prime_norm_reduced_equivalent_mlll_gram`, `lideal_create_mlll_gram`, `lideal_lideal_mul_reduced_mlll_gram`) | ✅ | `grep -n '_mlll_gram\b' src/quaternion/ref/generic/lll/lll_applications.c src/quaternion/ref/generic/lll/mlll_gram.c` |
| **MLLL_GRAM 핫 패스 라우팅 (alias)** | ✅ | `quaternion.h:777-791` 4개 함수 alias. CI `mlll-routing-check.yml` 양쪽 빌드 PASS |
| **`SQISIGN_USE_MLLL_GRAM` 매크로 사용처** | ✅ | `grep -rn 'SQISIGN_USE_MLLL_GRAM' src/ --include='*.c' --include='*.h'` 다수 (헤더 alias + ideal.c/lattice.c/mlll_benchmark.c 등 IMPL 가드) |
| **e2e 동치성 (양쪽 빌드 self-consistency)** | ✅ | ON 빌드의 signature / nistapi / threadsafety / SELFTEST / id2iso lvl1/3/5 PASS (commit `7e8f0eb`) |
| **fp 백엔드 (production lattice)** | 🟡 ibz fallback | ON 빌드는 `g_fp_mode=0` (ibz_t). fp는 단위 테스트/벤치 한정. Phase 4 task. |
| **KAT `.rsp.MLLL` 출처 = MLLL binary** | 🟡 | alias 방식이라 ON 빌드 `PQCgenKAT_sign` 자체가 자동 라우팅됨. 단 `KAT/*.rsp.MLLL` 3개 파일은 history 보관용 (Phase 3 이전 시점 출처). 회귀 검증에 쓰려면 ON 빌드로 재생성 필요. |
| **README/PLAN/commit msg 동기화** | ✅ | 2026-05-02 본 정정으로 동기화 |

→ **알고리즘 컴포넌트 수준: 100% 완료** (4개 paper algorithm 모두 MLLL_GRAM 변종 본체 + 단위 테스트 PASS)
→ **빌드 통합 수준: 매크로 alias 라우팅 + 양쪽 빌드 e2e PASS — fp 백엔드는 ibz fallback 상태**

## 마지막 검증 시점

- **2026-04-27 22:18 KST** (HNF 시점 baseline):
  - `build_default` (= `MLLL_FP_KIND=2 MLLL_FP_PREC=1024`, **`SQISIGN_USE_MLLL_GRAM=OFF`**)
  - `sqisign_test_scheme_lvl{1,3,5}` PASS
  - `sqisign_test_kat_lvl{1,3,5}` 100 iter PASS (50.34s / 81.94s / 112.21s)
  - 이 측정값은 HNF 핫 패스 + MLLL 보조 호출의 합.
- **2026-04-29 22:38 KST**: 외부 감사로 핫 패스 미통합 발견.
- **2026-04-29 23:0x KST (Phase 3-1)**: `quat_lideal_reduce_basis_mlll_gram` + `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` 본체 작성. `sqisign_test_mlll` 17/17 PASS.
- **2026-04-30 (Phase 3-2′ + 3-5, commits `2e4d178` … `c3e7b37`)**:
  - 매크로 alias 라우팅 (`quaternion.h:777-791`) + `sqisign_namespace.h` 등록 + `SQISIGN_MLLL_GRAM_IMPL` 가드 적용.
  - CI workflow `mlll-routing-check.yml` 신설 (matrix={OFF,ON}, ubuntu-latest, `nm` symbol-level 검증).
  - **양쪽 빌드 (`SQISIGN_USE_MLLL_GRAM=OFF`/`ON`) ctest 전 통과** — `signature` / `nistapi` / `threadsafety` / `SELFTEST` / `id2iso` lvl1/3/5 / `mlll` 단위 테스트.
  - **fp 백엔드 finding**: ON 빌드에서 production lattice gram bit-size가 749/1137/1507 (vs Phase 1 측정 budget 518/782/1026)로 fp 폭 부족 확인. 임시 결정으로 `g_fp_mode=0` (ibz_t fallback) default. fp 폭은 SQIsign protocol bound 기반 산정이 필요 (Phase 4).
- **2026-05-02 (본 정정)**: HEAD = `e82d86e`. README ↔ 코드 ↔ PLAN_KO.md 동기화 commit. 신규 검증 측정 없음 (e82d86e 시점 그대로 보고).

## 알려진 미통합/거짓 위험 항목

1. **이전 README 거짓 주장의 흔적**: 2026-04-29 이전 commit msg/README의 "MLLL 단독 기본 경로" 표현은 **알고리즘 한정 의미였음**. 본 정정으로 빌드 통합 상태가 분리 표기됨.
2. **KAT 파일명 호도**: `KAT/PQCsignKAT_*.rsp.MLLL` 3개 파일은 **Phase 3-2′ 통합 이전** 시점에 생성됨. alias 방식 도입 후로는 ON 빌드 `PQCgenKAT_sign_lvlN` 이 자동으로 MLLL 라우팅된 binary가 되지만, **현재 repo의 .rsp.MLLL 파일 자체는 그 binary로 재생성된 것이 아님**. 회귀 검증 기준선으로 쓰려면 ON 빌드로 재생성 후 commit.
3. **Lemma bound 측정값** (vec 259/391/513, gram 518/782/1026)은 OFF 빌드의 보조 호출에서 측정된 값. ON 빌드 production 호출 사슬에서는 더 큼 (749/1137/1507). MLLL 내부 settled state ≠ MLLL 입력 generator size — Phase 4 protocol-level analysis 대상.
4. **`g_fp_mode=0` default under `USE_MLLL_GRAM=ON`**: 페이퍼 contribution "compact 연산이 고정폭에 들어간다"의 빌드 통합 증명은 **production input 입력 폭 산정 후 재활성** 필요. 단위 테스트/벤치 (입력 bound 통제됨) 환경에서는 fp 경로가 여전히 검증되어 있음.
5. **commit `5ef9b76` / `d1073ee`** 의 표현은 알고리즘 한정 의미였음 — 본 README 표가 범위를 명시.

## Compact paper 매핑 (간단 요약)

| 페이퍼 | 구현 위치 | 상태 |
|---|---|---|
| Algorithm 1 MLLL 커널 | `src/quaternion/ref/generic/lll/mlll.c:170` `quat_mlll`, `mlll_gram.c` `quat_mlll_gram` | ✅ |
| Algorithm 2 CompactIdealMultiplication | `mlll.c:399` `quat_lattice_mul_mlll` (16-generator), `mlll_gram.c` `quat_lattice_mul_mlll_gram` | ✅ |
| Algorithm 3 RandomIdealGivenPrimeNorm | `mlll.c:451` `quat_lattice_add_mlll`, `mlll_gram.c:880` `quat_lideal_create_mlll_gram` | ✅ |
| Algorithm 4 RandomEquivalentPrimeIdeal | `lll_applications.c:213` `quat_lideal_prime_norm_reduced_equivalent_mlll_gram` | ✅ (Phase 3-1, 단위 테스트 5/5 PASS) |
| Lemma `mlll-bound` (정수 상한) | 측정값 ≤ 페이퍼 상한 | ✅ 3 레벨 (단, 측정 환경 = OFF 빌드 보조 호출. 본 README "알려진 미통합/거짓 위험 항목" 3 참조) |
| Appendix A.1 (μ/B real) | mpfr_t PREC=1024 | ✅ (단위 테스트/벤치 한정. ON 빌드 production 입력 폭 부족 — Phase 4 task) |
| **sign/keygen 핫 패스 통합 (정수 백엔드)** | `quaternion.h` 매크로 alias | ✅ 양쪽 빌드 e2e PASS |
| **sign/keygen 핫 패스 통합 (fp 백엔드)** | — | 🟡 ibz fallback (`g_fp_mode=0`), Phase 4 |

## 빌드 옵션

```bash
# 기본 빌드 (HNF 핫 패스 + MLLL 보조 호출)
cmake -DSQISIGN_BUILD_TYPE=ref ..

# MLLL_GRAM 핫 패스 라우팅 (alias 방식). 양쪽 빌드 e2e PASS.
# 현재 백엔드는 ibz fallback (`g_fp_mode=0`) — fp 백엔드 폭 산정은 Phase 4.
cmake -DSQISIGN_BUILD_TYPE=ref -DSQISIGN_USE_MLLL_GRAM=ON ..

# MLLL 백엔드/정밀도 (실험용 — 단위 테스트에서만 작동)
cmake -DMLLL_FP_KIND=2 -DMLLL_FP_PREC=1024 ..
#   MLLL_FP_KIND: 0=double, 1=mpf_t, 2=mpfr_t (default), 3=dpe_t
#   ⚠ KIND=0/3 (53 b mantissa)은 페이퍼 상한 (L1=506b 등) 수학적 미달
```

## 자세한 보고서

- `PRESENTATION_KO.md` — 사전 지식 0 사람 위한 11-Part 전체판 (735 줄)
- `src/quaternion/ref/generic/lll/PLAN_KO.md` — phase별 작업 계획 (시점 갱신 필요)
- `src/quaternion/ref/generic/lll/AUDIT_2026-04-19_KO.md` — 4월 19일 감사
- `CLAUDE.md` — 이 저장소 작업 시 Claude 세션 규칙

---

This library is a C implementation of SQIsign.

## Requirements

- CMake (version 3.13 or later)
- C11-compatible compiler
- GMP (version 6.0.0 or later)

### Pre-computation

The constant values in the `src/precomp` directory were generated using the
pre-computation scripts in the `scripts/precomp` directory. It is not necessary
to execute these scripts to compile the project. The scripts have the following
requirements:
- [two-isogenies](https://github.com/ThetaIsogenies/two-isogenies)
  (`Theta-SageMath` version).
- [deuring-2D](https://github.com/Jonathke/deuring-2D)

## Build

For a generic build
```
$ mkdir -p build
$ cd build
$ cmake -DSQISIGN_BUILD_TYPE=ref ..
$ make
$ make test
```

An optimized executable with debug code and assertions disabled can be built
replacing the `cmake` command above by
```
cmake -DSQISIGN_BUILD_TYPE=<ref/broadwell> -DCMAKE_BUILD_TYPE=Release ..
```

## Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

### SQISIGN_BUILD_TYPE

Specifies the build type for which SQIsign is built. The currently supported values are:
- `ref`: builds the plain reference implementation.
- `opt`: builds the optimized implementation which is the same as the reference
  implementation.
- `broadwell`: builds an additional optimized implementation targeting the Intel
  Broadwell architecture (and later). The optimizations are applied to the
  finite field arithmetic.

### GMP_LIBRARY

If set to `SYSTEM` (by default), the gmp library on the system is dynamically linked.

If set to `BUILD`, a custom gmp library is linked, which is built as part of the overall build process.
In this case, the following further options are available:
- `ENABLE_GMP_STATIC`: Does static linking against gmp. The default is `OFF`.
- `GMP_BUILD_CONFIG_ARGS`: Provides additional config arguments for the gmp build (for example `--disable-assembly`). By default, no config arguments are provided.

If set to `MINI`, the mini-gmp library is used, whose sources are included in the repository, in the folder `src/mini-gmp`. In this case, no copies of the full gmp library (system or custom-built) are required.

### ENABLE_SIGN

If set to `ON` (default), SQIsign is built with signature and verification functionality.
If set to `OFF`, SQIsign is built with verification functionality only.
In the latter case, GMP is no longer a dependency.

### CMAKE_BUILD_TYPE

Can be used to specify special build types. The options are:

- `Release`: Builds with optimizations enabled and assertions disabled.
- `Debug`: Builds with debug symbols.
- `ASAN`: Builds with AddressSanitizer memory error detector.
- `MSAN`: Builds with MemorySanitizer detector for uninitialized reads.
- `LSAN`: Builds with LeakSanitizer for run-time memory leak detection.
- `UBSAN`: Builds with UndefinedBehaviorSanitizer for undefined behavior detection.

The default build type uses the flags `-O3 -Wstrict-prototypes -Wno-error=strict-prototypes -fvisibility=hidden -Wno-error=implicit-function-declaration -Wno-error=attributes`. (Notice that assertions remain enabled in this configuration, which harms performance.)

## Test

In the build directory, run `make test` or `ctest`.

The test harness consists of the following units:

- KAT test: `SQIsign_<level>_KAT`- tests against the KAT files in the `KAT`
  directory.
- Self-tests: `SQIsign_<level>_SELFTEST` - runs random self-tests
  (key generation, signature and verification).
- Sub-library specific unit-tests.

Note that, `ctest` has a default timeout of 1500s, which is applied to all tests
except the KAT tests. To override the default timeout, run
`ctest --timeout <seconds>`.

## Known Answer Tests (KAT)

KAT are available in the `KAT` directory. They can be generated by running the
apps built in the `apps` directory:
```
apps/PQCgenKAT_sign_<level>
```

A successful execution will generate the `.req` and `.rsp` files.

A full KAT test is done as part of the test harness (see the [Test](#test)
section).

## Benchmarks

A benchmarking suite is built and can be executed with the following command:
```
apps/benchmark_<level> [--iterations=<iterations>]
```
where `<level>` specifies the SQIsign parameter set and `<iterations>` is the
number of iterations used for benchmarking; if the `--iterations` option is
omitted, a default of 10 iterations is used.

The benchmarks profile the key generation, signature and verification functions. The results are reported in CPU cycles if available on the host platform, and timing in nanoseconds otherwise.

## Examples

Example code that demonstrates how to use SQIsign with the NIST API is available
in `apps/example_nistapi.c`.

## Project Structure

The source code consists of a number of sub-libraries used to implement the
final SQIsign library:
- `common`: common code for hash function, seed expansion, PRNG, memory handling.
- `mp`: code for saturated-representation multiprecision arithmetic.
- `gf`: GF(p^2) and GF(p) arithmetic.
- `ec`: elliptic curves, isogenies and pairings. Everything that is purely
   finite-fieldy.
- `precomp`: constants and precomputed values.
- `quaternion`: quaternion orders and ideals.
- `hd`: code to compute (2,2)-isogenies in the theta model.
- `id2iso`: code for Ideal <-> Iso.
- `verification`: code for the verification protocol.
- `signature`: code for the key generation and signature protocols.

The dependencies are depicted below.
```
 ┌─┬──────────┬─┐        ┌─┬──────────┬─┐      ┌─┬──────────┬─┐
 │ ├──────────┤ │        │ ├──────────┤ │      │ ├──────────┤ │
 │ │  Keygen  │ │        │ │   Sign   │ │      │ │  Verify  │ │
 │ ├──────────┤ │        │ ├──────────┤ │      │ ├──────────┤ │
 └─┴────┬─────┴─┘        └─┴────┬─────┴─┘      └─┴────┬─────┴─┘
        │                       │                     │
        └──────────────────┐    │                     │
                           │    │                     │
┌─────────────────┐    ┌───▼────▼────────┐            │
│                 │    │                 │            │
│   Quaternions   ◄────┤  Ideal <-> Iso  ├────────┐   │
│                 │    │                 │        │   │
└────────┬────────┘    └────────┬────────┘        │   │
         │                      │                 │   │
         │                      │     ┌───────────────┘
         │                      │     │           │
┌────────▼────────┐    ┌────────▼─────▼──┐    ┌───▼────────────┐
│                 │    │                 │    │                │
│ Multiprecision  │    │       2D        ├────► Precomputation │
│ integers (GMP)  │    │    Isogenies    │    │                │
│                 │    │                 │    │                │
└─────────────────┘    └────────┬────────┘    └───▲────────────┘
                                │                 │
                                │                 │
                                │                 │
                       ┌────────▼────────┐        │
                       │                 │        │
                       │ Elliptic curves ├────────┘
                       │   & isogenies   │
                       │                 │
                       └──┬───────────┬──┘
                          │           │
                          │           │
                          │           │
              ┌───────────▼───┐   ┌───▼───────────┐
              │     GF(p)     │   │     Fixed     │
              │       &       │   │   precision   │
              │    GF(p^2)    │   │   integers    │
              └───────────────┘   └───────────────┘
```

## Cortex-M4 implementation

Verification routines are supported in 32-bit embedded architectures running on bare metal environments such as the ARM Cortex-M4, but they are not directly supported by the build system of the present repository. The [pqm4 project](https://github.com/mupq/pqm4) is supported for evaluating SQIsign verification in the ARM Cortex-M4.

pqm4 assumes that the full NIST API (keypair generation, signing and verification) is available. Since only verification is supported, the remaining routines are mocked and must meet certain constraints, such as sharing the public key, signing a prespecified message and being of a specific size used by the testing and benchmarking binaries of pqm4. Therefore, additional KATs must be generated specifically for pqm4, which is done by a dedicated KAT generator for pqm4, found in `apps/PQCgenKAT_sign_pqm4.c`.

A copy of the most recent version of pqm4 as of the round 2 submission deadline, including an implementation of SQIsign verification generated directly from this repository using the procedure explained next, is made available [here](https://github.com/SQISign/the-sqisign-pqm4), in the `sqisign` branch.

If changes are made to the library, the `scripts/gen_pqm4_sources.sh` shell script can be run, from the root folder of the repository, to generate a pqm4-compatible folder structure in `src/pqm4/sqisign_lvl{1,3,5}`, which can then be copied to the `crypto_sign` folder of pqm4. Note that the pqm4 KAT generator is automatically run by this script.

## Acknowledgements

The reference implementation for finite field arithemtic (i.e., `src/gf/ref`)
was generated using [modarith](https://github.com/mcarrickscott/modarith) by
Michael Scott.

## License

SQIsign is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

Third party code is used in some files:

- `src/common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `src/common/fips202.c`: CC0: Copyright (c) 2023, the PQClean team
- `src/common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>
- `src/common/broadwell/{aes_ni.c, vaes256_key_expansion.S}`: Apache-2.0: Copyright 2019 Amazon.com, Inc.
- `src/common/broadwell/ctr_drbg.c`: ISC: Copyright (c) 2017, Google Inc.
- `src/mini-gmp/mini-gmp.c` and `src/mini-gmp/mini-gmp.h`: LGPLv3: Copyright 1991-1997, 1999-2022 Free Software Foundation, Inc.
- `src/quaternion/ref/generic/dpe.h`: LGPLv3: Copyright (C) 2004-2024 Patrick Pelissier, Paul Zimmermann, LORIA/INRIA
- `apps/PQCgenKAT_sign.c`, `apps/PQCgenKAT_sign_pqm4.c`, `src/common/ref/randombytes_ctrdrbg.c`, `test/test_kat.c`: by NIST (Public Domain)
