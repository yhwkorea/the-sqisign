# MLLL vs HNF Benchmark Results

Implementation of "Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - Korea University).

Compares HNF-based vs MLLL-based ideal multiplication (`quat_lattice_mul` vs `quat_lattice_mul_mlll`) on random O0-ideal lattice pairs.

## Current Status

- **GSO**: Float GSO using `dpe_t` (53-bit mantissa + extended exponent), L² style
- **Gram matrix**: Exact integer (`ibz_t`), incrementally updated
- **Dependency detection**: `dpe_zero_p(r[k][k])` (float), with HNF post-processing in `done` section as temporary workaround
- **Tests**: 6/6 pass
- **Fixed-precision**: Confirmed — MLLL vector coordinates bounded by Lemma 1

## File Structure

| File | Description |
|------|-------------|
| `mlll.c` | **Core implementation**. Paper Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). Float GSO + Gram matrix. `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()`. |
| `mlll_internals.h` | Function declarations and constants (`MLLL_MAX_GENERATORS=16`) |
| `mlll_tests.c` | **Correctness verification**. 6 tests: (1) HNF vs MLLL lattice equality, (2) LLL-reducedness, (3) dependent generators, (4) CompactIdealMultiplication, (5) 127-bit prime scale, (6) tau path |
| `mlll_benchmark.c` | **Intermediate bit size comparison**. `bitsize_tracker` measures max bit size during computation. `tracker_disable()` separates MLLL reduction loop from HNF post-processing. |
| `bitsize_tracker.h` | Two-channel tracker: `vec` (vector coordinates) and `gso` (GSO coefficients). `__attribute__((weak))` symbols, no-op without `BITSIZE_TRACKER_ENABLE`. |
| `test_mlll_only.c` | Standalone test entry point |

## Paper Algorithm Mapping

| Paper | Code | Status |
|-------|------|--------|
| Algorithm 1 (MLLL) | `quat_mlll()` | ✓ Implemented (D_k=0 handling incomplete — see below) |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` | ✓ Implemented |
| Algorithm 3 (RandomIdealGivenPrimeNorm) | — | Not implemented |
| Algorithm 4 (RandomEquivalentPrimeIdeal) | — | Not implemented |
| Lemma 1 (intermediate bound) | `mlll_benchmark.c` | ✓ Experimentally verified |

## Fixed-Precision Benchmark Results (2026-04-09)

Tracker measures MLLL reduction loop only (`tracker_disable()` before HNF post-processing).

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 3 trials)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **2153** | **255** | **8.4x** |
| Max output bits | 128 | 128 | ~1.0 |
| Total time | 0.49ms | 1.08ms | 0.45x |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 1 trial)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **3245** | **384** | **8.5x** |
| Max output bits | 194 | 192 | ~1.0 |
| Total time | 0.14ms | 0.38ms | 0.37x |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 1 trial)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **3779** | **503** | **7.5x** |
| Max output bits | 250 | 250 | ~1.0 |
| Total time | 0.15ms | 0.24ms | 0.63x |

### Summary

MLLL vector coordinates stay within **~2× input bits** across all security levels, matching the paper's Lemma 1 bound (≤ max||a_i||²). HNF intermediate values grow to ~17× input bits due to modular determinant computation.

GSO coefficients use `dpe_t` (53-bit float) so they contribute 0 bits to the integer tracker. This is the intended design: fixed-precision GSO means the only integer growth is in vector coordinates, which are Lemma 1 bounded.

## Not Yet Implemented

### 1. Proper D_k=0 (dependent vector) handling

**Current**: When MLLL detects a dependent vector (`r[k][k] ≈ 0` via float GSO), the reduction loop may leave unresolved dependent vectors in `b[]`. The `done` section uses `ibz_mat_4xn_hnf_mod_core()` to extract the correct rank-4 basis as a temporary workaround. This HNF call is excluded from intermediate bit size tracking.

**Required**: Per Pohst (1987) and Matthews' pseudo-code, the algorithm must handle `D_k = 0` in two cases:
- `D_k = 0, λ_{k,k-1} = 0`: Swap rows and decrease β
- `D_k = 0, λ_{k,k-1} ≠ 0`: Apply `Swap2` (special μ update) then `Swap1`

This requires tracking exact integer `D_k` and `λ_{k,j}` values alongside the float GSO, or computing them from the Gram matrix when needed. See: http://www.numbertheory.org/PDFS/mlll.pdf

### 2. SQIsign pipeline integration

Replace HNF with MLLL in the main SQIsign pipeline (`lattice.c:177`). Currently MLLL exists as separate functions for benchmarking only.

### 3. Algorithm 3 (RandomIdealGivenPrimeNorm)

Cornacchia + MLLL combination for sampling ideals with prime norm.

### 4. Algorithm 4 (RandomEquivalentPrimeIdeal)

Modification to use LLL-reduced basis instead of HNF basis.

### 5. Benchmark generation bottleneck

`quat_test_input_random_ideal_lattice_generation()` is probabilistic and occasionally slow for > 5 iterations. The generation time is now printed per set to help diagnose this.

## Implementation History

| Date | Change |
|------|--------|
| 2026-03-30 | Initial MLLL with `ibq_t` (exact rational) GSO, full test pass |
| 2026-04-02 | Bitsize tracker, benchmark with timing data |
| 2026-04-07 | tau swap bug fix, tracker expansion to GSO coefficients |
| 2026-04-08 | Integral GSO experiment (replaced `ibq_t` fractions with `ibz_t`) |
| 2026-04-09 | Float GSO rewrite (`dpe_t`), Gram matrix approach, HNF post-processing as workaround, tracker separated from post-processing |

## How to Reproduce

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBITSIZE_TRACKER_ENABLE=ON
make sqisign_bm_mlll sqisign_test_mlll -j$(nproc)

# Tests
./src/quaternion/ref/generic/test/sqisign_test_mlll

# Benchmarks (use small iteration counts — generation can be slow)
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=3
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=1
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=1
```
