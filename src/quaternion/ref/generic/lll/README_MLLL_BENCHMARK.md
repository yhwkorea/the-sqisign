# MLLL vs HNF Benchmark Results

Implementation of "Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - Korea University).

Compares HNF-based vs MLLL-based ideal multiplication (`quat_lattice_mul` vs `quat_lattice_mul_mlll`) on random O0-ideal lattice pairs.

## Current Status

- **GSO**: Exact integer d[]/lam[][] (Cohen 2.6.7 style). d[k+1] = det(Gram(b[0..k])), lam[i][j] = d[j+1]*mu[i][j].
- **Gram matrix**: Exact integer (`ibz_t`), incrementally updated during size-reduce, recomputed after row shifts
- **Dependency detection**: `d[m+1] == 0` (exact integer test). Dependent vectors are eliminated via forced-swap cascade + REDI, no HNF post-processing needed.
- **Tests**: 6/6 pass
- **Fixed-precision**: Confirmed — MLLL vector coordinates bounded by ~3x input bits

## File Structure

| File | Description |
|------|-------------|
| `mlll.c` | **Core implementation**. Paper Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). Exact integer d[]/lam[][] with recompute-from-scratch after swaps. `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()`. |
| `mlll_internals.h` | Function declarations and constants (`MLLL_MAX_GENERATORS=16`) |
| `mlll_tests.c` | **Correctness verification**. 6 tests: (1) HNF vs MLLL lattice equality, (2) LLL-reducedness, (3) dependent generators, (4) CompactIdealMultiplication, (5) 127-bit prime scale, (6) tau path |
| `mlll_benchmark.c` | **Intermediate bit size comparison**. `bitsize_tracker` measures max bit size during computation. |
| `bitsize_tracker.h` | Two-channel tracker: `vec` (vector coordinates) and `gso` (GSO coefficients). `__attribute__((weak))` symbols, no-op without `BITSIZE_TRACKER_ENABLE`. |
| `test_mlll_only.c` | Standalone test entry point |

## Paper Algorithm Mapping

| Paper | Code | Status |
|-------|------|--------|
| Algorithm 1 (MLLL) | `quat_mlll()` | Done (exact integer d/lam, proper D_k=0 handling) |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` | Done |
| Algorithm 3 (RandomIdealGivenPrimeNorm) | — | Not implemented |
| Algorithm 4 (RandomEquivalentPrimeIdeal) | — | Not implemented |
| Lemma 1 (intermediate bound) | `mlll_benchmark.c` | Experimentally verified |

## D_k=0 Handling (Matthews/Pohst)

When `d[m+1] = 0` (dependent vector detected via exact integer test):

1. **Forced swap**: Exchange b[m] and b[m-1], bypass tau barrier (`m = max(m-1, 1)` instead of `m = max(m-1, tau)`)
2. **Cascade**: The dependent vector bubbles down toward position 1, where REDI against a shorter vector can eliminate it
3. **Recompute**: All d[]/lam[][] values are recomputed from the Gram matrix after each swap (O(beta^2), negligible for beta <= 16)
4. **Final reduction**: After dependent vector removal, full LLL reduction is re-run from m=1 to ensure the output is properly size-reduced

This replaces the previous HNF post-processing workaround with a mathematically correct approach.

## Exact Integer Benchmark Results (2026-04-09)

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 3 trials)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **2160** | **389** | **5.5x** |
| Max output bits | 128 | 127 | ~1.0 |
| Total time | 0.46ms | 8.17ms | 0.06x |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 1 trial)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **3238** | **571** | **5.7x** |
| Max output bits | 189 | 189 | ~1.0 |
| Total time | 0.14ms | 4.74ms | 0.03x |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 1 trial)

| Metric | HNF | MLLL vec | Improvement |
|--------|-----|----------|-------------|
| **Max intermediate bits** | **4307** | **758** | **5.7x** |
| Max output bits | 252 | 252 | ~1.0 |
| Total time | 0.15ms | 9.77ms | 0.02x |

### Summary

MLLL vector coordinates stay within **~3x input bits** across all security levels. The previous float GSO version achieved ~2x, but the exact integer version trades tighter bounds for correctness guarantees (proper D_k=0 handling without HNF fallback).

HNF intermediate values grow to ~17x input bits due to modular determinant computation, so the MLLL intermediate bit advantage is still **5.5-5.7x**.

GSO is now exact integer (d[]/lam[][]), so the "gso" tracker channel shows 0 (d/lam values are not separately tracked as they are internal to the algorithm).

### Performance Note

The recompute-from-scratch approach (recompute all d[]/lam[][] after every swap) is slower than the previous float GSO version. This is expected:
- Float GSO: O(1) per dpe_t operation, fewer swaps
- Exact integer: O(beta^2) per recompute, more forced swaps for dependent vector cascading

Future optimization: implement incremental Swap1/Swap2 formulas (Matthews) instead of full recomputation. This would reduce the per-swap cost from O(beta^2) to O(beta).

## Not Yet Implemented

### 1. SQIsign pipeline integration

Replace HNF with MLLL in the main SQIsign pipeline (`lattice.c:177`). Currently MLLL exists as separate functions for benchmarking only.

### 2. Incremental Swap1/Swap2 (performance optimization)

Replace recompute-from-scratch with Matthews' incremental update formulas for d[]/lam[][] during swaps. This would improve MLLL execution time from ~20-60x HNF to ~2-3x HNF (matching the previous float GSO version's speed while keeping exact integer correctness).

### 3. Algorithm 3 (RandomIdealGivenPrimeNorm)

Cornacchia + MLLL combination for sampling ideals with prime norm.

### 4. Algorithm 4 (RandomEquivalentPrimeIdeal)

Modification to use LLL-reduced basis instead of HNF basis.

## Implementation History

| Date | Change |
|------|--------|
| 2026-03-30 | Initial MLLL with `ibq_t` (exact rational) GSO, full test pass |
| 2026-04-02 | Bitsize tracker, benchmark with timing data |
| 2026-04-07 | tau swap bug fix, tracker expansion to GSO coefficients |
| 2026-04-08 | Integral GSO experiment (replaced `ibq_t` fractions with `ibz_t`) |
| 2026-04-09 | Float GSO rewrite (`dpe_t`), Gram matrix approach, HNF post-processing workaround |
| 2026-04-09 | **Exact integer d[]/lam[][] rewrite**. Proper D_k=0 handling (forced-swap cascade). HNF post-processing removed. Recompute-from-scratch for robustness. 6/6 tests pass. |

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
