# MLLL vs HNF Benchmark Results

Implementation of "Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - Korea University).

Compares HNF-based vs MLLL-based ideal multiplication (`quat_lattice_mul` vs `quat_lattice_mul_mlll`) on random O0-ideal lattice pairs.

## Current Status

The existing HNF-based ideal operations (`quat_lattice_mul`, `quat_lattice_add`) are **not replaced**. MLLL-based functions are implemented separately and compared against HNF on identical inputs to validate correctness and measure intermediate bit sizes.

## File Structure

| File | Description |
|------|-------------|
| `mlll.c` | **Core implementation**. Paper Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). Provides `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()`. Korean docstring describes Pohst 1987 algorithm and fixed-precision limitations. |
| `mlll_internals.h` | Function declarations and constants (`MLLL_MAX_GENERATORS=16`, etc.) |
| `mlll_tests.c` | **Correctness verification**. Runs both HNF and MLLL on same inputs, checks lattice equality via `quat_lattice_equal()`. 6 tests: (1) HNF vs MLLL result comparison, (2) LLL-reducedness of MLLL output, (3) linearly dependent generator handling, (4) full CompactIdealMultiplication flow, (5) realistic 127-bit prime scale, (6) tau > 1 swap path after remove_vector. |
| `mlll_benchmark.c` | **Intermediate bit size comparison**. Uses `bitsize_tracker` to record max bit size of all integers during HNF/MLLL execution. Flow: `tracker_reset()` → run operation → `tracker_get_max()`. |
| `bitsize_tracker.h` | Benchmark-only instrumentation. `tracker_update_vec4()` / `tracker_update_ibz()` calls inserted in `mlll.c` and `hnf.c` track intermediate integer sizes. Uses `__attribute__((weak))` symbols so any TU can link without a separate definition. No-op in normal builds (requires `BITSIZE_TRACKER_ENABLE`). |
| `test_mlll_only.c` | Standalone test entry point for MLLL tests only. |

## Paper Algorithm Mapping

| Paper | Code |
|-------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` in `mlll.c:192-586` |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` in `mlll.c:591-632` |
| Lemma 1 (intermediate bound ≤ max\|\|a_i\|\|²) | Experimentally verified by `mlll_benchmark.c` |

## Changes (2026-04-07)

- **Bug fix**: swap decrement condition changed from `m > 1` to `m > tau` — after `remove_vector` sets `tau = m+1`, the old condition could skip below `tau` and access invalid GS data.
- **Tracker coverage**: added `tracker_update_ibz()` for GS inner products, cross products, mu numerator/denominator, and B numerator/denominator during `compute_gs_single()`.
- **Weak symbols**: `bitsize_tracker.h` now uses `__attribute__((weak))` so linking works without `mlll_benchmark.c`.
- **Rank assertion**: `quat_lattice_add_mlll()` asserts `rank > 0 && rank <= 4`.
- **Korean docstring**: `mlll.c` header rewritten in Korean (Pohst 1987, fixed-precision limitations).
- **New tests**: `quat_test_mlll_realistic_scale` (127-bit prime), `quat_test_mlll_tau_path` (tau > 1 swap path).

## Not Yet Done

- Actual HNF → MLLL replacement in SQIsign pipeline (`lattice.c:177` still uses HNF)
- Algorithm 3 (RandomIdealGivenPrimeNorm)
- Algorithm 4 (RandomEquivalentPrimeIdeal) modification

## Key Metric: Maximum Intermediate Bit Size

The paper's core claim is that MLLL bounds intermediate integer sizes to `max ||a_i||^2` (input norm squared), while HNF can blow up far beyond that. This is what determines whether **fixed-precision** arithmetic is feasible.

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 10 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **2168** | **4513** | **2.082** |
| Avg intermediate bits | 2106 | 3596 | 1.707 |
| Max output bits | 128 | 127 | ~1.0 |
| Total time | 1.35 ms | 413.79 ms | x306 |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 5 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **3285** | **6881** | **2.095** |
| Avg intermediate bits | 2907 | 4788 | 1.647 |
| Max output bits | 194 | 193 | ~1.0 |
| Total time | 0.72 ms | 441.16 ms | x614 |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 3 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **4315** | **7699** | **1.784** |
| Avg intermediate bits | 3440 | 6931 | 2.015 |
| Max output bits | 254 | 253 | ~1.0 |
| Total time | 0.47 ms | 674.98 ms | x1449 |

## Analysis

### Intermediate Bit Size

**Important**: The tracker now covers both integer vector coordinates AND Gram-Schmidt rational coefficients (`ibq_t` numerators/denominators). Previous measurements (2026-04-02) only tracked vector coordinates, yielding misleadingly small MLLL intermediates (~255 bits). The updated tracker reveals the full picture:

- **HNF intermediate values grow to ~17x the input bit size** (e.g., 128-bit input -> 2168-bit intermediate). This is driven by the 4x4 determinant used as the modulus in HNF computation.
- **MLLL intermediate values grow to ~2x HNF** (e.g., 128-bit input -> 4513-bit intermediate). This is dominated by the exact rational Gram-Schmidt coefficients (`ibq_t` fractions), whose numerators/denominators grow without bound in this arbitrary-precision implementation.
- The paper's Lemma 1 bound (`max ||a_i||^2`) applies to **integer vector coordinates only**, not to GS coefficients. The vector coordinates do stay bounded as claimed (~255 bits for Level 1).

### Output Bit Size

Output basis element sizes are nearly identical between HNF and MLLL — both represent the same lattice, just in different forms (upper triangular vs LLL-reduced).

### Execution Time

MLLL is slower than HNF in this GMP-based (arbitrary precision) implementation because:
1. MLLL uses exact rational Gram-Schmidt (`ibq_t` fractions) which is expensive.
2. HNF benefits from modular reduction keeping intermediate values bounded at each step.

The speed comparison is not the point. The paper targets **fixed-precision** SQIsign implementations where:
- HNF's 2000-4000 bit intermediates (vector coordinates) **cannot fit** in the fixed-precision budget.
- MLLL's ~500 bit vector coordinates **can fit**, enabling fixed-precision computation.
- A fixed-precision MLLL would use integral GSO or L² representation instead of `ibq_t`, eliminating the GS coefficient blowup seen here.

## How to Reproduce

```bash
# Build
cd build_wsl
cmake .. -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release
make sqisign_bm_mlll -j4

# Run
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=10
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=5
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=3
```

## Date

2026-03-30 (initial benchmark), 2026-04-02 (full test pass, re-run with timing data), 2026-04-07 (tau bug fix, tracker expansion, new tests)
