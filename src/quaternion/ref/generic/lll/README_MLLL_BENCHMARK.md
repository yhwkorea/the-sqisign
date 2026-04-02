# MLLL vs HNF Benchmark Results

Implementation of "Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - Korea University).

Compares HNF-based vs MLLL-based ideal multiplication (`quat_lattice_mul` vs `quat_lattice_mul_mlll`) on random O0-ideal lattice pairs.

## Current Status

The existing HNF-based ideal operations (`quat_lattice_mul`, `quat_lattice_add`) are **not replaced**. MLLL-based functions are implemented separately and compared against HNF on identical inputs to validate correctness and measure intermediate bit sizes.

## File Structure

| File | Description |
|------|-------------|
| `mlll.c` | **Core implementation**. Paper Algorithm 1 (MLLL) + Algorithm 2 (CompactIdealMultiplication). Provides `quat_mlll()`, `quat_lattice_mul_mlll()`, `quat_lattice_add_mlll()`. |
| `mlll_internals.h` | Function declarations and constants (`MLLL_MAX_GENERATORS=16`, etc.) |
| `mlll_tests.c` | **Correctness verification**. Runs both HNF and MLLL on same inputs, checks lattice equality via `quat_lattice_equal()`. 4 tests: (1) HNF vs MLLL result comparison, (2) LLL-reducedness of MLLL output, (3) linearly dependent generator handling, (4) full CompactIdealMultiplication flow. |
| `mlll_benchmark.c` | **Intermediate bit size comparison**. Uses `bitsize_tracker` to record max bit size of all integers during HNF/MLLL execution. Flow: `tracker_reset()` → run operation → `tracker_get_max()`. |
| `bitsize_tracker.h` | Benchmark-only instrumentation. `tracker_update_vec4()` calls inserted in `mlll.c` and `hnf.c` track intermediate integer sizes. No-op in normal builds (requires `BITSIZE_TRACKER_ENABLE`). |
| `test_mlll_only.c` | Standalone test entry point for MLLL tests only. |

## Paper Algorithm Mapping

| Paper | Code |
|-------|------|
| Algorithm 1 (MLLL) | `quat_mlll()` in `mlll.c:177-571` |
| Algorithm 2 (CompactIdealMultiplication) | `quat_lattice_mul_mlll()` in `mlll.c:576-617` |
| Lemma 1 (intermediate bound ≤ max\|\|a_i\|\|²) | Experimentally verified by `mlll_benchmark.c` |

## Not Yet Done

- Actual HNF → MLLL replacement in SQIsign pipeline (`lattice.c:177` still uses HNF)
- Algorithm 3 (RandomIdealGivenPrimeNorm)
- Algorithm 4 (RandomEquivalentPrimeIdeal) modification

## Key Metric: Maximum Intermediate Bit Size

The paper's core claim is that MLLL bounds intermediate integer sizes to `max ||a_i||^2` (input norm squared), while HNF can blow up far beyond that. This is what determines whether **fixed-precision** arithmetic is feasible.

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 10 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **2162** | **255** | **0.118** |
| Avg intermediate bits | 2073 | 253 | 0.122 |
| Max output bits | 128 | 127 | ~1.0 |
| Total time | 1.16 ms | 4102 ms | x3527 |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 5 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **3290** | **388** | **0.118** |
| Avg intermediate bits | 3265 | 385 | 0.118 |
| Max output bits | 194 | 194 | ~1.0 |
| Total time | 0.74 ms | 10269 ms | x13914 |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 3 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **4325** | **510** | **0.118** |
| Avg intermediate bits | 3975 | 508 | 0.128 |
| Max output bits | 255 | 254 | ~1.0 |
| Total time | 0.51 ms | 6511 ms | x12743 |

## Analysis

### Intermediate Bit Size

- **HNF intermediate values grow to ~17x the input bit size** (e.g., 128-bit input -> 2163-bit intermediate). This is driven by the 4x4 determinant used as the modulus in HNF computation.
- **MLLL intermediate values stay at ~2x the input bit size** (e.g., 128-bit input -> 255-bit intermediate). This matches the paper's Lemma 1 bound of `max ||a_i||^2`.
- The ratio MLLL/HNF is consistently **~0.118** across all security levels.

### Output Bit Size

Output basis element sizes are nearly identical between HNF and MLLL — both represent the same lattice, just in different forms (upper triangular vs LLL-reduced).

### Execution Time

MLLL is significantly slower than HNF in this GMP-based (arbitrary precision) implementation because:
1. MLLL uses exact rational Gram-Schmidt (`ibq_t` fractions) which is expensive.
2. HNF benefits from modular reduction keeping intermediate values bounded at each step.

However, the speed comparison is not the point. The paper targets **fixed-precision** SQIsign implementations where:
- HNF's 2000-4000 bit intermediates **cannot fit** in the fixed-precision budget.
- MLLL's ~500 bit intermediates **can fit**, enabling fixed-precision computation.

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

2026-03-30 (initial benchmark), 2026-04-02 (full test pass, re-run with timing data)
