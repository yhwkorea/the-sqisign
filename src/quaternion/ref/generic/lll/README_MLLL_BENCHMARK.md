# MLLL vs HNF Benchmark Results

Implementation of "Compact Quaternion Algorithms for SQIsign" (Kim, Lee, Yoo - Korea University).

Compares HNF-based vs MLLL-based ideal multiplication (`quat_lattice_mul` vs `quat_lattice_mul_mlll`) on random O0-ideal lattice pairs.

## Key Metric: Maximum Intermediate Bit Size

The paper's core claim is that MLLL bounds intermediate integer sizes to `max ||a_i||^2` (input norm squared), while HNF can blow up far beyond that. This is what determines whether **fixed-precision** arithmetic is feasible.

### NIST Level 1 (p ~ 2^253, norm bitsize = 127, 10 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **2163** | **255** | **0.118** |
| Avg intermediate bits | 1759 | 253 | 0.144 |
| Max output bits | 128 | 127 | ~1.0 |

### NIST Level 3 (p ~ 2^381, norm bitsize = 193, 5 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **3279** | **387** | **0.118** |
| Avg intermediate bits | 3267 | 385 | 0.118 |
| Max output bits | 194 | 193 | ~1.0 |

### NIST Level 5 (p ~ 2^509, norm bitsize = 254, 3 trials)

| Metric | HNF | MLLL | Ratio (MLLL/HNF) |
|--------|-----|------|-------------------|
| **Max intermediate bits** | **4312** | **510** | **0.118** |
| Avg intermediate bits | 3971 | 508 | 0.128 |
| Max output bits | 255 | 254 | ~1.0 |

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

2026-03-30
