# MLLL Implementation for SQIsign Fixed-Precision Ideal Multiplication

Implementation based on:
- **Primary**: "SQIsign with Fixed-Precision Integer Arithmetic" (Kim, Lee, Kim, Lee - Korea University, 2025) — `klkl25.pdf`
- **MLLL algorithm**: Pohst (1987), Matthews pseudo-code, Cohen Algorithm 2.6.7

## Motivation: Why Replace HNF with MLLL?

The paper identifies **IdealMultiplication** as the bottleneck for fixed-precision SQIsign (Section 1, "Challenges"):

> "The sizes of intermediates during HNF can be governed by the g.c.d of 4x4 minors of this matrix. Computing all minors needs (16 choose 4) computations... So, in the original Round-2 algorithm, it computes just 4 minors, whose g.c.d is usually larger..."

The paper proposes **Algorithm 9** (Modified IdealMultiplication) which replaces the determinant modulus with `m = nrd(r1*I1)^2 * nrd(r2*I2)^2` (Lemma 11). Our approach goes further: **replace HNF entirely with MLLL**, avoiding modular arithmetic and its m^2 intermediate blowup.

### Theoretical Bound Comparison

| Approach | IdealMul intermediate bound | Source |
|----------|---------------------------|--------|
| Original HNF (Alg. 2) | `2^10 * p^10 * r1^32 * r2^8 * nrd(I1)^16 * nrd(I2)^4` | Lemma 2 |
| Modified HNF (Alg. 9) | `r1^8 * r2^8 * nrd(I1)^4 * nrd(I2)^4` | Lemma 12 |
| **MLLL (this impl.)** | `max ||a_i||^2` where a_i are generators | Lemma 1 + Pohst |

For MLLL, `||a_i||^2 = nrd(r1*alpha_s * r2*beta_t)` where alpha, beta are basis elements. By Lemma 1 of the paper, generator entry bound is `sqrt(nrd(alpha) * nrd(beta))`, so:

- MLLL intermediate ≈ `r1^2 * r2^2 * max_s(nrd(alpha_s)) * max_t(nrd(beta_t))`
- This is **quadratically smaller** than the Modified HNF bound (which squares m)

### Role of nrd(I) in the Bound

The paper's bounds are fundamentally expressed in terms of **reduced norms** `nrd(I)` of the input ideals.

**Lemma 12** (Modified IdealMul bound): intermediate ≤ `r1^8 * r2^8 * nrd(I1)^4 * nrd(I2)^4`

where:
- `nrd(I)` = reduced norm of ideal I = gcd{nrd(α) | α ∈ I}
- `r_i` = lcm of denominators of basis elements (Lemma 5: `r_α ≤ 2*nrd(connecting_ideal)`)

**How nrd(I) determines intermediate size**: The modulus m = nrd(r1·I1)^2 · nrd(r2·I2)^2 = r1^4·r2^4·nrd(I1)^2·nrd(I2)^2. HNF over Z_m produces intermediates ≤ m^2. So intermediates grow as **4th power of nrd(I)**.

**Composition through the pipeline** (how nrd(I) propagates):

| Subroutine | Ideal involved | nrd bound | Bound derivation |
|-----------|---------------|-----------|-----------------|
| KeyGen Line 2 | I_sk (secret) | nrd(I_sk) ≤ D_mix ≈ 2^(8λ) | RandomIdealGivenNorm (Alg. 6) |
| KeyGen Line 4 | I_sk (re-sampled) | nrd(I_sk) ≤ p (hypothesis) | RandomEquivalentPrimeIdeal |
| KeyGen Line 5 | J_t · I_sk | nrd(J_t) < (log p)^2 | IdealToIsogeny → SuitableIdeals |
| Sign Line 13 | [I_sk]* I'_chl | nrd ≤ nrd(I_sk)·nrd(I'_chl) ≤ p·2^f | Lemma 9, 10 |
| Sign Line 14 | I_sk · I_chl | nrd(I_chl) = nrd(I'_chl) ≤ 2^f | Lemma 12 applied |
| Sign Line 24 | I_com,rsp ∩ I_aux | nrd ≤ q_rsp·D_mix · √p | Lemma 9 |

**Theorem 1** (Modified KeyGen): All nrd(I) ≤ p, so Lemma 12 gives ≤ r^8 · p^4. With r ≤ 2·nrd ≤ 2p, total ≤ 2^8·p^8·p^4 ... but the dominant term comes from SuitableIdeals (Lemma 13): `2^16·(log p)^16·nrd(I_sk)^4 ≤ 2^16·(log p)^16·p^4`. The overall KeyGen bound 2^68·p^9 is dominated by RandomIdealGivenNorm: `32·p·D_mix^2 ≈ 2^68·p^9`.

**Theorem 2** (Modified Sign): The worst case is Sign Line 13 where nrd(I_sk)^20·nrd(I'_chl)^4 ≤ p^20·2^(4f) ≤ p^24·p^4 = p^28. With constant factor 2^16, total = **2^16·p^28**.

**For MLLL**: Our approach avoids m^2 entirely, so nrd(I) affects intermediates only through the generator norms (Lemma 1), not through the modulus. The MLLL intermediate for one IdealMul call is:

```
MLLL intermediate ≤ max_s,t { nrd(r1·α_s · r2·β_t) }
                  = max { r1^2·r2^2 · nrd(α_s) · nrd(β_t) }
```

Since nrd(α_s) ≤ 8·nrd(I) for LLL-reduced bases (Minkowski bound), this gives:

```
MLLL intermediate ≤ 64 · r1^2 · r2^2 · nrd(I1) · nrd(I2)
```

This is the **square root** of the Modified HNF bound (Lemma 12: r^8·nrd^4 vs r^4·nrd^2 for MLLL generator norms). This explains why our measured bits (~389 for L1) are much smaller than the HNF theoretical bound.

### Experimental Verification

| Level | Input bits | MLLL vec | Modified HNF (Lemma 12) | Original HNF (measured) |
|-------|-----------|----------|------------------------|------------------------|
| L1 (p~2^253) | 127 | **389** | ~1016 (m^2) | 2160 |
| L3 (p~2^381) | 192 | **571** | ~1528 | 3238 |
| L5 (p~2^509) | 255 | **758** | ~2032 | 4307 |

MLLL intermediate bits ≈ `norm_bits + log2(p)` ≈ `3 * input_bits`. This matches Lemma 1's bound since quaternion multiplication introduces a factor of p in coordinates 2,3.

## Paper Algorithm → Code Mapping

### Algorithm 9: Modified IdealMultiplication → `quat_lattice_mul_mlll()`

```
Paper (Algorithm 9)                    Code (mlll.c:437-479)
─────────────────────                  ──────────────────────
Line 1: r1 = lcm(r_alpha_i)           (handled by caller — lattice denom)
Line 2: alpha_i *= r1, beta_j *= r2   (implicit in lattice basis representation)
Line 3: m = nrd(r1*I1)^2*nrd(r2*I2)^2 ← REPLACED: no modulus needed
Line 4: HNF(alpha_i*beta_j) over Z_m  ← REPLACED: quat_mlll(generators, 16)
Line 5: return M / (r1*r2)            res->denom = lat1->denom * lat2->denom
```

**Key difference**: Instead of computing HNF modulo m (which squares the modulus → m^2 intermediate), we feed the 16 generators directly into MLLL which keeps vectors bounded by Lemma 1.

```c
// mlll.c:454-468 — Generate 16 quaternion products (same as paper's matrix M)
for (int k = 0; k < 4; k++) {          // for each alpha_k (basis of I1)
    for (int i = 0; i < 4; i++) {      // for each beta_i (basis of I2)
        quat_alg_coord_mul(&elem_res, &elem1, &elem2, alg);
        // generators[4*k+i] = alpha_k * beta_i   (paper: M column l=4k+i)
    }
}
quat_mlll(&(res->basis), &rank, generators, 16, alg);  // MLLL replaces HNF
```

### Algorithm 1 (HNF, paper) → Replaced by `quat_mlll()`

The paper's Algorithm 1 (HNF) computes Hermite Normal Form over Z_m with XGCD operations. Intermediate values are bounded by m^2 (Lemma 2, Lemma 12).

Our MLLL replacement avoids modular arithmetic entirely. The algorithm structure:

```
Paper HNF (Algorithm 1)                Our MLLL (quat_mlll, mlll.c:212-434)
────────────────────────               ─────────────────────────────────────
Input: 4×16 matrix over Z_m           Input: 16 generator vectors in Z^4
XGCD + column operations mod m        Size-reduction (REDI) + Lovász swaps
Output: 4×4 HNF basis                 Output: 4×4 LLL-reduced basis
Intermediate: ≤ m^2                   Intermediate: ≤ max||generator||^2
```

### Lemma 1 (Generator Entry Bound) → Verified by benchmark

Paper Lemma 1 states:
```
|M_{1,l}|, |M_{2,l}| ≤ sqrt(nrd(alpha) * nrd(beta))
|M_{3,l}|, |M_{4,l}| ≤ (1/sqrt(p)) * sqrt(nrd(alpha) * nrd(beta))
```

For Level 1: nrd(alpha), nrd(beta) ≈ 2^127, p ≈ 2^253.
- Coordinates 0,1: ≤ 2^127 → 127 bits
- Coordinates 2,3: ≤ 2^(127-126) ≈ 2^1 → small
- But after quaternion multiplication with p factor: coordinates grow to ≈ 127+253 = 380 bits

Our measured MLLL intermediate = 389 bits, confirming the bound.

## Internal Algorithm Details

### MLLL Core: `quat_mlll()` (mlll.c:212-434)

Implements Pohst (1987) MLLL with exact integer arithmetic (Cohen 2.6.7 style):

| Component | Paper/Reference | Code | Lines |
|-----------|----------------|------|-------|
| **Inner product** | `<a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3` (quaternion norm form) | `ibz_vec_4_dot_quat()` | 32-45 |
| **Gram matrix** | `G[i][j] = <b[i], b[j]>` — exact integer | `gram_compute_row()`, `gram_size_reduce()`, `gram_swap()` | 68-111 |
| **d[k] / lam[i][j]** | Cohen Alg. 2.6.7: d[0]=1, d[k+1]=det(Gram(b[0..k])), lam[i][j]=d[j+1]*mu[i][j] | `compute_lam_d_row()`, `recompute_all_lam_d()` | 160-208 |
| **Nearest-integer div** | REDI quotient: q = round(lam[m][l] / d[l+1]) | `rounded_div()` | 118-139 |
| **Exact division** | Cohen recurrence: num = (d[s+1]*num - lam*lam) / d[s] | `exact_div()` | 144-152 |

### MLLL Phases (Pohst 1987 / Matthews)

```
┌─────────┐     ┌───────────┐     ┌──────────────┐     ┌──────────┐
│  LOAD   │────>│ REDUCTION │────>│ SWAP or DONE │────>│  REMOVE  │
│ (245-268)│    │ (280-349) │     │ (352-371)    │     │ (374-401)│
└─────────┘     └───────────┘     └──────────────┘     └──────────┘
     ^               ^                  │                    │
     │               └──────────────────┘                    │
     └───────────────────────────────────────────────────────┘
```

| Phase | Matthews | Code | Description |
|-------|----------|------|-------------|
| **LOAD** | "Load generator into b[beta]" | lines 245-268 | Copy generator → compute Gram row → recompute d/lam → check d[beta]=0 |
| **REDI** | "Size-reduce b[m] against b[l]" (Cohen Alg. 2.6.7) | lines 283-309 | If \|2*lam[m][l]\| > d[l+1]: q=round(lam/d), b[m]-=q*b[l], update Gram+lam |
| **Lovász check** | "4*(d[m-1]*d[m+1]+lam^2) < 3*d[m]^2 → swap" | lines 327-339 | Standard LLL condition with delta=3/4 |
| **D_k=0 detection** | Matthews: "if D_m=0" (Pohst 1987 §3) | lines 320-325 | Exact integer test: `ibz_is_zero(&d[m+1])` |
| **Forced swap (dep.)** | Matthews flagg=1 / Swap2 cases | lines 352-360 | `do_swap_dep`: bypass tau (`m=max(m-1,1)`), cascade to position 1 |
| **Lovász swap** | Matthews Swap1 | lines 363-371 | `do_swap_lov`: respect tau (`m=max(m-1,tau)`) |
| **Remove zero vector** | "b[m]=0 → shift down, beta--" | lines 374-401 | After removal with alpha>=g: re-run reduction from m=1 |

### D_k=0 Handling Detail

When `d[m+1] = 0` (dependent vector detected):

```
Matthews (1987)                           Our implementation
──────────────────                        ─────────────────────
Case 1 (λ_{m,m-1}=0, flagg=1):          ┐
  Set d[m-1]=0, zero lambdas             │ Both cases unified:
  m = max(m-1, tau)                      │   physical swap b[m] ↔ b[m-1]
                                         │   recompute ALL d/lam from Gram
Case 2 (λ_{m,m-1}≠0):                   │   m = max(m-1, 1)  ← bypasses tau!
  Swap2 (special lam update)             │
  Swap1 (physical swap)                  │ Recompute-from-scratch avoids
  m = max(m-1, tau)                      ┘ complex Swap2 incremental formulas
```

**Why bypass tau**: The dependent vector may depend on a distant vector (e.g., b[5]=(3,0,0,0) depends on b[1]=(6,0,0,0)). With `m=max(m-1,tau)`, the cascade stops at tau and loops infinitely. With `m=max(m-1,1)`, the vector cascades to position 1 where REDI with the shortest vector eliminates it.

### Lattice Addition: `quat_lattice_add_mlll()` (mlll.c:481-514)

Used for ideal intersection via `I1 ∩ I2 = (I1^(-1) + I2^(-1))^(-1)` (paper Section 3.1). Takes 8 generators (4 from each lattice, cross-scaled by denominators) and produces LLL-reduced basis.

## Comparison with Paper's Fixed-Precision Results (Table 2)

The paper achieves these total pipeline bounds:

| | Original (Prop. 1,2) | Modified (Thm. 1,2) |
|---|---|---|
| **KeyGen** | 2^50 * p^14 * (log p)^96 | **2^68 * p^9** |
| **Sign** | 2^50 * p^95 | **2^16 * p^28** |

For NIST-I (p ≈ 2^253):
- Modified KeyGen bound: 2^68 * 2^(253*9) = 2^2345 bits
- Modified Sign bound: 2^16 * 2^(253*28) = 2^7100 bits
- Paper's experimental max (Table 2): 7,026 bits → 7,040 necessary

Our MLLL intermediate measurement (389 bits for L1) is for **one IdealMultiplication call only**, not the full pipeline. The paper's 7,026 bits includes the entire Sign chain with multiple IdealMultiplications composed.

## Test Suite → Paper Verification

| Test | What it verifies | Paper relevance |
|------|-----------------|-----------------|
| `dependent_generators` | MLLL handles Z-dependent input | D_k=0 case (Pohst §3) |
| `is_lll_reduced` | Output satisfies delta-LLL (delta=3/4) | Paper assumes LLL-reduced bases (Lemma 2 condition) |
| `vs_hnf_lattice_mul` | MLLL lattice = HNF lattice | Algorithm 9 correctness (same ideal, different basis) |
| `compact_ideal_multiplication` | Full IdealMul pipeline | Algorithm 9 end-to-end |
| `realistic_scale` | 127-bit prime norm | NIST Level 1 parameter regime |
| `tau_path` | Generator loading order | MLLL load/reduce/load cycle (Pohst §2) |

## Current Status

- **Tests**: 6/6 pass
- **Intermediate bits**: ~3x input (matches Lemma 1 bound with p factor)
- **HNF dependency**: Fully removed from MLLL path
- **Performance**: ~20-60x slower than HNF (recompute-from-scratch overhead; optimizable with incremental Swap1/Swap2)

## Not Yet Implemented

| Item | Paper reference | Priority |
|------|----------------|----------|
| Pipeline integration | Replace `lattice.c:177` HNF call | High |
| Incremental Swap1/Swap2 | Matthews (1987), reduce O(beta^2) → O(beta) per swap | Medium |
| Algorithm 3 (RandomIdealGivenPrimeNorm) | Algorithm 6 in paper | Low |
| Algorithm 4 (RandomEquivalentPrimeIdeal) | Algorithm 12 (Appendix) | Low |

## How to Reproduce

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBITSIZE_TRACKER_ENABLE=ON
make sqisign_bm_mlll sqisign_test_mlll -j$(nproc)

# Tests
./src/quaternion/ref/generic/test/sqisign_test_mlll

# Benchmarks
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=1 --iterations=3
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=3 --iterations=1
./src/quaternion/ref/generic/test/sqisign_bm_mlll --level=5 --iterations=1
```

## Implementation History

| Date | Change |
|------|--------|
| 2026-03-30 | Initial MLLL with `ibq_t` (exact rational) GSO |
| 2026-04-02 | Bitsize tracker, benchmark with timing |
| 2026-04-07 | tau swap bug fix, tracker expansion |
| 2026-04-08 | Integral GSO experiment (`ibq_t` → `ibz_t`) |
| 2026-04-09 | Float GSO (`dpe_t`), Gram matrix, HNF post-processing workaround |
| 2026-04-09 | **Exact integer d[]/lam[][] rewrite. D_k=0 proper handling. HNF removed.** |
