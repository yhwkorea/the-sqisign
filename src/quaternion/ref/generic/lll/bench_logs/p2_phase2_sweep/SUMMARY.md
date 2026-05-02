# Phase 2 Primary-Decision Sweep (2026-04-21)

Driver: `run_sweep.sh`
Iterations per combo: 30
Build: CMake Release (`-O3`), WSL Ubuntu on Windows host
Host: user's workstation (non-isolated; single-run, wall-clock)

Columns are GRAM total wall time (ms) across 30 iterations per combo.
MLLL column is the baseline ibz-MLLL reference (non-Gram path), included
for calibration. Alg3 MLLL field shows `?` because `bench_one_alg3` does
not invoke the non-Gram MLLL path — Alg3 is GRAM-only.

## Table 1 — trap-on (default, production safe)

Built with `FP_CHECK_VEC` / `FP_CHECK_GRAM` enabled at every mutation site.

| Level | Mode | Variant  | GRAM (ms) | MLLL (ms) |
| :---: | :--: | :------- | --------: | --------: |
| L1    | alg2 | baseline |     12.78 |    153.58 |
| L1    | alg2 | prealloc |     12.81 |    143.05 |
| L1    | alg2 | fp       |     67.10 |    152.34 |
| L1    | alg3 | baseline |      0.33 |         — |
| L1    | alg3 | prealloc |      0.46 |         — |
| L1    | alg3 | fp       |      0.60 |         — |
| L3    | alg2 | baseline |     20.09 |    362.43 |
| L3    | alg2 | prealloc |     19.45 |    339.87 |
| L3    | alg2 | fp       |     19.79 |    352.24 |
| L3    | alg3 | baseline |      0.34 |         — |
| L3    | alg3 | prealloc |      0.34 |         — |
| L3    | alg3 | fp       |      0.38 |         — |
| L5    | alg2 | baseline |     26.44 |    656.98 |
| L5    | alg2 | prealloc |     26.11 |    645.28 |
| L5    | alg2 | fp       |     26.25 |    640.00 |
| L5    | alg3 | baseline |      0.33 |         — |
| L5    | alg3 | prealloc |      0.34 |         — |
| L5    | alg3 | fp       |      0.33 |         — |

## Table 2 — trap-off (`MLLL_FP_NO_OVERFLOW_CHECK=1`)

Measurement-only rebuild to isolate overflow trap overhead. Trap on
production code should remain enabled — this table only answers the
question: *if trap were zero-cost, would fp be competitive?*

| Level | Mode | Variant  | GRAM (ms) | MLLL (ms) |
| :---: | :--: | :------- | --------: | --------: |
| L1    | alg2 | baseline |     12.97 |    149.75 |
| L1    | alg2 | prealloc |     13.19 |    151.14 |
| L1    | alg2 | fp       |     60.28 |    142.63 |
| L1    | alg3 | baseline |      0.35 |         — |
| L1    | alg3 | prealloc |      0.47 |         — |
| L1    | alg3 | fp       |      0.58 |         — |
| L3    | alg2 | baseline |     19.53 |    339.04 |
| L3    | alg2 | prealloc |     19.79 |    361.41 |
| L3    | alg2 | fp       |     19.78 |    351.20 |
| L3    | alg3 | baseline |      0.36 |         — |
| L3    | alg3 | prealloc |      0.34 |         — |
| L3    | alg3 | fp       |      0.33 |         — |
| L5    | alg2 | baseline |     25.90 |    635.08 |
| L5    | alg2 | prealloc |     26.45 |    649.65 |
| L5    | alg2 | fp       |     26.70 |    672.66 |
| L5    | alg3 | baseline |      0.34 |         — |
| L5    | alg3 | prealloc |      0.35 |         — |
| L5    | alg3 | fp       |      0.36 |         — |

## Ratios (fp / baseline, GRAM column)

Trap-on table:

| Level | Alg2 ratio | Alg3 ratio |
| :---: | ---------: | ---------: |
| L1    |   **5.25** |   **1.82** |
| L3    |       0.98 |       1.12 |
| L5    |       0.99 |       1.00 |

Trap-off table:

| Level | Alg2 ratio | Alg3 ratio |
| :---: | ---------: | ---------: |
| L1    |   **4.65** |   **1.66** |
| L3    |       1.01 |       0.92 |
| L5    |       1.03 |       1.06 |

## Findings

1. **Candidate C (fp path) is not competitive at L1** — the production-
   hottest level. 4.65× slower even with the overflow trap disabled.
   The trap itself accounts for only ~10% of the fp slowdown
   (67.10 → 60.28 = 6.82 ms of 54+ ms regression). The remaining 90%
   is the fixed-width schoolbook (5 × 5 / 5 × 9 limbs) losing to GMP's
   single-limb fast path that dominates L1 values. GMP internally
   dispatches to `mpn_mul_1` or a 1-limb special case when operands fit
   in one limb, which is the common case at L1 (vec values ≤ 259 bits ≈
   4-5 non-zero limbs but most products are much smaller after size-
   reduce tightens). Our schoolbook loop always does the full
   `nwords_vec × nwords_vec` limb work.

2. **Candidate B (prealloc) ≈ baseline across the board** — confirms
   what the Phase-2-B commit already concluded: `mpz_realloc2` is a
   seed hint, not a ceiling, so `ibz_mul` / `ibz_add` still realloc on
   demand and the hint buys nothing. Expected from the `mpz_realloc2`
   semantics in GMP docs.

3. **At L3/L5 fp roughly matches baseline** — recorded as a time
   trade-off, not a decision input.

## Decision

**Primary backend is the fp path** (`quat_mlll_gram_fp`).
`g_fp_mode` defaults to 1. The ibz body is kept as a regression
oracle, reachable via `--fp 0` / `set_fp_mode(0)`.

**Why not "baseline wins because fp is 5.25× slower at L1"?** Because
time was never a Phase 2 criterion. The four criteria that *were*
stated up front (see FIXED_PRECISION_DECISION_KO.md §5.1) are:

1. Heap-free / stack layout — fp ✓, ibz ✗.
2. Runtime Lemma-3 audit (overflow trap) — fp ✓, ibz ✗.
3. Paper contribution in code (KLKL25: "compact arithmetic fits in
   fixed widths") — fp ✓, ibz ✗.
4. CT-compatible layout (for eprint 2025/2192 future direction) —
   fp ✓, ibz ✗.

fp is the only candidate that satisfies all four. The 5.25× L1
regression is accepted; addressing it is a Phase-2-post redesign
task (inline `mpn_mul_1` for the 1-limb L1 case, or
`mul + narrow-sub` fusion) that reuses the fp types in place.

**Note on first-take reversal.** An earlier draft of this SUMMARY
(same day, AM) concluded the opposite — "baseline remains primary" —
on the basis of time ratios alone. That draft was overruled once we
re-read the decision criteria this document was supposed to
implement. The sweep numbers above are unchanged; only the
interpretation was corrected.

### What the ibz body is kept for

- **fp regression oracle.** The `gram_fp_equivalence` test compares
  fp output against ibz output bit-for-bit.
- **Width-unresolved fallback.** When `alg->p` doesn't map to
  L1/L3/L5, the dispatcher silently routes to ibz.

### What to remove/keep

- **Keep (production path)**: `quat_fixed_precision.[hc]`,
  `quat_mlll_gram_fp`, dispatcher, `set_fp_mode`,
  `gram_fp_equivalence` test, `--fp` bench flag, overflow trap.
- **Keep (oracle)**: `quat_mlll_gram_ibz`, prealloc scaffold,
  `quat_test_mlll_gram_prealloc_equivalence`.
- **Default**: `g_fp_mode = 1`. Production always goes through fp.
- **Remove** (future): nothing for now — ibz oracle cost is
  negligible. Re-evaluate if maintenance burden grows.
