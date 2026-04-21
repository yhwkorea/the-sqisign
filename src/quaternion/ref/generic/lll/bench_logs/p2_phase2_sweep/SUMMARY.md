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

3. **At L3/L5 fp roughly matches baseline** — but "matches" is not
   "beats". Neither candidate dominates, and we lose the L1
   production path. No level has a net win.

## Decision

**Primary backend remains the baseline ibz_t path** (the existing
`quat_mlll_gram_ibz` body). Neither B (prealloc) nor C (fp) yields a
competitive speedup at any security level; C actively regresses L1.

### What the fp scaffold is kept for

- **Upper-bound regression guard.** The `gram_fp_equivalence` test
  continues to prove the two paths produce identical lattices.
  Anything that breaks that test flags a correctness bug in either
  arithmetic path.
- **Runtime Lemma-3 audit.** The overflow trap runs whenever the fp
  path is opted into; a trip means our widths table or the Lemma 3
  argument is wrong.
- **Future redesign surface.** If a future improvement targets the
  size-reduce inner loop (e.g., replace schoolbook with `mpn_mul_1`
  inlined for the 1-limb L1 case, or fuse `mul + narrow-sub`), the
  fp machinery is the natural landing spot. The types and bridge are
  already in place.

### What to remove/keep

- **Keep**: `quat_fixed_precision.[hc]`, `quat_mlll_gram_fp`,
  dispatcher, `set_fp_mode`, `gram_fp_equivalence` test,
  `--fp` bench flag, overflow trap.
- **Default**: `g_fp_mode = 0`. Production never opts in.
- **Remove** (future): nothing for now — the scaffolding cost is
  negligible and the guard-rail value is real. Re-evaluate if the
  maintenance burden grows.
