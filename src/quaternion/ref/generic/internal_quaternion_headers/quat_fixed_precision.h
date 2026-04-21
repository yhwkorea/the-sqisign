#ifndef QUAT_FIXED_PRECISION_H
#define QUAT_FIXED_PRECISION_H

/**
 * @file quat_fixed_precision.h
 * @brief Phase 2 candidate C — stack-allocated fixed-precision integers
 *        for the Gram-based MLLL inner loop.
 *
 * Rationale:
 *   Phase 1 empirically established per-level bounds (L1 vec 259 / Gram 518,
 *   L3 vec 391 / Gram 782, L5 vec 513 / Gram 1026 bits; see §5.5 of
 *   FIXED_PRECISION_DECISION_KO.md). Candidate B (`mpz_realloc2` hint) failed
 *   to translate those bounds into heap-allocation reduction — GMP treats the
 *   hint only as a seed capacity, so `ibz_mul`/`ibz_add` still auto-realloc
 *   on demand. Candidate C replaces the `ibz_t` storage with raw digit_t
 *   arrays sized once, on the stack, at the top of `quat_mlll_gram`.
 *
 * Sizing policy:
 *   Arrays are dimensioned to the LARGEST level (L5). The runtime descriptor
 *   `quat_fp_widths_t` records the actual number of limbs in use for the
 *   current call — arithmetic routines operate on `nwords_*` not the physical
 *   array length. This keeps a single binary compatible with all three levels
 *   (matching the existing `mlll_gram_level_hints` pattern), at the cost of
 *   carrying up to ~200 bytes of unused stack at L1.
 *
 *   Margin (+64 bits / +1 limb) is included so size-reduce intermediate
 *   `X*G[i][j]` can fit in tmp before being subtracted back to the Gram
 *   bound. Overflow beyond the reserved margin is a program error and
 *   triggers `abort()` (see `quat_fp_vec_check_overflow` / gram variant).
 *
 * Representation:
 *   Two's complement, limb[0] = LSB, sign bit = MSB of the most significant
 *   in-use limb (limb[nwords_* - 1]). Multiplication internally splits sign
 *   and operates unsigned, then reapplies sign — this keeps `mp_mul` / future
 *   `mpn_mul` wrappers usable without re-implementing signed multiplication.
 *
 * Thread-safety:
 *   All structs are pass-by-pointer value types. No global state. Callers
 *   pass the widths descriptor explicitly.
 *
 * Scope:
 *   Internal-only header. Not exported via any public API. Consumed
 *   exclusively by `lll/mlll_gram.c` once Phase 2 C path lands.
 */

#include <quaternion.h>
#include <stdint.h>
#include "dpe.h"

/* digit_t is defined in src/gf/broadwell/lvl?/include/gf?????.h — pulled in
 * transitively via <quaternion.h> -> <sqisign.h> -> <tutil.h>. Asserted here
 * so this header fails fast if the chain ever changes. */
#ifndef RADIX
#  error "quat_fixed_precision.h requires RADIX (from tutil.h) to be defined"
#endif
#if RADIX != 64
#  error "quat_fixed_precision.h is only designed for 64-bit digit_t (RADIX=64)"
#endif

/* ---------- Per-level widths (from Phase 1, +64-bit margin, u64-aligned) ----- */

/* vec coord (Lemma 1 bound): L1 259, L3 391, L5 513 bits. */
#define NWORDS_QUAT_VEC_L1    5   /* 320 bits */
#define NWORDS_QUAT_VEC_L3    7   /* 448 bits */
#define NWORDS_QUAT_VEC_L5    9   /* 576 bits */
#define NWORDS_QUAT_VEC_MAX   NWORDS_QUAT_VEC_L5

/* Gram entry (Lemma 3 bound, 2 * max <a_k, a_k>): L1 518, L3 782, L5 1026. */
#define NWORDS_QUAT_GRAM_L1   9   /* 576 bits */
#define NWORDS_QUAT_GRAM_L3   13  /* 832 bits */
#define NWORDS_QUAT_GRAM_L5   17  /* 1088 bits */
#define NWORDS_QUAT_GRAM_MAX  NWORDS_QUAT_GRAM_L5

/* tmp = X * G[i][j] (vec limbs + Gram limbs + carry). */
#define NWORDS_QUAT_TMP_L1    (NWORDS_QUAT_VEC_L1 + NWORDS_QUAT_GRAM_L1 + 1)
#define NWORDS_QUAT_TMP_L3    (NWORDS_QUAT_VEC_L3 + NWORDS_QUAT_GRAM_L3 + 1)
#define NWORDS_QUAT_TMP_L5    (NWORDS_QUAT_VEC_L5 + NWORDS_QUAT_GRAM_L5 + 1)
#define NWORDS_QUAT_TMP_MAX   NWORDS_QUAT_TMP_L5

/* ---------- Storage types ---------- */

typedef struct {
    digit_t limbs[NWORDS_QUAT_VEC_MAX];
} quat_fp_vec_t;

typedef struct {
    digit_t limbs[NWORDS_QUAT_GRAM_MAX];
} quat_fp_gram_t;

typedef struct {
    digit_t limbs[NWORDS_QUAT_TMP_MAX];
} quat_fp_tmp_t;

/* ---------- Runtime widths descriptor ---------- */

typedef struct {
    unsigned nwords_vec;    /* L1 5 / L3 7 / L5 9 */
    unsigned nwords_gram;   /* L1 9 / L3 13 / L5 17 */
    unsigned nwords_tmp;    /* L1 15 / L3 21 / L5 27 */
} quat_fp_widths_t;

/**
 * Resolve per-level widths from `alg->p` bitsize. Mirrors the runtime
 * dispatch used by `mlll_gram_level_hints` (mlll_gram.c) so B and C paths
 * agree on the same tier for equivalence testing.
 *
 * @return 1 on match (L1/L3/L5), 0 otherwise (caller must fall back to
 * ibz_t path; widths left zero).
 */
int quat_fp_widths_from_alg(quat_fp_widths_t *out, const quat_alg_t *alg);

/* ---------- Element ops (declarations only; implemented in P2-C-gram) -------
 *
 * Naming convention:
 *   quat_fp_<storage>_<op>_... where storage in {vec, gram, tmp}.
 *
 * All ops take explicit widths from quat_fp_widths_t so a single binary
 * serves all levels. Widths are never read from the struct — they don't
 * fit there. That is deliberate: pass-by-value struct is cheap, and the
 * hot loop threads a single `widths` through anyway.
 */

/* Zeroing / copy / swap. */
void quat_fp_vec_set_zero(quat_fp_vec_t *v, const quat_fp_widths_t *w);
void quat_fp_gram_set_zero(quat_fp_gram_t *g, const quat_fp_widths_t *w);
void quat_fp_vec_copy(quat_fp_vec_t *dst, const quat_fp_vec_t *src,
                      const quat_fp_widths_t *w);
void quat_fp_gram_copy(quat_fp_gram_t *dst, const quat_fp_gram_t *src,
                       const quat_fp_widths_t *w);
void quat_fp_vec_swap(quat_fp_vec_t *a, quat_fp_vec_t *b,
                      const quat_fp_widths_t *w);
void quat_fp_gram_swap(quat_fp_gram_t *a, quat_fp_gram_t *b,
                       const quat_fp_widths_t *w);

/* Predicates. */
int quat_fp_vec_is_zero(const quat_fp_vec_t *v, const quat_fp_widths_t *w);
int quat_fp_gram_is_zero(const quat_fp_gram_t *g, const quat_fp_widths_t *w);
/* Bitsize helpers for the tracker (Phase 1 instrumentation compat). */
int quat_fp_vec_bitsize(const quat_fp_vec_t *v, const quat_fp_widths_t *w);
int quat_fp_gram_bitsize(const quat_fp_gram_t *g, const quat_fp_widths_t *w);

/* Add / sub. Two's complement with carry; MSB margin must survive. */
void quat_fp_vec_add(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                     const quat_fp_vec_t *b, const quat_fp_widths_t *w);
void quat_fp_vec_sub(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                     const quat_fp_vec_t *b, const quat_fp_widths_t *w);
void quat_fp_gram_add(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                      const quat_fp_gram_t *b, const quat_fp_widths_t *w);
void quat_fp_gram_sub(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                      const quat_fp_gram_t *b, const quat_fp_widths_t *w);

/* Multiplications (tmp = vec * vec, tmp = vec * gram). Sign handled
 * internally (sign-magnitude split → unsigned mul → sign combine). */
void quat_fp_tmp_mul_vec_vec(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                             const quat_fp_vec_t *b,
                             const quat_fp_widths_t *w);
void quat_fp_tmp_mul_vec_gram(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                              const quat_fp_gram_t *b,
                              const quat_fp_widths_t *w);

/* Narrowing: gram -= tmp (tmp must fit in gram width after the subtraction). */
void quat_fp_gram_sub_tmp(quat_fp_gram_t *r, const quat_fp_tmp_t *t,
                          const quat_fp_widths_t *w);
/* vec -= tmp (same discipline). */
void quat_fp_vec_sub_tmp(quat_fp_vec_t *r, const quat_fp_tmp_t *t,
                         const quat_fp_widths_t *w);

/* Overflow guard — fires abort() if the MSB + margin invariant is violated. */
void quat_fp_vec_check_overflow(const quat_fp_vec_t *v,
                                const quat_fp_widths_t *w,
                                const char *site);
void quat_fp_gram_check_overflow(const quat_fp_gram_t *g,
                                 const quat_fp_widths_t *w,
                                 const char *site);

/* ---------- Conversions (for B vs C equivalence oracle) ---------- */

/* ibz_t -> fixed. Returns 0 if src does not fit within the target width;
 * caller must treat that as a bug (not a normal condition post-Phase 1). */
int quat_fp_vec_set_ibz(quat_fp_vec_t *dst, const ibz_t *src,
                        const quat_fp_widths_t *w);
int quat_fp_gram_set_ibz(quat_fp_gram_t *dst, const ibz_t *src,
                         const quat_fp_widths_t *w);

/* fixed -> ibz_t (for tracker / equivalence dump). */
void quat_fp_vec_get_ibz(ibz_t *dst, const quat_fp_vec_t *src,
                         const quat_fp_widths_t *w);
void quat_fp_gram_get_ibz(ibz_t *dst, const quat_fp_gram_t *src,
                          const quat_fp_widths_t *w);

/* ---------- dpe bridge (for L² Cholesky computed from the fixed Gram) ------ */

/* Set a dpe from a fixed-width Gram entry (signed). */
void quat_fp_gram_to_dpe(dpe_t dst, const quat_fp_gram_t *src,
                         const quat_fp_widths_t *w);
/* Round a dpe to a signed integer X that fits in vec width; used for the
 * size-reduce multiplier. Returns 0 if it overflows the vec width. */
int quat_fp_vec_from_dpe_round(quat_fp_vec_t *dst, dpe_t src,
                               const quat_fp_widths_t *w);

#endif /* QUAT_FIXED_PRECISION_H */
