/**
 * @file quat_fixed_precision.c
 * @brief Phase 2 candidate C — stack fixed-precision integer ops.
 *
 * Status:
 *   P2-C-gram step 1 (2026-04-21): add/sub/bitsize + ibz_t bridge implemented.
 *   Remaining stubs (mul, narrow sub_tmp, overflow trap, dpe bridge) are
 *   gated behind fp_stub() so accidental early use is loud — we do not want
 *   a silent zero to impersonate a passing equivalence test during the
 *   `mlll_gram.c` cut-over in P2-C-gram step 3.
 *
 * Two's complement notes:
 *   - Storage is always `digit_t limbs[NWORDS]` with limbs[0] = LSB, MSB
 *     of limbs[nwords - 1] = sign bit. nwords comes from runtime
 *     `quat_fp_widths_t` so a single binary handles L1/L3/L5.
 *   - Add/sub are unsigned limb arithmetic with carry/borrow propagation;
 *     the final carry-out is discarded (2's comp semantics). Sign is
 *     implicit in the top bit; no separate tracking.
 *   - ibz_t → fixed rejects values whose magnitude bitsize is ≥ nwords*64
 *     (need 1 bit for sign). Phase 1 bounds leave a 64-bit margin on top
 *     so this check fires only on an upstream invariant violation.
 */

#include <quaternion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gmp.h>
#include "internal.h"
#include "quat_fixed_precision.h"

/* digit_t and mp_limb_t must be binary-compatible for mpz_import / getlimbn
 * to round-trip without conversion. Both are 64-bit on our target. */
_Static_assert(sizeof(digit_t) == sizeof(mp_limb_t),
               "digit_t must match GMP mp_limb_t size");

/* ---------- Unimplemented trap ---------- */

static void
fp_stub(const char *fname)
{
    fprintf(stderr,
            "quat_fixed_precision: %s is not implemented yet "
            "(Phase 2 P2-C-gram step 2+)\n",
            fname);
    abort();
}

/* ---------- Width dispatch ---------- */

int
quat_fp_widths_from_alg(quat_fp_widths_t *out, const quat_alg_t *alg)
{
    int p_bits = ibz_bitsize(&alg->p);
    if (p_bits <= 128) {
        out->nwords_vec  = NWORDS_QUAT_VEC_L1;
        out->nwords_gram = NWORDS_QUAT_GRAM_L1;
        out->nwords_tmp  = NWORDS_QUAT_TMP_L1;
        return 1;
    } else if (p_bits <= 200) {
        out->nwords_vec  = NWORDS_QUAT_VEC_L3;
        out->nwords_gram = NWORDS_QUAT_GRAM_L3;
        out->nwords_tmp  = NWORDS_QUAT_TMP_L3;
        return 1;
    } else if (p_bits <= 256) {
        out->nwords_vec  = NWORDS_QUAT_VEC_L5;
        out->nwords_gram = NWORDS_QUAT_GRAM_L5;
        out->nwords_tmp  = NWORDS_QUAT_TMP_L5;
        return 1;
    }
    out->nwords_vec = out->nwords_gram = out->nwords_tmp = 0;
    return 0;
}

/* ---------- Raw limb helpers (signed two's complement) ---------- */

/* In-place two's complement negate. */
static void
fp_limbs_negate(digit_t *x, unsigned nwords)
{
    digit_t carry = 1;
    for (unsigned i = 0; i < nwords; i++) {
        digit_t v = ~x[i] + carry;
        /* carry = 1 iff v == 0 after the add (overflow from ~x[i] + 1). */
        carry = (carry != 0 && v == 0) ? 1 : 0;
        x[i] = v;
    }
}

/* r = a + b (discarding final carry; 2's comp). r may alias a or b. */
static void
fp_limbs_add(digit_t *r, const digit_t *a, const digit_t *b, unsigned nwords)
{
    digit_t carry = 0;
    for (unsigned i = 0; i < nwords; i++) {
        digit_t ai = a[i], bi = b[i];
        digit_t s  = ai + carry;
        digit_t c1 = (s < ai) ? 1 : 0;  /* carry from + carry_in */
        digit_t t  = s + bi;
        digit_t c2 = (t < s) ? 1 : 0;   /* carry from + b */
        r[i] = t;
        carry = c1 | c2;
    }
}

/* r = a - b (discarding final borrow; 2's comp). r may alias a or b. */
static void
fp_limbs_sub(digit_t *r, const digit_t *a, const digit_t *b, unsigned nwords)
{
    digit_t borrow = 0;
    for (unsigned i = 0; i < nwords; i++) {
        digit_t ai = a[i], bi = b[i];
        digit_t s  = ai - borrow;
        digit_t b1 = (s > ai) ? 1 : 0;
        digit_t t  = s - bi;
        digit_t b2 = (t > s) ? 1 : 0;
        r[i] = t;
        borrow = b1 | b2;
    }
}

/* Set limbs from ibz_t (two's complement).
 * Returns 0 if |src| does not fit in nwords*64 - 1 bits (sign reserved). */
static int
fp_limbs_set_ibz(digit_t *limbs, unsigned nwords, const ibz_t *src)
{
    int sign = mpz_sgn(*src);
    int mag_bits = ibz_bitsize(src);

    /* mpz_sizeinbase(0, 2) = 1 by GMP convention, but a zero still fits. */
    if (sign != 0 && mag_bits >= (int)(nwords * 64))
        return 0;

    for (unsigned i = 0; i < nwords; i++)
        limbs[i] = 0;

    size_t nlimbs = mpz_size(*src);
    if (nlimbs > (size_t)nwords)
        return 0;  /* defensive; mag_bits check should have caught this */
    for (size_t i = 0; i < nlimbs; i++)
        limbs[i] = (digit_t)mpz_getlimbn(*src, (mp_size_t)i);

    if (sign < 0)
        fp_limbs_negate(limbs, nwords);

    return 1;
}

/* Get ibz_t from limbs (two's complement decode). */
static void
fp_limbs_get_ibz(ibz_t *dst, const digit_t *limbs, unsigned nwords)
{
    digit_t magnitude[NWORDS_QUAT_TMP_MAX] = {0};
    assert(nwords <= NWORDS_QUAT_TMP_MAX);

    int negative = (limbs[nwords - 1] >> 63) & 1;
    if (negative) {
        for (unsigned i = 0; i < nwords; i++)
            magnitude[i] = limbs[i];
        fp_limbs_negate(magnitude, nwords);
    } else {
        for (unsigned i = 0; i < nwords; i++)
            magnitude[i] = limbs[i];
    }

    mpz_import(*dst, nwords, -1 /* LSB first */, sizeof(digit_t),
               0 /* native endian */, 0 /* nails */, magnitude);
    if (negative)
        mpz_neg(*dst, *dst);
}

/* Bitsize of magnitude (matches `ibz_bitsize` = mpz_sizeinbase(_, 2)
 * semantics, which returns 1 for value 0). Cost: one negate + clz. */
static int
fp_limbs_bitsize(const digit_t *limbs, unsigned nwords)
{
    digit_t magnitude[NWORDS_QUAT_TMP_MAX] = {0};
    assert(nwords <= NWORDS_QUAT_TMP_MAX);

    int negative = (limbs[nwords - 1] >> 63) & 1;
    if (negative) {
        for (unsigned i = 0; i < nwords; i++)
            magnitude[i] = limbs[i];
        fp_limbs_negate(magnitude, nwords);
    } else {
        for (unsigned i = 0; i < nwords; i++)
            magnitude[i] = limbs[i];
    }

    for (int i = (int)nwords - 1; i >= 0; i--) {
        if (magnitude[i] != 0) {
            int bits = 64 - __builtin_clzll((unsigned long long)magnitude[i]);
            return i * 64 + bits;
        }
    }
    return 1;  /* zero → 1, matching mpz_sizeinbase(0, 2) */
}

/* ---------- Zero / copy / swap ---------- */

void
quat_fp_vec_set_zero(quat_fp_vec_t *v, const quat_fp_widths_t *w)
{
    memset(v->limbs, 0, w->nwords_vec * sizeof(digit_t));
}

void
quat_fp_gram_set_zero(quat_fp_gram_t *g, const quat_fp_widths_t *w)
{
    memset(g->limbs, 0, w->nwords_gram * sizeof(digit_t));
}

void
quat_fp_vec_copy(quat_fp_vec_t *dst, const quat_fp_vec_t *src,
                 const quat_fp_widths_t *w)
{
    memcpy(dst->limbs, src->limbs, w->nwords_vec * sizeof(digit_t));
}

void
quat_fp_gram_copy(quat_fp_gram_t *dst, const quat_fp_gram_t *src,
                  const quat_fp_widths_t *w)
{
    memcpy(dst->limbs, src->limbs, w->nwords_gram * sizeof(digit_t));
}

void
quat_fp_vec_swap(quat_fp_vec_t *a, quat_fp_vec_t *b,
                 const quat_fp_widths_t *w)
{
    for (unsigned i = 0; i < w->nwords_vec; i++) {
        digit_t t = a->limbs[i];
        a->limbs[i] = b->limbs[i];
        b->limbs[i] = t;
    }
}

void
quat_fp_gram_swap(quat_fp_gram_t *a, quat_fp_gram_t *b,
                  const quat_fp_widths_t *w)
{
    for (unsigned i = 0; i < w->nwords_gram; i++) {
        digit_t t = a->limbs[i];
        a->limbs[i] = b->limbs[i];
        b->limbs[i] = t;
    }
}

int
quat_fp_vec_is_zero(const quat_fp_vec_t *v, const quat_fp_widths_t *w)
{
    for (unsigned i = 0; i < w->nwords_vec; i++)
        if (v->limbs[i] != 0) return 0;
    return 1;
}

int
quat_fp_gram_is_zero(const quat_fp_gram_t *g, const quat_fp_widths_t *w)
{
    for (unsigned i = 0; i < w->nwords_gram; i++)
        if (g->limbs[i] != 0) return 0;
    return 1;
}

/* ---------- Bitsize ---------- */

int
quat_fp_vec_bitsize(const quat_fp_vec_t *v, const quat_fp_widths_t *w)
{
    return fp_limbs_bitsize(v->limbs, w->nwords_vec);
}

int
quat_fp_gram_bitsize(const quat_fp_gram_t *g, const quat_fp_widths_t *w)
{
    return fp_limbs_bitsize(g->limbs, w->nwords_gram);
}

/* ---------- Add / sub ---------- */

void
quat_fp_vec_add(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_add(r->limbs, a->limbs, b->limbs, w->nwords_vec);
}

void
quat_fp_vec_sub(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_sub(r->limbs, a->limbs, b->limbs, w->nwords_vec);
}

void
quat_fp_gram_add(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                 const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_add(r->limbs, a->limbs, b->limbs, w->nwords_gram);
}

void
quat_fp_gram_sub(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                 const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_sub(r->limbs, a->limbs, b->limbs, w->nwords_gram);
}

/* ---------- ibz_t bridge ---------- */

int
quat_fp_vec_set_ibz(quat_fp_vec_t *dst, const ibz_t *src,
                    const quat_fp_widths_t *w)
{
    return fp_limbs_set_ibz(dst->limbs, w->nwords_vec, src);
}

int
quat_fp_gram_set_ibz(quat_fp_gram_t *dst, const ibz_t *src,
                     const quat_fp_widths_t *w)
{
    return fp_limbs_set_ibz(dst->limbs, w->nwords_gram, src);
}

void
quat_fp_vec_get_ibz(ibz_t *dst, const quat_fp_vec_t *src,
                    const quat_fp_widths_t *w)
{
    fp_limbs_get_ibz(dst, src->limbs, w->nwords_vec);
}

void
quat_fp_gram_get_ibz(ibz_t *dst, const quat_fp_gram_t *src,
                     const quat_fp_widths_t *w)
{
    fp_limbs_get_ibz(dst, src->limbs, w->nwords_gram);
}

void
quat_fp_tmp_get_ibz(ibz_t *dst, const quat_fp_tmp_t *src,
                    const quat_fp_widths_t *w)
{
    fp_limbs_get_ibz(dst, src->limbs, w->nwords_tmp);
}

/* ---------- Multiplication (schoolbook, sign-magnitude internally) ----------
 *
 * The hot-loop product in `mlll_gram.c` is always `X * G[i][j]` or
 * `X * b[i][j]`, with operand widths known at compile time. 1st pass here is
 * a straight schoolbook `an × bn → rn` (rn == nwords_tmp) — no Karatsuba, no
 * mpn wrapping. Step 2b of P2-C-gram re-runs a micro-bench to decide whether
 * `mpn_mul` is worth the call overhead; for now the schoolbook keeps the
 * path dependency-free and easy to reason about during the mlll_gram cut-over.
 */

static void
fp_limbs_mul_unsigned(digit_t *r, unsigned rn,
                      const digit_t *a, unsigned an,
                      const digit_t *b, unsigned bn)
{
    assert(rn >= an + bn);
    for (unsigned i = 0; i < rn; i++) r[i] = 0;
    for (unsigned i = 0; i < an; i++) {
        digit_t ai = a[i];
        digit_t carry = 0;
        for (unsigned j = 0; j < bn; j++) {
            __uint128_t prod = (__uint128_t)ai * (__uint128_t)b[j];
            __uint128_t sum  = (__uint128_t)r[i + j] + prod + carry;
            r[i + j] = (digit_t)sum;
            carry    = (digit_t)(sum >> 64);
        }
        r[i + bn] += carry;
    }
}

/* Take magnitude of `src` (an signed 2's comp limbs) into `dst`; returns
 * sign (0 = non-negative, 1 = negative). `dst` must have capacity >= an. */
static int
fp_limbs_magnitude(digit_t *dst, const digit_t *src, unsigned an)
{
    int negative = (src[an - 1] >> 63) & 1;
    for (unsigned i = 0; i < an; i++) dst[i] = src[i];
    if (negative)
        fp_limbs_negate(dst, an);
    return negative;
}

/* r = a * b, signed. rn must hold the full product plus sign room. */
static void
fp_limbs_mul_signed(digit_t *r, unsigned rn,
                    const digit_t *a, unsigned an,
                    const digit_t *b, unsigned bn)
{
    digit_t mag_a[NWORDS_QUAT_TMP_MAX] = {0};
    digit_t mag_b[NWORDS_QUAT_TMP_MAX] = {0};
    assert(an <= NWORDS_QUAT_TMP_MAX && bn <= NWORDS_QUAT_TMP_MAX);

    int sa = fp_limbs_magnitude(mag_a, a, an);
    int sb = fp_limbs_magnitude(mag_b, b, bn);

    fp_limbs_mul_unsigned(r, rn, mag_a, an, mag_b, bn);

    if (sa ^ sb)
        fp_limbs_negate(r, rn);
}

void
quat_fp_tmp_mul_vec_vec(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                        const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_mul_signed(r->limbs, w->nwords_tmp,
                        a->limbs, w->nwords_vec,
                        b->limbs, w->nwords_vec);
}

void
quat_fp_tmp_mul_vec_gram(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                         const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    fp_limbs_mul_signed(r->limbs, w->nwords_tmp,
                        a->limbs, w->nwords_vec,
                        b->limbs, w->nwords_gram);
}

/* ---------- Narrow subtract: target -= tmp (low nwords_target limbs) ------
 *
 * Invariant (enforced by Lemma 3 in callers): the tmp value being subtracted
 * already fits within the narrower width after the subtraction — i.e. its
 * high limbs match the sign extension of the low part. If that invariant
 * is violated, the final result is wrong; P2-overflow adds a runtime trap
 * so bench/equiv runs surface the breach loudly.
 *
 * Mechanically this is the same as `fp_limbs_sub` with the narrow nwords;
 * low `nwords_target` limbs of tmp alias the narrow operand. */

void
quat_fp_vec_sub_tmp(quat_fp_vec_t *r, const quat_fp_tmp_t *t,
                    const quat_fp_widths_t *w)
{
    fp_limbs_sub(r->limbs, r->limbs, t->limbs, w->nwords_vec);
}

void
quat_fp_gram_sub_tmp(quat_fp_gram_t *r, const quat_fp_tmp_t *t,
                     const quat_fp_widths_t *w)
{
    fp_limbs_sub(r->limbs, r->limbs, t->limbs, w->nwords_gram);
}

void
quat_fp_vec_check_overflow(const quat_fp_vec_t *v,
                           const quat_fp_widths_t *w,
                           const char *site)
{
    (void)v; (void)w; (void)site; fp_stub(__func__);
}
void
quat_fp_gram_check_overflow(const quat_fp_gram_t *g,
                            const quat_fp_widths_t *w,
                            const char *site)
{
    (void)g; (void)w; (void)site; fp_stub(__func__);
}

void
quat_fp_gram_to_dpe(dpe_t dst, const quat_fp_gram_t *src,
                    const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__);
}
int
quat_fp_vec_from_dpe_round(quat_fp_vec_t *dst, dpe_t src,
                           const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__); return 0;
}
