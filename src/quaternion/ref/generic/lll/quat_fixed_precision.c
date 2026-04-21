/**
 * @file quat_fixed_precision.c
 * @brief Phase 2 candidate C — stack fixed-precision integer ops.
 *
 * Skeleton only (P2-C-types stage). Arithmetic bodies land in P2-C-gram.
 * Every unimplemented routine traps so that any accidental early use
 * during P2-C wiring is loud — we do not want silent zero results
 * to masquerade as passing equivalence tests.
 */

#include <quaternion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "internal.h"
#include "quat_fixed_precision.h"

/* Single entry point for unimplemented stubs. */
static void
fp_stub(const char *fname)
{
    fprintf(stderr,
            "quat_fixed_precision: %s is not implemented yet "
            "(Phase 2 P2-C-gram)\n",
            fname);
    abort();
}

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

/* ---------- Zeroing / copy / swap (implemented — trivial and needed for
 * basic smoke tests). The rest trap to fp_stub. ---------- */

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

/* ---------- Remaining ops: stubs (P2-C-gram) ---------- */

int
quat_fp_vec_bitsize(const quat_fp_vec_t *v, const quat_fp_widths_t *w)
{
    (void)v; (void)w; fp_stub(__func__); return 0;
}
int
quat_fp_gram_bitsize(const quat_fp_gram_t *g, const quat_fp_widths_t *w)
{
    (void)g; (void)w; fp_stub(__func__); return 0;
}

void
quat_fp_vec_add(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}
void
quat_fp_vec_sub(quat_fp_vec_t *r, const quat_fp_vec_t *a,
                const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}
void
quat_fp_gram_add(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                 const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}
void
quat_fp_gram_sub(quat_fp_gram_t *r, const quat_fp_gram_t *a,
                 const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}

void
quat_fp_tmp_mul_vec_vec(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                        const quat_fp_vec_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}
void
quat_fp_tmp_mul_vec_gram(quat_fp_tmp_t *r, const quat_fp_vec_t *a,
                         const quat_fp_gram_t *b, const quat_fp_widths_t *w)
{
    (void)r; (void)a; (void)b; (void)w; fp_stub(__func__);
}

void
quat_fp_gram_sub_tmp(quat_fp_gram_t *r, const quat_fp_tmp_t *t,
                     const quat_fp_widths_t *w)
{
    (void)r; (void)t; (void)w; fp_stub(__func__);
}
void
quat_fp_vec_sub_tmp(quat_fp_vec_t *r, const quat_fp_tmp_t *t,
                    const quat_fp_widths_t *w)
{
    (void)r; (void)t; (void)w; fp_stub(__func__);
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

int
quat_fp_vec_set_ibz(quat_fp_vec_t *dst, const ibz_t *src,
                    const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__); return 0;
}
int
quat_fp_gram_set_ibz(quat_fp_gram_t *dst, const ibz_t *src,
                     const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__); return 0;
}
void
quat_fp_vec_get_ibz(ibz_t *dst, const quat_fp_vec_t *src,
                    const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__);
}
void
quat_fp_gram_get_ibz(ibz_t *dst, const quat_fp_gram_t *src,
                     const quat_fp_widths_t *w)
{
    (void)dst; (void)src; (void)w; fp_stub(__func__);
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
