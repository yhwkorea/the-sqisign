/**
 * @file mlll_fp_ops.h
 * @brief MLLL μ/B 실수 백엔드 추상화 (D2a/b/c/d).
 *
 *  Kind:
 *    0 = D2a double        (53-bit mantissa, ±1023 exp)
 *                          ⚠ paper-strict 불가: 53b ≪ Lemma 3 bound (≥506b).
 *                          → μ round 부정확, 종속 검출 실패. 실험/벤치용만.
 *    1 = D2b mpf_t   (GMP, 임의 정밀도, truncation rounding)
 *                          PREC ≥ 2·log p 이면 paper-strict OK.
 *    2 = D2c mpfr_t  (MPFR, correctly-rounded) — **paper-strict 표준**
 *                          PREC ≥ 2·log p 필수 (L5: 1024b).
 *    3 = D2d dpe_t   (53-bit mantissa + int64 exp, fplll style)
 *                          ⚠ paper-strict 불가: mantissa 53b 그대로.
 *                          exp 외부화는 정밀도 부족을 해결하지 못함.
 *
 *  Select at compile time:  -DMLLL_FP_KIND=<0..3>
 *  Precision (mpf/mpfr):    -DMLLL_FP_PREC=<bits>   (default 1024 = L5 paper bound)
 *
 *  Default = KIND=2 PREC=1024 → Compact paper Alg 1 + Appendix A.1 paper-strict.
 *  실측: L1 G=506b, L3 G=769b, L5 G=1016b — paper Lemma mlll-bound와 ≤8b 일치.
 *
 *  Note: wrappers take non-const fp_t because mpz_t/mpf_t/dpe_t are array types
 *  and underlying libraries use non-const signatures. Callers must treat "input"
 *  fp_t args as logically const (wrappers do not mutate them).
 */
#ifndef MLLL_FP_OPS_H
#define MLLL_FP_OPS_H

#include <stdint.h>
#include <math.h>
#include <gmp.h>

#ifndef MLLL_FP_KIND
#define MLLL_FP_KIND 2  /* default: mpfr_t (paper-strict) */
#endif
#ifndef MLLL_FP_PREC
#define MLLL_FP_PREC 1024  /* default: covers L5 paper bound (~1010b) */
#endif

/* ========================================================================== */
#if   MLLL_FP_KIND == 0  /* ===== D2a: double ===== */

typedef double fp_t[1];

static inline void fp_init  (fp_t x)              { x[0] = 0.0; }
static inline void fp_clear (fp_t x)              { (void)x; }
static inline void fp_set   (fp_t x, fp_t y)      { x[0] = y[0]; }
static inline void fp_set_ui(fp_t x, unsigned long u){ x[0] = (double)u; }
static inline void fp_set_si(fp_t x, long s)      { x[0] = (double)s; }
static inline void fp_set_z (fp_t x, const mpz_t z){ x[0] = mpz_get_d(z); }
static inline void fp_add   (fp_t x, fp_t y, fp_t z){ x[0] = y[0] + z[0]; }
static inline void fp_sub   (fp_t x, fp_t y, fp_t z){ x[0] = y[0] - z[0]; }
static inline void fp_mul   (fp_t x, fp_t y, fp_t z){ x[0] = y[0] * z[0]; }
static inline void fp_div   (fp_t x, fp_t y, fp_t z){ x[0] = y[0] / z[0]; }
static inline void fp_neg   (fp_t x, fp_t y)      { x[0] = -y[0]; }
static inline void fp_abs   (fp_t x, fp_t y)      { x[0] = fabs(y[0]); }
static inline int  fp_cmp   (fp_t x, fp_t y)      {
    return (x[0] < y[0]) ? -1 : (x[0] > y[0]) ? 1 : 0;
}
static inline int  fp_sgn   (fp_t x)              {
    return (x[0] < 0.0) ? -1 : (x[0] > 0.0) ? 1 : 0;
}
static inline int  fp_is_zero(fp_t x)             { return x[0] == 0.0; }
static inline void fp_round_z(mpz_t z, fp_t x)  {
    double v = x[0];
    double r = (v >= 0.0) ? floor(v + 0.5) : -floor(-v + 0.5);
    if (isfinite(r)) mpz_set_d(z, r);
    else             mpz_set_ui(z, 0);
}
#define MLLL_FP_NAME "double"

/* ========================================================================== */
#elif MLLL_FP_KIND == 1  /* ===== D2b: mpf_t (GMP) ===== */

typedef mpf_t fp_t;

static inline void fp_init  (fp_t x)              { mpf_init2(x, MLLL_FP_PREC); }
static inline void fp_clear (fp_t x)              { mpf_clear(x); }
static inline void fp_set   (fp_t x, fp_t y)      { mpf_set(x, y); }
static inline void fp_set_ui(fp_t x, unsigned long u){ mpf_set_ui(x, u); }
static inline void fp_set_si(fp_t x, long s)      { mpf_set_si(x, s); }
static inline void fp_set_z (fp_t x, const mpz_t z){ mpf_set_z(x, z); }
static inline void fp_add   (fp_t x, fp_t y, fp_t z){ mpf_add(x, y, z); }
static inline void fp_sub   (fp_t x, fp_t y, fp_t z){ mpf_sub(x, y, z); }
static inline void fp_mul   (fp_t x, fp_t y, fp_t z){ mpf_mul(x, y, z); }
static inline void fp_div   (fp_t x, fp_t y, fp_t z){ mpf_div(x, y, z); }
static inline void fp_neg   (fp_t x, fp_t y)      { mpf_neg(x, y); }
static inline void fp_abs   (fp_t x, fp_t y)      { mpf_abs(x, y); }
static inline int  fp_cmp   (fp_t x, fp_t y)      { return mpf_cmp(x, y); }
static inline int  fp_sgn   (fp_t x)              { return mpf_sgn(x); }
static inline int  fp_is_zero(fp_t x)             { return mpf_sgn(x) == 0; }
static inline void fp_round_z(mpz_t z, fp_t x)  {
    mpf_t t, half;
    mpf_init2(t, MLLL_FP_PREC);
    mpf_init2(half, MLLL_FP_PREC);
    mpf_set_d(half, 0.5);
    if (mpf_sgn(x) >= 0) {
        mpf_add(t, x, half);
        mpf_floor(t, t);
    } else {
        mpf_sub(t, x, half);
        mpf_ceil(t, t);
    }
    mpz_set_f(z, t);
    mpf_clear(t);
    mpf_clear(half);
}
#define MLLL_FP_NAME "mpf_t"

/* ========================================================================== */
#elif MLLL_FP_KIND == 2  /* ===== D2c: mpfr_t ===== */

#include <mpfr.h>
typedef mpfr_t fp_t;

static inline void fp_init  (fp_t x)              { mpfr_init2(x, MLLL_FP_PREC); }
static inline void fp_clear (fp_t x)              { mpfr_clear(x); }
static inline void fp_set   (fp_t x, fp_t y)      { mpfr_set(x, y, MPFR_RNDN); }
static inline void fp_set_ui(fp_t x, unsigned long u){ mpfr_set_ui(x, u, MPFR_RNDN); }
static inline void fp_set_si(fp_t x, long s)      { mpfr_set_si(x, s, MPFR_RNDN); }
static inline void fp_set_z (fp_t x, const mpz_t z){ mpfr_set_z(x, z, MPFR_RNDN); }
static inline void fp_add   (fp_t x, fp_t y, fp_t z){ mpfr_add(x, y, z, MPFR_RNDN); }
static inline void fp_sub   (fp_t x, fp_t y, fp_t z){ mpfr_sub(x, y, z, MPFR_RNDN); }
static inline void fp_mul   (fp_t x, fp_t y, fp_t z){ mpfr_mul(x, y, z, MPFR_RNDN); }
static inline void fp_div   (fp_t x, fp_t y, fp_t z){ mpfr_div(x, y, z, MPFR_RNDN); }
static inline void fp_neg   (fp_t x, fp_t y)      { mpfr_neg(x, y, MPFR_RNDN); }
static inline void fp_abs   (fp_t x, fp_t y)      { mpfr_abs(x, y, MPFR_RNDN); }
static inline int  fp_cmp   (fp_t x, fp_t y)      { return mpfr_cmp(x, y); }
static inline int  fp_sgn   (fp_t x)              { return mpfr_sgn(x); }
static inline int  fp_is_zero(fp_t x)             { return mpfr_zero_p(x); }
static inline void fp_round_z(mpz_t z, fp_t x)  {
    mpfr_get_z(z, x, MPFR_RNDN);
}
#define MLLL_FP_NAME "mpfr_t"

/* ========================================================================== */
#elif MLLL_FP_KIND == 3  /* ===== D2d: dpe_t (double + exp) ===== */

#include "dpe.h"
typedef dpe_t fp_t;

static inline void fp_init  (fp_t x)              {
    dpe_init(x);
    DPE_MANT(x) = 0.0;
    DPE_EXP(x)  = DPE_EXPMIN;
}
static inline void fp_clear (fp_t x)              { dpe_clear(x); }
static inline void fp_set   (fp_t x, fp_t y)      {
    DPE_MANT(x) = DPE_MANT(y); DPE_EXP(x) = DPE_EXP(y);
}
static inline void fp_set_ui(fp_t x, unsigned long u){ dpe_set_ui(x, u); }
static inline void fp_set_si(fp_t x, long s)      { dpe_set_si(x, s); }
static inline void fp_set_z (fp_t x, const mpz_t z){
    /* dpe_set_z takes mpz_t (array), which decays to mpz_ptr. Cast const away. */
    dpe_set_z(x, (mpz_ptr)z);
}
static inline void fp_add   (fp_t x, fp_t y, fp_t z){ dpe_add(x, y, z); }
static inline void fp_sub   (fp_t x, fp_t y, fp_t z){ dpe_sub(x, y, z); }
static inline void fp_mul   (fp_t x, fp_t y, fp_t z){ dpe_mul(x, y, z); }
static inline void fp_div   (fp_t x, fp_t y, fp_t z){ dpe_div(x, y, z); }
static inline void fp_neg   (fp_t x, fp_t y)      { dpe_neg(x, y); }
static inline void fp_abs   (fp_t x, fp_t y)      { dpe_abs(x, y); }
static inline int  fp_cmp   (fp_t x, fp_t y)      { return dpe_cmp(x, y); }
static inline int  fp_sgn   (fp_t x)              { return DPE_SIGN(x); }
static inline int  fp_is_zero(fp_t x)             { return DPE_MANT(x) == 0.0; }
static inline void fp_round_z(mpz_t z, fp_t x)  {
    dpe_t r;
    dpe_init(r);
    dpe_round(r, x);
    dpe_get_z(z, r);
    dpe_clear(r);
}
#define MLLL_FP_NAME "dpe_t"

#else
#error "Unknown MLLL_FP_KIND (expected 0..3)"
#endif

#endif /* MLLL_FP_OPS_H */
