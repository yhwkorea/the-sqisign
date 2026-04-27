/**
 * @file mlll.c
 * @brief Modified LLL (MLLL) — Compact paper Algorithm 1, paper-strict + float μ/B.
 *
 * Paper-strict representation (D2 rewrite, 2026-04-24):
 *   b[i] ∈ ℤ^4           — basis vector (integer, Lemma 3 bound ‖a‖)
 *   G[i][j] = <b[i],b[j]> — Gram matrix (integer, Lemma 3 bound ‖a‖²)
 *   μ[i][j]               — real-valued (fp_t, 4 backends via mlll_fp_ops.h)
 *   B[j]    = ‖b[j]*‖²    — real-valued (fp_t)
 *
 * Backend selector (compile-time):
 *   -DMLLL_FP_KIND=0 → double           (53-bit, ±1023 exp; L5 overflow expected)
 *   -DMLLL_FP_KIND=1 → mpf_t (GMP)      (MLLL_FP_PREC bits, default 128)
 *   -DMLLL_FP_KIND=2 → mpfr_t           (correctly-rounded)
 *   -DMLLL_FP_KIND=3 → dpe_t            (double + int64 exp, fplll style) [default]
 *
 * Rationale: Paper Appendix A.1 explicitly states μ is "real-valued" and does
 * not introduce new unbounded integers. Previous integer-rational (D1) variant
 * inflated intermediate to ~6·p_bits due to 3-term product (μ²·B) blowup.
 * With float μ/B, the Lemma 3 bound ‖a‖² holds literally for (b, G, h).
 */

#include <quaternion.h>
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"
#include "mlll_fp_ops.h"

#define N MLLL_MAX_GENERATORS

/* ========== ibz vector helpers ========== */

/* Quaternion norm bilinear form: <a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3 */
static void
ibz_vec_4_dot(ibz_t *dot, const ibz_vec_4_t *a, const ibz_vec_4_t *b,
              const ibz_t *p)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_set(dot, 0);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&tmp, &((*a)[i]), &((*b)[i]));
        if (i >= 2)
            ibz_mul(&tmp, &tmp, p);
        ibz_add(dot, dot, &tmp);
    }
    ibz_finalize(&tmp);
}

static void
ibz_vec_4_sub_scalar_mul(ibz_vec_4_t *b, const ibz_t *t, const ibz_vec_4_t *a)
{
    ibz_t tmp;
    ibz_init(&tmp);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&tmp, t, &((*a)[i]));
        ibz_sub(&((*b)[i]), &((*b)[i]), &tmp);
    }
    ibz_finalize(&tmp);
}

static void
ibz_vec_4_swap(ibz_vec_4_t *a, ibz_vec_4_t *b)
{
    for (int i = 0; i < 4; i++)
        ibz_swap(&((*a)[i]), &((*b)[i]));
}

/* ========== Bitsize tracking ========== */

static int max_bits_vec = 0;
static int max_bits_gram = 0;
static int max_bits_intermediate = 0;
static int mlll_call_count = 0;

static inline void track_intermediate(const ibz_t *v) {
    int b = ibz_bitsize(v);
    if (b > max_bits_intermediate) max_bits_intermediate = b;
    tracker_update_gso_ibz(v);
}

static inline void track_gram(const ibz_t *v) {
    int b = ibz_bitsize(v);
    if (b > max_bits_gram) max_bits_gram = b;
}

void mlll_print_max_bits(void) {
    fprintf(stderr, "[MLLL bitsize, fp=%s prec=%d] calls=%d  "
            "vec=%d bits (%d u64)  G=%d bits (%d u64)  "
            "intermediate=%d bits (%d u64)\n",
            MLLL_FP_NAME, MLLL_FP_PREC,
            mlll_call_count,
            max_bits_vec, (max_bits_vec + 63) / 64,
            max_bits_gram, (max_bits_gram + 63) / 64,
            max_bits_intermediate, (max_bits_intermediate + 63) / 64);
}

/* ========== Gram / μ / B maintenance ========== */

static void
compute_gram_row(int idx, ibz_vec_4_t b[], ibz_t G[][N], const ibz_t *p)
{
    for (int j = 0; j <= idx; j++) {
        ibz_vec_4_dot(&G[idx][j], &b[idx], &b[j], p);
        track_gram(&G[idx][j]);
        if (j != idx)
            ibz_copy(&G[j][idx], &G[idx][j]);
    }
}

/* Compute μ[idx][j] for 0 ≤ j < idx, and B[idx].  Uses paper line 2 formula:
 *   μ[β][j] = (G[β][j] - Σ_{k<j} μ[β][k]·μ[j][k]·B[k]) / B[j]
 *   B[β]    = G[β][β] - Σ_{j<β} μ[β][j]²·B[j]
 */
static void
compute_mu_B_row(int idx, ibz_t G[][N], fp_t mu[][N], fp_t B[])
{
    fp_t acc, term, term2;
    fp_init(acc); fp_init(term); fp_init(term2);

    for (int j = 0; j < idx; j++) {
        fp_set_z(acc, G[idx][j]);  /* acc = G[idx][j] */

        for (int k = 0; k < j; k++) {
            if (fp_is_zero(B[k]))
                continue;
            /* term = μ[idx][k] · μ[j][k] · B[k] */
            fp_mul(term, mu[idx][k], mu[j][k]);
            fp_mul(term, term, B[k]);
            fp_sub(acc, acc, term);
        }

        if (fp_is_zero(B[j])) {
            fp_set_ui(mu[idx][j], 0);
        } else {
            fp_div(mu[idx][j], acc, B[j]);
        }
    }

    /* B[idx] = G[idx][idx] - Σ μ[idx][j]² · B[j] */
    fp_set_z(B[idx], G[idx][idx]);
    for (int j = 0; j < idx; j++) {
        if (fp_is_zero(B[j]))
            continue;
        fp_mul(term,  mu[idx][j], mu[idx][j]);
        fp_mul(term2, term, B[j]);
        fp_sub(B[idx], B[idx], term2);
    }

    fp_clear(acc); fp_clear(term); fp_clear(term2);
}

/* Recompute the entire GSO (rows 0..beta-1).  Also refills G. */
static void
recompute_all(int beta, ibz_vec_4_t b[], ibz_t G[][N],
              fp_t mu[][N], fp_t B[], const ibz_t *p)
{
    for (int i = 0; i < beta; i++)
        compute_gram_row(i, b, G, p);
    for (int i = 0; i < beta; i++)
        compute_mu_B_row(i, G, mu, B);
}

/* ========== MLLL core — paper Alg 1 ========== */

void
quat_mlll(ibz_mat_4x4_t *basis,
          int *rank,
          const ibz_vec_4_t *generators,
          int g,
          const quat_alg_t *alg)
{
    assert(g >= 1 && g <= N);

    ibz_vec_4_t b[N];
    ibz_t G[N][N];
    fp_t mu[N][N];
    fp_t B[N];

    for (int i = 0; i < g; i++) {
        ibz_vec_4_init(&b[i]);
        for (int j = 0; j < g; j++) {
            ibz_init(&G[i][j]);
            fp_init(mu[i][j]);
        }
        fp_init(B[i]);
    }

    ibz_t t;
    fp_t mu_abs, half, mu_sq, three_quarters, rhs, t_fp, tmp_fp;
    ibz_init(&t);
    fp_init(mu_abs); fp_init(half); fp_init(mu_sq);
    fp_init(three_quarters); fp_init(rhs); fp_init(t_fp); fp_init(tmp_fp);
    fp_set_ui(half, 1); { fp_t two; fp_init(two); fp_set_ui(two, 2); fp_div(half, half, two); fp_clear(two); }
    fp_set_ui(three_quarters, 3); { fp_t four; fp_init(four); fp_set_ui(four, 4); fp_div(three_quarters, three_quarters, four); fp_clear(four); }

    int alpha = 0, beta = 0, tau = 2, m, l;

    /* ===== LOAD ===== */
load:
    if (alpha >= g)
        goto start_reduce;

    while (alpha < g && ibz_vec_4_is_zero(&generators[alpha]))
        alpha++;
    if (alpha >= g)
        goto start_reduce;

    for (int i = 0; i < 4; i++)
        ibz_copy(&(b[beta][i]), &(generators[alpha][i]));
    alpha++;

    compute_gram_row(beta, b, G, &alg->p);
    compute_mu_B_row(beta, G, mu, B);
    beta++;

    if (!fp_is_zero(B[beta - 1]) && alpha < g)
        goto load;

start_reduce:
    if (beta <= 1)
        goto done;

    m = tau;
    if (m > beta) m = beta;
    if (m < 2)    m = 2;

    /* ===== REDUCTION ===== */
reduction:
    recompute_all(beta, b, G, mu, B, &alg->p);

    l = m - 1;

size_reduce:
    /* |μ[m-1][l-1]| > 1/2 ? */
    if (l >= 1 && !fp_is_zero(B[l - 1]) && !fp_is_zero(mu[m - 1][l - 1])) {
        fp_abs(mu_abs, mu[m - 1][l - 1]);
        if (fp_cmp(mu_abs, half) > 0) {
            /* t = round(μ[m-1][l-1]) */
            fp_round_z(t, mu[m - 1][l - 1]);

            if (!ibz_is_zero(&t)) {
                /* b[m-1] -= t · b[l-1] */
                ibz_vec_4_sub_scalar_mul(&b[m - 1], &t, &b[l - 1]);
                tracker_update_vec4(&b[m - 1]);

                /* Update G (ibz) for row/col m-1. */
                ibz_t gtmp;
                ibz_init(&gtmp);
                for (int j = 0; j < beta; j++) {
                    if (j == m - 1) continue;
                    ibz_mul(&gtmp, &t, &G[l - 1][j]);
                    track_intermediate(&gtmp);
                    ibz_sub(&G[m - 1][j], &G[m - 1][j], &gtmp);
                    ibz_copy(&G[j][m - 1], &G[m - 1][j]);
                    track_gram(&G[m - 1][j]);
                }
                /* Recompute G[m-1][m-1] from b to avoid 2t·G_ml + t²·G_ll blowup. */
                ibz_vec_4_dot(&G[m - 1][m - 1], &b[m - 1], &b[m - 1], &alg->p);
                track_gram(&G[m - 1][m - 1]);
                ibz_finalize(&gtmp);

                /* μ row update (paper line 9, float arithmetic):
                 *   μ[m-1][l-1] -= t
                 *   μ[m-1][j]   -= t·μ[l-1][j]  for j < l-1
                 */
                fp_set_z(t_fp, t);
                fp_sub(mu[m - 1][l - 1], mu[m - 1][l - 1], t_fp);
                for (int j = 0; j < l - 1; j++) {
                    fp_mul(tmp_fp, t_fp, mu[l - 1][j]);
                    fp_sub(mu[m - 1][j], mu[m - 1][j], tmp_fp);
                }
            }
        }
    }

    if (ibz_vec_4_is_zero(&b[m - 1]))
        goto remove_vector;

    /* Dependent but nonzero: swap down. */
    if (fp_is_zero(B[m - 1])) {
        if (m >= 2) {
            ibz_vec_4_swap(&b[m - 1], &b[m - 2]);
            if (m > 2)
                m = m - 1;
        }
        goto reduction;
    }

    if (l < m - 1)
        goto next_l;

    /* Lovász condition: B[m-1] < (3/4 - μ[m-1][m-2]²) · B[m-2] */
    if (m >= 2 && !fp_is_zero(B[m - 2])) {
        fp_mul(mu_sq, mu[m - 1][m - 2], mu[m - 1][m - 2]);
        fp_sub(rhs, three_quarters, mu_sq);
        fp_mul(rhs, rhs, B[m - 2]);

        if (fp_cmp(B[m - 1], rhs) < 0)
            goto do_swap;
    }

next_l:
    l = l - 1;
    if (l >= 1)
        goto size_reduce;
    m = m + 1;
    if (m > beta) {
        if (alpha < g) {
            tau = m;
            goto load;
        }
        goto done;
    }
    goto reduction;

    /* ===== SWAP ===== */
do_swap:
    ibz_vec_4_swap(&b[m - 1], &b[m - 2]);
    tracker_update_vec4(&b[m - 1]);
    tracker_update_vec4(&b[m - 2]);
    if (m > 2)
        m = m - 1;
    goto reduction;

    /* ===== REMOVE ===== */
remove_vector:
    for (int i = m; i <= beta - 1; i++) {
        for (int c = 0; c < 4; c++)
            ibz_copy(&(b[i - 1][c]), &(b[i][c]));
    }
    beta--;

    if (alpha >= g) {
        if (beta <= 1)
            goto done;
        recompute_all(beta, b, G, mu, B, &alg->p);
        m = 2;
        if (m > beta)
            goto done;
        goto reduction;
    }

    recompute_all(beta, b, G, mu, B, &alg->p);
    tau = m;
    if (tau < 2) tau = 2;
    goto load;

done:
    tracker_disable();

    mlll_call_count++;
    int nz_count = 0;
    for (int i = 0; i < beta; i++) {
        if (!ibz_vec_4_is_zero(&b[i])) nz_count++;
        for (int c = 0; c < 4; c++) {
            int bs = ibz_bitsize(&(b[i][c]));
            if (bs > max_bits_vec) max_bits_vec = bs;
        }
    }

    /* Initialize output basis */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&((*basis)[i][j]), 0);

    /* Paper Alg 1: MLLL must produce a clean rank-≤4 basis (all dependent
     * vectors size-reduced to zero and removed).  No post-processing. */
    assert(nz_count <= 4);
    *rank = 0;
    for (int i = 0; i < beta && *rank < 4; i++) {
        if (!ibz_vec_4_is_zero(&b[i])) {
            for (int j = 0; j < 4; j++)
                ibz_copy(&((*basis)[j][*rank]), &(b[i][j]));
            (*rank)++;
        }
    }

    /* Cleanup */
    for (int i = 0; i < g; i++) {
        ibz_vec_4_finalize(&b[i]);
        for (int j = 0; j < g; j++) {
            ibz_finalize(&G[i][j]);
            fp_clear(mu[i][j]);
        }
        fp_clear(B[i]);
    }
    ibz_finalize(&t);
    fp_clear(mu_abs); fp_clear(half); fp_clear(mu_sq);
    fp_clear(three_quarters); fp_clear(rhs); fp_clear(t_fp); fp_clear(tmp_fp);
}

/* ========== Lattice operations using MLLL ========== */

void
quat_lattice_mul_mlll(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg)
{
    ibz_vec_4_t elem1, elem2, elem_res;
    ibz_vec_4_t generators[16];
    int rank;

    ibz_vec_4_init(&elem1);
    ibz_vec_4_init(&elem2);
    ibz_vec_4_init(&elem_res);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_init(&generators[i]);

    ibz_mat_4x4_t red1, red2;
    ibz_mat_4x4_init(&red1);
    ibz_mat_4x4_init(&red2);

    quat_lattice_lll(&red1, lat1, alg);
    quat_lattice_lll(&red2, lat2, alg);

    for (int k = 0; k < 4; k++) {
        ibz_vec_4_copy_ibz(
            &elem1, &(red1[0][k]), &(red1[1][k]),
            &(red1[2][k]), &(red1[3][k]));
        for (int i = 0; i < 4; i++) {
            ibz_vec_4_copy_ibz(
                &elem2, &(red2[0][i]), &(red2[1][i]),
                &(red2[2][i]), &(red2[3][i]));
            quat_alg_coord_mul(&elem_res, &elem1, &elem2, alg);
            for (int j = 0; j < 4; j++)
                ibz_copy(&(generators[4 * k + i][j]), &(elem_res[j]));
        }
    }

    quat_mlll(&(res->basis), &rank, generators, 16, alg);
    assert(rank == 4);

    ibz_mul(&(res->denom), &(lat1->denom), &(lat2->denom));
    quat_lattice_reduce_denom(res, res);

    ibz_mat_4x4_finalize(&red1);
    ibz_mat_4x4_finalize(&red2);
    ibz_vec_4_finalize(&elem1);
    ibz_vec_4_finalize(&elem2);
    ibz_vec_4_finalize(&elem_res);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_finalize(&generators[i]);
}

void
quat_lattice_add_mlll(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg)
{
    ibz_vec_4_t generators[8];
    ibz_mat_4x4_t tmp;
    int rank;

    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);
    ibz_mat_4x4_init(&tmp);

    ibz_mat_4x4_scalar_mul(&tmp, &(lat1->denom), &(lat2->basis));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_copy(&(generators[j][i]), &(tmp[i][j]));

    ibz_mat_4x4_scalar_mul(&tmp, &(lat2->denom), &(lat1->basis));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_copy(&(generators[4 + j][i]), &(tmp[i][j]));

    quat_mlll(&(res->basis), &rank, generators, 8, alg);
    assert(rank > 0 && rank <= 4);

    ibz_mul(&(res->denom), &(lat1->denom), &(lat2->denom));
    quat_lattice_reduce_denom(res, res);

    ibz_mat_4x4_finalize(&tmp);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
}
