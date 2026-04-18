/**
 * @file mlll.c
 * @brief Modified LLL (MLLL) — Compact paper Algorithm 1
 *
 * Integral GSO (Cohen 2.6.3): d[k], Lambda[i][j] — all integer, no fractions.
 *   d[0] = 1
 *   Lambda[i][j] = d[j+1] * mu[i][j]  (integer)
 *   d[k+1] = d[k] * B[k]              (integer)
 * Vector coordinates b[i] stay bounded by max ||a_i||^2 (Lemma 3).
 */

#include <quaternion.h>
#include <stdio.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"

#define N MLLL_MAX_GENERATORS

/* ---------- helpers ---------- */

/* Quaternion norm bilinear form: <a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3
 * This matches quat_lll_bilinear used by the existing LLL and verification. */
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

/* (no row_gcd_reduce / B_reduce needed — integral GSO uses only integers) */

/* ========== Bitsize tracking (temporary) ========== */

static int max_bits_vec = 0;
static int max_bits_d = 0;
static int max_bits_lambda = 0;
static int max_bits_intermediate = 0;
static int mlll_call_count = 0;

static inline void track_intermediate(const ibz_t *v) {
    int b = ibz_bitsize(v);
    if (b > max_bits_intermediate) max_bits_intermediate = b;
    tracker_update_gso_ibz(v);
}

void mlll_print_max_bits(void) {
    fprintf(stderr, "[MLLL bitsize] calls=%d  vec=%d bits (%d u64)  "
            "d=%d bits (%d u64)  lambda=%d bits (%d u64)  "
            "intermediate=%d bits (%d u64)\n",
            mlll_call_count,
            max_bits_vec, (max_bits_vec + 63) / 64,
            max_bits_d, (max_bits_d + 63) / 64,
            max_bits_lambda, (max_bits_lambda + 63) / 64,
            max_bits_intermediate, (max_bits_intermediate + 63) / 64);
}

/* ========== Integral GSO row computation (Cohen 2.6.3) ========== */

/**
 * @brief Compute Lambda[idx][j] for j < idx, and d[idx+1].
 *
 * Cohen's Algorithm 2.6.3 (integral GSO, incremental update):
 *   d[0] = 1,  d[k] = det(Gram of b[0]..b[k-1])
 *   Lambda[i][j] = d[j+1] * mu[i][j]  (integer)
 *
 * The Lambda and d values are computed incrementally:
 *   Lambda[idx][j] starts as <b[idx], b[j]>, then transformed for k=0..j-1
 *   d[idx+1] starts as <b[idx], b[idx]>, then updated for each j
 *
 * d[idx+1] = 0 means b[idx] is linearly dependent on b[0]..b[idx-1].
 */
static void
compute_gso_row(int idx, const ibz_vec_4_t *b,
                ibz_t Lambda[][N], ibz_t d[], const ibz_t *p)
{
    ibz_t tmp, prod, rem;
    ibz_init(&tmp); ibz_init(&prod); ibz_init(&rem);

    /* d[idx+1] starts as <b[idx], b[idx]> */
    ibz_vec_4_dot(&d[idx + 1], &b[idx], &b[idx], p);
    track_intermediate(&d[idx + 1]);

    for (int j = 0; j < idx; j++) {
        /* Lambda[idx][j] starts as <b[idx], b[j]> */
        ibz_vec_4_dot(&Lambda[idx][j], &b[idx], &b[j], p);

        /* Transform: Lambda[idx][j] = (d[k+1]*Lambda[idx][j] - Lambda[idx][k]*Lambda[j][k]) / d[k]
         * for k = 0 to j-1 */
        for (int k = 0; k < j; k++) {
            ibz_mul(&tmp, &d[k + 1], &Lambda[idx][j]);
            ibz_mul(&prod, &Lambda[idx][k], &Lambda[j][k]);
            ibz_sub(&tmp, &tmp, &prod);
            track_intermediate(&tmp);
            ibz_div(&Lambda[idx][j], &rem, &tmp, &d[k]);
            assert(ibz_is_zero(&rem));
        }
        track_intermediate(&Lambda[idx][j]);

        /* Update d[idx+1] = (d[j+1]*d[idx+1] - Lambda[idx][j]^2) / d[j] */
        if (!ibz_is_zero(&d[j])) {
            ibz_mul(&tmp, &d[j + 1], &d[idx + 1]);
            ibz_mul(&prod, &Lambda[idx][j], &Lambda[idx][j]);
            ibz_sub(&tmp, &tmp, &prod);
            track_intermediate(&tmp);
            ibz_div(&d[idx + 1], &rem, &tmp, &d[j]);
            assert(ibz_is_zero(&rem));
        }
        /* If d[j]=0, b[j] was dependent — shouldn't happen if b[0..idx-1] are independent */
        track_intermediate(&d[idx + 1]);
    }

    ibz_finalize(&tmp); ibz_finalize(&prod); ibz_finalize(&rem);
}

/* ========== MLLL core (integral GSO) ========== */

void
quat_mlll(ibz_mat_4x4_t *basis,
          int *rank,
          const ibz_vec_4_t *generators,
          int g,
          const quat_alg_t *alg)
{
    assert(g >= 1 && g <= N);

    ibz_vec_4_t b[N];
    ibz_t Lambda[N][N];  /* Lambda[i][j] = d[j+1] * mu[i][j], integer */
    ibz_t d[N + 1];      /* d[0]=1, d[k+1] = d[k]*B[k], integer */

    for (int i = 0; i < g; i++) {
        ibz_vec_4_init(&b[i]);
        for (int j = 0; j < g; j++)
            ibz_init(&Lambda[i][j]);
    }
    for (int i = 0; i <= g; i++) {
        ibz_init(&d[i]);
    }
    ibz_set(&d[0], 1);

    ibz_t t, tmp1, tmp2, rem;
    ibz_init(&t);
    ibz_init(&tmp1);
    ibz_init(&tmp2);
    ibz_init(&rem);

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

    compute_gso_row(beta, b, Lambda, d, &alg->p);
    beta++;

    if (!ibz_is_zero(&d[beta]) && alpha < g)
        goto load;

start_reduce:
    if (beta <= 1)
        goto done;

    m = tau;
    if (m > beta)
        m = beta;
    if (m < 2)
        m = 2;

    /* ===== REDUCTION (1-indexed m, 0-indexed arrays) ===== */
reduction:
    /* Recompute GSO from 0 to m-1 to ensure all Lambda/d are fresh.
     * Necessary because swaps (both Lovász and dependent-vector) can leave
     * d values stale. Safe because b[0]..b[m-2] are always independent
     * (dependent vectors are at position m-1, not earlier).
     * Cost: O(m^3) per step, acceptable for dimension 4 / beta ≤ 16. */
    ibz_set(&d[0], 1);
    for (int _i = 0; _i < m; _i++)
        compute_gso_row(_i, b, Lambda, d, &alg->p);

    l = m - 1;

size_reduce:
    /* |mu[m][l]| > 1/2 ?
     * mu[m-1][l-1] = Lambda[m-1][l-1] / d[l]
     * |Lambda[m-1][l-1]| * 2 > d[l] ?
     */
    if (l >= 1 && !ibz_is_zero(&d[l])) {
        ibz_abs(&tmp1, &Lambda[m - 1][l - 1]);
        ibz_add(&tmp1, &tmp1, &tmp1);
        if (ibz_cmp(&tmp1, &d[l]) > 0) {
            /* t = round(Lambda[m-1][l-1] / d[l]) */
            ibz_add(&tmp1, &Lambda[m - 1][l - 1], &Lambda[m - 1][l - 1]);
            ibz_add(&tmp1, &tmp1, &d[l]);
            ibz_add(&tmp2, &d[l], &d[l]);
            ibz_div_floor(&t, &tmp1, &tmp1, &tmp2);

            if (!ibz_is_zero(&t)) {
                /* b[m-1] -= t * b[l-1] */
                ibz_vec_4_sub_scalar_mul(&b[m - 1], &t, &b[l - 1]);
                tracker_update_vec4(&b[m - 1]);

                /* Lambda[m-1][j] -= t * Lambda[l-1][j]  for j < l-1 */
                for (int j = 0; j < l - 1; j++) {
                    ibz_mul(&tmp1, &t, &Lambda[l - 1][j]);
                    track_intermediate(&tmp1);
                    ibz_sub(&Lambda[m - 1][j], &Lambda[m - 1][j], &tmp1);
                    track_intermediate(&Lambda[m - 1][j]);
                }

                /* Lambda[m-1][l-1] -= t * d[l] */
                ibz_mul(&tmp1, &t, &d[l]);
                track_intermediate(&tmp1);
                ibz_sub(&Lambda[m - 1][l - 1], &Lambda[m - 1][l - 1], &tmp1);
                track_intermediate(&Lambda[m - 1][l - 1]);
            }
        }
    }

    /* Zero vector: remove immediately */
    if (ibz_vec_4_is_zero(&b[m - 1]))
        goto remove_vector;

    /* Real-dependent but nonzero (d[m]=0, b[m-1]≠0):
     * b[m-1] is in the real span of b[0]..b[m-2] but NOT the integer span.
     * Swap it down one position. The next reduction iteration recomputes GSO
     * at the new position, where the vector may become independent or
     * continue swapping until it reaches position 1 and size-reduction
     * zeros it out (Euclidean-algorithm-like convergence). */
    if (ibz_is_zero(&d[m])) {
        if (m >= 2) {
            ibz_vec_4_swap(&b[m - 1], &b[m - 2]);
            if (m > 2)
                m = m - 1;
        }
        goto reduction;
    }

    if (l < m - 1)
        goto next_l;

    /* Lovász condition (integral form, Cohen 2.6.3):
     * For 0-based vectors (m-2, m-1):
     *   4*(d[m]*d[m-2] + Lambda[m-1][m-2]^2) < 3*d[m-1]^2
     */
    if (m >= 2 && !ibz_is_zero(&d[m - 1])) {
        ibz_t lhs, rhs, lam_sq;
        ibz_init(&lhs); ibz_init(&rhs); ibz_init(&lam_sq);

        /* LHS = 4 * (d[m] * d[m-2] + Lambda[m-1][m-2]^2) */
        ibz_mul(&lhs, &d[m], &d[m - 2]);
        ibz_mul(&lam_sq, &Lambda[m - 1][m - 2], &Lambda[m - 1][m - 2]);
        track_intermediate(&lam_sq);
        ibz_add(&lhs, &lhs, &lam_sq);
        ibz_add(&lhs, &lhs, &lhs);
        ibz_add(&lhs, &lhs, &lhs);  /* *4 */
        track_intermediate(&lhs);

        /* RHS = 3 * d[m-1]^2 */
        ibz_mul(&rhs, &d[m - 1], &d[m - 1]);
        ibz_mul(&tmp1, &rhs, &ibz_const_three);
        track_intermediate(&tmp1);

        int need_swap = (ibz_cmp(&lhs, &tmp1) < 0);

        ibz_finalize(&lhs); ibz_finalize(&rhs); ibz_finalize(&lam_sq);

        if (need_swap)
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

    /* ===== SWAP (NTL-style incremental update) ===== */
do_swap:
    {
        /* 0-based: swap b[k] and b[k-1], where k = m-1 */
        int k = m - 1;
        ibz_t lam, t1_s, t2_s, tmp_a, tmp_b;
        ibz_init(&lam); ibz_init(&t1_s); ibz_init(&t2_s);
        ibz_init(&tmp_a); ibz_init(&tmp_b);

        ibz_copy(&lam, &Lambda[k][k - 1]);

        /* Update Lambda[i][k-1] and Lambda[i][k] for i > k */
        for (int i = k + 1; i < beta; i++) {
            /* t1 = (Lambda[i][k-1]*lam + Lambda[i][k]*d[k-1]) / d[k] */
            ibz_mul(&tmp_a, &Lambda[i][k - 1], &lam);
            ibz_mul(&tmp_b, &Lambda[i][k], &d[k - 1]);
            ibz_add(&tmp_a, &tmp_a, &tmp_b);
            track_intermediate(&tmp_a);
            ibz_div(&t1_s, &rem, &tmp_a, &d[k]);
            assert(ibz_is_zero(&rem));

            /* t2 = (Lambda[i][k-1]*d[k+1] - Lambda[i][k]*lam) / d[k] */
            ibz_mul(&tmp_a, &Lambda[i][k - 1], &d[k + 1]);
            ibz_mul(&tmp_b, &Lambda[i][k], &lam);
            ibz_sub(&tmp_a, &tmp_a, &tmp_b);
            track_intermediate(&tmp_a);
            ibz_div(&t2_s, &rem, &tmp_a, &d[k]);
            assert(ibz_is_zero(&rem));

            ibz_copy(&Lambda[i][k - 1], &t1_s);
            ibz_copy(&Lambda[i][k], &t2_s);
            track_intermediate(&Lambda[i][k - 1]);
            track_intermediate(&Lambda[i][k]);
        }

        /* Swap Lambda[k][j] and Lambda[k-1][j] for j < k-1 */
        for (int j = 0; j < k - 1; j++)
            ibz_swap(&Lambda[k][j], &Lambda[k - 1][j]);

        /* d[k]_new = (d[k+1]*d[k-1] + lam^2) / d[k]_old */
        ibz_mul(&tmp_a, &d[k + 1], &d[k - 1]);
        ibz_mul(&tmp_b, &lam, &lam);
        ibz_add(&tmp_a, &tmp_a, &tmp_b);
        ibz_div(&d[k], &rem, &tmp_a, &d[k]);
        assert(ibz_is_zero(&rem));
        track_intermediate(&d[k]);

        /* Swap the vectors */
        ibz_vec_4_swap(&b[k], &b[k - 1]);
        tracker_update_vec4(&b[k]);
        tracker_update_vec4(&b[k - 1]);

        ibz_finalize(&lam); ibz_finalize(&t1_s); ibz_finalize(&t2_s);
        ibz_finalize(&tmp_a); ibz_finalize(&tmp_b);
    }

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
        ibz_set(&d[0], 1);
        for (int i = 0; i < beta; i++)
            compute_gso_row(i, b, Lambda, d, &alg->p);
        m = 2;
        if (m > beta)
            goto done;
        goto reduction;
    }

    ibz_set(&d[0], 1);
    for (int i = 0; i < beta; i++)
        compute_gso_row(i, b, Lambda, d, &alg->p);
    tau = m;
    if (tau < 2)
        tau = 2;
    goto load;

done:
    tracker_disable();

    /* Track bitsizes */
    mlll_call_count++;
    for (int i = 0; i < beta; i++) {
        for (int c = 0; c < 4; c++) {
            int bs = ibz_bitsize(&(b[i][c]));
            if (bs > max_bits_vec) max_bits_vec = bs;
        }
        for (int j = 0; j < i; j++) {
            int bs = ibz_bitsize(&Lambda[i][j]);
            if (bs > max_bits_lambda) max_bits_lambda = bs;
        }
    }
    for (int i = 0; i <= beta; i++) {
        int bs = ibz_bitsize(&d[i]);
        if (bs > max_bits_d) max_bits_d = bs;
    }

    /* Extract result */
    *rank = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&((*basis)[i][j]), 0);

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
        for (int j = 0; j < g; j++)
            ibz_finalize(&Lambda[i][j]);
    }
    for (int i = 0; i <= g; i++)
        ibz_finalize(&d[i]);
    ibz_finalize(&t);
    ibz_finalize(&tmp1);
    ibz_finalize(&tmp2);
    ibz_finalize(&rem);
}

/* ========== Lattice operations using MLLL ========== */

void
quat_lattice_mul_mlll(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg)
{
    /* Compact paper Algorithm 2: CompactIdealMultiplication
     * 1. Multiply out denominators so all vectors are integer
     * 2. LLL-reduce each basis
     * 3. Form 16 products (integer quaternion multiplication)
     * 4. MLLL on integer generators
     * 5. Result denominator = r1 * r2
     */
    ibz_vec_4_t elem1, elem2, elem_res;
    ibz_vec_4_t generators[16];
    int rank;

    ibz_vec_4_init(&elem1);
    ibz_vec_4_init(&elem2);
    ibz_vec_4_init(&elem_res);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_init(&generators[i]);

    /* Step 1-2: Multiply out denominators (r1=lat1->denom, r2=lat2->denom)
     * and LLL-reduce each basis.
     * quat_lattice_lll takes a lattice (basis+denom) and returns LLL-reduced basis.
     * The basis already includes the denom scaling implicitly,
     * so the reduced basis columns are integer coordinates with denom factored in.
     */
    ibz_mat_4x4_t red1, red2;
    ibz_mat_4x4_init(&red1);
    ibz_mat_4x4_init(&red2);

    quat_lattice_lll(&red1, lat1, alg);
    quat_lattice_lll(&red2, lat2, alg);

    /* Step 3: M = (alpha_i * beta_j) for 1<=i,j<=4
     * red1/red2 columns are the LLL-reduced bases (integer coords, denom separate).
     * We use the integer coordinates directly — denom is handled at the end.
     */
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

    /* Step 4: MLLL on integer generators */
    quat_mlll(&(res->basis), &rank, generators, 16, alg);
    assert(rank == 4);

    /* Step 5: Result denominator = r1 * r2 */
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
