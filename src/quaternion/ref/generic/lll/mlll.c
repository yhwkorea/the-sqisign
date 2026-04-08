/**
 * @file mlll.c
 * @brief Modified LLL (MLLL) algorithm — exact integer d/λ version
 *
 * Pohst (1987) MLLL with exact integer d[] and λ[][] (Cohen 2.6.7 style).
 * Gram matrix is exact integer, incrementally updated during size-reduce/swap.
 * After each swap, d/lam are recomputed from scratch for robustness.
 *
 * - Vector coordinates b[i]: ibz_t (Lemma 1 bounded)
 * - Gram matrix G[i][j] = <b[i], b[j]>: ibz_t (exact)
 * - d[k]: Gram determinant of b[0..k-1], ibz_t (exact)
 * - lam[i][j] = d[j+1]*μ[i][j]: ibz_t (exact)
 * - Dependency detection: d[k+1] == 0 (exact integer test)
 *
 * No HNF post-processing needed — dependent vectors are eliminated through
 * exact size-reduction and swaps.
 */

#include <quaternion.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"

#define N MLLL_MAX_GENERATORS

/* Access symmetric Gram entry: G[max(i,j)][min(i,j)] */
#define GRAM(G, i, j) ((i) < (j) ? &(G)[(j)][(i)] : &(G)[(i)][(j)])

/* ---------- helpers ---------- */

static void
ibz_vec_4_dot_quat(ibz_t *dot, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const ibz_t *p)
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

/* ---------- Gram matrix operations ---------- */

static void
gram_compute_row(int idx, int beta, ibz_t G[][N],
                 const ibz_vec_4_t *b, const ibz_t *p)
{
    (void)beta;
    for (int j = 0; j <= idx; j++)
        ibz_vec_4_dot_quat(&G[idx][j], &b[idx], &b[j], p);
}

/**
 * @brief Update Gram matrix after b[k] -= X * b[l]
 *
 * Two-step update: G[k][k] -= X*G[k][l], then G[k][j] -= X*G[l][j] for all j.
 * Since l < k, the G[k][l] entry is updated before G[k][k] in the loop,
 * yielding the correct result: G[k][k] = old - 2X*G[k][l] + X²*G[l][l].
 */
static void
gram_size_reduce(int k, int l, const ibz_t *X, int beta,
                 ibz_t G[][N])
{
    ibz_t tmp;
    ibz_init(&tmp);

    ibz_mul(&tmp, X, GRAM(G, k, l));
    ibz_sub(&G[k][k], &G[k][k], &tmp);

    for (int j = 0; j < beta; j++) {
        ibz_mul(&tmp, X, GRAM(G, l, j));
        ibz_sub(GRAM(G, k, j), GRAM(G, k, j), &tmp);
    }

    ibz_finalize(&tmp);
}

static void
gram_swap(int k, int beta, ibz_t G[][N])
{
    for (int j = 0; j < beta; j++) {
        if (j == k || j == k - 1)
            continue;
        ibz_swap(GRAM(G, k, j), GRAM(G, k - 1, j));
    }
    ibz_swap(&G[k][k], &G[k - 1][k - 1]);
}

/* ---------- exact integer helpers ---------- */

/**
 * @brief Nearest-integer division: q = round(a/b)
 */
static void
rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t r, sign_q, abs_b;
    ibz_init(&r);
    ibz_init(&sign_q);
    ibz_init(&abs_b);

    ibz_abs(&abs_b, b);
    ibz_mul(&sign_q, a, b);
    ibz_div(q, &r, a, b);
    ibz_abs(&r, &r);
    ibz_add(&r, &r, &r);
    ibz_set(&sign_q,
            (1 - 2 * (ibz_cmp(&sign_q, &ibz_const_zero) < 0)) *
            (ibz_cmp(&r, &abs_b) > 0));
    ibz_add(q, q, &sign_q);

    ibz_finalize(&r);
    ibz_finalize(&sign_q);
    ibz_finalize(&abs_b);
}

/**
 * @brief Exact integer division: q = a/b (asserts remainder == 0)
 */
static void
exact_div(ibz_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t r;
    ibz_init(&r);
    ibz_div(q, &r, a, b);
    assert(ibz_is_zero(&r));
    ibz_finalize(&r);
}

/**
 * @brief Compute lam[k][0..k-1] and d[k+1] from Gram matrix.
 *
 * Cohen Algorithm 2.6.3 recurrence. When d[s] = 0 (R-dependent vectors),
 * the numerator is provably also 0 — we handle this by setting num = 0.
 */
static void
compute_lam_d_row(int k, ibz_t G[][N], ibz_t lam[][N], ibz_t *d)
{
    ibz_t num, t1, t2;
    ibz_init(&num);
    ibz_init(&t1);
    ibz_init(&t2);

    for (int j = 0; j <= k; j++) {
        ibz_copy(&num, GRAM(G, k, j));

        for (int s = 0; s < j; s++) {
            ibz_mul(&t1, &d[s + 1], &num);
            if (j < k)
                ibz_mul(&t2, &lam[k][s], &lam[j][s]);
            else
                ibz_mul(&t2, &lam[k][s], &lam[k][s]);
            ibz_sub(&t1, &t1, &t2);

            if (ibz_is_zero(&d[s])) {
                /* d[s]=0: vectors 0..s-1 are R-dependent.
                 * By algebraic structure, numerator is also 0. */
                ibz_set(&num, 0);
            } else {
                exact_div(&num, &t1, &d[s]);
            }
        }

        if (j < k)
            ibz_copy(&lam[k][j], &num);
        else
            ibz_copy(&d[k + 1], &num);
    }

    ibz_finalize(&num);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
}

/**
 * @brief Recompute all d[] and lam[][] from scratch using current Gram matrix.
 */
static void
recompute_all_lam_d(int beta, ibz_t G[][N], ibz_t lam[][N], ibz_t *d)
{
    ibz_set(&d[0], 1);
    for (int k = 0; k < beta; k++)
        compute_lam_d_row(k, G, lam, d);
}

/* ========== MLLL core ========== */

void
quat_mlll(ibz_mat_4x4_t *basis,
          int *rank,
          const ibz_vec_4_t *generators,
          int g,
          const quat_alg_t *alg)
{
    assert(g >= 1 && g <= N);

    ibz_vec_4_t b[N];
    ibz_t G[N][N];      /* Gram matrix (lower triangular, exact) */
    ibz_t d[N + 1];     /* d[0]=1, d[k+1] = det(Gram(b[0..k])) */
    ibz_t lam[N][N];    /* lam[i][j] = d[j+1]*μ[i][j], j < i */

    for (int i = 0; i < g; i++) {
        ibz_vec_4_init(&b[i]);
        for (int j = 0; j < g; j++) {
            ibz_init(&G[i][j]);
            ibz_init(&lam[i][j]);
        }
    }
    for (int i = 0; i <= g; i++)
        ibz_init(&d[i]);
    ibz_set(&d[0], 1);

    ibz_t X, tmp, tmp2;
    ibz_init(&X);
    ibz_init(&tmp);
    ibz_init(&tmp2);

    int alpha = 0, beta = 0, tau = 1, m, l;

    /* ===== LOAD ===== */
load:
    while (alpha < g) {
        if (ibz_vec_4_is_zero(&generators[alpha])) {
            alpha++;
            continue;
        }

        for (int i = 0; i < 4; i++)
            ibz_copy(&(b[beta][i]), &(generators[alpha][i]));
        alpha++;

        /* Compute Gram row for new vector */
        gram_compute_row(beta, beta + 1, G, b, &alg->p);
        beta++;

        /* Recompute ALL d/lam from scratch (handles stale values from
         * prior flagg=1 swaps and natural R-dependencies) */
        recompute_all_lam_d(beta, G, lam, d);

        if (!ibz_is_zero(&d[beta]) && alpha < g)
            continue;   /* independent vector, load more */

        break;  /* dependent or last generator → reduce */
    }

    if (beta <= 1)
        goto done;

    m = tau;
    if (m >= beta)
        m = beta - 1;
    if (m < 1)
        m = 1;

    /* ===== REDUCTION ===== */
reduction:
    l = m - 1;

size_reduce:
    /* REDI(m, l): exact size-reduce b[m] against b[l] */
    if (l >= 0 && !ibz_is_zero(&d[l + 1])) {
        /* Check |2*lam[m][l]| > d[l+1] */
        ibz_add(&tmp, &lam[m][l], &lam[m][l]);
        ibz_abs(&tmp, &tmp);
        if (ibz_cmp(&tmp, &d[l + 1]) > 0) {
            rounded_div(&X, &lam[m][l], &d[l + 1]);

            if (!ibz_is_zero(&X)) {
                /* b[m] -= X * b[l] */
                ibz_vec_4_sub_scalar_mul(&b[m], &X, &b[l]);
                tracker_update_vec4(&b[m]);

                /* Update Gram matrix */
                gram_size_reduce(m, l, &X, beta, G);

                /* Update lam: lam[m][l] -= X*d[l+1], lam[m][j] -= X*lam[l][j] */
                ibz_mul(&tmp, &X, &d[l + 1]);
                ibz_sub(&lam[m][l], &lam[m][l], &tmp);
                for (int j = 0; j < l; j++) {
                    ibz_mul(&tmp, &X, &lam[l][j]);
                    ibz_sub(&lam[m][j], &lam[m][j], &tmp);
                }
            }
        }
    }

    if (ibz_vec_4_is_zero(&b[m]))
        goto remove_vector;

    l--;
    if (l >= 0)
        goto size_reduce;

    /* === All size-reductions for b[m] done. Check swap conditions. === */

    /* Case 1: d[m+1]=0 — dependent vector (Matthews D_k=0).
     * Force swap to cascade dependent vector toward position 1,
     * where REDI against a shorter vector can eliminate it.
     * Must bypass tau barrier (unlike Lovász swaps). */
    if (ibz_is_zero(&d[m + 1]))
        goto do_swap_dep;

    /* Case 2: Lovász condition.
     * Swap if 4*(d[m-1]*d[m+1] + lam[m][m-1]²) < 3*d[m]² */
    ibz_mul(&tmp, &d[m - 1], &d[m + 1]);
    ibz_mul(&tmp2, &lam[m][m - 1], &lam[m][m - 1]);
    ibz_add(&tmp, &tmp, &tmp2);
    ibz_add(&tmp, &tmp, &tmp);    /* ×2 */
    ibz_add(&tmp, &tmp, &tmp);    /* ×4 */

    ibz_mul(&tmp2, &d[m], &d[m]);
    ibz_mul(&tmp2, &tmp2, &ibz_const_three);

    if (ibz_cmp(&tmp, &tmp2) < 0)
        goto do_swap_lov;

    m++;
    if (m >= beta) {
        if (alpha < g) {
            tau = beta;
            goto load;
        }
        goto done;
    }
    goto reduction;

    /* ===== SWAP (dependent vector) ===== */
do_swap_dep:
    ibz_vec_4_swap(&b[m], &b[m - 1]);
    tracker_update_vec4(&b[m]);
    tracker_update_vec4(&b[m - 1]);
    gram_swap(m, beta, G);
    recompute_all_lam_d(beta, G, lam, d);
    if (m > 1)
        m--;
    goto reduction;

    /* ===== SWAP (Lovász) ===== */
do_swap_lov:
    ibz_vec_4_swap(&b[m], &b[m - 1]);
    tracker_update_vec4(&b[m]);
    tracker_update_vec4(&b[m - 1]);
    gram_swap(m, beta, G);
    recompute_all_lam_d(beta, G, lam, d);
    if (m > tau)
        m--;
    goto reduction;

    /* ===== REMOVE ===== */
remove_vector:
    /* Shift b[] down from position m */
    for (int i = m; i < beta - 1; i++) {
        for (int c = 0; c < 4; c++)
            ibz_copy(&(b[i][c]), &(b[i + 1][c]));
    }
    beta--;

    /* Recompute Gram for shifted positions */
    for (int i = m; i < beta; i++)
        gram_compute_row(i, beta, G, b, &alg->p);

    /* Recompute all d/lam from scratch */
    recompute_all_lam_d(beta, G, lam, d);

    if (alpha >= g) {
        /* All generators loaded. Re-run reduction from m=1 to
         * ensure full LLL-reducedness after forced-swap cascades. */
        if (beta <= 1)
            goto done;
        m = 1;
        goto reduction;
    }

    tau = m + 1;
    if (tau < 1)
        tau = 1;
    goto load;

done:
    /* Stop tracking: reduction loop is done */
    tracker_disable();

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
        for (int j = 0; j < g; j++) {
            ibz_finalize(&G[i][j]);
            ibz_finalize(&lam[i][j]);
        }
    }
    for (int i = 0; i <= g; i++)
        ibz_finalize(&d[i]);
    ibz_finalize(&X);
    ibz_finalize(&tmp);
    ibz_finalize(&tmp2);
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

    for (int k = 0; k < 4; k++) {
        ibz_vec_4_copy_ibz(
            &elem1, &(lat1->basis[0][k]), &(lat1->basis[1][k]),
            &(lat1->basis[2][k]), &(lat1->basis[3][k]));
        for (int i = 0; i < 4; i++) {
            ibz_vec_4_copy_ibz(
                &elem2, &(lat2->basis[0][i]), &(lat2->basis[1][i]),
                &(lat2->basis[2][i]), &(lat2->basis[3][i]));
            quat_alg_coord_mul(&elem_res, &elem1, &elem2, alg);
            for (int j = 0; j < 4; j++)
                ibz_copy(&(generators[4 * k + i][j]), &(elem_res[j]));
        }
    }

    quat_mlll(&(res->basis), &rank, generators, 16, alg);
    assert(rank == 4);

    ibz_mul(&(res->denom), &(lat1->denom), &(lat2->denom));
    quat_lattice_reduce_denom(res, res);

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
