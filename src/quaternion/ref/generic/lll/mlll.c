/**
 * @file mlll.c
 * @brief Modified LLL (MLLL) 알고리즘 구현 — Float GSO + Gram matrix 버전
 *
 * Pohst(1987)의 MLLL에 L² 스타일 float GSO를 적용.
 * KLKL25와 동일 접근: Gram matrix는 정수로 정확하게 유지,
 * GSO(Cholesky) 계수는 dpe_t로 근사 계산.
 *
 * - 벡터 좌표 b[i]: ibz_t (Lemma 1 bounded)
 * - Gram matrix G[i][j] = <b[i], b[j]>: ibz_t (exact)
 * - GSO 계수 r, u: dpe_t (53-bit float + extended exponent)
 * - Dependency detection: b[m] = 0 (exact integer check)
 *
 * Gram matrix는 size-reduce/swap 시 점진적 업데이트.
 * Cholesky는 Gram matrix에서 dpe_t로 빠르게 재계산.
 */

#include <quaternion.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"
#include "dpe.h"

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

/**
 * @brief Compute Gram row for position idx: G[idx][j] = <b[idx], b[j]> for j <= idx
 */
static void
gram_compute_row(int idx, int beta, ibz_t G[][N],
                 const ibz_vec_4_t *b, const ibz_t *p)
{
    for (int j = 0; j <= idx; j++)
        ibz_vec_4_dot_quat(&G[idx][j], &b[idx], &b[j], p);
}

/**
 * @brief Update Gram matrix after b[k] -= X * b[l]
 *
 * G[k][j] = <b[k]-X*b[l], b[j]> = G[k][j] - X*G[l][j]  for all j
 * G[k][k] needs special handling since both sides changed.
 */
static void
gram_size_reduce(int k, int l, const ibz_t *X, int beta,
                 ibz_t G[][N])
{
    ibz_t tmp;
    ibz_init(&tmp);

    /* G[k][k] = G[k][k] - 2*X*G[k][l] + X^2*G[l][l] */
    /* But easier: first update G[k][j] for all j != k, then G[k][k] */

    /* Step 1: G[k][k] -= X * G[k][l]  (partial, for the "left" b[k] change) */
    ibz_mul(&tmp, X, GRAM(G, k, l));
    ibz_sub(&G[k][k], &G[k][k], &tmp);

    /* Step 2: G[k][j] -= X * G[l][j]  for all j (including j=k after step 1) */
    for (int j = 0; j < beta; j++) {
        ibz_mul(&tmp, X, GRAM(G, l, j));
        ibz_sub(GRAM(G, k, j), GRAM(G, k, j), &tmp);
    }

    ibz_finalize(&tmp);
}

/**
 * @brief Swap rows/columns k and k-1 in the Gram matrix
 */
static void
gram_swap(int k, int beta, ibz_t G[][N])
{
    /* Swap G[k][·] and G[k-1][·] symmetrically */
    for (int j = 0; j < beta; j++) {
        if (j == k || j == k - 1)
            continue;
        ibz_swap(GRAM(G, k, j), GRAM(G, k - 1, j));
    }
    /* Swap diagonal */
    ibz_swap(&G[k][k], &G[k - 1][k - 1]);
    /* G[k][k-1] stays the same (it's symmetric) */
}

/* ---------- Cholesky from Gram (dpe_t) ---------- */

/**
 * @brief Compute Cholesky row idx from Gram matrix.
 *
 * r[idx][j] = G[idx][j] - sum_{k<j} r[idx][k]*u[j][k]
 * u[idx][j] = r[idx][j] / r[j][j]
 * r[idx][idx] = ||b*[idx]||²
 */
static void
cholesky_row(int idx,
             ibz_t G[][N],
             dpe_t r[][N],
             dpe_t u[][N])
{
    dpe_t tmpF;
    dpe_init(tmpF);

    for (int j = 0; j <= idx; j++) {
        dpe_set_z(r[idx][j], *GRAM(G, idx, j));
        for (int k = 0; k < j; k++) {
            dpe_mul(tmpF, r[idx][k], u[j][k]);
            dpe_sub(r[idx][j], r[idx][j], tmpF);
        }
        if (j < idx) {
            if (dpe_zero_p(r[j][j]))
                dpe_set_d(u[idx][j], 0.0);
            else
                dpe_div(u[idx][j], r[idx][j], r[j][j]);
        }
    }

    dpe_clear(tmpF);
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
    ibz_t G[N][N]; /* Gram matrix (lower triangular) */
    dpe_t r[N][N]; /* Cholesky r */
    dpe_t u[N][N]; /* Cholesky u = mu */

    for (int i = 0; i < g; i++) {
        ibz_vec_4_init(&b[i]);
        for (int j = 0; j < g; j++) {
            ibz_init(&G[i][j]);
            dpe_init(r[i][j]);
            dpe_init(u[i][j]);
        }
    }

    dpe_t Xf, tmpF, delta_bar;
    ibz_t X, tmpI;
    dpe_init(Xf);
    dpe_init(tmpF);
    dpe_init(delta_bar);
    ibz_init(&X);
    ibz_init(&tmpI);

    dpe_set_d(delta_bar, 0.75);

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

        /* Compute Gram row and Cholesky */
        gram_compute_row(beta, beta + 1, G, b, &alg->p);
        cholesky_row(beta, G, r, u);
        beta++;

        if (!dpe_zero_p(r[beta - 1][beta - 1]) && alpha < g)
            continue;

        break;
    }

    if (beta <= 1)
        goto done;

    m = tau;
    if (m >= beta)
        m = beta - 1;  /* ensure last loaded vector gets reduced */
    if (m < 1)
        m = 1;

    /* ===== REDUCTION ===== */
reduction:
    l = m - 1;

size_reduce:
    if (l >= 0 && !dpe_zero_p(r[l][l])) {
        /* Recompute Cholesky row m from Gram */
        cholesky_row(m, G, r, u);

        if (dpe_cmp_d(u[m][l], 0.5) >= 0 || dpe_cmp_d(u[m][l], -0.5) <= 0) {
            dpe_set(Xf, u[m][l]);
            dpe_round(Xf, Xf);
            dpe_get_z(X, Xf);

            if (!ibz_is_zero(&X)) {
                ibz_vec_4_sub_scalar_mul(&b[m], &X, &b[l]);
                tracker_update_vec4(&b[m]);
                gram_size_reduce(m, l, &X, beta, G);

                /* Update u for j < l */
                for (int j = 0; j < l; j++) {
                    dpe_mul(tmpF, Xf, u[l][j]);
                    dpe_sub(u[m][j], u[m][j], tmpF);
                }
                dpe_sub(u[m][l], u[m][l], Xf);
            }
        }
    }

    if (ibz_vec_4_is_zero(&b[m]))
        goto remove_vector;

    /* Lovász: r[m][m] < (delta - u[m][m-1]^2) * r[m-1][m-1] */
    cholesky_row(m, G, r, u);
    if (m > 0 && !dpe_zero_p(r[m - 1][m - 1])) {
        dpe_t lhs, rhs, mu_sq;
        dpe_init(lhs);
        dpe_init(rhs);
        dpe_init(mu_sq);

        dpe_set(lhs, r[m][m]);
        dpe_mul(mu_sq, u[m][m - 1], u[m][m - 1]);
        dpe_sub(rhs, delta_bar, mu_sq);
        dpe_mul(rhs, rhs, r[m - 1][m - 1]);

        int need_swap = (dpe_cmp(lhs, rhs) < 0);

        dpe_clear(lhs);
        dpe_clear(rhs);
        dpe_clear(mu_sq);

        if (need_swap)
            goto do_swap;
    }

    l--;
    if (l >= 0)
        goto size_reduce;
    m++;
    if (m >= beta) {
        if (alpha < g) {
            tau = beta;
            goto load;
        }
        goto done;
    }
    goto reduction;

    /* ===== SWAP ===== */
do_swap:
    ibz_vec_4_swap(&b[m], &b[m - 1]);
    tracker_update_vec4(&b[m]);
    tracker_update_vec4(&b[m - 1]);
    gram_swap(m, beta, G);

    /* Recompute Cholesky from m-1 */
    for (int i = m - 1; i < beta; i++)
        cholesky_row(i, G, r, u);

    if (m > tau)
        m--;
    goto reduction;

    /* ===== REMOVE ===== */
remove_vector:
    for (int i = m; i < beta - 1; i++) {
        for (int c = 0; c < 4; c++)
            ibz_copy(&(b[i][c]), &(b[i + 1][c]));
    }
    beta--;

    if (alpha >= g)
        goto done;

    /* Recompute Gram and Cholesky from m */
    for (int i = m; i < beta; i++) {
        gram_compute_row(i, beta, G, b, &alg->p);
        cholesky_row(i, G, r, u);
    }

    tau = m + 1;
    if (tau < 1)
        tau = 1;
    goto load;

done:
    assert(alpha >= g);

    /* Use HNF to extract the correct rank-4 basis from all b[] vectors,
     * then LLL-reduce it with L². */
    {
        /* Compute modulus: determinant of any non-singular 4x4 submatrix */
        ibz_t mod;
        ibz_init(&mod);

        /* Try first 4 non-zero vectors */
        ibz_mat_4x4_t tmpmat;
        ibz_mat_4x4_init(&tmpmat);
        int cnt = 0;
        for (int i = 0; i < beta && cnt < 4; i++) {
            if (!ibz_vec_4_is_zero(&b[i])) {
                for (int j = 0; j < 4; j++)
                    ibz_copy(&tmpmat[j][cnt], &(b[i][j]));
                cnt++;
            }
        }
        if (cnt == 4)
            ibz_mat_4x4_inv_with_det_as_denom(NULL, &mod, &tmpmat);
        ibz_abs(&mod, &mod);

        /* If det=0, try all subsets until we find a non-singular one */
        if (ibz_is_zero(&mod) && beta >= 4) {
            for (int a = 0; a < beta - 3 && ibz_is_zero(&mod); a++)
                for (int c = a + 1; c < beta - 2 && ibz_is_zero(&mod); c++)
                    for (int d = c + 1; d < beta - 1 && ibz_is_zero(&mod); d++)
                        for (int e = d + 1; e < beta && ibz_is_zero(&mod); e++) {
                            int idx[4] = {a, c, d, e};
                            for (int ii = 0; ii < 4; ii++)
                                for (int jj = 0; jj < 4; jj++)
                                    ibz_copy(&tmpmat[jj][ii], &(b[idx[ii]][jj]));
                            ibz_mat_4x4_inv_with_det_as_denom(NULL, &mod, &tmpmat);
                            ibz_abs(&mod, &mod);
                        }
        }
        ibz_mat_4x4_finalize(&tmpmat);

        if (!ibz_is_zero(&mod) && beta > 4) {
            /* HNF to get exact rank-4 basis */
            ibz_mat_4xn_hnf_mod_core(basis, beta, (const ibz_vec_4_t *)b, &mod);
            *rank = 4;
        } else if (!ibz_is_zero(&mod)) {
            /* beta <= 4: just copy */
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
        } else {
            /* Fallback: just copy first non-zero vectors */
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
        }
        ibz_finalize(&mod);
    }

    /* LLL-reduce the rank-4 basis using the proven L² algorithm */
    if (*rank == 4) {
        ibz_mat_4x4_t Gfinal;
        ibz_mat_4x4_init(&Gfinal);
        ibz_t tmp2;
        ibz_init(&tmp2);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j <= i; j++) {
                ibz_set(&Gfinal[i][j], 0);
                for (int k = 0; k < 4; k++) {
                    ibz_mul(&tmp2, &(*basis)[k][i], &(*basis)[k][j]);
                    if (k >= 2)
                        ibz_mul(&tmp2, &tmp2, &alg->p);
                    ibz_add(&Gfinal[i][j], &Gfinal[i][j], &tmp2);
                }
                ibz_mul(&Gfinal[i][j], &Gfinal[i][j], &ibz_const_two);
            }
            for (int j = i + 1; j < 4; j++)
                ibz_copy(&Gfinal[i][j], &Gfinal[j][i]);
        }
        ibz_finalize(&tmp2);
        quat_lll_core(&Gfinal, basis);
        ibz_mat_4x4_finalize(&Gfinal);
    }

    /* Cleanup */
    for (int i = 0; i < g; i++) {
        ibz_vec_4_finalize(&b[i]);
        for (int j = 0; j < g; j++) {
            ibz_finalize(&G[i][j]);
            dpe_clear(r[i][j]);
            dpe_clear(u[i][j]);
        }
    }
    dpe_clear(Xf);
    dpe_clear(tmpF);
    dpe_clear(delta_bar);
    ibz_finalize(&X);
    ibz_finalize(&tmpI);
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
