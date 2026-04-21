/**
 * @file mlll_gram.c
 * @brief Gram-based Modified LLL — intermediate integers bounded by Lemma 3.
 *
 * Combines:
 *   - L2 algorithm (l2.c): integer Gram matrix G[i][j] = <b_i, b_j>,
 *     floating-point (dpe) Cholesky/GSO r, u, lovasz.
 *   - Pohst (1987): generator loading loop + dependent-vector removal.
 *
 * Key invariant:
 *   All integer objects manipulated (basis coords b[i], Gram entries G[i][j])
 *   have bitsize at most max_{k} <a_k, a_k>. This matches Lemma 3 of the
 *   "Compact Quaternion Algorithms for SQIsign" paper.
 */

#include <quaternion.h>
#include <stdio.h>
#include <assert.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"
#include "dpe.h"

#define N MLLL_MAX_GENERATORS

/* Symmetric Gram access: store lower triangular (i >= j). */
#define G_SYM(G, i, j) ((i) >= (j) ? &(G)[(i)][(j)] : &(G)[(j)][(i)])

/* ---------- Phase 2 candidate B: mpz_realloc2 prealloc scaffold ---------- */

/* Process-global mode flag. 0 = baseline, 1 = prealloc hints on.
 * Set via quat_mlll_gram_set_prealloc_mode() by benchmarks/tests. */
static int g_prealloc_mode = 0;

void
quat_mlll_gram_set_prealloc_mode(int mode)
{
    g_prealloc_mode = (mode != 0) ? 1 : 0;
}

int
quat_mlll_gram_get_prealloc_mode(void)
{
    return g_prealloc_mode;
}

/* Level hints: map p bitsize to (vec, Gram) widths measured in Phase 1
 * (10k trials, C1-fixed). Widths rounded up to u64 limb boundary with
 * +64 bit margin beyond the max-observed value — tracked in
 * project_sqisign_mlll.md B*_{1,3,5} table.
 *
 * Returns 1 if level matched, 0 otherwise (no prealloc done in that case).
 */
static int
mlll_gram_level_hints(const quat_alg_t *alg, int *vec_bits, int *gram_bits)
{
    int p_bits = ibz_bitsize(&alg->p);
    if (p_bits <= 128) {             /* L1: p ~127 bits, vec 259, Gram 518 */
        *vec_bits  = 320;            /* 5 u64 */
        *gram_bits = 576;            /* 9 u64 */
        return 1;
    } else if (p_bits <= 200) {      /* L3: p ~193 bits, vec 391, Gram 782 */
        *vec_bits  = 448;            /* 7 u64 */
        *gram_bits = 832;            /* 13 u64 */
        return 1;
    } else if (p_bits <= 256) {      /* L5: p ~254 bits, vec 513, Gram 1026 */
        *vec_bits  = 576;            /* 9 u64 */
        *gram_bits = 1088;           /* 17 u64 */
        return 1;
    }
    *vec_bits = 0;
    *gram_bits = 0;
    return 0;
}

/* Prealloc hint on all buffers used by quat_mlll_gram. Called right after
 * the ibz_init/dpe_init loop when g_prealloc_mode == 1. The b[] and G[][]
 * matrices hold Lemma 1 / Lemma 3 bounded values; tmp holds transient
 * products X*b[i][j] or X*G[i][j] of up to (vec+gram) bits. */
static void
mlll_gram_prealloc_buffers(ibz_vec_4_t b[N], ibz_t G[N][N],
                           ibz_t *X, ibz_t *tmp,
                           int vec_bits, int gram_bits)
{
    mp_bitcnt_t vw = (mp_bitcnt_t)vec_bits;
    mp_bitcnt_t gw = (mp_bitcnt_t)gram_bits;
    mp_bitcnt_t tw = (mp_bitcnt_t)(vec_bits + gram_bits + 64);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < 4; j++)
            mpz_realloc2(b[i][j], vw);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            mpz_realloc2(G[i][j], gw);
    mpz_realloc2(*X, vw);
    mpz_realloc2(*tmp, tw);
}

/* ---------- helpers ---------- */

/* Quaternion norm bilinear form: <a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3 */
static void
vec4_dot_p(ibz_t *dot, const ibz_vec_4_t *a, const ibz_vec_4_t *b,
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

/* Fill G[row][j] for j=0..row (inclusive) from b[]. */
static void
gram_fill_row(ibz_t G[N][N], const ibz_vec_4_t b[N],
              int row, const ibz_t *p)
{
    for (int j = 0; j <= row; j++) {
        vec4_dot_p(&G[row][j], &b[row], &b[j], p);
        tracker_update_gso_ibz(&G[row][j]);
    }
}

/* Swap row/col a and b in the symmetric Gram (lower-triangular storage). */
static void
gram_swap(ibz_t G[N][N], int a, int b, int beta)
{
    if (a == b) return;
    if (a > b) { int t = a; a = b; b = t; }
    /* Now a < b. */
    ibz_swap(&G[a][a], &G[b][b]);  /* diagonals */
    /* G[b][a] stays (it's <b_a_new, b_b_new> = <old_b, old_a>, same value). */
    for (int k = 0; k < beta; k++) {
        if (k == a || k == b) continue;
        ibz_swap(G_SYM(G, a, k), G_SYM(G, b, k));
    }
}

/* ---------- main ---------- */

void
quat_mlll_gram(ibz_mat_4x4_t *basis,
               int *rank,
               const ibz_vec_4_t *generators,
               int g,
               const quat_alg_t *alg)
{
    assert(g >= 1 && g <= N);

    ibz_vec_4_t b[N];
    ibz_t G[N][N];            /* symmetric, lower-triangular storage */
    dpe_t r[N][N], u[N][N];   /* Cholesky r, GSO coefficients u (mu) */
    dpe_t lovasz[N];
    dpe_t delta_bar, Xf, tmpF;
    ibz_t X, tmp;

    for (int i = 0; i < N; i++) {
        ibz_vec_4_init(&b[i]);
        for (int j = 0; j < N; j++) {
            ibz_init(&G[i][j]);
        }
        for (int j = 0; j <= i; j++) {
            dpe_init(r[i][j]);
            dpe_init(u[i][j]);
        }
        dpe_init(lovasz[i]);
    }
    ibz_init(&X);
    ibz_init(&tmp);
    dpe_init(delta_bar); dpe_set_d(delta_bar, DELTABAR);
    dpe_init(Xf);
    dpe_init(tmpF);

    /* Phase 2 candidate B scaffold: prealloc hint for ibz_t buffers. */
    if (g_prealloc_mode) {
        int vec_hint_bits = 0, gram_hint_bits = 0;
        if (mlll_gram_level_hints(alg, &vec_hint_bits, &gram_hint_bits))
            mlll_gram_prealloc_buffers(b, G, &X, &tmp,
                                       vec_hint_bits, gram_hint_bits);
    }

    int alpha = 0, beta = 0, kappa = 0;

    /* Load first non-zero generator. */
    while (alpha < g && ibz_vec_4_is_zero(&generators[alpha]))
        alpha++;
    if (alpha >= g) {
        *rank = 0;
        goto cleanup;
    }
    for (int j = 0; j < 4; j++)
        ibz_copy(&b[0][j], &generators[alpha][j]);
    alpha++;
    gram_fill_row(G, b, 0, &alg->p);
    tracker_update_vec4(&b[0]);
    dpe_set_z(r[0][0], G[0][0]);
    beta = 1;
    kappa = 1;

    /* Iteration caps: guard against dpe precision loss causing infinite loops. */
    int outer_iter = 0;
    while (1) {
        if (++outer_iter >= 100000) {
            fprintf(stderr, "mlll_gram: outer iteration cap exceeded\n");
            abort();
        }
        /* Load next generator if we've consumed the current working set. */
        if (kappa >= beta) {
            while (alpha < g && ibz_vec_4_is_zero(&generators[alpha]))
                alpha++;
            if (alpha >= g) break;
            for (int j = 0; j < 4; j++)
                ibz_copy(&b[beta][j], &generators[alpha][j]);
            alpha++;
            gram_fill_row(G, b, beta, &alg->p);
            tracker_update_vec4(&b[beta]);
            beta++;
            /* kappa == beta-1 now is the position to process. */
        }

        /* ----- Size-reduce b[kappa] (L2 style with Gram updates). ----- */
        int done = 0;
        int size_iter = 0;
        while (!done) {
            if (++size_iter >= 64) {
                fprintf(stderr, "mlll_gram: size-reduce iteration cap exceeded\n");
                abort();
            }
            /* Recompute row kappa of Cholesky from Gram. */
            for (int j = 0; j <= kappa; j++) {
                dpe_set_z(r[kappa][j], G[kappa][j]);
                for (int k = 0; k < j; k++) {
                    dpe_mul(tmpF, r[kappa][k], u[j][k]);
                    dpe_sub(r[kappa][j], r[kappa][j], tmpF);
                }
                if (j < kappa)
                    dpe_div(u[kappa][j], r[kappa][j], r[j][j]);
            }

            done = 1;
            for (int i = kappa - 1; i >= 0; i--) {
                if (dpe_cmp_d(u[kappa][i], ETABAR) > 0 ||
                    dpe_cmp_d(u[kappa][i], -ETABAR) < 0) {
                    done = 0;
                    dpe_set(Xf, u[kappa][i]);
                    dpe_round(Xf, Xf);
                    dpe_get_z(X, Xf);

                    /* b[kappa] -= X * b[i] */
                    for (int j = 0; j < 4; j++) {
                        ibz_mul(&tmp, &X, &b[i][j]);
                        ibz_sub(&b[kappa][j], &b[kappa][j], &tmp);
                    }
                    tracker_update_vec4(&b[kappa]);

                    /* Gram update (same identity as l2.c):
                     *   first  G[kappa][kappa] -= X * G[kappa][i]
                     *   then   G[kappa][j]     -= X * G[i][j]  for all j in [0, beta).
                     * The second loop at j=i absorbs the remaining X^2*G[i][i] term,
                     * and at j=kappa adds another -X*G[i][kappa] matching the first.
                     */
                    ibz_mul(&tmp, &X, G_SYM(G, kappa, i));
                    ibz_sub(&G[kappa][kappa], &G[kappa][kappa], &tmp);
                    for (int j = 0; j < beta; j++) {
                        ibz_t *gkj = G_SYM(G, kappa, j);
                        ibz_mul(&tmp, &X, G_SYM(G, i, j));
                        ibz_sub(gkj, gkj, &tmp);
                        tracker_update_gso_ibz(gkj);
                    }

                    /* u[kappa][j] -= X * u[i][j] for j < i */
                    for (int j = 0; j < i; j++) {
                        dpe_mul(tmpF, Xf, u[i][j]);
                        dpe_sub(u[kappa][j], u[kappa][j], tmpF);
                    }
                }
            }
        }

        /* ----- Dependent-vector detection. ----- */
        /* Quaternion norm is positive definite, so <b_k, b_k> == 0 iff b_k == 0.
         * After size-reduce, dependent vectors collapse to zero. */
        if (ibz_is_zero(&G[kappa][kappa])) {
            /* REMOVE b[kappa]: shift b[kappa+1..beta-1] down, likewise Gram rows.
             * Simplest: do repeated adjacent swaps then drop last. */
            for (int i = kappa; i + 1 < beta; i++) {
                for (int c = 0; c < 4; c++)
                    ibz_swap(&b[i][c], &b[i + 1][c]);
                gram_swap(G, i, i + 1, beta);
            }
            beta--;
            /* Stay at kappa; it now holds what was b[kappa+1]. */
            if (kappa >= beta) {
                /* Dropped vector was last; try to load more. */
                continue;
            }
            /* Otherwise, the row at kappa is unchanged in position but its
             * Gram entries against lower rows might have shifted. Since we
             * used adjacent swaps that preserve Gram correctness, nothing to
             * recompute. But Cholesky row at kappa needs redo next iter. */
            continue;
        }

        /* ----- Lovasz check (same as l2.c). ----- */
        dpe_set_z(lovasz[0], G[kappa][kappa]);
        for (int i = 1; i < kappa; i++) {
            dpe_mul(tmpF, u[kappa][i - 1], r[kappa][i - 1]);
            dpe_sub(lovasz[i], lovasz[i - 1], tmpF);
        }
        int swap;
        for (swap = kappa; swap > 0; swap--) {
            dpe_mul(tmpF, delta_bar, r[swap - 1][swap - 1]);
            if (dpe_cmp(tmpF, lovasz[swap - 1]) < 0)
                break;
        }

        if (swap != kappa) {
            /* Move b[kappa] to position `swap` via successive adjacent swaps
             * (j = kappa, kappa-1, ..., swap+1). This matches l2.c semantics
             * even though our Gram is NxN (not 4x4). */
            for (int j = kappa; j > swap; j--) {
                for (int c = 0; c < 4; c++)
                    ibz_swap(&b[j][c], &b[j - 1][c]);
                gram_swap(G, j - 1, j, beta);
            }
            /* Copy r, u row to new position `swap`. Entries at indices >= swap
             * will be recomputed when those kappa values are revisited. */
            for (int i = 0; i < swap; i++) {
                dpe_set(u[swap][i], u[kappa][i]);
                dpe_set(r[swap][i], r[kappa][i]);
            }
            dpe_set(r[swap][swap], lovasz[swap]);
            kappa = swap;
        }

        kappa++;
    }

    *rank = beta;

    /* Extract result: columns of basis = b[0..beta-1] (truncate or zero). */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&((*basis)[i][j]), 0);
    int out = 0;
    for (int i = 0; i < beta && out < 4; i++) {
        if (!ibz_vec_4_is_zero(&b[i])) {
            for (int j = 0; j < 4; j++)
                ibz_copy(&((*basis)[j][out]), &b[i][j]);
            out++;
        }
    }
    *rank = out;

cleanup:
    for (int i = 0; i < N; i++) {
        ibz_vec_4_finalize(&b[i]);
        for (int j = 0; j < N; j++) {
            ibz_finalize(&G[i][j]);
        }
        for (int j = 0; j <= i; j++) {
            dpe_clear(r[i][j]);
            dpe_clear(u[i][j]);
        }
        dpe_clear(lovasz[i]);
    }
    ibz_finalize(&X);
    ibz_finalize(&tmp);
    dpe_clear(delta_bar);
    dpe_clear(Xf);
    dpe_clear(tmpF);
}

/* ========== Lattice operations using Gram-based MLLL ========== */

void
quat_lattice_mul_mlll_gram(quat_lattice_t *res,
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

    quat_mlll_gram(&(res->basis), &rank, generators, 16, alg);
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
quat_lideal_create_mlll_gram(quat_left_ideal_t *lideal,
                             const quat_alg_elem_t *x,
                             const ibz_t *norm_N,
                             const quat_lattice_t *order,
                             const quat_alg_t *alg)
{
    assert(quat_order_is_maximal(order, alg));
    assert(!quat_alg_elem_is_zero(x));

    quat_lattice_t ON;
    quat_lattice_init(&ON);

    quat_lideal_create_principal(lideal, x, order, alg);

    ibz_mat_4x4_scalar_mul(&ON.basis, norm_N, &order->basis);
    ibz_copy(&ON.denom, &order->denom);

    quat_lattice_add_mlll_gram(&lideal->lattice, &lideal->lattice, &ON, alg);
    lideal->parent_order = order;
    quat_lideal_norm(lideal);

    quat_lattice_finalize(&ON);
}

void
quat_lattice_add_mlll_gram(quat_lattice_t *res,
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

    quat_mlll_gram(&(res->basis), &rank, generators, 8, alg);
    assert(rank > 0 && rank <= 4);

    ibz_mul(&(res->denom), &(lat1->denom), &(lat2->denom));
    quat_lattice_reduce_denom(res, res);

    ibz_mat_4x4_finalize(&tmp);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
}
