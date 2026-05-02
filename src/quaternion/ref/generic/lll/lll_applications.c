#include <quaternion.h>
#include <internal.h>
#include "lll_internals.h"
#include "mlll_internals.h"

void
quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced,
                         ibz_mat_4x4_t *gram,
                         const quat_left_ideal_t *lideal,
                         const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal->parent_order), alg));
    ibz_t gram_corrector;
    ibz_init(&gram_corrector);
    ibz_mul(&gram_corrector, &(lideal->lattice.denom), &(lideal->lattice.denom));
    quat_lideal_class_gram(gram, lideal, alg);
    ibz_mat_4x4_copy(reduced, &(lideal->lattice.basis));
    quat_lll_core(gram, reduced);
    ibz_mat_4x4_scalar_mul(gram, &gram_corrector, gram);
    for (int i = 0; i < 4; i++) {
        ibz_div_2exp(&((*gram)[i][i]), &((*gram)[i][i]), 1);
        for (int j = i + 1; j < 4; j++) {
            ibz_set(&((*gram)[i][j]), 0);
        }
    }
    ibz_finalize(&gram_corrector);
}

/* ---------------------------------------------------------------------------
 * MLLL_GRAM variant of quat_lideal_reduce_basis
 *
 * Same external contract as quat_lideal_reduce_basis (Phase 3-1). Internally
 * loads the 4 columns of `lideal->lattice.basis` as MLLL generators and
 * drives the reduction with `quat_mlll_gram` (Lemma 3 width). The output
 * `gram` matches the HNF version's contract bit-for-bit: lower-triangular
 * Gram of the reduced columns divided by `lideal->norm`, with the diagonal
 * already halved (so diagonal entries are nrd(b_i)/norm) and the strict
 * upper triangle zeroed.
 * --------------------------------------------------------------------------- */
void
quat_lideal_reduce_basis_mlll_gram(ibz_mat_4x4_t *reduced,
                                   ibz_mat_4x4_t *gram,
                                   const quat_left_ideal_t *lideal,
                                   const quat_alg_t *alg)
{
    assert(quat_order_is_maximal((lideal->parent_order), alg));

    /* 4 columns of stored basis -> MLLL generators */
    ibz_vec_4_t generators[4];
    for (int i = 0; i < 4; i++)
        ibz_vec_4_init(&generators[i]);
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
            ibz_copy(&(generators[j][i]), &((lideal->lattice.basis)[i][j]));
        }
    }

    int rank = 0;
    quat_mlll_gram(reduced, &rank, generators, 4, alg);
    assert(rank == 4);

    /* Compute Gram of `reduced` via quat_lattice_gram on a synthetic lattice
     * (denom irrelevant: quat_lattice_gram uses only the basis matrix entries). */
    quat_lattice_t tmp_lat;
    quat_lattice_init(&tmp_lat);
    ibz_mat_4x4_copy(&(tmp_lat.basis), reduced);
    ibz_set(&(tmp_lat.denom), 1);

    quat_lattice_gram(gram, &tmp_lat, alg);

    /* quat_lattice_gram returns 2*nrd_form on the basis matrix entries, i.e.
     * 2 * nrd_form(reduced) = 2 * d^2 * nrd_form(actual basis). HNF version
     * ends with 2 * nrd_form(reduced) / norm before halving the diagonal.
     * Match that contract: divide by `lideal->norm`. */
    ibz_t rmd;
    ibz_init(&rmd);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <= i; j++) {
            ibz_div(&((*gram)[i][j]), &rmd, &((*gram)[i][j]), &(lideal->norm));
            assert(ibz_is_zero(&rmd));
        }
    }
    /* halve diagonal */
    for (int i = 0; i < 4; i++) {
        ibz_div_2exp(&((*gram)[i][i]), &((*gram)[i][i]), 1);
    }
    /* zero strict upper-triangular */
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            ibz_set(&((*gram)[i][j]), 0);
        }
    }

    ibz_finalize(&rmd);
    quat_lattice_finalize(&tmp_lat);
    for (int i = 0; i < 4; i++)
        ibz_vec_4_finalize(&generators[i]);
}

void
quat_lideal_lideal_mul_reduced(quat_left_ideal_t *prod,
                               ibz_mat_4x4_t *gram,
                               const quat_left_ideal_t *lideal1,
                               const quat_left_ideal_t *lideal2,
                               const quat_alg_t *alg)
{
    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);

    quat_lattice_mul(&(prod->lattice), &(lideal1->lattice), &(lideal2->lattice), alg);
    prod->parent_order = lideal1->parent_order;
    quat_lideal_norm(prod);
    quat_lideal_reduce_basis(&red, gram, prod, alg);
    ibz_mat_4x4_copy(&(prod->lattice.basis), &red);

    ibz_mat_4x4_finalize(&red);
}

int
quat_lideal_prime_norm_reduced_equivalent(quat_left_ideal_t *lideal,
                                          const quat_alg_t *alg,
                                          const int primality_num_iter,
                                          const int equiv_bound_coeff)
{
    ibz_mat_4x4_t gram, red;
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);

    int found = 0;

    // computing the reduced basis
    quat_lideal_reduce_basis(&red, &gram, lideal, alg);

    quat_alg_elem_t new_alpha;
    quat_alg_elem_init(&new_alpha);
    ibz_t tmp, remainder, adjusted_norm;
    ibz_init(&tmp);
    ibz_init(&remainder);
    ibz_init(&adjusted_norm);

    ibz_mul(&adjusted_norm, &lideal->lattice.denom, &lideal->lattice.denom);

    int ctr = 0;

    // equiv_num_iter = (2 * equiv_bound_coeff + 1)^4
    assert(equiv_bound_coeff < (1 << 20));
    int equiv_num_iter = (2 * equiv_bound_coeff + 1);
    equiv_num_iter = equiv_num_iter * equiv_num_iter;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;

    while (!found && ctr < equiv_num_iter) {
        ctr++;
        // we select our linear combination at random
        ibz_rand_interval_minm_m(&new_alpha.coord[0], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[1], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[2], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[3], equiv_bound_coeff);

        // computation of the norm of the vector sampled
        quat_qf_eval(&tmp, &gram, &new_alpha.coord);

        // compute the norm of the equivalent ideal
        // can be improved by removing the power of two first and the odd part only if the trial
        // division failed (this should always be called on an ideal of norm 2^x * N for some
        // big prime N )
        ibz_div(&tmp, &remainder, &tmp, &adjusted_norm);

        // debug : check that the remainder is zero
        assert(ibz_is_zero(&remainder));

        // pseudo-primality test
        if (ibz_probab_prime(&tmp, primality_num_iter)) {

            // computes the generator using a matrix multiplication
            ibz_mat_4x4_eval(&new_alpha.coord, &red, &new_alpha.coord);
            ibz_copy(&new_alpha.denom, &lideal->lattice.denom);
            assert(quat_lattice_contains(NULL, &lideal->lattice, &new_alpha));

            quat_alg_conj(&new_alpha, &new_alpha);
            ibz_mul(&new_alpha.denom, &new_alpha.denom, &lideal->norm);
            quat_lideal_mul(lideal, lideal, &new_alpha, alg);
            assert(ibz_probab_prime(&lideal->norm, primality_num_iter));

            found = 1;
            break;
        }
    }
    assert(found);

    ibz_finalize(&tmp);
    ibz_finalize(&remainder);
    ibz_finalize(&adjusted_norm);
    quat_alg_elem_finalize(&new_alpha);

    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);

    return found;
}

/* ---------------------------------------------------------------------------
 * MLLL_GRAM variant of quat_lideal_prime_norm_reduced_equivalent
 *
 * Same external contract and return semantics as the HNF version. Identical
 * random-sampling loop and primality handshake; the only structural change is
 * that the initial reduced basis is computed via
 * `quat_lideal_reduce_basis_mlll_gram` instead of `quat_lideal_reduce_basis`.
 *
 * Drop-in replacement for the HNF version in keygen / sign Alg 4 hot paths
 * once the dispatcher routing (Phase 3-2) is wired up.
 * --------------------------------------------------------------------------- */
int
quat_lideal_prime_norm_reduced_equivalent_mlll_gram(quat_left_ideal_t *lideal,
                                                    const quat_alg_t *alg,
                                                    const int primality_num_iter,
                                                    const int equiv_bound_coeff)
{
    ibz_mat_4x4_t gram, red;
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);

    int found = 0;

    /* MLLL_GRAM-driven reduced basis (only line that differs from HNF) */
    quat_lideal_reduce_basis_mlll_gram(&red, &gram, lideal, alg);

    quat_alg_elem_t new_alpha;
    quat_alg_elem_init(&new_alpha);
    ibz_t tmp, remainder, adjusted_norm;
    ibz_init(&tmp);
    ibz_init(&remainder);
    ibz_init(&adjusted_norm);

    ibz_mul(&adjusted_norm, &lideal->lattice.denom, &lideal->lattice.denom);

    int ctr = 0;

    /* equiv_num_iter = (2 * equiv_bound_coeff + 1)^4 */
    assert(equiv_bound_coeff < (1 << 20));
    int equiv_num_iter = (2 * equiv_bound_coeff + 1);
    equiv_num_iter = equiv_num_iter * equiv_num_iter;
    equiv_num_iter = equiv_num_iter * equiv_num_iter;

    while (!found && ctr < equiv_num_iter) {
        ctr++;
        ibz_rand_interval_minm_m(&new_alpha.coord[0], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[1], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[2], equiv_bound_coeff);
        ibz_rand_interval_minm_m(&new_alpha.coord[3], equiv_bound_coeff);

        quat_qf_eval(&tmp, &gram, &new_alpha.coord);

        ibz_div(&tmp, &remainder, &tmp, &adjusted_norm);
        assert(ibz_is_zero(&remainder));

        if (ibz_probab_prime(&tmp, primality_num_iter)) {
            ibz_mat_4x4_eval(&new_alpha.coord, &red, &new_alpha.coord);
            ibz_copy(&new_alpha.denom, &lideal->lattice.denom);
            assert(quat_lattice_contains(NULL, &lideal->lattice, &new_alpha));

            quat_alg_conj(&new_alpha, &new_alpha);
            ibz_mul(&new_alpha.denom, &new_alpha.denom, &lideal->norm);
            quat_lideal_mul(lideal, lideal, &new_alpha, alg);
            assert(ibz_probab_prime(&lideal->norm, primality_num_iter));

            found = 1;
            break;
        }
    }
    assert(found);

    ibz_finalize(&tmp);
    ibz_finalize(&remainder);
    ibz_finalize(&adjusted_norm);
    quat_alg_elem_finalize(&new_alpha);

    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);

    return found;
}
