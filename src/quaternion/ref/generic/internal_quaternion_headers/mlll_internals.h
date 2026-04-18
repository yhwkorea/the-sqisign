#ifndef MLLL_INTERNALS_H
#define MLLL_INTERNALS_H

/** @file
 *
 * @brief Declarations for Modified LLL (Pohst 1987) algorithm
 *
 * MLLL takes a generating set of a lattice (possibly linearly dependent)
 * and outputs an LLL-reduced basis. Used to replace HNF in ideal operations.
 */

#include <quaternion.h>

/** @brief Maximum number of generators MLLL can accept */
#define MLLL_MAX_GENERATORS 16

/** @brief Lattice dimension (quaternion algebra = 4) */
#define MLLL_DIM 4

/**
 * @brief Modified LLL algorithm (Pohst 1987)
 *
 * Given g generating vectors of a lattice in Z^4, compute an LLL-reduced
 * basis of the lattice they span. Handles linear dependencies automatically.
 *
 * @param basis Output: 4x4 matrix whose columns are the LLL-reduced basis.
 *              Only the first `rank` columns are meaningful.
 * @param rank Output: rank of the lattice (number of basis vectors found)
 * @param generators Array of g vectors in Z^4 (each is ibz_vec_4_t)
 * @param g Number of generators (must be <= MLLL_MAX_GENERATORS)
 * @param alg The quaternion algebra (used for norm computation)
 */
void quat_mlll(ibz_mat_4x4_t *basis,
               int *rank,
               const ibz_vec_4_t *generators,
               int g,
               const quat_alg_t *alg);

/**
 * @brief Gram-based MLLL (L2-style: integer Gram matrix + dpe float GSO)
 *
 * Same interface as quat_mlll but keeps intermediate integers bounded by
 * max ||a_i||^2 (Lemma 3 of "Compact Quaternion Algorithms for SQIsign").
 * Built on top of the Nguyen-Stehle L2 algorithm (see l2.c) with Pohst's
 * generator-loading and dependent-vector removal added.
 */
void quat_mlll_gram(ibz_mat_4x4_t *basis,
                    int *rank,
                    const ibz_vec_4_t *generators,
                    int g,
                    const quat_alg_t *alg);

/**
 * @brief Lattice multiplication using MLLL instead of HNF
 *
 * Replaces quat_lattice_mul with MLLL-based approach.
 * Computes res = lat1 * lat2 using MLLL on 16 product generators.
 *
 * @param res Output: product lattice
 * @param lat1 First lattice
 * @param lat2 Second lattice
 * @param alg The quaternion algebra
 */
void quat_lattice_mul_mlll(quat_lattice_t *res,
                           const quat_lattice_t *lat1,
                           const quat_lattice_t *lat2,
                           const quat_alg_t *alg);

/**
 * @brief Lattice addition using MLLL instead of HNF
 *
 * Replaces quat_lattice_add with MLLL-based approach.
 *
 * @param res Output: sum lattice
 * @param lat1 First lattice
 * @param lat2 Second lattice
 */
void quat_lattice_add_mlll(quat_lattice_t *res,
                           const quat_lattice_t *lat1,
                           const quat_lattice_t *lat2,
                           const quat_alg_t *alg);

/**
 * @brief Lattice multiplication using Gram-based MLLL
 *
 * Same as quat_lattice_mul_mlll but uses quat_mlll_gram internally,
 * giving the Lemma 3 intermediate bitsize bound.
 */
void quat_lattice_mul_mlll_gram(quat_lattice_t *res,
                                const quat_lattice_t *lat1,
                                const quat_lattice_t *lat2,
                                const quat_alg_t *alg);

/**
 * @brief Lattice addition using Gram-based MLLL
 *
 * Same as quat_lattice_add_mlll but uses quat_mlll_gram internally.
 */
void quat_lattice_add_mlll_gram(quat_lattice_t *res,
                                const quat_lattice_t *lat1,
                                const quat_lattice_t *lat2,
                                const quat_alg_t *alg);

/**
 * @brief Left ideal creation using Gram-based MLLL for lattice add
 *
 * Same as quat_lideal_create but uses quat_lattice_add_mlll_gram internally.
 * This is the "compact" path of Algorithm 3 (RandomIdealGivenPrimeNorm):
 * the final "O·x + O·N" step becomes MLLL on 8 integer generators.
 */
void quat_lideal_create_mlll_gram(quat_left_ideal_t *lideal,
                                  const quat_alg_elem_t *x,
                                  const ibz_t *norm_N,
                                  const quat_lattice_t *order,
                                  const quat_alg_t *alg);

#endif
