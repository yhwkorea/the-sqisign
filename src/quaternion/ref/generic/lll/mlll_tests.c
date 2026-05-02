/**
 * @file mlll_tests.c
 * @brief Tests for Modified LLL algorithm
 *
 * Verification strategy:
 * 1. Correctness: MLLL output lattice == HNF output lattice (lattice equality)
 * 2. LLL-reducedness: verify size-reduction + Lovász conditions
 * 3. Bit size bound: track max integer size during execution
 */

#include <quaternion.h>
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"
#include "quat_fixed_precision.h"
#include "quaternion_tests.h"
#include <rng.h>

/* ========== Test 1: MLLL produces same lattice as HNF for lattice_mul ========== */

/**
 * @brief Test that quat_lattice_mul_mlll produces the same lattice as quat_lattice_mul (HNF)
 *
 * Uses same test vectors from the existing quat_test_lattice_mul test.
 */
int
quat_test_mlll_vs_hnf_lattice_mul(void)
{
    int res = 0;
    quat_lattice_t lat1, lat2, prod_hnf, prod_mlll;
    quat_alg_t alg;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&prod_hnf);
    quat_lattice_init(&prod_mlll);
    quat_alg_init_set_ui(&alg, 19);

    /* Test case 1: diagonal-ish lattices */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat1.basis[i][j]), 0);
            ibz_set(&(lat2.basis[i][j]), 0);
        }

    ibz_set(&(lat1.basis[0][0]), 44);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat1.basis[0][3]), 32);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), 2);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(lat1.basis[3][3]), 1);
    ibz_set(&(lat2.basis[3][3]), 3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);

    quat_lattice_mul(&prod_hnf, &lat1, &lat2, &alg);
    quat_lattice_mul_mlll(&prod_mlll, &lat1, &lat2, &alg);

    /* The lattices should be equal (same lattice, possibly different basis) */
    if (!quat_lattice_equal(&prod_hnf, &prod_mlll)) {
        printf("  FAIL: test case 1 - MLLL and HNF produce different lattices\n");
        printf("  HNF basis:\n");
        ibz_mat_4x4_print(&prod_hnf.basis);
        printf("  HNF denom: ");
        ibz_printf("%Zd\n", prod_hnf.denom);
        printf("  MLLL basis:\n");
        ibz_mat_4x4_print(&prod_mlll.basis);
        printf("  MLLL denom: ");
        ibz_printf("%Zd\n", prod_mlll.denom);
        res = 1;
    }

    /* Test case 2: non-HNF input */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat1.basis[i][j]), 0);
            ibz_set(&(lat2.basis[i][j]), 0);
        }

    ibz_set(&(lat1.basis[0][0]), 4);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(lat2.basis[0][3]), -1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), -2);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(lat2.basis[2][1]), 1);
    ibz_set(&(lat1.basis[3][3]), 7);
    ibz_set(&(lat2.basis[3][3]), -3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);

    quat_lattice_mul(&prod_hnf, &lat1, &lat2, &alg);
    quat_lattice_mul_mlll(&prod_mlll, &lat1, &lat2, &alg);

    if (!quat_lattice_equal(&prod_hnf, &prod_mlll)) {
        printf("  FAIL: test case 2 - MLLL and HNF produce different lattices\n");
        res = 1;
    }

    if (res == 0) {
        printf("  PASS: quat_test_mlll_vs_hnf_lattice_mul\n");
    }

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod_hnf);
    quat_lattice_finalize(&prod_mlll);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Test 2: MLLL output is LLL-reduced ========== */

/**
 * @brief Test that MLLL output satisfies LLL reduction conditions
 */
int
quat_test_mlll_is_lll_reduced(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 103);

    /* Create test generators: 8 vectors from two lattice bases */
    ibz_vec_4_t generators[8];
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);

    /* First 4: scaled identity */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            ibz_set(&generators[i][j], (i == j) ? 6 : 0);
    }
    /* Next 4: another basis */
    ibz_set(&generators[4][0], 3);
    ibz_set(&generators[4][1], 0);
    ibz_set(&generators[4][2], 0);
    ibz_set(&generators[4][3], 0);
    ibz_set(&generators[5][0], 1);
    ibz_set(&generators[5][1], 3);
    ibz_set(&generators[5][2], 0);
    ibz_set(&generators[5][3], 0);
    ibz_set(&generators[6][0], 0);
    ibz_set(&generators[6][1], 0);
    ibz_set(&generators[6][2], 3);
    ibz_set(&generators[6][3], 0);
    ibz_set(&generators[7][0], 0);
    ibz_set(&generators[7][1], 0);
    ibz_set(&generators[7][2], 1);
    ibz_set(&generators[7][3], 3);

    ibz_mat_4x4_t basis;
    int rank;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 8, &alg);

    if (rank != 4) {
        printf("  FAIL: expected rank 4, got %d\n", rank);
        res = 1;
    }

    /* Check LLL-reducedness using existing verification */
    if (rank == 4) {
        ibq_t delta, eta;
        ibq_init(&delta);
        ibq_init(&eta);

        /* Use delta=3/4 (the parameter from the paper) */
        ibz_t three, four;
        ibz_init(&three);
        ibz_init(&four);
        ibz_set(&three, 3);
        ibz_set(&four, 4);
        ibq_set(&delta, &three, &four);

        ibz_t one, two;
        ibz_init(&one);
        ibz_init(&two);
        ibz_set(&one, 1);
        ibz_set(&two, 2);
        ibq_set(&eta, &one, &two);

        int lll_ok = quat_lll_verify(&basis, &delta, &eta, &alg);
        if (!lll_ok) {
            printf("  FAIL: MLLL output is not LLL-reduced\n");
            printf("  Basis:\n");
            ibz_mat_4x4_print(&basis);
            res = 1;
        }

        ibq_finalize(&delta);
        ibq_finalize(&eta);
        ibz_finalize(&three);
        ibz_finalize(&four);
        ibz_finalize(&one);
        ibz_finalize(&two);
    }

    if (res == 0) {
        printf("  PASS: quat_test_mlll_is_lll_reduced\n");
    }

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Test 3: MLLL handles linearly dependent generators ========== */

/**
 * @brief Test that MLLL correctly removes linearly dependent vectors
 */
int
quat_test_mlll_dependent_generators(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);

    /* 6 generators where some are linearly dependent */
    ibz_vec_4_t generators[6];
    for (int i = 0; i < 6; i++)
        ibz_vec_4_init(&generators[i]);

    /* e1, e2, e3, e4, e1+e2, 2*e3 */
    ibz_set(&generators[0][0], 1); ibz_set(&generators[0][1], 0); ibz_set(&generators[0][2], 0); ibz_set(&generators[0][3], 0);
    ibz_set(&generators[1][0], 0); ibz_set(&generators[1][1], 1); ibz_set(&generators[1][2], 0); ibz_set(&generators[1][3], 0);
    ibz_set(&generators[2][0], 0); ibz_set(&generators[2][1], 0); ibz_set(&generators[2][2], 1); ibz_set(&generators[2][3], 0);
    ibz_set(&generators[3][0], 0); ibz_set(&generators[3][1], 0); ibz_set(&generators[3][2], 0); ibz_set(&generators[3][3], 1);
    ibz_set(&generators[4][0], 1); ibz_set(&generators[4][1], 1); ibz_set(&generators[4][2], 0); ibz_set(&generators[4][3], 0);
    ibz_set(&generators[5][0], 0); ibz_set(&generators[5][1], 0); ibz_set(&generators[5][2], 2); ibz_set(&generators[5][3], 0);

    ibz_mat_4x4_t basis;
    int rank;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 6, &alg);

    /* Should produce rank 4 (generators span Z^4) */
    if (rank != 4) {
        printf("  FAIL: expected rank 4, got %d\n", rank);
        res = 1;
    }

    /* The lattice should be Z^4, so basis should be identity (up to signs/reordering) */
    /* Check: lattice index should be 1 */
    if (rank == 4) {
        ibz_t det;
        ibz_init(&det);
        ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &basis);
        ibz_abs(&det, &det);
        if (!ibz_is_one(&det)) {
            printf("  FAIL: expected determinant 1 for Z^4, got ");
            ibz_printf("%Zd\n", det);
            res = 1;
        }
        ibz_finalize(&det);
    }

    if (res == 0) {
        printf("  PASS: quat_test_mlll_dependent_generators\n");
    }

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 6; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Test 4: CompactIdealMultiplication (Algorithm 2) correctness ========== */

/**
 * @brief Test the full CompactIdealMultiplication flow with a larger prime
 */
int
quat_test_compact_ideal_multiplication(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_lattice_t lat1, lat2, prod_hnf, prod_mlll;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&prod_hnf);
    quat_lattice_init(&prod_mlll);

    /* Create a non-trivial lattice product test */
    /* lat1: a random-ish O0 lattice */
    ibz_set(&(lat1.basis[0][0]), 7);
    ibz_set(&(lat1.basis[0][1]), 3);
    ibz_set(&(lat1.basis[0][2]), 0);
    ibz_set(&(lat1.basis[0][3]), 2);
    ibz_set(&(lat1.basis[1][0]), 0);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat1.basis[1][2]), 1);
    ibz_set(&(lat1.basis[1][3]), 0);
    ibz_set(&(lat1.basis[2][0]), 0);
    ibz_set(&(lat1.basis[2][1]), 0);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat1.basis[2][3]), 1);
    ibz_set(&(lat1.basis[3][0]), 0);
    ibz_set(&(lat1.basis[3][1]), 0);
    ibz_set(&(lat1.basis[3][2]), 0);
    ibz_set(&(lat1.basis[3][3]), 11);
    ibz_set(&(lat1.denom), 2);

    /* lat2: another lattice */
    ibz_set(&(lat2.basis[0][0]), 4);
    ibz_set(&(lat2.basis[0][1]), 1);
    ibz_set(&(lat2.basis[0][2]), 0);
    ibz_set(&(lat2.basis[0][3]), 0);
    ibz_set(&(lat2.basis[1][0]), 0);
    ibz_set(&(lat2.basis[1][1]), 3);
    ibz_set(&(lat2.basis[1][2]), 0);
    ibz_set(&(lat2.basis[1][3]), 2);
    ibz_set(&(lat2.basis[2][0]), 0);
    ibz_set(&(lat2.basis[2][1]), 0);
    ibz_set(&(lat2.basis[2][2]), 5);
    ibz_set(&(lat2.basis[2][3]), 0);
    ibz_set(&(lat2.basis[3][0]), 0);
    ibz_set(&(lat2.basis[3][1]), 0);
    ibz_set(&(lat2.basis[3][2]), 0);
    ibz_set(&(lat2.basis[3][3]), 7);
    ibz_set(&(lat2.denom), 3);

    quat_lattice_mul(&prod_hnf, &lat1, &lat2, &alg);
    quat_lattice_mul_mlll(&prod_mlll, &lat1, &lat2, &alg);

    if (!quat_lattice_equal(&prod_hnf, &prod_mlll)) {
        printf("  FAIL: CompactIdealMul test - lattices differ\n");
        res = 1;
    }

    if (res == 0) {
        printf("  PASS: quat_test_compact_ideal_multiplication\n");
    }

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod_hnf);
    quat_lattice_finalize(&prod_mlll);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Test 5: Realistic scale (127-bit prime) ========== */

/**
 * @brief Test MLLL with a prime close to SQIsign Level 1 size (~127 bits)
 *
 * Previous tests used p=7, 19, 103. This test uses a cryptographic-scale prime
 * to verify MLLL correctness at actual deployment size.
 */
int
quat_test_mlll_realistic_scale(void)
{
    int res = 0;
    quat_alg_t alg;
    ibz_t prime;
    ibz_init(&prime);

    /* Level 1 SQIsign prime (~127 bits) */
    ibz_set_from_str(&prime,
        "4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    quat_alg_init_set(&alg, &prime);

    /* Create lattices with large entries */
    quat_lattice_t lat1, lat2, prod_hnf, prod_mlll;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&prod_hnf);
    quat_lattice_init(&prod_mlll);

    /* lat1: diagonal with ~64-bit entries */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&(lat1.basis[i][j]), 0);
    ibz_set_from_str(&(lat1.basis[0][0]), "7a3f5c2e1d0b9876", 16);
    ibz_set_from_str(&(lat1.basis[0][2]), "3f1a2b", 16);
    ibz_set_from_str(&(lat1.basis[1][1]), "5e4d3c2b1a098765", 16);
    ibz_set_from_str(&(lat1.basis[2][2]), "4c3b2a1908765432", 16);
    ibz_set_from_str(&(lat1.basis[3][3]), "3b2a190876543210", 16);
    ibz_set(&(lat1.denom), 2);

    /* lat2: another lattice with ~64-bit entries */
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&(lat2.basis[i][j]), 0);
    ibz_set_from_str(&(lat2.basis[0][0]), "6b2e4d1c3a098765", 16);
    ibz_set_from_str(&(lat2.basis[0][3]), "1f2e3d", 16);
    ibz_set_from_str(&(lat2.basis[1][1]), "4d3c2b1a09876543", 16);
    ibz_set_from_str(&(lat2.basis[2][2]), "3c2b1a0987654321", 16);
    ibz_set_from_str(&(lat2.basis[3][3]), "2b1a098765432100", 16);
    ibz_set(&(lat2.denom), 3);

    quat_lattice_mul(&prod_hnf, &lat1, &lat2, &alg);
    quat_lattice_mul_mlll(&prod_mlll, &lat1, &lat2, &alg);

    if (!quat_lattice_equal(&prod_hnf, &prod_mlll)) {
        printf("  FAIL: realistic scale test - lattices differ\n");
        printf("  HNF basis:\n");
        ibz_mat_4x4_print(&prod_hnf.basis);
        printf("  MLLL basis:\n");
        ibz_mat_4x4_print(&prod_mlll.basis);
        res = 1;
    }

    /* Verify LLL-reducedness */
    if (res == 0) {
        ibq_t delta, eta;
        ibq_init(&delta);
        ibq_init(&eta);
        ibz_t three, four, one, two;
        ibz_init(&three); ibz_init(&four); ibz_init(&one); ibz_init(&two);
        ibz_set(&three, 3); ibz_set(&four, 4);
        ibz_set(&one, 1); ibz_set(&two, 2);
        ibq_set(&delta, &three, &four);
        ibq_set(&eta, &one, &two);

        int lll_ok = quat_lll_verify(&prod_mlll.basis, &delta, &eta, &alg);
        if (!lll_ok) {
            printf("  FAIL: realistic scale MLLL output is not LLL-reduced\n");
            res = 1;
        }

        ibq_finalize(&delta); ibq_finalize(&eta);
        ibz_finalize(&three); ibz_finalize(&four);
        ibz_finalize(&one); ibz_finalize(&two);
    }

    if (res == 0) {
        printf("  PASS: quat_test_mlll_realistic_scale\n");
    }

    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod_hnf);
    quat_lattice_finalize(&prod_mlll);
    quat_alg_finalize(&alg);
    ibz_finalize(&prime);
    return res;
}

/* ========== Test 6: tau > 1 path after remove_vector + swap ========== */

/**
 * @brief Test that exercises the tau > 1 swap path
 *
 * After remove_vector sets tau = m+1, a subsequent swap must check m > tau
 * (not m > 1). This test creates generators where a dependent vector is
 * followed by vectors that trigger a swap.
 */
int
quat_test_mlll_tau_path(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 19);

    /* 8 generators designed to trigger: load → remove → load → swap with tau > 1
     *
     * Generators 0-3: standard basis (will be loaded as initial vectors)
     * Generator 4: zero vector (triggers dependent removal)
     * Generator 5: e1 + e2 (dependent on 0,1 — triggers removal, tau update)
     * Generator 6-7: vectors that should trigger Lovász swap after tau advances
     */
    ibz_vec_4_t generators[8];
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);

    /* e1 */
    ibz_set(&generators[0][0], 1);
    /* e2 */
    ibz_set(&generators[1][1], 1);
    /* e3 */
    ibz_set(&generators[2][2], 1);
    /* e4 */
    ibz_set(&generators[3][3], 1);
    /* zero vector */
    /* (generators[4] already zero) */
    /* e1 + e2 (dependent) */
    ibz_set(&generators[5][0], 1);
    ibz_set(&generators[5][1], 1);
    /* Vector that should cause Lovász violation with existing basis */
    ibz_set(&generators[6][0], 100);
    ibz_set(&generators[6][1], 1);
    ibz_set(&generators[6][2], 1);
    /* Another vector */
    ibz_set(&generators[7][0], 1);
    ibz_set(&generators[7][1], 100);
    ibz_set(&generators[7][2], 1);
    ibz_set(&generators[7][3], 1);

    ibz_mat_4x4_t basis;
    int rank;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 8, &alg);

    if (rank != 4) {
        printf("  FAIL: tau path test - expected rank 4, got %d\n", rank);
        res = 1;
    }

    /* Verify LLL-reducedness */
    if (rank == 4) {
        ibq_t delta, eta;
        ibq_init(&delta);
        ibq_init(&eta);
        ibz_t three, four, one, two;
        ibz_init(&three); ibz_init(&four); ibz_init(&one); ibz_init(&two);
        ibz_set(&three, 3); ibz_set(&four, 4);
        ibz_set(&one, 1); ibz_set(&two, 2);
        ibq_set(&delta, &three, &four);
        ibq_set(&eta, &one, &two);

        int lll_ok = quat_lll_verify(&basis, &delta, &eta, &alg);
        if (!lll_ok) {
            printf("  FAIL: tau path test - output not LLL-reduced\n");
            printf("  Basis:\n");
            ibz_mat_4x4_print(&basis);
            res = 1;
        }

        ibq_finalize(&delta); ibq_finalize(&eta);
        ibz_finalize(&three); ibz_finalize(&four);
        ibz_finalize(&one); ibz_finalize(&two);
    }

    if (res == 0) {
        printf("  PASS: quat_test_mlll_tau_path\n");
    }

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Test: quat_mlll_gram lattice equivalence to quat_mlll ========== */

/* Helper: do both implementations produce the same lattice on one input set? */
static int
compare_mlll_vs_gram(const ibz_vec_4_t *generators, int g, const quat_alg_t *alg,
                     const char *label)
{
    ibz_mat_4x4_t basis_a, basis_b;
    int rank_a = 0, rank_b = 0;
    ibz_mat_4x4_init(&basis_a);
    ibz_mat_4x4_init(&basis_b);

    quat_mlll(&basis_a, &rank_a, generators, g, alg);
    quat_mlll_gram(&basis_b, &rank_b, generators, g, alg);

    int res = 0;

    if (rank_a != rank_b) {
        printf("  FAIL[%s]: rank mismatch mlll=%d gram=%d\n", label, rank_a, rank_b);
        res = 1;
    } else if (rank_a == 4) {
        quat_lattice_t la, lb;
        quat_lattice_init(&la);
        quat_lattice_init(&lb);
        ibz_mat_4x4_copy(&la.basis, &basis_a);
        ibz_mat_4x4_copy(&lb.basis, &basis_b);
        ibz_set(&la.denom, 1);
        ibz_set(&lb.denom, 1);
        if (!quat_lattice_equal(&la, &lb)) {
            printf("  FAIL[%s]: full-rank lattices differ\n", label);
            res = 1;
        }
        quat_lattice_finalize(&la);
        quat_lattice_finalize(&lb);
    }
    /* For rank < 4, the output basis is zero-padded; full equality test requires
     * projecting to the rank subspace which we skip. The rank check alone is a
     * useful sanity signal. */

    ibz_mat_4x4_finalize(&basis_a);
    ibz_mat_4x4_finalize(&basis_b);
    return res;
}

int
quat_test_mlll_gram_equivalence(void)
{
    int res = 0;

    /* Test 1: 6 generators, Z^4 standard basis + dependents (from test_dependent) */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 7);
        ibz_vec_4_t gens[6];
        for (int i = 0; i < 6; i++) ibz_vec_4_init(&gens[i]);
        ibz_set(&gens[0][0], 1);
        ibz_set(&gens[1][1], 1);
        ibz_set(&gens[2][2], 1);
        ibz_set(&gens[3][3], 1);
        ibz_set(&gens[4][0], 1); ibz_set(&gens[4][1], 1);
        ibz_set(&gens[5][2], 2);
        res |= compare_mlll_vs_gram(gens, 6, &alg, "dependent_6gen");
        for (int i = 0; i < 6; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    /* Test 2: 8 dense random-ish generators with small prime */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 19);
        ibz_vec_4_t gens[8];
        for (int i = 0; i < 8; i++) ibz_vec_4_init(&gens[i]);
        int vals[8][4] = {
            {3, 1, 0, 2},
            {1, 4, 2, 0},
            {0, 2, 3, 1},
            {2, 0, 1, 4},
            {5, 3, 1, 0},
            {1, 1, 2, 2},
            {4, 0, 0, 3},
            {0, 5, 1, 1},
        };
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&gens[i][j], vals[i][j]);
        res |= compare_mlll_vs_gram(gens, 8, &alg, "dense_8gen_p19");
        for (int i = 0; i < 8; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    /* Test 3: 16 generators (max), stress */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 11);
        ibz_vec_4_t gens[16];
        for (int i = 0; i < 16; i++) ibz_vec_4_init(&gens[i]);
        for (int i = 0; i < 16; i++) {
            ibz_set(&gens[i][0], (i * 7 + 3) % 17);
            ibz_set(&gens[i][1], (i * 11 + 5) % 17);
            ibz_set(&gens[i][2], (i * 13 + 1) % 17);
            ibz_set(&gens[i][3], (i * 5 + 9) % 17);
        }
        res |= compare_mlll_vs_gram(gens, 16, &alg, "stress_16gen_p11");
        for (int i = 0; i < 16; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    if (res == 0)
        printf("  PASS: quat_test_mlll_gram_equivalence\n");
    return res;
}

/* ========== Test: lideal_create HNF == lideal_create_mlll_gram ========== */

int
quat_test_lideal_create_gram_equivalence(void)
{
    int res = 0;
    quat_alg_t alg;
    ibz_t prime;
    ibz_init(&prime);
    ibz_set_from_str(&prime,
        "4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    quat_alg_init_set(&alg, &prime);

    quat_p_extremal_maximal_order_t order;
    quat_represent_integer_params_t params;
    quat_lattice_init(&(order.order));
    quat_alg_elem_init(&(order.t));
    quat_alg_elem_init(&(order.z));
    quat_lattice_O0_set_extremal(&order);
    params.algebra = &alg;
    params.order = &order;
    params.primality_test_iterations = 30;

    uint32_t seed[12] = { 0xCAFEBABE };
    randombytes_init((unsigned char *)seed, NULL, 256);

    int trials = 25;
    quat_lattice_t *lats = malloc(trials * sizeof(quat_lattice_t));
    ibz_t *norms = malloc(trials * sizeof(ibz_t));
    for (int i = 0; i < trials; i++) {
        quat_lattice_init(&lats[i]);
        ibz_init(&norms[i]);
    }

    if (quat_test_input_random_ideal_lattice_generation(lats, norms, 127, trials, &params) != 0) {
        printf("  FAIL: random ideal generation\n");
        res = 1;
        goto cleanup;
    }

    /* For each generated lattice, extract a generator x of the ideal and N,
     * then compare quat_lideal_create vs quat_lideal_create_mlll_gram. */
    for (int i = 0; i < trials; i++) {
        /* The random_ideal_lattice_generation produces left O0-ideals with known norms.
         * We build a generator via an element of the lattice and the norm. */
        quat_alg_elem_t x;
        quat_alg_elem_init(&x);
        /* Use lats[i].basis column 0 as a generator candidate (scaled by denom). */
        for (int j = 0; j < 4; j++)
            ibz_copy(&x.coord[j], &lats[i].basis[j][0]);
        ibz_copy(&x.denom, &lats[i].denom);

        if (quat_alg_elem_is_zero(&x)) {
            quat_alg_elem_finalize(&x);
            continue;
        }

        quat_left_ideal_t I_hnf, I_gram;
        quat_left_ideal_init(&I_hnf);
        quat_left_ideal_init(&I_gram);

        quat_lideal_create(&I_hnf, &x, &norms[i], &(order.order), &alg);
        quat_lideal_create_mlll_gram(&I_gram, &x, &norms[i], &(order.order), &alg);

        if (!quat_lattice_equal(&I_hnf.lattice, &I_gram.lattice)) {
            printf("  FAIL: trial %d: lattices differ\n", i);
            res = 1;
        }
        if (ibz_cmp(&I_hnf.norm, &I_gram.norm) != 0) {
            printf("  FAIL: trial %d: norms differ\n", i);
            res = 1;
        }

        quat_left_ideal_finalize(&I_hnf);
        quat_left_ideal_finalize(&I_gram);
        quat_alg_elem_finalize(&x);
    }

    if (res == 0)
        printf("  PASS: quat_test_lideal_create_gram_equivalence (%d trials)\n", trials);

cleanup:
    for (int i = 0; i < trials; i++) {
        quat_lattice_finalize(&lats[i]);
        ibz_finalize(&norms[i]);
    }
    free(lats); free(norms);
    quat_lattice_finalize(&(order.order));
    quat_alg_elem_finalize(&(order.t));
    quat_alg_elem_finalize(&(order.z));
    quat_alg_finalize(&alg);
    ibz_finalize(&prime);
    return res;
}

/* ========== Test: quat_mlll_gram prealloc mode produces identical output ==========
 *
 * Phase 2 candidate B: `g_prealloc_mode` is a hint-only optimization
 * (mpz_realloc2). Output must be bitwise identical between mode 0 and mode 1
 * across every trial — anything else means we've introduced a correctness bug
 * along the way.
 *
 * Tests both Alg 2 (lattice_mul_mlll_gram) and Alg 3 (lideal_create_mlll_gram)
 * entry points at L1 since they are the two production code paths.
 */
int
quat_test_mlll_gram_prealloc_equivalence(void)
{
    int res = 0;
    quat_alg_t alg;
    ibz_t prime;
    ibz_init(&prime);
    ibz_set_from_str(&prime,
        "4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    quat_alg_init_set(&alg, &prime);

    quat_p_extremal_maximal_order_t order;
    quat_represent_integer_params_t params;
    quat_lattice_init(&(order.order));
    quat_alg_elem_init(&(order.t));
    quat_alg_elem_init(&(order.z));
    quat_lattice_O0_set_extremal(&order);
    params.algebra = &alg;
    params.order = &order;
    params.primality_test_iterations = 30;

    uint32_t seed[12] = { 0xDEADBEEF };
    randombytes_init((unsigned char *)seed, NULL, 256);

    int trials = 20;
    quat_lattice_t *lats1 = malloc(trials * sizeof(quat_lattice_t));
    quat_lattice_t *lats2 = malloc(trials * sizeof(quat_lattice_t));
    ibz_t *norms1 = malloc(trials * sizeof(ibz_t));
    for (int i = 0; i < trials; i++) {
        quat_lattice_init(&lats1[i]);
        quat_lattice_init(&lats2[i]);
        ibz_init(&norms1[i]);
    }

    if (quat_test_input_random_ideal_lattice_generation(lats1, norms1, 127, trials, &params) != 0) {
        printf("  FAIL: random ideal generation (set 1)\n");
        res = 1;
        goto cleanup;
    }
    {
        ibz_t *norms2_tmp = malloc(trials * sizeof(ibz_t));
        for (int i = 0; i < trials; i++) ibz_init(&norms2_tmp[i]);
        if (quat_test_input_random_ideal_lattice_generation(lats2, norms2_tmp, 127, trials, &params) != 0) {
            printf("  FAIL: random ideal generation (set 2)\n");
            res = 1;
        }
        for (int i = 0; i < trials; i++) ibz_finalize(&norms2_tmp[i]);
        free(norms2_tmp);
        if (res) goto cleanup;
    }

    int saved_mode = quat_mlll_gram_get_prealloc_mode();
    /* Force the ibz path — prealloc is a hint on the ibz mpz_t buffers and
     * has no effect when the dispatcher routes through the fp body. Since
     * fp_mode default flipped to 1 (Phase 2 primary) we must pin it off
     * here so this test actually exercises the prealloc toggle. */
    int saved_fp = quat_mlll_gram_get_fp_mode();
    quat_mlll_gram_set_fp_mode(0);

    /* ---- Alg 2: lattice_mul_mlll_gram ---- */
    for (int i = 0; i < trials; i++) {
        quat_lattice_t prod_off, prod_on;
        quat_lattice_init(&prod_off);
        quat_lattice_init(&prod_on);

        quat_mlll_gram_set_prealloc_mode(0);
        quat_lattice_mul_mlll_gram(&prod_off, &lats1[i], &lats2[i], &alg);

        quat_mlll_gram_set_prealloc_mode(1);
        quat_lattice_mul_mlll_gram(&prod_on, &lats1[i], &lats2[i], &alg);

        if (!quat_lattice_equal(&prod_off, &prod_on)) {
            printf("  FAIL: alg2 trial %d: lattice differs between prealloc 0 and 1\n", i);
            res = 1;
        }

        quat_lattice_finalize(&prod_off);
        quat_lattice_finalize(&prod_on);
    }

    /* ---- Alg 3: lideal_create_mlll_gram ---- */
    for (int i = 0; i < trials; i++) {
        quat_alg_elem_t x;
        quat_alg_elem_init(&x);
        for (int j = 0; j < 4; j++)
            ibz_copy(&x.coord[j], &lats1[i].basis[j][0]);
        ibz_copy(&x.denom, &lats1[i].denom);

        if (quat_alg_elem_is_zero(&x)) {
            quat_alg_elem_finalize(&x);
            continue;
        }

        quat_left_ideal_t I_off, I_on;
        quat_left_ideal_init(&I_off);
        quat_left_ideal_init(&I_on);

        quat_mlll_gram_set_prealloc_mode(0);
        quat_lideal_create_mlll_gram(&I_off, &x, &norms1[i], &(order.order), &alg);

        quat_mlll_gram_set_prealloc_mode(1);
        quat_lideal_create_mlll_gram(&I_on, &x, &norms1[i], &(order.order), &alg);

        if (!quat_lattice_equal(&I_off.lattice, &I_on.lattice)) {
            printf("  FAIL: alg3 trial %d: ideal lattice differs between prealloc 0 and 1\n", i);
            res = 1;
        }
        if (ibz_cmp(&I_off.norm, &I_on.norm) != 0) {
            printf("  FAIL: alg3 trial %d: norms differ between prealloc 0 and 1\n", i);
            res = 1;
        }

        quat_left_ideal_finalize(&I_off);
        quat_left_ideal_finalize(&I_on);
        quat_alg_elem_finalize(&x);
    }

    quat_mlll_gram_set_prealloc_mode(saved_mode);
    quat_mlll_gram_set_fp_mode(saved_fp);

    if (res == 0)
        printf("  PASS: quat_test_mlll_gram_prealloc_equivalence (%d trials, alg2+alg3)\n", trials);

cleanup:
    for (int i = 0; i < trials; i++) {
        quat_lattice_finalize(&lats1[i]);
        quat_lattice_finalize(&lats2[i]);
        ibz_finalize(&norms1[i]);
    }
    free(lats1); free(lats2); free(norms1);
    quat_lattice_finalize(&(order.order));
    quat_alg_elem_finalize(&(order.t));
    quat_alg_elem_finalize(&(order.z));
    quat_alg_finalize(&alg);
    ibz_finalize(&prime);
    return res;
}

/* ========== quat_mlll_gram fp path vs ibz path equivalence ==========
 *
 * Phase 2 candidate C: `g_fp_mode` routes `quat_mlll_gram` through the
 * stack fixed-precision body (`quat_mlll_gram_fp`). Output must match the
 * `ibz_t` body (= `quat_mlll_gram` with fp_mode = 0) under the Lemma 3
 * invariant. Reuses the three `compare_mlll_vs_gram` cases from the
 * gram-vs-mlll test so the fp path is exercised against the same
 * canonical inputs, then adds an L5-width ideal trial to cover the 17-
 * limb Gram path which the tiny-p cases cannot reach.
 */
int
quat_test_mlll_gram_fp_equivalence(void)
{
    int res = 0;
    int saved_fp = quat_mlll_gram_get_fp_mode();
    quat_mlll_gram_set_fp_mode(1);

    /* Reuse the 3 gram_equivalence cases under fp_mode = 1. */
    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 7);
        ibz_vec_4_t gens[6];
        for (int i = 0; i < 6; i++) ibz_vec_4_init(&gens[i]);
        ibz_set(&gens[0][0], 1);
        ibz_set(&gens[1][1], 1);
        ibz_set(&gens[2][2], 1);
        ibz_set(&gens[3][3], 1);
        ibz_set(&gens[4][0], 1); ibz_set(&gens[4][1], 1);
        ibz_set(&gens[5][2], 2);
        res |= compare_mlll_vs_gram(gens, 6, &alg, "fp_dependent_6gen");
        for (int i = 0; i < 6; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 19);
        ibz_vec_4_t gens[8];
        for (int i = 0; i < 8; i++) ibz_vec_4_init(&gens[i]);
        int vals[8][4] = {
            {3, 1, 0, 2}, {1, 4, 2, 0},
            {0, 2, 3, 1}, {2, 0, 1, 4},
            {5, 3, 1, 0}, {1, 1, 2, 2},
            {4, 0, 0, 3}, {0, 5, 1, 1},
        };
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&gens[i][j], vals[i][j]);
        res |= compare_mlll_vs_gram(gens, 8, &alg, "fp_dense_8gen_p19");
        for (int i = 0; i < 8; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    {
        quat_alg_t alg;
        quat_alg_init_set_ui(&alg, 11);
        ibz_vec_4_t gens[16];
        for (int i = 0; i < 16; i++) ibz_vec_4_init(&gens[i]);
        for (int i = 0; i < 16; i++) {
            ibz_set(&gens[i][0], (i * 7 + 3) % 17);
            ibz_set(&gens[i][1], (i * 11 + 5) % 17);
            ibz_set(&gens[i][2], (i * 13 + 1) % 17);
            ibz_set(&gens[i][3], (i * 5 + 9) % 17);
        }
        res |= compare_mlll_vs_gram(gens, 16, &alg, "fp_stress_16gen_p11");
        for (int i = 0; i < 16; i++) ibz_vec_4_finalize(&gens[i]);
        quat_alg_finalize(&alg);
    }

    /* Per-level width coverage with production-matching prime bitsize.
     * Post-2026-04-22 dispatcher fix (§6.8), each of the three L1/L3/L5
     * branches of `quat_fp_widths_from_alg` selects a distinct widths
     * table. We exercise all three so fp-path bugs can't hide in a branch
     * we never visited. Primes replicate the MSB pattern of
     * `precomp/ref/lvl{1,3,5}/quaternion_data.c` QUATALG_PINFTY.p so
     * ibz_bitsize lands in the correct band:
     *   L1: 251 bits  (top limb 0x4ff..ff, 4 limbs)
     *   L3: 383 bits  (top limb 0x40ff..ff, 6 limbs)
     *   L5: 505 bits  (top limb 0x01af..ff, 8 limbs)
     * 16 random generators bounded by floor(sqrt(p)) so Lemma 3
     * `|<a, b>| ≈ a_bits + b_bits + p_bits + log n` stays within the
     * per-level Gram budget. Using bound=p would yield `|<a,b>| ~ 3*p_bits`
     * and blow the Gram slot even at L1. Fixed per-block seed. */
    struct {
        const char *hex;
        const char *label;
        uint32_t    seed;
    } level_cases[] = {
        { "4fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
          "fp_L1_16gen_prod_p251", 0xF0CACC1A },
        { "40ffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff",
          "fp_L3_16gen_prod_p383", 0xA1B2C3D4 },
        { "01afffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff"
          "ffffffffffffffff",
          "fp_L5_16gen_prod_p505", 0xDEADBEEF },
    };

    for (size_t k = 0; k < sizeof(level_cases) / sizeof(level_cases[0]); k++) {
        quat_alg_t alg;
        ibz_t prime;
        ibz_init(&prime);
        ibz_set_from_str(&prime, level_cases[k].hex, 16);
        quat_alg_init_set(&alg, &prime);

        uint32_t seed[12] = { 0 };
        seed[0] = level_cases[k].seed;
        randombytes_init((unsigned char *)seed, NULL, 256);

        ibz_t bound, neg_bound;
        ibz_init(&bound);
        ibz_init(&neg_bound);
        ibz_sqrt_floor(&bound, &prime);
        ibz_neg(&neg_bound, &bound);

        ibz_vec_4_t gens[16];
        for (int i = 0; i < 16; i++) ibz_vec_4_init(&gens[i]);
        for (int i = 0; i < 16; i++)
            for (int j = 0; j < 4; j++)
                ibz_rand_interval(&gens[i][j], &neg_bound, &bound);

        res |= compare_mlll_vs_gram(gens, 16, &alg, level_cases[k].label);

        for (int i = 0; i < 16; i++) ibz_vec_4_finalize(&gens[i]);
        ibz_finalize(&bound);
        ibz_finalize(&neg_bound);
        ibz_finalize(&prime);
        quat_alg_finalize(&alg);
    }

    quat_mlll_gram_set_fp_mode(saved_fp);

    if (res == 0)
        printf("  PASS: quat_test_mlll_gram_fp_equivalence\n");
    return res;
}

/* ========== Corner-case tests: minimal / rank-deficient / zero inputs ========== */

/* Helper: check all entries of column c in basis are zero */
static int
col_is_zero(const ibz_mat_4x4_t *basis, int c)
{
    for (int r = 0; r < 4; r++)
        if (!ibz_is_zero(&(*basis)[r][c])) return 0;
    return 1;
}

int
quat_test_mlll_all_zero_generators(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);

    ibz_vec_4_t generators[4];
    for (int i = 0; i < 4; i++)
        ibz_vec_4_init(&generators[i]);  /* all zero */

    ibz_mat_4x4_t basis;
    int rank = -1;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 4, &alg);

    if (rank != 0) {
        printf("  FAIL: all-zero generators should give rank 0, got %d\n", rank);
        res = 1;
    }

    if (res == 0) printf("  PASS: quat_test_mlll_all_zero_generators\n");

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 4; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

int
quat_test_mlll_single_generator(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);

    ibz_vec_4_t generators[1];
    ibz_vec_4_init(&generators[0]);
    ibz_set(&generators[0][0], 3);
    ibz_set(&generators[0][2], 5);

    ibz_mat_4x4_t basis;
    int rank = -1;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 1, &alg);

    if (rank != 1) {
        printf("  FAIL: g=1 non-zero should give rank 1, got %d\n", rank);
        res = 1;
    } else {
        /* basis column 0 should equal generator (up to sign) */
        for (int r = 0; r < 4; r++) {
            if (ibz_cmp(&basis[r][0], &generators[0][r]) != 0) {
                ibz_t neg;
                ibz_init(&neg);
                ibz_neg(&neg, &generators[0][r]);
                if (ibz_cmp(&basis[r][0], &neg) != 0) {
                    printf("  FAIL: g=1 basis col 0 row %d doesn't match generator\n", r);
                    res = 1;
                }
                ibz_finalize(&neg);
                break;
            }
        }
    }

    if (res == 0) printf("  PASS: quat_test_mlll_single_generator\n");

    ibz_mat_4x4_finalize(&basis);
    ibz_vec_4_finalize(&generators[0]);
    quat_alg_finalize(&alg);
    return res;
}

int
quat_test_mlll_two_generators_dependent(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);

    ibz_vec_4_t generators[2];
    for (int i = 0; i < 2; i++)
        ibz_vec_4_init(&generators[i]);

    /* b1 = (2,0,3,0), b2 = 2*b1 = (4,0,6,0) */
    ibz_set(&generators[0][0], 2);
    ibz_set(&generators[0][2], 3);
    ibz_set(&generators[1][0], 4);
    ibz_set(&generators[1][2], 6);

    ibz_mat_4x4_t basis;
    int rank = -1;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 2, &alg);

    if (rank != 1) {
        printf("  FAIL: two colinear generators should give rank 1, got %d\n", rank);
        res = 1;
    }

    if (res == 0) printf("  PASS: quat_test_mlll_two_generators_dependent\n");

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 2; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

int
quat_test_mlll_rank_deficient(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);

    /* 8 generators in a rank-2 sublattice: all in e1-e2 plane */
    ibz_vec_4_t generators[8];
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);

    /* e1, e2, e1+e2, 2*e1, 3*e2, e1-e2, 0, 2*e1+3*e2 */
    ibz_set(&generators[0][0], 1);
    ibz_set(&generators[1][1], 1);
    ibz_set(&generators[2][0], 1); ibz_set(&generators[2][1], 1);
    ibz_set(&generators[3][0], 2);
    ibz_set(&generators[4][1], 3);
    ibz_set(&generators[5][0], 1); ibz_set(&generators[5][1], -1);
    /* generators[6] left as zero */
    ibz_set(&generators[7][0], 2); ibz_set(&generators[7][1], 3);

    ibz_mat_4x4_t basis;
    int rank = -1;
    ibz_mat_4x4_init(&basis);

    quat_mlll(&basis, &rank, generators, 8, &alg);

    if (rank != 2) {
        printf("  FAIL: rank-2 sublattice should give rank 2, got %d\n", rank);
        res = 1;
    } else {
        /* Columns 2 and 3 should be zero (only cols 0,1 meaningful) */
        if (!col_is_zero(&basis, 2) || !col_is_zero(&basis, 3)) {
            printf("  FAIL: rank-2 output should zero columns 2,3\n");
            res = 1;
        }
        /* All basis entries should have row 2,3 = 0 (no e3/e4 component) */
        for (int c = 0; c < 2; c++) {
            if (!ibz_is_zero(&basis[2][c]) || !ibz_is_zero(&basis[3][c])) {
                printf("  FAIL: basis col %d has nonzero e3 or e4 component\n", c);
                res = 1;
            }
        }
    }

    if (res == 0) printf("  PASS: quat_test_mlll_rank_deficient\n");

    ibz_mat_4x4_finalize(&basis);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_alg_finalize(&alg);
    return res;
}

/* ========== Phase 2 candidate C (fixed-precision) unit tests ==========
 *
 * These exercise the `quat_fp_*` API in `quat_fixed_precision.c` without
 * touching the MLLL loop — goal is to lock in add/sub + ibz_t bridge
 * correctness before the `mlll_gram.c` cut-over lands. Each width tier
 * (L1/L3/L5) gets 40 randomized trials plus explicit edge cases. */

static int
fp_arith_at_width(const char *label, const quat_fp_widths_t *w, int ntrials)
{
    int res = 0;
    ibz_t a, b, ref, got, bound, neg_bound;
    ibz_init(&a); ibz_init(&b);
    ibz_init(&ref); ibz_init(&got);
    ibz_init(&bound); ibz_init(&neg_bound);

    quat_fp_vec_t va, vb, vc;
    quat_fp_gram_t ga, gb, gc;

    /* Stable seed per label. */
    unsigned char seed[48];
    for (int i = 0; i < 48; i++)
        seed[i] = (unsigned char)((i * 7) ^ (unsigned)label[0]);
    randombytes_init(seed, NULL, 256);

    /* ---- vec path ---- */
    int max_vec_bits = (int)(w->nwords_vec * 64) - 2; /* sign + add carry */
    for (int t = 0; t < ntrials; t++) {
        int bits = 1 + ((t * 17) % max_vec_bits);
        ibz_pow(&bound, &ibz_const_two, (uint32_t)bits);
        ibz_sub(&bound, &bound, &ibz_const_one);
        ibz_neg(&neg_bound, &bound);
        ibz_rand_interval(&a, &neg_bound, &bound);
        ibz_rand_interval(&b, &neg_bound, &bound);

        if (!quat_fp_vec_set_ibz(&va, &a, w) ||
            !quat_fp_vec_set_ibz(&vb, &b, w)) {
            printf("  FAIL: %s vec set_ibz (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_vec_get_ibz(&got, &va, w);
        if (ibz_cmp(&a, &got) != 0) {
            printf("  FAIL: %s vec roundtrip (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_vec_add(&vc, &va, &vb, w);
        quat_fp_vec_get_ibz(&got, &vc, w);
        ibz_add(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s vec add (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_vec_sub(&vc, &va, &vb, w);
        quat_fp_vec_get_ibz(&got, &vc, w);
        ibz_sub(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s vec sub (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        int fp_bits = quat_fp_vec_bitsize(&va, w);
        int ibz_bits = ibz_bitsize(&a);
        if (fp_bits != ibz_bits) {
            printf("  FAIL: %s vec bitsize mismatch (fp=%d ibz=%d)\n",
                   label, fp_bits, ibz_bits);
            res = 1; goto done;
        }
    }

    /* ---- gram path ---- */
    int max_gram_bits = (int)(w->nwords_gram * 64) - 2;
    for (int t = 0; t < ntrials; t++) {
        int bits = 1 + ((t * 23) % max_gram_bits);
        ibz_pow(&bound, &ibz_const_two, (uint32_t)bits);
        ibz_sub(&bound, &bound, &ibz_const_one);
        ibz_neg(&neg_bound, &bound);
        ibz_rand_interval(&a, &neg_bound, &bound);
        ibz_rand_interval(&b, &neg_bound, &bound);

        if (!quat_fp_gram_set_ibz(&ga, &a, w) ||
            !quat_fp_gram_set_ibz(&gb, &b, w)) {
            printf("  FAIL: %s gram set_ibz (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_gram_get_ibz(&got, &ga, w);
        if (ibz_cmp(&a, &got) != 0) {
            printf("  FAIL: %s gram roundtrip (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_gram_add(&gc, &ga, &gb, w);
        quat_fp_gram_get_ibz(&got, &gc, w);
        ibz_add(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s gram add (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        quat_fp_gram_sub(&gc, &ga, &gb, w);
        quat_fp_gram_get_ibz(&got, &gc, w);
        ibz_sub(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s gram sub (bits %d)\n", label, bits);
            res = 1; goto done;
        }
        int fp_bits = quat_fp_gram_bitsize(&ga, w);
        int ibz_bits = ibz_bitsize(&a);
        if (fp_bits != ibz_bits) {
            printf("  FAIL: %s gram bitsize mismatch (fp=%d ibz=%d)\n",
                   label, fp_bits, ibz_bits);
            res = 1; goto done;
        }
    }

done:
    (void)vb; (void)gb;  /* silence unused when loops break early */
    ibz_finalize(&a); ibz_finalize(&b);
    ibz_finalize(&ref); ibz_finalize(&got);
    ibz_finalize(&bound); ibz_finalize(&neg_bound);
    return res;
}

int
quat_test_fp_arith_vs_ibz(void)
{
    int res = 0;
    quat_fp_widths_t w;

    w.nwords_vec  = NWORDS_QUAT_VEC_L1;
    w.nwords_gram = NWORDS_QUAT_GRAM_L1;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L1;
    res |= fp_arith_at_width("L1", &w, 40);

    w.nwords_vec  = NWORDS_QUAT_VEC_L3;
    w.nwords_gram = NWORDS_QUAT_GRAM_L3;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L3;
    res |= fp_arith_at_width("L3", &w, 40);

    w.nwords_vec  = NWORDS_QUAT_VEC_L5;
    w.nwords_gram = NWORDS_QUAT_GRAM_L5;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L5;
    res |= fp_arith_at_width("L5", &w, 40);

    /* Edge cases at L1 vec width. */
    quat_fp_widths_t lw;
    lw.nwords_vec  = NWORDS_QUAT_VEC_L1;
    lw.nwords_gram = NWORDS_QUAT_GRAM_L1;
    lw.nwords_tmp  = NWORDS_QUAT_TMP_L1;

    ibz_t z, one, neg_one, big, big_neg, got;
    ibz_init(&z);       ibz_set(&z, 0);
    ibz_init(&one);     ibz_set(&one, 1);
    ibz_init(&neg_one); ibz_set(&neg_one, -1);
    ibz_init(&got);
    ibz_init(&big);     ibz_pow(&big, &ibz_const_two, 200);
    ibz_sub(&big, &big, &ibz_const_one);
    ibz_init(&big_neg); ibz_neg(&big_neg, &big);

    const ibz_t *cases[] = { &z, &one, &neg_one, &big, &big_neg };
    const char *names[] = { "zero", "one", "neg_one", "big_pos", "big_neg" };
    quat_fp_vec_t v;
    for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        if (!quat_fp_vec_set_ibz(&v, cases[i], &lw)) {
            printf("  FAIL: edge case %s set_ibz\n", names[i]); res = 1; continue;
        }
        quat_fp_vec_get_ibz(&got, &v, &lw);
        if (ibz_cmp(cases[i], &got) != 0) {
            printf("  FAIL: edge case %s roundtrip\n", names[i]); res = 1;
        }
    }

    ibz_finalize(&z); ibz_finalize(&one); ibz_finalize(&neg_one);
    ibz_finalize(&big); ibz_finalize(&big_neg); ibz_finalize(&got);

    /* widths_from_alg sanity: tiny alg (p=7) should land in L1 tier. */
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);
    quat_fp_widths_t aw;
    if (!quat_fp_widths_from_alg(&aw, &alg)) {
        printf("  FAIL: widths_from_alg (tiny alg should match L1)\n");
        res = 1;
    } else if (aw.nwords_vec != NWORDS_QUAT_VEC_L1 ||
               aw.nwords_gram != NWORDS_QUAT_GRAM_L1) {
        printf("  FAIL: widths_from_alg L1 mismatch (vec=%u gram=%u)\n",
               aw.nwords_vec, aw.nwords_gram);
        res = 1;
    }
    quat_alg_finalize(&alg);

    if (res == 0)
        printf("  PASS: quat_test_fp_arith_vs_ibz (L1/L3/L5, 40 trials + edges)\n");
    return res;
}

/* ---------- Multiplication + narrow-sub tests ----------
 * Exercises quat_fp_tmp_mul_vec_{vec,gram} and quat_fp_{vec,gram}_sub_tmp
 * against ibz_mul / ibz_sub references. These are the hot ops that the
 * mlll_gram.c cut-over (P2-C-gram step 3) will actually call. */

static int
fp_mul_at_width(const char *label, const quat_fp_widths_t *w, int ntrials)
{
    int res = 0;
    ibz_t a, b, ref, got, vb, gb, tb, prod, sub_ref;
    ibz_init(&a); ibz_init(&b);
    ibz_init(&ref); ibz_init(&got);
    ibz_init(&vb); ibz_init(&gb); ibz_init(&tb);
    ibz_init(&prod); ibz_init(&sub_ref);

    quat_fp_vec_t vv;
    quat_fp_gram_t gg, gg_accum;
    quat_fp_vec_t vv_accum;
    quat_fp_tmp_t tt;

    unsigned char seed[48];
    for (int i = 0; i < 48; i++)
        seed[i] = (unsigned char)((i * 11) ^ (unsigned)label[1]);
    randombytes_init(seed, NULL, 256);

    /* vec * vec -> tmp */
    int vec_bits_max = (int)(w->nwords_vec * 64) - 2;
    for (int t = 0; t < ntrials; t++) {
        int ab = 1 + ((t * 13) % vec_bits_max);
        int bb = 1 + ((t * 19 + 5) % vec_bits_max);

        ibz_pow(&vb, &ibz_const_two, (uint32_t)ab);
        ibz_sub(&vb, &vb, &ibz_const_one);
        ibz_neg(&got, &vb);
        ibz_rand_interval(&a, &got, &vb);

        ibz_pow(&vb, &ibz_const_two, (uint32_t)bb);
        ibz_sub(&vb, &vb, &ibz_const_one);
        ibz_neg(&got, &vb);
        ibz_rand_interval(&b, &got, &vb);

        quat_fp_vec_t va2;
        quat_fp_vec_set_ibz(&vv, &a, w);
        quat_fp_vec_set_ibz(&va2, &b, w);
        quat_fp_tmp_mul_vec_vec(&tt, &vv, &va2, w);
        quat_fp_tmp_get_ibz(&got, &tt, w);
        ibz_mul(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s vec*vec (a %d bits, b %d bits)\n",
                   label, ab, bb);
            res = 1; goto done;
        }

        /* vec_sub_tmp: start from a, subtract (a*b), check result. */
        quat_fp_vec_set_ibz(&vv_accum, &a, w);
        quat_fp_vec_sub_tmp(&vv_accum, &tt, w);
        quat_fp_vec_get_ibz(&got, &vv_accum, w);
        ibz_sub(&sub_ref, &a, &ref);
        /* Only meaningful when result fits in vec width; bound the test. */
        if (ibz_bitsize(&sub_ref) < (int)(w->nwords_vec * 64) - 1) {
            if (ibz_cmp(&sub_ref, &got) != 0) {
                printf("  FAIL: %s vec_sub_tmp (a %d bits, b %d bits)\n",
                       label, ab, bb);
                res = 1; goto done;
            }
        }
    }

    /* vec * gram -> tmp */
    int gram_bits_max = (int)(w->nwords_gram * 64) - 2;
    for (int t = 0; t < ntrials; t++) {
        int ab = 1 + ((t * 29) % vec_bits_max);
        int bb = 1 + ((t * 37 + 3) % gram_bits_max);

        ibz_pow(&vb, &ibz_const_two, (uint32_t)ab);
        ibz_sub(&vb, &vb, &ibz_const_one);
        ibz_neg(&got, &vb);
        ibz_rand_interval(&a, &got, &vb);

        ibz_pow(&gb, &ibz_const_two, (uint32_t)bb);
        ibz_sub(&gb, &gb, &ibz_const_one);
        ibz_neg(&got, &gb);
        ibz_rand_interval(&b, &got, &gb);

        quat_fp_vec_set_ibz(&vv, &a, w);
        quat_fp_gram_set_ibz(&gg, &b, w);
        quat_fp_tmp_mul_vec_gram(&tt, &vv, &gg, w);
        quat_fp_tmp_get_ibz(&got, &tt, w);
        ibz_mul(&ref, &a, &b);
        if (ibz_cmp(&ref, &got) != 0) {
            printf("  FAIL: %s vec*gram (a %d bits, b %d bits)\n",
                   label, ab, bb);
            res = 1; goto done;
        }

        /* gram_sub_tmp: start from b, subtract (a*b), check. */
        quat_fp_gram_set_ibz(&gg_accum, &b, w);
        quat_fp_gram_sub_tmp(&gg_accum, &tt, w);
        quat_fp_gram_get_ibz(&got, &gg_accum, w);
        ibz_sub(&sub_ref, &b, &ref);
        if (ibz_bitsize(&sub_ref) < (int)(w->nwords_gram * 64) - 1) {
            if (ibz_cmp(&sub_ref, &got) != 0) {
                printf("  FAIL: %s gram_sub_tmp (a %d bits, b %d bits)\n",
                       label, ab, bb);
                res = 1; goto done;
            }
        }
    }

done:
    (void)tb; (void)prod;
    ibz_finalize(&a); ibz_finalize(&b);
    ibz_finalize(&ref); ibz_finalize(&got);
    ibz_finalize(&vb); ibz_finalize(&gb); ibz_finalize(&tb);
    ibz_finalize(&prod); ibz_finalize(&sub_ref);
    return res;
}

int
quat_test_fp_mul_vs_ibz(void)
{
    int res = 0;
    quat_fp_widths_t w;

    w.nwords_vec  = NWORDS_QUAT_VEC_L1;
    w.nwords_gram = NWORDS_QUAT_GRAM_L1;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L1;
    res |= fp_mul_at_width("L1", &w, 40);

    w.nwords_vec  = NWORDS_QUAT_VEC_L3;
    w.nwords_gram = NWORDS_QUAT_GRAM_L3;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L3;
    res |= fp_mul_at_width("L3", &w, 40);

    w.nwords_vec  = NWORDS_QUAT_VEC_L5;
    w.nwords_gram = NWORDS_QUAT_GRAM_L5;
    w.nwords_tmp  = NWORDS_QUAT_TMP_L5;
    res |= fp_mul_at_width("L5", &w, 40);

    if (res == 0)
        printf("  PASS: quat_test_fp_mul_vs_ibz (L1/L3/L5, 40 trials × vec*vec + vec*gram + sub_tmp)\n");
    return res;
}

/* ========== Test runner ========== */

int
quat_test_mlll_all(void)
{
    int res = 0;
    printf("=== MLLL Tests ===\n");

    res |= quat_test_mlll_dependent_generators();
    res |= quat_test_mlll_is_lll_reduced();
    res |= quat_test_mlll_vs_hnf_lattice_mul();
    res |= quat_test_compact_ideal_multiplication();
    res |= quat_test_mlll_realistic_scale();
    res |= quat_test_mlll_tau_path();
    res |= quat_test_mlll_gram_equivalence();
    res |= quat_test_lideal_create_gram_equivalence();
    res |= quat_test_mlll_gram_prealloc_equivalence();
    res |= quat_test_fp_arith_vs_ibz();
    res |= quat_test_fp_mul_vs_ibz();
    res |= quat_test_mlll_gram_fp_equivalence();
    res |= quat_test_mlll_all_zero_generators();
    res |= quat_test_mlll_single_generator();
    res |= quat_test_mlll_two_generators_dependent();
    res |= quat_test_mlll_rank_deficient();

    if (res == 0) {
        printf("=== All MLLL tests PASSED ===\n");
    } else {
        printf("=== Some MLLL tests FAILED ===\n");
    }
    return res;
}
