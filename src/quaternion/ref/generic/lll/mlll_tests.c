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
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"

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

    if (res == 0) {
        printf("=== All MLLL tests PASSED ===\n");
    } else {
        printf("=== Some MLLL tests FAILED ===\n");
    }
    return res;
}
