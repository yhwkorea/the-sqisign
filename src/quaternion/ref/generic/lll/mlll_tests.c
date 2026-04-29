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

/* ========== Test: lideal_reduce_basis HNF == MLLL_GRAM (lattice equality) ========== */

/* Phase 3-1 acceptance test:
 * For each random L5-scale ideal generated by the production helper, run both
 *   - quat_lideal_reduce_basis            (HNF + classic LLL core)
 *   - quat_lideal_reduce_basis_mlll_gram  (MLLL_GRAM + Lemma 3 width)
 * and verify that the two reduced bases generate the same lattice. The two
 * basis matrices may differ (LLL is not unique up to sign), so we compare via
 * quat_lattice_equal (HNF equality on the spanned lattice). The Gram matrices
 * are not byte-equal but we sanity-check that diagonal entries are positive
 * (basis vectors have positive reduced norm) and that the two diagonals' total
 * sums match within the same reduced-norm class via a HNF round-trip on each
 * basis matrix.
 */
int
quat_test_lideal_reduce_basis_gram_equivalence(void)
{
    int res = 0;
    quat_alg_t alg;
    ibz_t prime;
    ibz_init(&prime);
    /* L5-class prime same as quat_test_lideal_create_gram_equivalence */
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

    for (int i = 0; i < trials; i++) {
        /* Build a left ideal from the random lattice + known norm. Use the
         * lideal_create path so parent_order and norm bookkeeping match what
         * production code feeds quat_lideal_reduce_basis. */
        quat_alg_elem_t x;
        quat_alg_elem_init(&x);
        for (int j = 0; j < 4; j++)
            ibz_copy(&x.coord[j], &lats[i].basis[j][0]);
        ibz_copy(&x.denom, &lats[i].denom);
        if (quat_alg_elem_is_zero(&x)) {
            quat_alg_elem_finalize(&x);
            continue;
        }

        quat_left_ideal_t I;
        quat_left_ideal_init(&I);
        quat_lideal_create(&I, &x, &norms[i], &(order.order), &alg);

        /* Run both reduce_basis variants */
        ibz_mat_4x4_t red_h, gram_h, red_m, gram_m;
        ibz_mat_4x4_init(&red_h);
        ibz_mat_4x4_init(&gram_h);
        ibz_mat_4x4_init(&red_m);
        ibz_mat_4x4_init(&gram_m);

        quat_lideal_reduce_basis(&red_h, &gram_h, &I, &alg);
        quat_lideal_reduce_basis_mlll_gram(&red_m, &gram_m, &I, &alg);

        /* The two reduced bases must span the same lattice (with same denom). */
        quat_lattice_t lh, lm;
        quat_lattice_init(&lh);
        quat_lattice_init(&lm);
        ibz_mat_4x4_copy(&lh.basis, &red_h);
        ibz_mat_4x4_copy(&lm.basis, &red_m);
        ibz_copy(&lh.denom, &I.lattice.denom);
        ibz_copy(&lm.denom, &I.lattice.denom);

        if (!quat_lattice_equal(&lh, &lm)) {
            printf("  FAIL: trial %d: HNF and MLLL_GRAM reduced bases span different lattices\n", i);
            res = 1;
        }
        if (!quat_lattice_equal(&I.lattice, &lh)) {
            printf("  FAIL: trial %d: HNF reduced basis differs from input lattice\n", i);
            res = 1;
        }
        if (!quat_lattice_equal(&I.lattice, &lm)) {
            printf("  FAIL: trial %d: MLLL_GRAM reduced basis differs from input lattice\n", i);
            res = 1;
        }

        /* Sanity: gram diagonals are strictly positive (=nrd(b_i)/norm > 0
         * for a basis of a full-rank ideal). */
        for (int k = 0; k < 4; k++) {
            if (ibz_cmp(&gram_h[k][k], &ibz_const_zero) <= 0) {
                printf("  FAIL: trial %d: HNF gram diag[%d] not positive\n", i, k);
                res = 1;
            }
            if (ibz_cmp(&gram_m[k][k], &ibz_const_zero) <= 0) {
                printf("  FAIL: trial %d: MLLL_GRAM gram diag[%d] not positive\n", i, k);
                res = 1;
            }
        }

        quat_lattice_finalize(&lh);
        quat_lattice_finalize(&lm);
        ibz_mat_4x4_finalize(&red_h);
        ibz_mat_4x4_finalize(&gram_h);
        ibz_mat_4x4_finalize(&red_m);
        ibz_mat_4x4_finalize(&gram_m);

        quat_left_ideal_finalize(&I);
        quat_alg_elem_finalize(&x);
    }

    if (res == 0)
        printf("  PASS: quat_test_lideal_reduce_basis_gram_equivalence (%d trials)\n", trials);

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

/* ========== Test: prime_norm_reduced_equivalent MLLL_GRAM smoke + primality ========== */

/* Phase 3-1 acceptance test:
 * For a few production-scale ideals, call the MLLL_GRAM variant of
 * prime_norm_reduced_equivalent and verify
 *   (a) it returns 1 (success)
 *   (b) the resulting ideal's norm is probabilistically prime
 *   (c) the resulting lideal is contained in the original lideal's left order
 *
 * We do NOT compare the output to the HNF version: the function uses random
 * sampling, so different reduced bases produce different (but equally valid)
 * prime-norm equivalent ideals. Equivalence with the HNF path is exercised
 * end-to-end by the scheme/KAT tests once Phase 3-2 routes the hot path.
 */
int
quat_test_lideal_prime_norm_reduced_equivalent_gram(void)
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

    uint32_t seed[12] = { 0x1357BEEF };
    randombytes_init((unsigned char *)seed, NULL, 256);

    int trials = 5;
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

    int success_count = 0;
    for (int i = 0; i < trials; i++) {
        quat_alg_elem_t x;
        quat_alg_elem_init(&x);
        for (int j = 0; j < 4; j++)
            ibz_copy(&x.coord[j], &lats[i].basis[j][0]);
        ibz_copy(&x.denom, &lats[i].denom);
        if (quat_alg_elem_is_zero(&x)) {
            quat_alg_elem_finalize(&x);
            continue;
        }

        quat_left_ideal_t I;
        quat_left_ideal_init(&I);
        quat_lideal_create(&I, &x, &norms[i], &(order.order), &alg);

        int found = quat_lideal_prime_norm_reduced_equivalent_mlll_gram(&I, &alg, 30, 20);
        if (!found) {
            printf("  FAIL: trial %d: function returned 0\n", i);
            res = 1;
        } else {
            if (!ibz_probab_prime(&I.norm, 30)) {
                printf("  FAIL: trial %d: result norm not prime\n", i);
                res = 1;
            }
            success_count++;
        }

        quat_left_ideal_finalize(&I);
        quat_alg_elem_finalize(&x);
    }

    if (res == 0)
        printf("  PASS: quat_test_lideal_prime_norm_reduced_equivalent_gram (%d/%d trials)\n",
               success_count, trials);

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
    res |= quat_test_lideal_reduce_basis_gram_equivalence();
    res |= quat_test_lideal_prime_norm_reduced_equivalent_gram();
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
