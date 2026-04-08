/**
 * @file mlll_benchmark.c
 * @brief Benchmark: HNF vs MLLL for ideal multiplication
 *
 * Compares:
 *   1. Maximum intermediate bit size during computation (fixed-precision requirement)
 *   2. Maximum bit size of output basis elements
 *   3. Execution time (wall clock)
 *
 * Uses the same random ideal lattice generation as lll_benchmarks.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <inttypes.h>

#include <quaternion.h>
#include <internal.h>
#include "lll_internals.h"
#include "mlll_internals.h"
#include "bitsize_tracker.h"
#include "quaternion_tests.h"
#include <rng.h>

/* Global tracker state — provided by weak symbols in bitsize_tracker.h */

/* ---------- helpers ---------- */

static int
lattice_max_bitsize(const quat_lattice_t *lat)
{
    int maxbits = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            int b = ibz_bitsize(&(lat->basis[i][j]));
            if (b > maxbits) maxbits = b;
        }
    return maxbits;
}

static int
generators_max_bitsize(const quat_lattice_t *lat1,
                       const quat_lattice_t *lat2)
{
    int maxbits = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            int b1 = ibz_bitsize(&(lat1->basis[i][j]));
            int b2 = ibz_bitsize(&(lat2->basis[i][j]));
            if (b1 > maxbits) maxbits = b1;
            if (b2 > maxbits) maxbits = b2;
        }
    return maxbits;
}

/**
 * @brief Run the HNF vs MLLL comparison for one pair of lattices
 */
static void
bench_one_mul(const quat_lattice_t *lat1,
              const quat_lattice_t *lat2,
              const quat_alg_t *alg,
              int *hnf_out_bits,
              int *mlll_out_bits,
              int *hnf_inter_bits,
              int *mlll_inter_bits,
              int *input_bits,
              double *hnf_time_ms,
              double *mlll_time_ms)
{
    quat_lattice_t prod_hnf, prod_mlll;
    quat_lattice_init(&prod_hnf);
    quat_lattice_init(&prod_mlll);
    clock_t t0, t1;

    *input_bits = generators_max_bitsize(lat1, lat2);

    /* HNF with tracking */
    tracker_reset();
    t0 = clock();
    quat_lattice_mul(&prod_hnf, lat1, lat2, alg);
    t1 = clock();
    *hnf_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    *hnf_out_bits = lattice_max_bitsize(&prod_hnf);
    *hnf_inter_bits = tracker_get_max();
    tracker_disable();

    /* MLLL with tracking */
    tracker_reset();
    t0 = clock();
    quat_lattice_mul_mlll(&prod_mlll, lat1, lat2, alg);
    t1 = clock();
    *mlll_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    *mlll_out_bits = lattice_max_bitsize(&prod_mlll);
    *mlll_inter_bits = tracker_get_max();
    tracker_disable();

    quat_lattice_finalize(&prod_hnf);
    quat_lattice_finalize(&prod_mlll);
}

/* ---------- main ---------- */

int
main(int argc, char *argv[])
{
    int iterations = 20;
    int norm_bitsize = 127;
    int level = 1;

    for (int i = 1; i < argc; i++) {
        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) continue;
        if (sscanf(argv[i], "--level=%d", &level) == 1) continue;
        if (strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s [--level=1|3|5] [--iterations=N]\n", argv[0]);
            return 0;
        }
    }

    ibz_t prime;
    ibz_init(&prime);

    if (level == 5) {
        ibz_set_from_str(&prime,
            "1afffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
        norm_bitsize = 254;
    } else if (level == 3) {
        ibz_set_from_str(&prime,
            "40fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
            "fffffffffffffffffffffff", 16);
        norm_bitsize = 193;
    } else {
        ibz_set_from_str(&prime,
            "4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
        norm_bitsize = 127;
    }

    printf("=== MLLL vs HNF Intermediate Bit Size Benchmark ===\n");
    printf("Level: %d, Norm bitsize: %d, Iterations: %d\n\n", level, norm_bitsize, iterations);

    quat_alg_t alg;
    quat_p_extremal_maximal_order_t order;
    quat_represent_integer_params_t params;
    quat_alg_init_set(&alg, &prime);
    quat_lattice_init(&(order.order));
    quat_alg_elem_init(&(order.t));
    quat_alg_elem_init(&(order.z));
    quat_lattice_O0_set_extremal(&order);
    params.algebra = &alg;
    params.order = &order;
    params.primality_test_iterations = 30;

    uint32_t seed[12] = { 0x12345678 };
    randombytes_init((unsigned char *)seed, NULL, 256);

    quat_lattice_t *lats1 = malloc(iterations * sizeof(quat_lattice_t));
    quat_lattice_t *lats2 = malloc(iterations * sizeof(quat_lattice_t));
    ibz_t *norms1 = malloc(iterations * sizeof(ibz_t));
    ibz_t *norms2 = malloc(iterations * sizeof(ibz_t));

    for (int i = 0; i < iterations; i++) {
        quat_lattice_init(&lats1[i]);
        quat_lattice_init(&lats2[i]);
        ibz_init(&norms1[i]);
        ibz_init(&norms2[i]);
    }

    printf("Generating %d random ideal lattice pairs (bitsize=%d)...\n", iterations, norm_bitsize);
    int ret1 = quat_test_input_random_ideal_lattice_generation(lats1, norms1, norm_bitsize, iterations, &params);
    int ret2 = quat_test_input_random_ideal_lattice_generation(lats2, norms2, norm_bitsize, iterations, &params);
    if (ret1 || ret2) {
        printf("ERROR: random lattice generation failed\n");
        return 1;
    }
    printf("Generation done.\n\n");

    /* Column headers */
    printf("%-5s | %5s | %12s %12s | %12s %12s | %8s %8s\n",
           "Trial", "Input", "HNF inter", "MLLL inter", "HNF out", "MLLL out", "HNF ms", "MLLL ms");
    printf("------+-------+---------------------------+---------------------------+-------------------\n");

    long sum_hnf_inter = 0, sum_mlll_inter = 0;
    long sum_hnf_out = 0, sum_mlll_out = 0;
    int max_hnf_inter = 0, max_mlll_inter = 0;
    int max_hnf_out = 0, max_mlll_out = 0;
    double total_hnf_ms = 0, total_mlll_ms = 0;

    for (int i = 0; i < iterations; i++) {
        int hnf_out, mlll_out, hnf_inter, mlll_inter, input_bits;
        double hnf_ms, mlll_ms;

        bench_one_mul(&lats1[i], &lats2[i], &alg,
                      &hnf_out, &mlll_out,
                      &hnf_inter, &mlll_inter,
                      &input_bits,
                      &hnf_ms, &mlll_ms);

        printf("%-5d | %5d | %12d %12d | %12d %12d | %8.2f %8.2f\n",
               i, input_bits, hnf_inter, mlll_inter, hnf_out, mlll_out, hnf_ms, mlll_ms);

        sum_hnf_inter += hnf_inter;
        sum_mlll_inter += mlll_inter;
        sum_hnf_out += hnf_out;
        sum_mlll_out += mlll_out;
        if (hnf_inter > max_hnf_inter) max_hnf_inter = hnf_inter;
        if (mlll_inter > max_mlll_inter) max_mlll_inter = mlll_inter;
        if (hnf_out > max_hnf_out) max_hnf_out = hnf_out;
        if (mlll_out > max_mlll_out) max_mlll_out = mlll_out;
        total_hnf_ms += hnf_ms;
        total_mlll_ms += mlll_ms;
    }

    printf("------+-------+---------------------------+---------------------------+-------------------\n");
    printf("\n=== Summary (Intermediate = max bit size during computation) ===\n");
    printf("Avg intermediate bits:  HNF=%ld  MLLL=%ld  (MLLL/HNF=%.3f)\n",
           sum_hnf_inter / iterations, sum_mlll_inter / iterations,
           sum_hnf_inter > 0 ? (double)sum_mlll_inter / sum_hnf_inter : 0.0);
    printf("Max intermediate bits:  HNF=%d  MLLL=%d  (MLLL/HNF=%.3f)\n",
           max_hnf_inter, max_mlll_inter,
           max_hnf_inter > 0 ? (double)max_mlll_inter / max_hnf_inter : 0.0);
    printf("Avg output bits:        HNF=%ld  MLLL=%ld\n",
           sum_hnf_out / iterations, sum_mlll_out / iterations);
    printf("Max output bits:        HNF=%d  MLLL=%d\n", max_hnf_out, max_mlll_out);
    printf("Total time:             HNF=%.2fms  MLLL=%.2fms  (MLLL/HNF=%.1f)\n",
           total_hnf_ms, total_mlll_ms,
           total_hnf_ms > 0 ? total_mlll_ms / total_hnf_ms : 0.0);

    /* Cleanup */
    for (int i = 0; i < iterations; i++) {
        quat_lattice_finalize(&lats1[i]);
        quat_lattice_finalize(&lats2[i]);
        ibz_finalize(&norms1[i]);
        ibz_finalize(&norms2[i]);
    }
    free(lats1); free(lats2);
    free(norms1); free(norms2);
    quat_lattice_finalize(&(order.order));
    quat_alg_elem_finalize(&(order.t));
    quat_alg_elem_finalize(&(order.z));
    quat_alg_finalize(&alg);
    ibz_finalize(&prime);
    return 0;
}
