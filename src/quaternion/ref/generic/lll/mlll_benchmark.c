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

#define SQISIGN_MLLL_GRAM_IMPL
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
typedef struct {
    int input_bits;
    int hnf_out_bits, mlll_out_bits, gram_out_bits;
    int hnf_inter_bits;
    int mlll_inter_bits, mlll_vec_bits, mlll_gso_bits;
    int gram_inter_bits, gram_vec_bits, gram_gso_bits;
    double hnf_time_ms, mlll_time_ms, gram_time_ms;
} bench_result_t;

static void
bench_one_mul(const quat_lattice_t *lat1,
              const quat_lattice_t *lat2,
              const quat_alg_t *alg,
              bench_result_t *out)
{
    quat_lattice_t prod_hnf, prod_mlll, prod_gram;
    quat_lattice_init(&prod_hnf);
    quat_lattice_init(&prod_mlll);
    quat_lattice_init(&prod_gram);
    clock_t t0, t1;

    out->input_bits = generators_max_bitsize(lat1, lat2);

    /* HNF with tracking */
    tracker_reset();
    t0 = clock();
    quat_lattice_mul(&prod_hnf, lat1, lat2, alg);
    t1 = clock();
    out->hnf_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->hnf_out_bits = lattice_max_bitsize(&prod_hnf);
    out->hnf_inter_bits = tracker_get_max();
    tracker_disable();

    /* MLLL (Cohen integral GSO) with tracking */
    tracker_reset();
    t0 = clock();
    quat_lattice_mul_mlll(&prod_mlll, lat1, lat2, alg);
    t1 = clock();
    out->mlll_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->mlll_out_bits = lattice_max_bitsize(&prod_mlll);
    out->mlll_inter_bits = tracker_get_max();
    out->mlll_vec_bits = tracker_get_vec_max();
    out->mlll_gso_bits = tracker_get_gso_max();
    tracker_disable();

    /* MLLL Gram (L2-style) with tracking */
    tracker_reset();
    t0 = clock();
    quat_lattice_mul_mlll_gram(&prod_gram, lat1, lat2, alg);
    t1 = clock();
    out->gram_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->gram_out_bits = lattice_max_bitsize(&prod_gram);
    out->gram_inter_bits = tracker_get_max();
    out->gram_vec_bits = tracker_get_vec_max();
    out->gram_gso_bits = tracker_get_gso_max();
    tracker_disable();

    quat_lattice_finalize(&prod_hnf);
    quat_lattice_finalize(&prod_mlll);
    quat_lattice_finalize(&prod_gram);
}

/* ---------- Alg 3 bench: quat_lideal_create HNF vs Gram ---------- */

typedef struct {
    int input_bits;
    int hnf_out_bits, gram_out_bits;
    int hnf_inter_bits;
    int gram_inter_bits, gram_vec_bits, gram_gso_bits;
    double hnf_time_ms, gram_time_ms;
} alg3_result_t;

static void
bench_one_alg3(const quat_lattice_t *src_lat,
               const ibz_t *norm,
               const quat_lattice_t *order,
               const quat_alg_t *alg,
               alg3_result_t *out)
{
    quat_alg_elem_t x;
    quat_alg_elem_init(&x);
    for (int j = 0; j < 4; j++)
        ibz_copy(&x.coord[j], &src_lat->basis[j][0]);
    ibz_copy(&x.denom, &src_lat->denom);

    out->input_bits = 0;
    for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&x.coord[j]);
        if (b > out->input_bits) out->input_bits = b;
    }
    int nb = ibz_bitsize(norm);
    if (nb > out->input_bits) out->input_bits = nb;

    quat_left_ideal_t I_hnf, I_gram;
    quat_left_ideal_init(&I_hnf);
    quat_left_ideal_init(&I_gram);

    clock_t t0, t1;

    /* HNF path */
    tracker_reset();
    t0 = clock();
    quat_lideal_create(&I_hnf, &x, norm, order, alg);
    t1 = clock();
    out->hnf_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->hnf_out_bits = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&I_hnf.lattice.basis[i][j]);
        if (b > out->hnf_out_bits) out->hnf_out_bits = b;
    }
    out->hnf_inter_bits = tracker_get_max();
    tracker_disable();

    /* Gram path */
    tracker_reset();
    t0 = clock();
    quat_lideal_create_mlll_gram(&I_gram, &x, norm, order, alg);
    t1 = clock();
    out->gram_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->gram_out_bits = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&I_gram.lattice.basis[i][j]);
        if (b > out->gram_out_bits) out->gram_out_bits = b;
    }
    out->gram_inter_bits = tracker_get_max();
    out->gram_vec_bits = tracker_get_vec_max();
    out->gram_gso_bits = tracker_get_gso_max();
    tracker_disable();

    quat_left_ideal_finalize(&I_hnf);
    quat_left_ideal_finalize(&I_gram);
    quat_alg_elem_finalize(&x);
}

/* ---------- Alg 4 bench: quat_lideal_prime_norm_reduced_equivalent HNF vs Gram ----------
 *
 * Note: per paper Sampling §04 / PLAN_KO §P3.2-b, Alg 4 itself does not call
 * MLLL — it uses a reduced basis (LLL). The HNF↔MLLL_GRAM split here lives
 * inside `quat_lideal_reduce_basis(_mlll_gram)` which is invoked at the very
 * top of each variant. So this bench measures the cost/quality difference of
 * basis reduction as it propagates through Alg 4's prime-norm search loop.
 */

typedef struct {
    int input_bits;
    int hnf_out_bits, gram_out_bits;
    int hnf_norm_bits, gram_norm_bits;
    int hnf_inter_bits;
    int gram_inter_bits, gram_vec_bits, gram_gso_bits;
    double hnf_time_ms, gram_time_ms;
    int hnf_found, gram_found;
} alg4_result_t;

static void
bench_one_alg4(const quat_lattice_t *src_lat,
               const ibz_t *norm,
               const quat_lattice_t *order,
               const quat_alg_t *alg,
               alg4_result_t *out)
{
    quat_alg_elem_t x;
    quat_alg_elem_init(&x);
    for (int j = 0; j < 4; j++)
        ibz_copy(&x.coord[j], &src_lat->basis[j][0]);
    ibz_copy(&x.denom, &src_lat->denom);

    quat_left_ideal_t I_hnf, I_gram;
    quat_left_ideal_init(&I_hnf);
    quat_left_ideal_init(&I_gram);

    /* Build identical input ideals for both paths (Alg 3 = quat_lideal_create). */
    quat_lideal_create(&I_hnf, &x, norm, order, alg);
    quat_lideal_create(&I_gram, &x, norm, order, alg);

    out->input_bits = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&I_hnf.lattice.basis[i][j]);
        if (b > out->input_bits) out->input_bits = b;
    }

    clock_t t0, t1;

    /* HNF path: in-place, returns 1 on success. */
    tracker_reset();
    t0 = clock();
    out->hnf_found = quat_lideal_prime_norm_reduced_equivalent(&I_hnf, alg, 30, 20);
    t1 = clock();
    out->hnf_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->hnf_out_bits = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&I_hnf.lattice.basis[i][j]);
        if (b > out->hnf_out_bits) out->hnf_out_bits = b;
    }
    out->hnf_norm_bits = ibz_bitsize(&I_hnf.norm);
    out->hnf_inter_bits = tracker_get_max();
    tracker_disable();

    /* GRAM path: same contract. */
    tracker_reset();
    t0 = clock();
    out->gram_found = quat_lideal_prime_norm_reduced_equivalent_mlll_gram(&I_gram, alg, 30, 20);
    t1 = clock();
    out->gram_time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    out->gram_out_bits = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int b = ibz_bitsize(&I_gram.lattice.basis[i][j]);
        if (b > out->gram_out_bits) out->gram_out_bits = b;
    }
    out->gram_norm_bits = ibz_bitsize(&I_gram.norm);
    out->gram_inter_bits = tracker_get_max();
    out->gram_vec_bits = tracker_get_vec_max();
    out->gram_gso_bits = tracker_get_gso_max();
    tracker_disable();

    quat_left_ideal_finalize(&I_hnf);
    quat_left_ideal_finalize(&I_gram);
    quat_alg_elem_finalize(&x);
}

/* ---------- main ---------- */

int
main(int argc, char *argv[])
{
    int iterations = 20;
    int norm_bitsize = 127;
    int level = 1;
    const char *mode = "alg2";
    int prealloc = 0;
    int fp = 0;
    int ret_code = 0;

    for (int i = 1; i < argc; i++) {
        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) continue;
        if (sscanf(argv[i], "--level=%d", &level) == 1) continue;
        if (strncmp(argv[i], "--mode=", 7) == 0) { mode = argv[i] + 7; continue; }
        if (strcmp(argv[i], "--prealloc") == 0) { prealloc = 1; continue; }
        if (strcmp(argv[i], "--fp") == 0) { fp = 1; continue; }
        if (strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s [--level=1|3|5] [--iterations=N] "
                   "[--mode=alg2|alg3|alg4] [--prealloc] [--fp]\n", argv[0]);
            printf("  --mode=alg2: CompactIdealMultiplication (lattice_mul HNF vs MLLL vs MLLL_GRAM).\n");
            printf("  --mode=alg3: RandomIdealGivenPrimeNorm (lideal_create HNF vs MLLL_GRAM).\n");
            printf("  --mode=alg4: RandomEquivalentPrimeIdeal (prime_norm_reduced_equivalent HNF vs MLLL_GRAM).\n");
            printf("  --prealloc: enable Phase 2 candidate B mpz_realloc2 hints on the GRAM path.\n");
            printf("  --fp:       enable Phase 2 candidate C fixed-precision body on the GRAM path.\n");
            printf("  (--prealloc and --fp are independent; fp takes precedence when both set.)\n");
            return 0;
        }
    }

    if (prealloc && fp) {
        fprintf(stderr, "WARN: both --prealloc and --fp set; fp path wins "
                "(prealloc hints do not apply to the fp body).\n");
    }

    if (strcmp(mode, "alg2") != 0 && strcmp(mode, "alg3") != 0 && strcmp(mode, "alg4") != 0) {
        fprintf(stderr, "ERROR: unknown --mode=%s (expected alg2, alg3, or alg4)\n", mode);
        return 1;
    }
    if (iterations <= 0) {
        fprintf(stderr, "ERROR: --iterations must be positive (got %d)\n", iterations);
        return 1;
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

    quat_mlll_gram_set_prealloc_mode(prealloc);
    quat_mlll_gram_set_fp_mode(fp);

    printf("=== MLLL vs HNF Intermediate Bit Size Benchmark ===\n");
    printf("Level: %d, Norm bitsize: %d, Iterations: %d\n", level, norm_bitsize, iterations);
    printf("GRAM prealloc (Phase 2 candidate B): %s\n",
           prealloc ? "ON" : "OFF");
    printf("GRAM fp       (Phase 2 candidate C): %s\n\n",
           fp ? "ON" : "OFF");

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
    if (!lats1 || !lats2 || !norms1 || !norms2) {
        fprintf(stderr, "ERROR: malloc failed for iterations=%d\n", iterations);
        free(lats1); free(lats2); free(norms1); free(norms2);
        quat_lattice_finalize(&(order.order));
        quat_alg_elem_finalize(&(order.t));
        quat_alg_elem_finalize(&(order.z));
        quat_alg_finalize(&alg);
        ibz_finalize(&prime);
        return 1;
    }

    for (int i = 0; i < iterations; i++) {
        quat_lattice_init(&lats1[i]);
        quat_lattice_init(&lats2[i]);
        ibz_init(&norms1[i]);
        ibz_init(&norms2[i]);
    }

    printf("Generating %d random ideal lattice pairs (bitsize=%d)...\n", iterations, norm_bitsize);
    int is_alg3 = (strcmp(mode, "alg3") == 0);
    int is_alg4 = (strcmp(mode, "alg4") == 0);
    int needs_set2 = !is_alg3 && !is_alg4;
    printf("Mode: %s\n", mode);

    clock_t gen_t0 = clock();
    int ret1 = quat_test_input_random_ideal_lattice_generation(lats1, norms1, norm_bitsize, iterations, &params);
    clock_t gen_t1 = clock();
    printf("  Set 1: %.2f ms\n", (double)(gen_t1 - gen_t0) * 1000.0 / CLOCKS_PER_SEC);
    int ret2 = 0;
    if (needs_set2) {
        gen_t0 = clock();
        ret2 = quat_test_input_random_ideal_lattice_generation(lats2, norms2, norm_bitsize, iterations, &params);
        gen_t1 = clock();
        printf("  Set 2: %.2f ms\n", (double)(gen_t1 - gen_t0) * 1000.0 / CLOCKS_PER_SEC);
    }
    if (ret1 || ret2) {
        fprintf(stderr, "ERROR: random lattice generation failed\n");
        ret_code = 1;
        goto done;
    }
    printf("Generation done.\n\n");

    if (is_alg4) {
        /* ---- Alg 4 bench ---- */
        printf("%-4s | %5s | %7s | %7s %7s %7s | %7s %7s | %5s %5s | %4s %4s | %6s %6s\n",
               "T", "Inp", "HNFint",
               "GRMtot", "GRMv", "GRMg",
               "HNFout", "GRMout",
               "HNFn", "GRMn",
               "Hok", "Gok",
               "HNFms", "GRMms");
        printf("-----+-------+---------+-------------------------+-----------------+-----------+---------+---------------\n");

        long sum_hnf_inter = 0, sum_gram_inter = 0, sum_gram_vec = 0, sum_gram_gso = 0;
        long sum_hnf_out = 0, sum_gram_out = 0;
        long sum_hnf_norm = 0, sum_gram_norm = 0;
        int max_hnf_inter = 0, max_gram_inter = 0, max_gram_vec = 0, max_gram_gso = 0;
        int max_hnf_out = 0, max_gram_out = 0;
        int max_hnf_norm = 0, max_gram_norm = 0;
        int hnf_success = 0, gram_success = 0;
        double total_hnf_ms = 0, total_gram_ms = 0;

        for (int i = 0; i < iterations; i++) {
            alg4_result_t r;
            bench_one_alg4(&lats1[i], &norms1[i], &(order.order), &alg, &r);

            if (iterations <= 50) {
                printf("%-4d | %5d | %7d | %7d %7d %7d | %7d %7d | %5d %5d | %4d %4d | %6.1f %6.1f\n",
                       i, r.input_bits, r.hnf_inter_bits,
                       r.gram_inter_bits, r.gram_vec_bits, r.gram_gso_bits,
                       r.hnf_out_bits, r.gram_out_bits,
                       r.hnf_norm_bits, r.gram_norm_bits,
                       r.hnf_found, r.gram_found,
                       r.hnf_time_ms, r.gram_time_ms);
            }

            sum_hnf_inter += r.hnf_inter_bits;
            sum_gram_inter += r.gram_inter_bits;
            sum_gram_vec += r.gram_vec_bits;
            sum_gram_gso += r.gram_gso_bits;
            sum_hnf_out += r.hnf_out_bits;
            sum_gram_out += r.gram_out_bits;
            sum_hnf_norm += r.hnf_norm_bits;
            sum_gram_norm += r.gram_norm_bits;
            if (r.hnf_inter_bits > max_hnf_inter) max_hnf_inter = r.hnf_inter_bits;
            if (r.gram_inter_bits > max_gram_inter) max_gram_inter = r.gram_inter_bits;
            if (r.gram_vec_bits > max_gram_vec) max_gram_vec = r.gram_vec_bits;
            if (r.gram_gso_bits > max_gram_gso) max_gram_gso = r.gram_gso_bits;
            if (r.hnf_out_bits > max_hnf_out) max_hnf_out = r.hnf_out_bits;
            if (r.gram_out_bits > max_gram_out) max_gram_out = r.gram_out_bits;
            if (r.hnf_norm_bits > max_hnf_norm) max_hnf_norm = r.hnf_norm_bits;
            if (r.gram_norm_bits > max_gram_norm) max_gram_norm = r.gram_norm_bits;
            hnf_success += (r.hnf_found ? 1 : 0);
            gram_success += (r.gram_found ? 1 : 0);
            total_hnf_ms += r.hnf_time_ms;
            total_gram_ms += r.gram_time_ms;
        }

        printf("\n=== Alg 4 Summary ===\n");
        printf("HNF  intermediate:  avg=%ld  max=%d\n",
               sum_hnf_inter / iterations, max_hnf_inter);
        printf("GRAM total:         avg=%ld  max=%d\n",
               sum_gram_inter / iterations, max_gram_inter);
        printf("GRAM vec coords:    avg=%ld  max=%d  (Lemma 1)\n",
               sum_gram_vec / iterations, max_gram_vec);
        printf("GRAM Gram entries:  avg=%ld  max=%d  (Lemma 3)\n",
               sum_gram_gso / iterations, max_gram_gso);
        printf("Output basis bits:  HNF avg=%ld max=%d  GRAM avg=%ld max=%d\n",
               sum_hnf_out / iterations, max_hnf_out,
               sum_gram_out / iterations, max_gram_out);
        printf("Output norm bits:   HNF avg=%ld max=%d  GRAM avg=%ld max=%d\n",
               sum_hnf_norm / iterations, max_hnf_norm,
               sum_gram_norm / iterations, max_gram_norm);
        printf("Success rate:       HNF %d/%d  GRAM %d/%d\n",
               hnf_success, iterations, gram_success, iterations);
        printf("Total time:         HNF=%.2fms  GRAM=%.2fms  (GRAM/HNF=%.2f)\n",
               total_hnf_ms, total_gram_ms,
               total_hnf_ms > 0 ? total_gram_ms / total_hnf_ms : 0.0);
        goto done;
    }

    if (is_alg3) {
        /* ---- Alg 3 bench ---- */
        printf("%-4s | %5s | %7s | %7s %7s %7s | %7s %7s | %6s %6s\n",
               "T", "Inp", "HNFint",
               "GRMtot", "GRMv", "GRMg",
               "HNFout", "GRMout",
               "HNFms", "GRMms");
        printf("-----+-------+---------+-------------------------+-----------------+---------------\n");

        long sum_hnf_inter = 0, sum_gram_inter = 0, sum_gram_vec = 0, sum_gram_gso = 0;
        long sum_hnf_out = 0, sum_gram_out = 0;
        int max_hnf_inter = 0, max_gram_inter = 0, max_gram_vec = 0, max_gram_gso = 0;
        int max_hnf_out = 0, max_gram_out = 0;
        double total_hnf_ms = 0, total_gram_ms = 0;

        for (int i = 0; i < iterations; i++) {
            alg3_result_t r;
            bench_one_alg3(&lats1[i], &norms1[i], &(order.order), &alg, &r);

            if (iterations <= 50) {
                printf("%-4d | %5d | %7d | %7d %7d %7d | %7d %7d | %6.1f %6.1f\n",
                       i, r.input_bits, r.hnf_inter_bits,
                       r.gram_inter_bits, r.gram_vec_bits, r.gram_gso_bits,
                       r.hnf_out_bits, r.gram_out_bits,
                       r.hnf_time_ms, r.gram_time_ms);
            }

            sum_hnf_inter += r.hnf_inter_bits;
            sum_gram_inter += r.gram_inter_bits;
            sum_gram_vec += r.gram_vec_bits;
            sum_gram_gso += r.gram_gso_bits;
            sum_hnf_out += r.hnf_out_bits;
            sum_gram_out += r.gram_out_bits;
            if (r.hnf_inter_bits > max_hnf_inter) max_hnf_inter = r.hnf_inter_bits;
            if (r.gram_inter_bits > max_gram_inter) max_gram_inter = r.gram_inter_bits;
            if (r.gram_vec_bits > max_gram_vec) max_gram_vec = r.gram_vec_bits;
            if (r.gram_gso_bits > max_gram_gso) max_gram_gso = r.gram_gso_bits;
            if (r.hnf_out_bits > max_hnf_out) max_hnf_out = r.hnf_out_bits;
            if (r.gram_out_bits > max_gram_out) max_gram_out = r.gram_out_bits;
            total_hnf_ms += r.hnf_time_ms;
            total_gram_ms += r.gram_time_ms;
        }

        printf("\n=== Alg 3 Summary ===\n");
        printf("HNF  intermediate:  avg=%ld  max=%d\n",
               sum_hnf_inter / iterations, max_hnf_inter);
        printf("GRAM total:         avg=%ld  max=%d\n",
               sum_gram_inter / iterations, max_gram_inter);
        printf("GRAM vec coords:    avg=%ld  max=%d  (Lemma 1)\n",
               sum_gram_vec / iterations, max_gram_vec);
        printf("GRAM Gram entries:  avg=%ld  max=%d  (Lemma 3)\n",
               sum_gram_gso / iterations, max_gram_gso);
        printf("Output bits:        HNF avg=%ld max=%d  GRAM avg=%ld max=%d\n",
               sum_hnf_out / iterations, max_hnf_out,
               sum_gram_out / iterations, max_gram_out);
        printf("Total time:         HNF=%.2fms  GRAM=%.2fms  (GRAM/HNF=%.2f)\n",
               total_hnf_ms, total_gram_ms,
               total_hnf_ms > 0 ? total_gram_ms / total_hnf_ms : 0.0);
        goto done;
    }

    /* Column headers */
    printf("%-4s | %5s | %7s | %7s %7s %7s | %7s %7s %7s | %6s %6s %6s\n",
           "T", "Inp", "HNFint",
           "MLLLtot", "MLLLv", "MLLLg",
           "GRMtot", "GRMv", "GRMg",
           "HNFms", "MLLLms", "GRMms");
    printf("-----+-------+---------+-------------------------+-------------------------+--------------------\n");

    long sum_hnf_inter = 0, sum_mlll_inter = 0, sum_gram_inter = 0;
    long sum_mlll_vec = 0, sum_mlll_gso = 0;
    long sum_gram_vec = 0, sum_gram_gso = 0;
    long sum_hnf_out = 0, sum_mlll_out = 0, sum_gram_out = 0;
    int max_hnf_inter = 0, max_mlll_inter = 0, max_gram_inter = 0;
    int max_mlll_vec = 0, max_mlll_gso = 0;
    int max_gram_vec = 0, max_gram_gso = 0;
    int max_hnf_out = 0, max_mlll_out = 0, max_gram_out = 0;
    double total_hnf_ms = 0, total_mlll_ms = 0, total_gram_ms = 0;

    for (int i = 0; i < iterations; i++) {
        bench_result_t r;
        bench_one_mul(&lats1[i], &lats2[i], &alg, &r);

        printf("%-4d | %5d | %7d | %7d %7d %7d | %7d %7d %7d | %6.1f %6.1f %6.1f\n",
               i, r.input_bits, r.hnf_inter_bits,
               r.mlll_inter_bits, r.mlll_vec_bits, r.mlll_gso_bits,
               r.gram_inter_bits, r.gram_vec_bits, r.gram_gso_bits,
               r.hnf_time_ms, r.mlll_time_ms, r.gram_time_ms);

        sum_hnf_inter += r.hnf_inter_bits;
        sum_mlll_inter += r.mlll_inter_bits;
        sum_mlll_vec += r.mlll_vec_bits;
        sum_mlll_gso += r.mlll_gso_bits;
        sum_gram_inter += r.gram_inter_bits;
        sum_gram_vec += r.gram_vec_bits;
        sum_gram_gso += r.gram_gso_bits;
        sum_hnf_out += r.hnf_out_bits;
        sum_mlll_out += r.mlll_out_bits;
        sum_gram_out += r.gram_out_bits;
        if (r.hnf_inter_bits > max_hnf_inter) max_hnf_inter = r.hnf_inter_bits;
        if (r.mlll_inter_bits > max_mlll_inter) max_mlll_inter = r.mlll_inter_bits;
        if (r.mlll_vec_bits > max_mlll_vec) max_mlll_vec = r.mlll_vec_bits;
        if (r.mlll_gso_bits > max_mlll_gso) max_mlll_gso = r.mlll_gso_bits;
        if (r.gram_inter_bits > max_gram_inter) max_gram_inter = r.gram_inter_bits;
        if (r.gram_vec_bits > max_gram_vec) max_gram_vec = r.gram_vec_bits;
        if (r.gram_gso_bits > max_gram_gso) max_gram_gso = r.gram_gso_bits;
        if (r.hnf_out_bits > max_hnf_out) max_hnf_out = r.hnf_out_bits;
        if (r.mlll_out_bits > max_mlll_out) max_mlll_out = r.mlll_out_bits;
        if (r.gram_out_bits > max_gram_out) max_gram_out = r.gram_out_bits;
        total_hnf_ms += r.hnf_time_ms;
        total_mlll_ms += r.mlll_time_ms;
        total_gram_ms += r.gram_time_ms;
    }

    printf("-----+-------+---------+-------------------------+-------------------------+--------------------\n");
    printf("\n=== Summary ===\n");
    printf("HNF  intermediate:  avg=%ld  max=%d\n",
           sum_hnf_inter / iterations, max_hnf_inter);
    printf("MLLL total:         avg=%ld  max=%d\n",
           sum_mlll_inter / iterations, max_mlll_inter);
    printf("MLLL vec coords:    avg=%ld  max=%d  (Lemma 1 bounded)\n",
           sum_mlll_vec / iterations, max_mlll_vec);
    printf("MLLL GSO coeffs:    avg=%ld  max=%d  (integral GSO overhead)\n",
           sum_mlll_gso / iterations, max_mlll_gso);
    printf("GRAM total:         avg=%ld  max=%d\n",
           sum_gram_inter / iterations, max_gram_inter);
    printf("GRAM vec coords:    avg=%ld  max=%d  (Lemma 1 bounded)\n",
           sum_gram_vec / iterations, max_gram_vec);
    printf("GRAM Gram entries:  avg=%ld  max=%d  (Lemma 3 bound target)\n",
           sum_gram_gso / iterations, max_gram_gso);
    printf("Output bits:        HNF avg=%ld max=%d  MLLL avg=%ld max=%d  GRAM avg=%ld max=%d\n",
           sum_hnf_out / iterations, max_hnf_out,
           sum_mlll_out / iterations, max_mlll_out,
           sum_gram_out / iterations, max_gram_out);
    printf("Total time:         HNF=%.2fms  MLLL=%.2fms  GRAM=%.2fms  (GRAM/MLLL=%.2f)\n",
           total_hnf_ms, total_mlll_ms, total_gram_ms,
           total_mlll_ms > 0 ? total_gram_ms / total_mlll_ms : 0.0);

done:
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
    return ret_code;
}
