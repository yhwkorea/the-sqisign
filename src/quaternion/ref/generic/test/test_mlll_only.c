#include <stdio.h>
#include "quaternion_tests.h"
#include "mlll_internals.h"

/* Run individual tests to isolate which one hangs */
int main(void)
{
    printf("Starting MLLL test: dependent_generators...\n");
    fflush(stdout);
    int res = quat_test_mlll_dependent_generators();
    printf("dependent_generators done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: is_lll_reduced...\n");
    fflush(stdout);
    res |= quat_test_mlll_is_lll_reduced();
    printf("is_lll_reduced done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: vs_hnf_lattice_mul...\n");
    fflush(stdout);
    res |= quat_test_mlll_vs_hnf_lattice_mul();
    printf("vs_hnf_lattice_mul done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: compact_ideal_multiplication...\n");
    fflush(stdout);
    res |= quat_test_compact_ideal_multiplication();
    printf("compact_ideal_multiplication done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: realistic_scale...\n");
    fflush(stdout);
    res |= quat_test_mlll_realistic_scale();
    printf("realistic_scale done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: tau_path...\n");
    fflush(stdout);
    res |= quat_test_mlll_tau_path();
    printf("tau_path done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: gram_equivalence...\n");
    fflush(stdout);
    res |= quat_test_mlll_gram_equivalence();
    printf("gram_equivalence done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: lideal_create_gram_equivalence...\n");
    fflush(stdout);
    res |= quat_test_lideal_create_gram_equivalence();
    printf("lideal_create_gram_equivalence done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: gram_prealloc_equivalence...\n");
    fflush(stdout);
    res |= quat_test_mlll_gram_prealloc_equivalence();
    printf("gram_prealloc_equivalence done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: all_zero_generators...\n");
    fflush(stdout);
    res |= quat_test_mlll_all_zero_generators();
    printf("all_zero_generators done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: single_generator...\n");
    fflush(stdout);
    res |= quat_test_mlll_single_generator();
    printf("single_generator done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: two_generators_dependent...\n");
    fflush(stdout);
    res |= quat_test_mlll_two_generators_dependent();
    printf("two_generators_dependent done (res=%d)\n", res);
    fflush(stdout);

    printf("Starting MLLL test: rank_deficient...\n");
    fflush(stdout);
    res |= quat_test_mlll_rank_deficient();
    printf("rank_deficient done (res=%d)\n", res);
    fflush(stdout);

    return res;
}
