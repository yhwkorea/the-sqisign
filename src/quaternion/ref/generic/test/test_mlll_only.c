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

    return res;
}
