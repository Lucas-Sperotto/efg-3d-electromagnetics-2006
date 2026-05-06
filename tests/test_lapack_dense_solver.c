#include "lapack_dense_solver.h"
#include "lagrange_system.h"
#include "sparse_matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    exit(1);
}

static void expect_close(const char *name,
                         double actual,
                         double expected,
                         double tolerance)
{
    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr,
                "FAIL %s: expected %.17g got %.17g\n",
                name, expected, actual);
        exit(1);
    }
}

static void test_csr_to_dense_column_major(void)
{
    SparseCOO coo;
    SparseCSR csr;
    double dense[9];

    if (sparse_coo_create(&coo, 3, 3, 6) != SPARSE_OK) {
        fail("coo create failed");
    }

    sparse_coo_add(&coo, 0, 0, 1.0);
    sparse_coo_add(&coo, 0, 2, 2.0);
    sparse_coo_add(&coo, 1, 1, 3.0);
    sparse_coo_add(&coo, 2, 0, 4.0);
    sparse_coo_add(&coo, 2, 2, 5.0);

    if (sparse_coo_to_csr(&coo, &csr) != SPARSE_OK) {
        fail("csr conversion failed");
    }

    if (sparse_csr_to_dense_column_major(&csr, dense, 3)
        != LAPACK_DENSE_SOLVER_OK) {
        fail("CSR -> dense column-major conversion failed");
    }

    expect_close("A(0,0)", dense[0 + 0 * 3], 1.0, 0.0);
    expect_close("A(1,0)", dense[1 + 0 * 3], 0.0, 0.0);
    expect_close("A(2,0)", dense[2 + 0 * 3], 4.0, 0.0);
    expect_close("A(0,1)", dense[0 + 1 * 3], 0.0, 0.0);
    expect_close("A(1,1)", dense[1 + 1 * 3], 3.0, 0.0);
    expect_close("A(2,1)", dense[2 + 1 * 3], 0.0, 0.0);
    expect_close("A(0,2)", dense[0 + 2 * 3], 2.0, 0.0);
    expect_close("A(1,2)", dense[1 + 2 * 3], 0.0, 0.0);
    expect_close("A(2,2)", dense[2 + 2 * 3], 5.0, 0.0);

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);

    printf("PASS test_csr_to_dense_column_major\n");
}

static void test_lapack_dense_manual_5x5(void)
{
    enum { N = 3, M = 2, TOTAL = 5 };

    const double K_data[N * N] = {
         2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0
    };
    const double F_data[N] = {0.0, 0.0, 0.0};
    const double G_data[M * N] = {
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0
    };
    const double q_data[M] = {0.0, 10.0};
    const double expected[TOTAL] = {0.0, 5.0, 10.0, 5.0, -15.0};

    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;
    double b_aug[TOTAL];
    double x[TOTAL];
    LapackDenseSolveResult result;
    int status;

    sparse_coo_create(&K_coo, N, N, N * N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (fabs(K_data[r * N + c]) > 0.0) {
                sparse_coo_add(&K_coo, r, c, K_data[r * N + c]);
            }
        }
    }

    sparse_coo_create(&A_aug_coo, TOTAL, TOTAL,
                      K_coo.count + 2 * M * N + 4);
    if (lagrange_assemble_augmented_sparse_coo(&K_coo, F_data, N,
                                               G_data, q_data, M,
                                               &A_aug_coo, b_aug)
        != LAGRANGE_SYSTEM_OK) {
        fail("augmented sparse assembly failed");
    }
    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    memset(x, 0, sizeof(x));
    status = solve_augmented_lapack_dense(&A_aug_csr, b_aug, x, &result);

#ifdef EFG_HAVE_LAPACK
    if (status != LAPACK_DENSE_SOLVER_OK) {
        fprintf(stderr, "FAIL LAPACK solve status=%d info=%d\n",
                status, result.lapack_info);
        exit(1);
    }
    for (int i = 0; i < TOTAL; ++i) {
        expect_close("LAPACK x", x[i], expected[i], 1.0e-10);
    }
    if (result.residual_final > 1.0e-10) {
        fprintf(stderr,
                "FAIL LAPACK residual %.3e exceeds tolerance\n",
                result.residual_final);
        exit(1);
    }
    printf("PASS test_lapack_dense_manual_5x5: "
           "residual=%.3e conversion=%.3e solve=%.3e\n",
           result.residual_final,
           result.conversion_time_s,
           result.solve_time_s);
#else
    if (status != LAPACK_DENSE_SOLVER_UNAVAILABLE || result.available != 0) {
        fail("LAPACK unavailable path returned unexpected status");
    }
    printf("PASS test_lapack_dense_manual_5x5: LAPACK unavailable\n");
#endif

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);
}

int main(void)
{
    test_csr_to_dense_column_major();
    test_lapack_dense_manual_5x5();
    return 0;
}
