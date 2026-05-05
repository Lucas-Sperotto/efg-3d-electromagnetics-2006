#include "lagrange_system.h"
#include "sparse_matrix.h"
#include "global_stiffness.h"
#include "dirichlet.h"
#include "cube_problem.h"
#include "support.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    exit(1);
}

/*
 * Dense matvec: y = A * x, where A is stored row-major in a flat array of
 * size n*n.
 */
static void dense_matvec(const double *A, const double *x, double *y, int n)
{
    for (int row = 0; row < n; ++row) {
        double sum = 0.0;
        for (int col = 0; col < n; ++col) {
            sum += A[row * n + col] * x[col];
        }
        y[row] = sum;
    }
}

/*
 * Compare sparse and dense augmented systems using the manual 3-node /
 * 2-constraint system from test_lagrange_system.c.
 *
 *   K = tridiag(2,-1,-1,2,-1,0,-1,2)  (3×3 symmetric)
 *   G = [[1,0,0],[0,0,1]]
 *   F = 0, q = [0, 10]
 *
 * The test builds K as a SparseCOO, calls both assemblers, converts the
 * sparse result to CSR, then compares A_aug * e_i and b_aug column by column.
 */
static void test_sparse_aug_equals_dense_manual(void)
{
    enum { N = 3, M = 2, TOTAL = 5 };
    const double TOL = 1e-12;

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

    double A_dense[TOTAL * TOTAL];
    double b_dense[TOTAL];
    double b_sparse[TOTAL];

    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;

    /* Assemble dense augmented system. */
    if (lagrange_assemble_augmented_system_dense(K_data, F_data, N,
                                                  G_data, q_data, M,
                                                  A_dense, b_dense)
        != LAGRANGE_SYSTEM_OK) {
        fail("dense augmented assembly failed");
    }

    /* Build K as SparseCOO (insert all non-zero entries). */
    sparse_coo_create(&K_coo, N, N, N * N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (fabs(K_data[r * N + c]) > 0.0) {
                sparse_coo_add(&K_coo, r, c, K_data[r * N + c]);
            }
        }
    }

    /* Assemble sparse augmented COO. */
    sparse_coo_create(&A_aug_coo, TOTAL, TOTAL,
                      K_coo.count + 2 * M * N + 4);
    if (lagrange_assemble_augmented_sparse_coo(&K_coo, F_data, N,
                                               G_data, q_data, M,
                                               &A_aug_coo, b_sparse)
        != LAGRANGE_SYSTEM_OK) {
        fail("sparse augmented assembly failed");
    }

    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    /* Compare b_aug. */
    for (int i = 0; i < TOTAL; ++i) {
        if (fabs(b_dense[i] - b_sparse[i]) > TOL) {
            fprintf(stderr,
                    "FAIL b_aug[%d]: dense=%.15e sparse=%.15e\n",
                    i, b_dense[i], b_sparse[i]);
            exit(1);
        }
    }

    /* Compare A_aug column by column via unit-vector matvec. */
    double e[TOTAL];
    double y_dense[TOTAL];
    double y_sparse[TOTAL];

    for (int col = 0; col < TOTAL; ++col) {
        memset(e, 0, sizeof(e));
        e[col] = 1.0;

        dense_matvec(A_dense, e, y_dense, TOTAL);
        sparse_csr_matvec(&A_aug_csr, e, y_sparse);

        for (int row = 0; row < TOTAL; ++row) {
            if (fabs(y_dense[row] - y_sparse[row]) > TOL) {
                fprintf(stderr,
                        "FAIL A_aug[%d][%d]: dense=%.15e sparse=%.15e diff=%.3e\n",
                        row, col, y_dense[row], y_sparse[row],
                        fabs(y_dense[row] - y_sparse[row]));
                exit(1);
            }
        }
    }

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);

    printf("PASS test_sparse_aug_equals_dense_manual: "
           "total_size=%d nnz=%d\n", TOTAL, A_aug_csr.nnz);
}

/*
 * End-to-end test using the 3×3×3 regular cube problem.
 *
 * Assembles both dense and sparse augmented systems and compares them via
 * matvec with each standard basis vector.
 */
static void test_sparse_aug_equals_dense_3x3x3(void)
{
    enum { CELLS = 3, MAX_NODES = 27, MAX_CONSTRAINTS = 400 };
    const double L   = 10.0;
    const double V0  = 10.0;
    const double TOL = 1e-10;

    Node3D        nodes[MAX_NODES];
    DirichletPoint dp[MAX_CONSTRAINTS];
    int node_count    = 0;
    int point_count   = 0;
    int total_size;

    double *K_dense    = NULL;
    double *G_dense    = NULL;
    double *F          = NULL;
    double *q          = NULL;
    double *A_dense    = NULL;
    double *b_dense    = NULL;
    double *b_sparse   = NULL;
    double *e_vec      = NULL;
    double *y_dense    = NULL;
    double *y_sparse   = NULL;

    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;

    cube_generate_regular_nodes(L, CELLS, CELLS, CELLS,
                                nodes, MAX_NODES, &node_count);
    support_assign_article_default(nodes, node_count);

    cube_generate_dirichlet_points_from_nodes(L, V0,
                                              nodes, node_count,
                                              dp, MAX_CONSTRAINTS,
                                              &point_count);

    total_size = node_count + point_count;

    K_dense  = calloc((size_t)(node_count * node_count), sizeof(double));
    G_dense  = calloc((size_t)(point_count * node_count), sizeof(double));
    F        = calloc((size_t)node_count,  sizeof(double));
    q        = calloc((size_t)point_count, sizeof(double));
    A_dense  = calloc((size_t)(total_size * total_size), sizeof(double));
    b_dense  = calloc((size_t)total_size, sizeof(double));
    b_sparse = calloc((size_t)total_size, sizeof(double));
    e_vec    = calloc((size_t)total_size, sizeof(double));
    y_dense  = calloc((size_t)total_size, sizeof(double));
    y_sparse = calloc((size_t)total_size, sizeof(double));

    if (!K_dense || !G_dense || !F || !q || !A_dense ||
        !b_dense || !b_sparse || !e_vec || !y_dense || !y_sparse) {
        fail("allocation failed");
    }

    /* Assemble dense K. */
    if (global_stiffness_assemble_dense(nodes, node_count,
                                         0.0, L, 0.0, L, 0.0, L,
                                         CELLS, CELLS, CELLS,
                                         K_dense) != GLOBAL_STIFFNESS_OK) {
        fail("dense K assembly failed");
    }

    /* Assemble G (Dirichlet constraints). */
    if (dirichlet_assemble_constraints_dense(nodes, node_count,
                                              dp, point_count,
                                              G_dense, q) != DIRICHLET_OK) {
        fail("dirichlet assembly failed");
    }

    /* Dense augmented system. */
    if (lagrange_assemble_augmented_system_dense(K_dense, F, node_count,
                                                  G_dense, q, point_count,
                                                  A_dense, b_dense)
        != LAGRANGE_SYSTEM_OK) {
        fail("dense augmented assembly failed");
    }

    /* Sparse K. */
    sparse_coo_create(&K_coo, node_count, node_count,
                      node_count * node_count * 4);
    if (global_stiffness_assemble_sparse_coo(nodes, node_count,
                                              0.0, L, 0.0, L, 0.0, L,
                                              CELLS, CELLS, CELLS,
                                              &K_coo) != GLOBAL_STIFFNESS_OK) {
        fail("sparse K assembly failed");
    }

    /* Sparse augmented system. */
    sparse_coo_create(&A_aug_coo, total_size, total_size,
                      K_coo.count + 2 * point_count * node_count + 4);
    if (lagrange_assemble_augmented_sparse_coo(&K_coo, F, node_count,
                                               G_dense, q, point_count,
                                               &A_aug_coo, b_sparse)
        != LAGRANGE_SYSTEM_OK) {
        fail("sparse augmented assembly failed");
    }

    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    /* Compare b_aug. */
    for (int i = 0; i < total_size; ++i) {
        if (fabs(b_dense[i] - b_sparse[i]) > TOL) {
            fprintf(stderr, "FAIL b_aug[%d]: dense=%.15e sparse=%.15e\n",
                    i, b_dense[i], b_sparse[i]);
            exit(1);
        }
    }

    /* Compare A_aug column by column. */
    for (int col = 0; col < total_size; ++col) {
        memset(e_vec, 0, (size_t)total_size * sizeof(double));
        e_vec[col] = 1.0;

        dense_matvec(A_dense, e_vec, y_dense, total_size);
        sparse_csr_matvec(&A_aug_csr, e_vec, y_sparse);

        for (int row = 0; row < total_size; ++row) {
            if (fabs(y_dense[row] - y_sparse[row]) > TOL) {
                fprintf(stderr,
                        "FAIL A_aug[%d][%d]: dense=%.15e sparse=%.15e diff=%.3e\n",
                        row, col, y_dense[row], y_sparse[row],
                        fabs(y_dense[row] - y_sparse[row]));
                exit(1);
            }
        }
    }

    printf("PASS test_sparse_aug_equals_dense_3x3x3: "
           "nodes=%d constraints=%d total=%d nnz=%d\n",
           node_count, point_count, total_size, A_aug_csr.nnz);

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);

    free(K_dense);
    free(G_dense);
    free(F);
    free(q);
    free(A_dense);
    free(b_dense);
    free(b_sparse);
    free(e_vec);
    free(y_dense);
    free(y_sparse);
}

int main(void)
{
    test_sparse_aug_equals_dense_manual();
    test_sparse_aug_equals_dense_3x3x3();
    return 0;
}
