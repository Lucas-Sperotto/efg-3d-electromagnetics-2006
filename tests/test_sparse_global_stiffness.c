#include "global_stiffness.h"
#include "sparse_matrix.h"
#include "cube_problem.h"
#include "support.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    exit(1);
}

/*
 * For a 3×3×3 node cloud and 3×3×3 integration cells the sparse and dense K
 * matrices must agree to double precision (entries matched via CSR matvec).
 *
 * The comparison is done via y_dense = K_dense * e_i and y_sparse = K_csr * e_i
 * for each standard basis vector e_i, verifying ||y_dense - y_sparse|| < tol.
 */
static void test_sparse_equals_dense_3x3x3(void)
{
    enum { N = 27, CELLS = 3 };
    const double L = 10.0;
    const double TOL = 1e-10;

    Node3D nodes[N];
    double K_dense[N * N];
    SparseCOO coo;
    SparseCSR csr;
    int node_count = 0;

    cube_generate_regular_nodes(L, CELLS, CELLS, CELLS, nodes, N, &node_count);
    support_assign_article_default(nodes, node_count);

    /* Assemble dense K. */
    int rc = global_stiffness_assemble_dense(nodes, node_count,
                                              0.0, L, 0.0, L, 0.0, L,
                                              CELLS, CELLS, CELLS,
                                              K_dense);
    if (rc != GLOBAL_STIFFNESS_OK) {
        fail("dense assembly failed");
    }

    /* Assemble sparse K. */
    int initial_cap = node_count * node_count * 8; /* generous initial capacity */
    sparse_coo_create(&coo, node_count, node_count, initial_cap);

    rc = global_stiffness_assemble_sparse_coo(nodes, node_count,
                                               0.0, L, 0.0, L, 0.0, L,
                                               CELLS, CELLS, CELLS,
                                               &coo);
    if (rc != GLOBAL_STIFFNESS_OK) {
        fail("sparse assembly failed");
    }

    sparse_coo_to_csr(&coo, &csr);

    /* Compare column by column: K_dense * e_i vs K_csr * e_i. */
    double e[N];
    double y_dense[N];
    double y_sparse[N];

    for (int col = 0; col < node_count; ++col) {
        /* Unit vector e_col. */
        for (int k = 0; k < node_count; ++k) {
            e[k] = 0.0;
        }
        e[col] = 1.0;

        /* Dense matvec: y_dense = K_dense * e. */
        for (int row = 0; row < node_count; ++row) {
            double sum = 0.0;
            for (int k = 0; k < node_count; ++k) {
                sum += K_dense[row * node_count + k] * e[k];
            }
            y_dense[row] = sum;
        }

        /* Sparse matvec. */
        sparse_csr_matvec(&csr, e, y_sparse);

        /* Compare. */
        for (int row = 0; row < node_count; ++row) {
            if (fabs(y_dense[row] - y_sparse[row]) > TOL) {
                fprintf(stderr,
                        "FAIL sparse_K != dense_K at K[%d][%d]: "
                        "dense=%.15e sparse=%.15e diff=%.3e\n",
                        row, col, y_dense[row], y_sparse[row],
                        fabs(y_dense[row] - y_sparse[row]));
                exit(1);
            }
        }
    }

    printf("PASS test_sparse_equals_dense_3x3x3: nnz=%d (dense=%d)\n",
           csr.nnz, node_count * node_count);

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
}

int main(void)
{
    test_sparse_equals_dense_3x3x3();
    return 0;
}
