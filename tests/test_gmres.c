#include "efg/gmres.h"
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
 * Solve the 5×5 augmented system from the manual Lagrange test with GMRES
 * and compare the solution to the known exact result.
 *
 * Known solution (from test_lagrange_solve.c): [0.0, 5.0, 10.0, 5.0, -15.0]
 */
static void test_gmres_manual_5x5(void)
{
    enum { N = 3, M = 2, TOTAL = 5 };
    const double TOL_GMRES  = 1e-10;
    const double TOL_SOL    = 1e-8;

    const double K_data[N * N] = {
         2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0
    };
    const double F_data[N]     = {0.0, 0.0, 0.0};
    const double G_data[M * N] = {
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0
    };
    const double q_data[M] = {0.0, 10.0};

    const double expected[TOTAL] = {0.0, 5.0, 10.0, 5.0, -15.0};

    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;
    GmresResult res;

    double b_aug[TOTAL];
    double x[TOTAL];

    /* Build K_coo */
    sparse_coo_create(&K_coo, N, N, N * N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (fabs(K_data[r * N + c]) > 0.0) {
                sparse_coo_add(&K_coo, r, c, K_data[r * N + c]);
            }
        }
    }

    /* Assemble sparse augmented system */
    sparse_coo_create(&A_aug_coo, TOTAL, TOTAL,
                      K_coo.count + 2 * M * N + 4);
    if (lagrange_assemble_augmented_sparse_coo(&K_coo, F_data, N,
                                               G_data, q_data, M,
                                               &A_aug_coo, b_aug)
        != LAGRANGE_SYSTEM_OK) {
        fail("augmented sparse assembly failed");
    }

    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    /* GMRES solve — zero initial guess */
    memset(x, 0, sizeof(x));
    if (gmres_solve(&A_aug_csr, b_aug, x,
                    TOL_GMRES, 50, 10, NULL, &res) != GMRES_OK) {
        fail("gmres_solve returned error");
    }

    if (!res.converged) {
        fprintf(stderr,
                "FAIL: GMRES did not converge in %d iterations "
                "(r_init=%.3e r_final=%.3e)\n",
                res.iterations, res.residual_init, res.residual_final);
        exit(1);
    }

    for (int i = 0; i < TOTAL; ++i) {
        if (fabs(x[i] - expected[i]) > TOL_SOL) {
            fprintf(stderr,
                    "FAIL solution[%d]: expected %.10g got %.10g\n",
                    i, expected[i], x[i]);
            exit(1);
        }
    }

    printf("PASS test_gmres_manual_5x5: "
           "iters=%d r_init=%.3e r_final=%.3e\n",
           res.iterations, res.residual_init, res.residual_final);

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);
}

/*
 * Solve the 3×3×3 cube augmented system with GMRES and verify the residual
 * falls below tolerance. Compare the solution norm against the dense solver
 * to confirm they agree within 1e-4 relative tolerance.
 */
static void test_gmres_vs_dense_3x3x3(void)
{
    enum { CELLS = 3, MAX_NODES = 27, MAX_CONSTRAINTS = 400 };
    const double L        = 10.0;
    const double V0       = 10.0;
    const double TOL_GMRES = 1e-9;
    const double TOL_CMP   = 1e-4;

    Node3D        nodes[MAX_NODES];
    DirichletPoint dp[MAX_CONSTRAINTS];
    int node_count  = 0;
    int point_count = 0;
    int total_size;

    double *K_dense  = NULL;
    double *G_dense  = NULL;
    double *F        = NULL;
    double *q        = NULL;
    double *A_dense  = NULL;
    double *b_dense  = NULL;
    double *x_dense  = NULL;
    double *x_gmres  = NULL;

    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;
    GmresResult res;

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
    b_dense  = calloc((size_t)total_size,  sizeof(double));
    x_dense  = calloc((size_t)total_size,  sizeof(double));
    x_gmres  = calloc((size_t)total_size,  sizeof(double));

    if (!K_dense || !G_dense || !F || !q ||
        !A_dense || !b_dense || !x_dense || !x_gmres) {
        fail("allocation failed");
    }

    if (global_stiffness_assemble_dense(nodes, node_count,
                                         0.0, L, 0.0, L, 0.0, L,
                                         CELLS, CELLS, CELLS,
                                         K_dense) != GLOBAL_STIFFNESS_OK) {
        fail("dense K assembly failed");
    }

    if (dirichlet_assemble_constraints_dense(nodes, node_count,
                                              dp, point_count,
                                              G_dense, q) != DIRICHLET_OK) {
        fail("dirichlet assembly failed");
    }

    if (lagrange_assemble_augmented_system_dense(K_dense, F, node_count,
                                                  G_dense, q, point_count,
                                                  A_dense, b_dense)
        != LAGRANGE_SYSTEM_OK) {
        fail("dense augmented assembly failed");
    }

    /* Dense solve (using Gauss elimination on A_dense / b_dense copy) */
    {
        double *A_copy = malloc((size_t)(total_size * total_size) * sizeof(double));
        double *b_copy = malloc((size_t)total_size * sizeof(double));
        if (!A_copy || !b_copy) fail("dense copy allocation failed");
        memcpy(A_copy, A_dense, (size_t)(total_size * total_size) * sizeof(double));
        memcpy(b_copy, b_dense, (size_t)total_size * sizeof(double));

        /* Forward elimination with partial pivoting */
        for (int col = 0; col < total_size; ++col) {
            /* Find pivot */
            int pivot = col;
            for (int row = col + 1; row < total_size; ++row) {
                if (fabs(A_copy[row * total_size + col]) >
                    fabs(A_copy[pivot * total_size + col])) {
                    pivot = row;
                }
            }
            /* Swap rows */
            if (pivot != col) {
                for (int k = 0; k < total_size; ++k) {
                    double tmp = A_copy[col * total_size + k];
                    A_copy[col * total_size + k] = A_copy[pivot * total_size + k];
                    A_copy[pivot * total_size + k] = tmp;
                }
                double tmp = b_copy[col]; b_copy[col] = b_copy[pivot]; b_copy[pivot] = tmp;
            }
            double diag = A_copy[col * total_size + col];
            if (fabs(diag) < 1e-30) continue;
            for (int row = col + 1; row < total_size; ++row) {
                double factor = A_copy[row * total_size + col] / diag;
                for (int k = col; k < total_size; ++k) {
                    A_copy[row * total_size + k] -= factor * A_copy[col * total_size + k];
                }
                b_copy[row] -= factor * b_copy[col];
            }
        }
        /* Back substitution */
        for (int row = total_size - 1; row >= 0; --row) {
            double s = b_copy[row];
            for (int k = row + 1; k < total_size; ++k) {
                s -= A_copy[row * total_size + k] * x_dense[k];
            }
            double diag = A_copy[row * total_size + row];
            x_dense[row] = (fabs(diag) > 1e-30) ? s / diag : 0.0;
        }
        free(A_copy);
        free(b_copy);
    }

    /* Sparse augmented system */
    sparse_coo_create(&K_coo, node_count, node_count,
                      node_count * node_count * 4);
    if (global_stiffness_assemble_sparse_coo(nodes, node_count,
                                              0.0, L, 0.0, L, 0.0, L,
                                              CELLS, CELLS, CELLS,
                                              &K_coo) != GLOBAL_STIFFNESS_OK) {
        fail("sparse K assembly failed");
    }

    sparse_coo_create(&A_aug_coo, total_size, total_size,
                      K_coo.count + 2 * point_count * node_count + 4);
    if (lagrange_assemble_augmented_sparse_coo(&K_coo, F, node_count,
                                               G_dense, q, point_count,
                                               &A_aug_coo, b_dense)
        != LAGRANGE_SYSTEM_OK) {
        fail("sparse augmented assembly failed");
    }
    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    /* GMRES solve */
    memset(x_gmres, 0, (size_t)total_size * sizeof(double));
    if (gmres_solve(&A_aug_csr, b_dense, x_gmres,
                    TOL_GMRES, 2000, total_size, NULL, &res) != GMRES_OK) {
        fail("gmres_solve returned error");
    }

    if (!res.converged) {
        fprintf(stderr,
                "FAIL: GMRES did not converge in %d iterations "
                "(r_init=%.3e r_final=%.3e)\n",
                res.iterations, res.residual_init, res.residual_final);
        exit(1);
    }

    /* Compare GMRES solution with dense solution */
    {
        double diff_sq = 0.0;
        double ref_sq  = 0.0;
        for (int i = 0; i < total_size; ++i) {
            double d = x_gmres[i] - x_dense[i];
            diff_sq += d * d;
            ref_sq  += x_dense[i] * x_dense[i];
        }
        double rel_err = (ref_sq > 1e-30) ? sqrt(diff_sq / ref_sq) : sqrt(diff_sq);
        if (rel_err > TOL_CMP) {
            fprintf(stderr,
                    "FAIL GMRES vs dense relative error = %.3e (tol=%.3e)\n",
                    rel_err, TOL_CMP);
            exit(1);
        }
        printf("PASS test_gmres_vs_dense_3x3x3: nodes=%d constraints=%d "
               "iters=%d r_final=%.3e rel_err=%.3e\n",
               node_count, point_count,
               res.iterations, res.residual_final, rel_err);
    }

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);

    free(K_dense); free(G_dense); free(F); free(q);
    free(A_dense); free(b_dense); free(x_dense); free(x_gmres);
}

/*
 * Verify that Jacobi preconditioning (a) gives the same solution as without
 * preconditioner and (b) converges in fewer or equal iterations on a system
 * with large diagonal variation.
 *
 * System: 4×4 diagonal-dominant with entries scaled by 1..4 so the diagonal
 * varies from 20 to 80.  Without preconditioning the condition number is
 * proportional to max/min diagonal ≈ 4.  With Jacobi (D⁻¹A) it is 1.
 * GMRES with restart=4 (full Arnoldi) should converge in ≤ 2 iterations
 * with Jacobi and possibly more without.
 */
static void test_gmres_jacobi_precond(void)
{
    enum { N = 4 };
    const double TOL_GMRES = 1e-10;
    const double TOL_SOL   = 1e-8;

    /* A is diagonally dominant with diagonal 20, 40, 60, 80 */
    const double A_data[N * N] = {
        20.0, -1.0, -1.0, -1.0,
        -1.0, 40.0, -1.0, -1.0,
        -1.0, -1.0, 60.0, -1.0,
        -1.0, -1.0, -1.0, 80.0
    };
    const double b[N] = {1.0, 2.0, 3.0, 4.0};

    SparseCOO coo;
    SparseCSR csr;
    GmresResult res_no, res_jac;
    double x_no[N], x_jac[N];
    double diag_inv[N];

    sparse_coo_create(&coo, N, N, N * N);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (fabs(A_data[r * N + c]) > 0.0) {
                sparse_coo_add(&coo, r, c, A_data[r * N + c]);
            }
        }
    }
    sparse_coo_to_csr(&coo, &csr);

    /* Build Jacobi preconditioner: 1/diag */
    for (int i = 0; i < N; ++i) {
        double d = 0.0;
        for (int p = csr.row_ptr[i]; p < csr.row_ptr[i + 1]; ++p) {
            if (csr.col_ind[p] == i) { d = csr.values[p]; break; }
        }
        diag_inv[i] = (fabs(d) > 1e-15) ? 1.0 / d : 1.0;
    }

    /* Solve without preconditioner */
    memset(x_no, 0, sizeof(x_no));
    if (gmres_solve(&csr, b, x_no, TOL_GMRES, 100, N, NULL, &res_no)
            != GMRES_OK) {
        fail("gmres_solve (no precond) returned error");
    }
    if (!res_no.converged) fail("GMRES (no precond) did not converge");

    /* Solve with Jacobi */
    memset(x_jac, 0, sizeof(x_jac));
    if (gmres_solve(&csr, b, x_jac, TOL_GMRES, 100, N, diag_inv, &res_jac)
            != GMRES_OK) {
        fail("gmres_solve (Jacobi) returned error");
    }
    if (!res_jac.converged) fail("GMRES (Jacobi) did not converge");

    /* Both solutions must agree */
    for (int i = 0; i < N; ++i) {
        if (fabs(x_no[i] - x_jac[i]) > TOL_SOL) {
            fprintf(stderr,
                    "FAIL solution[%d]: no_precond=%.10g jacobi=%.10g\n",
                    i, x_no[i], x_jac[i]);
            exit(1);
        }
    }

    /* Jacobi must not need more iterations than unpreconditioned */
    if (res_jac.iterations > res_no.iterations) {
        fprintf(stderr,
                "FAIL: Jacobi used more iterations (%d) than no-precond (%d)\n",
                res_jac.iterations, res_no.iterations);
        exit(1);
    }

    printf("PASS test_gmres_jacobi_precond: "
           "no_precond iters=%d  jacobi iters=%d  "
           "r_final_no=%.3e  r_final_jac=%.3e\n",
           res_no.iterations, res_jac.iterations,
           res_no.residual_final, res_jac.residual_final);

    sparse_coo_destroy(&coo);
    sparse_csr_destroy(&csr);
}

int main(void)
{
    test_gmres_manual_5x5();
    test_gmres_vs_dense_3x3x3();
    test_gmres_jacobi_precond();
    return 0;
}
