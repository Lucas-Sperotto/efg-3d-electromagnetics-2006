#include "analytical.h"
#include "cube_problem.h"
#include "dense_solver.h"
#include "dirichlet.h"
#include "efg/gmres.h"
#include "error_norm.h"
#include "global_stiffness.h"
#include "lagrange_system.h"
#include "mls.h"
#include "mls_diagnostic.h"
#include "nodes.h"
#include "sparse_matrix.h"
#include "support.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define ANALYTICAL_TERMS 25

static double elapsed_s(clock_t start)
{
    return (double)(clock() - start) / (double)CLOCKS_PER_SEC;
}

/*
 * Sparse EFG pipeline for one problem configuration.
 *
 * Runs: node generation → MLS diagnostic → dense K/G assembly →
 * K_dense→K_coo conversion (filtering exact zeros) → sparse augmented
 * system → GMRES → error vs analytical solution.
 *
 * When run_dense_cmp != 0 the dense solver is also run and the two
 * solutions are compared.  Use only for small cases (total_size < ~400).
 *
 * Returns 0 on success, 1 on any fatal error.
 */
static int run_case(const char *label,
                    double L, double V0,
                    int nx, int ny, int nz,
                    int nx_cells, int ny_cells, int nz_cells,
                    double gmres_tol, int gmres_restart, int gmres_max_iter,
                    int sample_n,
                    int run_dense_cmp)
{
    int node_count      = 0;
    int point_count     = 0;
    int total_size      = 0;
    int mls_eval_fail   = 0;
    int coo_count_aug   = 0;
    double t_asm = 0.0, t_sol = 0.0;
    clock_t t_start;

    Node3D         *nodes     = NULL;
    DirichletPoint *dp        = NULL;
    double         *K_dense   = NULL;
    double         *G_dense   = NULL;
    double         *F         = NULL;
    double         *q_vec     = NULL;
    double         *b_aug     = NULL;
    double         *x_sol     = NULL;
    MlsShapeValue  *shape_buf = NULL;

    SparseCOO K_coo     = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCOO A_aug_coo = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCSR A_aug_csr = {NULL, NULL, NULL, 0, 0, 0};
    GmresResult       gmres_res  = {0, 0, 0.0, 0.0};
    MlsConnectivityStats diag    = {0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};

    double dense_rel_diff = -1.0;  /* < 0 means comparison was not run */
    int rc = 1;

    const int expected_nodes       = cube_regular_node_count(nx, ny, nz);
    const int expected_constraints = cube_dirichlet_point_count(nx, ny, nz);

    /* ---------------------------------------------------------------- alloc */

    nodes     = malloc((size_t)expected_nodes * sizeof(*nodes));
    dp        = malloc((size_t)expected_constraints * sizeof(*dp));
    K_dense   = calloc((size_t)expected_nodes * (size_t)expected_nodes,
                       sizeof(double));
    G_dense   = calloc((size_t)expected_constraints * (size_t)expected_nodes,
                       sizeof(double));
    F         = calloc((size_t)expected_nodes,       sizeof(double));
    q_vec     = calloc((size_t)expected_constraints, sizeof(double));
    shape_buf = malloc((size_t)expected_nodes * sizeof(*shape_buf));

    if (!nodes || !dp || !K_dense || !G_dense || !F || !q_vec || !shape_buf) {
        fprintf(stderr, "[%s] allocation failed\n", label);
        goto done;
    }

    /* ------------------------------------------------------- nodes & support */

    if (cube_generate_regular_nodes(L, nx, ny, nz,
                                    nodes, expected_nodes,
                                    &node_count) != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[%s] node generation failed\n", label);
        goto done;
    }

    if (support_assign_article_default(nodes, node_count) != SUPPORT_OK) {
        fprintf(stderr, "[%s] support assignment failed\n", label);
        goto done;
    }

    /* ---------------------------------------------------- Dirichlet points */

    if (cube_generate_dirichlet_points_from_nodes(L, V0,
                                                  nodes, node_count,
                                                  dp, expected_constraints,
                                                  &point_count)
        != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[%s] Dirichlet generation failed\n", label);
        goto done;
    }

    total_size = node_count + point_count;

    /* -------------------------------------------------- MLS diagnostic */
    /*
     * Evaluate connectivity at Gauss cell centres.  For nx_cells cells over
     * [0,L], the centres lie at (i + 0.5) * L / nx_cells for i = 0..nx_cells-1.
     * mls_connectivity_stats_uniform_grid places nx points uniformly in
     * [xmin, xmax], so setting xmin = half_cell and xmax = L - half_cell
     * gives exactly the cell centres.
     */
    {
        double hx = 0.5 * L / (double)nx_cells;
        double hy = 0.5 * L / (double)ny_cells;
        double hz = 0.5 * L / (double)nz_cells;
        mls_connectivity_stats_uniform_grid(nodes, node_count,
                                            hx, L - hx, nx_cells,
                                            hy, L - hy, ny_cells,
                                            hz, L - hz, nz_cells,
                                            &diag);
    }

    /* ----------------------------------------- dense K and G assembly */

    t_start = clock();

    if (global_stiffness_assemble_dense(nodes, node_count,
                                        0.0, L, 0.0, L, 0.0, L,
                                        nx_cells, ny_cells, nz_cells,
                                        K_dense) != GLOBAL_STIFFNESS_OK) {
        fprintf(stderr, "[%s] K assembly failed\n", label);
        goto done;
    }

    if (dirichlet_assemble_constraints_dense(nodes, node_count,
                                             dp, point_count,
                                             G_dense, q_vec) != DIRICHLET_OK) {
        fprintf(stderr, "[%s] G assembly failed\n", label);
        goto done;
    }

    /* --------------------------------- K_dense → K_coo (exact zeros dropped) */
    /*
     * Using the dense K as intermediate avoids the O(cells × Gauss² × n_active²)
     * memory that the direct sparse assembler would require before compression.
     * For n ≤ ~2000 nodes the dense K is small (≤ 32 MB) and the conversion
     * loop is negligible.
     */
    {
        int k_cap = node_count * 30 + 16;
        if (sparse_coo_create(&K_coo, node_count, node_count, k_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] K_coo alloc failed\n", label);
            goto done;
        }
        for (int r = 0; r < node_count; ++r) {
            for (int c = 0; c < node_count; ++c) {
                double v = K_dense[r * node_count + c];
                if (v != 0.0) {
                    sparse_coo_add(&K_coo, r, c, v);
                }
            }
        }
    }

    /* -------------------- augmented COO: K block + G / G^T (zeros dropped) */
    /*
     * The G matrix is calloc-initialised; entries for nodes outside the
     * support of each Dirichlet point are exactly 0.0 and are skipped here.
     * This keeps the COO sparse: only the ~n_active ≈ 20-30 non-zero entries
     * per constraint row enter the sparse structure.
     */
    {
        int aug_cap = K_coo.count + point_count * 60 + 16;
        if (sparse_coo_create(&A_aug_coo, total_size, total_size, aug_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] A_aug_coo alloc failed\n", label);
            goto done;
        }

        /* K block (rows and cols 0..node_count-1) */
        for (int e = 0; e < K_coo.count; ++e) {
            sparse_coo_add(&A_aug_coo,
                           K_coo.row[e], K_coo.col[e], K_coo.val[e]);
        }

        /* G and G^T blocks */
        for (int ci = 0; ci < point_count; ++ci) {
            for (int ni = 0; ni < node_count; ++ni) {
                double val = G_dense[ci * node_count + ni];
                if (val != 0.0) {
                    /* G^T: A_aug[ni, node_count+ci] */
                    sparse_coo_add(&A_aug_coo, ni, node_count + ci, val);
                    /* G:   A_aug[node_count+ci, ni] */
                    sparse_coo_add(&A_aug_coo, node_count + ci, ni, val);
                }
            }
        }
    }

    /* b_aug: F (zeros) followed by q */
    b_aug = calloc((size_t)total_size, sizeof(double));
    if (!b_aug) { fprintf(stderr, "[%s] b_aug alloc failed\n", label); goto done; }
    for (int ci = 0; ci < point_count; ++ci) {
        b_aug[node_count + ci] = q_vec[ci];
    }

    coo_count_aug = A_aug_coo.count;
    sparse_coo_to_csr(&A_aug_coo, &A_aug_csr);

    t_asm = elapsed_s(t_start);

    /* ------------------------------------------------------------ GMRES */

    x_sol = calloc((size_t)total_size, sizeof(double));
    if (!x_sol) { fprintf(stderr, "[%s] x_sol alloc failed\n", label); goto done; }

    t_start = clock();

    if (gmres_solve(&A_aug_csr, b_aug, x_sol,
                    gmres_tol, gmres_max_iter, gmres_restart,
                    &gmres_res) != GMRES_OK) {
        fprintf(stderr, "[%s] gmres_solve returned error\n", label);
        goto done;
    }

    t_sol = elapsed_s(t_start);

    /* --------------------------------- dense comparison (small cases only) */

    if (run_dense_cmp) {
        double *A_aug_d = calloc((size_t)total_size * (size_t)total_size,
                                 sizeof(double));
        double *b_aug_d = calloc((size_t)total_size, sizeof(double));
        double *x_d     = calloc((size_t)total_size, sizeof(double));

        if (A_aug_d && b_aug_d && x_d) {
            lagrange_assemble_augmented_system_dense(K_dense, F, node_count,
                                                     G_dense, q_vec,
                                                     point_count,
                                                     A_aug_d, b_aug_d);
            if (dense_solve(A_aug_d, b_aug_d, total_size, x_d)
                == DENSE_SOLVER_OK) {
                double diff2 = 0.0, ref2 = 0.0;
                for (int i = 0; i < total_size; ++i) {
                    double d = x_sol[i] - x_d[i];
                    diff2 += d * d;
                    ref2  += x_d[i] * x_d[i];
                }
                dense_rel_diff = (ref2 > 0.0) ? sqrt(diff2 / ref2)
                                               : sqrt(diff2);
            }
        }

        free(A_aug_d);
        free(b_aug_d);
        free(x_d);
    }

    /* ----------------------------------------------- error vs analytical */

    {
        double global_err2 = 0.0, global_ref2 = 0.0;
        double intrnl_err2 = 0.0, intrnl_ref2 = 0.0;
        double max_abs_err = 0.0;
        int valid_all = 0, valid_int = 0;

        for (int ix = 0; ix < sample_n; ++ix) {
            double sx = (double)ix / (double)(sample_n - 1) * L;

            for (int iy = 0; iy < sample_n; ++iy) {
                double sy = (double)iy / (double)(sample_n - 1) * L;

                for (int iz = 0; iz < sample_n; ++iz) {
                    double sz = (double)iz / (double)(sample_n - 1) * L;
                    int val_count = 0;
                    double v_num  = 0.0;

                    int mls_st = mls_linear3d_shape_functions(
                                     nodes, node_count, sx, sy, sz,
                                     shape_buf, node_count, &val_count);
                    if (mls_st != MLS_OK) {
                        ++mls_eval_fail;
                        continue;
                    }

                    for (int k = 0; k < val_count; ++k) {
                        v_num += shape_buf[k].phi
                                 * x_sol[shape_buf[k].node_index];
                    }

                    double v_exact = analytical_potential_cube(
                                         sx, sy, sz, L, V0,
                                         ANALYTICAL_TERMS);
                    double d = v_num - v_exact;
                    double ae = fabs(d);

                    if (ae > max_abs_err) max_abs_err = ae;
                    global_err2 += d * d;
                    global_ref2 += v_exact * v_exact;
                    ++valid_all;

                    int interior = (ix > 0 && ix < sample_n - 1 &&
                                    iy > 0 && iy < sample_n - 1 &&
                                    iz > 0 && iz < sample_n - 1);
                    if (interior) {
                        intrnl_err2 += d * d;
                        intrnl_ref2 += v_exact * v_exact;
                        ++valid_int;
                    }
                }
            }
        }

        double rel_global = (global_ref2 > 0.0)
                            ? sqrt(global_err2 / global_ref2) : 0.0;
        double rel_intrnl = (intrnl_ref2 > 0.0)
                            ? sqrt(intrnl_err2 / intrnl_ref2) : 0.0;

        /* ------------------------------------------------- print report */

        printf("\n");
        printf("=== %s ===\n", label);
        printf("\n");

        printf("--- Problem ---\n");
        printf("  nodes:                  %d\n",     node_count);
        printf("  constraints:            %d\n",     point_count);
        printf("  augmented size:         %d\n",     total_size);
        printf("  node grid:              %dx%dx%d\n", nx, ny, nz);
        printf("  integration cells:      %dx%dx%d\n",
               nx_cells, ny_cells, nz_cells);
        printf("\n");

        printf("--- MLS diagnostic  (cell-centre grid %dx%dx%d = %d pts) ---\n",
               nx_cells, ny_cells, nz_cells, diag.n_total);
        printf("  invalid (active < 4):   %d\n", diag.n_invalid);
        printf("  suspect (active < 8):   %d\n", diag.n_suspect);
        printf("  moment matrix failures: %d\n", diag.n_moment_fail);
        printf("  min active nodes:       %d\n", diag.min_nodes);
        printf("  max active nodes:       %d\n", diag.max_nodes);
        printf("  mean active nodes:      %.2f\n", diag.mean_nodes);
        printf("  max cond estimate:      %.3e\n", diag.max_cond);
        printf("  mean cond estimate:     %.3e\n", diag.mean_cond);
        printf("\n");

        printf("--- Sparse assembly ---\n");
        printf("  K nnz  (|v| > 0):       %d\n", K_coo.count);
        printf("  A_aug COO entries:       %d\n", coo_count_aug);
        printf("  A_aug CSR nnz:           %d\n", A_aug_csr.nnz);
        printf("  assembly time:           %.3f s\n", t_asm);
        printf("\n");

        printf("--- GMRES ---\n");
        printf("  restart:                %d\n",    gmres_restart);
        printf("  max iter:               %d\n",    gmres_max_iter);
        printf("  tolerance:              %.1e\n",  gmres_tol);
        printf("  iterations:             %d\n",    gmres_res.iterations);
        printf("  residual init:          %.6e\n",  gmres_res.residual_init);
        printf("  residual final:         %.6e\n",  gmres_res.residual_final);
        printf("  rel residual:           %.3e\n",
               (gmres_res.residual_init > 0.0)
               ? gmres_res.residual_final / gmres_res.residual_init
               : 0.0);
        printf("  converged:              %s\n",
               gmres_res.converged ? "YES" : "NO");
        printf("  solution time:          %.3f s\n", t_sol);
        if (dense_rel_diff >= 0.0) {
            printf("  GMRES vs dense diff:    %.3e\n", dense_rel_diff);
        }
        printf("\n");

        printf("--- Errors  (sample grid %dx%dx%d,"
               " valid=%d interior=%d) ---\n",
               sample_n, sample_n, sample_n, valid_all, valid_int);
        printf("  max abs error:          %.6e\n", max_abs_err);
        printf("  rel error (global):     %.6e\n", rel_global);
        printf("  rel error (interior):   %.6e\n", rel_intrnl);
        if (mls_eval_fail > 0) {
            printf("  MLS eval failures:      %d\n", mls_eval_fail);
        }
        printf("\n");
    }

    rc = 0;

done:
    free(nodes);
    free(dp);
    free(K_dense);
    free(G_dense);
    free(F);
    free(q_vec);
    free(b_aug);
    free(x_sol);
    free(shape_buf);

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);

    return rc;
}

/*
 * Sparse EFG cube solver.
 *
 * Runs two cases:
 *   1. 5x5x5 sanity check  — compared against the dense solver.
 *   2. 11x11x11 main case  — GMRES only (dense solve is O(n^3) too slow).
 */
int main(void)
{
    int status = 0;

    /* --- Case 1: sanity check --- */
    status |= run_case(
        "5x5x5 sanity check (dense comparison enabled)",
        /*L=*/10.0, /*V0=*/10.0,
        /*nx,ny,nz=*/5, 5, 5,
        /*cells=*/5, 5, 5,
        /*gmres_tol=*/1e-9,
        /*restart=*/30,
        /*max_iter=*/2000,
        /*sample_n=*/11,
        /*dense_cmp=*/1);

    /* --- Case 2: 11x11x11 --- */
    status |= run_case(
        "11x11x11  (15x15x15 integration cells, GMRES only)",
        /*L=*/10.0, /*V0=*/10.0,
        /*nx,ny,nz=*/11, 11, 11,
        /*cells=*/15, 15, 15,
        /*gmres_tol=*/1e-9,
        /*restart=*/200,
        /*max_iter=*/10000,
        /*sample_n=*/11,
        /*dense_cmp=*/0);

    return status;
}
