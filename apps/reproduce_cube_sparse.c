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
#define REPORT_CSV_PATH "data/output/reproduce_cube_sparse_report.csv"
#define PLANE15_CSV_PATH "data/output/cube_plane_x_5_33_refine15.csv"

typedef struct CubeSparseReport {
    const char *label;
    int nodes;
    int constraints;
    int augmented_size;
    int K_nnz;
    int A_aug_nnz;
    int support_min;
    int support_max;
    int support_lt_4;
    int support_lt_8;
    int mls_failures;
    int gmres_converged;
    int gmres_iterations;
    double support_mean;
    double max_cond_estimate;
    double residual_initial;
    double residual_final;
    double relative_error_global;
    double relative_error_interior;
    double max_abs_error;
    double assembly_time_s;
    double solve_time_s;
} CubeSparseReport;

typedef enum CubeSparseSelection {
    CUBE_SPARSE_ALL,
    CUBE_SPARSE_SANITY,
    CUBE_SPARSE_TARGET,
    CUBE_SPARSE_REFINE13,
    CUBE_SPARSE_REFINE15,
    CUBE_SPARSE_PLANE15
} CubeSparseSelection;

typedef struct PlaneExportConfig {
    int enabled;
    double x;
    int y_samples;
    int z_samples;
    const char *path;
} PlaneExportConfig;

static double elapsed_s(clock_t start)
{
    return (double)(clock() - start) / (double)CLOCKS_PER_SEC;
}

static void report_init(CubeSparseReport *report, const char *label)
{
    report->label = label;
    report->nodes = -1;
    report->constraints = -1;
    report->augmented_size = -1;
    report->K_nnz = -1;
    report->A_aug_nnz = -1;
    report->support_min = -1;
    report->support_max = -1;
    report->support_lt_4 = -1;
    report->support_lt_8 = -1;
    report->mls_failures = -1;
    report->gmres_converged = 0;
    report->gmres_iterations = -1;
    report->support_mean = NAN;
    report->max_cond_estimate = NAN;
    report->residual_initial = NAN;
    report->residual_final = NAN;
    report->relative_error_global = NAN;
    report->relative_error_interior = NAN;
    report->max_abs_error = NAN;
    report->assembly_time_s = NAN;
    report->solve_time_s = NAN;
}

static void print_required_report(const CubeSparseReport *report)
{
    printf("--- Required report ---\n");
    printf("label: %s\n", report->label);
    printf("nodes: %d\n", report->nodes);
    printf("constraints: %d\n", report->constraints);
    printf("augmented_size: %d\n", report->augmented_size);
    printf("K_nnz: %d\n", report->K_nnz);
    printf("A_aug_nnz: %d\n", report->A_aug_nnz);
    printf("support_min: %d\n", report->support_min);
    printf("support_mean: %.6f\n", report->support_mean);
    printf("support_max: %d\n", report->support_max);
    printf("support_lt_4: %d\n", report->support_lt_4);
    printf("support_lt_8: %d\n", report->support_lt_8);
    printf("mls_failures: %d\n", report->mls_failures);
    printf("max_cond_estimate: %.6e\n", report->max_cond_estimate);
    printf("gmres_converged: %s\n", report->gmres_converged ? "YES" : "NO");
    printf("gmres_iterations: %d\n", report->gmres_iterations);
    printf("residual_initial: %.6e\n", report->residual_initial);
    printf("residual_final: %.6e\n", report->residual_final);
    printf("relative_error_global: %.6e\n", report->relative_error_global);
    printf("relative_error_interior: %.6e\n", report->relative_error_interior);
    printf("max_abs_error: %.6e\n", report->max_abs_error);
    printf("assembly_time_s: %.6f\n", report->assembly_time_s);
    printf("solve_time_s: %.6f\n", report->solve_time_s);
    printf("\n");
}

static int write_report_csv_header(void)
{
    FILE *file = fopen(REPORT_CSV_PATH, "w");

    if (file == NULL) {
        fprintf(stderr, "Could not open report CSV: %s\n", REPORT_CSV_PATH);
        return 1;
    }

    fprintf(file,
            "label,nodes,constraints,augmented_size,K_nnz,A_aug_nnz,"
            "support_min,support_mean,support_max,support_lt_4,"
            "support_lt_8,mls_failures,max_cond_estimate,"
            "gmres_converged,gmres_iterations,"
            "residual_initial,residual_final,relative_error_global,"
            "relative_error_interior,max_abs_error,assembly_time_s,"
            "solve_time_s\n");

    if (fclose(file) != 0) {
        fprintf(stderr, "Could not close report CSV: %s\n", REPORT_CSV_PATH);
        return 1;
    }

    return 0;
}

static int append_required_report_csv(const CubeSparseReport *report)
{
    FILE *file = fopen(REPORT_CSV_PATH, "a");

    if (file == NULL) {
        fprintf(stderr, "Could not append report CSV: %s\n", REPORT_CSV_PATH);
        return 1;
    }

    fprintf(file,
            "\"%s\",%d,%d,%d,%d,%d,%d,%.17g,%d,%d,%d,%d,%.17g,%s,%d,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            report->label,
            report->nodes,
            report->constraints,
            report->augmented_size,
            report->K_nnz,
            report->A_aug_nnz,
            report->support_min,
            report->support_mean,
            report->support_max,
            report->support_lt_4,
            report->support_lt_8,
            report->mls_failures,
            report->max_cond_estimate,
            report->gmres_converged ? "YES" : "NO",
            report->gmres_iterations,
            report->residual_initial,
            report->residual_final,
            report->relative_error_global,
            report->relative_error_interior,
            report->max_abs_error,
            report->assembly_time_s,
            report->solve_time_s);

    if (fclose(file) != 0) {
        fprintf(stderr, "Could not close report CSV: %s\n", REPORT_CSV_PATH);
        return 1;
    }

    return 0;
}

static int emit_required_report(const CubeSparseReport *report)
{
    print_required_report(report);
    return append_required_report_csv(report);
}

static int export_plane_csv(const PlaneExportConfig *config,
                            const Node3D *nodes,
                            int node_count,
                            const double *solution,
                            MlsShapeValue *shape_buf,
                            double L,
                            double V0)
{
    FILE *file = NULL;
    int mls_failures = 0;
    int valid_points = 0;
    int total_points = 0;
    double err2 = 0.0;
    double ref2 = 0.0;
    double max_abs_error = 0.0;
    double sum_abs_error = 0.0;
    double v_num_min = INFINITY;
    double v_num_max = -INFINITY;
    double v_exact_min = INFINITY;
    double v_exact_max = -INFINITY;

    if (config == NULL || nodes == NULL || solution == NULL ||
        shape_buf == NULL || node_count <= 0 ||
        config->y_samples < 2 || config->z_samples < 2 ||
        config->path == NULL) {
        return 1;
    }

    file = fopen(config->path, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open plane CSV: %s\n", config->path);
        return 1;
    }

    fprintf(file, "x,y,z,V_num,V_exact,abs_error\n");

    for (int iy = 0; iy < config->y_samples; ++iy) {
        double y = L * (double)iy / (double)(config->y_samples - 1);

        for (int iz = 0; iz < config->z_samples; ++iz) {
            double z = L * (double)iz / (double)(config->z_samples - 1);
            double v_exact = analytical_potential_cube(config->x, y, z,
                                                       L, V0,
                                                       ANALYTICAL_TERMS);
            int val_count = 0;
            double v_num = 0.0;
            double abs_error = NAN;
            int mls_status;

            if (v_exact < v_exact_min) v_exact_min = v_exact;
            if (v_exact > v_exact_max) v_exact_max = v_exact;

            mls_status = mls_linear3d_shape_functions(nodes, node_count,
                                                       config->x, y, z,
                                                       shape_buf, node_count,
                                                       &val_count);
            if (mls_status != MLS_OK) {
                ++mls_failures;
                fprintf(file, "%.17g,%.17g,%.17g,nan,%.17g,nan\n",
                        config->x, y, z, v_exact);
                ++total_points;
                continue;
            }

            for (int k = 0; k < val_count; ++k) {
                v_num += shape_buf[k].phi
                         * solution[shape_buf[k].node_index];
            }

            abs_error = fabs(v_num - v_exact);
            err2 += (v_num - v_exact) * (v_num - v_exact);
            ref2 += v_exact * v_exact;
            sum_abs_error += abs_error;
            if (abs_error > max_abs_error) max_abs_error = abs_error;
            if (v_num < v_num_min) v_num_min = v_num;
            if (v_num > v_num_max) v_num_max = v_num;

            fprintf(file, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    config->x, y, z, v_num, v_exact, abs_error);

            ++valid_points;
            ++total_points;
        }
    }

    if (fclose(file) != 0) {
        fprintf(stderr, "Could not close plane CSV: %s\n", config->path);
        return 1;
    }

    printf("--- Plane export ---\n");
    printf("  csv:                    %s\n", config->path);
    printf("  x plane:                %.17g\n", config->x);
    printf("  y samples:              %d\n", config->y_samples);
    printf("  z samples:              %d\n", config->z_samples);
    printf("  exported points:        %d\n", total_points);
    printf("  valid points:           %d\n", valid_points);
    printf("  MLS failures:           %d\n", mls_failures);
    printf("  max abs error:          %.6e\n", max_abs_error);
    printf("  mean abs error:         %.6e\n",
           (valid_points > 0) ? sum_abs_error / (double)valid_points : NAN);
    printf("  relative error plane:   %.6e\n",
           (ref2 > 0.0) ? sqrt(err2 / ref2) : NAN);
    printf("  V_num min:              %.6e\n",
           (valid_points > 0) ? v_num_min : NAN);
    printf("  V_num max:              %.6e\n",
           (valid_points > 0) ? v_num_max : NAN);
    printf("  V_exact min:            %.6e\n", v_exact_min);
    printf("  V_exact max:            %.6e\n", v_exact_max);
    printf("\n");

    if (mls_failures > 0) {
        fprintf(stderr, "STOP: plane MLS failures = %d (> 0)\n",
                mls_failures);
        return 1;
    }

    return 0;
}

static void print_usage(const char *program_name)
{
    printf("Usage:\n");
    printf("  %s\n", program_name);
    printf("  %s --case sanity\n", program_name);
    printf("  %s --case target\n", program_name);
    printf("  %s --case refine13\n", program_name);
    printf("  %s --case refine15\n", program_name);
    printf("  %s --case plane15\n", program_name);
    printf("  %s --case all\n", program_name);
    printf("\n");
    printf("Cases:\n");
    printf("  sanity  regular 5x5x5 nodes, 5x5x5 integration cells, dense comparison\n");
    printf("  target  regular 11x11x11 nodes, 15x15x15 integration cells, GMRES\n");
    printf("  refine13 regular 13x13x13 nodes, 15x15x15 integration cells, GMRES\n");
    printf("  refine15 regular 15x15x15 nodes, 15x15x15 integration cells, GMRES\n");
    printf("  plane15 solve refine15 and export x=5.33 plane CSV\n");
    printf("  all     run sanity, target, refine13, and refine15 (default)\n");
}

static int parse_args(int argc, char **argv, CubeSparseSelection *selection)
{
    *selection = CUBE_SPARSE_ALL;

    if (argc == 1) {
        return 0;
    }

    if (argc == 2 &&
        (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        print_usage(argv[0]);
        return 1;
    }

    if (argc == 3 && strcmp(argv[1], "--case") == 0) {
        if (strcmp(argv[2], "sanity") == 0) {
            *selection = CUBE_SPARSE_SANITY;
            return 0;
        }
        if (strcmp(argv[2], "target") == 0) {
            *selection = CUBE_SPARSE_TARGET;
            return 0;
        }
        if (strcmp(argv[2], "refine13") == 0) {
            *selection = CUBE_SPARSE_REFINE13;
            return 0;
        }
        if (strcmp(argv[2], "refine15") == 0) {
            *selection = CUBE_SPARSE_REFINE15;
            return 0;
        }
        if (strcmp(argv[2], "plane15") == 0) {
            *selection = CUBE_SPARSE_PLANE15;
            return 0;
        }
        if (strcmp(argv[2], "all") == 0) {
            *selection = CUBE_SPARSE_ALL;
            return 0;
        }
    }

    fprintf(stderr, "Invalid command-line arguments.\n");
    print_usage(argv[0]);
    return 2;
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
                    int run_dense_cmp,
                    const PlaneExportConfig *plane_export)
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
    CubeSparseReport report;

    double dense_rel_diff = -1.0;  /* < 0 means comparison was not run */
    int rc = 1;

    const int expected_nodes       = cube_regular_node_count(nx, ny, nz);
    const int expected_constraints = cube_dirichlet_point_count(nx, ny, nz);

    report_init(&report, label);

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
    report.nodes          = node_count;
    report.constraints    = point_count;
    report.augmented_size = total_size;

    /* -------------------------------------------------- MLS diagnostic */
    /*
     * Evaluate connectivity at the nx_cells × ny_cells × nz_cells Gauss cell
     * centres. mls_connectivity_stats_uniform_grid with [xmin=0, xmax=L, n]
     * places points at (ix + 0.5) * L / n, which are exactly the cell centres.
     *
     * The MLS linear basis p^T = [1, x, y, z] has 4 terms. For the moment
     * matrix A(x) to be invertible, each evaluation point x must have at
     * least 4 active nodes in its support. This is a critical condition.
     *
     * - active_nodes < 4 (support_lt_4 > 0): fatal, A(x) is singular.
     * - active_nodes < 8 (support_lt_8 > 0): warning, A(x) may be ill-
     *   conditioned, leading to numerical instability.
     */
    mls_connectivity_stats_uniform_grid(nodes, node_count,
                                        0.0, L, nx_cells,
                                        0.0, L, ny_cells,
                                        0.0, L, nz_cells,
                                        &diag);
    report.support_min   = diag.min_nodes;
    report.support_mean  = diag.mean_nodes;
    report.support_max   = diag.max_nodes;
    report.support_lt_4  = diag.n_invalid;
    report.support_lt_8  = diag.n_invalid + diag.n_suspect;
    report.mls_failures  = diag.n_moment_fail;
    report.max_cond_estimate = diag.max_cond;

    if (report.support_lt_4 > 0) {
        fprintf(stderr, "[%s] STOP: support_lt_4 = %d (> 0)\n",
                label, report.support_lt_4);
        emit_required_report(&report);
        goto done;
    }

    if (report.mls_failures > 0) {
        fprintf(stderr, "[%s] STOP: mls_failures = %d (> 0)\n",
                label, report.mls_failures);
        emit_required_report(&report);
        goto done;
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
        int k_cap = node_count * 100 + 16;
        if (sparse_coo_create(&K_coo, node_count, node_count, k_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] K_coo alloc failed\n", label);
            goto done;
        }
        for (int r = 0; r < node_count; ++r) {
            for (int c = 0; c < node_count; ++c) {
                double v = K_dense[r * node_count + c];
                if (v != 0.0) {
                    if (sparse_coo_add(&K_coo, r, c, v) != SPARSE_OK) {
                        fprintf(stderr, "[%s] K_coo append failed\n", label);
                        goto done;
                    }
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
        /*
         * The augmented system [K G^T; G 0] arises from imposing Dirichlet
         * boundary conditions via Lagrange multipliers. This increases the
         * system size but avoids empirical penalty parameters. The resulting
         * matrix is symmetric but not positive-definite, requiring a solver
         * like GMRES.
         */
        int aug_cap = K_coo.count + point_count * 60 + 16;
        if (sparse_coo_create(&A_aug_coo, total_size, total_size, aug_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] A_aug_coo alloc failed\n", label);
            goto done;
        }

        /* K block (rows and cols 0..node_count-1) */
        for (int e = 0; e < K_coo.count; ++e) {
            if (sparse_coo_add(&A_aug_coo,
                               K_coo.row[e], K_coo.col[e], K_coo.val[e])
                != SPARSE_OK) {
                fprintf(stderr, "[%s] A_aug K-block append failed\n", label);
                goto done;
            }
        }

        /* G and G^T blocks */
        for (int ci = 0; ci < point_count; ++ci) {
            for (int ni = 0; ni < node_count; ++ni) {
                double val = G_dense[ci * node_count + ni];
                if (val != 0.0) {
                    /* G^T: A_aug[ni, node_count+ci] */
                    if (sparse_coo_add(&A_aug_coo, ni, node_count + ci, val)
                        != SPARSE_OK) {
                        fprintf(stderr, "[%s] A_aug G^T append failed\n", label);
                        goto done;
                    }
                    /* G:   A_aug[node_count+ci, ni] */
                    if (sparse_coo_add(&A_aug_coo, node_count + ci, ni, val)
                        != SPARSE_OK) {
                        fprintf(stderr, "[%s] A_aug G append failed\n", label);
                        goto done;
                    }
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
    if (sparse_coo_to_csr(&A_aug_coo, &A_aug_csr) != SPARSE_OK) {
        fprintf(stderr, "[%s] A_aug CSR conversion failed\n", label);
        goto done;
    }

    t_asm = elapsed_s(t_start);
    report.K_nnz          = K_coo.count;
    report.A_aug_nnz      = A_aug_csr.nnz;
    report.assembly_time_s = t_asm;

    /* ------------------------------------------------------------ GMRES */

    x_sol = calloc((size_t)total_size, sizeof(double));
    if (!x_sol) { fprintf(stderr, "[%s] x_sol alloc failed\n", label); goto done; }

    t_start = clock();

    /*
     * GMRES is used because the augmented system with Lagrange multipliers
     * is symmetric but indefinite (due to the zero block), so standard
     * conjugate gradient methods cannot be used.
     */
    if (gmres_solve(&A_aug_csr, b_aug, x_sol,
                    gmres_tol, gmres_max_iter, gmres_restart,
                    &gmres_res) != GMRES_OK) {
        fprintf(stderr, "[%s] gmres_solve returned error\n", label);
        goto done;
    }

    t_sol = elapsed_s(t_start);
    report.gmres_converged = gmres_res.converged;
    report.gmres_iterations = gmres_res.iterations;
    report.residual_initial = gmres_res.residual_init;
    report.residual_final   = gmres_res.residual_final;
    report.solve_time_s     = t_sol;

    if (!gmres_res.converged) {
        fprintf(stderr, "[%s] STOP: GMRES did not converge\n", label);
        emit_required_report(&report);
        goto done;
    }

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

    if (plane_export != NULL && plane_export->enabled) {
        if (export_plane_csv(plane_export, nodes, node_count, x_sol,
                             shape_buf, L, V0) != 0) {
            goto done;
        }
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
        report.max_abs_error = max_abs_err;
        report.relative_error_global = rel_global;
        report.relative_error_interior = rel_intrnl;
        report.mls_failures += mls_eval_fail;

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
        /*
         * The max absolute error is expected to be ~V0 due to the singularity
         * at the top edges (e.g., at (0, y, L), V=0 from the side face and
         * V=V0 from the top face). The analytical solution resolves this to 0,
         * while the Lagrange multiplier formulation gives priority to the top
         * face, resulting in V=V0. The interior error is the relevant metric.
         */
        printf("  max abs error:          %.6e\n", max_abs_err);
        printf("  rel error (global):     %.6e\n", rel_global);
        printf("  rel error (interior):   %.6e\n", rel_intrnl);
        if (mls_eval_fail > 0) {
            printf("  MLS eval failures:      %d\n", mls_eval_fail);
        }
        printf("\n");

        if (emit_required_report(&report) != 0) {
            goto done;
        }

        if (report.mls_failures > 0) {
            fprintf(stderr, "[%s] STOP: mls_failures = %d (> 0)\n",
                    label, report.mls_failures);
            goto done;
        }
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
 * Runs two cases by default:
 *   1. 5x5x5 sanity check  — compared against the dense solver.
 *   2. 11x11x11 target     — GMRES only (dense solve is O(n^3) too slow).
 *   3. 13x13x13 refinement — GMRES only, larger restart.
 *   4. 15x15x15 refinement — GMRES only, larger restart.
 * The plane15 case reuses the 15x15x15 setup and exports x = 5.33.
 */
int main(int argc, char **argv)
{
    CubeSparseSelection selection;
    int status = 0;
    int parse_status = parse_args(argc, argv, &selection);

    if (parse_status == 1) {
        return 0;
    }
    if (parse_status != 0) {
        return 1;
    }

    if (write_report_csv_header() != 0) {
        return 1;
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_SANITY) {
        status |= run_case(
            "5x5x5 sanity check (dense comparison enabled)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/5, 5, 5,
            /*cells=*/5, 5, 5,
            /*gmres_tol=*/1e-9,
            /*restart=*/30,
            /*max_iter=*/2000,
            /*sample_n=*/11,
            /*dense_cmp=*/1,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_TARGET) {
        status |= run_case(
            "11x11x11 target (15x15x15 integration cells, GMRES only)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/11, 11, 11,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/200,
            /*max_iter=*/10000,
            /*sample_n=*/11,
            /*dense_cmp=*/0,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_REFINE13) {
        status |= run_case(
            "13x13x13 refine13 (15x15x15 integration cells, GMRES only)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/13, 13, 13,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            /*dense_cmp=*/0,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_REFINE15) {
        status |= run_case(
            "15x15x15 refine15 (15x15x15 integration cells, GMRES only)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/15, 15, 15,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            /*dense_cmp=*/0,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_PLANE15) {
        const PlaneExportConfig plane_export = {
            1,
            5.33,
            101,
            101,
            PLANE15_CSV_PATH
        };
        status |= run_case(
            "15x15x15 plane15 (export x=5.33 plane)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/15, 15, 15,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            /*dense_cmp=*/0,
            &plane_export);
    }

    return status;
}
