#include "analytical.h"
#include "cube_problem.h"
#include "dense_solver.h"
#include "dirichlet.h"
#include "efg/gmres.h"
#include "error_norm.h"
#include "gauss.h"
#include "global_stiffness.h"
#include "lagrange_system.h"
#include "lapack_dense_solver.h"
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
#define NONUNIFORM_CSV_PATH "data/output/cube_plane_x_5_33_nonuniform.csv"
#define NONUNIFORM_REFINE_CSV_PATH "data/output/cube_plane_x_5_33_nonuniform_refine.csv"
#define INTEGRATION_NORM_STUDY_CSV_PATH "data/output/integration_order_norm_sensitivity_refine15.csv"
#define INTEGRATION_SOLUTION_STUDY_CSV_PATH "data/output/integration_order_solution_sensitivity_refine15.csv"
#define ERROR_REGION_STUDY_CSV_PATH "data/output/error_region_integral_refine15.csv"
#define ERROR_REGION_COUNT 8

typedef struct CubeSparseReport {
    const char *label;
    const char *solver;
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
    int lapack_available;
    int lapack_info;
    double support_mean;
    double max_cond_estimate;
    double residual_initial;
    double residual_final;
    double lapack_conversion_time_s;
    double lapack_solve_time_s;
    double lapack_residual_final;
    double gmres_lapack_rel_diff;
    double relative_error_global;
    double relative_error_discrete_global;
    double relative_error_interior;
    double plane_relative_error;
    double plane_max_abs_error;
    double plane_mean_abs_error;
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
    CUBE_SPARSE_PLANE15,
    CUBE_SPARSE_NONUNIFORM,
    CUBE_SPARSE_NONUNIFORM_REFINE
} CubeSparseSelection;

typedef enum CubeSparseSolver {
    CUBE_SOLVER_GMRES,
    CUBE_SOLVER_LAPACK_DENSE
} CubeSparseSolver;

typedef enum CubeIntegrationStudy {
    CUBE_INTEGRATION_STUDY_NONE,
    CUBE_INTEGRATION_STUDY_NORM,
    CUBE_INTEGRATION_STUDY_ASSEMBLY
} CubeIntegrationStudy;

typedef struct PlaneExportConfig {
    int enabled;
    double x;
    int y_samples;
    int z_samples;
    const char *path;
} PlaneExportConfig;

typedef struct PlaneExportMetrics {
    int valid_points;
    int mls_failures;
    double relative_error;
    double max_abs_error;
    double mean_abs_error;
} PlaneExportMetrics;

typedef struct CubeAnalyticalContext {
    double L;
    double V0;
    int terms;
} CubeAnalyticalContext;

typedef struct DiscreteErrorMetrics {
    double relative_global;
    double relative_interior;
    double max_abs_error;
    int mls_failures;
    int valid_global;
    int valid_interior;
} DiscreteErrorMetrics;

typedef struct ErrorRegionAccumulator {
    const char *name;
    int points_or_cells_used;
    int mls_failures;
    double numerator_integral;
    double denominator_integral;
    double abs_error_integral;
    double volume_integral;
    double relative_error_integral;
    double mean_abs_error_integral;
    double max_abs_error_sampled;
} ErrorRegionAccumulator;

typedef struct StudyRegularSolution {
    int node_count;
    int point_count;
    int total_size;
    int coo_count_aug;
    int assembly_order;
    int mls_failures_assembly;
    double assembly_time_s;
    double solve_time_s;
    Node3D *nodes;
    DirichletPoint *dp;
    double *K_dense;
    double *G_dense;
    double *F;
    double *q_vec;
    double *b_aug;
    double *x_sol;
    double *jacobi_precond;
    SparseCOO K_coo;
    SparseCOO A_aug_coo;
    SparseCSR A_aug_csr;
    GmresResult gmres_res;
    MlsConnectivityStats diag;
} StudyRegularSolution;

static double elapsed_s(clock_t start)
{
    return (double)(clock() - start) / (double)CLOCKS_PER_SEC;
}

static const char *solver_name(CubeSparseSolver solver)
{
    switch (solver) {
    case CUBE_SOLVER_GMRES:
        return "gmres";
    case CUBE_SOLVER_LAPACK_DENSE:
        return "lapack-dense";
    }
    return "unknown";
}

static double vector_norm2(const double *x, int n)
{
    double sum = 0.0;

    if (x == NULL || n <= 0) {
        return NAN;
    }

    for (int i = 0; i < n; ++i) {
        sum += x[i] * x[i];
    }

    return sqrt(sum);
}

static double relative_solution_difference(const double *candidate,
                                           const double *reference,
                                           int n)
{
    double diff2 = 0.0;
    double ref2 = 0.0;

    if (candidate == NULL || reference == NULL || n <= 0) {
        return NAN;
    }

    for (int i = 0; i < n; ++i) {
        const double d = candidate[i] - reference[i];
        diff2 += d * d;
        ref2 += reference[i] * reference[i];
    }

    return (ref2 > 0.0) ? sqrt(diff2 / ref2) : sqrt(diff2);
}

static double cube_exact_potential_callback(double x,
                                            double y,
                                            double z,
                                            void *context)
{
    const CubeAnalyticalContext *ctx = (const CubeAnalyticalContext *)context;

    return analytical_potential_cube(x, y, z,
                                     ctx->L, ctx->V0, ctx->terms);
}

static void report_init(CubeSparseReport *report,
                        const char *label,
                        CubeSparseSolver solver)
{
    report->label = label;
    report->solver = solver_name(solver);
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
#ifdef EFG_HAVE_LAPACK
    report->lapack_available = 1;
#else
    report->lapack_available = 0;
#endif
    report->lapack_info = -1;
    report->support_mean = NAN;
    report->max_cond_estimate = NAN;
    report->residual_initial = NAN;
    report->residual_final = NAN;
    report->lapack_conversion_time_s = NAN;
    report->lapack_solve_time_s = NAN;
    report->lapack_residual_final = NAN;
    report->gmres_lapack_rel_diff = NAN;
    report->relative_error_global = NAN;
    report->relative_error_discrete_global = NAN;
    report->relative_error_interior = NAN;
    report->plane_relative_error = NAN;
    report->plane_max_abs_error = NAN;
    report->plane_mean_abs_error = NAN;
    report->max_abs_error = NAN;
    report->assembly_time_s = NAN;
    report->solve_time_s = NAN;
}

static void print_required_report(const CubeSparseReport *report)
{
    printf("--- Required report ---\n");
    printf("label: %s\n", report->label);
    printf("solver: %s\n", report->solver);
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
    printf("lapack_available: %s\n", report->lapack_available ? "YES" : "NO");
    printf("lapack_info: %d\n", report->lapack_info);
    printf("residual_initial: %.6e\n", report->residual_initial);
    printf("residual_final: %.6e\n", report->residual_final);
    printf("lapack_conversion_time_s: %.6f\n", report->lapack_conversion_time_s);
    printf("lapack_solve_time_s: %.6f\n", report->lapack_solve_time_s);
    printf("lapack_residual_final: %.6e\n", report->lapack_residual_final);
    printf("gmres_lapack_rel_diff: %.6e\n", report->gmres_lapack_rel_diff);
    printf("relative_error_global: %.6e\n", report->relative_error_global);
    printf("relative_error_discrete_global: %.6e\n",
           report->relative_error_discrete_global);
    printf("relative_error_interior: %.6e\n", report->relative_error_interior);
    printf("plane_relative_error: %.6e\n", report->plane_relative_error);
    printf("plane_max_abs_error: %.6e\n", report->plane_max_abs_error);
    printf("plane_mean_abs_error: %.6e\n", report->plane_mean_abs_error);
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
            "label,solver,nodes,constraints,augmented_size,K_nnz,A_aug_nnz,"
            "support_min,support_mean,support_max,support_lt_4,"
            "support_lt_8,mls_failures,max_cond_estimate,"
            "gmres_converged,gmres_iterations,lapack_available,lapack_info,"
            "residual_initial,residual_final,lapack_conversion_time_s,"
            "lapack_solve_time_s,lapack_residual_final,"
            "gmres_lapack_rel_diff,relative_error_global,"
            "relative_error_discrete_global,relative_error_interior,"
            "plane_relative_error,"
            "plane_max_abs_error,plane_mean_abs_error,"
            "max_abs_error,assembly_time_s,"
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
            "\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%.17g,%d,%d,%d,%d,%.17g,%s,%d,"
            "%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            report->label,
            report->solver,
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
            report->lapack_available ? "YES" : "NO",
            report->lapack_info,
            report->residual_initial,
            report->residual_final,
            report->lapack_conversion_time_s,
            report->lapack_solve_time_s,
            report->lapack_residual_final,
            report->gmres_lapack_rel_diff,
            report->relative_error_global,
            report->relative_error_discrete_global,
            report->relative_error_interior,
            report->plane_relative_error,
            report->plane_max_abs_error,
            report->plane_mean_abs_error,
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
                            double V0,
                            PlaneExportMetrics *metrics)
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
    /*
     * The largest plane errors remain close to the upper side edges, where
     * the analytical sine series is difficult to approximate with a finite
     * number of terms and the MLS field crosses a steep boundary layer. The
     * visualization script also reports interior-plane metrics, which are more
     * representative of the reproduced contour figure.
     */
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

    if (metrics != NULL) {
        metrics->valid_points = valid_points;
        metrics->mls_failures = mls_failures;
        metrics->max_abs_error = max_abs_error;
        metrics->mean_abs_error =
            (valid_points > 0) ? sum_abs_error / (double)valid_points : NAN;
        metrics->relative_error = (ref2 > 0.0) ? sqrt(err2 / ref2) : NAN;
    }

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
    printf("  %s --case nonuniform\n", program_name);
    printf("  %s --case nonuniform_refine\n", program_name);
    printf("  %s --case refine15 --solver gmres\n", program_name);
    printf("  %s --case refine15 --solver lapack-dense\n", program_name);
    printf("  %s --case refine15 --integration-study norm\n", program_name);
    printf("  %s --case refine15 --integration-study assembly\n", program_name);
    printf("  %s --case refine15 --error-region-study\n", program_name);
    printf("  %s --case all\n", program_name);
    printf("\n");
    printf("Cases:\n");
    printf("  sanity     regular 5x5x5 nodes, 5x5x5 integration cells, dense comparison\n");
    printf("  target     regular 11x11x11 nodes, 15x15x15 integration cells\n");
    printf("  refine13   regular 13x13x13 nodes, 15x15x15 integration cells\n");
    printf("  refine15   regular 15x15x15 nodes, 15x15x15 integration cells\n");
    printf("  plane15    solve refine15 and export x=5.33 plane CSV\n");
    printf("  nonuniform non-uniform cloud (Fig. 2), export x=5.33 plane CSV\n");
    printf("  nonuniform_refine non-uniform cloud base=13 top=15, export plane CSV\n");
    printf("  all        run sanity, target, refine13, and refine15 (default)\n");
    printf("\n");
    printf("Solvers:\n");
    printf("  gmres        restarted GMRES with current Jacobi path (default)\n");
    printf("  lapack-dense optional LAPACK dgesv validation solver\n");
    printf("\n");
    printf("Integration studies:\n");
    printf("  norm      solve refine15 once with assembly order 2, then evaluate Eq.16 with norm orders 1..8\n");
    printf("  assembly  solve refine15 with assembly orders 1..5, then evaluate Eq.16 with norm order 6\n");
    printf("\n");
    printf("Error-region study:\n");
    printf("  --error-region-study  solve refine15 once and evaluate regional Eq.16 integrals with norm order 6\n");
}

static int parse_args(int argc,
                      char **argv,
                      CubeSparseSelection *selection,
                      CubeSparseSolver *solver,
                      CubeIntegrationStudy *integration_study,
                      int *error_region_study)
{
    *selection = CUBE_SPARSE_ALL;
    *solver = CUBE_SOLVER_GMRES;
    *integration_study = CUBE_INTEGRATION_STUDY_NONE;
    *error_region_study = 0;

    if (argc == 1) {
        return 0;
    }

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 1;
        }

        if (strcmp(argv[i], "--case") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --case.\n");
                print_usage(argv[0]);
                return 2;
            }
            ++i;
            if (strcmp(argv[i], "sanity") == 0) {
                *selection = CUBE_SPARSE_SANITY;
                continue;
            }
            if (strcmp(argv[i], "target") == 0) {
                *selection = CUBE_SPARSE_TARGET;
                continue;
            }
            if (strcmp(argv[i], "refine13") == 0) {
                *selection = CUBE_SPARSE_REFINE13;
                continue;
            }
            if (strcmp(argv[i], "refine15") == 0) {
                *selection = CUBE_SPARSE_REFINE15;
                continue;
            }
            if (strcmp(argv[i], "plane15") == 0) {
                *selection = CUBE_SPARSE_PLANE15;
                continue;
            }
            if (strcmp(argv[i], "nonuniform") == 0) {
                *selection = CUBE_SPARSE_NONUNIFORM;
                continue;
            }
            if (strcmp(argv[i], "nonuniform_refine") == 0) {
                *selection = CUBE_SPARSE_NONUNIFORM_REFINE;
                continue;
            }
            if (strcmp(argv[i], "all") == 0) {
                *selection = CUBE_SPARSE_ALL;
                continue;
            }

            fprintf(stderr, "Invalid --case value: %s\n", argv[i]);
            print_usage(argv[0]);
            return 2;
        }

        if (strcmp(argv[i], "--solver") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --solver.\n");
                print_usage(argv[0]);
                return 2;
            }
            ++i;
            if (strcmp(argv[i], "gmres") == 0) {
                *solver = CUBE_SOLVER_GMRES;
                continue;
            }
            if (strcmp(argv[i], "lapack-dense") == 0) {
                *solver = CUBE_SOLVER_LAPACK_DENSE;
                continue;
            }

            fprintf(stderr, "Invalid --solver value: %s\n", argv[i]);
            print_usage(argv[0]);
            return 2;
        }

        if (strcmp(argv[i], "--integration-study") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --integration-study.\n");
                print_usage(argv[0]);
                return 2;
            }
            ++i;
            if (strcmp(argv[i], "norm") == 0) {
                *integration_study = CUBE_INTEGRATION_STUDY_NORM;
                continue;
            }
            if (strcmp(argv[i], "assembly") == 0) {
                *integration_study = CUBE_INTEGRATION_STUDY_ASSEMBLY;
                continue;
            }

            fprintf(stderr, "Invalid --integration-study value: %s\n", argv[i]);
            print_usage(argv[0]);
            return 2;
        }

        if (strcmp(argv[i], "--error-region-study") == 0) {
            *error_region_study = 1;
            continue;
        }

        fprintf(stderr, "Invalid command-line argument: %s\n", argv[i]);
        print_usage(argv[0]);
        return 2;
    }

    return 0;
}

/*
 * Build left Jacobi preconditioner: diag_inv[i] = 1/A_ii.
 * Rows with zero diagonal (the Lagrange multiplier block) get diag_inv[i]=1.
 */
static void build_jacobi_precond(const SparseCSR *A, double *diag_inv)
{
    for (int i = 0; i < A->nrows; ++i) {
        double d = 0.0;
        for (int p = A->row_ptr[i]; p < A->row_ptr[i + 1]; ++p) {
            if (A->col_ind[p] == i) { d = A->values[p]; break; }
        }
        diag_inv[i] = (fabs(d) > 1e-15) ? 1.0 / d : 1.0;
    }
}

static void study_regular_solution_init(StudyRegularSolution *solution)
{
    if (solution == NULL) {
        return;
    }

    memset(solution, 0, sizeof(*solution));
    solution->K_coo = (SparseCOO){NULL, NULL, NULL, 0, 0, 0, 0};
    solution->A_aug_coo = (SparseCOO){NULL, NULL, NULL, 0, 0, 0, 0};
    solution->A_aug_csr = (SparseCSR){NULL, NULL, NULL, 0, 0, 0};
    solution->gmres_res = (GmresResult){0, 0, NAN, NAN};
    solution->diag = (MlsConnectivityStats){0, 0, 0, 0, 0, 0,
                                            NAN, NAN, NAN};
}

static void study_regular_solution_destroy(StudyRegularSolution *solution)
{
    if (solution == NULL) {
        return;
    }

    free(solution->nodes);
    free(solution->dp);
    free(solution->K_dense);
    free(solution->G_dense);
    free(solution->F);
    free(solution->q_vec);
    free(solution->b_aug);
    free(solution->x_sol);
    free(solution->jacobi_precond);
    sparse_coo_destroy(&solution->K_coo);
    sparse_coo_destroy(&solution->A_aug_coo);
    sparse_csr_destroy(&solution->A_aug_csr);
    study_regular_solution_init(solution);
}

static int dense_stiffness_to_coo(const double *K_dense,
                                  int node_count,
                                  SparseCOO *K_coo)
{
    int k_cap;

    if (K_dense == NULL || K_coo == NULL || node_count <= 0) {
        return 1;
    }

    k_cap = node_count * 100 + 16;
    if (sparse_coo_create(K_coo, node_count, node_count, k_cap)
        != SPARSE_OK) {
        return 1;
    }

    for (int r = 0; r < node_count; ++r) {
        for (int c = 0; c < node_count; ++c) {
            const double v = K_dense[r * node_count + c];
            if (v != 0.0) {
                if (sparse_coo_add(K_coo, r, c, v) != SPARSE_OK) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

static int assemble_augmented_from_dense_constraints(const SparseCOO *K_coo,
                                                     const double *G_dense,
                                                     int node_count,
                                                     int point_count,
                                                     SparseCOO *A_aug_coo)
{
    const int total_size = node_count + point_count;
    int aug_cap;

    if (K_coo == NULL || G_dense == NULL || A_aug_coo == NULL ||
        node_count <= 0 || point_count <= 0) {
        return 1;
    }

    aug_cap = K_coo->count + point_count * 60 + 16;
    if (sparse_coo_create(A_aug_coo, total_size, total_size, aug_cap)
        != SPARSE_OK) {
        return 1;
    }

    for (int e = 0; e < K_coo->count; ++e) {
        if (sparse_coo_add(A_aug_coo,
                           K_coo->row[e], K_coo->col[e], K_coo->val[e])
            != SPARSE_OK) {
            return 1;
        }
    }

    for (int ci = 0; ci < point_count; ++ci) {
        for (int ni = 0; ni < node_count; ++ni) {
            const double val = G_dense[ci * node_count + ni];
            if (val != 0.0) {
                if (sparse_coo_add(A_aug_coo, ni, node_count + ci, val)
                    != SPARSE_OK) {
                    return 1;
                }
                if (sparse_coo_add(A_aug_coo, node_count + ci, ni, val)
                    != SPARSE_OK) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

static int compute_discrete_error_metrics(const Node3D *nodes,
                                          int node_count,
                                          const double *solution,
                                          double L,
                                          double V0,
                                          int sample_n,
                                          DiscreteErrorMetrics *metrics)
{
    MlsShapeValue *shape_buf = NULL;
    double global_err2 = 0.0;
    double global_ref2 = 0.0;
    double interior_err2 = 0.0;
    double interior_ref2 = 0.0;
    double max_abs_error = 0.0;
    int mls_failures = 0;
    int valid_global = 0;
    int valid_interior = 0;

    if (nodes == NULL || solution == NULL || metrics == NULL ||
        node_count <= 0 || L <= 0.0 || sample_n < 2) {
        return 1;
    }

    shape_buf = malloc((size_t)node_count * sizeof(shape_buf[0]));
    if (shape_buf == NULL) {
        return 1;
    }

    for (int ix = 0; ix < sample_n; ++ix) {
        const double sx = (double)ix / (double)(sample_n - 1) * L;

        for (int iy = 0; iy < sample_n; ++iy) {
            const double sy = (double)iy / (double)(sample_n - 1) * L;

            for (int iz = 0; iz < sample_n; ++iz) {
                const double sz = (double)iz / (double)(sample_n - 1) * L;
                const double v_exact =
                    analytical_potential_cube(sx, sy, sz, L, V0,
                                              ANALYTICAL_TERMS);
                double v_num = 0.0;
                int value_count = 0;
                const int mls_status =
                    mls_linear3d_shape_functions(nodes, node_count,
                                                 sx, sy, sz,
                                                 shape_buf, node_count,
                                                 &value_count);

                if (mls_status != MLS_OK) {
                    ++mls_failures;
                    continue;
                }

                for (int k = 0; k < value_count; ++k) {
                    v_num += shape_buf[k].phi
                             * solution[shape_buf[k].node_index];
                }

                {
                    const double diff = v_num - v_exact;
                    const double abs_error = fabs(diff);
                    const int interior =
                        (ix > 0 && ix < sample_n - 1 &&
                         iy > 0 && iy < sample_n - 1 &&
                         iz > 0 && iz < sample_n - 1);

                    if (abs_error > max_abs_error) {
                        max_abs_error = abs_error;
                    }
                    global_err2 += diff * diff;
                    global_ref2 += v_exact * v_exact;
                    ++valid_global;

                    if (interior) {
                        interior_err2 += diff * diff;
                        interior_ref2 += v_exact * v_exact;
                        ++valid_interior;
                    }
                }
            }
        }
    }

    free(shape_buf);

    metrics->relative_global =
        (global_ref2 > 0.0) ? sqrt(global_err2 / global_ref2) : NAN;
    metrics->relative_interior =
        (interior_ref2 > 0.0) ? sqrt(interior_err2 / interior_ref2) : NAN;
    metrics->max_abs_error = max_abs_error;
    metrics->mls_failures = mls_failures;
    metrics->valid_global = valid_global;
    metrics->valid_interior = valid_interior;

    return 0;
}

static int solve_refine15_for_integration_study(int assembly_order,
                                                StudyRegularSolution *solution)
{
    const double L = 10.0;
    const double V0 = 10.0;
    const int nx = 15;
    const int ny = 15;
    const int nz = 15;
    const int nx_cells = 15;
    const int ny_cells = 15;
    const int nz_cells = 15;
    const int gmres_restart = 300;
    const int gmres_max_iter = 20000;
    const double gmres_tol = 1.0e-9;
    const int expected_nodes = cube_regular_node_count(nx, ny, nz);
    const int expected_constraints = cube_dirichlet_point_count(nx, ny, nz);
    clock_t t_start;

    if (solution == NULL ||
        assembly_order < GAUSS_LEGENDRE_MIN_ORDER ||
        assembly_order > GAUSS_LEGENDRE_MAX_ORDER ||
        expected_nodes <= 0 || expected_constraints <= 0) {
        return 1;
    }

    study_regular_solution_init(solution);
    solution->assembly_order = assembly_order;

    solution->nodes = malloc((size_t)expected_nodes * sizeof(solution->nodes[0]));
    solution->dp = malloc((size_t)expected_constraints * sizeof(solution->dp[0]));
    solution->K_dense = calloc((size_t)expected_nodes * (size_t)expected_nodes,
                               sizeof(solution->K_dense[0]));
    solution->G_dense =
        calloc((size_t)expected_constraints * (size_t)expected_nodes,
               sizeof(solution->G_dense[0]));
    solution->F = calloc((size_t)expected_nodes, sizeof(solution->F[0]));
    solution->q_vec =
        calloc((size_t)expected_constraints, sizeof(solution->q_vec[0]));

    if (solution->nodes == NULL || solution->dp == NULL ||
        solution->K_dense == NULL || solution->G_dense == NULL ||
        solution->F == NULL || solution->q_vec == NULL) {
        fprintf(stderr, "[integration-study] allocation failed\n");
        return 1;
    }

    if (cube_generate_regular_nodes(L, nx, ny, nz,
                                    solution->nodes, expected_nodes,
                                    &solution->node_count)
        != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[integration-study] node generation failed\n");
        return 1;
    }

    if (support_assign_article_default(solution->nodes,
                                       solution->node_count) != SUPPORT_OK) {
        fprintf(stderr, "[integration-study] support assignment failed\n");
        return 1;
    }

    if (cube_generate_dirichlet_points_from_nodes(L, V0,
                                                  solution->nodes,
                                                  solution->node_count,
                                                  solution->dp,
                                                  expected_constraints,
                                                  &solution->point_count)
        != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[integration-study] Dirichlet generation failed\n");
        return 1;
    }

    solution->total_size = solution->node_count + solution->point_count;

    mls_connectivity_stats_uniform_grid(solution->nodes,
                                        solution->node_count,
                                        0.0, L, nx_cells,
                                        0.0, L, ny_cells,
                                        0.0, L, nz_cells,
                                        &solution->diag);
    solution->mls_failures_assembly = solution->diag.n_moment_fail;

    if (solution->diag.n_invalid > 0 ||
        solution->diag.n_moment_fail > 0) {
        fprintf(stderr,
                "[integration-study] support/MLS diagnostic failed "
                "(active<4=%d, MLS failures=%d)\n",
                solution->diag.n_invalid, solution->diag.n_moment_fail);
        return 1;
    }

    t_start = clock();

    if (global_stiffness_assemble_dense_with_order(
            solution->nodes, solution->node_count,
            0.0, L, 0.0, L, 0.0, L,
            nx_cells, ny_cells, nz_cells,
            assembly_order,
            solution->K_dense) != GLOBAL_STIFFNESS_OK) {
        fprintf(stderr,
                "[integration-study] K assembly failed for order %d\n",
                assembly_order);
        return 1;
    }

    if (dirichlet_assemble_constraints_dense(solution->nodes,
                                             solution->node_count,
                                             solution->dp,
                                             solution->point_count,
                                             solution->G_dense,
                                             solution->q_vec)
        != DIRICHLET_OK) {
        fprintf(stderr, "[integration-study] G assembly failed\n");
        return 1;
    }

    if (dense_stiffness_to_coo(solution->K_dense,
                               solution->node_count,
                               &solution->K_coo) != 0) {
        fprintf(stderr, "[integration-study] K COO conversion failed\n");
        return 1;
    }

    if (assemble_augmented_from_dense_constraints(&solution->K_coo,
                                                  solution->G_dense,
                                                  solution->node_count,
                                                  solution->point_count,
                                                  &solution->A_aug_coo)
        != 0) {
        fprintf(stderr, "[integration-study] augmented COO assembly failed\n");
        return 1;
    }

    solution->b_aug = calloc((size_t)solution->total_size,
                             sizeof(solution->b_aug[0]));
    if (solution->b_aug == NULL) {
        fprintf(stderr, "[integration-study] b_aug allocation failed\n");
        return 1;
    }
    for (int ci = 0; ci < solution->point_count; ++ci) {
        solution->b_aug[solution->node_count + ci] = solution->q_vec[ci];
    }

    solution->coo_count_aug = solution->A_aug_coo.count;
    if (sparse_coo_to_csr(&solution->A_aug_coo,
                          &solution->A_aug_csr) != SPARSE_OK) {
        fprintf(stderr, "[integration-study] augmented CSR conversion failed\n");
        return 1;
    }

    solution->assembly_time_s = elapsed_s(t_start);

    solution->x_sol = calloc((size_t)solution->total_size,
                             sizeof(solution->x_sol[0]));
    solution->jacobi_precond =
        malloc((size_t)solution->total_size
               * sizeof(solution->jacobi_precond[0]));
    if (solution->x_sol == NULL || solution->jacobi_precond == NULL) {
        fprintf(stderr, "[integration-study] solver allocation failed\n");
        return 1;
    }

    build_jacobi_precond(&solution->A_aug_csr, solution->jacobi_precond);

    t_start = clock();
    if (gmres_solve(&solution->A_aug_csr,
                    solution->b_aug,
                    solution->x_sol,
                    gmres_tol,
                    gmres_max_iter,
                    gmres_restart,
                    solution->jacobi_precond,
                    &solution->gmres_res) != GMRES_OK) {
        fprintf(stderr,
                "[integration-study] GMRES returned error for assembly order %d\n",
                assembly_order);
        return 1;
    }
    solution->solve_time_s = elapsed_s(t_start);

    if (!solution->gmres_res.converged) {
        fprintf(stderr,
                "[integration-study] WARNING: GMRES did not converge for "
                "assembly order %d (iters=%d residual=%.6e)\n",
                assembly_order,
                solution->gmres_res.iterations,
                solution->gmres_res.residual_final);
    }

    return 0;
}

static int run_integration_study_norm_refine15(void)
{
    const double L = 10.0;
    const double V0 = 10.0;
    const int nx_cells = 15;
    const int ny_cells = 15;
    const int nz_cells = 15;
    const int assembly_order = 2;
    StudyRegularSolution solution;
    CubeAnalyticalContext exact_ctx = {L, V0, ANALYTICAL_TERMS};
    FILE *file = NULL;
    int rc = 1;

    study_regular_solution_init(&solution);

    printf("\n=== Integration-order norm study: refine15 ===\n");
    printf("assembly order: %d\n", assembly_order);
    printf("norm orders:    1..8\n\n");

    if (solve_refine15_for_integration_study(assembly_order, &solution) != 0) {
        goto done;
    }

    file = fopen(INTEGRATION_NORM_STUDY_CSV_PATH, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open CSV: %s\n",
                INTEGRATION_NORM_STUDY_CSV_PATH);
        goto done;
    }

    fprintf(file,
            "case_name,node_grid,integration_cells,assembly_order,"
            "norm_order,relative_error_global_integral,numerator_integral,"
            "denominator_integral,mls_failures_norm,norm_eval_time_s\n");

    for (int norm_order = 1; norm_order <= 8; ++norm_order) {
        RelativeErrorDomainResult domain_error = {NAN, 0.0, 0.0, 0, 0};
        clock_t t_start = clock();
        const int status =
            relative_error_norm_domain_integral_with_order(
                solution.nodes, solution.node_count, solution.x_sol,
                0.0, L, 0.0, L, 0.0, L,
                nx_cells, ny_cells, nz_cells,
                norm_order,
                cube_exact_potential_callback,
                &exact_ctx,
                &domain_error);
        const double norm_time = elapsed_s(t_start);

        fprintf(file,
                "\"refine15\",\"15x15x15\",\"15x15x15\",%d,%d,"
                "%.17g,%.17g,%.17g,%d,%.17g\n",
                assembly_order,
                norm_order,
                domain_error.error,
                domain_error.numerator_integral,
                domain_error.denominator_integral,
                domain_error.mls_failures,
                norm_time);

        printf("norm_order=%d rel_error=%.6e numerator=%.6e "
               "denominator=%.6e mls_failures=%d time=%.3f s\n",
               norm_order,
               domain_error.error,
               domain_error.numerator_integral,
               domain_error.denominator_integral,
               domain_error.mls_failures,
               norm_time);

        if (status != RELATIVE_ERROR_NORM_OK) {
            fprintf(stderr,
                    "STOP: Eq.16 norm study failed at norm_order=%d "
                    "(status=%d, mls_failures=%d)\n",
                    norm_order, status, domain_error.mls_failures);
            goto done;
        }
    }

    if (fclose(file) != 0) {
        file = NULL;
        fprintf(stderr, "Could not close CSV: %s\n",
                INTEGRATION_NORM_STUDY_CSV_PATH);
        goto done;
    }
    file = NULL;

    printf("\nCSV written: %s\n\n", INTEGRATION_NORM_STUDY_CSV_PATH);
    rc = 0;

done:
    if (file != NULL) {
        fclose(file);
    }
    study_regular_solution_destroy(&solution);
    return rc;
}

static int run_integration_study_assembly_refine15(void)
{
    const double L = 10.0;
    const double V0 = 10.0;
    const int nx_cells = 15;
    const int ny_cells = 15;
    const int nz_cells = 15;
    const int norm_order_ref = 6;
    const int sample_n = 11;
    CubeAnalyticalContext exact_ctx = {L, V0, ANALYTICAL_TERMS};
    FILE *file = NULL;
    int rc = 1;

    printf("\n=== Integration-order assembly study: refine15 ===\n");
    printf("assembly orders: 1..5\n");
    printf("norm order ref:  %d\n\n", norm_order_ref);

    file = fopen(INTEGRATION_SOLUTION_STUDY_CSV_PATH, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open CSV: %s\n",
                INTEGRATION_SOLUTION_STUDY_CSV_PATH);
        return 1;
    }

    fprintf(file,
            "case_name,node_grid,integration_cells,assembly_order,"
            "norm_order_ref,nodes,constraints,augmented_size,K_nnz,"
            "A_aug_nnz,support_min,support_mean,support_max,support_lt_4,"
            "support_lt_8,mls_failures_assembly,mls_failures_norm,"
            "gmres_converged,gmres_iterations,residual_initial,"
            "residual_final,relative_residual,"
            "relative_error_global_integral,relative_error_discrete_global,"
            "relative_error_interior,max_abs_error,assembly_time_s,"
            "solve_time_s,norm_eval_time_s\n");

    for (int assembly_order = 1; assembly_order <= 5; ++assembly_order) {
        StudyRegularSolution solution;
        RelativeErrorDomainResult domain_error = {NAN, 0.0, 0.0, 0, 0};
        DiscreteErrorMetrics discrete_metrics = {NAN, NAN, NAN, 0, 0, 0};
        double norm_time = NAN;
        double relative_residual = NAN;
        int norm_status = RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
        int solve_status;

        study_regular_solution_init(&solution);

        printf("assembly_order=%d: assembling and solving...\n",
               assembly_order);
        fflush(stdout);

        solve_status =
            solve_refine15_for_integration_study(assembly_order, &solution);

        if (solve_status == 0) {
            clock_t t_start = clock();
            norm_status =
                relative_error_norm_domain_integral_with_order(
                    solution.nodes, solution.node_count, solution.x_sol,
                    0.0, L, 0.0, L, 0.0, L,
                    nx_cells, ny_cells, nz_cells,
                    norm_order_ref,
                    cube_exact_potential_callback,
                    &exact_ctx,
                    &domain_error);
            norm_time = elapsed_s(t_start);

            if (compute_discrete_error_metrics(solution.nodes,
                                               solution.node_count,
                                               solution.x_sol,
                                               L, V0, sample_n,
                                               &discrete_metrics) != 0) {
                fprintf(stderr,
                        "[integration-study] discrete metrics failed for "
                        "assembly_order=%d\n",
                        assembly_order);
                study_regular_solution_destroy(&solution);
                goto done;
            }

            relative_residual =
                (solution.gmres_res.residual_init > 0.0)
                ? solution.gmres_res.residual_final
                  / solution.gmres_res.residual_init
                : NAN;
        }

        fprintf(file,
                "\"refine15\",\"15x15x15\",\"15x15x15\",%d,%d,"
                "%d,%d,%d,%d,%d,%d,%.17g,%d,%d,%d,%d,%d,%s,%d,"
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%.17g\n",
                assembly_order,
                norm_order_ref,
                solution.node_count,
                solution.point_count,
                solution.total_size,
                solution.K_coo.count,
                solution.A_aug_csr.nnz,
                solution.diag.min_nodes,
                solution.diag.mean_nodes,
                solution.diag.max_nodes,
                solution.diag.n_invalid,
                solution.diag.n_invalid + solution.diag.n_suspect,
                solution.mls_failures_assembly,
                domain_error.mls_failures,
                solution.gmres_res.converged ? "YES" : "NO",
                solution.gmres_res.iterations,
                solution.gmres_res.residual_init,
                solution.gmres_res.residual_final,
                relative_residual,
                domain_error.error,
                discrete_metrics.relative_global,
                discrete_metrics.relative_interior,
                discrete_metrics.max_abs_error,
                solution.assembly_time_s,
                solution.solve_time_s,
                norm_time);

        printf("assembly_order=%d converged=%s iters=%d "
               "rel_residual=%.3e rel_error_integral=%.6e "
               "discrete_global=%.6e mls_norm=%d assembly=%.3f s "
               "solve=%.3f s norm=%.3f s\n",
               assembly_order,
               solution.gmres_res.converged ? "YES" : "NO",
               solution.gmres_res.iterations,
               relative_residual,
               domain_error.error,
               discrete_metrics.relative_global,
               domain_error.mls_failures,
               solution.assembly_time_s,
               solution.solve_time_s,
               norm_time);

        if (solve_status != 0 ||
            norm_status != RELATIVE_ERROR_NORM_OK ||
            domain_error.mls_failures > 0) {
            fprintf(stderr,
                    "STOP: assembly integration study failed at "
                    "assembly_order=%d (solve_status=%d, norm_status=%d, "
                    "mls_failures_norm=%d)\n",
                    assembly_order,
                    solve_status,
                    norm_status,
                    domain_error.mls_failures);
            study_regular_solution_destroy(&solution);
            goto done;
        }

        study_regular_solution_destroy(&solution);
    }

    if (fclose(file) != 0) {
        file = NULL;
        fprintf(stderr, "Could not close CSV: %s\n",
                INTEGRATION_SOLUTION_STUDY_CSV_PATH);
        goto done;
    }
    file = NULL;

    printf("\nCSV written: %s\n\n", INTEGRATION_SOLUTION_STUDY_CSV_PATH);
    rc = 0;

done:
    if (file != NULL) {
        fclose(file);
    }
    return rc;
}

static void error_regions_init(ErrorRegionAccumulator regions[ERROR_REGION_COUNT])
{
    static const char *names[ERROR_REGION_COUNT] = {
        "full_domain",
        "open_interior",
        "core_delta_0_5",
        "core_delta_1_0",
        "central_box",
        "upper_layer",
        "upper_edge_bands",
        "upper_face_interior"
    };

    for (int i = 0; i < ERROR_REGION_COUNT; ++i) {
        regions[i].name = names[i];
        regions[i].points_or_cells_used = 0;
        regions[i].mls_failures = 0;
        regions[i].numerator_integral = 0.0;
        regions[i].denominator_integral = 0.0;
        regions[i].abs_error_integral = 0.0;
        regions[i].volume_integral = 0.0;
        regions[i].relative_error_integral = NAN;
        regions[i].mean_abs_error_integral = NAN;
        regions[i].max_abs_error_sampled = NAN;
    }
}

static int error_region_contains(int region_index,
                                 double x,
                                 double y,
                                 double z,
                                 double L)
{
    switch (region_index) {
    case 0:
        return x >= 0.0 && x <= L &&
               y >= 0.0 && y <= L &&
               z >= 0.0 && z <= L;
    case 1:
        return x > 0.0 && x < L &&
               y > 0.0 && y < L &&
               z > 0.0 && z < L;
    case 2:
        return x >= 0.5 && x <= 9.5 &&
               y >= 0.5 && y <= 9.5 &&
               z >= 0.5 && z <= 9.5;
    case 3:
        return x >= 1.0 && x <= 9.0 &&
               y >= 1.0 && y <= 9.0 &&
               z >= 1.0 && z <= 9.0;
    case 4:
        return x >= 2.0 && x <= 8.0 &&
               y >= 2.0 && y <= 8.0 &&
               z >= 2.0 && z <= 8.0;
    case 5:
        return z >= 9.0 && z <= 10.0 &&
               x >= 0.0 && x <= L &&
               y >= 0.0 && y <= L;
    case 6:
        return z >= 9.0 && z <= 10.0 &&
               (x <= 1.0 || x >= 9.0 || y <= 1.0 || y >= 9.0);
    case 7:
        return z >= 9.0 && z <= 10.0 &&
               x >= 1.0 && x <= 9.0 &&
               y >= 1.0 && y <= 9.0;
    default:
        return 0;
    }
}

/*
 * Regional Eq. (16) integral for the refine15 error study.
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: relative L2 error restricted to named subregions of
 * the cube. V_EFG(x_gp) is evaluated once per Gauss point by the MLS expansion
 * sum_I Phi_I(x_gp) u_I, then accumulated in every region whose indicator
 * contains that quadrature point.
 */
static int evaluate_error_regions_refine15(
    const Node3D *nodes,
    int node_count,
    const double *solution,
    double L,
    double V0,
    int nx_cells,
    int ny_cells,
    int nz_cells,
    int norm_order,
    ErrorRegionAccumulator regions[ERROR_REGION_COUNT])
{
    MlsShapeValue *shape_values = NULL;
    const double dx = L / (double)nx_cells;
    const double dy = L / (double)ny_cells;
    const double dz = L / (double)nz_cells;

    if (nodes == NULL || solution == NULL || regions == NULL ||
        node_count <= 0 || L <= 0.0 ||
        nx_cells <= 0 || ny_cells <= 0 || nz_cells <= 0 ||
        norm_order < GAUSS_LEGENDRE_MIN_ORDER ||
        norm_order > GAUSS_LEGENDRE_MAX_ORDER) {
        return RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
    }

    shape_values = malloc((size_t)node_count * sizeof(shape_values[0]));
    if (shape_values == NULL) {
        return RELATIVE_ERROR_NORM_ALLOCATION_FAILED;
    }

    for (int ix = 0; ix < nx_cells; ++ix) {
        const double cell_xmin = (double)ix * dx;
        const double cell_xmax = cell_xmin + dx;

        for (int iy = 0; iy < ny_cells; ++iy) {
            const double cell_ymin = (double)iy * dy;
            const double cell_ymax = cell_ymin + dy;

            for (int iz = 0; iz < nz_cells; ++iz) {
                GaussPoint3D gauss_points[GAUSS_LEGENDRE_MAX_CUBE_POINTS];
                const double cell_zmin = (double)iz * dz;
                const double cell_zmax = cell_zmin + dz;
                int gauss_count = 0;

                if (gauss_legendre_cube(norm_order,
                                        cell_xmin, cell_xmax,
                                        cell_ymin, cell_ymax,
                                        cell_zmin, cell_zmax,
                                        gauss_points,
                                        GAUSS_LEGENDRE_MAX_CUBE_POINTS,
                                        &gauss_count) != GAUSS_OK) {
                    free(shape_values);
                    return RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
                }

                for (int gp = 0; gp < gauss_count; ++gp) {
                    const double x = gauss_points[gp].x;
                    const double y = gauss_points[gp].y;
                    const double z = gauss_points[gp].z;
                    const double weight = gauss_points[gp].weight;
                    const double v_exact =
                        analytical_potential_cube(x, y, z, L, V0,
                                                  ANALYTICAL_TERMS);
                    int value_count = 0;
                    double v_efg = 0.0;
                    const int mls_status =
                        mls_linear3d_shape_functions(nodes,
                                                     node_count,
                                                     x, y, z,
                                                     shape_values,
                                                     node_count,
                                                     &value_count);

                    if (mls_status != MLS_OK) {
                        for (int r = 0; r < ERROR_REGION_COUNT; ++r) {
                            if (error_region_contains(r, x, y, z, L)) {
                                ++regions[r].points_or_cells_used;
                                ++regions[r].mls_failures;
                            }
                        }
                        free(shape_values);
                        return RELATIVE_ERROR_NORM_MLS_FAILED;
                    }

                    for (int k = 0; k < value_count; ++k) {
                        v_efg += shape_values[k].phi
                                 * solution[shape_values[k].node_index];
                    }

                    {
                        const double diff = v_efg - v_exact;
                        const double abs_error = fabs(diff);

                        for (int r = 0; r < ERROR_REGION_COUNT; ++r) {
                            if (error_region_contains(r, x, y, z, L)) {
                                ++regions[r].points_or_cells_used;
                                regions[r].numerator_integral +=
                                    weight * diff * diff;
                                regions[r].denominator_integral +=
                                    weight * v_exact * v_exact;
                                regions[r].abs_error_integral +=
                                    weight * abs_error;
                                regions[r].volume_integral += weight;
                                if (isnan(regions[r].max_abs_error_sampled) ||
                                    abs_error >
                                    regions[r].max_abs_error_sampled) {
                                    regions[r].max_abs_error_sampled =
                                        abs_error;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    free(shape_values);

    for (int r = 0; r < ERROR_REGION_COUNT; ++r) {
        if (regions[r].denominator_integral <= 0.0 ||
            regions[r].volume_integral <= 0.0) {
            return RELATIVE_ERROR_NORM_ZERO_DENOMINATOR;
        }
        regions[r].relative_error_integral =
            sqrt(regions[r].numerator_integral /
                 regions[r].denominator_integral);
        regions[r].mean_abs_error_integral =
            regions[r].abs_error_integral / regions[r].volume_integral;
    }

    return RELATIVE_ERROR_NORM_OK;
}

static int run_error_region_study_refine15(void)
{
    const double L = 10.0;
    const double V0 = 10.0;
    const int nx_cells = 15;
    const int ny_cells = 15;
    const int nz_cells = 15;
    const int assembly_order = 2;
    const int norm_order = 6;
    StudyRegularSolution solution;
    ErrorRegionAccumulator regions[ERROR_REGION_COUNT];
    FILE *file = NULL;
    clock_t t_start;
    double eval_time_s = NAN;
    int eval_status;
    int rc = 1;

    study_regular_solution_init(&solution);
    error_regions_init(regions);

    printf("\n=== Error-region integral study: refine15 ===\n");
    printf("node grid:       15x15x15\n");
    printf("cells:           15x15x15\n");
    printf("assembly order:  %d\n", assembly_order);
    printf("norm order:      %d\n", norm_order);
    printf("solver:          gmres\n\n");

    if (solve_refine15_for_integration_study(assembly_order, &solution) != 0) {
        goto done;
    }

    if (!solution.gmres_res.converged) {
        fprintf(stderr,
                "STOP: GMRES did not converge for error-region study\n");
        goto done;
    }

    t_start = clock();
    eval_status = evaluate_error_regions_refine15(
        solution.nodes, solution.node_count, solution.x_sol,
        L, V0,
        nx_cells, ny_cells, nz_cells,
        norm_order,
        regions);
    eval_time_s = elapsed_s(t_start);

    if (eval_status != RELATIVE_ERROR_NORM_OK) {
        fprintf(stderr,
                "STOP: regional Eq.16 evaluation failed (status=%d)\n",
                eval_status);
        goto done;
    }

    file = fopen(ERROR_REGION_STUDY_CSV_PATH, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open CSV: %s\n",
                ERROR_REGION_STUDY_CSV_PATH);
        goto done;
    }

    fprintf(file,
            "case_name,region_name,norm_order,points_or_cells_used,"
            "numerator_integral,denominator_integral,"
            "relative_error_integral,mean_abs_error_integral,"
            "max_abs_error_sampled,mls_failures,eval_time_s\n");

    printf("region,points,relative_error,numerator,mean_abs,max_abs,"
           "mls_failures,numerator_fraction_full\n");
    for (int r = 0; r < ERROR_REGION_COUNT; ++r) {
        const double numerator_fraction =
            (regions[0].numerator_integral > 0.0)
            ? regions[r].numerator_integral / regions[0].numerator_integral
            : NAN;

        fprintf(file,
                "\"refine15\",\"%s\",%d,%d,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%d,%.17g\n",
                regions[r].name,
                norm_order,
                regions[r].points_or_cells_used,
                regions[r].numerator_integral,
                regions[r].denominator_integral,
                regions[r].relative_error_integral,
                regions[r].mean_abs_error_integral,
                regions[r].max_abs_error_sampled,
                regions[r].mls_failures,
                eval_time_s);

        printf("%s,%d,%.6e,%.6e,%.6e,%.6e,%d,%.3f\n",
               regions[r].name,
               regions[r].points_or_cells_used,
               regions[r].relative_error_integral,
               regions[r].numerator_integral,
               regions[r].mean_abs_error_integral,
               regions[r].max_abs_error_sampled,
               regions[r].mls_failures,
               numerator_fraction);
    }

    {
        int dominant = 2;
        for (int r = 3; r < ERROR_REGION_COUNT; ++r) {
            if (regions[r].numerator_integral >
                regions[dominant].numerator_integral) {
                dominant = r;
            }
        }
        printf("\nDominant non-global region by numerator: %s "
               "(%.2f%% of full-domain numerator)\n",
               regions[dominant].name,
               100.0 * regions[dominant].numerator_integral
               / regions[0].numerator_integral);
        printf("Upper split: edge bands %.2f%%, upper face interior %.2f%% "
               "of full-domain numerator\n",
               100.0 * regions[6].numerator_integral
               / regions[0].numerator_integral,
               100.0 * regions[7].numerator_integral
               / regions[0].numerator_integral);
    }

    if (fclose(file) != 0) {
        file = NULL;
        fprintf(stderr, "Could not close CSV: %s\n",
                ERROR_REGION_STUDY_CSV_PATH);
        goto done;
    }
    file = NULL;

    for (int r = 0; r < ERROR_REGION_COUNT; ++r) {
        if (regions[r].mls_failures != 0) {
            fprintf(stderr,
                    "STOP: region %s has MLS failures = %d\n",
                    regions[r].name,
                    regions[r].mls_failures);
            goto done;
        }
    }

    printf("\nCSV written: %s\n\n", ERROR_REGION_STUDY_CSV_PATH);
    rc = 0;

done:
    if (file != NULL) {
        fclose(file);
    }
    study_regular_solution_destroy(&solution);
    return rc;
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
                    CubeSparseSolver solver,
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

    Node3D         *nodes           = NULL;
    DirichletPoint *dp              = NULL;
    double         *K_dense         = NULL;
    double         *G_dense         = NULL;
    double         *F               = NULL;
    double         *q_vec           = NULL;
    double         *b_aug           = NULL;
    double         *x_sol           = NULL;
    double         *x_gmres_cmp     = NULL;
    double         *jacobi_precond  = NULL;
    MlsShapeValue  *shape_buf       = NULL;

    SparseCOO K_coo     = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCOO A_aug_coo = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCSR A_aug_csr = {NULL, NULL, NULL, 0, 0, 0};
    GmresResult       gmres_res  = {0, 0, 0.0, 0.0};
    LapackDenseSolveResult lapack_res = {0, -1, NAN, NAN, NAN};
    MlsConnectivityStats diag    = {0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    CubeSparseReport report;
    PlaneExportMetrics plane_metrics = {0, 0, NAN, NAN, NAN};
    RelativeErrorDomainResult domain_error = {NAN, 0.0, 0.0, 0, 0};
    CubeAnalyticalContext exact_ctx = {L, V0, ANALYTICAL_TERMS};

    double dense_rel_diff = -1.0;  /* < 0 means comparison was not run */
    int rc = 1;

    const int expected_nodes       = cube_regular_node_count(nx, ny, nz);
    const int expected_constraints = cube_dirichlet_point_count(nx, ny, nz);

    report_init(&report, label, solver);

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

    /* ------------------------------------------------------------ solver */

    x_sol = calloc((size_t)total_size, sizeof(double));
    if (!x_sol) { fprintf(stderr, "[%s] x_sol alloc failed\n", label); goto done; }

    jacobi_precond = malloc((size_t)total_size * sizeof(double));
    if (!jacobi_precond) { fprintf(stderr, "[%s] jacobi alloc failed\n", label); goto done; }
    build_jacobi_precond(&A_aug_csr, jacobi_precond);

    if (solver == CUBE_SOLVER_GMRES) {
        t_start = clock();

        /*
         * GMRES is used because the augmented system with Lagrange multipliers
         * is symmetric but indefinite (due to the zero block), so standard
         * conjugate gradient methods cannot be used.
         * Left Jacobi preconditioner applied: scales each row by 1/A_ii.
         */
        if (gmres_solve(&A_aug_csr, b_aug, x_sol,
                        gmres_tol, gmres_max_iter, gmres_restart,
                        jacobi_precond, &gmres_res) != GMRES_OK) {
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
    } else {
        int lapack_status;

        report.residual_initial = vector_norm2(b_aug, total_size);

        lapack_status = solve_augmented_lapack_dense(&A_aug_csr, b_aug,
                                                     x_sol, &lapack_res);
        report.lapack_available = lapack_res.available;
        report.lapack_info = lapack_res.lapack_info;
        report.lapack_conversion_time_s = lapack_res.conversion_time_s;
        report.lapack_solve_time_s = lapack_res.solve_time_s;
        report.lapack_residual_final = lapack_res.residual_final;
        report.residual_final = lapack_res.residual_final;
        report.solve_time_s = lapack_res.conversion_time_s
                              + lapack_res.solve_time_s;

        if (lapack_status != LAPACK_DENSE_SOLVER_OK) {
            fprintf(stderr,
                    "[%s] STOP: lapack-dense solver failed "
                    "(status=%d info=%d available=%d)\n",
                    label, lapack_status, lapack_res.lapack_info,
                    lapack_res.available);
            emit_required_report(&report);
            goto done;
        }

        x_gmres_cmp = calloc((size_t)total_size, sizeof(double));
        if (!x_gmres_cmp) {
            fprintf(stderr, "[%s] x_gmres_cmp alloc failed\n", label);
            goto done;
        }

        t_start = clock();
        if (gmres_solve(&A_aug_csr, b_aug, x_gmres_cmp,
                        gmres_tol, gmres_max_iter, gmres_restart,
                        jacobi_precond, &gmres_res) != GMRES_OK) {
            fprintf(stderr, "[%s] gmres_solve comparison returned error\n",
                    label);
            goto done;
        }
        t_sol = elapsed_s(t_start);
        report.gmres_converged = gmres_res.converged;
        report.gmres_iterations = gmres_res.iterations;
        report.gmres_lapack_rel_diff =
            relative_solution_difference(x_gmres_cmp, x_sol, total_size);

        if (!gmres_res.converged) {
            fprintf(stderr,
                    "[%s] WARNING: GMRES comparison did not converge "
                    "(iters=%d residual=%.6e)\n",
                    label, gmres_res.iterations, gmres_res.residual_final);
        }
        if (report.gmres_lapack_rel_diff > 1.0e-6) {
            fprintf(stderr,
                    "[%s] WARNING: GMRES vs LAPACK relative difference is "
                    "%.6e; possible convergence criterion or solver bug.\n",
                    label, report.gmres_lapack_rel_diff);
        }
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
                             shape_buf, L, V0, &plane_metrics) != 0) {
            goto done;
        }
        report.plane_relative_error = plane_metrics.relative_error;
        report.plane_max_abs_error = plane_metrics.max_abs_error;
        report.plane_mean_abs_error = plane_metrics.mean_abs_error;
    }

    if (relative_error_norm_domain_integral(nodes, node_count, x_sol,
                                            0.0, L, 0.0, L, 0.0, L,
                                            nx_cells, ny_cells, nz_cells,
                                            cube_exact_potential_callback,
                                            &exact_ctx,
                                            &domain_error)
        != RELATIVE_ERROR_NORM_OK) {
        report.mls_failures += domain_error.mls_failures;
        fprintf(stderr,
                "[%s] STOP: Eq.16 domain error integration failed "
                "(mls_failures=%d)\n",
                label, domain_error.mls_failures);
        emit_required_report(&report);
        goto done;
    }
    report.relative_error_global = domain_error.error;

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
        report.relative_error_discrete_global = rel_global;
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
        printf("  solver:                 %s\n",     solver_name(solver));
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

        if (solver == CUBE_SOLVER_GMRES) {
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
            printf("  solution time:          %.3f s\n", report.solve_time_s);
        } else {
            printf("--- LAPACK dense ---\n");
            printf("  available:              %s\n",
                   lapack_res.available ? "YES" : "NO");
            printf("  dgesv info:             %d\n", lapack_res.lapack_info);
            printf("  conversion time:        %.3f s\n",
                   lapack_res.conversion_time_s);
            printf("  factor/solve time:      %.3f s\n",
                   lapack_res.solve_time_s);
            printf("  residual final:         %.6e\n",
                   lapack_res.residual_final);
            printf("  total solution time:    %.3f s\n", report.solve_time_s);
            printf("\n");

            printf("--- GMRES comparison ---\n");
            printf("  restart:                %d\n",    gmres_restart);
            printf("  max iter:               %d\n",    gmres_max_iter);
            printf("  tolerance:              %.1e\n",  gmres_tol);
            printf("  iterations:             %d\n",    gmres_res.iterations);
            printf("  residual init:          %.6e\n",  gmres_res.residual_init);
            printf("  residual final:         %.6e\n",  gmres_res.residual_final);
            printf("  converged:              %s\n",
                   gmres_res.converged ? "YES" : "NO");
            printf("  GMRES vs LAPACK diff:   %.6e\n",
                   report.gmres_lapack_rel_diff);
            printf("  GMRES comparison time:  %.3f s\n", t_sol);
        }
        if (dense_rel_diff >= 0.0) {
            printf("  selected vs dense diff: %.3e\n", dense_rel_diff);
        }
        printf("\n");

        printf("--- Errors  (sample grid %dx%dx%d,"
               " valid=%d interior=%d) ---\n",
               sample_n, sample_n, sample_n, valid_all, valid_int);
        /*
         * The boundary-compatible Dirichlet policy grounds the upper edges and
         * applies V0 only on the open top face. The global max error is still
         * useful for regression, but the interior relative error is the main
         * quality metric for the smooth part of the field.
         */
        printf("  Eq.16 quadrature pts:   %d\n",
               domain_error.quadrature_points);
        printf("  Eq.16 numerator int:    %.6e\n",
               domain_error.numerator_integral);
        printf("  Eq.16 denominator int:  %.6e\n",
               domain_error.denominator_integral);
        printf("  max abs error:          %.6e\n", max_abs_err);
        printf("  rel error Eq.16 domain: %.6e\n", domain_error.error);
        printf("  rel error discrete global:   %.6e\n", rel_global);
        printf("  rel error discrete interior: %.6e\n", rel_intrnl);
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
    free(x_gmres_cmp);
    free(jacobi_precond);
    free(shape_buf);

    sparse_coo_destroy(&K_coo);
    sparse_coo_destroy(&A_aug_coo);
    sparse_csr_destroy(&A_aug_csr);

    return rc;
}

/*
 * Non-uniform EFG pipeline using cube_generate_article_cloud.
 *
 * Identical to run_case except for node generation: uses a coarse base_n^3
 * regular grid plus fine interior-in-xy slices near z = L (upper zone).
 */
static int run_case_nonuniform(const char *label,
                               double L, double V0,
                               int base_n, int top_n,
                               int n_extra_slices, double z_frac,
                               int nx_cells, int ny_cells, int nz_cells,
                               double gmres_tol, int gmres_restart,
                               int gmres_max_iter,
                               int sample_n,
                               CubeSparseSolver solver,
                               const PlaneExportConfig *plane_export)
{
    int node_count      = 0;
    int point_count     = 0;
    int total_size      = 0;
    int mls_eval_fail   = 0;
    int coo_count_aug   = 0;
    double t_asm = 0.0, t_sol = 0.0;
    clock_t t_start;

    Node3D         *nodes          = NULL;
    DirichletPoint *dp             = NULL;
    double         *K_dense        = NULL;
    double         *G_dense        = NULL;
    double         *F              = NULL;
    double         *q_vec          = NULL;
    double         *b_aug          = NULL;
    double         *x_sol          = NULL;
    double         *x_gmres_cmp    = NULL;
    double         *jacobi_precond = NULL;
    MlsShapeValue  *shape_buf      = NULL;

    SparseCOO K_coo     = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCOO A_aug_coo = {NULL, NULL, NULL, 0, 0, 0, 0};
    SparseCSR A_aug_csr = {NULL, NULL, NULL, 0, 0, 0};
    GmresResult       gmres_res  = {0, 0, 0.0, 0.0};
    LapackDenseSolveResult lapack_res = {0, -1, NAN, NAN, NAN};
    MlsConnectivityStats diag    = {0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    CubeSparseReport report;
    PlaneExportMetrics plane_metrics = {0, 0, NAN, NAN, NAN};
    RelativeErrorDomainResult domain_error = {NAN, 0.0, 0.0, 0, 0};
    CubeAnalyticalContext exact_ctx = {L, V0, ANALYTICAL_TERMS};

    int rc = 1;

    const int max_nodes = cube_article_cloud_max_nodes(base_n, top_n,
                                                       n_extra_slices);
    if (max_nodes < 0) {
        fprintf(stderr, "[%s] invalid cloud parameters\n", label);
        return 1;
    }
    /* Conservative upper bound for constraints: at most all nodes on surface */
    const int max_constraints = max_nodes;

    report_init(&report, label, solver);

    /* ---------------------------------------------------------------- alloc */

    nodes     = malloc((size_t)max_nodes * sizeof(*nodes));
    dp        = malloc((size_t)max_constraints * sizeof(*dp));
    K_dense   = calloc((size_t)max_nodes * (size_t)max_nodes, sizeof(double));
    G_dense   = calloc((size_t)max_constraints * (size_t)max_nodes,
                       sizeof(double));
    F         = calloc((size_t)max_nodes, sizeof(double));
    q_vec     = calloc((size_t)max_constraints, sizeof(double));
    shape_buf = malloc((size_t)max_nodes * sizeof(*shape_buf));

    if (!nodes || !dp || !K_dense || !G_dense || !F || !q_vec || !shape_buf) {
        fprintf(stderr, "[%s] allocation failed\n", label);
        goto nu_done;
    }

    /* ------------------------------------------------------- cloud & support */

    if (cube_generate_article_cloud(L, base_n, top_n, n_extra_slices, z_frac,
                                    nodes, max_nodes,
                                    &node_count) != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[%s] non-uniform cloud generation failed\n", label);
        goto nu_done;
    }

    if (support_assign_article_default(nodes, node_count) != SUPPORT_OK) {
        fprintf(stderr, "[%s] support assignment failed\n", label);
        goto nu_done;
    }

    /* ---------------------------------------------------- Dirichlet points */

    if (cube_generate_dirichlet_points_from_nodes(L, V0,
                                                  nodes, node_count,
                                                  dp, max_constraints,
                                                  &point_count)
        != CUBE_PROBLEM_OK) {
        fprintf(stderr, "[%s] Dirichlet generation failed\n", label);
        goto nu_done;
    }

    total_size = node_count + point_count;
    report.nodes          = node_count;
    report.constraints    = point_count;
    report.augmented_size = total_size;

    /* -------------------------------------------------- MLS diagnostic */

    mls_connectivity_stats_uniform_grid(nodes, node_count,
                                        0.0, L, nx_cells,
                                        0.0, L, ny_cells,
                                        0.0, L, nz_cells,
                                        &diag);
    report.support_min       = diag.min_nodes;
    report.support_mean      = diag.mean_nodes;
    report.support_max       = diag.max_nodes;
    report.support_lt_4      = diag.n_invalid;
    report.support_lt_8      = diag.n_invalid + diag.n_suspect;
    report.mls_failures      = diag.n_moment_fail;
    report.max_cond_estimate = diag.max_cond;

    if (report.support_lt_4 > 0) {
        fprintf(stderr, "[%s] STOP: support_lt_4 = %d (> 0)\n",
                label, report.support_lt_4);
        emit_required_report(&report);
        goto nu_done;
    }

    if (report.mls_failures > 0) {
        fprintf(stderr, "[%s] STOP: mls_failures = %d (> 0)\n",
                label, report.mls_failures);
        emit_required_report(&report);
        goto nu_done;
    }

    /* ----------------------------------------- dense K and G assembly */

    t_start = clock();

    if (global_stiffness_assemble_dense(nodes, node_count,
                                        0.0, L, 0.0, L, 0.0, L,
                                        nx_cells, ny_cells, nz_cells,
                                        K_dense) != GLOBAL_STIFFNESS_OK) {
        fprintf(stderr, "[%s] K assembly failed\n", label);
        goto nu_done;
    }

    if (dirichlet_assemble_constraints_dense(nodes, node_count,
                                             dp, point_count,
                                             G_dense, q_vec) != DIRICHLET_OK) {
        fprintf(stderr, "[%s] G assembly failed\n", label);
        goto nu_done;
    }

    /* --------------------------------- K_dense → K_coo (exact zeros dropped) */

    {
        int k_cap = node_count * 100 + 16;
        if (sparse_coo_create(&K_coo, node_count, node_count, k_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] K_coo alloc failed\n", label);
            goto nu_done;
        }
        for (int r = 0; r < node_count; ++r) {
            for (int c = 0; c < node_count; ++c) {
                double v = K_dense[r * node_count + c];
                if (v != 0.0) {
                    if (sparse_coo_add(&K_coo, r, c, v) != SPARSE_OK) {
                        fprintf(stderr, "[%s] K_coo append failed\n", label);
                        goto nu_done;
                    }
                }
            }
        }
    }

    /* -------------------- augmented COO: K block + G / G^T (zeros dropped) */

    {
        int aug_cap = K_coo.count + point_count * 60 + 16;
        if (sparse_coo_create(&A_aug_coo, total_size, total_size, aug_cap)
            != SPARSE_OK) {
            fprintf(stderr, "[%s] A_aug_coo alloc failed\n", label);
            goto nu_done;
        }

        for (int e = 0; e < K_coo.count; ++e) {
            if (sparse_coo_add(&A_aug_coo,
                               K_coo.row[e], K_coo.col[e], K_coo.val[e])
                != SPARSE_OK) {
                fprintf(stderr, "[%s] A_aug K-block append failed\n", label);
                goto nu_done;
            }
        }

        for (int ci = 0; ci < point_count; ++ci) {
            for (int ni = 0; ni < node_count; ++ni) {
                double val = G_dense[ci * node_count + ni];
                if (val != 0.0) {
                    if (sparse_coo_add(&A_aug_coo, ni, node_count + ci, val)
                        != SPARSE_OK) {
                        fprintf(stderr,
                                "[%s] A_aug G^T append failed\n", label);
                        goto nu_done;
                    }
                    if (sparse_coo_add(&A_aug_coo, node_count + ci, ni, val)
                        != SPARSE_OK) {
                        fprintf(stderr, "[%s] A_aug G append failed\n", label);
                        goto nu_done;
                    }
                }
            }
        }
    }

    /* b_aug: F (zeros) followed by q */
    b_aug = calloc((size_t)total_size, sizeof(double));
    if (!b_aug) {
        fprintf(stderr, "[%s] b_aug alloc failed\n", label);
        goto nu_done;
    }
    for (int ci = 0; ci < point_count; ++ci) {
        b_aug[node_count + ci] = q_vec[ci];
    }

    coo_count_aug = A_aug_coo.count;
    if (sparse_coo_to_csr(&A_aug_coo, &A_aug_csr) != SPARSE_OK) {
        fprintf(stderr, "[%s] A_aug CSR conversion failed\n", label);
        goto nu_done;
    }

    t_asm = elapsed_s(t_start);
    report.K_nnz           = K_coo.count;
    report.A_aug_nnz       = A_aug_csr.nnz;
    report.assembly_time_s = t_asm;

    /* ------------------------------------------------------------ solver */

    x_sol = calloc((size_t)total_size, sizeof(double));
    if (!x_sol) {
        fprintf(stderr, "[%s] x_sol alloc failed\n", label);
        goto nu_done;
    }

    jacobi_precond = malloc((size_t)total_size * sizeof(double));
    if (!jacobi_precond) {
        fprintf(stderr, "[%s] jacobi alloc failed\n", label);
        goto nu_done;
    }
    build_jacobi_precond(&A_aug_csr, jacobi_precond);

    if (solver == CUBE_SOLVER_GMRES) {
        t_start = clock();

        if (gmres_solve(&A_aug_csr, b_aug, x_sol,
                        gmres_tol, gmres_max_iter, gmres_restart,
                        jacobi_precond, &gmres_res) != GMRES_OK) {
            fprintf(stderr, "[%s] gmres_solve returned error\n", label);
            goto nu_done;
        }

        t_sol = elapsed_s(t_start);
        report.gmres_converged   = gmres_res.converged;
        report.gmres_iterations  = gmres_res.iterations;
        report.residual_initial  = gmres_res.residual_init;
        report.residual_final    = gmres_res.residual_final;
        report.solve_time_s      = t_sol;

        if (!gmres_res.converged) {
            fprintf(stderr, "[%s] STOP: GMRES did not converge\n", label);
            emit_required_report(&report);
            goto nu_done;
        }
    } else {
        int lapack_status;

        report.residual_initial = vector_norm2(b_aug, total_size);

        lapack_status = solve_augmented_lapack_dense(&A_aug_csr, b_aug,
                                                     x_sol, &lapack_res);
        report.lapack_available = lapack_res.available;
        report.lapack_info = lapack_res.lapack_info;
        report.lapack_conversion_time_s = lapack_res.conversion_time_s;
        report.lapack_solve_time_s = lapack_res.solve_time_s;
        report.lapack_residual_final = lapack_res.residual_final;
        report.residual_final = lapack_res.residual_final;
        report.solve_time_s = lapack_res.conversion_time_s
                              + lapack_res.solve_time_s;

        if (lapack_status != LAPACK_DENSE_SOLVER_OK) {
            fprintf(stderr,
                    "[%s] STOP: lapack-dense solver failed "
                    "(status=%d info=%d available=%d)\n",
                    label, lapack_status, lapack_res.lapack_info,
                    lapack_res.available);
            emit_required_report(&report);
            goto nu_done;
        }

        x_gmres_cmp = calloc((size_t)total_size, sizeof(double));
        if (!x_gmres_cmp) {
            fprintf(stderr, "[%s] x_gmres_cmp alloc failed\n", label);
            goto nu_done;
        }

        t_start = clock();
        if (gmres_solve(&A_aug_csr, b_aug, x_gmres_cmp,
                        gmres_tol, gmres_max_iter, gmres_restart,
                        jacobi_precond, &gmres_res) != GMRES_OK) {
            fprintf(stderr, "[%s] gmres_solve comparison returned error\n",
                    label);
            goto nu_done;
        }
        t_sol = elapsed_s(t_start);
        report.gmres_converged = gmres_res.converged;
        report.gmres_iterations = gmres_res.iterations;
        report.gmres_lapack_rel_diff =
            relative_solution_difference(x_gmres_cmp, x_sol, total_size);

        if (!gmres_res.converged) {
            fprintf(stderr,
                    "[%s] WARNING: GMRES comparison did not converge "
                    "(iters=%d residual=%.6e)\n",
                    label, gmres_res.iterations, gmres_res.residual_final);
        }
        if (report.gmres_lapack_rel_diff > 1.0e-6) {
            fprintf(stderr,
                    "[%s] WARNING: GMRES vs LAPACK relative difference is "
                    "%.6e; possible convergence criterion or solver bug.\n",
                    label, report.gmres_lapack_rel_diff);
        }
    }

    if (plane_export != NULL && plane_export->enabled) {
        if (export_plane_csv(plane_export, nodes, node_count, x_sol,
                             shape_buf, L, V0, &plane_metrics) != 0) {
            goto nu_done;
        }
        report.plane_relative_error = plane_metrics.relative_error;
        report.plane_max_abs_error = plane_metrics.max_abs_error;
        report.plane_mean_abs_error = plane_metrics.mean_abs_error;
    }

    if (relative_error_norm_domain_integral(nodes, node_count, x_sol,
                                            0.0, L, 0.0, L, 0.0, L,
                                            nx_cells, ny_cells, nz_cells,
                                            cube_exact_potential_callback,
                                            &exact_ctx,
                                            &domain_error)
        != RELATIVE_ERROR_NORM_OK) {
        report.mls_failures += domain_error.mls_failures;
        fprintf(stderr,
                "[%s] STOP: Eq.16 domain error integration failed "
                "(mls_failures=%d)\n",
                label, domain_error.mls_failures);
        emit_required_report(&report);
        goto nu_done;
    }
    report.relative_error_global = domain_error.error;

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
                                         sx, sy, sz, L, V0, ANALYTICAL_TERMS);
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
        report.max_abs_error            = max_abs_err;
        report.relative_error_discrete_global = rel_global;
        report.relative_error_interior  = rel_intrnl;
        report.mls_failures            += mls_eval_fail;

        printf("\n");
        printf("=== %s ===\n", label);
        printf("\n");

        printf("--- Problem ---\n");
        printf("  nodes:                  %d\n", node_count);
        printf("  constraints:            %d\n", point_count);
        printf("  augmented size:         %d\n", total_size);
        printf("  solver:                 %s\n", solver_name(solver));
        printf("  cloud: base=%dx%dx%d top=%dx%d slices=%d z_frac=%.2f\n",
               base_n, base_n, base_n, top_n, top_n, n_extra_slices, z_frac);
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

        if (solver == CUBE_SOLVER_GMRES) {
            printf("--- GMRES ---\n");
            printf("  restart:                %d\n",   gmres_restart);
            printf("  max iter:               %d\n",   gmres_max_iter);
            printf("  tolerance:              %.1e\n", gmres_tol);
            printf("  iterations:             %d\n",   gmres_res.iterations);
            printf("  residual init:          %.6e\n", gmres_res.residual_init);
            printf("  residual final:         %.6e\n", gmres_res.residual_final);
            printf("  rel residual:           %.3e\n",
                   (gmres_res.residual_init > 0.0)
                   ? gmres_res.residual_final / gmres_res.residual_init : 0.0);
            printf("  converged:              %s\n",
                   gmres_res.converged ? "YES" : "NO");
            printf("  solution time:          %.3f s\n", report.solve_time_s);
        } else {
            printf("--- LAPACK dense ---\n");
            printf("  available:              %s\n",
                   lapack_res.available ? "YES" : "NO");
            printf("  dgesv info:             %d\n", lapack_res.lapack_info);
            printf("  conversion time:        %.3f s\n",
                   lapack_res.conversion_time_s);
            printf("  factor/solve time:      %.3f s\n",
                   lapack_res.solve_time_s);
            printf("  residual final:         %.6e\n",
                   lapack_res.residual_final);
            printf("  total solution time:    %.3f s\n", report.solve_time_s);
            printf("\n");

            printf("--- GMRES comparison ---\n");
            printf("  restart:                %d\n",   gmres_restart);
            printf("  max iter:               %d\n",   gmres_max_iter);
            printf("  tolerance:              %.1e\n", gmres_tol);
            printf("  iterations:             %d\n",   gmres_res.iterations);
            printf("  residual init:          %.6e\n", gmres_res.residual_init);
            printf("  residual final:         %.6e\n", gmres_res.residual_final);
            printf("  converged:              %s\n",
                   gmres_res.converged ? "YES" : "NO");
            printf("  GMRES vs LAPACK diff:   %.6e\n",
                   report.gmres_lapack_rel_diff);
            printf("  GMRES comparison time:  %.3f s\n", t_sol);
        }
        printf("\n");

        printf("--- Errors  (sample grid %dx%dx%d,"
               " valid=%d interior=%d) ---\n",
               sample_n, sample_n, sample_n, valid_all, valid_int);
        printf("  Eq.16 quadrature pts:   %d\n",
               domain_error.quadrature_points);
        printf("  Eq.16 numerator int:    %.6e\n",
               domain_error.numerator_integral);
        printf("  Eq.16 denominator int:  %.6e\n",
               domain_error.denominator_integral);
        printf("  max abs error:          %.6e\n", max_abs_err);
        printf("  rel error Eq.16 domain: %.6e\n", domain_error.error);
        printf("  rel error discrete global:   %.6e\n", rel_global);
        printf("  rel error discrete interior: %.6e\n", rel_intrnl);
        if (mls_eval_fail > 0) {
            printf("  MLS eval failures:      %d\n", mls_eval_fail);
        }
        printf("\n");

        if (emit_required_report(&report) != 0) {
            goto nu_done;
        }

        if (report.mls_failures > 0) {
            fprintf(stderr, "[%s] STOP: mls_failures = %d (> 0)\n",
                    label, report.mls_failures);
            goto nu_done;
        }
    }

    rc = 0;

nu_done:
    free(nodes);
    free(dp);
    free(K_dense);
    free(G_dense);
    free(F);
    free(q_vec);
    free(b_aug);
    free(x_sol);
    free(x_gmres_cmp);
    free(jacobi_precond);
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
 * The nonuniform case uses cube_generate_article_cloud (Fig. 2 inspired).
 */
int main(int argc, char **argv)
{
    CubeSparseSelection selection;
    CubeSparseSolver solver;
    CubeIntegrationStudy integration_study;
    int error_region_study;
    int status = 0;
    int parse_status = parse_args(argc, argv, &selection, &solver,
                                  &integration_study,
                                  &error_region_study);

    if (parse_status == 1) {
        return 0;
    }
    if (parse_status != 0) {
        return 1;
    }

    if (integration_study != CUBE_INTEGRATION_STUDY_NONE &&
        error_region_study) {
        fprintf(stderr,
                "Choose either --integration-study or --error-region-study, "
                "not both.\n");
        return 1;
    }

    if (integration_study != CUBE_INTEGRATION_STUDY_NONE) {
        if (selection != CUBE_SPARSE_REFINE15 &&
            selection != CUBE_SPARSE_ALL) {
            fprintf(stderr,
                    "--integration-study is defined only for --case refine15.\n");
            return 1;
        }
        if (solver != CUBE_SOLVER_GMRES) {
            fprintf(stderr,
                    "--integration-study uses the GMRES reference path; "
                    "ignoring --solver %s.\n",
                    solver_name(solver));
        }
        if (integration_study == CUBE_INTEGRATION_STUDY_NORM) {
            return run_integration_study_norm_refine15();
        }
        if (integration_study == CUBE_INTEGRATION_STUDY_ASSEMBLY) {
            return run_integration_study_assembly_refine15();
        }
    }

    if (error_region_study) {
        if (selection != CUBE_SPARSE_REFINE15 &&
            selection != CUBE_SPARSE_ALL) {
            fprintf(stderr,
                    "--error-region-study is defined only for --case refine15.\n");
            return 1;
        }
        if (solver != CUBE_SOLVER_GMRES) {
            fprintf(stderr,
                    "--error-region-study uses the GMRES reference path; "
                    "ignoring --solver %s.\n",
                    solver_name(solver));
        }
        return run_error_region_study_refine15();
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
            solver,
            /*dense_cmp=*/1,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_TARGET) {
        status |= run_case(
            "11x11x11 target (15x15x15 integration cells)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/11, 11, 11,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/200,
            /*max_iter=*/10000,
            /*sample_n=*/11,
            solver,
            /*dense_cmp=*/0,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_REFINE13) {
        status |= run_case(
            "13x13x13 refine13 (15x15x15 integration cells)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/13, 13, 13,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            solver,
            /*dense_cmp=*/0,
            /*plane_export=*/NULL);
    }

    if (selection == CUBE_SPARSE_ALL ||
        selection == CUBE_SPARSE_REFINE15) {
        status |= run_case(
            "15x15x15 refine15 (15x15x15 integration cells)",
            /*L=*/10.0, /*V0=*/10.0,
            /*nx,ny,nz=*/15, 15, 15,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            solver,
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
            solver,
            /*dense_cmp=*/0,
            &plane_export);
    }

    if (selection == CUBE_SPARSE_NONUNIFORM) {
        const PlaneExportConfig plane_export = {
            1,
            5.33,
            101,
            101,
            NONUNIFORM_CSV_PATH
        };
        status |= run_case_nonuniform(
            "nonuniform cloud (base=11 top=13 slices=4 z_frac=0.30)",
            /*L=*/10.0, /*V0=*/10.0,
            /*base_n=*/11, /*top_n=*/13,
            /*n_extra_slices=*/4, /*z_frac=*/0.30,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            solver,
            &plane_export);
    }

    if (selection == CUBE_SPARSE_NONUNIFORM_REFINE) {
        const PlaneExportConfig plane_export = {
            1,
            5.33,
            101,
            101,
            NONUNIFORM_REFINE_CSV_PATH
        };
        status |= run_case_nonuniform(
            "nonuniform refine cloud (base=13 top=15 slices=4 z_frac=0.30)",
            /*L=*/10.0, /*V0=*/10.0,
            /*base_n=*/13, /*top_n=*/15,
            /*n_extra_slices=*/4, /*z_frac=*/0.30,
            /*cells=*/15, 15, 15,
            /*gmres_tol=*/1e-9,
            /*restart=*/300,
            /*max_iter=*/20000,
            /*sample_n=*/11,
            solver,
            &plane_export);
    }

    return status;
}
