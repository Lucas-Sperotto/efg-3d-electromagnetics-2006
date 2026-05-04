#include "analytical.h"
#include "cube_problem.h"
#include "dense_solver.h"
#include "dirichlet.h"
#include "error_norm.h"
#include "global_stiffness.h"
#include "lagrange_system.h"
#include "mls.h"
#include "support.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define DEFAULT_L 10.0
#define DEFAULT_V0 10.0
#define DEFAULT_NX 3
#define DEFAULT_NY 3
#define DEFAULT_NZ 3
#define DEFAULT_NX_CELLS 2
#define DEFAULT_NY_CELLS 2
#define DEFAULT_NZ_CELLS 2
#define DEFAULT_SAMPLE_NX 5
#define DEFAULT_SAMPLE_NY 5
#define DEFAULT_SAMPLE_NZ 5
#define DEFAULT_OUTPUT_PATH "data/output/cube_smoke_solution.csv"
#define ANALYTICAL_TERMS 25

typedef struct CubeSmokeConfig {
    double L;
    double V0;
    int nx;
    int ny;
    int nz;
    int nx_cells;
    int ny_cells;
    int nz_cells;
    int sample_nx;
    int sample_ny;
    int sample_nz;
    const char *output_path;
} CubeSmokeConfig;

static int ensure_directory(const char *path)
{
    if (mkdir(path, 0777) == 0) {
        return 0;
    }

    if (errno == EEXIST) {
        return 0;
    }

    fprintf(stderr, "Could not create directory %s: %s\n", path, strerror(errno));
    return 1;
}

static int ensure_output_directory(void)
{
    if (ensure_directory("data") != 0) {
        return 1;
    }

    return ensure_directory("data/output");
}

static void print_usage(const char *program_name)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s\n", program_name);
    fprintf(stderr,
            "  %s nx ny nz nx_cells ny_cells nz_cells sample_nx sample_ny sample_nz output.csv\n",
            program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Constraints:\n");
    fprintf(stderr, "  nx, ny, nz >= 3\n");
    fprintf(stderr, "  nx_cells, ny_cells, nz_cells >= 1\n");
    fprintf(stderr, "  sample_nx, sample_ny, sample_nz >= 2\n");
    fprintf(stderr, "  output.csv must be a non-empty path\n");
}

static int parse_int_at_least(const char *text, int minimum, const char *name, int *value)
{
    char *end = NULL;
    long parsed = 0;

    errno = 0;
    parsed = strtol(text, &end, 10);

    if (end == text || *end != '\0' || errno != 0 || parsed < minimum || parsed > INT_MAX) {
        fprintf(stderr, "Invalid %s: %s\n", name, text);
        return 1;
    }

    *value = (int)parsed;
    return 0;
}

static int parse_config(int argc, char **argv, CubeSmokeConfig *config)
{
    *config = (CubeSmokeConfig){
        DEFAULT_L,
        DEFAULT_V0,
        DEFAULT_NX,
        DEFAULT_NY,
        DEFAULT_NZ,
        DEFAULT_NX_CELLS,
        DEFAULT_NY_CELLS,
        DEFAULT_NZ_CELLS,
        DEFAULT_SAMPLE_NX,
        DEFAULT_SAMPLE_NY,
        DEFAULT_SAMPLE_NZ,
        DEFAULT_OUTPUT_PATH
    };

    if (argc == 1) {
        return 0;
    }

    if (argc != 11) {
        print_usage(argv[0]);
        return 1;
    }

    if (parse_int_at_least(argv[1], 3, "nx", &config->nx) != 0 ||
        parse_int_at_least(argv[2], 3, "ny", &config->ny) != 0 ||
        parse_int_at_least(argv[3], 3, "nz", &config->nz) != 0 ||
        parse_int_at_least(argv[4], 1, "nx_cells", &config->nx_cells) != 0 ||
        parse_int_at_least(argv[5], 1, "ny_cells", &config->ny_cells) != 0 ||
        parse_int_at_least(argv[6], 1, "nz_cells", &config->nz_cells) != 0 ||
        parse_int_at_least(argv[7], 2, "sample_nx", &config->sample_nx) != 0 ||
        parse_int_at_least(argv[8], 2, "sample_ny", &config->sample_ny) != 0 ||
        parse_int_at_least(argv[9], 2, "sample_nz", &config->sample_nz) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    if (argv[10][0] == '\0') {
        fprintf(stderr, "Invalid output path: empty string\n");
        print_usage(argv[0]);
        return 1;
    }

    config->output_path = argv[10];
    return 0;
}

static void print_config(const CubeSmokeConfig *config)
{
    printf("Cube smoke EFG configuration:\n");
    printf("L: %.17g\n", config->L);
    printf("V0: %.17g\n", config->V0);
    printf("nx ny nz: %d %d %d\n", config->nx, config->ny, config->nz);
    printf("nx_cells ny_cells nz_cells: %d %d %d\n",
           config->nx_cells,
           config->ny_cells,
           config->nz_cells);
    printf("sample_nx sample_ny sample_nz: %d %d %d\n",
           config->sample_nx,
           config->sample_ny,
           config->sample_nz);
    printf("output path: %s\n", config->output_path);
}

static int checked_product3(int a, int b, int c, int *result)
{
    const long long product = (long long)a * (long long)b * (long long)c;

    if (product <= 0 || product > (long long)INT_MAX) {
        return 1;
    }

    *result = (int)product;
    return 0;
}

static int check_status(const char *step, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s failed: expected code %d, got %d\n", step, expected, actual);
        return 1;
    }

    return 0;
}

static int nearest_node_index(const Node3D *nodes,
                              int node_count,
                              double x,
                              double y,
                              double z)
{
    int nearest_index = 0;
    double nearest_distance_squared = INFINITY;

    for (int i = 0; i < node_count; ++i) {
        const double dx = x - nodes[i].x;
        const double dy = y - nodes[i].y;
        const double dz = z - nodes[i].z;
        const double distance_squared = (dx * dx) + (dy * dy) + (dz * dz);

        if (distance_squared < nearest_distance_squared) {
            nearest_distance_squared = distance_squared;
            nearest_index = i;
        }
    }

    return nearest_index;
}

static int evaluate_mls_potential(const Node3D *nodes,
                                  int node_count,
                                  const double *solution,
                                  MlsShapeValue *shape_values,
                                  double x,
                                  double y,
                                  double z,
                                  double *potential)
{
    int value_count = 0;
    double reconstructed = 0.0;
    const int status = mls_linear3d_shape_functions(nodes,
                                                    node_count,
                                                    x,
                                                    y,
                                                    z,
                                                    shape_values,
                                                    node_count,
                                                    &value_count);

    if (status != MLS_OK) {
        fprintf(stderr,
                "MLS evaluation failed at (%.17g, %.17g, %.17g): code %d\n",
                x,
                y,
                z,
                status);
        return 1;
    }

    for (int i = 0; i < value_count; ++i) {
        reconstructed += shape_values[i].phi * solution[shape_values[i].node_index];
    }

    *potential = reconstructed;
    return 0;
}

/*
 * Didactic smoke reproduction of the electrostatic cube pipeline.
 *
 * The default run preserves the original tiny 3x3x3 node cloud and dense
 * direct solver. Optional command-line parameters let the same smoke pipeline
 * be inspected at other small sizes without recompilation. This is not a
 * quantitative reproduction of the article's final figures and it does not
 * implement GMRES.
 */
int main(int argc, char **argv)
{
    CubeSmokeConfig config;
    Node3D *nodes = NULL;
    DirichletPoint *points = NULL;
    MlsShapeValue *shape_values = NULL;
    double *K = NULL;
    double *G = NULL;
    double *q = NULL;
    double *F = NULL;
    double *A_aug = NULL;
    double *b_aug = NULL;
    double *solution = NULL;
    double *numerical_samples = NULL;
    double *exact_samples = NULL;
    double max_abs_error = 0.0;
    double relative_error = 0.0;
    int expected_node_count = 0;
    int node_count = 0;
    int expected_constraint_count = 0;
    int constraint_count = 0;
    int total_size = 0;
    int sample_count = 0;
    int sample_index = 0;
    int exit_code = 1;
    FILE *file = NULL;

    if (parse_config(argc, argv, &config) != 0) {
        return 1;
    }

    print_config(&config);

    expected_node_count = cube_regular_node_count(config.nx, config.ny, config.nz);
    expected_constraint_count = cube_dirichlet_point_count(config.nx, config.ny, config.nz);

    if (expected_node_count <= 0 || expected_constraint_count <= 0 ||
        checked_product3(config.sample_nx, config.sample_ny, config.sample_nz, &sample_count) != 0) {
        fprintf(stderr, "Invalid derived problem size.\n");
        return 1;
    }

    if (expected_node_count > INT_MAX - expected_constraint_count) {
        fprintf(stderr, "Augmented system size is too large.\n");
        return 1;
    }

    total_size = expected_node_count + expected_constraint_count;

    nodes = malloc((size_t)expected_node_count * sizeof(nodes[0]));
    points = malloc((size_t)expected_constraint_count * sizeof(points[0]));
    shape_values = malloc((size_t)expected_node_count * sizeof(shape_values[0]));
    K = malloc((size_t)expected_node_count * (size_t)expected_node_count * sizeof(K[0]));
    G = malloc((size_t)expected_constraint_count * (size_t)expected_node_count * sizeof(G[0]));
    q = malloc((size_t)expected_constraint_count * sizeof(q[0]));
    F = malloc((size_t)expected_node_count * sizeof(F[0]));
    A_aug = malloc((size_t)total_size * (size_t)total_size * sizeof(A_aug[0]));
    b_aug = malloc((size_t)total_size * sizeof(b_aug[0]));
    solution = malloc((size_t)total_size * sizeof(solution[0]));
    numerical_samples = malloc((size_t)sample_count * sizeof(numerical_samples[0]));
    exact_samples = malloc((size_t)sample_count * sizeof(exact_samples[0]));

    if (nodes == NULL || points == NULL || shape_values == NULL || K == NULL ||
        G == NULL || q == NULL || F == NULL || A_aug == NULL || b_aug == NULL ||
        solution == NULL || numerical_samples == NULL || exact_samples == NULL) {
        fprintf(stderr, "Could not allocate dense smoke problem arrays.\n");
        goto cleanup;
    }

    if (ensure_output_directory() != 0) {
        goto cleanup;
    }

    if (check_status("cube_generate_regular_nodes",
                     cube_generate_regular_nodes(config.L,
                                                 config.nx,
                                                 config.ny,
                                                 config.nz,
                                                 nodes,
                                                 expected_node_count,
                                                 &node_count),
                     CUBE_PROBLEM_OK) != 0) {
        goto cleanup;
    }

    if (node_count != expected_node_count) {
        fprintf(stderr, "Unexpected node count: %d\n", node_count);
        goto cleanup;
    }

    if (check_status("support_assign_article_default",
                     support_assign_article_default(nodes, node_count),
                     SUPPORT_OK) != 0) {
        goto cleanup;
    }

    if (check_status("cube_generate_dirichlet_points",
                     cube_generate_dirichlet_points(config.L,
                                                    config.V0,
                                                    config.nx,
                                                    config.ny,
                                                    config.nz,
                                                    points,
                                                    expected_constraint_count,
                                                    &constraint_count),
                     CUBE_PROBLEM_OK) != 0) {
        goto cleanup;
    }

    if (constraint_count != expected_constraint_count) {
        fprintf(stderr, "Unexpected constraint count: %d\n", constraint_count);
        goto cleanup;
    }

    if (check_status("global_stiffness_assemble_dense",
                     global_stiffness_assemble_dense(nodes,
                                                     node_count,
                                                     0.0,
                                                     config.L,
                                                     0.0,
                                                     config.L,
                                                     0.0,
                                                     config.L,
                                                     config.nx_cells,
                                                     config.ny_cells,
                                                     config.nz_cells,
                                                     K),
                     GLOBAL_STIFFNESS_OK) != 0) {
        goto cleanup;
    }

    for (int i = 0; i < node_count; ++i) {
        F[i] = 0.0;
    }

    if (check_status("dirichlet_assemble_constraints_dense",
                     dirichlet_assemble_constraints_dense(nodes,
                                                          node_count,
                                                          points,
                                                          constraint_count,
                                                          G,
                                                          q),
                     DIRICHLET_OK) != 0) {
        goto cleanup;
    }

    if (check_status("lagrange_assemble_augmented_system_dense",
                     lagrange_assemble_augmented_system_dense(K,
                                                              F,
                                                              node_count,
                                                              G,
                                                              q,
                                                              constraint_count,
                                                              A_aug,
                                                              b_aug),
                     LAGRANGE_SYSTEM_OK) != 0) {
        goto cleanup;
    }

    if (check_status("dense_solve",
                     dense_solve(A_aug, b_aug, total_size, solution),
                     DENSE_SOLVER_OK) != 0) {
        goto cleanup;
    }

    file = fopen(config.output_path, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open output CSV: %s\n", config.output_path);
        goto cleanup;
    }

    fprintf(file, "x,y,z,V_eval,V_exact,abs_error,nearest_node_index,V_coeff_nearest\n");

    for (int ix = 0; ix < config.sample_nx; ++ix) {
        const double x = (config.L * (double)ix) / (double)(config.sample_nx - 1);

        for (int iy = 0; iy < config.sample_ny; ++iy) {
            const double y = (config.L * (double)iy) / (double)(config.sample_ny - 1);

            for (int iz = 0; iz < config.sample_nz; ++iz) {
                const double z = (config.L * (double)iz) / (double)(config.sample_nz - 1);
                const int nearest_index = nearest_node_index(nodes, node_count, x, y, z);
                double reconstructed = 0.0;
                double exact = 0.0;
                double abs_error = 0.0;

                if (evaluate_mls_potential(nodes,
                                           node_count,
                                           solution,
                                           shape_values,
                                           x,
                                           y,
                                           z,
                                           &reconstructed) != 0) {
                    goto cleanup;
                }

                exact = analytical_potential_cube(x,
                                                  y,
                                                  z,
                                                  config.L,
                                                  config.V0,
                                                  ANALYTICAL_TERMS);
                abs_error = fabs(reconstructed - exact);

                numerical_samples[sample_index] = reconstructed;
                exact_samples[sample_index] = exact;

                if (abs_error > max_abs_error) {
                    max_abs_error = abs_error;
                }

                fprintf(file,
                        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g\n",
                        x,
                        y,
                        z,
                        reconstructed,
                        exact,
                        abs_error,
                        nearest_index,
                        solution[nearest_index]);

                ++sample_index;
            }
        }
    }

    if (fclose(file) != 0) {
        file = NULL;
        fprintf(stderr, "Could not close output CSV: %s\n", config.output_path);
        goto cleanup;
    }
    file = NULL;

    if (check_status("relative_error_norm",
                     relative_error_norm(numerical_samples,
                                         exact_samples,
                                         sample_count,
                                         &relative_error),
                     RELATIVE_ERROR_NORM_OK) != 0) {
        goto cleanup;
    }

    printf("Cube smoke EFG run completed.\n");
    printf("nodes: %d\n", node_count);
    printf("constraints: %d\n", constraint_count);
    printf("augmented system size: %d\n", total_size);
    printf("sample points: %d\n", sample_count);
    printf("max abs error: %.17g\n", max_abs_error);
    printf("relative discrete error: %.17g\n", relative_error);
    printf("csv: %s\n", config.output_path);

    exit_code = 0;

cleanup:
    if (file != NULL) {
        fclose(file);
    }

    free(nodes);
    free(points);
    free(shape_values);
    free(K);
    free(G);
    free(q);
    free(F);
    free(A_aug);
    free(b_aug);
    free(solution);
    free(numerical_samples);
    free(exact_samples);

    return exit_code;
}
