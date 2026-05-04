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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define CUBE_L 10.0
#define CUBE_V0 10.0
#define NODE_NX 3
#define NODE_NY 3
#define NODE_NZ 3
#define NODE_COUNT 27
#define CONSTRAINT_COUNT 26
#define TOTAL_SIZE (NODE_COUNT + CONSTRAINT_COUNT)
#define NODE_MATRIX_SIZE (NODE_COUNT * NODE_COUNT)
#define CONSTRAINT_MATRIX_SIZE (CONSTRAINT_COUNT * NODE_COUNT)
#define AUG_MATRIX_SIZE (TOTAL_SIZE * TOTAL_SIZE)
#define CELL_COUNT 2
#define SAMPLE_NX 5
#define SAMPLE_NY 5
#define SAMPLE_NZ 5
#define SAMPLE_COUNT (SAMPLE_NX * SAMPLE_NY * SAMPLE_NZ)
#define ANALYTICAL_TERMS 25

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

static int check_status(const char *step, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s failed: expected code %d, got %d\n", step, expected, actual);
        return 1;
    }

    return 0;
}

static int nearest_node_index(const Node3D nodes[NODE_COUNT], double x, double y, double z)
{
    int nearest_index = 0;
    double nearest_distance_squared = INFINITY;

    for (int i = 0; i < NODE_COUNT; ++i) {
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

static int evaluate_mls_potential(const Node3D nodes[NODE_COUNT],
                                  const double solution[TOTAL_SIZE],
                                  double x,
                                  double y,
                                  double z,
                                  double *potential)
{
    MlsShapeValue values[NODE_COUNT];
    int value_count = 0;
    double reconstructed = 0.0;
    const int status = mls_linear3d_shape_functions(nodes,
                                                    NODE_COUNT,
                                                    x,
                                                    y,
                                                    z,
                                                    values,
                                                    NODE_COUNT,
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
        reconstructed += values[i].phi * solution[values[i].node_index];
    }

    *potential = reconstructed;
    return 0;
}

/*
 * Didactic smoke reproduction of the electrostatic cube pipeline.
 *
 * This app intentionally uses a tiny 3x3x3 node cloud and a dense direct
 * solver. It exports MLS-reconstructed values for later inspection, but it is
 * not a quantitative reproduction of the article's final figures and it does
 * not implement GMRES.
 */
int main(void)
{
    const char *output_path = "data/output/cube_smoke_solution.csv";
    Node3D nodes[NODE_COUNT];
    DirichletPoint points[CONSTRAINT_COUNT];
    double K[NODE_MATRIX_SIZE];
    double G[CONSTRAINT_MATRIX_SIZE];
    double q[CONSTRAINT_COUNT];
    double F[NODE_COUNT];
    double A_aug[AUG_MATRIX_SIZE];
    double b_aug[TOTAL_SIZE];
    double solution[TOTAL_SIZE];
    double numerical_samples[SAMPLE_COUNT];
    double exact_samples[SAMPLE_COUNT];
    double max_abs_error = 0.0;
    double relative_error = 0.0;
    int node_count = 0;
    int constraint_count = 0;
    int sample_index = 0;
    FILE *file = NULL;

    if (ensure_output_directory() != 0) {
        return 1;
    }

    if (check_status("cube_generate_regular_nodes",
                     cube_generate_regular_nodes(CUBE_L,
                                                 NODE_NX,
                                                 NODE_NY,
                                                 NODE_NZ,
                                                 nodes,
                                                 NODE_COUNT,
                                                 &node_count),
                     CUBE_PROBLEM_OK) != 0) {
        return 1;
    }

    if (node_count != NODE_COUNT) {
        fprintf(stderr, "Unexpected node count: %d\n", node_count);
        return 1;
    }

    if (check_status("support_assign_article_default",
                     support_assign_article_default(nodes, node_count),
                     SUPPORT_OK) != 0) {
        return 1;
    }

    if (check_status("cube_generate_dirichlet_points",
                     cube_generate_dirichlet_points(CUBE_L,
                                                    CUBE_V0,
                                                    NODE_NX,
                                                    NODE_NY,
                                                    NODE_NZ,
                                                    points,
                                                    CONSTRAINT_COUNT,
                                                    &constraint_count),
                     CUBE_PROBLEM_OK) != 0) {
        return 1;
    }

    if (constraint_count != CONSTRAINT_COUNT) {
        fprintf(stderr, "Unexpected constraint count: %d\n", constraint_count);
        return 1;
    }

    if (check_status("global_stiffness_assemble_dense",
                     global_stiffness_assemble_dense(nodes,
                                                     node_count,
                                                     0.0,
                                                     CUBE_L,
                                                     0.0,
                                                     CUBE_L,
                                                     0.0,
                                                     CUBE_L,
                                                     CELL_COUNT,
                                                     CELL_COUNT,
                                                     CELL_COUNT,
                                                     K),
                     GLOBAL_STIFFNESS_OK) != 0) {
        return 1;
    }

    for (int i = 0; i < NODE_COUNT; ++i) {
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
        return 1;
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
        return 1;
    }

    if (check_status("dense_solve",
                     dense_solve(A_aug, b_aug, TOTAL_SIZE, solution),
                     DENSE_SOLVER_OK) != 0) {
        return 1;
    }

    file = fopen(output_path, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open output CSV: %s\n", output_path);
        return 1;
    }

    fprintf(file, "x,y,z,V_eval,V_exact,abs_error,nearest_node_index,V_coeff_nearest\n");

    for (int ix = 0; ix < SAMPLE_NX; ++ix) {
        const double x = (CUBE_L * (double)ix) / (double)(SAMPLE_NX - 1);

        for (int iy = 0; iy < SAMPLE_NY; ++iy) {
            const double y = (CUBE_L * (double)iy) / (double)(SAMPLE_NY - 1);

            for (int iz = 0; iz < SAMPLE_NZ; ++iz) {
                const double z = (CUBE_L * (double)iz) / (double)(SAMPLE_NZ - 1);
                const int nearest_index = nearest_node_index(nodes, x, y, z);
                double reconstructed = 0.0;
                double exact = 0.0;
                double abs_error = 0.0;

                if (evaluate_mls_potential(nodes, solution, x, y, z, &reconstructed) != 0) {
                    fclose(file);
                    return 1;
                }

                exact = analytical_potential_cube(x,
                                                  y,
                                                  z,
                                                  CUBE_L,
                                                  CUBE_V0,
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
        fprintf(stderr, "Could not close output CSV: %s\n", output_path);
        return 1;
    }

    if (check_status("relative_error_norm",
                     relative_error_norm(numerical_samples,
                                         exact_samples,
                                         SAMPLE_COUNT,
                                         &relative_error),
                     RELATIVE_ERROR_NORM_OK) != 0) {
        return 1;
    }

    printf("Cube smoke EFG run completed.\n");
    printf("nodes: %d\n", node_count);
    printf("constraints: %d\n", constraint_count);
    printf("augmented system size: %d\n", TOTAL_SIZE);
    printf("sample points: %d\n", SAMPLE_COUNT);
    printf("max abs error: %.17g\n", max_abs_error);
    printf("relative discrete error: %.17g\n", relative_error);
    printf("csv: %s\n", output_path);

    return 0;
}
