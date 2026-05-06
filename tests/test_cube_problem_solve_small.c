#include "analytical.h"
#include "cube_problem.h"
#include "dense_solver.h"
#include "dirichlet.h"
#include "error_norm.h"
#include "global_stiffness.h"
#include "lagrange_system.h"
#include "support.h"

#include <math.h>
#include <stdio.h>

#define TEST_L 10.0
#define TEST_V0 10.0
#define TEST_NX 3
#define TEST_NY 3
#define TEST_NZ 3
#define TEST_NODE_COUNT 27
#define TEST_CONSTRAINT_COUNT 26
#define TEST_TOTAL_SIZE (TEST_NODE_COUNT + TEST_CONSTRAINT_COUNT)
#define TEST_NODE_MATRIX_SIZE (TEST_NODE_COUNT * TEST_NODE_COUNT)
#define TEST_CONSTRAINT_MATRIX_SIZE (TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT)
#define TEST_AUG_MATRIX_SIZE (TEST_TOTAL_SIZE * TEST_TOTAL_SIZE)
#define TEST_ANALYTICAL_TERMS 15
#define TEST_RELATIVE_ERROR_LIMIT 2.0

static int expect_code(const char *name, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s: expected code %d, got %d\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int expect_close(const char *name, double actual, double expected, double tolerance)
{
    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int same_point(double ax,
                      double ay,
                      double az,
                      double bx,
                      double by,
                      double bz)
{
    const double tolerance = 1.0e-10;

    return fabs(ax - bx) <= tolerance &&
           fabs(ay - by) <= tolerance &&
           fabs(az - bz) <= tolerance;
}

static int find_node_index(const Node3D nodes[TEST_NODE_COUNT],
                           int node_count,
                           double x,
                           double y,
                           double z)
{
    for (int i = 0; i < node_count; ++i) {
        if (same_point(nodes[i].x, nodes[i].y, nodes[i].z, x, y, z)) {
            return i;
        }
    }

    return -1;
}

static double apply_constraint_row(const double G[TEST_CONSTRAINT_MATRIX_SIZE],
                                   const double V[TEST_NODE_COUNT],
                                   int row)
{
    double value = 0.0;

    for (int node = 0; node < TEST_NODE_COUNT; ++node) {
        value += G[(row * TEST_NODE_COUNT) + node] * V[node];
    }

    return value;
}

static void multiply_augmented(const double A[TEST_AUG_MATRIX_SIZE],
                               const double x[TEST_TOTAL_SIZE],
                               double result[TEST_TOTAL_SIZE])
{
    for (int row = 0; row < TEST_TOTAL_SIZE; ++row) {
        result[row] = 0.0;

        for (int col = 0; col < TEST_TOTAL_SIZE; ++col) {
            result[row] += A[(row * TEST_TOTAL_SIZE) + col] * x[col];
        }
    }
}

/*
 * Small integrated smoke test for the electrostatic cube pipeline.
 *
 * This is not the final reproduction of Figures 3, 4, or 5 of the article.
 * It only checks that the small regular cube geometry, support assignment,
 * stiffness assembly, Dirichlet constraints, Lagrange augmentation, dense
 * direct solve, analytical reference, and relative error helper interoperate.
 * GMRES remains intentionally outside this test.
 */
static int test_small_cube_pipeline(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    DirichletPoint points[TEST_CONSTRAINT_COUNT];
    double K[TEST_NODE_MATRIX_SIZE];
    double G[TEST_CONSTRAINT_MATRIX_SIZE];
    double q[TEST_CONSTRAINT_COUNT];
    double F[TEST_NODE_COUNT];
    double A_aug[TEST_AUG_MATRIX_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    double solution[TEST_TOTAL_SIZE];
    double residual[TEST_TOTAL_SIZE];
    double exact_at_nodes[TEST_NODE_COUNT];
    double interior_numerical[1];
    double interior_exact[1];
    double relative_error = 0.0;
    int node_count = 0;
    int point_count = 0;
    int center_index = -1;
    int failures = 0;

    failures += expect_code("generate cube nodes",
                            cube_generate_regular_nodes(TEST_L,
                                                        TEST_NX,
                                                        TEST_NY,
                                                        TEST_NZ,
                                                        nodes,
                                                        TEST_NODE_COUNT,
                                                        &node_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("generated node count", node_count, TEST_NODE_COUNT);
    if (failures != 0) {
        return failures;
    }

    failures += expect_code("assign article support",
                            support_assign_article_default(nodes, node_count),
                            SUPPORT_OK);
    if (failures != 0) {
        return failures;
    }

    for (int node = 0; node < node_count; ++node) {
        if (nodes[node].support_radius <= 0.0) {
            fprintf(stderr, "node %d has non-positive support radius\n", node);
            ++failures;
        }
    }

    failures += expect_code("generate Dirichlet points",
                            cube_generate_dirichlet_points(TEST_L,
                                                           TEST_V0,
                                                           TEST_NX,
                                                           TEST_NY,
                                                           TEST_NZ,
                                                           points,
                                                           TEST_CONSTRAINT_COUNT,
                                                           &point_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("generated Dirichlet point count",
                            point_count,
                            TEST_CONSTRAINT_COUNT);

    for (int i = 0; i < TEST_NODE_COUNT; ++i) {
        F[i] = 0.0;
    }

    failures += expect_code("global stiffness",
                            global_stiffness_assemble_dense(nodes,
                                                            node_count,
                                                            0.0,
                                                            TEST_L,
                                                            0.0,
                                                            TEST_L,
                                                            0.0,
                                                            TEST_L,
                                                            2,
                                                            2,
                                                            2,
                                                            K),
                            GLOBAL_STIFFNESS_OK);
    failures += expect_code("Dirichlet constraints",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 node_count,
                                                                 points,
                                                                 point_count,
                                                                 G,
                                                                 q),
                            DIRICHLET_OK);
    failures += expect_code("Lagrange augmented system",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     node_count,
                                                                     G,
                                                                     q,
                                                                     point_count,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_OK);
    if (failures != 0) {
        return failures;
    }

    failures += expect_code("dense solve",
                            dense_solve(A_aug, b_aug, TEST_TOTAL_SIZE, solution),
                            DENSE_SOLVER_OK);
    if (failures != 0) {
        return failures;
    }

    for (int i = 0; i < TEST_TOTAL_SIZE; ++i) {
        if (!isfinite(solution[i])) {
            fprintf(stderr, "solution[%d] is not finite\n", i);
            ++failures;
        }
    }

    /*
     * The first node_count unknowns are kept finite here, but this initial
     * smoke test does not enforce a strict [0,V0] bound on them. In EFG the
     * unknown coefficients are not Kronecker-nodal values, and this very coarse
     * 3x3x3 system shows small boundary overshoots while still satisfying
     * G * V = q and the augmented residual.
     */

    for (int constraint = 0; constraint < point_count; ++constraint) {
        char label[96];
        const double constrained_value = apply_constraint_row(G, solution, constraint);

        snprintf(label, sizeof(label), "Dirichlet constraint[%d]", constraint);
        failures += expect_close(label, constrained_value, q[constraint], 1.0e-7);
    }

    multiply_augmented(A_aug, solution, residual);
    for (int row = 0; row < TEST_TOTAL_SIZE; ++row) {
        char label[96];

        snprintf(label, sizeof(label), "augmented residual[%d]", row);
        failures += expect_close(label, residual[row], b_aug[row], 1.0e-7);
    }

    for (int node = 0; node < node_count; ++node) {
        exact_at_nodes[node] = analytical_potential_cube(nodes[node].x,
                                                        nodes[node].y,
                                                        nodes[node].z,
                                                        TEST_L,
                                                        TEST_V0,
                                                        TEST_ANALYTICAL_TERMS);
        if (!isfinite(exact_at_nodes[node])) {
            fprintf(stderr, "analytical value at node %d is not finite\n", node);
            ++failures;
        }
    }

    center_index = find_node_index(nodes, node_count, 5.0, 5.0, 5.0);
    if (center_index < 0) {
        fprintf(stderr, "center node (5,5,5) was not generated\n");
        return failures + 1;
    }

    interior_numerical[0] = solution[center_index];
    interior_exact[0] = exact_at_nodes[center_index];

    failures += expect_code("interior relative error",
                            relative_error_norm(interior_numerical,
                                                interior_exact,
                                                1,
                                                &relative_error),
                            RELATIVE_ERROR_NORM_OK);
    if (!isfinite(relative_error)) {
        fprintf(stderr, "relative error is not finite\n");
        ++failures;
    } else if (relative_error >= TEST_RELATIVE_ERROR_LIMIT) {
        fprintf(stderr, "relative error is too large: %.17g\n", relative_error);
        ++failures;
    }

    if (!isfinite(interior_numerical[0]) || !isfinite(interior_exact[0])) {
        fprintf(stderr, "center-node numerical or analytical value is not finite\n");
        ++failures;
    } else {
        failures += expect_close("center node loose analytical comparison",
                                 interior_numerical[0],
                                 interior_exact[0],
                                 5.0);
    }

    /*
     * The full-node analytical vector is computed above, but the relative norm
     * is evaluated on the single interior node for this first smoke test. With
     * this tiny 3x3x3 grid, boundary values are too coarse to serve as a
     * meaningful convergence metric. The single interior sample keeps the test
     * focused on pipeline interoperability rather than boundary singularities.
     */

    return failures;
}

int main(void)
{
    const int failures = test_small_cube_pipeline();

    return failures == 0 ? 0 : 1;
}
