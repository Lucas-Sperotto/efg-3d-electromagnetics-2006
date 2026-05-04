#include "dense_solver.h"
#include "dirichlet.h"
#include "global_stiffness.h"
#include "lagrange_system.h"

#include <math.h>
#include <stdio.h>

#define TEST_NODE_COUNT 8
#define TEST_CONSTRAINT_COUNT 8
#define TEST_TOTAL_SIZE (TEST_NODE_COUNT + TEST_CONSTRAINT_COUNT)
#define TEST_NODE_MATRIX_SIZE (TEST_NODE_COUNT * TEST_NODE_COUNT)
#define TEST_AUG_MATRIX_SIZE (TEST_TOTAL_SIZE * TEST_TOTAL_SIZE)

static int expect_code(const char *name, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s: expected code %d, got %d\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int expect_close(const char *name, double actual, double expected)
{
    const double tolerance = 1.0e-8;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static void make_unit_cube_nodes(Node3D nodes[TEST_NODE_COUNT])
{
    const double support_radius = 2.0;

    nodes[0] = (Node3D){0.0, 0.0, 0.0, support_radius};
    nodes[1] = (Node3D){1.0, 0.0, 0.0, support_radius};
    nodes[2] = (Node3D){0.0, 1.0, 0.0, support_radius};
    nodes[3] = (Node3D){1.0, 1.0, 0.0, support_radius};
    nodes[4] = (Node3D){0.0, 0.0, 1.0, support_radius};
    nodes[5] = (Node3D){1.0, 0.0, 1.0, support_radius};
    nodes[6] = (Node3D){0.0, 1.0, 1.0, support_radius};
    nodes[7] = (Node3D){1.0, 1.0, 1.0, support_radius};
}

static void make_vertex_dirichlet_points(const Node3D nodes[TEST_NODE_COUNT],
                                         DirichletPoint points[TEST_CONSTRAINT_COUNT])
{
    for (int i = 0; i < TEST_CONSTRAINT_COUNT; ++i) {
        points[i] = (DirichletPoint){
            nodes[i].x,
            nodes[i].y,
            nodes[i].z,
            10.0 * nodes[i].z
        };
    }
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

static double apply_constraint_row(const double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT],
                                   const double V[TEST_NODE_COUNT],
                                   int row)
{
    double value = 0.0;

    for (int node = 0; node < TEST_NODE_COUNT; ++node) {
        value += G[(row * TEST_NODE_COUNT) + node] * V[node];
    }

    return value;
}

/*
 * Artificial linear patch test for the integrated didactic EFG pipeline.
 *
 * This is not the final reproduction of the article's cube problem. It uses
 * the manufactured linear field V(x,y,z) = 10 z only to validate that dense
 * stiffness assembly, MLS-based Dirichlet constraints, Lagrange augmentation,
 * and the dense direct solver work together on a small case. GMRES is still
 * intentionally outside this test.
 */
static int run_linear_patch_test(int cell_count)
{
    Node3D nodes[TEST_NODE_COUNT];
    DirichletPoint points[TEST_CONSTRAINT_COUNT];
    double K[TEST_NODE_MATRIX_SIZE];
    double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_CONSTRAINT_COUNT];
    double F[TEST_NODE_COUNT];
    double A_aug[TEST_AUG_MATRIX_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    double solution[TEST_TOTAL_SIZE];
    double residual[TEST_TOTAL_SIZE];
    int failures = 0;

    make_unit_cube_nodes(nodes);
    make_vertex_dirichlet_points(nodes, points);

    for (int i = 0; i < TEST_NODE_COUNT; ++i) {
        F[i] = 0.0;
    }

    failures += expect_code("global stiffness",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            cell_count,
                                                            cell_count,
                                                            cell_count,
                                                            K),
                            GLOBAL_STIFFNESS_OK);
    failures += expect_code("Dirichlet constraints",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_CONSTRAINT_COUNT,
                                                                 G,
                                                                 q),
                            DIRICHLET_OK);
    failures += expect_code("Lagrange augmented system",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
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
            fprintf(stderr, "cell_count=%d: solution[%d] is not finite\n", cell_count, i);
            ++failures;
        }
    }

    for (int node = 0; node < TEST_NODE_COUNT; ++node) {
        char label[80];
        const double expected_value = 10.0 * nodes[node].z;

        snprintf(label, sizeof(label), "cell_count=%d nodal V[%d]", cell_count, node);
        failures += expect_close(label, solution[node], expected_value);
    }

    for (int constraint = 0; constraint < TEST_CONSTRAINT_COUNT; ++constraint) {
        char label[96];
        const double constraint_value = apply_constraint_row(G, solution, constraint);

        snprintf(label,
                 sizeof(label),
                 "cell_count=%d constraint[%d]",
                 cell_count,
                 constraint);
        failures += expect_close(label, constraint_value, q[constraint]);
    }

    multiply_augmented(A_aug, solution, residual);
    for (int row = 0; row < TEST_TOTAL_SIZE; ++row) {
        char label[96];

        snprintf(label, sizeof(label), "cell_count=%d augmented residual[%d]", cell_count, row);
        failures += expect_close(label, residual[row], b_aug[row]);
    }

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += run_linear_patch_test(1);
    failures += run_linear_patch_test(2);

    return failures == 0 ? 0 : 1;
}
