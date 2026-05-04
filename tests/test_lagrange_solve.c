#include "dense_solver.h"
#include "lagrange_system.h"

#include <math.h>
#include <stdio.h>

#define TEST_NODE_COUNT 3
#define TEST_CONSTRAINT_COUNT 2
#define TEST_TOTAL_SIZE (TEST_NODE_COUNT + TEST_CONSTRAINT_COUNT)
#define TEST_AUG_SIZE (TEST_TOTAL_SIZE * TEST_TOTAL_SIZE)

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
    const double tolerance = 1.0e-10;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static void make_manual_blocks(double K[TEST_NODE_COUNT * TEST_NODE_COUNT],
                               double F[TEST_NODE_COUNT],
                               double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT],
                               double q[TEST_CONSTRAINT_COUNT])
{
    const double K_values[TEST_NODE_COUNT * TEST_NODE_COUNT] = {
        2.0, -1.0, 0.0,
        -1.0, 2.0, -1.0,
        0.0, -1.0, 2.0
    };
    const double G_values[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT] = {
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0
    };

    for (int i = 0; i < TEST_NODE_COUNT * TEST_NODE_COUNT; ++i) {
        K[i] = K_values[i];
    }

    for (int i = 0; i < TEST_NODE_COUNT; ++i) {
        F[i] = 0.0;
    }

    for (int i = 0; i < TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT; ++i) {
        G[i] = G_values[i];
    }

    q[0] = 0.0;
    q[1] = 10.0;
}

static void multiply_dense(const double A[TEST_AUG_SIZE],
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
 * Minimal integrated check for the direct didactic path:
 * lagrange_assemble_augmented_system_dense() -> dense_solve().
 *
 * This test intentionally uses manual K, F, G and q blocks. It does not touch
 * the real EFG stiffness assembly, MLS, Dirichlet point assembly, GMRES, or
 * analytical solution code.
 */
static int test_lagrange_assembly_and_dense_solve(void)
{
    double K[TEST_NODE_COUNT * TEST_NODE_COUNT];
    double F[TEST_NODE_COUNT];
    double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_CONSTRAINT_COUNT];
    double A_aug[TEST_AUG_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    double solution[TEST_TOTAL_SIZE];
    double residual_check[TEST_TOTAL_SIZE];
    const double expected_solution[TEST_TOTAL_SIZE] = {
        0.0, 5.0, 10.0, 5.0, -15.0
    };
    int failures = 0;

    make_manual_blocks(K, F, G, q);

    failures += expect_code("assemble augmented system",
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

    multiply_dense(A_aug, solution, residual_check);
    for (int i = 0; i < TEST_TOTAL_SIZE; ++i) {
        failures += expect_close("A_aug*x residual", residual_check[i], b_aug[i]);
        failures += expect_close("expected solution", solution[i], expected_solution[i]);
    }

    return failures;
}

int main(void)
{
    const int failures = test_lagrange_assembly_and_dense_solve();

    return failures == 0 ? 0 : 1;
}
