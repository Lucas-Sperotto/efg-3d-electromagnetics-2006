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
    const double tolerance = 1.0e-12;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static void make_manual_system(double K[TEST_NODE_COUNT * TEST_NODE_COUNT],
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

static int test_manual_augmented_system(void)
{
    double K[TEST_NODE_COUNT * TEST_NODE_COUNT];
    double F[TEST_NODE_COUNT];
    double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_CONSTRAINT_COUNT];
    double A_aug[TEST_AUG_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    const double expected_A[TEST_AUG_SIZE] = {
        2.0, -1.0, 0.0, 1.0, 0.0,
        -1.0, 2.0, -1.0, 0.0, 0.0,
        0.0, -1.0, 2.0, 0.0, 1.0,
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0
    };
    const double expected_b[TEST_TOTAL_SIZE] = {
        0.0, 0.0, 0.0, 0.0, 10.0
    };
    int failures = 0;

    make_manual_system(K, F, G, q);

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

    for (int i = 0; i < TEST_AUG_SIZE; ++i) {
        failures += expect_close("A_aug manual value", A_aug[i], expected_A[i]);
    }

    for (int i = 0; i < TEST_TOTAL_SIZE; ++i) {
        failures += expect_close("b_aug manual value", b_aug[i], expected_b[i]);
    }

    return failures;
}

static int test_symmetry_and_zero_block(void)
{
    double K[TEST_NODE_COUNT * TEST_NODE_COUNT];
    double F[TEST_NODE_COUNT];
    double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_CONSTRAINT_COUNT];
    double A_aug[TEST_AUG_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    int failures = 0;

    make_manual_system(K, F, G, q);

    failures += expect_code("assemble for structural checks",
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

    for (int row = 0; row < TEST_TOTAL_SIZE; ++row) {
        for (int col = 0; col < TEST_TOTAL_SIZE; ++col) {
            failures += expect_close("A_aug symmetry",
                                     A_aug[(row * TEST_TOTAL_SIZE) + col],
                                     A_aug[(col * TEST_TOTAL_SIZE) + row]);
        }
    }

    for (int row = TEST_NODE_COUNT; row < TEST_TOTAL_SIZE; ++row) {
        for (int col = TEST_NODE_COUNT; col < TEST_TOTAL_SIZE; ++col) {
            failures += expect_close("lower-right zero block",
                                     A_aug[(row * TEST_TOTAL_SIZE) + col],
                                     0.0);
        }
    }

    return failures;
}

static int test_error_paths(void)
{
    double K[TEST_NODE_COUNT * TEST_NODE_COUNT];
    double F[TEST_NODE_COUNT];
    double G[TEST_CONSTRAINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_CONSTRAINT_COUNT];
    double A_aug[TEST_AUG_SIZE];
    double b_aug[TEST_TOTAL_SIZE];
    int failures = 0;

    make_manual_system(K, F, G, q);

    failures += expect_code("null K",
                            lagrange_assemble_augmented_system_dense(NULL,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("null F",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     NULL,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("null G",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     NULL,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("null q",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     NULL,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("null A_aug",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     NULL,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("null b_aug",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     NULL),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("invalid node_count",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     0,
                                                                     G,
                                                                     q,
                                                                     TEST_CONSTRAINT_COUNT,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);
    failures += expect_code("invalid constraint_count",
                            lagrange_assemble_augmented_system_dense(K,
                                                                     F,
                                                                     TEST_NODE_COUNT,
                                                                     G,
                                                                     q,
                                                                     0,
                                                                     A_aug,
                                                                     b_aug),
                            LAGRANGE_SYSTEM_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_manual_augmented_system();
    failures += test_symmetry_and_zero_block();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
