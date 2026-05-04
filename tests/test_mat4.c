#include "mat4.h"

#include <math.h>
#include <stdio.h>

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
    const double tolerance = 1.0e-11;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int expect_vector_close(const char *name, const double actual[4], const double expected[4])
{
    int failures = 0;

    for (int i = 0; i < 4; ++i) {
        char label[64];

        snprintf(label, sizeof(label), "%s[%d]", name, i);
        failures += expect_close(label, actual[i], expected[i]);
    }

    return failures;
}

static int expect_matrix_unchanged(const double actual[4][4], const double expected[4][4])
{
    int failures = 0;

    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            char label[64];

            snprintf(label, sizeof(label), "A unchanged[%d][%d]", row, col);
            failures += expect_close(label, actual[row][col], expected[row][col]);
        }
    }

    return failures;
}

static int test_zero_identity_and_copy(void)
{
    double zero[4][4];
    double identity[4][4];
    double copy[4][4];
    int failures = 0;

    mat4_zero(zero);
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            failures += expect_close("zero", zero[row][col], 0.0);
        }
    }

    mat4_identity(identity);
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            failures += expect_close("identity", identity[row][col], row == col ? 1.0 : 0.0);
        }
    }

    mat4_copy((const double(*)[4])identity, copy);
    failures += expect_matrix_unchanged((const double(*)[4])copy, (const double(*)[4])identity);

    return failures;
}

static int test_identity_solve(void)
{
    double A[4][4];
    const double b[4] = {3.0, -2.0, 7.0, 0.5};
    const double expected[4] = {3.0, -2.0, 7.0, 0.5};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    mat4_identity(A);
    failures += expect_code("identity solve code", mat4_solve((const double(*)[4])A, b, x), MAT4_OK);
    failures += expect_vector_close("identity solve", x, expected);

    return failures;
}

static int test_diagonal_solve(void)
{
    double A[4][4];
    const double b[4] = {2.0, 6.0, -12.0, 20.0};
    const double expected[4] = {1.0, 2.0, -3.0, 4.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    mat4_zero(A);
    A[0][0] = 2.0;
    A[1][1] = 3.0;
    A[2][2] = 4.0;
    A[3][3] = 5.0;

    failures += expect_code("diagonal solve code", mat4_solve((const double(*)[4])A, b, x), MAT4_OK);
    failures += expect_vector_close("diagonal solve", x, expected);

    return failures;
}

static int test_known_4x4_solve(void)
{
    const double A[4][4] = {
        {4.0, 1.0, 2.0, 0.0},
        {1.0, 3.0, 0.0, 1.0},
        {2.0, 0.0, 5.0, 1.0},
        {0.0, 1.0, 1.0, 2.0},
    };
    const double b[4] = {12.0, 11.0, 21.0, 13.0};
    const double expected[4] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    failures += expect_code("known solve code", mat4_solve(A, b, x), MAT4_OK);
    failures += expect_vector_close("known solve", x, expected);

    return failures;
}

static int test_singular_and_near_singular(void)
{
    const double singular[4][4] = {
        {1.0, 2.0, 3.0, 4.0},
        {1.0, 2.0, 3.0, 4.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
    };
    const double near_singular[4][4] = {
        {1.0e-14, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
    };
    const double b[4] = {1.0, 1.0, 1.0, 1.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    failures += expect_code("singular solve", mat4_solve(singular, b, x), MAT4_SINGULAR);
    failures += expect_code("near singular solve",
                            mat4_solve(near_singular, b, x),
                            MAT4_NEAR_SINGULAR);

    return failures;
}

static int test_inputs_are_not_modified(void)
{
    const double A_original[4][4] = {
        {2.0, 1.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0},
        {0.0, 1.0, 2.0, 1.0},
        {0.0, 0.0, 1.0, 2.0},
    };
    const double b_original[4] = {1.0, 2.0, 3.0, 4.0};
    double A[4][4];
    double b[4] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    mat4_copy(A_original, A);

    failures += expect_code("preserve solve code", mat4_solve((const double(*)[4])A, b, x), MAT4_OK);
    failures += expect_matrix_unchanged((const double(*)[4])A, A_original);
    failures += expect_vector_close("b unchanged", b, b_original);

    return failures;
}

static int test_invalid_arguments(void)
{
    double A[4][4];
    double b[4] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    mat4_identity(A);
    failures += expect_code("null A", mat4_solve(NULL, b, x), MAT4_INVALID_ARGUMENT);
    failures += expect_code("null b", mat4_solve((const double(*)[4])A, NULL, x), MAT4_INVALID_ARGUMENT);
    failures += expect_code("null x", mat4_solve((const double(*)[4])A, b, NULL), MAT4_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_zero_identity_and_copy();
    failures += test_identity_solve();
    failures += test_diagonal_solve();
    failures += test_known_4x4_solve();
    failures += test_singular_and_near_singular();
    failures += test_inputs_are_not_modified();
    failures += test_invalid_arguments();

    return failures == 0 ? 0 : 1;
}
