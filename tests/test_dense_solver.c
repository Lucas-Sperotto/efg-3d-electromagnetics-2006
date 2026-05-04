#include "dense_solver.h"

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

static int expect_close(const char *name, double actual, double expected, double tolerance)
{
    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static void multiply_dense(const double *A, const double *x, int n, double *result)
{
    for (int row = 0; row < n; ++row) {
        result[row] = 0.0;
        for (int col = 0; col < n; ++col) {
            result[row] += A[(row * n) + col] * x[col];
        }
    }
}

static int expect_residual(const char *name,
                           const double *A,
                           const double *x,
                           const double *b,
                           int n)
{
    double result[5];
    int failures = 0;

    multiply_dense(A, x, n, result);

    for (int i = 0; i < n; ++i) {
        failures += expect_close(name, result[i], b[i], 1.0e-10);
    }

    return failures;
}

static int test_one_by_one(void)
{
    const double A[1] = {4.0};
    const double b[1] = {8.0};
    double x[1] = {0.0};
    int failures = 0;

    failures += expect_code("1x1 solve", dense_solve(A, b, 1, x), DENSE_SOLVER_OK);
    failures += expect_close("1x1 x[0]", x[0], 2.0, 1.0e-12);

    return failures;
}

static int test_two_by_two(void)
{
    const double A[4] = {
        3.0, 2.0,
        1.0, 2.0
    };
    const double b[2] = {5.0, 5.0};
    double x[2] = {0.0, 0.0};
    int failures = 0;

    failures += expect_code("2x2 solve", dense_solve(A, b, 2, x), DENSE_SOLVER_OK);
    failures += expect_close("2x2 x[0]", x[0], 0.0, 1.0e-12);
    failures += expect_close("2x2 x[1]", x[1], 2.5, 1.0e-12);

    return failures;
}

static int test_three_by_three(void)
{
    const double A[9] = {
        4.0, -1.0, 0.0,
        -1.0, 4.0, -1.0,
        0.0, -1.0, 3.0
    };
    const double b[3] = {2.0, 4.0, 7.0};
    double x[3] = {0.0, 0.0, 0.0};
    int failures = 0;

    failures += expect_code("3x3 solve", dense_solve(A, b, 3, x), DENSE_SOLVER_OK);
    failures += expect_close("3x3 x[0]", x[0], 1.0, 1.0e-12);
    failures += expect_close("3x3 x[1]", x[1], 2.0, 1.0e-12);
    failures += expect_close("3x3 x[2]", x[2], 3.0, 1.0e-12);

    return failures;
}

static int test_lagrange_five_by_five_residual(void)
{
    const double A[25] = {
        2.0, -1.0, 0.0, 1.0, 0.0,
        -1.0, 2.0, -1.0, 0.0, 0.0,
        0.0, -1.0, 2.0, 0.0, 1.0,
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0
    };
    const double b[5] = {0.0, 0.0, 0.0, 0.0, 10.0};
    double x[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int failures = 0;

    failures += expect_code("5x5 Lagrange solve", dense_solve(A, b, 5, x), DENSE_SOLVER_OK);
    if (failures != 0) {
        return failures;
    }

    failures += expect_residual("5x5 Lagrange residual", A, x, b, 5);

    return failures;
}

static int test_singular_and_near_singular(void)
{
    const double singular_A[4] = {
        1.0, 2.0,
        2.0, 4.0
    };
    const double singular_b[2] = {3.0, 6.0};
    const double near_A[1] = {1.0e-15};
    const double near_b[1] = {1.0};
    double x[2] = {0.0, 0.0};
    int failures = 0;

    failures += expect_code("singular matrix",
                            dense_solve(singular_A, singular_b, 2, x),
                            DENSE_SOLVER_SINGULAR);
    failures += expect_code("near singular matrix",
                            dense_solve(near_A, near_b, 1, x),
                            DENSE_SOLVER_NEAR_SINGULAR);

    return failures;
}

static int test_invalid_arguments(void)
{
    const double A[1] = {1.0};
    const double b[1] = {1.0};
    double x[1] = {0.0};
    int failures = 0;

    failures += expect_code("null A", dense_solve(NULL, b, 1, x), DENSE_SOLVER_INVALID_ARGUMENT);
    failures += expect_code("null b", dense_solve(A, NULL, 1, x), DENSE_SOLVER_INVALID_ARGUMENT);
    failures += expect_code("null x", dense_solve(A, b, 1, NULL), DENSE_SOLVER_INVALID_ARGUMENT);
    failures += expect_code("invalid n", dense_solve(A, b, 0, x), DENSE_SOLVER_INVALID_ARGUMENT);

    return failures;
}

static int test_preserves_inputs(void)
{
    double A[9] = {
        4.0, -1.0, 0.0,
        -1.0, 4.0, -1.0,
        0.0, -1.0, 3.0
    };
    double b[3] = {2.0, 4.0, 7.0};
    const double original_A[9] = {
        4.0, -1.0, 0.0,
        -1.0, 4.0, -1.0,
        0.0, -1.0, 3.0
    };
    const double original_b[3] = {2.0, 4.0, 7.0};
    double x[3] = {0.0, 0.0, 0.0};
    int failures = 0;

    failures += expect_code("preservation solve", dense_solve(A, b, 3, x), DENSE_SOLVER_OK);

    for (int i = 0; i < 9; ++i) {
        failures += expect_close("A preserved", A[i], original_A[i], 0.0);
    }

    for (int i = 0; i < 3; ++i) {
        failures += expect_close("b preserved", b[i], original_b[i], 0.0);
    }

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_one_by_one();
    failures += test_two_by_two();
    failures += test_three_by_three();
    failures += test_lagrange_five_by_five_residual();
    failures += test_singular_and_near_singular();
    failures += test_invalid_arguments();
    failures += test_preserves_inputs();

    return failures == 0 ? 0 : 1;
}
