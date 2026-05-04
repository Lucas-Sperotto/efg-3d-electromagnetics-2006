#include "error_norm.h"

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
    const double tolerance = 1.0e-12;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int test_equal_vectors_have_zero_error(void)
{
    const double numerical[] = {1.0, 2.0, 3.0};
    const double exact[] = {1.0, 2.0, 3.0};
    double error = -1.0;
    int failures = 0;

    failures += expect_code("equal vectors code",
                            relative_error_norm(numerical, exact, 3, &error),
                            RELATIVE_ERROR_NORM_OK);
    failures += expect_close("equal vectors error", error, 0.0);

    return failures;
}

static int test_known_small_vector_error(void)
{
    const double numerical[] = {2.0, 2.0, 2.0};
    const double exact[] = {1.0, 2.0, 2.0};
    double error = -1.0;
    int failures = 0;

    failures += expect_code("known vector code",
                            relative_error_norm(numerical, exact, 3, &error),
                            RELATIVE_ERROR_NORM_OK);
    failures += expect_close("known vector error", error, 1.0 / 3.0);

    return failures;
}

static int test_zero_denominator(void)
{
    const double numerical[] = {1.0, -1.0};
    const double exact[] = {0.0, 0.0};
    double error = 42.0;
    int failures = 0;

    failures += expect_code("zero denominator code",
                            relative_error_norm(numerical, exact, 2, &error),
                            RELATIVE_ERROR_NORM_ZERO_DENOMINATOR);
    failures += expect_close("zero denominator leaves output unchanged", error, 42.0);

    return failures;
}

static int test_invalid_pointers(void)
{
    const double numerical[] = {1.0};
    const double exact[] = {1.0};
    double error = 0.0;
    int failures = 0;

    failures += expect_code("null numerical",
                            relative_error_norm(NULL, exact, 1, &error),
                            RELATIVE_ERROR_NORM_INVALID_ARGUMENT);
    failures += expect_code("null exact",
                            relative_error_norm(numerical, NULL, 1, &error),
                            RELATIVE_ERROR_NORM_INVALID_ARGUMENT);
    failures += expect_code("null error",
                            relative_error_norm(numerical, exact, 1, NULL),
                            RELATIVE_ERROR_NORM_INVALID_ARGUMENT);

    return failures;
}

static int test_invalid_n(void)
{
    const double numerical[] = {1.0};
    const double exact[] = {1.0};
    double error = 0.0;
    int failures = 0;

    failures += expect_code("n zero",
                            relative_error_norm(numerical, exact, 0, &error),
                            RELATIVE_ERROR_NORM_INVALID_ARGUMENT);
    failures += expect_code("n negative",
                            relative_error_norm(numerical, exact, -1, &error),
                            RELATIVE_ERROR_NORM_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_equal_vectors_have_zero_error();
    failures += test_known_small_vector_error();
    failures += test_zero_denominator();
    failures += test_invalid_pointers();
    failures += test_invalid_n();

    return failures == 0 ? 0 : 1;
}
