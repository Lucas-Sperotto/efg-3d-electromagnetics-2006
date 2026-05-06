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

static int expect_close_tol(const char *name,
                            double actual,
                            double expected,
                            double tolerance)
{
    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g (tol %.3e)\n",
                name, expected, actual, tolerance);
        return 1;
    }

    return 0;
}

static double exact_one(double x, double y, double z, void *context)
{
    (void)x;
    (void)y;
    (void)z;
    (void)context;
    return 1.0;
}

static double exact_zero(double x, double y, double z, void *context)
{
    (void)x;
    (void)y;
    (void)z;
    (void)context;
    return 0.0;
}

static double exact_linear(double x, double y, double z, void *context)
{
    (void)context;
    return 1.0 + (2.0 * x) - (3.0 * y) + (4.0 * z);
}

static void make_unit_cube_nodes(Node3D nodes[8], double support_radius)
{
    int index = 0;

    for (int ix = 0; ix <= 1; ++ix) {
        for (int iy = 0; iy <= 1; ++iy) {
            for (int iz = 0; iz <= 1; ++iz) {
                nodes[index].x = (double)ix;
                nodes[index].y = (double)iy;
                nodes[index].z = (double)iz;
                nodes[index].support_radius = support_radius;
                ++index;
            }
        }
    }
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

static int test_domain_integral_constant_exact_match(void)
{
    Node3D nodes[8];
    double solution[8];
    RelativeErrorDomainResult result;
    int failures = 0;

    make_unit_cube_nodes(nodes, 2.0);
    for (int i = 0; i < 8; ++i) {
        solution[i] = 1.0;
    }

    failures += expect_code("domain constant match code",
                            relative_error_norm_domain_integral(
                                nodes, 8, solution,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                1, 1, 1,
                                exact_one, NULL, &result),
                            RELATIVE_ERROR_NORM_OK);
    failures += expect_close_tol("domain constant match error",
                                 result.error, 0.0, 1.0e-12);
    failures += expect_close_tol("domain constant match numerator",
                                 result.numerator_integral, 0.0, 1.0e-12);
    failures += expect_close_tol("domain constant match denominator",
                                 result.denominator_integral, 1.0, 1.0e-12);

    return failures;
}

static int test_domain_integral_constant_unit_error(void)
{
    Node3D nodes[8];
    double solution[8];
    RelativeErrorDomainResult result;
    int failures = 0;

    make_unit_cube_nodes(nodes, 2.0);
    for (int i = 0; i < 8; ++i) {
        solution[i] = 2.0;
    }

    failures += expect_code("domain constant unit error code",
                            relative_error_norm_domain_integral(
                                nodes, 8, solution,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                1, 1, 1,
                                exact_one, NULL, &result),
                            RELATIVE_ERROR_NORM_OK);
    failures += expect_close_tol("domain constant unit error",
                                 result.error, 1.0, 1.0e-12);

    return failures;
}

static int test_domain_integral_linear_reproduction(void)
{
    Node3D nodes[8];
    double solution[8];
    RelativeErrorDomainResult result;
    int failures = 0;

    make_unit_cube_nodes(nodes, 2.0);
    for (int i = 0; i < 8; ++i) {
        solution[i] = exact_linear(nodes[i].x, nodes[i].y, nodes[i].z, NULL);
    }

    failures += expect_code("domain linear reproduction code",
                            relative_error_norm_domain_integral(
                                nodes, 8, solution,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                2, 2, 2,
                                exact_linear, NULL, &result),
                            RELATIVE_ERROR_NORM_OK);
    failures += expect_close_tol("domain linear reproduction error",
                                 result.error, 0.0, 1.0e-10);
    if (result.quadrature_points != 64) {
        fprintf(stderr, "domain linear reproduction: expected 64 points, got %d\n",
                result.quadrature_points);
        ++failures;
    }

    return failures;
}

static int test_domain_integral_zero_denominator(void)
{
    Node3D nodes[8];
    double solution[8];
    RelativeErrorDomainResult result;
    int failures = 0;

    make_unit_cube_nodes(nodes, 2.0);
    for (int i = 0; i < 8; ++i) {
        solution[i] = 1.0;
    }

    failures += expect_code("domain zero denominator code",
                            relative_error_norm_domain_integral(
                                nodes, 8, solution,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                1, 1, 1,
                                exact_zero, NULL, &result),
                            RELATIVE_ERROR_NORM_ZERO_DENOMINATOR);

    return failures;
}

static int test_domain_integral_mls_failure(void)
{
    Node3D node = {0.0, 0.0, 0.0, 1.0};
    double solution = 1.0;
    RelativeErrorDomainResult result;
    int failures = 0;

    failures += expect_code("domain MLS failure code",
                            relative_error_norm_domain_integral(
                                &node, 1, &solution,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                1, 1, 1,
                                exact_one, NULL, &result),
                            RELATIVE_ERROR_NORM_MLS_FAILED);
    if (result.mls_failures != 1) {
        fprintf(stderr, "domain MLS failure: expected 1 failure, got %d\n",
                result.mls_failures);
        ++failures;
    }

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
    failures += test_domain_integral_constant_exact_match();
    failures += test_domain_integral_constant_unit_error();
    failures += test_domain_integral_linear_reproduction();
    failures += test_domain_integral_zero_denominator();
    failures += test_domain_integral_mls_failure();

    return failures == 0 ? 0 : 1;
}
