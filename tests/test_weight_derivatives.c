#include "weight.h"

#include <math.h>
#include <stdio.h>

static int expect_close(const char *name, double actual, double expected, double tolerance)
{
    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static double finite_difference(double (*weight)(double), double r)
{
    const double h = 1.0e-6;

    return (weight(r + h) - weight(r - h)) / (2.0 * h);
}

static int test_cubic_exact_values(void)
{
    int failures = 0;

    failures += expect_close("cubic r < 0", weight_cubic_spline_derivative(-0.25), 0.0, 1.0e-12);
    failures += expect_close("cubic r = 0", weight_cubic_spline_derivative(0.0), 0.0, 1.0e-12);
    failures += expect_close("cubic first interval", weight_cubic_spline_derivative(0.25), -1.25, 1.0e-12);
    failures += expect_close("cubic change point", weight_cubic_spline_derivative(0.5), -1.0, 1.0e-12);
    failures += expect_close("cubic second interval", weight_cubic_spline_derivative(0.75), -0.25, 1.0e-12);
    failures += expect_close("cubic r = 1", weight_cubic_spline_derivative(1.0), 0.0, 1.0e-12);
    failures += expect_close("cubic r > 1", weight_cubic_spline_derivative(1.25), 0.0, 1.0e-12);

    return failures;
}

static int test_quadratic_exact_values(void)
{
    int failures = 0;

    failures += expect_close("quadratic r < 0", weight_quadratic_derivative(-0.25), 0.0, 1.0e-12);
    failures += expect_close("quadratic r = 0", weight_quadratic_derivative(0.0), 0.0, 1.0e-12);
    failures += expect_close("quadratic first interval", weight_quadratic_derivative(0.25), -1.0, 1.0e-12);
    failures += expect_close("quadratic change point", weight_quadratic_derivative(0.5), -2.0, 1.0e-12);
    failures += expect_close("quadratic second interval", weight_quadratic_derivative(0.75), -1.0, 1.0e-12);
    failures += expect_close("quadratic r = 1", weight_quadratic_derivative(1.0), 0.0, 1.0e-12);
    failures += expect_close("quadratic r > 1", weight_quadratic_derivative(1.25), 0.0, 1.0e-12);

    return failures;
}

static int test_cubic_finite_difference(void)
{
    int failures = 0;
    const double points[] = {0.10, 0.25, 0.70, 0.90};

    for (int i = 0; i < 4; ++i) {
        const double r = points[i];

        failures += expect_close("cubic finite difference",
                                 weight_cubic_spline_derivative(r),
                                 finite_difference(weight_cubic_spline, r),
                                 1.0e-8);
    }

    return failures;
}

static int test_quadratic_finite_difference(void)
{
    int failures = 0;
    const double points[] = {0.10, 0.25, 0.70, 0.90};

    for (int i = 0; i < 4; ++i) {
        const double r = points[i];

        failures += expect_close("quadratic finite difference",
                                 weight_quadratic_derivative(r),
                                 finite_difference(weight_quadratic, r),
                                 1.0e-8);
    }

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_cubic_exact_values();
    failures += test_quadratic_exact_values();
    failures += test_cubic_finite_difference();
    failures += test_quadratic_finite_difference();

    return failures == 0 ? 0 : 1;
}
