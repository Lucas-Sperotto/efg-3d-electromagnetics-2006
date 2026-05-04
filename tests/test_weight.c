#include "weight.h"

#include <stdio.h>

static double abs_double(double value)
{
    return value < 0.0 ? -value : value;
}

static int almost_equal(double actual, double expected)
{
    const double tolerance = 1.0e-12;
    return abs_double(actual - expected) <= tolerance;
}

static int expect_weight(const char *name, double actual, double expected)
{
    if (!almost_equal(actual, expected)) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

int main(void)
{
    int failures = 0;

    failures += expect_weight("cubic r < 0", weight_cubic_spline(-0.25), 0.0);
    failures += expect_weight("cubic r = 0", weight_cubic_spline(0.0), 2.0 / 3.0);
    failures += expect_weight("cubic r = 0.25", weight_cubic_spline(0.25), 23.0 / 48.0);
    failures += expect_weight("cubic r = 0.75", weight_cubic_spline(0.75), 1.0 / 48.0);
    failures += expect_weight("cubic r = 1", weight_cubic_spline(1.0), 0.0);
    failures += expect_weight("cubic r > 1", weight_cubic_spline(1.25), 0.0);

    failures += expect_weight("quadratic r < 0", weight_quadratic(-0.25), 0.0);
    failures += expect_weight("quadratic r = 0", weight_quadratic(0.0), 1.0);
    failures += expect_weight("quadratic r = 0.25", weight_quadratic(0.25), 7.0 / 8.0);
    failures += expect_weight("quadratic r = 0.75", weight_quadratic(0.75), 1.0 / 8.0);
    failures += expect_weight("quadratic r = 1", weight_quadratic(1.0), 0.0);
    failures += expect_weight("quadratic r > 1", weight_quadratic(1.25), 0.0);

    return failures == 0 ? 0 : 1;
}
