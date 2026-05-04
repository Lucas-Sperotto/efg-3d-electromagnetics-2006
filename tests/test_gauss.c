#include "gauss.h"

#include <math.h>
#include <stdio.h>

typedef double (*GaussTestFunction)(double x, double y, double z);

static double function_one(double x, double y, double z)
{
    (void)x;
    (void)y;
    (void)z;
    return 1.0;
}

static double function_x(double x, double y, double z)
{
    (void)y;
    (void)z;
    return x;
}

static double function_y(double x, double y, double z)
{
    (void)x;
    (void)z;
    return y;
}

static double function_z(double x, double y, double z)
{
    (void)x;
    (void)y;
    return z;
}

static double function_x_times_y(double x, double y, double z)
{
    (void)z;
    return x * y;
}

static double integrate_unit_cube(GaussTestFunction function)
{
    GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS];
    double integral = 0.0;

    gauss_legendre_order2_cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, points);

    for (int i = 0; i < GAUSS_LEGENDRE_ORDER2_CUBE_POINTS; ++i) {
        integral += points[i].weight * function(points[i].x, points[i].y, points[i].z);
    }

    return integral;
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

int main(void)
{
    int failures = 0;

    failures += expect_close("integral 1", integrate_unit_cube(function_one), 1.0);
    failures += expect_close("integral x", integrate_unit_cube(function_x), 0.5);
    failures += expect_close("integral y", integrate_unit_cube(function_y), 0.5);
    failures += expect_close("integral z", integrate_unit_cube(function_z), 0.5);
    failures += expect_close("integral x*y", integrate_unit_cube(function_x_times_y), 0.25);

    return failures == 0 ? 0 : 1;
}
