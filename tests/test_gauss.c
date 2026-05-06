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

static double integrate_cell_order(int order,
                                   double xmin,
                                   double xmax,
                                   double ymin,
                                   double ymax,
                                   double zmin,
                                   double zmax,
                                   GaussTestFunction function)
{
    GaussPoint3D points[GAUSS_LEGENDRE_MAX_CUBE_POINTS];
    int point_count = 0;
    double integral = 0.0;

    if (gauss_legendre_cube(order, xmin, xmax, ymin, ymax, zmin, zmax,
                            points, GAUSS_LEGENDRE_MAX_CUBE_POINTS,
                            &point_count) != GAUSS_OK) {
        return NAN;
    }

    for (int i = 0; i < point_count; ++i) {
        integral += points[i].weight
                    * function(points[i].x, points[i].y, points[i].z);
    }

    return integral;
}

static double function_linear(double x, double y, double z)
{
    return 2.0 * x - 3.0 * y + 4.0 * z + 5.0;
}

static double function_x_squared(double x, double y, double z)
{
    (void)y;
    (void)z;
    return x * x;
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

static int test_weight_sum_1d_orders_1_to_8(void)
{
    int failures = 0;

    for (int order = 1; order <= 8; ++order) {
        double points[GAUSS_LEGENDRE_MAX_ORDER];
        double weights[GAUSS_LEGENDRE_MAX_ORDER];
        double sum = 0.0;

        if (gauss_legendre_rule_1d(order, points, weights) != GAUSS_OK) {
            fprintf(stderr, "gauss 1D order %d failed\n", order);
            ++failures;
            continue;
        }
        for (int i = 0; i < order; ++i) {
            sum += weights[i];
        }
        failures += expect_close("1D weight sum", sum, 2.0);
    }

    return failures;
}

static int test_weight_sum_reference_cube_orders_1_to_8(void)
{
    int failures = 0;

    for (int order = 1; order <= 8; ++order) {
        GaussPoint3D points[GAUSS_LEGENDRE_MAX_CUBE_POINTS];
        int point_count = 0;
        double sum = 0.0;

        if (gauss_legendre_cube(order, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
                                points, GAUSS_LEGENDRE_MAX_CUBE_POINTS,
                                &point_count) != GAUSS_OK) {
            fprintf(stderr, "gauss cube order %d failed\n", order);
            ++failures;
            continue;
        }
        if (point_count != order * order * order) {
            fprintf(stderr, "gauss cube order %d expected %d points, got %d\n",
                    order, order * order * order, point_count);
            ++failures;
        }
        for (int i = 0; i < point_count; ++i) {
            sum += points[i].weight;
        }
        failures += expect_close("3D reference weight sum", sum, 8.0);
    }

    return failures;
}

static int test_physical_cell_constant_integrates_to_volume(void)
{
    const double xmin = 1.0;
    const double xmax = 3.0;
    const double ymin = -2.0;
    const double ymax = 1.0;
    const double zmin = 0.0;
    const double zmax = 4.0;
    const double volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);

    return expect_close("physical constant volume",
                        integrate_cell_order(3, xmin, xmax, ymin, ymax,
                                             zmin, zmax, function_one),
                        volume);
}

static int test_linear_integral_exact_order_one(void)
{
    const double xmin = 1.0;
    const double xmax = 3.0;
    const double ymin = -2.0;
    const double ymax = 1.0;
    const double zmin = 0.0;
    const double zmax = 4.0;
    const double volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    const double xc = 0.5 * (xmin + xmax);
    const double yc = 0.5 * (ymin + ymax);
    const double zc = 0.5 * (zmin + zmax);
    const double exact = volume * function_linear(xc, yc, zc);

    return expect_close("linear physical cell order 1",
                        integrate_cell_order(1, xmin, xmax, ymin, ymax,
                                             zmin, zmax, function_linear),
                        exact);
}

static int test_quadratic_improves_with_order_two(void)
{
    const double exact = 8.0 / 3.0;
    const double order1 = integrate_cell_order(1, -1.0, 1.0, -1.0, 1.0,
                                               -1.0, 1.0,
                                               function_x_squared);
    const double order2 = integrate_cell_order(2, -1.0, 1.0, -1.0, 1.0,
                                               -1.0, 1.0,
                                               function_x_squared);
    const double error1 = fabs(order1 - exact);
    const double error2 = fabs(order2 - exact);
    int failures = 0;

    if (!(error2 < error1)) {
        fprintf(stderr,
                "quadratic improvement: expected order2 error %.17g < order1 %.17g\n",
                error2, error1);
        ++failures;
    }
    failures += expect_close("quadratic exact order 2", order2, exact);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_weight_sum_1d_orders_1_to_8();
    failures += test_weight_sum_reference_cube_orders_1_to_8();
    failures += test_physical_cell_constant_integrates_to_volume();
    failures += test_linear_integral_exact_order_one();
    failures += test_quadratic_improves_with_order_two();
    failures += expect_close("integral 1", integrate_unit_cube(function_one), 1.0);
    failures += expect_close("integral x", integrate_unit_cube(function_x), 0.5);
    failures += expect_close("integral y", integrate_unit_cube(function_y), 0.5);
    failures += expect_close("integral z", integrate_unit_cube(function_z), 0.5);
    failures += expect_close("integral x*y", integrate_unit_cube(function_x_times_y), 0.25);

    return failures == 0 ? 0 : 1;
}
