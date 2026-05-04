#include "weight.h"

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

static double tensor_weight_value(double x, double y, double z)
{
    double value = -1.0;

    if (weight_cubic_tensor3d(x, y, z, 0.0, 0.0, 0.0, 1.0, &value) != WEIGHT_OK) {
        return -1.0;
    }

    return value;
}

static int test_weight_values(void)
{
    double value = -1.0;
    int failures = 0;
    const double center_expected = (2.0 / 3.0) * (2.0 / 3.0) * (2.0 / 3.0);
    const double internal_expected =
        weight_cubic_spline(0.25) * weight_cubic_spline(0.25) * weight_cubic_spline(0.25);

    failures += expect_code("center code",
                            weight_cubic_tensor3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, &value),
                            WEIGHT_OK);
    failures += expect_close("center weight", value, center_expected, 1.0e-12);

    failures += expect_code("outside x code",
                            weight_cubic_tensor3d(1.25, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, &value),
                            WEIGHT_OK);
    failures += expect_close("outside x weight", value, 0.0, 1.0e-12);

    failures += expect_code("outside y code",
                            weight_cubic_tensor3d(0.0, -1.25, 0.0, 0.0, 0.0, 0.0, 1.0, &value),
                            WEIGHT_OK);
    failures += expect_close("outside y weight", value, 0.0, 1.0e-12);

    failures += expect_code("outside z code",
                            weight_cubic_tensor3d(0.0, 0.0, 1.25, 0.0, 0.0, 0.0, 1.0, &value),
                            WEIGHT_OK);
    failures += expect_close("outside z weight", value, 0.0, 1.0e-12);

    failures += expect_code("internal code",
                            weight_cubic_tensor3d(0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 1.0, &value),
                            WEIGHT_OK);
    failures += expect_close("internal weight", value, internal_expected, 1.0e-12);

    return failures;
}

static int test_center_gradient_uses_zero_sign(void)
{
    double dw_dx = 1.0;
    double dw_dy = 1.0;
    double dw_dz = 1.0;
    int failures = 0;

    failures += expect_code("center gradient code",
                            weight_cubic_tensor3d_gradient(0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           1.0,
                                                           &dw_dx,
                                                           &dw_dy,
                                                           &dw_dz),
                            WEIGHT_OK);
    failures += expect_close("center dw_dx", dw_dx, 0.0, 1.0e-12);
    failures += expect_close("center dw_dy", dw_dy, 0.0, 1.0e-12);
    failures += expect_close("center dw_dz", dw_dz, 0.0, 1.0e-12);

    return failures;
}

static int test_gradient_finite_difference(void)
{
    const double h = 1.0e-6;
    const double x = 0.23;
    const double y = 0.37;
    const double z = 0.72;
    double dw_dx = 0.0;
    double dw_dy = 0.0;
    double dw_dz = 0.0;
    double fd_x;
    double fd_y;
    double fd_z;
    int failures = 0;

    failures += expect_code("gradient code",
                            weight_cubic_tensor3d_gradient(x,
                                                           y,
                                                           z,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           1.0,
                                                           &dw_dx,
                                                           &dw_dy,
                                                           &dw_dz),
                            WEIGHT_OK);

    fd_x = (tensor_weight_value(x + h, y, z) - tensor_weight_value(x - h, y, z)) / (2.0 * h);
    fd_y = (tensor_weight_value(x, y + h, z) - tensor_weight_value(x, y - h, z)) / (2.0 * h);
    fd_z = (tensor_weight_value(x, y, z + h) - tensor_weight_value(x, y, z - h)) / (2.0 * h);

    failures += expect_close("finite difference dw_dx", dw_dx, fd_x, 1.0e-8);
    failures += expect_close("finite difference dw_dy", dw_dy, fd_y, 1.0e-8);
    failures += expect_close("finite difference dw_dz", dw_dz, fd_z, 1.0e-8);

    return failures;
}

static int test_invalid_arguments(void)
{
    double value = 0.0;
    double dw_dx = 0.0;
    double dw_dy = 0.0;
    double dw_dz = 0.0;
    int failures = 0;

    failures += expect_code("weight d zero",
                            weight_cubic_tensor3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &value),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("weight d negative",
                            weight_cubic_tensor3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, &value),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("weight null output",
                            weight_cubic_tensor3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, NULL),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("gradient d zero",
                            weight_cubic_tensor3d_gradient(0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           &dw_dx,
                                                           &dw_dy,
                                                           &dw_dz),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("gradient null dx",
                            weight_cubic_tensor3d_gradient(0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           1.0,
                                                           NULL,
                                                           &dw_dy,
                                                           &dw_dz),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("gradient null dy",
                            weight_cubic_tensor3d_gradient(0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           1.0,
                                                           &dw_dx,
                                                           NULL,
                                                           &dw_dz),
                            WEIGHT_INVALID_ARGUMENT);
    failures += expect_code("gradient null dz",
                            weight_cubic_tensor3d_gradient(0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           0.0,
                                                           1.0,
                                                           &dw_dx,
                                                           &dw_dy,
                                                           NULL),
                            WEIGHT_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_weight_values();
    failures += test_center_gradient_uses_zero_sign();
    failures += test_gradient_finite_difference();
    failures += test_invalid_arguments();

    return failures == 0 ? 0 : 1;
}
