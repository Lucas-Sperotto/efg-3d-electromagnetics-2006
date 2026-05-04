#include "weight.h"

#include <math.h>
#include <stddef.h>

static int sign_of_difference(double difference)
{
    if (difference > 0.0) {
        return 1;
    }

    if (difference < 0.0) {
        return -1;
    }

    return 0;
}

/*
 * Cubic spline weight function.
 *
 * Article equation: Eq. (10), cited in docs/03_resultados_numericos.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.2.
 * Mathematical meaning: compact-support scalar weight w(r) for MLS shape
 * functions in the EFG method.
 *
 * r is already normalized. Negative values are outside the physical meaning of
 * a normalized distance and are treated as outside the compact support.
 */
double weight_cubic_spline(double r)
{
    if (r < 0.0 || r > 1.0) {
        return 0.0;
    }

    if (r <= 0.5) {
        return (2.0 / 3.0) - (4.0 * r * r) + (4.0 * r * r * r);
    }

    return (4.0 / 3.0) - (4.0 * r) + (4.0 * r * r) - ((4.0 / 3.0) * r * r * r);
}

/*
 * Derivative of the cubic spline weight with respect to normalized r.
 *
 * This scalar derivative will be used later when computing grad Phi_I in the
 * MLS implementation. Tensor-product 3-D derivatives, shape-function
 * derivatives, stiffness assembly, Lagrange terms, and GMRES are intentionally
 * out of scope for this step.
 */
double weight_cubic_spline_derivative(double r)
{
    if (r < 0.0 || r > 1.0) {
        return 0.0;
    }

    if (r <= 0.5) {
        return (-8.0 * r) + (12.0 * r * r);
    }

    return -4.0 + (8.0 * r) - (4.0 * r * r);
}

/*
 * Quadratic weight function.
 *
 * Article equation: Eq. (11), cited in docs/03_resultados_numericos.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.2.
 * Mathematical meaning: comparison compact-support scalar weight w(r) for MLS
 * shape functions in the EFG method.
 *
 * r is already normalized. Negative values are outside the physical meaning of
 * a normalized distance and are treated as outside the compact support.
 */
double weight_quadratic(double r)
{
    if (r < 0.0 || r > 1.0) {
        return 0.0;
    }

    if (r <= 0.5) {
        return 1.0 - (2.0 * r * r);
    }

    return 2.0 * (1.0 - r) * (1.0 - r);
}

/*
 * Derivative of the quadratic weight with respect to normalized r.
 *
 * This scalar derivative will be used later when computing grad Phi_I in the
 * MLS implementation. Tensor-product 3-D derivatives, shape-function
 * derivatives, stiffness assembly, Lagrange terms, and GMRES are intentionally
 * out of scope for this step.
 */
double weight_quadratic_derivative(double r)
{
    if (r < 0.0 || r > 1.0) {
        return 0.0;
    }

    if (r <= 0.5) {
        return -4.0 * r;
    }

    return -4.0 * (1.0 - r);
}

int weight_cubic_tensor3d(double x,
                          double y,
                          double z,
                          double xI,
                          double yI,
                          double zI,
                          double d,
                          double *w)
{
    const double dx = fabs(x - xI);
    const double dy = fabs(y - yI);
    const double dz = fabs(z - zI);

    if (d <= 0.0 || w == NULL) {
        return WEIGHT_INVALID_ARGUMENT;
    }

    if (dx > d || dy > d || dz > d) {
        *w = 0.0;
        return WEIGHT_OK;
    }

    *w = weight_cubic_spline(dx / d) *
         weight_cubic_spline(dy / d) *
         weight_cubic_spline(dz / d);

    return WEIGHT_OK;
}

/*
 * Physical gradient of the cubic tensor-product weight.
 *
 * These helpers are for the future computation of grad Phi_I in the MLS linear
 * 3-D implementation. Tensorial derivatives of Phi_I, stiffness assembly,
 * Lagrange terms, and GMRES are intentionally not implemented here.
 */
int weight_cubic_tensor3d_gradient(double x,
                                   double y,
                                   double z,
                                   double xI,
                                   double yI,
                                   double zI,
                                   double d,
                                   double *dw_dx,
                                   double *dw_dy,
                                   double *dw_dz)
{
    const double x_difference = x - xI;
    const double y_difference = y - yI;
    const double z_difference = z - zI;
    const double dx = fabs(x_difference);
    const double dy = fabs(y_difference);
    const double dz = fabs(z_difference);
    double wx;
    double wy;
    double wz;
    double dw_drx;
    double dw_dry;
    double dw_drz;

    if (d <= 0.0 || dw_dx == NULL || dw_dy == NULL || dw_dz == NULL) {
        return WEIGHT_INVALID_ARGUMENT;
    }

    if (dx > d || dy > d || dz > d) {
        *dw_dx = 0.0;
        *dw_dy = 0.0;
        *dw_dz = 0.0;
        return WEIGHT_OK;
    }

    wx = weight_cubic_spline(dx / d);
    wy = weight_cubic_spline(dy / d);
    wz = weight_cubic_spline(dz / d);
    dw_drx = weight_cubic_spline_derivative(dx / d);
    dw_dry = weight_cubic_spline_derivative(dy / d);
    dw_drz = weight_cubic_spline_derivative(dz / d);

    *dw_dx = dw_drx * (double)sign_of_difference(x_difference) * wy * wz / d;
    *dw_dy = dw_dry * (double)sign_of_difference(y_difference) * wx * wz / d;
    *dw_dz = dw_drz * (double)sign_of_difference(z_difference) * wx * wy / d;

    return WEIGHT_OK;
}

/*
 * Link-time smoke symbol for the weight-function module.
 */
int efg_weight_module_ready(void)
{
    return 1;
}
