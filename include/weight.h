#ifndef WEIGHT_H
#define WEIGHT_H

/*
 * Weight functions used by the Element-Free Galerkin method.
 *
 * The argument r is a normalized distance. This module implements only the
 * scalar functions w(r); it does not define the 3-D influence domain yet.
 */

#define WEIGHT_OK 0
#define WEIGHT_INVALID_ARGUMENT 1

/*
 * Cubic spline weight function.
 *
 * Article equation: Eq. (10), cited in docs/03_resultados_numericos.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.2.
 * Mathematical meaning: compact-support scalar weight w(r) for MLS shape
 * functions in the EFG method.
 */
double weight_cubic_spline(double r);

/*
 * Derivative of the cubic spline weight with respect to normalized r.
 *
 * This scalar derivative will be used later when computing grad Phi_I in the
 * MLS implementation. It does not implement tensor-product 3-D derivatives or
 * shape-function derivatives yet.
 */
double weight_cubic_spline_derivative(double r);

/*
 * Quadratic weight function.
 *
 * Article equation: Eq. (11), cited in docs/03_resultados_numericos.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.2.
 * Mathematical meaning: comparison compact-support scalar weight w(r) for MLS
 * shape functions in the EFG method.
 */
double weight_quadratic(double r);

/*
 * Derivative of the quadratic weight with respect to normalized r.
 *
 * This scalar derivative will be used later when computing grad Phi_I in the
 * MLS implementation. It does not implement tensor-product 3-D derivatives or
 * shape-function derivatives yet.
 */
double weight_quadratic_derivative(double r);

/*
 * Cubic tensor-product weight for a cubic influence domain.
 *
 * This helper will be used later when computing grad Phi_I in the MLS linear
 * 3-D implementation. It does not implement shape-function derivatives,
 * stiffness assembly, Lagrange terms, or GMRES.
 */
int weight_cubic_tensor3d(double x,
                          double y,
                          double z,
                          double xI,
                          double yI,
                          double zI,
                          double d,
                          double *w);

/*
 * Physical gradient of the cubic tensor-product weight.
 *
 * The derivatives are with respect to physical x, y, and z. The scalar weight
 * derivative is still taken with respect to normalized r.
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
                                   double *dw_dz);

#endif
