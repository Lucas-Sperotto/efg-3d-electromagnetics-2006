#ifndef GAUSS_H
#define GAUSS_H

/*
 * Gauss-Legendre quadrature helpers.
 *
 * The article uses Gauss quadrature of order 2 in cubic integration cells
 * (docs/03_resultados_numericos.md, sections 3.3 and 3.5). This module keeps
 * that order as the default while also exposing orders 1..8 for controlled
 * quadrature studies. It does not implement MLS, assembly, GMRES, or the EFG
 * solver.
 */

#define GAUSS_OK 0
#define GAUSS_INVALID_ARGUMENT 1

#define GAUSS_LEGENDRE_MIN_ORDER 1
#define GAUSS_LEGENDRE_MAX_ORDER 8
#define GAUSS_LEGENDRE_MAX_CUBE_POINTS 512

#define GAUSS_LEGENDRE_ORDER2_1D_POINTS 2
#define GAUSS_LEGENDRE_ORDER2_CUBE_POINTS 8

typedef struct GaussPoint1D {
    double xi;
    double weight;
} GaussPoint1D;

typedef struct GaussPoint3D {
    double x;
    double y;
    double z;
    double weight;
} GaussPoint3D;

/*
 * Fill the Gauss-Legendre points and weights on the reference interval
 * [-1, 1] for order 1..8.
 */
int gauss_legendre_rule_1d(int order, double *points, double *weights);

/*
 * Generate tensor-product Gauss-Legendre points for one rectangular cell.
 *
 * The returned coordinates are physical coordinates and the weights include
 * the mapping Jacobian from [-1,1]^3 to the requested cell.
 */
int gauss_legendre_cube(int order,
                        double x_min,
                        double x_max,
                        double y_min,
                        double y_max,
                        double z_min,
                        double z_max,
                        GaussPoint3D *points,
                        int max_points,
                        int *point_count);

/*
 * Fill the two Gauss-Legendre points and weights on the reference interval
 * [-1, 1].
 */
void gauss_legendre_order2_1d(GaussPoint1D points[GAUSS_LEGENDRE_ORDER2_1D_POINTS]);

/*
 * Generate the 8 quadrature points for one rectangular cubic cell.
 *
 * Each 3-D point is obtained by tensor product of the order-2 rule in x, y,
 * and z, then mapped from [-1, 1]^3 to the requested cell.
 */
void gauss_legendre_order2_cube(double x_min,
                                double x_max,
                                double y_min,
                                double y_max,
                                double z_min,
                                double z_max,
                                GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS]);

#endif
