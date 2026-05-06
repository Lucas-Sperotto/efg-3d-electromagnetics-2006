#include "gauss.h"

#include <math.h>
#include <stddef.h>

static const double gauss_nodes[GAUSS_LEGENDRE_MAX_ORDER + 1][GAUSS_LEGENDRE_MAX_ORDER] = {
    {0.0},
    {0.0},
    {-0.57735026918962576451, 0.57735026918962576451},
    {-0.77459666924148337704, 0.0, 0.77459666924148337704},
    {-0.86113631159405257522, -0.33998104358485626480,
      0.33998104358485626480,  0.86113631159405257522},
    {-0.90617984593866399280, -0.53846931010568309104, 0.0,
      0.53846931010568309104,  0.90617984593866399280},
    {-0.93246951420315202781, -0.66120938646626451366,
     -0.23861918608319690863,  0.23861918608319690863,
      0.66120938646626451366,  0.93246951420315202781},
    {-0.94910791234275852453, -0.74153118559939443986,
     -0.40584515137739716691,  0.0,
      0.40584515137739716691,  0.74153118559939443986,
      0.94910791234275852453},
    {-0.96028985649753623168, -0.79666647741362673959,
     -0.52553240991632898582, -0.18343464249564980494,
      0.18343464249564980494,  0.52553240991632898582,
      0.79666647741362673959,  0.96028985649753623168}
};

static const double gauss_weights[GAUSS_LEGENDRE_MAX_ORDER + 1][GAUSS_LEGENDRE_MAX_ORDER] = {
    {0.0},
    {2.0},
    {1.0, 1.0},
    {0.55555555555555555556, 0.88888888888888888889,
     0.55555555555555555556},
    {0.34785484513745385737, 0.65214515486254614263,
     0.65214515486254614263, 0.34785484513745385737},
    {0.23692688505618908751, 0.47862867049936646804,
     0.56888888888888888889, 0.47862867049936646804,
     0.23692688505618908751},
    {0.17132449237917034504, 0.36076157304813860757,
     0.46791393457269104739, 0.46791393457269104739,
     0.36076157304813860757, 0.17132449237917034504},
    {0.12948496616886969327, 0.27970539148927666790,
     0.38183005050511894495, 0.41795918367346938776,
     0.38183005050511894495, 0.27970539148927666790,
     0.12948496616886969327},
    {0.10122853629037625915, 0.22238103445337447054,
     0.31370664587788728734, 0.36268378337836198297,
     0.36268378337836198297, 0.31370664587788728734,
     0.22238103445337447054, 0.10122853629037625915}
};

int gauss_legendre_rule_1d(int order, double *points, double *weights)
{
    if (points == NULL || weights == NULL ||
        order < GAUSS_LEGENDRE_MIN_ORDER ||
        order > GAUSS_LEGENDRE_MAX_ORDER) {
        return GAUSS_INVALID_ARGUMENT;
    }

    for (int i = 0; i < order; ++i) {
        points[i] = gauss_nodes[order][i];
        weights[i] = gauss_weights[order][i];
    }

    return GAUSS_OK;
}

int gauss_legendre_cube(int order,
                        double x_min,
                        double x_max,
                        double y_min,
                        double y_max,
                        double z_min,
                        double z_max,
                        GaussPoint3D *points,
                        int max_points,
                        int *point_count)
{
    double one_d_points[GAUSS_LEGENDRE_MAX_ORDER];
    double one_d_weights[GAUSS_LEGENDRE_MAX_ORDER];
    const double x_center = 0.5 * (x_min + x_max);
    const double y_center = 0.5 * (y_min + y_max);
    const double z_center = 0.5 * (z_min + z_max);
    const double x_half_length = 0.5 * (x_max - x_min);
    const double y_half_length = 0.5 * (y_max - y_min);
    const double z_half_length = 0.5 * (z_max - z_min);
    const double jacobian = x_half_length * y_half_length * z_half_length;
    const int required_points = order * order * order;
    int index = 0;

    if (points == NULL || point_count == NULL ||
        order < GAUSS_LEGENDRE_MIN_ORDER ||
        order > GAUSS_LEGENDRE_MAX_ORDER ||
        max_points < required_points ||
        x_min >= x_max || y_min >= y_max || z_min >= z_max) {
        return GAUSS_INVALID_ARGUMENT;
    }

    if (gauss_legendre_rule_1d(order, one_d_points, one_d_weights)
        != GAUSS_OK) {
        return GAUSS_INVALID_ARGUMENT;
    }

    /*
     * Tensor product of the 1-D rule. A rule of order n has n^3 quadrature
     * points in each physical integration cell.
     */
    for (int ix = 0; ix < order; ++ix) {
        for (int iy = 0; iy < order; ++iy) {
            for (int iz = 0; iz < order; ++iz) {
                points[index].x = x_center + (x_half_length * one_d_points[ix]);
                points[index].y = y_center + (y_half_length * one_d_points[iy]);
                points[index].z = z_center + (z_half_length * one_d_points[iz]);
                points[index].weight =
                    one_d_weights[ix] *
                    one_d_weights[iy] *
                    one_d_weights[iz] *
                    jacobian;
                ++index;
            }
        }
    }

    *point_count = index;
    return GAUSS_OK;
}

void gauss_legendre_order2_1d(GaussPoint1D points[GAUSS_LEGENDRE_ORDER2_1D_POINTS])
{
    double xi[GAUSS_LEGENDRE_ORDER2_1D_POINTS];
    double weight[GAUSS_LEGENDRE_ORDER2_1D_POINTS];

    /*
     * Order-2 Gauss-Legendre quadrature on [-1, 1].
     *
     * The two weights are 1.0. The article later uses this rule in cubic
     * integration cells (docs/03_resultados_numericos.md, sections 3.3 and
     * 3.5).
     */
    if (points == NULL ||
        gauss_legendre_rule_1d(2, xi, weight) != GAUSS_OK) {
        return;
    }

    for (int i = 0; i < GAUSS_LEGENDRE_ORDER2_1D_POINTS; ++i) {
        points[i].xi = xi[i];
        points[i].weight = weight[i];
    }
}

void gauss_legendre_order2_cube(double x_min,
                                double x_max,
                                double y_min,
                                double y_max,
                                double z_min,
                                double z_max,
                                GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS])
{
    int point_count = 0;

    (void)gauss_legendre_cube(2, x_min, x_max, y_min, y_max, z_min, z_max,
                              points, GAUSS_LEGENDRE_ORDER2_CUBE_POINTS,
                              &point_count);
}

/*
 * Link-time smoke symbol for numerical integration.
 */
int efg_gauss_module_ready(void)
{
    return 1;
}
