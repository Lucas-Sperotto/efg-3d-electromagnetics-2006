#include "gauss.h"

#include <math.h>

void gauss_legendre_order2_1d(GaussPoint1D points[GAUSS_LEGENDRE_ORDER2_1D_POINTS])
{
    const double point = 1.0 / sqrt(3.0);

    /*
     * Order-2 Gauss-Legendre quadrature on [-1, 1].
     *
     * The two weights are 1.0. The article later uses this rule in cubic
     * integration cells (docs/03_resultados_numericos.md, sections 3.3 and
     * 3.5).
     */
    points[0].xi = -point;
    points[0].weight = 1.0;
    points[1].xi = point;
    points[1].weight = 1.0;
}

void gauss_legendre_order2_cube(double x_min,
                                double x_max,
                                double y_min,
                                double y_max,
                                double z_min,
                                double z_max,
                                GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS])
{
    GaussPoint1D one_dimensional[GAUSS_LEGENDRE_ORDER2_1D_POINTS];
    const double x_center = 0.5 * (x_min + x_max);
    const double y_center = 0.5 * (y_min + y_max);
    const double z_center = 0.5 * (z_min + z_max);
    const double x_half_length = 0.5 * (x_max - x_min);
    const double y_half_length = 0.5 * (y_max - y_min);
    const double z_half_length = 0.5 * (z_max - z_min);
    const double jacobian = x_half_length * y_half_length * z_half_length;
    int index = 0;

    gauss_legendre_order2_1d(one_dimensional);

    /*
     * Tensor product of the 1-D rule. With two points in each direction, one
     * cubic cell has 2 x 2 x 2 = 8 quadrature points.
     */
    for (int ix = 0; ix < GAUSS_LEGENDRE_ORDER2_1D_POINTS; ++ix) {
        for (int iy = 0; iy < GAUSS_LEGENDRE_ORDER2_1D_POINTS; ++iy) {
            for (int iz = 0; iz < GAUSS_LEGENDRE_ORDER2_1D_POINTS; ++iz) {
                points[index].x = x_center + (x_half_length * one_dimensional[ix].xi);
                points[index].y = y_center + (y_half_length * one_dimensional[iy].xi);
                points[index].z = z_center + (z_half_length * one_dimensional[iz].xi);
                points[index].weight =
                    one_dimensional[ix].weight *
                    one_dimensional[iy].weight *
                    one_dimensional[iz].weight *
                    jacobian;
                ++index;
            }
        }
    }
}

/*
 * Link-time smoke symbol for numerical integration.
 */
int efg_gauss_module_ready(void)
{
    return 1;
}
