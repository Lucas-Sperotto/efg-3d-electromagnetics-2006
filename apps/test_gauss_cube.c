#include "gauss.h"

#include <stdio.h>

/*
 * Small demonstration app for the order-2 Gauss rule in one unit cube.
 *
 * The article uses order-2 Gauss quadrature in cubic cells. This program only
 * prints the 8 points and checks simple polynomial integrals; it does not run
 * MLS, assembly, GMRES, or the EFG solver.
 */

static void print_integrals(const GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS])
{
    double integral_one = 0.0;
    double integral_x = 0.0;
    double integral_y = 0.0;
    double integral_z = 0.0;
    double integral_x_y = 0.0;

    for (int i = 0; i < GAUSS_LEGENDRE_ORDER2_CUBE_POINTS; ++i) {
        integral_one += points[i].weight;
        integral_x += points[i].weight * points[i].x;
        integral_y += points[i].weight * points[i].y;
        integral_z += points[i].weight * points[i].z;
        integral_x_y += points[i].weight * points[i].x * points[i].y;
    }

    printf("integral_1,%.17g\n", integral_one);
    printf("integral_x,%.17g\n", integral_x);
    printf("integral_y,%.17g\n", integral_y);
    printf("integral_z,%.17g\n", integral_z);
    printf("integral_x_y,%.17g\n", integral_x_y);
}

int main(void)
{
    GaussPoint3D points[GAUSS_LEGENDRE_ORDER2_CUBE_POINTS];

    gauss_legendre_order2_cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, points);

    printf("index,x,y,z,weight\n");
    for (int i = 0; i < GAUSS_LEGENDRE_ORDER2_CUBE_POINTS; ++i) {
        printf("%d,%.17g,%.17g,%.17g,%.17g\n",
               i,
               points[i].x,
               points[i].y,
               points[i].z,
               points[i].weight);
    }

    print_integrals(points);

    return 0;
}
