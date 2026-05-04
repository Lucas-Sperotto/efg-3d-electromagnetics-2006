#ifndef ANALYTICAL_H
#define ANALYTICAL_H

/*
 * Analytical solution for the electrostatic cube benchmark.
 *
 * The function evaluates only the closed-form reference solution used to
 * compare numerical methods. It does not implement EFG, FEM, node generation,
 * MLS shape functions, or numerical assembly.
 */

/*
 * Evaluate the truncated analytical potential in the cubic box.
 *
 * Article equations: Eqs. (13)-(15).
 * Formula reference: docs/03_resultados_numericos.md, section 3.1.
 * Mathematical meaning: Fourier-series potential for the cube benchmark with
 * V = 0 on the lateral/bottom grounded walls and V = V0 on the top wall.
 *
 * max_terms is the number of odd terms used in each direction. For example,
 * max_terms = 3 uses indices 1, 3, and 5 for both m and n.
 */
double analytical_potential_cube(double x,
                                 double y,
                                 double z,
                                 double L,
                                 double V0,
                                 int max_terms);

#endif
