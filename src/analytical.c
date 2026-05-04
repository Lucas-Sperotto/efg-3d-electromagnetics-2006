#include "analytical.h"

#include <math.h>

static const double EFG_PI = 3.141592653589793238462643383279502884;

/*
 * Evaluate the truncated analytical potential in the cubic box.
 *
 * Article equations: Eqs. (13)-(15).
 * Formula reference: docs/03_resultados_numericos.md, section 3.1.
 * Mathematical meaning: Fourier-series potential for the cube benchmark with
 * V = 0 on the lateral/bottom grounded walls and V = V0 on the top wall.
 *
 * The article prints the particular case L = 10 and V0 = 10, where the
 * coefficient becomes 160. The expression below keeps the same equations in
 * parameterized form: 16 * V0 replaces 160 and L replaces 10.
 */
double analytical_potential_cube(double x,
                                 double y,
                                 double z,
                                 double L,
                                 double V0,
                                 int max_terms)
{
    double potential = 0.0;

    if (L <= 0.0 || max_terms <= 0) {
        return 0.0;
    }

    for (int n_index = 0; n_index < max_terms; ++n_index) {
        const int n = (2 * n_index) + 1;
        const double n_factor = ((double)n * EFG_PI) / L;
        const double sin_nx = sin(n_factor * x);

        for (int m_index = 0; m_index < max_terms; ++m_index) {
            const int m = (2 * m_index) + 1;
            const double m_factor = ((double)m * EFG_PI) / L;
            const double k_mn = sqrt((n_factor * n_factor) + (m_factor * m_factor));
            const double denominator =
                (double)m * (double)n * EFG_PI * EFG_PI * sinh(k_mn * L);
            const double a_mn = (16.0 * V0) / denominator;

            potential += a_mn * sin_nx * sin(m_factor * y) * sinh(k_mn * z);
        }
    }

    return potential;
}

/*
 * Link-time smoke symbol for the analytical benchmark.
 */
int efg_analytical_module_ready(void)
{
    return 1;
}
