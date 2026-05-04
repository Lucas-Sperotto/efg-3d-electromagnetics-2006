#include "weight.h"

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
 * Link-time smoke symbol for the weight-function module.
 */
int efg_weight_module_ready(void)
{
    return 1;
}
