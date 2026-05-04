#include "error_norm.h"

#include <math.h>
#include <stddef.h>

/*
 * Compute the discrete relative error norm.
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: relative norm used to compare a numerical potential
 * field against the analytical reference solution.
 *
 * This is the discrete version requested for the current step:
 *
 * e = sqrt(sum_i (numerical_i - exact_i)^2 / sum_i exact_i^2)
 */
int relative_error_norm(const double *numerical, const double *exact, int n, double *error)
{
    double numerator = 0.0;
    double denominator = 0.0;

    if (n <= 0 || numerical == NULL || exact == NULL || error == NULL) {
        return RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
    }

    for (int i = 0; i < n; ++i) {
        const double difference = numerical[i] - exact[i];

        numerator += difference * difference;
        denominator += exact[i] * exact[i];
    }

    if (denominator == 0.0) {
        return RELATIVE_ERROR_NORM_ZERO_DENOMINATOR;
    }

    *error = sqrt(numerator / denominator);
    return RELATIVE_ERROR_NORM_OK;
}

/*
 * Link-time smoke symbol for error metrics.
 */
int efg_error_norm_module_ready(void)
{
    return 1;
}
