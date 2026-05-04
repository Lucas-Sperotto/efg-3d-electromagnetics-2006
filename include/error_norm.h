#ifndef ERROR_NORM_H
#define ERROR_NORM_H

/*
 * Discrete relative error norm for comparing numerical and exact values.
 *
 * This module implements only the finite-vector version of the norm. It does
 * not perform continuous integration over the domain in this step.
 */

#define RELATIVE_ERROR_NORM_OK 0
#define RELATIVE_ERROR_NORM_INVALID_ARGUMENT 1
#define RELATIVE_ERROR_NORM_ZERO_DENOMINATOR 2

/*
 * Compute:
 *
 * e = sqrt(sum_i (numerical_i - exact_i)^2 / sum_i exact_i^2)
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: relative norm used to compare a numerical potential
 * field against the analytical reference solution.
 */
int relative_error_norm(const double *numerical, const double *exact, int n, double *error);

#endif
