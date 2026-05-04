#ifndef WEIGHT_H
#define WEIGHT_H

/*
 * Weight functions used by the Element-Free Galerkin method.
 *
 * The argument r is a normalized distance. This module implements only the
 * scalar functions w(r); it does not define the 3-D influence domain yet.
 */

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
 * Quadratic weight function.
 *
 * Article equation: Eq. (11), cited in docs/03_resultados_numericos.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.2.
 * Mathematical meaning: comparison compact-support scalar weight w(r) for MLS
 * shape functions in the EFG method.
 */
double weight_quadratic(double r);

#endif
