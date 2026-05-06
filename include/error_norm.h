#ifndef ERROR_NORM_H
#define ERROR_NORM_H

#include "nodes.h"

/*
 * Relative error norms for comparing numerical and exact potentials.
 *
 * The primary article metric is the domain-integrated Eq. (16). The discrete
 * finite-vector helper remains available for regression diagnostics and older
 * smoke tests.
 */

#define RELATIVE_ERROR_NORM_OK 0
#define RELATIVE_ERROR_NORM_INVALID_ARGUMENT 1
#define RELATIVE_ERROR_NORM_ZERO_DENOMINATOR 2
#define RELATIVE_ERROR_NORM_ALLOCATION_FAILED 3
#define RELATIVE_ERROR_NORM_MLS_FAILED 4

typedef double (*RelativeErrorExactPotential)(double x,
                                              double y,
                                              double z,
                                              void *context);

typedef struct RelativeErrorDomainResult {
    double error;
    double numerator_integral;
    double denominator_integral;
    int quadrature_points;
    int mls_failures;
} RelativeErrorDomainResult;

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

/*
 * Compute the domain-integrated relative error from Eq. (16):
 *
 * e = sqrt( int_Omega (V_EFG - V_exact)^2 dOmega
 *          / int_Omega V_exact^2 dOmega )
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: relative L2 error of the EFG potential field in the
 * whole domain, approximated with the same Gauss 2x2x2 cell quadrature used
 * in the stiffness assembly.
 */
int relative_error_norm_domain_integral(const Node3D *nodes,
                                        int node_count,
                                        const double *solution,
                                        double xmin,
                                        double xmax,
                                        double ymin,
                                        double ymax,
                                        double zmin,
                                        double zmax,
                                        int nx_cells,
                                        int ny_cells,
                                        int nz_cells,
                                        RelativeErrorExactPotential exact,
                                        void *exact_context,
                                        RelativeErrorDomainResult *result);

/*
 * Same Eq. (16) integral, with an explicit Gauss-Legendre order for the norm
 * evaluation. This separates the error-norm quadrature from the K assembly
 * quadrature in integration-order sensitivity studies.
 */
int relative_error_norm_domain_integral_with_order(
    const Node3D *nodes,
    int node_count,
    const double *solution,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zmin,
    double zmax,
    int nx_cells,
    int ny_cells,
    int nz_cells,
    int quadrature_order,
    RelativeErrorExactPotential exact,
    void *exact_context,
    RelativeErrorDomainResult *result);

#endif
