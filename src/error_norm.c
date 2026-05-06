#include "error_norm.h"

#include "gauss.h"
#include "mls.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

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

static void domain_result_init(RelativeErrorDomainResult *result)
{
    if (result == NULL) {
        return;
    }

    result->error = NAN;
    result->numerator_integral = 0.0;
    result->denominator_integral = 0.0;
    result->quadrature_points = 0;
    result->mls_failures = 0;
}

/*
 * Compute the domain integral version of Eq. (16).
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: continuous-domain relative L2 error of the EFG
 * potential, approximated by summing the Gauss weights over all integration
 * cells used by the assembled problem.
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
                                        RelativeErrorDomainResult *result)
{
    return relative_error_norm_domain_integral_with_order(
        nodes, node_count, solution,
        xmin, xmax, ymin, ymax, zmin, zmax,
        nx_cells, ny_cells, nz_cells,
        2,
        exact, exact_context, result);
}

/*
 * Compute the domain integral version of Eq. (16) with explicit norm
 * quadrature order.
 *
 * Article equation: Eq. (16), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: continuous-domain relative L2 error of the EFG
 * potential. The EFG potential V_EFG(x_gp) is approximated at each quadrature
 * point by the MLS basis expansion sum_I Phi_I(x_gp) u_I.
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
    RelativeErrorDomainResult *result)
{
    MlsShapeValue *shape_values = NULL;
    double dx;
    double dy;
    double dz;

    if (nodes == NULL || solution == NULL || exact == NULL || result == NULL ||
        node_count <= 0 || xmin >= xmax || ymin >= ymax || zmin >= zmax ||
        nx_cells <= 0 || ny_cells <= 0 || nz_cells <= 0 ||
        quadrature_order < GAUSS_LEGENDRE_MIN_ORDER ||
        quadrature_order > GAUSS_LEGENDRE_MAX_ORDER) {
        return RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
    }

    domain_result_init(result);

    shape_values = malloc((size_t)node_count * sizeof(shape_values[0]));
    if (shape_values == NULL) {
        return RELATIVE_ERROR_NORM_ALLOCATION_FAILED;
    }

    dx = (xmax - xmin) / (double)nx_cells;
    dy = (ymax - ymin) / (double)ny_cells;
    dz = (zmax - zmin) / (double)nz_cells;

    for (int ix = 0; ix < nx_cells; ++ix) {
        const double cell_xmin = xmin + ((double)ix * dx);
        const double cell_xmax = cell_xmin + dx;

        for (int iy = 0; iy < ny_cells; ++iy) {
            const double cell_ymin = ymin + ((double)iy * dy);
            const double cell_ymax = cell_ymin + dy;

            for (int iz = 0; iz < nz_cells; ++iz) {
                GaussPoint3D gauss_points[GAUSS_LEGENDRE_MAX_CUBE_POINTS];
                const double cell_zmin = zmin + ((double)iz * dz);
                const double cell_zmax = cell_zmin + dz;
                int gauss_count = 0;

                if (gauss_legendre_cube(quadrature_order,
                                        cell_xmin, cell_xmax,
                                        cell_ymin, cell_ymax,
                                        cell_zmin, cell_zmax,
                                        gauss_points,
                                        GAUSS_LEGENDRE_MAX_CUBE_POINTS,
                                        &gauss_count) != GAUSS_OK) {
                    free(shape_values);
                    return RELATIVE_ERROR_NORM_INVALID_ARGUMENT;
                }

                for (int gp = 0; gp < gauss_count; ++gp) {
                    int value_count = 0;
                    double v_efg = 0.0;
                    const double x = gauss_points[gp].x;
                    const double y = gauss_points[gp].y;
                    const double z = gauss_points[gp].z;
                    const double weight = gauss_points[gp].weight;
                    const double v_exact = exact(x, y, z, exact_context);
                    const int mls_status =
                        mls_linear3d_shape_functions(nodes,
                                                     node_count,
                                                     x, y, z,
                                                     shape_values,
                                                     node_count,
                                                     &value_count);

                    ++result->quadrature_points;

                    if (mls_status != MLS_OK) {
                        ++result->mls_failures;
                        free(shape_values);
                        return RELATIVE_ERROR_NORM_MLS_FAILED;
                    }

                    for (int k = 0; k < value_count; ++k) {
                        v_efg += shape_values[k].phi
                                 * solution[shape_values[k].node_index];
                    }

                    {
                        const double diff = v_efg - v_exact;
                        result->numerator_integral += weight * diff * diff;
                        result->denominator_integral += weight * v_exact * v_exact;
                    }
                }
            }
        }
    }

    free(shape_values);

    if (result->denominator_integral == 0.0) {
        return RELATIVE_ERROR_NORM_ZERO_DENOMINATOR;
    }

    result->error = sqrt(result->numerator_integral /
                         result->denominator_integral);
    return RELATIVE_ERROR_NORM_OK;
}

/*
 * Link-time smoke symbol for error metrics.
 */
int efg_error_norm_module_ready(void)
{
    return 1;
}
