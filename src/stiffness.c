#include "stiffness.h"

#include "gauss.h"
#include "mls.h"

#include <stddef.h>
#include <stdlib.h>

/*
 * Assemble local stiffness entries for one axis-aligned cubic cell.
 *
 * Mathematical reference: docs/05_mapa_equacoes_codigo.md maps the stiffness
 * matrix K_ij to the weak-form term integral grad(Phi_i) dot grad(Phi_j).
 *
 * Boundary-condition reference: docs/06_condicoes_de_contorno_e_lagrange.md
 * explains that Dirichlet/Lagrange treatment is separate. This function does
 * not apply boundary conditions, Lagrange multipliers, or GMRES.
 *
 * Convention note: gauss_legendre_cube() returns physical coordinates and
 * physical weights. Its weights already include the cell Jacobian
 * (xmax-xmin)(ymax-ymin)(zmax-zmin)/8, so this function does not multiply by a
 * second Jacobian.
 *
 * Entries are intentionally not consolidated. Repeated (row_node, col_node)
 * pairs can appear once per active pair per Gauss point; callers may consolidate
 * them later when a global matrix representation exists.
 */
int stiffness_assemble_cell(const Node3D *nodes,
                            int node_count,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            double zmin,
                            double zmax,
                            StiffnessEntry *entries,
                            int max_entries,
                            int *entry_count)
{
    return stiffness_assemble_cell_with_order(nodes, node_count,
                                              xmin, xmax, ymin, ymax,
                                              zmin, zmax,
                                              2,
                                              entries, max_entries,
                                              entry_count);
}

/*
 * Assemble local stiffness entries for Eq. (8).
 *
 * Article equation: Eq. (8), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: local contribution to K_ij from
 * grad(Phi_i) dot grad(Phi_j), evaluated by Gauss-Legendre quadrature of the
 * requested order in the physical integration cell.
 */
int stiffness_assemble_cell_with_order(const Node3D *nodes,
                                       int node_count,
                                       double xmin,
                                       double xmax,
                                       double ymin,
                                       double ymax,
                                       double zmin,
                                       double zmax,
                                       int quadrature_order,
                                       StiffnessEntry *entries,
                                       int max_entries,
                                       int *entry_count)
{
    GaussPoint3D gauss_points[GAUSS_LEGENDRE_MAX_CUBE_POINTS];
    MlsShapeGradientValue *shape_values = NULL;
    int gauss_count = 0;

    if (nodes == NULL || entries == NULL || entry_count == NULL ||
        node_count <= 0 || max_entries <= 0 ||
        quadrature_order < GAUSS_LEGENDRE_MIN_ORDER ||
        quadrature_order > GAUSS_LEGENDRE_MAX_ORDER) {
        return STIFFNESS_INVALID_ARGUMENT;
    }

    *entry_count = 0;

    if (xmin >= xmax || ymin >= ymax || zmin >= zmax) {
        return STIFFNESS_INVALID_CELL;
    }

    shape_values = malloc((size_t)node_count * sizeof(shape_values[0]));
    if (shape_values == NULL) {
        return STIFFNESS_INVALID_ARGUMENT;
    }

    if (gauss_legendre_cube(quadrature_order, xmin, xmax, ymin, ymax,
                            zmin, zmax, gauss_points,
                            GAUSS_LEGENDRE_MAX_CUBE_POINTS,
                            &gauss_count) != GAUSS_OK) {
        free(shape_values);
        return STIFFNESS_INVALID_CELL;
    }

    for (int gp = 0; gp < gauss_count; ++gp) {
        int value_count = 0;
        const int mls_status =
            mls_linear3d_shape_functions_with_gradients(nodes,
                                                        node_count,
                                                        gauss_points[gp].x,
                                                        gauss_points[gp].y,
                                                        gauss_points[gp].z,
                                                        shape_values,
                                                        node_count,
                                                        &value_count);

        if (mls_status != MLS_OK) {
            free(shape_values);
            *entry_count = 0;
            return STIFFNESS_MLS_FAILED;
        }

        if (*entry_count + (value_count * value_count) > max_entries) {
            free(shape_values);
            return STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL;
        }

        for (int i = 0; i < value_count; ++i) {
            for (int j = 0; j < value_count; ++j) {
                const double gradient_dot =
                    (shape_values[i].dphi_dx * shape_values[j].dphi_dx) +
                    (shape_values[i].dphi_dy * shape_values[j].dphi_dy) +
                    (shape_values[i].dphi_dz * shape_values[j].dphi_dz);

                entries[*entry_count].row_node = shape_values[i].node_index;
                entries[*entry_count].col_node = shape_values[j].node_index;
                entries[*entry_count].value = gauss_points[gp].weight * gradient_dot;
                ++(*entry_count);
            }
        }
    }

    free(shape_values);
    return STIFFNESS_OK;
}
