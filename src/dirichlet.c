#include "dirichlet.h"

#include "mls.h"

#include <stddef.h>
#include <stdlib.h>

static void zero_constraints(double *G, double *q, int point_count, int node_count)
{
    for (int i = 0; i < point_count * node_count; ++i) {
        G[i] = 0.0;
    }

    for (int i = 0; i < point_count; ++i) {
        q[i] = 0.0;
    }
}

/*
 * Assemble dense Dirichlet constraints for EFG Lagrange multipliers.
 *
 * Mathematical reference: docs/05_mapa_equacoes_codigo.md maps the boundary
 * block G and vector b/q to the Dirichlet constraint terms.
 *
 * Boundary-condition reference: docs/06_condicoes_de_contorno_e_lagrange.md
 * explains that EFG shape functions do not have the Kronecker-delta property,
 * so Dirichlet values are imposed through constraints rather than direct nodal
 * assignment.
 *
 * Each Dirichlet point x_b creates one row:
 *
 *     G[row, I] = Phi_I(x_b)
 *     q[row] = Vd(x_b)
 *
 * The augmented system, Lagrange multiplier unknowns, and GMRES solve are
 * intentionally outside this step.
 */
int dirichlet_assemble_constraints_dense(const Node3D *nodes,
                                         int node_count,
                                         const DirichletPoint *points,
                                         int point_count,
                                         double *G,
                                         double *q)
{
    MlsShapeValue *shape_values = NULL;

    if (nodes == NULL || points == NULL || G == NULL || q == NULL ||
        node_count <= 0 || point_count <= 0) {
        return DIRICHLET_INVALID_ARGUMENT;
    }

    zero_constraints(G, q, point_count, node_count);

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].support_radius <= 0.0) {
            return DIRICHLET_INVALID_SUPPORT_RADIUS;
        }
    }

    shape_values = malloc((size_t)node_count * sizeof(shape_values[0]));
    if (shape_values == NULL) {
        return DIRICHLET_ALLOCATION_FAILED;
    }

    for (int point_index = 0; point_index < point_count; ++point_index) {
        int value_count = 0;
        const DirichletPoint *point = &points[point_index];
        const int mls_status = mls_linear3d_shape_functions(nodes,
                                                           node_count,
                                                           point->x,
                                                           point->y,
                                                           point->z,
                                                           shape_values,
                                                           node_count,
                                                           &value_count);

        if (mls_status != MLS_OK) {
            free(shape_values);
            zero_constraints(G, q, point_count, node_count);
            return DIRICHLET_MLS_FAILED;
        }

        for (int value_index = 0; value_index < value_count; ++value_index) {
            const int node_index = shape_values[value_index].node_index;

            G[(point_index * node_count) + node_index] = shape_values[value_index].phi;
        }

        q[point_index] = point->value;
    }

    free(shape_values);
    return DIRICHLET_OK;
}
