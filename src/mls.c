#include "mls.h"

#include "mat4.h"
#include "weight.h"

#include <math.h>
#include <stddef.h>

static void linear_basis(double x, double y, double z, double p[4])
{
    p[0] = 1.0;
    p[1] = x;
    p[2] = y;
    p[3] = z;
}

static int node_contains_point_cubic_support(const Node3D *node, double x, double y, double z)
{
    return fabs(x - node->x) <= node->support_radius &&
           fabs(y - node->y) <= node->support_radius &&
           fabs(z - node->z) <= node->support_radius;
}

static double tensor_cubic_spline_weight(const Node3D *node, double x, double y, double z)
{
    const double rx = fabs(x - node->x) / node->support_radius;
    const double ry = fabs(y - node->y) / node->support_radius;
    const double rz = fabs(z - node->z) / node->support_radius;

    return weight_cubic_spline(rx) * weight_cubic_spline(ry) * weight_cubic_spline(rz);
}

static int map_mat4_status(int status)
{
    if (status == MAT4_SINGULAR) {
        return MLS_SINGULAR_MOMENT_MATRIX;
    }

    if (status == MAT4_NEAR_SINGULAR) {
        return MLS_NEAR_SINGULAR_MOMENT_MATRIX;
    }

    return MLS_INVALID_ARGUMENT;
}

/*
 * Compute linear 3-D MLS shape functions Phi_I(x), without derivatives.
 *
 * Mathematical reference: docs/05_mapa_equacoes_codigo.md maps the MLS
 * expression Phi_I(x) = p^T(x) A^{-1}(x) B_I(x) to this module.
 *
 * Implementation reference: docs/07_notas_sobre_esflib.md describes the
 * sequence used by ESFLib, but this code is an independent didactic
 * implementation and does not copy ESFLib code.
 */
int mls_linear3d_shape_functions(const Node3D *nodes,
                                 int node_count,
                                 double x,
                                 double y,
                                 double z,
                                 MlsShapeValue *values,
                                 int max_values,
                                 int *value_count)
{
    double moment[4][4];
    double p_eval[4];
    int support_count = 0;

    if (nodes == NULL || values == NULL || value_count == NULL ||
        node_count <= 0 || max_values <= 0) {
        return MLS_INVALID_ARGUMENT;
    }

    *value_count = 0;

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].support_radius <= 0.0) {
            return MLS_INVALID_SUPPORT_RADIUS;
        }
    }

    for (int i = 0; i < node_count; ++i) {
        if (node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            ++support_count;
        }
    }

    if (support_count < 4) {
        return MLS_NOT_ENOUGH_SUPPORT_NODES;
    }

    if (support_count > max_values) {
        return MLS_OUTPUT_CAPACITY_TOO_SMALL;
    }

    mat4_zero(moment);
    linear_basis(x, y, z, p_eval);

    for (int i = 0; i < node_count; ++i) {
        double p_node[4];
        const double weight = tensor_cubic_spline_weight(&nodes[i], x, y, z);

        if (!node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            continue;
        }

        linear_basis(nodes[i].x, nodes[i].y, nodes[i].z, p_node);

        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                moment[row][col] += weight * p_node[row] * p_node[col];
            }
        }
    }

    for (int i = 0; i < node_count; ++i) {
        double p_node[4];
        double weighted_basis[4];
        double solved[4];
        double phi = 0.0;
        int solve_status;
        const double weight = tensor_cubic_spline_weight(&nodes[i], x, y, z);

        if (!node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            continue;
        }

        linear_basis(nodes[i].x, nodes[i].y, nodes[i].z, p_node);

        for (int component = 0; component < 4; ++component) {
            weighted_basis[component] = weight * p_node[component];
        }

        solve_status = mat4_solve((const double(*)[4])moment, weighted_basis, solved);
        if (solve_status != MAT4_OK) {
            *value_count = 0;
            return map_mat4_status(solve_status);
        }

        for (int component = 0; component < 4; ++component) {
            phi += p_eval[component] * solved[component];
        }

        values[*value_count].node_index = i;
        values[*value_count].phi = phi;
        ++(*value_count);
    }

    return MLS_OK;
}

/*
 * Link-time smoke symbol for MLS shape functions.
 */
int efg_mls_module_ready(void)
{
    return 1;
}
