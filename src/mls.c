#include "mls.h"

#include "mat4.h"
#include "weight.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

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
    double weight = 0.0;

    if (weight_cubic_tensor3d(x,
                              y,
                              z,
                              node->x,
                              node->y,
                              node->z,
                              node->support_radius,
                              &weight) != WEIGHT_OK) {
        return 0.0;
    }

    return weight;
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

static void add_weighted_outer_product(double matrix[4][4], double scale, const double p_node[4])
{
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            matrix[row][col] += scale * p_node[row] * p_node[col];
        }
    }
}

static void mat4_vector_product(const double A[4][4], const double x[4], double result[4])
{
    for (int row = 0; row < 4; ++row) {
        result[row] = 0.0;
        for (int col = 0; col < 4; ++col) {
            result[row] += A[row][col] * x[col];
        }
    }
}

static double dot4(const double a[4], const double b[4])
{
    double result = 0.0;

    for (int i = 0; i < 4; ++i) {
        result += a[i] * b[i];
    }

    return result;
}

static int validate_mls_inputs(const Node3D *nodes,
                               int node_count,
                               const void *values,
                               int max_values,
                               int *value_count)
{
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

    return MLS_OK;
}

static int count_support_nodes(const Node3D *nodes,
                               int node_count,
                               double x,
                               double y,
                               double z)
{
    int support_count = 0;

    for (int i = 0; i < node_count; ++i) {
        if (node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            ++support_count;
        }
    }

    return support_count;
}

static int validate_support_capacity(int support_count, int max_values)
{
    if (support_count < 4) {
        return MLS_NOT_ENOUGH_SUPPORT_NODES;
    }

    if (support_count > max_values) {
        return MLS_OUTPUT_CAPACITY_TOO_SMALL;
    }

    return MLS_OK;
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
    MlsShapeGradientValue *gradient_values = NULL;
    int status;
    int gradient_count = 0;

    if (max_values <= 0) {
        return MLS_INVALID_ARGUMENT;
    }

    gradient_values = malloc((size_t)max_values * sizeof(gradient_values[0]));
    if (gradient_values == NULL) {
        return MLS_INVALID_ARGUMENT;
    }

    status = mls_linear3d_shape_functions_with_gradients(nodes,
                                                         node_count,
                                                         x,
                                                         y,
                                                         z,
                                                         gradient_values,
                                                         max_values,
                                                         &gradient_count);
    if (status != MLS_OK) {
        if (value_count != NULL) {
            *value_count = 0;
        }
        free(gradient_values);
        return status;
    }

    for (int i = 0; i < gradient_count; ++i) {
        values[i].node_index = gradient_values[i].node_index;
        values[i].phi = gradient_values[i].phi;
    }

    *value_count = gradient_count;
    free(gradient_values);
    return MLS_OK;
}

/*
 * Compute linear 3-D MLS shape functions and gradients.
 *
 * Mathematical reference: docs/05_mapa_equacoes_codigo.md maps MLS
 * construction to this module. The gradient follows the differentiated
 * q_I = A^{-1} B_I relation used for MLS shape functions.
 *
 * Implementation reference: docs/07_notas_sobre_esflib.md is used only as a
 * conceptual checklist. This code is independent and does not copy ESFLib.
 */
int mls_linear3d_shape_functions_with_gradients(const Node3D *nodes,
                                                int node_count,
                                                double x,
                                                double y,
                                                double z,
                                                MlsShapeGradientValue *values,
                                                int max_values,
                                                int *value_count)
{
    double moment[4][4];
    double dmoment_dx[4][4];
    double dmoment_dy[4][4];
    double dmoment_dz[4][4];
    double p_eval[4];
    int support_count;
    int status;

    status = validate_mls_inputs(nodes, node_count, values, max_values, value_count);
    if (status != MLS_OK) {
        return status;
    }

    support_count = count_support_nodes(nodes, node_count, x, y, z);
    status = validate_support_capacity(support_count, max_values);
    if (status != MLS_OK) {
        return status;
    }

    mat4_zero(moment);
    mat4_zero(dmoment_dx);
    mat4_zero(dmoment_dy);
    mat4_zero(dmoment_dz);
    linear_basis(x, y, z, p_eval);

    for (int i = 0; i < node_count; ++i) {
        double p_node[4];
        double weight;
        double dw_dx;
        double dw_dy;
        double dw_dz;

        if (!node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            continue;
        }

        linear_basis(nodes[i].x, nodes[i].y, nodes[i].z, p_node);
        weight = tensor_cubic_spline_weight(&nodes[i], x, y, z);
        weight_cubic_tensor3d_gradient(x,
                                       y,
                                       z,
                                       nodes[i].x,
                                       nodes[i].y,
                                       nodes[i].z,
                                       nodes[i].support_radius,
                                       &dw_dx,
                                       &dw_dy,
                                       &dw_dz);

        add_weighted_outer_product(moment, weight, p_node);
        add_weighted_outer_product(dmoment_dx, dw_dx, p_node);
        add_weighted_outer_product(dmoment_dy, dw_dy, p_node);
        add_weighted_outer_product(dmoment_dz, dw_dz, p_node);
    }

    for (int i = 0; i < node_count; ++i) {
        double p_node[4];
        double weighted_basis[4];
        double dweighted_basis_dx[4];
        double dweighted_basis_dy[4];
        double dweighted_basis_dz[4];
        double q[4];
        double dmoment_q_dx[4];
        double dmoment_q_dy[4];
        double dmoment_q_dz[4];
        double rhs_dx[4];
        double rhs_dy[4];
        double rhs_dz[4];
        double r_dx[4];
        double r_dy[4];
        double r_dz[4];
        int solve_status;
        double weight;
        double dw_dx;
        double dw_dy;
        double dw_dz;

        if (!node_contains_point_cubic_support(&nodes[i], x, y, z)) {
            continue;
        }

        linear_basis(nodes[i].x, nodes[i].y, nodes[i].z, p_node);
        weight = tensor_cubic_spline_weight(&nodes[i], x, y, z);
        weight_cubic_tensor3d_gradient(x,
                                       y,
                                       z,
                                       nodes[i].x,
                                       nodes[i].y,
                                       nodes[i].z,
                                       nodes[i].support_radius,
                                       &dw_dx,
                                       &dw_dy,
                                       &dw_dz);

        for (int component = 0; component < 4; ++component) {
            weighted_basis[component] = weight * p_node[component];
            dweighted_basis_dx[component] = dw_dx * p_node[component];
            dweighted_basis_dy[component] = dw_dy * p_node[component];
            dweighted_basis_dz[component] = dw_dz * p_node[component];
        }

        solve_status = mat4_solve((const double(*)[4])moment, weighted_basis, q);
        if (solve_status != MAT4_OK) {
            *value_count = 0;
            return map_mat4_status(solve_status);
        }

        mat4_vector_product((const double(*)[4])dmoment_dx, q, dmoment_q_dx);
        mat4_vector_product((const double(*)[4])dmoment_dy, q, dmoment_q_dy);
        mat4_vector_product((const double(*)[4])dmoment_dz, q, dmoment_q_dz);

        for (int component = 0; component < 4; ++component) {
            rhs_dx[component] = dweighted_basis_dx[component] - dmoment_q_dx[component];
            rhs_dy[component] = dweighted_basis_dy[component] - dmoment_q_dy[component];
            rhs_dz[component] = dweighted_basis_dz[component] - dmoment_q_dz[component];
        }

        solve_status = mat4_solve((const double(*)[4])moment, rhs_dx, r_dx);
        if (solve_status != MAT4_OK) {
            *value_count = 0;
            return map_mat4_status(solve_status);
        }

        solve_status = mat4_solve((const double(*)[4])moment, rhs_dy, r_dy);
        if (solve_status != MAT4_OK) {
            *value_count = 0;
            return map_mat4_status(solve_status);
        }

        solve_status = mat4_solve((const double(*)[4])moment, rhs_dz, r_dz);
        if (solve_status != MAT4_OK) {
            *value_count = 0;
            return map_mat4_status(solve_status);
        }

        values[*value_count].node_index = i;
        values[*value_count].phi = dot4(p_eval, q);
        values[*value_count].dphi_dx = q[1] + dot4(p_eval, r_dx);
        values[*value_count].dphi_dy = q[2] + dot4(p_eval, r_dy);
        values[*value_count].dphi_dz = q[3] + dot4(p_eval, r_dz);
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
