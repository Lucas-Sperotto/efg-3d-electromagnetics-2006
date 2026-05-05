#include "cube_problem.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>

static int valid_grid_dimensions(int nx, int ny, int nz)
{
    return nx >= 2 && ny >= 2 && nz >= 2;
}

static int product3_to_int(int a, int b, int c)
{
    const long long result = (long long)a * (long long)b * (long long)c;

    if (result > (long long)INT_MAX) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    return (int)result;
}

static double grid_coordinate(double L, int index, int count)
{
    return ((double)index * L) / (double)(count - 1);
}

static double cube_tolerance(double L)
{
    return 1.0e-12 + (fabs(L) * 1.0e-12);
}

static int same_coordinate(double a, double b, double tolerance)
{
    return fabs(a - b) <= tolerance;
}

static int same_node_position(const Node3D *node,
                              double x,
                              double y,
                              double z,
                              double tolerance)
{
    return same_coordinate(node->x, x, tolerance) &&
           same_coordinate(node->y, y, tolerance) &&
           same_coordinate(node->z, z, tolerance);
}

static int append_unique_node(Node3D *nodes,
                              int max_nodes,
                              int *node_count,
                              double x,
                              double y,
                              double z,
                              double L)
{
    const double tolerance = cube_tolerance(L);

    if (x < -tolerance || x > L + tolerance ||
        y < -tolerance || y > L + tolerance ||
        z < -tolerance || z > L + tolerance) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    if (same_coordinate(x, 0.0, tolerance)) {
        x = 0.0;
    } else if (same_coordinate(x, L, tolerance)) {
        x = L;
    }

    if (same_coordinate(y, 0.0, tolerance)) {
        y = 0.0;
    } else if (same_coordinate(y, L, tolerance)) {
        y = L;
    }

    if (same_coordinate(z, 0.0, tolerance)) {
        z = 0.0;
    } else if (same_coordinate(z, L, tolerance)) {
        z = L;
    }

    for (int i = 0; i < *node_count; ++i) {
        if (same_node_position(&nodes[i], x, y, z, tolerance)) {
            return CUBE_PROBLEM_OK;
        }
    }

    if (*node_count >= max_nodes) {
        return CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL;
    }

    nodes[*node_count] = (Node3D){x, y, z, 0.0};
    ++(*node_count);

    return CUBE_PROBLEM_OK;
}

static double top_refined_layer_coordinate(double L, int layer, int refine_n)
{
    /*
     * Extra layers occupy the upper quarter of the cube. This keeps the
     * construction deliberately simple while concentrating points where the
     * boundary changes from grounded lateral faces to the top electrode.
     */
    return L - ((0.25 * L * (double)layer) / (double)(refine_n - 1));
}

int cube_regular_node_count(int nx, int ny, int nz)
{
    if (!valid_grid_dimensions(nx, ny, nz)) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    return product3_to_int(nx, ny, nz);
}

/*
 * Generate a regular node cloud on [0,L]^3.
 *
 * The support radius is deliberately left as zero. The article-inspired
 * support parameter d_mI should be assigned later by the support module, after
 * the geometry exists.
 */
int cube_generate_regular_nodes(double L,
                                int nx,
                                int ny,
                                int nz,
                                Node3D *nodes,
                                int max_nodes,
                                int *node_count)
{
    int required_count;
    int output_index = 0;

    if (L <= 0.0 || nodes == NULL || node_count == NULL ||
        !valid_grid_dimensions(nx, ny, nz)) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    required_count = cube_regular_node_count(nx, ny, nz);
    if (required_count < 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    if (max_nodes < required_count) {
        return CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL;
    }

    for (int i = 0; i < nx; ++i) {
        const double x = grid_coordinate(L, i, nx);

        for (int j = 0; j < ny; ++j) {
            const double y = grid_coordinate(L, j, ny);

            for (int k = 0; k < nz; ++k) {
                const double z = grid_coordinate(L, k, nz);

                nodes[output_index] = (Node3D){x, y, z, 0.0};
                ++output_index;
            }
        }
    }

    *node_count = output_index;
    return CUBE_PROBLEM_OK;
}

/*
 * Generate a first didactic non-uniform cloud inspired by the article's
 * denser distribution near upper edges and corners.
 *
 * The function starts from the regular base_n^3 grid, then adds a refined
 * top slab and additional strips close to the upper lateral edges. Duplicates
 * are removed by coordinate comparison with a small tolerance. As with the
 * regular generator, support_radius is left at zero so the support module can
 * assign d_mI later.
 */
int cube_generate_top_refined_nodes(double L,
                                    int base_n,
                                    int refine_n,
                                    Node3D *nodes,
                                    int max_nodes,
                                    int *node_count)
{
    int status = CUBE_PROBLEM_OK;
    int output_count = 0;

    if (L <= 0.0 || nodes == NULL || node_count == NULL ||
        base_n < 2 || refine_n < 2 || max_nodes <= 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    for (int i = 0; i < base_n; ++i) {
        const double x = grid_coordinate(L, i, base_n);

        for (int j = 0; j < base_n; ++j) {
            const double y = grid_coordinate(L, j, base_n);

            for (int k = 0; k < base_n; ++k) {
                const double z = grid_coordinate(L, k, base_n);

                status = append_unique_node(nodes, max_nodes, &output_count, x, y, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }
            }
        }
    }

    for (int layer = 0; layer < refine_n; ++layer) {
        const double z = top_refined_layer_coordinate(L, layer, refine_n);

        for (int ix = 0; ix < refine_n; ++ix) {
            const double x = grid_coordinate(L, ix, refine_n);

            for (int iy = 0; iy < refine_n; ++iy) {
                const double y = grid_coordinate(L, iy, refine_n);

                status = append_unique_node(nodes, max_nodes, &output_count, x, y, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }
            }
        }
    }

    for (int layer = 0; layer < refine_n; ++layer) {
        const double z = top_refined_layer_coordinate(L, layer, refine_n);

        for (int offset_index = 0; offset_index < refine_n; ++offset_index) {
            const double edge_offset =
                (0.25 * L * (double)offset_index) / (double)(refine_n - 1);

            for (int t_index = 0; t_index < refine_n; ++t_index) {
                const double t = grid_coordinate(L, t_index, refine_n);
                const double lower_offset = edge_offset;
                const double upper_offset = L - edge_offset;

                status = append_unique_node(nodes, max_nodes, &output_count, t, lower_offset, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }

                status = append_unique_node(nodes, max_nodes, &output_count, t, upper_offset, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }

                status = append_unique_node(nodes, max_nodes, &output_count, lower_offset, t, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }

                status = append_unique_node(nodes, max_nodes, &output_count, upper_offset, t, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }
            }
        }
    }

    *node_count = output_count;
    return CUBE_PROBLEM_OK;
}

int cube_dirichlet_point_count(int nx, int ny, int nz)
{
    int total_count;
    int interior_count;

    if (!valid_grid_dimensions(nx, ny, nz)) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    total_count = product3_to_int(nx, ny, nz);
    interior_count = product3_to_int(nx - 2, ny - 2, nz - 2);

    if (total_count < 0 || interior_count < 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    return total_count - interior_count;
}

/*
 * Generate unique Dirichlet points on the cube boundary.
 *
 * Each surface grid point is emitted once, including edges and corners. The
 * top face z = L has priority, so upper corners receive V0 even though they
 * also belong to lateral faces.
 */
int cube_generate_dirichlet_points(double L,
                                   double V0,
                                   int nx,
                                   int ny,
                                   int nz,
                                   DirichletPoint *points,
                                   int max_points,
                                   int *point_count)
{
    int required_count;
    int output_index = 0;

    if (L <= 0.0 || points == NULL || point_count == NULL ||
        !valid_grid_dimensions(nx, ny, nz)) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    required_count = cube_dirichlet_point_count(nx, ny, nz);
    if (required_count < 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    if (max_points < required_count) {
        return CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL;
    }

    for (int i = 0; i < nx; ++i) {
        const double x = grid_coordinate(L, i, nx);

        for (int j = 0; j < ny; ++j) {
            const double y = grid_coordinate(L, j, ny);

            for (int k = 0; k < nz; ++k) {
                const int on_top_face = k == nz - 1;
                const int on_boundary = i == 0 || i == nx - 1 ||
                                        j == 0 || j == ny - 1 ||
                                        k == 0 || on_top_face;

                if (on_boundary) {
                    const double z = grid_coordinate(L, k, nz);
                    const double value = on_top_face ? V0 : 0.0;

                    points[output_index] = (DirichletPoint){x, y, z, value};
                    ++output_index;
                }
            }
        }
    }

    *point_count = output_index;
    return CUBE_PROBLEM_OK;
}
