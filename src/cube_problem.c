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

static int same_dirichlet_position(const DirichletPoint *point,
                                   double x,
                                   double y,
                                   double z,
                                   double tolerance)
{
    return same_coordinate(point->x, x, tolerance) &&
           same_coordinate(point->y, y, tolerance) &&
           same_coordinate(point->z, z, tolerance);
}

static int coordinate_on_boundary(double coordinate, double L, double tolerance)
{
    return same_coordinate(coordinate, 0.0, tolerance) ||
           same_coordinate(coordinate, L, tolerance);
}

static int node_on_cube_surface(const Node3D *node, double L, double tolerance)
{
    return coordinate_on_boundary(node->x, L, tolerance) ||
           coordinate_on_boundary(node->y, L, tolerance) ||
           coordinate_on_boundary(node->z, L, tolerance);
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

static int append_unique_dirichlet_point(DirichletPoint *points,
                                         int max_points,
                                         int *point_count,
                                         double x,
                                         double y,
                                         double z,
                                         double value,
                                         double L)
{
    const double tolerance = cube_tolerance(L);

    if (x < -tolerance || x > L + tolerance ||
        y < -tolerance || y > L + tolerance ||
        z < -tolerance || z > L + tolerance) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    for (int i = 0; i < *point_count; ++i) {
        if (same_dirichlet_position(&points[i], x, y, z, tolerance)) {
            return CUBE_PROBLEM_OK;
        }
    }

    if (*point_count >= max_points) {
        return CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL;
    }

    points[*point_count] = (DirichletPoint){x, y, z, value};
    ++(*point_count);

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
 * The function starts from the regular base_n^3 grid, then adds a refined top
 * face and layered upper-edge lines. Duplicates are removed by coordinate
 * comparison with a small tolerance. As with the regular generator,
 * support_radius is left at zero so the support module can assign d_mI later.
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

    for (int layer = 0; layer <= 1; ++layer) {
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

    for (int layer = 2; layer < refine_n; ++layer) {
        const double z = top_refined_layer_coordinate(L, layer, refine_n);
        const double edge_offset =
            (0.25 * L * (double)layer) / (double)(refine_n - 1);
        const double lower_offset = edge_offset;
        const double upper_offset = L - edge_offset;

        for (int t_index = 0; t_index < refine_n; ++t_index) {
            const double t = grid_coordinate(L, t_index, refine_n);

            status = append_unique_node(nodes, max_nodes, &output_count, t, 0.0, z, L);
            if (status != CUBE_PROBLEM_OK) {
                return status;
            }

            status = append_unique_node(nodes, max_nodes, &output_count, t, L, z, L);
            if (status != CUBE_PROBLEM_OK) {
                return status;
            }

            status = append_unique_node(nodes, max_nodes, &output_count, 0.0, t, z, L);
            if (status != CUBE_PROBLEM_OK) {
                return status;
            }

            status = append_unique_node(nodes, max_nodes, &output_count, L, t, z, L);
            if (status != CUBE_PROBLEM_OK) {
                return status;
            }

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

/*
 * Generate Dirichlet points from an arbitrary cube node cloud.
 *
 * This variant is intended for non-uniform clouds, where a regular
 * nx-by-ny-by-nz surface count is not available. Every unique node on the cube
 * surface produces one Dirichlet point. The top face z = L has priority over
 * lateral and bottom faces, matching the electrostatic cube boundary
 * conditions.
 */
int cube_article_cloud_max_nodes(int base_n, int top_n, int n_extra_slices)
{
    if (base_n < 2 || top_n < 2 || n_extra_slices < 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    const long long base  = (long long)base_n * base_n * base_n;
    const long long face  = (long long)top_n * top_n;
    const long long inner = (long long)(top_n - 2) * (top_n - 2);
    const long long total = base + face + inner * n_extra_slices;

    if (total > (long long)INT_MAX) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    return (int)total;
}

/*
 * Generate a non-uniform node cloud inspired by Fig. 2 of the article.
 *
 * Builds a regular base_n^3 coarse grid, then overlays a full top_n x top_n
 * face at z = L (all Dirichlet V0) and n_extra_slices horizontal layers of
 * (top_n-2) x (top_n-2) interior-in-xy nodes at z-values uniformly spaced
 * within the top z_frac fraction of the cube. Only interior nodes are added
 * in the extra slices so that the Dirichlet set is not overwhelmed.
 */
int cube_generate_article_cloud(double L,
                                int base_n,
                                int top_n,
                                int n_extra_slices,
                                double z_frac,
                                Node3D *nodes,
                                int max_nodes,
                                int *node_count)
{
    int status = CUBE_PROBLEM_OK;
    int output_count = 0;

    if (L <= 0.0 || nodes == NULL || node_count == NULL ||
        base_n < 2 || top_n < 4 || n_extra_slices < 0 ||
        z_frac <= 0.0 || z_frac > 1.0 || max_nodes <= 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    /* Step 1: regular base_n^3 grid */
    for (int i = 0; i < base_n; ++i) {
        const double x = grid_coordinate(L, i, base_n);

        for (int j = 0; j < base_n; ++j) {
            const double y = grid_coordinate(L, j, base_n);

            for (int k = 0; k < base_n; ++k) {
                const double z = grid_coordinate(L, k, base_n);

                status = append_unique_node(nodes, max_nodes, &output_count,
                                            x, y, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }
            }
        }
    }

    /* Step 2: full top_n x top_n face at z = L (all become V0 Dirichlet) */
    for (int ix = 0; ix < top_n; ++ix) {
        const double x = grid_coordinate(L, ix, top_n);

        for (int iy = 0; iy < top_n; ++iy) {
            const double y = grid_coordinate(L, iy, top_n);

            status = append_unique_node(nodes, max_nodes, &output_count,
                                        x, y, L, L);
            if (status != CUBE_PROBLEM_OK) {
                return status;
            }
        }
    }

    /* Step 3: interior-in-xy fine slices in the top z_frac zone */
    for (int s = 1; s <= n_extra_slices; ++s) {
        const double z = L - (double)s * (L * z_frac) / (double)n_extra_slices;

        for (int ix = 1; ix <= top_n - 2; ++ix) {
            const double x = grid_coordinate(L, ix, top_n);

            for (int iy = 1; iy <= top_n - 2; ++iy) {
                const double y = grid_coordinate(L, iy, top_n);

                status = append_unique_node(nodes, max_nodes, &output_count,
                                            x, y, z, L);
                if (status != CUBE_PROBLEM_OK) {
                    return status;
                }
            }
        }
    }

    *node_count = output_count;
    return CUBE_PROBLEM_OK;
}

int cube_generate_dirichlet_points_from_nodes(double L,
                                              double V0,
                                              const Node3D *nodes,
                                              int node_count,
                                              DirichletPoint *points,
                                              int max_points,
                                              int *point_count)
{
    int output_count = 0;
    const double tolerance = cube_tolerance(L);

    if (L <= 0.0 || nodes == NULL || points == NULL || point_count == NULL ||
        node_count <= 0 || max_points <= 0) {
        return CUBE_PROBLEM_INVALID_ARGUMENT;
    }

    for (int i = 0; i < node_count; ++i) {
        int status;
        double value = 0.0;

        if (nodes[i].x < -tolerance || nodes[i].x > L + tolerance ||
            nodes[i].y < -tolerance || nodes[i].y > L + tolerance ||
            nodes[i].z < -tolerance || nodes[i].z > L + tolerance) {
            return CUBE_PROBLEM_INVALID_ARGUMENT;
        }

        if (!node_on_cube_surface(&nodes[i], L, tolerance)) {
            continue;
        }

        if (same_coordinate(nodes[i].z, L, tolerance)) {
            value = V0;
        }

        status = append_unique_dirichlet_point(points,
                                               max_points,
                                               &output_count,
                                               nodes[i].x,
                                               nodes[i].y,
                                               nodes[i].z,
                                               value,
                                               L);
        if (status != CUBE_PROBLEM_OK) {
            return status;
        }
    }

    *point_count = output_count;
    return CUBE_PROBLEM_OK;
}
