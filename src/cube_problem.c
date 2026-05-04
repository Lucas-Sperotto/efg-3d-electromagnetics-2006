#include "cube_problem.h"

#include <limits.h>
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
