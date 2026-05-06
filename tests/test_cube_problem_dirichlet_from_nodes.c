#include "cube_problem.h"

#include <math.h>
#include <stdio.h>

#define TEST_L 10.0
#define TEST_V0 10.0
#define TEST_TOLERANCE 1.0e-10

static int expect_code(const char *name, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s: expected code %d, got %d\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int same_point(const DirichletPoint *a, const DirichletPoint *b)
{
    return fabs(a->x - b->x) <= TEST_TOLERANCE &&
           fabs(a->y - b->y) <= TEST_TOLERANCE &&
           fabs(a->z - b->z) <= TEST_TOLERANCE;
}

static int same_coordinates(double ax,
                            double ay,
                            double az,
                            double bx,
                            double by,
                            double bz)
{
    return fabs(ax - bx) <= TEST_TOLERANCE &&
           fabs(ay - by) <= TEST_TOLERANCE &&
           fabs(az - bz) <= TEST_TOLERANCE;
}

static int point_on_surface(const DirichletPoint *point)
{
    return fabs(point->x) <= TEST_TOLERANCE ||
           fabs(point->x - TEST_L) <= TEST_TOLERANCE ||
           fabs(point->y) <= TEST_TOLERANCE ||
           fabs(point->y - TEST_L) <= TEST_TOLERANCE ||
           fabs(point->z) <= TEST_TOLERANCE ||
           fabs(point->z - TEST_L) <= TEST_TOLERANCE;
}

static int point_inside_cube(const DirichletPoint *point)
{
    return point->x >= -TEST_TOLERANCE && point->x <= TEST_L + TEST_TOLERANCE &&
           point->y >= -TEST_TOLERANCE && point->y <= TEST_L + TEST_TOLERANCE &&
           point->z >= -TEST_TOLERANCE && point->z <= TEST_L + TEST_TOLERANCE;
}

static int validate_dirichlet_values(const DirichletPoint *points, int point_count)
{
    int failures = 0;

    for (int i = 0; i < point_count; ++i) {
        if (!point_inside_cube(&points[i])) {
            fprintf(stderr, "Dirichlet point %d is outside [0,L]^3\n", i);
            ++failures;
        }

        if (!point_on_surface(&points[i])) {
            fprintf(stderr, "Dirichlet point %d is not on the cube surface\n", i);
            ++failures;
        }

        const int on_open_top = fabs(points[i].z - TEST_L) <= TEST_TOLERANCE &&
                                fabs(points[i].x) > TEST_TOLERANCE &&
                                fabs(points[i].x - TEST_L) > TEST_TOLERANCE &&
                                fabs(points[i].y) > TEST_TOLERANCE &&
                                fabs(points[i].y - TEST_L) > TEST_TOLERANCE;

        if (on_open_top) {
            if (fabs(points[i].value - TEST_V0) > TEST_TOLERANCE) {
                fprintf(stderr, "open top point %d should have value V0\n", i);
                ++failures;
            }
        } else if (fabs(points[i].value) > TEST_TOLERANCE) {
            fprintf(stderr, "grounded boundary point %d should have value zero\n", i);
            ++failures;
        }
    }

    return failures;
}

static int validate_no_duplicates(const DirichletPoint *points, int point_count)
{
    int failures = 0;

    for (int i = 0; i < point_count; ++i) {
        for (int j = i + 1; j < point_count; ++j) {
            if (same_point(&points[i], &points[j])) {
                fprintf(stderr, "duplicate Dirichlet points at indices %d and %d\n", i, j);
                ++failures;
            }
        }
    }

    return failures;
}

static int find_matching_point(const DirichletPoint *points,
                               int point_count,
                               const DirichletPoint *target)
{
    for (int i = 0; i < point_count; ++i) {
        if (same_point(&points[i], target) &&
            fabs(points[i].value - target->value) <= TEST_TOLERANCE) {
            return 1;
        }
    }

    return 0;
}

static int test_regular_grid_matches_existing_generator(void)
{
    Node3D nodes[27];
    DirichletPoint from_grid[26];
    DirichletPoint from_nodes[26];
    int node_count = 0;
    int grid_point_count = 0;
    int node_point_count = 0;
    int failures = 0;

    failures += expect_code("generate regular nodes",
                            cube_generate_regular_nodes(TEST_L, 3, 3, 3, nodes, 27, &node_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("generate grid Dirichlet",
                            cube_generate_dirichlet_points(TEST_L,
                                                           TEST_V0,
                                                           3,
                                                           3,
                                                           3,
                                                           from_grid,
                                                           26,
                                                           &grid_point_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("generate Dirichlet from regular nodes",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      node_count,
                                                                      from_nodes,
                                                                      26,
                                                                      &node_point_count),
                            CUBE_PROBLEM_OK);

    if (failures != 0) {
        return failures;
    }

    if (grid_point_count != 26 || node_point_count != 26) {
        fprintf(stderr,
                "expected both regular Dirichlet generators to produce 26 points, got %d and %d\n",
                grid_point_count,
                node_point_count);
        ++failures;
    }

    failures += validate_dirichlet_values(from_nodes, node_point_count);
    failures += validate_no_duplicates(from_nodes, node_point_count);

    for (int i = 0; i < grid_point_count; ++i) {
        if (!find_matching_point(from_nodes, node_point_count, &from_grid[i])) {
            fprintf(stderr,
                    "missing regular Dirichlet point (%.17g, %.17g, %.17g)\n",
                    from_grid[i].x,
                    from_grid[i].y,
                    from_grid[i].z);
            ++failures;
        }
    }

    return failures;
}

static int test_top_refined_nodes_generate_surface_dirichlet(void)
{
    Node3D nodes[1024];
    DirichletPoint points[1024];
    int node_count = 0;
    int point_count = 0;
    int top_corner_count = 0;
    int failures = 0;

    failures += expect_code("generate top-refined nodes",
                            cube_generate_top_refined_nodes(TEST_L, 3, 5, nodes, 1024, &node_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("generate Dirichlet from top-refined nodes",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      node_count,
                                                                      points,
                                                                      1024,
                                                                      &point_count),
                            CUBE_PROBLEM_OK);

    if (failures != 0) {
        return failures;
    }

    if (point_count <= 0) {
        fprintf(stderr, "expected at least one Dirichlet point from non-uniform cloud\n");
        ++failures;
    }

    failures += validate_dirichlet_values(points, point_count);
    failures += validate_no_duplicates(points, point_count);

    for (int i = 0; i < point_count; ++i) {
        if (same_coordinates(points[i].x, points[i].y, points[i].z, 0.0, 0.0, TEST_L) ||
            same_coordinates(points[i].x, points[i].y, points[i].z, TEST_L, 0.0, TEST_L) ||
            same_coordinates(points[i].x, points[i].y, points[i].z, 0.0, TEST_L, TEST_L) ||
            same_coordinates(points[i].x, points[i].y, points[i].z, TEST_L, TEST_L, TEST_L)) {
            if (fabs(points[i].value) > TEST_TOLERANCE) {
                fprintf(stderr, "upper corner %d should remain grounded\n", i);
                ++failures;
            }
            ++top_corner_count;
        }
    }

    if (top_corner_count != 4) {
        fprintf(stderr, "expected 4 upper Dirichlet corners, got %d\n", top_corner_count);
        ++failures;
    }

    return failures;
}

static int test_duplicate_nodes_are_merged(void)
{
    const Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {TEST_L, TEST_L, TEST_L, 0.0},
        {TEST_L, TEST_L, TEST_L, 0.0},
        {0.5 * TEST_L, 0.5 * TEST_L, 0.5 * TEST_L, 0.0}
    };
    DirichletPoint points[4];
    int point_count = 0;
    int failures = 0;

    failures += expect_code("generate from duplicate nodes",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      5,
                                                                      points,
                                                                      4,
                                                                      &point_count),
                            CUBE_PROBLEM_OK);

    if (point_count != 2) {
        fprintf(stderr, "expected duplicate surface nodes to merge into 2 points, got %d\n",
                point_count);
        ++failures;
    }

    failures += validate_no_duplicates(points, point_count);
    failures += validate_dirichlet_values(points, point_count);

    return failures;
}

static int test_error_paths(void)
{
    const Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {TEST_L, TEST_L, TEST_L, 0.0}
    };
    const Node3D invalid_nodes[] = {
        {-1.0, 0.0, 0.0, 0.0}
    };
    DirichletPoint points[2];
    int point_count = 0;
    int failures = 0;

    failures += expect_code("invalid L",
                            cube_generate_dirichlet_points_from_nodes(0.0,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      2,
                                                                      points,
                                                                      2,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null nodes",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      NULL,
                                                                      2,
                                                                      points,
                                                                      2,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null points",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      2,
                                                                      NULL,
                                                                      2,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null point_count",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      2,
                                                                      points,
                                                                      2,
                                                                      NULL),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid node_count",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      0,
                                                                      points,
                                                                      2,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid max_points",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      2,
                                                                      points,
                                                                      0,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("insufficient capacity",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      nodes,
                                                                      2,
                                                                      points,
                                                                      1,
                                                                      &point_count),
                            CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL);
    failures += expect_code("node outside cube",
                            cube_generate_dirichlet_points_from_nodes(TEST_L,
                                                                      TEST_V0,
                                                                      invalid_nodes,
                                                                      1,
                                                                      points,
                                                                      2,
                                                                      &point_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_regular_grid_matches_existing_generator();
    failures += test_top_refined_nodes_generate_surface_dirichlet();
    failures += test_duplicate_nodes_are_merged();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
