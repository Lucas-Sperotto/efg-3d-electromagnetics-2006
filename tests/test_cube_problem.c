#include "cube_problem.h"

#include <math.h>
#include <stdio.h>

#define TOLERANCE 1.0e-12

static int expect_code(const char *name, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s: expected code %d, got %d\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int expect_close(const char *name, double actual, double expected)
{
    if (fabs(actual - expected) > TOLERANCE) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int same_point(double ax,
                      double ay,
                      double az,
                      double bx,
                      double by,
                      double bz)
{
    return fabs(ax - bx) <= TOLERANCE &&
           fabs(ay - by) <= TOLERANCE &&
           fabs(az - bz) <= TOLERANCE;
}

static int find_node(const Node3D *nodes,
                     int node_count,
                     double x,
                     double y,
                     double z)
{
    for (int i = 0; i < node_count; ++i) {
        if (same_point(nodes[i].x, nodes[i].y, nodes[i].z, x, y, z)) {
            return 1;
        }
    }

    return 0;
}

static int find_dirichlet_point(const DirichletPoint *points,
                                int point_count,
                                double x,
                                double y,
                                double z)
{
    for (int i = 0; i < point_count; ++i) {
        if (same_point(points[i].x, points[i].y, points[i].z, x, y, z)) {
            return 1;
        }
    }

    return 0;
}

static int test_counts(void)
{
    int failures = 0;

    failures += expect_code("2x2x2 node count", cube_regular_node_count(2, 2, 2), 8);
    failures += expect_code("3x3x3 node count", cube_regular_node_count(3, 3, 3), 27);
    failures += expect_code("2x2x2 Dirichlet count", cube_dirichlet_point_count(2, 2, 2), 8);
    failures += expect_code("3x3x3 Dirichlet count", cube_dirichlet_point_count(3, 3, 3), 26);

    return failures;
}

static int test_generate_vertices(void)
{
    Node3D nodes[8];
    int node_count = 0;
    int failures = 0;

    failures += expect_code("generate 2x2x2 nodes",
                            cube_generate_regular_nodes(10.0, 2, 2, 2, nodes, 8, &node_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("2x2x2 generated count", node_count, 8);

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].x < -TOLERANCE || nodes[i].x > 10.0 + TOLERANCE ||
            nodes[i].y < -TOLERANCE || nodes[i].y > 10.0 + TOLERANCE ||
            nodes[i].z < -TOLERANCE || nodes[i].z > 10.0 + TOLERANCE) {
            fprintf(stderr, "node %d is outside [0,L]^3\n", i);
            ++failures;
        }

        failures += expect_close("initial support radius", nodes[i].support_radius, 0.0);
    }

    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                if (!find_node(nodes, node_count, 10.0 * ix, 10.0 * iy, 10.0 * iz)) {
                    fprintf(stderr, "missing vertex (%g,%g,%g)\n",
                            10.0 * ix,
                            10.0 * iy,
                            10.0 * iz);
                    ++failures;
                }
            }
        }
    }

    return failures;
}

static int test_generate_three_by_three_nodes(void)
{
    Node3D nodes[27];
    int node_count = 0;
    int failures = 0;

    failures += expect_code("generate 3x3x3 nodes",
                            cube_generate_regular_nodes(10.0, 3, 3, 3, nodes, 27, &node_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("3x3x3 generated count", node_count, 27);

    if (!find_node(nodes, node_count, 5.0, 5.0, 5.0)) {
        fprintf(stderr, "missing center node (5,5,5)\n");
        ++failures;
    }

    return failures;
}

static int test_generate_dirichlet_points(void)
{
    DirichletPoint points[26];
    int point_count = 0;
    int failures = 0;

    failures += expect_code("generate 3x3x3 Dirichlet points",
                            cube_generate_dirichlet_points(10.0,
                                                           10.0,
                                                           3,
                                                           3,
                                                           3,
                                                           points,
                                                           26,
                                                           &point_count),
                            CUBE_PROBLEM_OK);
    failures += expect_code("3x3x3 Dirichlet generated count", point_count, 26);

    for (int a = 0; a < point_count; ++a) {
        const int on_top = fabs(points[a].z - 10.0) <= TOLERANCE;
        const int on_surface = fabs(points[a].x) <= TOLERANCE ||
                               fabs(points[a].x - 10.0) <= TOLERANCE ||
                               fabs(points[a].y) <= TOLERANCE ||
                               fabs(points[a].y - 10.0) <= TOLERANCE ||
                               fabs(points[a].z) <= TOLERANCE ||
                               on_top;

        if (!on_surface) {
            fprintf(stderr, "Dirichlet point %d is not on the surface\n", a);
            ++failures;
        }

        if (on_top) {
            failures += expect_close("top face prescribed value", points[a].value, 10.0);
        } else {
            failures += expect_close("non-top surface prescribed value", points[a].value, 0.0);
        }

        for (int b = a + 1; b < point_count; ++b) {
            if (same_point(points[a].x,
                           points[a].y,
                           points[a].z,
                           points[b].x,
                           points[b].y,
                           points[b].z)) {
                fprintf(stderr, "duplicate Dirichlet point at indices %d and %d\n", a, b);
                ++failures;
            }
        }
    }

    if (find_dirichlet_point(points, point_count, 5.0, 5.0, 5.0)) {
        fprintf(stderr, "internal point (5,5,5) should not be a Dirichlet point\n");
        ++failures;
    }

    return failures;
}

static int test_error_paths(void)
{
    Node3D nodes[8];
    DirichletPoint points[8];
    int count = 0;
    int failures = 0;

    failures += expect_code("invalid node count grid",
                            cube_regular_node_count(1, 2, 2),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid Dirichlet count grid",
                            cube_dirichlet_point_count(2, 1, 2),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid L for nodes",
                            cube_generate_regular_nodes(0.0, 2, 2, 2, nodes, 8, &count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null node output",
                            cube_generate_regular_nodes(10.0, 2, 2, 2, NULL, 8, &count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null node_count output",
                            cube_generate_regular_nodes(10.0, 2, 2, 2, nodes, 8, NULL),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("insufficient node capacity",
                            cube_generate_regular_nodes(10.0, 2, 2, 2, nodes, 7, &count),
                            CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL);
    failures += expect_code("invalid L for Dirichlet",
                            cube_generate_dirichlet_points(0.0,
                                                           10.0,
                                                           2,
                                                           2,
                                                           2,
                                                           points,
                                                           8,
                                                           &count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null Dirichlet output",
                            cube_generate_dirichlet_points(10.0,
                                                           10.0,
                                                           2,
                                                           2,
                                                           2,
                                                           NULL,
                                                           8,
                                                           &count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null Dirichlet count output",
                            cube_generate_dirichlet_points(10.0,
                                                           10.0,
                                                           2,
                                                           2,
                                                           2,
                                                           points,
                                                           8,
                                                           NULL),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("insufficient Dirichlet capacity",
                            cube_generate_dirichlet_points(10.0,
                                                           10.0,
                                                           2,
                                                           2,
                                                           2,
                                                           points,
                                                           7,
                                                           &count),
                            CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_counts();
    failures += test_generate_vertices();
    failures += test_generate_three_by_three_nodes();
    failures += test_generate_dirichlet_points();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
