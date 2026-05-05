#include "cube_problem.h"

#include <math.h>
#include <stdio.h>

#define TEST_L 10.0
#define TEST_BASE_N 3
#define TEST_REFINE_N 5
#define TEST_MAX_NODES 1024
#define TEST_TOLERANCE 1.0e-10

static int expect_code(const char *name, int actual, int expected)
{
    if (actual != expected) {
        fprintf(stderr, "%s: expected code %d, got %d\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int same_point(const Node3D *a, const Node3D *b)
{
    return fabs(a->x - b->x) <= TEST_TOLERANCE &&
           fabs(a->y - b->y) <= TEST_TOLERANCE &&
           fabs(a->z - b->z) <= TEST_TOLERANCE;
}

static int node_near(const Node3D *node, double x, double y, double z, double radius)
{
    const double dx = node->x - x;
    const double dy = node->y - y;
    const double dz = node->z - z;
    const double distance_squared = (dx * dx) + (dy * dy) + (dz * dz);

    return distance_squared <= radius * radius;
}

static int count_nodes_near_z(const Node3D *nodes, int node_count, double z, double tolerance)
{
    int count = 0;

    for (int i = 0; i < node_count; ++i) {
        if (fabs(nodes[i].z - z) <= tolerance) {
            ++count;
        }
    }

    return count;
}

static int test_top_refined_cloud_properties(void)
{
    Node3D nodes[TEST_MAX_NODES];
    int node_count = 0;
    const int base_count = TEST_BASE_N * TEST_BASE_N * TEST_BASE_N;
    int top_like_count = 0;
    int middle_like_count = 0;
    int upper_corner_hits = 0;
    int near_upper_corner_hits = 0;
    int failures = 0;

    failures += expect_code("generate top-refined nodes",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            TEST_BASE_N,
                                                            TEST_REFINE_N,
                                                            nodes,
                                                            TEST_MAX_NODES,
                                                            &node_count),
                            CUBE_PROBLEM_OK);
    if (failures != 0) {
        return failures;
    }

    if (node_count <= base_count) {
        fprintf(stderr, "expected more than %d base nodes, got %d\n", base_count, node_count);
        ++failures;
    }

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].x < -TEST_TOLERANCE || nodes[i].x > TEST_L + TEST_TOLERANCE ||
            nodes[i].y < -TEST_TOLERANCE || nodes[i].y > TEST_L + TEST_TOLERANCE ||
            nodes[i].z < -TEST_TOLERANCE || nodes[i].z > TEST_L + TEST_TOLERANCE) {
            fprintf(stderr, "node %d is outside [0,L]^3\n", i);
            ++failures;
        }

        if (fabs(nodes[i].support_radius) > TEST_TOLERANCE) {
            fprintf(stderr, "node %d support_radius should start at zero\n", i);
            ++failures;
        }

        if (nodes[i].z >= 0.75 * TEST_L - TEST_TOLERANCE) {
            ++top_like_count;
        }

        if (fabs(nodes[i].z - (0.5 * TEST_L)) <= TEST_TOLERANCE) {
            ++middle_like_count;
        }

        if (node_near(&nodes[i], 0.0, 0.0, TEST_L, TEST_TOLERANCE) ||
            node_near(&nodes[i], TEST_L, 0.0, TEST_L, TEST_TOLERANCE) ||
            node_near(&nodes[i], 0.0, TEST_L, TEST_L, TEST_TOLERANCE) ||
            node_near(&nodes[i], TEST_L, TEST_L, TEST_L, TEST_TOLERANCE)) {
            ++upper_corner_hits;
        }

        if (node_near(&nodes[i], 0.0, 0.0, TEST_L, 2.0) ||
            node_near(&nodes[i], TEST_L, 0.0, TEST_L, 2.0) ||
            node_near(&nodes[i], 0.0, TEST_L, TEST_L, 2.0) ||
            node_near(&nodes[i], TEST_L, TEST_L, TEST_L, 2.0)) {
            ++near_upper_corner_hits;
        }

        for (int j = i + 1; j < node_count; ++j) {
            if (same_point(&nodes[i], &nodes[j])) {
                fprintf(stderr, "duplicate node at indices %d and %d\n", i, j);
                ++failures;
            }
        }
    }

    if (top_like_count <= middle_like_count) {
        fprintf(stderr,
                "expected more top-refined nodes than middle nodes, got top=%d middle=%d\n",
                top_like_count,
                middle_like_count);
        ++failures;
    }

    if (upper_corner_hits != 4) {
        fprintf(stderr, "expected exactly 4 upper corner nodes, got %d\n", upper_corner_hits);
        ++failures;
    }

    if (near_upper_corner_hits <= upper_corner_hits) {
        fprintf(stderr,
                "expected additional nodes near upper corners, got %d total hits\n",
                near_upper_corner_hits);
        ++failures;
    }

    if (count_nodes_near_z(nodes, node_count, TEST_L, TEST_TOLERANCE) <=
        count_nodes_near_z(nodes, node_count, 0.5 * TEST_L, TEST_TOLERANCE)) {
        fprintf(stderr, "expected denser z=L layer than z=L/2 layer\n");
        ++failures;
    }

    return failures;
}

static int test_error_paths(void)
{
    Node3D nodes[TEST_MAX_NODES];
    int node_count = 0;
    int failures = 0;

    failures += expect_code("invalid L",
                            cube_generate_top_refined_nodes(0.0,
                                                            TEST_BASE_N,
                                                            TEST_REFINE_N,
                                                            nodes,
                                                            TEST_MAX_NODES,
                                                            &node_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid base_n",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            1,
                                                            TEST_REFINE_N,
                                                            nodes,
                                                            TEST_MAX_NODES,
                                                            &node_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("invalid refine_n",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            TEST_BASE_N,
                                                            1,
                                                            nodes,
                                                            TEST_MAX_NODES,
                                                            &node_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null nodes",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            TEST_BASE_N,
                                                            TEST_REFINE_N,
                                                            NULL,
                                                            TEST_MAX_NODES,
                                                            &node_count),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("null node_count",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            TEST_BASE_N,
                                                            TEST_REFINE_N,
                                                            nodes,
                                                            TEST_MAX_NODES,
                                                            NULL),
                            CUBE_PROBLEM_INVALID_ARGUMENT);
    failures += expect_code("insufficient capacity",
                            cube_generate_top_refined_nodes(TEST_L,
                                                            TEST_BASE_N,
                                                            TEST_REFINE_N,
                                                            nodes,
                                                            4,
                                                            &node_count),
                            CUBE_PROBLEM_OUTPUT_CAPACITY_TOO_SMALL);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_top_refined_cloud_properties();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
