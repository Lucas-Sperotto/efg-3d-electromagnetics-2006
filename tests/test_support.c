#include "support.h"

#include <math.h>
#include <stdio.h>

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
    const double tolerance = 1.0e-12;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int test_distance(void)
{
    const Node3D a = node3d_make(0.0, 0.0, 0.0);
    const Node3D b = node3d_make(1.0, 2.0, 2.0);

    return expect_close("euclidean distance", node3d_distance(&a, &b), 3.0);
}

static int test_kth_neighbor_distance(void)
{
    const Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0, 0.0},
        {3.0, 0.0, 0.0, 0.0},
    };
    double distance = 0.0;
    int failures = 0;

    failures += expect_code("third neighbor status",
                            support_kth_neighbor_distance(nodes, 4, 0, 3, &distance),
                            SUPPORT_OK);
    failures += expect_close("third neighbor distance", distance, 3.0);

    return failures;
}

static int test_article_default_support_radius(void)
{
    Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0, 0.0},
        {3.0, 0.0, 0.0, 0.0},
    };
    int failures = 0;

    failures += expect_code("article support status",
                            support_assign_article_default(nodes, 4),
                            SUPPORT_OK);
    failures += expect_close("node 0 support", nodes[0].support_radius, 4.5);
    failures += expect_close("node 1 support", nodes[1].support_radius, 3.0);
    failures += expect_close("node 2 support", nodes[2].support_radius, 3.0);
    failures += expect_close("node 3 support", nodes[3].support_radius, 4.5);

    return failures;
}

static int test_not_enough_neighbors(void)
{
    Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0, 0.0},
    };

    return expect_code("not enough neighbors",
                       support_assign_article_default(nodes, 3),
                       SUPPORT_NOT_ENOUGH_NEIGHBORS);
}

static int test_invalid_arguments(void)
{
    Node3D nodes[] = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0, 0.0},
        {3.0, 0.0, 0.0, 0.0},
    };
    double distance = 0.0;
    int failures = 0;

    failures += expect_code("null nodes",
                            support_kth_neighbor_distance(NULL, 4, 0, 3, &distance),
                            SUPPORT_INVALID_ARGUMENT);
    failures += expect_code("null distance",
                            support_kth_neighbor_distance(nodes, 4, 0, 3, NULL),
                            SUPPORT_INVALID_ARGUMENT);
    failures += expect_code("bad index",
                            support_kth_neighbor_distance(nodes, 4, 4, 3, &distance),
                            SUPPORT_INVALID_ARGUMENT);
    failures += expect_code("bad k",
                            support_kth_neighbor_distance(nodes, 4, 0, 0, &distance),
                            SUPPORT_INVALID_ARGUMENT);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_distance();
    failures += test_kth_neighbor_distance();
    failures += test_article_default_support_radius();
    failures += test_not_enough_neighbors();
    failures += test_invalid_arguments();

    return failures == 0 ? 0 : 1;
}
