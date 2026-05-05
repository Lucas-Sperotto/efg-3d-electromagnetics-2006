#include "mls_diagnostic.h"
#include "cube_problem.h"
#include "support.h"

#include <stdio.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ helpers */

static void fail(const char *msg)
{
    fprintf(stderr, "FAIL: %s\n", msg);
    exit(1);
}

/*
 * Allocate and fill a regular node cloud for an NxNxN grid on [0,L]^3,
 * assign support radii using the article default (k=3, scale=1.5).
 * Returns the node count.
 */
static int make_regular_nodes(double L, int n_per_axis, Node3D *nodes, int max_nodes)
{
    int count = 0;
    int rc = cube_generate_regular_nodes(L, n_per_axis, n_per_axis, n_per_axis,
                                         nodes, max_nodes, &count);
    if (rc != CUBE_PROBLEM_OK || count != n_per_axis * n_per_axis * n_per_axis) {
        fail("make_regular_nodes: generation failed");
    }
    support_assign_article_default(nodes, count);
    return count;
}

/* ------------------------------------------------------------------ tests */

/*
 * A single isolated node must report INVALID (active_nodes = 1 < 4).
 */
static void test_single_node_is_invalid(void)
{
    Node3D node = node3d_make(5.0, 5.0, 5.0);
    MlsPointDiagnostic d;

    node.support_radius = 3.0;
    mls_diagnostic_at_point(&node, 1, 5.0, 5.0, 5.0, &d);

    if (d.status != MLS_POINT_INVALID) {
        fail("single node: expected INVALID status");
    }
    if (d.active_nodes != 1) {
        fail("single node: expected active_nodes = 1");
    }

    printf("PASS test_single_node_is_invalid: active=%d\n", d.active_nodes);
}

/*
 * A 3×3×3 regular grid evaluated at 3×3×3 cell centres must have no invalid
 * points (every cell centre sees the nearby nodes).
 */
static void test_3x3x3_no_invalid_points(void)
{
    Node3D nodes[27];
    MlsConnectivityStats stats;
    int n;

    n = make_regular_nodes(10.0, 3, nodes, 27);

    mls_connectivity_stats_uniform_grid(nodes, n,
                                        0.0, 10.0, 3,
                                        0.0, 10.0, 3,
                                        0.0, 10.0, 3,
                                        &stats);

    if (stats.n_total != 27) {
        fail("3x3x3: expected 27 evaluation points");
    }
    if (stats.n_invalid > 0) {
        fprintf(stderr, "FAIL 3x3x3: %d invalid points (min_nodes=%d)\n",
                stats.n_invalid, stats.min_nodes);
        exit(1);
    }

    printf("PASS test_3x3x3_no_invalid_points: "
           "min=%d max=%d mean=%.1f invalid=%d suspect=%d "
           "moment_fail=%d max_cond=%.2e\n",
           stats.min_nodes, stats.max_nodes, stats.mean_nodes,
           stats.n_invalid, stats.n_suspect, stats.n_moment_fail,
           stats.max_cond);
}

/*
 * A 5×5×5 regular grid must also have no invalid points at cell centres.
 */
static void test_5x5x5_no_invalid_points(void)
{
    Node3D nodes[125];
    MlsConnectivityStats stats;
    int n;

    n = make_regular_nodes(10.0, 5, nodes, 125);

    mls_connectivity_stats_uniform_grid(nodes, n,
                                        0.0, 10.0, 5,
                                        0.0, 10.0, 5,
                                        0.0, 10.0, 5,
                                        &stats);

    if (stats.n_invalid > 0) {
        fprintf(stderr, "FAIL 5x5x5: %d invalid points\n", stats.n_invalid);
        exit(1);
    }

    printf("PASS test_5x5x5_no_invalid_points: "
           "min=%d max=%d mean=%.1f max_cond=%.2e\n",
           stats.min_nodes, stats.max_nodes, stats.mean_nodes, stats.max_cond);
}

/*
 * stats_from_points must produce the same counts as the equivalent grid call.
 */
static void test_stats_from_points_matches_grid(void)
{
    Node3D nodes[27];
    MlsConnectivityStats sg;
    MlsConnectivityStats sp;
    double pts[27 * 3];
    int k = 0;
    int n;

    n = make_regular_nodes(10.0, 3, nodes, 27);

    for (int ix = 0; ix < 3; ++ix) {
        for (int iy = 0; iy < 3; ++iy) {
            for (int iz = 0; iz < 3; ++iz) {
                pts[k++] = (ix + 0.5) * 10.0 / 3.0;
                pts[k++] = (iy + 0.5) * 10.0 / 3.0;
                pts[k++] = (iz + 0.5) * 10.0 / 3.0;
            }
        }
    }

    mls_connectivity_stats_uniform_grid(nodes, n,
                                        0.0, 10.0, 3,
                                        0.0, 10.0, 3,
                                        0.0, 10.0, 3, &sg);
    mls_connectivity_stats_from_points(nodes, n, pts, 27, &sp);

    if (sg.n_invalid != sp.n_invalid ||
        sg.min_nodes != sp.min_nodes ||
        sg.max_nodes != sp.max_nodes) {
        fail("stats_from_points does not match uniform_grid");
    }

    printf("PASS test_stats_from_points_matches_grid\n");
}

/* ------------------------------------------------------------------ main */

int main(void)
{
    test_single_node_is_invalid();
    test_3x3x3_no_invalid_points();
    test_5x5x5_no_invalid_points();
    test_stats_from_points_matches_grid();
    return 0;
}
