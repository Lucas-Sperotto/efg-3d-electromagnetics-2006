#include "dirichlet.h"

#include <math.h>
#include <stdio.h>

#define TEST_NODE_COUNT 8
#define TEST_POINT_COUNT 3

static void make_unit_cube_nodes(Node3D nodes[TEST_NODE_COUNT])
{
    const double support_radius = 2.0;

    nodes[0] = (Node3D){0.0, 0.0, 0.0, support_radius};
    nodes[1] = (Node3D){1.0, 0.0, 0.0, support_radius};
    nodes[2] = (Node3D){0.0, 1.0, 0.0, support_radius};
    nodes[3] = (Node3D){1.0, 1.0, 0.0, support_radius};
    nodes[4] = (Node3D){0.0, 0.0, 1.0, support_radius};
    nodes[5] = (Node3D){1.0, 0.0, 1.0, support_radius};
    nodes[6] = (Node3D){0.0, 1.0, 1.0, support_radius};
    nodes[7] = (Node3D){1.0, 1.0, 1.0, support_radius};
}

static void make_dirichlet_points(DirichletPoint points[TEST_POINT_COUNT])
{
    points[0] = (DirichletPoint){0.5, 0.5, 0.0, 0.0};
    points[1] = (DirichletPoint){0.5, 0.5, 1.0, 10.0};
    points[2] = (DirichletPoint){0.0, 0.5, 0.5, 0.0};
}

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
    const double tolerance = 1.0e-9;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static double apply_constraint_row(const double *G, const double *nodal_values, int row)
{
    double value = 0.0;

    for (int node = 0; node < TEST_NODE_COUNT; ++node) {
        value += G[(row * TEST_NODE_COUNT) + node] * nodal_values[node];
    }

    return value;
}

static int test_constraint_reproduction(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    DirichletPoint points[TEST_POINT_COUNT];
    double G[TEST_POINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_POINT_COUNT];
    double constant_values[TEST_NODE_COUNT];
    double x_values[TEST_NODE_COUNT];
    double y_values[TEST_NODE_COUNT];
    double z_values[TEST_NODE_COUNT];
    int failures = 0;

    make_unit_cube_nodes(nodes);
    make_dirichlet_points(points);

    for (int i = 0; i < TEST_NODE_COUNT; ++i) {
        constant_values[i] = 3.25;
        x_values[i] = nodes[i].x;
        y_values[i] = nodes[i].y;
        z_values[i] = nodes[i].z;
    }

    failures += expect_code("assemble constraints",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_POINT_COUNT,
                                                                 G,
                                                                 q),
                            DIRICHLET_OK);
    if (failures != 0) {
        return failures;
    }

    for (int point = 0; point < TEST_POINT_COUNT; ++point) {
        double row_sum = 0.0;
        char label[64];

        snprintf(label, sizeof(label), "q[%d]", point);
        failures += expect_close(label, q[point], points[point].value);

        for (int node = 0; node < TEST_NODE_COUNT; ++node) {
            row_sum += G[(point * TEST_NODE_COUNT) + node];
        }

        snprintf(label, sizeof(label), "row sum[%d]", point);
        failures += expect_close(label, row_sum, 1.0);

        snprintf(label, sizeof(label), "constant reproduction[%d]", point);
        failures += expect_close(label, apply_constraint_row(G, constant_values, point), 3.25);

        snprintf(label, sizeof(label), "x reproduction[%d]", point);
        failures += expect_close(label, apply_constraint_row(G, x_values, point), points[point].x);

        snprintf(label, sizeof(label), "y reproduction[%d]", point);
        failures += expect_close(label, apply_constraint_row(G, y_values, point), points[point].y);

        snprintf(label, sizeof(label), "z reproduction[%d]", point);
        failures += expect_close(label, apply_constraint_row(G, z_values, point), points[point].z);
    }

    return failures;
}

static int test_error_paths(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    DirichletPoint points[TEST_POINT_COUNT];
    double G[TEST_POINT_COUNT * TEST_NODE_COUNT];
    double q[TEST_POINT_COUNT];
    int failures = 0;

    make_unit_cube_nodes(nodes);
    make_dirichlet_points(points);

    failures += expect_code("null nodes",
                            dirichlet_assemble_constraints_dense(NULL,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_POINT_COUNT,
                                                                 G,
                                                                 q),
                            DIRICHLET_INVALID_ARGUMENT);
    failures += expect_code("null points",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 NULL,
                                                                 TEST_POINT_COUNT,
                                                                 G,
                                                                 q),
                            DIRICHLET_INVALID_ARGUMENT);
    failures += expect_code("null G",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_POINT_COUNT,
                                                                 NULL,
                                                                 q),
                            DIRICHLET_INVALID_ARGUMENT);
    failures += expect_code("null q",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_POINT_COUNT,
                                                                 G,
                                                                 NULL),
                            DIRICHLET_INVALID_ARGUMENT);

    nodes[0].support_radius = 0.0;
    failures += expect_code("invalid support",
                            dirichlet_assemble_constraints_dense(nodes,
                                                                 TEST_NODE_COUNT,
                                                                 points,
                                                                 TEST_POINT_COUNT,
                                                                 G,
                                                                 q),
                            DIRICHLET_INVALID_SUPPORT_RADIUS);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_constraint_reproduction();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
