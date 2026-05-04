#include "mls.h"

#include <math.h>
#include <stdio.h>

static void make_unit_cube_nodes(Node3D nodes[8])
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

static int expect_close(const char *name, double actual, double expected)
{
    const double tolerance = 1.0e-10;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int test_partition_at(double x, double y, double z)
{
    Node3D nodes[8];
    MlsShapeValue values[8];
    int value_count = 0;
    double sum_phi = 0.0;
    int status;

    make_unit_cube_nodes(nodes);

    status = mls_linear3d_shape_functions(nodes, 8, x, y, z, values, 8, &value_count);
    if (status != MLS_OK) {
        fprintf(stderr, "mls failed at (%.17g, %.17g, %.17g): %d\n", x, y, z, status);
        return 1;
    }

    for (int i = 0; i < value_count; ++i) {
        sum_phi += values[i].phi;
    }

    return expect_close("partition of unity", sum_phi, 1.0);
}

int main(void)
{
    int failures = 0;

    failures += test_partition_at(0.5, 0.5, 0.5);
    failures += test_partition_at(0.25, 0.5, 0.75);
    failures += test_partition_at(0.8, 0.2, 0.4);

    return failures == 0 ? 0 : 1;
}
