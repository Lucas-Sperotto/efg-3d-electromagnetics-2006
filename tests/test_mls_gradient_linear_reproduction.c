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
    const double tolerance = 1.0e-9;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int test_gradient_linear_reproduction_at(double x, double y, double z)
{
    Node3D nodes[8];
    MlsShapeGradientValue values[8];
    int value_count = 0;
    double reproduced_x = 0.0;
    double reproduced_y = 0.0;
    double reproduced_z = 0.0;
    double grad_x[3] = {0.0, 0.0, 0.0};
    double grad_y[3] = {0.0, 0.0, 0.0};
    double grad_z[3] = {0.0, 0.0, 0.0};
    int failures = 0;

    make_unit_cube_nodes(nodes);

    if (mls_linear3d_shape_functions_with_gradients(nodes,
                                                    8,
                                                    x,
                                                    y,
                                                    z,
                                                    values,
                                                    8,
                                                    &value_count) != MLS_OK) {
        fprintf(stderr, "mls gradient failed at (%.17g, %.17g, %.17g)\n", x, y, z);
        return 1;
    }

    for (int i = 0; i < value_count; ++i) {
        const Node3D *node = &nodes[values[i].node_index];

        reproduced_x += values[i].phi * node->x;
        reproduced_y += values[i].phi * node->y;
        reproduced_z += values[i].phi * node->z;

        grad_x[0] += node->x * values[i].dphi_dx;
        grad_x[1] += node->x * values[i].dphi_dy;
        grad_x[2] += node->x * values[i].dphi_dz;

        grad_y[0] += node->y * values[i].dphi_dx;
        grad_y[1] += node->y * values[i].dphi_dy;
        grad_y[2] += node->y * values[i].dphi_dz;

        grad_z[0] += node->z * values[i].dphi_dx;
        grad_z[1] += node->z * values[i].dphi_dy;
        grad_z[2] += node->z * values[i].dphi_dz;
    }

    failures += expect_close("reproduce x", reproduced_x, x);
    failures += expect_close("reproduce y", reproduced_y, y);
    failures += expect_close("reproduce z", reproduced_z, z);

    failures += expect_close("grad x wrt x", grad_x[0], 1.0);
    failures += expect_close("grad x wrt y", grad_x[1], 0.0);
    failures += expect_close("grad x wrt z", grad_x[2], 0.0);

    failures += expect_close("grad y wrt x", grad_y[0], 0.0);
    failures += expect_close("grad y wrt y", grad_y[1], 1.0);
    failures += expect_close("grad y wrt z", grad_y[2], 0.0);

    failures += expect_close("grad z wrt x", grad_z[0], 0.0);
    failures += expect_close("grad z wrt y", grad_z[1], 0.0);
    failures += expect_close("grad z wrt z", grad_z[2], 1.0);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_gradient_linear_reproduction_at(0.5, 0.5, 0.5);
    failures += test_gradient_linear_reproduction_at(0.25, 0.5, 0.75);
    failures += test_gradient_linear_reproduction_at(0.8, 0.2, 0.4);

    return failures == 0 ? 0 : 1;
}
