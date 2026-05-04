#include "global_stiffness.h"

#include "stiffness.h"

#include <math.h>
#include <stdio.h>

#define TEST_NODE_COUNT 8
#define TEST_MATRIX_SIZE (TEST_NODE_COUNT * TEST_NODE_COUNT)
#define TEST_MAX_LOCAL_ENTRIES 512

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
    const double tolerance = 1.0e-8;

    if (fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%s: expected %.17g, got %.17g\n", name, expected, actual);
        return 1;
    }

    return 0;
}

static int validate_dense_matrix(const double K[TEST_MATRIX_SIZE])
{
    int failures = 0;

    for (int row = 0; row < TEST_NODE_COUNT; ++row) {
        double row_sum = 0.0;

        if (!isfinite(K[(row * TEST_NODE_COUNT) + row])) {
            fprintf(stderr, "diagonal %d is not finite\n", row);
            ++failures;
        }

        if (K[(row * TEST_NODE_COUNT) + row] < -1.0e-8) {
            fprintf(stderr, "diagonal %d is negative: %.17g\n",
                    row,
                    K[(row * TEST_NODE_COUNT) + row]);
            ++failures;
        }

        for (int col = 0; col < TEST_NODE_COUNT; ++col) {
            char label[64];
            const double value = K[(row * TEST_NODE_COUNT) + col];

            if (!isfinite(value)) {
                fprintf(stderr, "K[%d][%d] is not finite\n", row, col);
                ++failures;
            }

            row_sum += value;
            snprintf(label, sizeof(label), "global symmetry[%d][%d]", row, col);
            failures += expect_close(label, value, K[(col * TEST_NODE_COUNT) + row]);
        }

        failures += expect_close("global constant-field row sum", row_sum, 0.0);
    }

    return failures;
}

static void consolidate_local_entries(const StiffnessEntry *entries,
                                      int entry_count,
                                      double K[TEST_MATRIX_SIZE])
{
    for (int i = 0; i < TEST_MATRIX_SIZE; ++i) {
        K[i] = 0.0;
    }

    for (int entry = 0; entry < entry_count; ++entry) {
        K[(entries[entry].row_node * TEST_NODE_COUNT) + entries[entry].col_node] +=
            entries[entry].value;
    }
}

static int test_one_cell_matches_local_consolidation(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    StiffnessEntry entries[TEST_MAX_LOCAL_ENTRIES];
    double global_K[TEST_MATRIX_SIZE];
    double local_K[TEST_MATRIX_SIZE];
    int local_count = 0;
    int failures = 0;

    make_unit_cube_nodes(nodes);

    failures += expect_code("global one-cell assemble",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            1,
                                                            1,
                                                            1,
                                                            global_K),
                            GLOBAL_STIFFNESS_OK);

    failures += expect_code("local one-cell assemble",
                            stiffness_assemble_cell(nodes,
                                                    TEST_NODE_COUNT,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    entries,
                                                    TEST_MAX_LOCAL_ENTRIES,
                                                    &local_count),
                            STIFFNESS_OK);

    if (failures != 0) {
        return failures;
    }

    consolidate_local_entries(entries, local_count, local_K);

    for (int i = 0; i < TEST_MATRIX_SIZE; ++i) {
        failures += expect_close("global equals local", global_K[i], local_K[i]);
    }

    failures += validate_dense_matrix(global_K);

    return failures;
}

static int test_two_by_two_by_two_cells(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    double K[TEST_MATRIX_SIZE];
    int failures = 0;

    make_unit_cube_nodes(nodes);

    failures += expect_code("global 2x2x2 assemble",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            2,
                                                            2,
                                                            2,
                                                            K),
                            GLOBAL_STIFFNESS_OK);

    if (failures != 0) {
        return failures;
    }

    failures += validate_dense_matrix(K);

    return failures;
}

static int test_error_paths(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    double K[TEST_MATRIX_SIZE];
    int failures = 0;

    make_unit_cube_nodes(nodes);

    failures += expect_code("null nodes",
                            global_stiffness_assemble_dense(NULL,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            1,
                                                            1,
                                                            1,
                                                            K),
                            GLOBAL_STIFFNESS_INVALID_ARGUMENT);

    failures += expect_code("invalid domain",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            1.0,
                                                            0.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            1,
                                                            1,
                                                            1,
                                                            K),
                            GLOBAL_STIFFNESS_INVALID_DOMAIN);

    failures += expect_code("invalid cells",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0,
                                                            1,
                                                            1,
                                                            K),
                            GLOBAL_STIFFNESS_INVALID_CELLS);

    nodes[0].support_radius = 0.0;
    failures += expect_code("invalid support",
                            global_stiffness_assemble_dense(nodes,
                                                            TEST_NODE_COUNT,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            0.0,
                                                            1.0,
                                                            1,
                                                            1,
                                                            1,
                                                            K),
                            GLOBAL_STIFFNESS_INVALID_SUPPORT_RADIUS);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_one_cell_matches_local_consolidation();
    failures += test_two_by_two_by_two_cells();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
