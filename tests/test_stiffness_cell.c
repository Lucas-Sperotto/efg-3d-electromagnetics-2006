#include "stiffness.h"

#include <math.h>
#include <stdio.h>

#define TEST_NODE_COUNT 8
#define TEST_MAX_ENTRIES 512

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

static void consolidate_entries(const StiffnessEntry *entries,
                                int entry_count,
                                double matrix[TEST_NODE_COUNT][TEST_NODE_COUNT])
{
    for (int row = 0; row < TEST_NODE_COUNT; ++row) {
        for (int col = 0; col < TEST_NODE_COUNT; ++col) {
            matrix[row][col] = 0.0;
        }
    }

    for (int entry = 0; entry < entry_count; ++entry) {
        matrix[entries[entry].row_node][entries[entry].col_node] += entries[entry].value;
    }
}

static int validate_entries(const StiffnessEntry *entries, int entry_count)
{
    int failures = 0;

    if (entry_count <= 0) {
        fprintf(stderr, "expected at least one stiffness entry\n");
        return 1;
    }

    for (int i = 0; i < entry_count; ++i) {
        if (entries[i].row_node < 0 || entries[i].row_node >= TEST_NODE_COUNT ||
            entries[i].col_node < 0 || entries[i].col_node >= TEST_NODE_COUNT) {
            fprintf(stderr, "entry %d has invalid node pair (%d, %d)\n",
                    i,
                    entries[i].row_node,
                    entries[i].col_node);
            ++failures;
        }

        if (!isfinite(entries[i].value)) {
            fprintf(stderr, "entry %d is not finite\n", i);
            ++failures;
        }
    }

    return failures;
}

static int validate_matrix_properties(const double matrix[TEST_NODE_COUNT][TEST_NODE_COUNT])
{
    const double tolerance = 1.0e-8;
    int failures = 0;

    for (int row = 0; row < TEST_NODE_COUNT; ++row) {
        double row_sum = 0.0;

        if (matrix[row][row] < -tolerance) {
            fprintf(stderr, "diagonal %d is negative: %.17g\n", row, matrix[row][row]);
            ++failures;
        }

        for (int col = 0; col < TEST_NODE_COUNT; ++col) {
            char label[64];

            row_sum += matrix[row][col];
            snprintf(label, sizeof(label), "symmetry[%d][%d]", row, col);
            failures += expect_close(label, matrix[row][col], matrix[col][row]);
        }

        failures += expect_close("constant-field row sum", row_sum, 0.0);
    }

    return failures;
}

static int test_unit_cube_cell(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    StiffnessEntry entries[TEST_MAX_ENTRIES];
    double matrix[TEST_NODE_COUNT][TEST_NODE_COUNT];
    int entry_count = 0;
    int failures = 0;

    make_unit_cube_nodes(nodes);

    failures += expect_code("assemble unit cube",
                            stiffness_assemble_cell(nodes,
                                                    TEST_NODE_COUNT,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    entries,
                                                    TEST_MAX_ENTRIES,
                                                    &entry_count),
                            STIFFNESS_OK);

    if (failures != 0) {
        return failures;
    }

    failures += validate_entries(entries, entry_count);
    consolidate_entries(entries, entry_count, matrix);
    failures += validate_matrix_properties((const double(*)[TEST_NODE_COUNT])matrix);

    return failures;
}

static int test_error_paths(void)
{
    Node3D nodes[TEST_NODE_COUNT];
    StiffnessEntry entries[TEST_MAX_ENTRIES];
    int entry_count = 0;
    int failures = 0;

    make_unit_cube_nodes(nodes);

    failures += expect_code("null nodes",
                            stiffness_assemble_cell(NULL,
                                                    TEST_NODE_COUNT,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    entries,
                                                    TEST_MAX_ENTRIES,
                                                    &entry_count),
                            STIFFNESS_INVALID_ARGUMENT);

    failures += expect_code("invalid cell",
                            stiffness_assemble_cell(nodes,
                                                    TEST_NODE_COUNT,
                                                    1.0,
                                                    0.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    entries,
                                                    TEST_MAX_ENTRIES,
                                                    &entry_count),
                            STIFFNESS_INVALID_CELL);

    failures += expect_code("insufficient capacity",
                            stiffness_assemble_cell(nodes,
                                                    TEST_NODE_COUNT,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    1.0,
                                                    entries,
                                                    1,
                                                    &entry_count),
                            STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL);

    return failures;
}

int main(void)
{
    int failures = 0;

    failures += test_unit_cube_cell();
    failures += test_error_paths();

    return failures == 0 ? 0 : 1;
}
