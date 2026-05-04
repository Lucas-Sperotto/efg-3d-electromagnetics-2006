#include "analytical.h"

#include <stdio.h>
#include <stdlib.h>

/*
 * Generate the analytical yz slice used as a reference for Fig. 3.
 *
 * This app writes a CSV file for x = 5.33 using only the analytical solution
 * from docs/03_resultados_numericos.md, section 3.1. It intentionally does not
 * implement EFG, assemble matrices, or solve any linear system.
 */

static int parse_int_at_least(const char *text, int minimum, int fallback)
{
    char *end = NULL;
    const long value = strtol(text, &end, 10);

    if (end == text || *end != '\0' || value < minimum || value > 1000000L) {
        return fallback;
    }

    return (int)value;
}

int main(int argc, char **argv)
{
    const char *output_path = "data/output/analytical_slice_x5_33.csv";
    const double L = 10.0;
    const double V0 = 10.0;
    const double x_slice = 5.33;
    int grid_points = 101;
    int max_terms = 25;

    if (argc > 1) {
        output_path = argv[1];
    }

    if (argc > 2) {
        grid_points = parse_int_at_least(argv[2], 2, grid_points);
    }

    if (argc > 3) {
        max_terms = parse_int_at_least(argv[3], 1, max_terms);
    }

    FILE *file = fopen(output_path, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open output CSV: %s\n", output_path);
        return 1;
    }

    fprintf(file, "x,y,z,potential\n");

    for (int z_index = 0; z_index < grid_points; ++z_index) {
        const double z = (L * (double)z_index) / (double)(grid_points - 1);

        for (int y_index = 0; y_index < grid_points; ++y_index) {
            const double y = (L * (double)y_index) / (double)(grid_points - 1);
            const double potential =
                analytical_potential_cube(x_slice, y, z, L, V0, max_terms);

            fprintf(file, "%.17g,%.17g,%.17g,%.17g\n", x_slice, y, z, potential);
        }
    }

    if (fclose(file) != 0) {
        fprintf(stderr, "Could not close output CSV: %s\n", output_path);
        return 1;
    }

    return 0;
}
