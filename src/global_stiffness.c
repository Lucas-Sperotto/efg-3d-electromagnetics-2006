#include "global_stiffness.h"

#include "gauss.h"
#include "stiffness.h"

#include <stddef.h>
#include <stdlib.h>

static void zero_dense_matrix(double *K, int node_count)
{
    for (int i = 0; i < node_count * node_count; ++i) {
        K[i] = 0.0;
    }
}

/*
 * Assemble a dense global stiffness matrix from regular cubic cells.
 *
 * This is a didactic dense assembly before sparse data structures. It sums the
 * local entries returned by stiffness_assemble_cell(), including repeated
 * entries from Gauss points. It assumes epsilon = 1 and does not apply
 * Dirichlet conditions, Lagrange multipliers, load vectors, or GMRES.
 */
int global_stiffness_assemble_dense(const Node3D *nodes,
                                    int node_count,
                                    double xmin,
                                    double xmax,
                                    double ymin,
                                    double ymax,
                                    double zmin,
                                    double zmax,
                                    int nx_cells,
                                    int ny_cells,
                                    int nz_cells,
                                    double *K)
{
    int max_local_entries;
    StiffnessEntry *local_entries = NULL;
    double dx;
    double dy;
    double dz;

    if (nodes == NULL || K == NULL || node_count <= 0) {
        return GLOBAL_STIFFNESS_INVALID_ARGUMENT;
    }

    zero_dense_matrix(K, node_count);

    if (xmin >= xmax || ymin >= ymax || zmin >= zmax) {
        return GLOBAL_STIFFNESS_INVALID_DOMAIN;
    }

    if (nx_cells <= 0 || ny_cells <= 0 || nz_cells <= 0) {
        return GLOBAL_STIFFNESS_INVALID_CELLS;
    }

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].support_radius <= 0.0) {
            return GLOBAL_STIFFNESS_INVALID_SUPPORT_RADIUS;
        }
    }

    max_local_entries = GAUSS_LEGENDRE_ORDER2_CUBE_POINTS * node_count * node_count;
    dx = (xmax - xmin) / (double)nx_cells;
    dy = (ymax - ymin) / (double)ny_cells;
    dz = (zmax - zmin) / (double)nz_cells;

    local_entries = malloc((size_t)max_local_entries * sizeof(local_entries[0]));
    if (local_entries == NULL) {
        return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
    }

    for (int ix = 0; ix < nx_cells; ++ix) {
        const double cell_xmin = xmin + ((double)ix * dx);
        const double cell_xmax = cell_xmin + dx;

        for (int iy = 0; iy < ny_cells; ++iy) {
            const double cell_ymin = ymin + ((double)iy * dy);
            const double cell_ymax = cell_ymin + dy;

            for (int iz = 0; iz < nz_cells; ++iz) {
                int local_count = 0;
                const double cell_zmin = zmin + ((double)iz * dz);
                const double cell_zmax = cell_zmin + dz;
                const int local_status = stiffness_assemble_cell(nodes,
                                                                 node_count,
                                                                 cell_xmin,
                                                                 cell_xmax,
                                                                 cell_ymin,
                                                                 cell_ymax,
                                                                 cell_zmin,
                                                                 cell_zmax,
                                                                 local_entries,
                                                                 max_local_entries,
                                                                 &local_count);

                if (local_status != STIFFNESS_OK) {
                    free(local_entries);
                    zero_dense_matrix(K, node_count);
                    return GLOBAL_STIFFNESS_LOCAL_ASSEMBLY_FAILED;
                }

                for (int entry = 0; entry < local_count; ++entry) {
                    const int row = local_entries[entry].row_node;
                    const int col = local_entries[entry].col_node;

                    K[(row * node_count) + col] += local_entries[entry].value;
                }
            }
        }
    }

    free(local_entries);
    return GLOBAL_STIFFNESS_OK;
}
