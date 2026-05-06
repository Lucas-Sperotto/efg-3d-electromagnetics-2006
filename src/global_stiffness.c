#include "global_stiffness.h"

#include "gauss.h"
#include "sparse_matrix.h"
#include "stiffness.h"

#include <stddef.h>
#include <stdlib.h>

static void zero_dense_matrix(double *K, int node_count)
{
    for (int i = 0; i < node_count * node_count; ++i) {
        K[i] = 0.0;
    }
}

static int initial_local_entry_capacity(int node_count, int quadrature_order)
{
    const int qpts = quadrature_order * quadrature_order * quadrature_order;
    const int active_cap = (node_count < 64) ? node_count : 64;

    return qpts * active_cap * active_cap;
}

static int grow_local_entries(StiffnessEntry **entries, int *capacity)
{
    StiffnessEntry *grown = NULL;
    int new_capacity;

    if (entries == NULL || *entries == NULL || capacity == NULL ||
        *capacity <= 0) {
        return GLOBAL_STIFFNESS_INVALID_ARGUMENT;
    }

    new_capacity = *capacity * 2;
    if (new_capacity <= *capacity) {
        return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
    }

    grown = realloc(*entries, (size_t)new_capacity * sizeof((*entries)[0]));
    if (grown == NULL) {
        return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
    }

    *entries = grown;
    *capacity = new_capacity;
    return GLOBAL_STIFFNESS_OK;
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
    return global_stiffness_assemble_dense_with_order(nodes, node_count,
                                                      xmin, xmax,
                                                      ymin, ymax,
                                                      zmin, zmax,
                                                      nx_cells, ny_cells,
                                                      nz_cells,
                                                      2,
                                                      K);
}

/*
 * Assemble the dense global stiffness matrix for Eq. (8).
 *
 * Article equation: Eq. (8), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: global K_ij obtained by summing all cell integrals of
 * grad(Phi_i) dot grad(Phi_j), with an explicit Gauss order for integration
 * sensitivity studies.
 */
int global_stiffness_assemble_dense_with_order(const Node3D *nodes,
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
                                               int quadrature_order,
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

    if (quadrature_order < GAUSS_LEGENDRE_MIN_ORDER ||
        quadrature_order > GAUSS_LEGENDRE_MAX_ORDER) {
        return GLOBAL_STIFFNESS_INVALID_ARGUMENT;
    }

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].support_radius <= 0.0) {
            return GLOBAL_STIFFNESS_INVALID_SUPPORT_RADIUS;
        }
    }

    max_local_entries = initial_local_entry_capacity(node_count,
                                                     quadrature_order);
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
                int local_status;
                const double cell_zmin = zmin + ((double)iz * dz);
                const double cell_zmax = cell_zmin + dz;

                do {
                    local_status =
                        stiffness_assemble_cell_with_order(nodes,
                                                           node_count,
                                                           cell_xmin,
                                                           cell_xmax,
                                                           cell_ymin,
                                                           cell_ymax,
                                                           cell_zmin,
                                                           cell_zmax,
                                                           quadrature_order,
                                                           local_entries,
                                                           max_local_entries,
                                                           &local_count);
                    if (local_status == STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL) {
                        if (grow_local_entries(&local_entries,
                                               &max_local_entries)
                            != GLOBAL_STIFFNESS_OK) {
                            free(local_entries);
                            zero_dense_matrix(K, node_count);
                            return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
                        }
                    }
                } while (local_status == STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL);

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

/*
 * Sparse COO variant: reuses the same cell-loop logic; each local entry is
 * appended to the COO instead of being scatter-added into a dense array.
 * Duplicate (row,col) pairs accumulate naturally and are summed at CSR
 * conversion time (sparse_coo_to_csr).
 */
int global_stiffness_assemble_sparse_coo(const Node3D *nodes,
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
                                          SparseCOO *coo)
{
    return global_stiffness_assemble_sparse_coo_with_order(nodes, node_count,
                                                           xmin, xmax,
                                                           ymin, ymax,
                                                           zmin, zmax,
                                                           nx_cells, ny_cells,
                                                           nz_cells,
                                                           2,
                                                           coo);
}

/*
 * Assemble the sparse COO global stiffness matrix for Eq. (8).
 *
 * Article equation: Eq. (8), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/03_resultados_numericos.md, section 3.3.
 * Mathematical meaning: sparse accumulation of the same global K_ij integral,
 * leaving duplicate COO entries to be summed during CSR conversion.
 */
int global_stiffness_assemble_sparse_coo_with_order(const Node3D *nodes,
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
                                                    int quadrature_order,
                                                    SparseCOO *coo)
{
    int max_local_entries;
    StiffnessEntry *local_entries = NULL;
    double dx;
    double dy;
    double dz;

    if (nodes == NULL || coo == NULL || node_count <= 0) {
        return GLOBAL_STIFFNESS_INVALID_ARGUMENT;
    }

    if (xmin >= xmax || ymin >= ymax || zmin >= zmax) {
        return GLOBAL_STIFFNESS_INVALID_DOMAIN;
    }

    if (nx_cells <= 0 || ny_cells <= 0 || nz_cells <= 0) {
        return GLOBAL_STIFFNESS_INVALID_CELLS;
    }

    if (quadrature_order < GAUSS_LEGENDRE_MIN_ORDER ||
        quadrature_order > GAUSS_LEGENDRE_MAX_ORDER) {
        return GLOBAL_STIFFNESS_INVALID_ARGUMENT;
    }

    for (int i = 0; i < node_count; ++i) {
        if (nodes[i].support_radius <= 0.0) {
            return GLOBAL_STIFFNESS_INVALID_SUPPORT_RADIUS;
        }
    }

    max_local_entries = initial_local_entry_capacity(node_count,
                                                     quadrature_order);
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
                int local_status;
                const double cell_zmin = zmin + ((double)iz * dz);
                const double cell_zmax = cell_zmin + dz;

                do {
                    local_status =
                        stiffness_assemble_cell_with_order(nodes,
                                                           node_count,
                                                           cell_xmin,
                                                           cell_xmax,
                                                           cell_ymin,
                                                           cell_ymax,
                                                           cell_zmin,
                                                           cell_zmax,
                                                           quadrature_order,
                                                           local_entries,
                                                           max_local_entries,
                                                           &local_count);
                    if (local_status == STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL) {
                        if (grow_local_entries(&local_entries,
                                               &max_local_entries)
                            != GLOBAL_STIFFNESS_OK) {
                            free(local_entries);
                            return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
                        }
                    }
                } while (local_status == STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL);

                if (local_status != STIFFNESS_OK) {
                    free(local_entries);
                    return GLOBAL_STIFFNESS_LOCAL_ASSEMBLY_FAILED;
                }

                for (int entry = 0; entry < local_count; ++entry) {
                    const int row = local_entries[entry].row_node;
                    const int col = local_entries[entry].col_node;
                    const int rc  = sparse_coo_add(coo, row, col,
                                                   local_entries[entry].value);
                    if (rc != SPARSE_OK) {
                        free(local_entries);
                        return GLOBAL_STIFFNESS_ALLOCATION_FAILED;
                    }
                }
            }
        }
    }

    free(local_entries);
    return GLOBAL_STIFFNESS_OK;
}
