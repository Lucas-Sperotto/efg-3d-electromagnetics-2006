#ifndef GLOBAL_STIFFNESS_H
#define GLOBAL_STIFFNESS_H

#include "nodes.h"

/*
 * Dense global stiffness assembly.
 *
 * This module is intentionally didactic and dense. It is useful for validating
 * the local EFG stiffness contributions before introducing sparse storage,
 * Dirichlet constraints, Lagrange multipliers, right-hand sides, or GMRES.
 */

#define GLOBAL_STIFFNESS_OK 0
#define GLOBAL_STIFFNESS_INVALID_ARGUMENT 1
#define GLOBAL_STIFFNESS_INVALID_DOMAIN 2
#define GLOBAL_STIFFNESS_INVALID_CELLS 3
#define GLOBAL_STIFFNESS_INVALID_SUPPORT_RADIUS 4
#define GLOBAL_STIFFNESS_ALLOCATION_FAILED 5
#define GLOBAL_STIFFNESS_LOCAL_ASSEMBLY_FAILED 6

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
                                    double *K);

#endif
