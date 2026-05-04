#ifndef STIFFNESS_H
#define STIFFNESS_H

#include "nodes.h"

/*
 * Didactic local stiffness assembly for one cubic integration cell.
 *
 * This is the first local contribution associated with the weak form of the
 * Laplace/Poisson problem. It computes only grad(Phi_i) dot grad(Phi_j) with
 * epsilon = 1. It does not assemble a global matrix, impose Dirichlet
 * conditions, use Lagrange multipliers, or solve with GMRES.
 */

#define STIFFNESS_OK 0
#define STIFFNESS_INVALID_ARGUMENT 1
#define STIFFNESS_INVALID_CELL 2
#define STIFFNESS_MLS_FAILED 3
#define STIFFNESS_OUTPUT_CAPACITY_TOO_SMALL 4

typedef struct StiffnessEntry {
    int row_node;
    int col_node;
    double value;
} StiffnessEntry;

int stiffness_assemble_cell(const Node3D *nodes,
                            int node_count,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            double zmin,
                            double zmax,
                            StiffnessEntry *entries,
                            int max_entries,
                            int *entry_count);

#endif
