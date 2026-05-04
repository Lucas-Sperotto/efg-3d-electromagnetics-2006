#ifndef DIRICHLET_H
#define DIRICHLET_H

#include "nodes.h"

/*
 * Dense Dirichlet constraint assembly for EFG.
 *
 * This module builds the G matrix and q vector used by Lagrange multiplier
 * constraints. It does not build the augmented system [K G^T; G 0], does not
 * apply boundary conditions directly, and does not call GMRES.
 */

#define DIRICHLET_OK 0
#define DIRICHLET_INVALID_ARGUMENT 1
#define DIRICHLET_INVALID_SUPPORT_RADIUS 2
#define DIRICHLET_ALLOCATION_FAILED 3
#define DIRICHLET_MLS_FAILED 4

typedef struct DirichletPoint {
    double x;
    double y;
    double z;
    double value;
} DirichletPoint;

int dirichlet_assemble_constraints_dense(const Node3D *nodes,
                                         int node_count,
                                         const DirichletPoint *points,
                                         int point_count,
                                         double *G,
                                         double *q);

#endif
