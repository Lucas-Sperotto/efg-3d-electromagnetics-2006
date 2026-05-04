#ifndef MLS_H
#define MLS_H

#include "nodes.h"

/*
 * Didactic Moving Least Squares helpers.
 *
 * This first MLS step computes only the scalar shape-function values Phi_I(x)
 * for a linear 3-D basis p^T = [1, x, y, z]. It does not compute derivatives,
 * assembly terms, Lagrange multipliers, or GMRES systems.
 */

#define MLS_OK 0
#define MLS_INVALID_ARGUMENT 1
#define MLS_INVALID_SUPPORT_RADIUS 2
#define MLS_NOT_ENOUGH_SUPPORT_NODES 3
#define MLS_OUTPUT_CAPACITY_TOO_SMALL 4
#define MLS_SINGULAR_MOMENT_MATRIX 5
#define MLS_NEAR_SINGULAR_MOMENT_MATRIX 6

typedef struct MlsShapeValue {
    int node_index;
    double phi;
} MlsShapeValue;

int mls_linear3d_shape_functions(const Node3D *nodes,
                                 int node_count,
                                 double x,
                                 double y,
                                 double z,
                                 MlsShapeValue *values,
                                 int max_values,
                                 int *value_count);

#endif
