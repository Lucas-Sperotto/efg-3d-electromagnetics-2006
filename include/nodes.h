#ifndef NODES_H
#define NODES_H

/*
 * Basic 3-D node used by the didactic EFG implementation.
 *
 * The support_radius field stores d_mI, the radius/half-size parameter used by
 * the weight-function support. The article defines it from local node density:
 * three nearest neighbors and a factor of 1.5
 * (docs/03_resultados_numericos.md, section 3.2).
 */
typedef struct Node3D {
    double x;
    double y;
    double z;
    double support_radius;
} Node3D;

Node3D node3d_make(double x, double y, double z);

double node3d_distance(const Node3D *a, const Node3D *b);

#endif
