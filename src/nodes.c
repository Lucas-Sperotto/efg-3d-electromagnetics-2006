#include "nodes.h"

#include <math.h>

Node3D node3d_make(double x, double y, double z)
{
    Node3D node;

    node.x = x;
    node.y = y;
    node.z = z;
    node.support_radius = 0.0;

    return node;
}

double node3d_distance(const Node3D *a, const Node3D *b)
{
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;

    return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

/*
 * Link-time smoke symbol for nodal data.
 */
int efg_nodes_module_ready(void)
{
    return 1;
}
