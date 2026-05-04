#ifndef SUPPORT_H
#define SUPPORT_H

#include "nodes.h"

/*
 * Didactic support-radius helpers.
 *
 * This step intentionally uses a simple O(n^2) search over all nodes. The
 * article mentions efficient point-location structures elsewhere, but here we
 * avoid k-d trees, octrees, and optimizations until the baseline is tested.
 */

#define SUPPORT_OK 0
#define SUPPORT_INVALID_ARGUMENT 1
#define SUPPORT_NOT_ENOUGH_NEIGHBORS 2

int support_kth_neighbor_distance(const Node3D *nodes,
                                  int count,
                                  int node_index,
                                  int k,
                                  double *distance);

int support_assign_from_kth_neighbor(Node3D *nodes, int count, int k, double scale);

int support_assign_article_default(Node3D *nodes, int count);

#endif
