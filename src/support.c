#include "support.h"

#include <stddef.h>
#include <stdlib.h>

/*
 * Find the k-th nearest neighbor distance using a plain O(n^2)-style scan.
 *
 * The article uses the three nearest neighbors and sets d_mI to 1.5 times the
 * distance to the third neighbor (docs/03_resultados_numericos.md, section
 * 3.2). This implementation keeps the search intentionally simple and
 * didactic; no k-d tree or octree is used in this step.
 */
int support_kth_neighbor_distance(const Node3D *nodes,
                                  int count,
                                  int node_index,
                                  int k,
                                  double *distance)
{
    double *nearest_distances = NULL;
    int status = SUPPORT_OK;

    if (nodes == NULL || distance == NULL || count <= 1 || node_index < 0 ||
        node_index >= count || k <= 0) {
        return SUPPORT_INVALID_ARGUMENT;
    }

    if (count - 1 < k) {
        return SUPPORT_NOT_ENOUGH_NEIGHBORS;
    }

    nearest_distances = malloc((size_t)k * sizeof(nearest_distances[0]));
    if (nearest_distances == NULL) {
        return SUPPORT_INVALID_ARGUMENT;
    }

    for (int i = 0; i < k; ++i) {
        nearest_distances[i] = -1.0;
    }

    for (int i = 0; i < count; ++i) {
        double candidate_distance;
        int insert_at = -1;

        if (i == node_index) {
            continue;
        }

        candidate_distance = node3d_distance(&nodes[node_index], &nodes[i]);

        for (int slot = 0; slot < k; ++slot) {
            if (nearest_distances[slot] < 0.0 || candidate_distance < nearest_distances[slot]) {
                insert_at = slot;
                break;
            }
        }

        if (insert_at >= 0) {
            for (int slot = k - 1; slot > insert_at; --slot) {
                nearest_distances[slot] = nearest_distances[slot - 1];
            }
            nearest_distances[insert_at] = candidate_distance;
        }
    }

    if (nearest_distances[k - 1] < 0.0) {
        status = SUPPORT_NOT_ENOUGH_NEIGHBORS;
    } else {
        *distance = nearest_distances[k - 1];
    }

    free(nearest_distances);
    return status;
}

int support_assign_from_kth_neighbor(Node3D *nodes, int count, int k, double scale)
{
    if (nodes == NULL || count <= 1 || k <= 0 || scale <= 0.0) {
        return SUPPORT_INVALID_ARGUMENT;
    }

    for (int i = 0; i < count; ++i) {
        double kth_distance = 0.0;
        const int status = support_kth_neighbor_distance(nodes, count, i, k, &kth_distance);

        if (status != SUPPORT_OK) {
            return status;
        }

        nodes[i].support_radius = scale * kth_distance;
    }

    return SUPPORT_OK;
}

int support_assign_article_default(Node3D *nodes, int count)
{
    /*
     * Article default from docs/03_resultados_numericos.md, section 3.2:
     * k = 3 nearest neighbors and d_mI = 1.5 times the third-neighbor distance.
     */
    return support_assign_from_kth_neighbor(nodes, count, 3, 1.5);
}
