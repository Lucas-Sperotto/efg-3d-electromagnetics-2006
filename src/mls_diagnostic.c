#include "mls_diagnostic.h"

#include "mat4.h"
#include "weight.h"

#include <math.h>
#include <stddef.h>

/* Cubic-box containment test identical to the one in mls.c. */
static int node_in_cubic_support(const Node3D *node, double x, double y, double z)
{
    return fabs(x - node->x) <= node->support_radius &&
           fabs(y - node->y) <= node->support_radius &&
           fabs(z - node->z) <= node->support_radius;
}

static void linear_basis_at(double x, double y, double z, double p[4])
{
    p[0] = 1.0;
    p[1] = x;
    p[2] = y;
    p[3] = z;
}

/*
 * Build moment matrix A(x) = Σ_I w_I p(x_I) p^T(x_I) at evaluation point
 * (ex, ey, ez) using the cubic tensor-product weight. Only nodes whose cubic
 * support covers (ex,ey,ez) contribute.
 *
 * Returns the number of active nodes.
 */
static int build_moment_matrix(const Node3D *nodes,
                                int           node_count,
                                double        ex,
                                double        ey,
                                double        ez,
                                double        moment[4][4])
{
    int active = 0;

    mat4_zero(moment);

    for (int i = 0; i < node_count; ++i) {
        double p_node[4];
        double w = 0.0;

        if (!node_in_cubic_support(&nodes[i], ex, ey, ez)) {
            continue;
        }

        ++active;

        if (weight_cubic_tensor3d(ex, ey, ez,
                                  nodes[i].x, nodes[i].y, nodes[i].z,
                                  nodes[i].support_radius, &w) != WEIGHT_OK) {
            continue;
        }

        linear_basis_at(nodes[i].x, nodes[i].y, nodes[i].z, p_node);

        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                moment[r][c] += w * p_node[r] * p_node[c];
            }
        }
    }

    return active;
}

int mls_diagnostic_at_point(const Node3D       *nodes,
                              int                 node_count,
                              double              x,
                              double              y,
                              double              z,
                              MlsPointDiagnostic *out)
{
    double moment[4][4];
    int active;

    if (nodes == NULL || out == NULL || node_count <= 0) {
        return MLS_DIAGNOSTIC_INVALID_ARG;
    }

    out->x = x;
    out->y = y;
    out->z = z;
    out->cond_estimate = 0.0;
    out->moment_status = MAT4_OK;

    active = build_moment_matrix(nodes, node_count, x, y, z, moment);
    out->active_nodes = active;

    if (active < MLS_MIN_SUPPORT_NODES) {
        out->status = MLS_POINT_INVALID;
        out->moment_status = MAT4_SINGULAR; /* cannot form a valid A */
        return MLS_DIAGNOSTIC_OK;
    }

    /* Estimate condition number from max/min pivot ratio. */
    out->moment_status = mat4_cond_estimate((const double (*)[4])moment, &out->cond_estimate);

    if (active < MLS_SUSPECT_SUPPORT_NODES) {
        out->status = MLS_POINT_SUSPECT;
    } else {
        out->status = MLS_POINT_VALID;
    }

    return MLS_DIAGNOSTIC_OK;
}

/* ------------------------------------------------------------------ */

static void stats_init(MlsConnectivityStats *s)
{
    s->n_total       = 0;
    s->n_invalid     = 0;
    s->n_suspect     = 0;
    s->n_moment_fail = 0;
    s->min_nodes     = -1;
    s->max_nodes     = 0;
    s->mean_nodes    = 0.0;
    s->max_cond      = 0.0;
    s->mean_cond     = 0.0;
}

static void stats_accumulate(MlsConnectivityStats *s, const MlsPointDiagnostic *d)
{
    ++s->n_total;

    if (d->status == MLS_POINT_INVALID) {
        ++s->n_invalid;
    } else if (d->status == MLS_POINT_SUSPECT) {
        ++s->n_suspect;
    }

    if (d->moment_status != MAT4_OK) {
        ++s->n_moment_fail;
    }

    if (s->min_nodes < 0 || d->active_nodes < s->min_nodes) {
        s->min_nodes = d->active_nodes;
    }
    if (d->active_nodes > s->max_nodes) {
        s->max_nodes = d->active_nodes;
    }

    /* Welford running mean for node count. */
    s->mean_nodes += (d->active_nodes - s->mean_nodes) / (double)s->n_total;

    if (d->cond_estimate > s->max_cond) {
        s->max_cond = d->cond_estimate;
    }
    /* Running mean for condition (all points, including suspect). */
    s->mean_cond += (d->cond_estimate - s->mean_cond) / (double)s->n_total;
}

int mls_connectivity_stats_from_points(const Node3D         *nodes,
                                        int                   node_count,
                                        const double         *eval_pts,
                                        int                   n_eval,
                                        MlsConnectivityStats *stats)
{
    if (nodes == NULL || eval_pts == NULL || stats == NULL ||
        node_count <= 0 || n_eval <= 0) {
        return MLS_DIAGNOSTIC_INVALID_ARG;
    }

    stats_init(stats);

    for (int i = 0; i < n_eval; ++i) {
        MlsPointDiagnostic d;
        double px = eval_pts[3 * i + 0];
        double py = eval_pts[3 * i + 1];
        double pz = eval_pts[3 * i + 2];

        mls_diagnostic_at_point(nodes, node_count, px, py, pz, &d);
        stats_accumulate(stats, &d);
    }

    return MLS_DIAGNOSTIC_OK;
}

int mls_connectivity_stats_uniform_grid(const Node3D         *nodes,
                                         int                   node_count,
                                         double xmin, double xmax, int nx,
                                         double ymin, double ymax, int ny,
                                         double zmin, double zmax, int nz,
                                         MlsConnectivityStats *stats)
{
    if (nodes == NULL || stats == NULL || node_count <= 0 ||
        nx <= 0 || ny <= 0 || nz <= 0 ||
        xmin >= xmax || ymin >= ymax || zmin >= zmax) {
        return MLS_DIAGNOSTIC_INVALID_ARG;
    }

    stats_init(stats);

    for (int ix = 0; ix < nx; ++ix) {
        double px = xmin + (ix + 0.5) * (xmax - xmin) / (double)nx;
        for (int iy = 0; iy < ny; ++iy) {
            double py = ymin + (iy + 0.5) * (ymax - ymin) / (double)ny;
            for (int iz = 0; iz < nz; ++iz) {
                double pz = zmin + (iz + 0.5) * (zmax - zmin) / (double)nz;
                MlsPointDiagnostic d;

                mls_diagnostic_at_point(nodes, node_count, px, py, pz, &d);
                stats_accumulate(stats, &d);
            }
        }
    }

    return MLS_DIAGNOSTIC_OK;
}
