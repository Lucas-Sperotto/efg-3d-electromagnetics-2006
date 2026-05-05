#ifndef MLS_DIAGNOSTIC_H
#define MLS_DIAGNOSTIC_H

#include "nodes.h"

/*
 * Connectivity and moment-matrix diagnostic for MLS shape functions.
 *
 * Implements the diagnostics described in docs/TODO.md sections 3 and 4:
 *   - count active nodes in the cubic support of each evaluation point;
 *   - detect points with fewer than MLS_MIN_SUPPORT_NODES (= 4) active nodes;
 *   - estimate the condition number of A(x) via max/min pivot ratio (see mat4.h);
 *   - aggregate statistics over a set of evaluation points.
 *
 * The minimum cardinality for a linear 3-D basis p^T = [1, x, y, z] is 4.
 * Points with active_nodes < 4 are classified as INVALID.
 * Points with active_nodes < MLS_SUSPECT_SUPPORT_NODES (= 8) are classified
 * as SUSPECT following the guideline in TODO.md section 3.3.
 */

#define MLS_DIAGNOSTIC_OK              0
#define MLS_DIAGNOSTIC_INVALID_ARG     1

/* Minimum nodes required for the 4-term linear 3-D basis. */
#define MLS_MIN_SUPPORT_NODES    4

/* Threshold below which a point is considered numerically suspect. */
#define MLS_SUSPECT_SUPPORT_NODES 8

typedef enum MlsPointStatus {
    MLS_POINT_VALID   = 0, /* active_nodes >= MLS_MIN_SUPPORT_NODES */
    MLS_POINT_SUSPECT = 1, /* active_nodes in [MLS_MIN_SUPPORT_NODES, MLS_SUSPECT_SUPPORT_NODES) */
    MLS_POINT_INVALID = 2  /* active_nodes < MLS_MIN_SUPPORT_NODES */
} MlsPointStatus;

typedef struct MlsPointDiagnostic {
    double x;
    double y;
    double z;
    int    active_nodes;   /* number of nodes whose cubic support covers (x,y,z) */
    double cond_estimate;  /* max/min pivot ratio of A(x); 0.0 if A is singular */
    int    moment_status;  /* MAT4_OK, MAT4_SINGULAR, or MAT4_NEAR_SINGULAR */
    MlsPointStatus status; /* VALID / SUSPECT / INVALID */
} MlsPointDiagnostic;

typedef struct MlsConnectivityStats {
    int    n_total;        /* total points evaluated */
    int    n_invalid;      /* points with active_nodes < MLS_MIN_SUPPORT_NODES */
    int    n_suspect;      /* points with active_nodes in suspect range */
    int    n_moment_fail;  /* points where A(x) is singular or near-singular */
    int    min_nodes;      /* minimum active_nodes seen */
    int    max_nodes;      /* maximum active_nodes seen */
    double mean_nodes;     /* average active_nodes */
    double max_cond;       /* maximum condition estimate seen */
    double mean_cond;      /* mean condition estimate (valid points only) */
} MlsConnectivityStats;

/*
 * Evaluate MLS connectivity and condition diagnostics at a single point.
 *
 * Counts the active nodes in the cubic support of (x,y,z), builds A(x), and
 * estimates its condition number. The result is written to *out.
 *
 * Returns MLS_DIAGNOSTIC_OK on success or MLS_DIAGNOSTIC_INVALID_ARG if any
 * pointer is NULL or node_count <= 0.
 */
int mls_diagnostic_at_point(const Node3D *nodes,
                             int           node_count,
                             double        x,
                             double        y,
                             double        z,
                             MlsPointDiagnostic *out);

/*
 * Aggregate connectivity statistics over an arbitrary set of evaluation points.
 *
 * eval_pts is an array of (3 * n_eval) doubles: {x0,y0,z0, x1,y1,z1, ...}.
 * The result is written to *stats.
 *
 * Returns MLS_DIAGNOSTIC_OK on success.
 */
int mls_connectivity_stats_from_points(const Node3D *nodes,
                                       int           node_count,
                                       const double *eval_pts,
                                       int           n_eval,
                                       MlsConnectivityStats *stats);

/*
 * Aggregate connectivity statistics over a regular nx×ny×nz grid inside the
 * axis-aligned box [xmin,xmax] × [ymin,ymax] × [zmin,zmax].
 *
 * Returns MLS_DIAGNOSTIC_OK on success.
 */
int mls_connectivity_stats_uniform_grid(const Node3D *nodes,
                                        int           node_count,
                                        double xmin, double xmax, int nx,
                                        double ymin, double ymax, int ny,
                                        double zmin, double zmax, int nz,
                                        MlsConnectivityStats *stats);

#endif
