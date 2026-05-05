#ifndef EFG_GMRES_H
#define EFG_GMRES_H

#include "sparse_matrix.h"

/*
 * Restarted GMRES solver for the augmented Lagrange system [K G^T; G 0].
 *
 * The system is symmetric indefinite; GMRES is used because the standard
 * CG variants require positive definiteness.
 */

#define GMRES_OK                 0
#define GMRES_INVALID_ARGUMENT   1
#define GMRES_ALLOCATION_FAILED  2

typedef struct GmresResult {
    int    converged;       /* 1 if converged within tolerance, 0 otherwise */
    int    iterations;      /* total Arnoldi steps performed */
    double residual_init;   /* ||b - A*x0|| at start */
    double residual_final;  /* ||b - A*x||  at exit */
} GmresResult;

/*
 * Solve A*x = b using restarted GMRES(restart).
 *
 * x   — initial guess on input, approximate solution on output.
 * tol — relative tolerance: convergence when ||r_k|| / ||r_0|| < tol.
 * max_iter — maximum total Arnoldi steps (across all restarts).
 * restart  — Krylov subspace dimension per restart cycle (e.g. 30).
 *
 * Returns GMRES_OK on both convergence and non-convergence; inspect
 * result->converged to distinguish the two cases.
 */
int gmres_solve(const SparseCSR *A,
                const double    *b,
                double          *x,
                double           tol,
                int              max_iter,
                int              restart,
                GmresResult     *result);

/* Link-time smoke symbol — kept for backward compatibility. */
int efg_gmres_module_ready(void);

#endif
