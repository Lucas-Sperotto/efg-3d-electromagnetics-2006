#ifndef DENSE_SOLVER_H
#define DENSE_SOLVER_H

/*
 * Didactic dense linear solver for small validation systems.
 *
 * The article ultimately motivates iterative solvers such as GMRES for larger
 * systems. This module is intentionally simpler: it solves dense A x = b
 * systems so early EFG assemblies can be checked on small examples first.
 */

#define DENSE_SOLVER_OK 0
#define DENSE_SOLVER_INVALID_ARGUMENT 1
#define DENSE_SOLVER_SINGULAR 2
#define DENSE_SOLVER_NEAR_SINGULAR 3
#define DENSE_SOLVER_ALLOCATION_FAILED 4

int dense_solve(const double *A, const double *b, int n, double *x);

#endif
