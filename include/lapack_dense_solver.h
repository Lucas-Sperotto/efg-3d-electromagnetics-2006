#ifndef LAPACK_DENSE_SOLVER_H
#define LAPACK_DENSE_SOLVER_H

#include "sparse_matrix.h"

/*
 * Optional LAPACK dense validation solver for the augmented Lagrange system.
 *
 * The core didactic solver remains GMRES. This module is a validation path:
 * it converts the already validated CSR matrix to LAPACK's dense column-major
 * storage and calls dgesv when LAPACK is available at CMake configuration time.
 */

#define LAPACK_DENSE_SOLVER_OK                 0
#define LAPACK_DENSE_SOLVER_INVALID_ARGUMENT   1
#define LAPACK_DENSE_SOLVER_ALLOCATION_FAILED  2
#define LAPACK_DENSE_SOLVER_UNAVAILABLE        3
#define LAPACK_DENSE_SOLVER_SINGULAR           4
#define LAPACK_DENSE_SOLVER_LAPACK_ERROR       5

typedef struct LapackDenseSolveResult {
    int available;              /* 1 when compiled with LAPACK, 0 otherwise */
    int lapack_info;            /* dgesv INFO value; 0 means success */
    double conversion_time_s;   /* CSR -> dense column-major time */
    double solve_time_s;        /* dgesv factorization/solve time */
    double residual_final;      /* ||b - A*x|| after the LAPACK solve */
} LapackDenseSolveResult;

int sparse_csr_to_dense_column_major(const SparseCSR *csr,
                                     double *dense,
                                     int leading_dim);

/*
 * Solve the augmented Lagrange system A*x = b with LAPACK dgesv.
 *
 * Article equation: Eq. (7), mapped in docs/05_mapa_equacoes_codigo.md.
 * Formula reference: docs/02_formulacao_matematica.md, section 2.1, and
 * docs/06_condicoes_de_contorno_e_lagrange.md.
 * Mathematical meaning: direct dense validation of the saddle-point system
 * [K G^T; G 0] that imposes Dirichlet conditions through Lagrange
 * multipliers.
 */
int solve_augmented_lapack_dense(const SparseCSR *A,
                                 const double *b,
                                 double *x,
                                 LapackDenseSolveResult *result);

#endif
