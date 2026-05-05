#ifndef LAGRANGE_SYSTEM_H
#define LAGRANGE_SYSTEM_H

#include "sparse_matrix.h"

/*
 * Augmented Lagrange system assembly — dense and sparse variants.
 *
 * Both functions arrange EFG blocks into [K G^T; G 0]. They do not solve
 * the system and do not call GMRES.
 */

#define LAGRANGE_SYSTEM_OK                 0
#define LAGRANGE_SYSTEM_INVALID_ARGUMENT   1
#define LAGRANGE_SYSTEM_ALLOCATION_FAILED  2

int lagrange_assemble_augmented_system_dense(const double *K,
                                             const double *F,
                                             int node_count,
                                             const double *G,
                                             const double *q,
                                             int constraint_count,
                                             double *A_aug,
                                             double *b_aug);

/*
 * Sparse COO variant of the augmented system assembler.
 *
 * K_coo must already contain the assembled stiffness matrix (from
 * global_stiffness_assemble_sparse_coo). A_aug_coo must have been
 * initialised with sparse_coo_create() with dimensions
 * (node_count + constraint_count) x (node_count + constraint_count).
 * b_aug must point to an array of length (node_count + constraint_count).
 *
 * Duplicate (row,col) entries are allowed in A_aug_coo; they will be summed
 * when the COO is converted to CSR.
 */
int lagrange_assemble_augmented_sparse_coo(const SparseCOO *K_coo,
                                           const double *F,
                                           int node_count,
                                           const double *G,
                                           const double *q,
                                           int constraint_count,
                                           SparseCOO *A_aug_coo,
                                           double *b_aug);

#endif
