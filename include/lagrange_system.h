#ifndef LAGRANGE_SYSTEM_H
#define LAGRANGE_SYSTEM_H

/*
 * Dense augmented Lagrange system assembly.
 *
 * This module only arranges the already assembled EFG blocks into
 * [K G^T; G 0]. It does not solve the system and does not call GMRES.
 */

#define LAGRANGE_SYSTEM_OK 0
#define LAGRANGE_SYSTEM_INVALID_ARGUMENT 1

int lagrange_assemble_augmented_system_dense(const double *K,
                                             const double *F,
                                             int node_count,
                                             const double *G,
                                             const double *q,
                                             int constraint_count,
                                             double *A_aug,
                                             double *b_aug);

#endif
