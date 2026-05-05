#include "lagrange_system.h"

#include <stddef.h>
#include <stdlib.h>

static void zero_augmented_system(double *A_aug, double *b_aug, int total_size)
{
    for (int i = 0; i < total_size * total_size; ++i) {
        A_aug[i] = 0.0;
    }

    for (int i = 0; i < total_size; ++i) {
        b_aug[i] = 0.0;
    }
}

/*
 * Assemble the dense augmented system used to impose Dirichlet conditions by
 * Lagrange multipliers in EFG.
 *
 * The block layout follows docs/05_mapa_equacoes_codigo.md and
 * docs/06_condicoes_de_contorno_e_lagrange.md:
 *
 *     [ K  G^T ] [ V      ] = [ F ]
 *     [ G   0  ] [ lambda ]   [ q ]
 *
 * This function performs only the dense block copy. It does not apply a
 * solver, does not form GMRES data structures, and does not modify K, F, G,
 * or q.
 */
int lagrange_assemble_augmented_system_dense(const double *K,
                                             const double *F,
                                             int node_count,
                                             const double *G,
                                             const double *q,
                                             int constraint_count,
                                             double *A_aug,
                                             double *b_aug)
{
    int total_size;

    if (K == NULL || F == NULL || G == NULL || q == NULL ||
        A_aug == NULL || b_aug == NULL ||
        node_count <= 0 || constraint_count <= 0) {
        return LAGRANGE_SYSTEM_INVALID_ARGUMENT;
    }

    total_size = node_count + constraint_count;
    zero_augmented_system(A_aug, b_aug, total_size);

    for (int row = 0; row < node_count; ++row) {
        b_aug[row] = F[row];

        for (int col = 0; col < node_count; ++col) {
            A_aug[(row * total_size) + col] = K[(row * node_count) + col];
        }
    }

    for (int constraint = 0; constraint < constraint_count; ++constraint) {
        const int aug_row = node_count + constraint;

        b_aug[aug_row] = q[constraint];

        for (int node = 0; node < node_count; ++node) {
            const double value = G[(constraint * node_count) + node];

            A_aug[(node * total_size) + aug_row] = value;
            A_aug[(aug_row * total_size) + node] = value;
        }
    }

    return LAGRANGE_SYSTEM_OK;
}

/*
 * Sparse COO variant of the augmented Lagrange system assembler.
 *
 * The K block is copied from K_coo (indices are unchanged — it occupies
 * rows/cols 0..node_count-1 in the augmented system). The G and G^T blocks
 * are inserted from the dense G matrix with column offset node_count.
 * Duplicate (row,col) entries are accumulated and summed at CSR conversion.
 */
int lagrange_assemble_augmented_sparse_coo(const SparseCOO *K_coo,
                                           const double *F,
                                           int node_count,
                                           const double *G,
                                           const double *q,
                                           int constraint_count,
                                           SparseCOO *A_aug_coo,
                                           double *b_aug)
{
    int e;
    int c;
    int n;

    if (K_coo == NULL || F == NULL || G == NULL || q == NULL ||
        A_aug_coo == NULL || b_aug == NULL ||
        node_count <= 0 || constraint_count <= 0) {
        return LAGRANGE_SYSTEM_INVALID_ARGUMENT;
    }

    /* b_aug = [F; q] */
    for (int i = 0; i < node_count; ++i) {
        b_aug[i] = F[i];
    }
    for (c = 0; c < constraint_count; ++c) {
        b_aug[node_count + c] = q[c];
    }

    /* K block: copy COO entries verbatim (already in rows/cols 0..node_count-1) */
    for (e = 0; e < K_coo->count; ++e) {
        const int rc = sparse_coo_add(A_aug_coo,
                                      K_coo->row[e],
                                      K_coo->col[e],
                                      K_coo->val[e]);
        if (rc != SPARSE_OK) {
            return LAGRANGE_SYSTEM_ALLOCATION_FAILED;
        }
    }

    /* G and G^T blocks */
    for (c = 0; c < constraint_count; ++c) {
        for (n = 0; n < node_count; ++n) {
            const double val = G[c * node_count + n];
            int rc;

            /* G^T: A_aug[n, node_count+c] */
            rc = sparse_coo_add(A_aug_coo, n, node_count + c, val);
            if (rc != SPARSE_OK) {
                return LAGRANGE_SYSTEM_ALLOCATION_FAILED;
            }

            /* G: A_aug[node_count+c, n] */
            rc = sparse_coo_add(A_aug_coo, node_count + c, n, val);
            if (rc != SPARSE_OK) {
                return LAGRANGE_SYSTEM_ALLOCATION_FAILED;
            }
        }
    }

    return LAGRANGE_SYSTEM_OK;
}
