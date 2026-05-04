#include "lagrange_system.h"

#include <stddef.h>

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
