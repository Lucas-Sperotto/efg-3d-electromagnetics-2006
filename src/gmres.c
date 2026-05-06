#include "efg/gmres.h"
#include "sparse_matrix.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int efg_gmres_module_ready(void)
{
    return 1;
}

/* ---------------------------------------------------------------- helpers */

static double dot_n(const double *a, const double *b, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        s += a[i] * b[i];
    }
    return s;
}

static double norm2_n(const double *a, int n)
{
    return sqrt(dot_n(a, a, n));
}

/* ---------------------------------------------------------------- gmres_solve */

/*
 * Restarted GMRES(restart).
 *
 * Memory layout:
 *   V  — Krylov basis: row j is the j-th basis vector, V[j*n + k]
 *   H  — upper Hessenberg: H[i*restart + j] = H_{ij}  (i rows, j cols)
 *   g  — current RHS of the projected LS problem after Givens rotations
 *   cs, sn — Givens rotation cosines and sines
 */
int gmres_solve(const SparseCSR *A,
                const double    *b,
                double          *x,
                double           tol,
                int              max_iter,
                int              restart,
                const double    *diag_precond,
                GmresResult     *result)
{
    const double eps = 1.0e-30;
    int n;
    int total_iter;
    int converged;
    double beta0;
    double beta;

    double *V  = NULL;
    double *H  = NULL;
    double *g  = NULL;
    double *cs = NULL;
    double *sn = NULL;
    double *w  = NULL;
    double *r  = NULL;
    double *y  = NULL;

    if (A == NULL || b == NULL || x == NULL || result == NULL ||
        tol <= 0.0 || max_iter <= 0 || restart <= 0) {
        return GMRES_INVALID_ARGUMENT;
    }

    n = A->nrows;

    V  = malloc((size_t)(restart + 1) * (size_t)n * sizeof(double));
    H  = calloc((size_t)(restart + 1) * (size_t)restart, sizeof(double));
    g  = calloc((size_t)(restart + 1), sizeof(double));
    cs = malloc((size_t)restart * sizeof(double));
    sn = malloc((size_t)restart * sizeof(double));
    w  = malloc((size_t)n * sizeof(double));
    r  = malloc((size_t)n * sizeof(double));
    y  = malloc((size_t)restart * sizeof(double));

    if (!V || !H || !g || !cs || !sn || !w || !r || !y) {
        free(V); free(H); free(g);
        free(cs); free(sn); free(w); free(r); free(y);
        return GMRES_ALLOCATION_FAILED;
    }

    /* Initial residual r = b - A*x (true residual, for reporting) */
    sparse_csr_matvec(A, x, r);
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - r[i];
    }
    result->residual_init = norm2_n(r, n);

    /* Apply left preconditioner M⁻¹ = diag(diag_precond) if supplied */
    if (diag_precond) {
        for (int i = 0; i < n; ++i) r[i] *= diag_precond[i];
    }
    beta0 = norm2_n(r, n);

    total_iter = 0;
    converged  = 0;

    if (beta0 < eps) {
        converged = 1;
        goto done;
    }

    while (total_iter < max_iter && !converged) {

        /* Restart: recompute residual and initialise Krylov space */
        sparse_csr_matvec(A, x, r);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - r[i];
        }
        if (diag_precond) {
            for (int i = 0; i < n; ++i) r[i] *= diag_precond[i];
        }
        beta = norm2_n(r, n);

        if (beta / beta0 < tol) {
            converged = 1;
            break;
        }

        /* V[0] = r / beta */
        for (int k = 0; k < n; ++k) {
            V[k] = r[k] / beta;
        }

        /* g = beta * e_1 */
        memset(g, 0, (size_t)(restart + 1) * sizeof(double));
        g[0] = beta;

        /* Clear Hessenberg for this restart */
        memset(H, 0, (size_t)(restart + 1) * (size_t)restart * sizeof(double));

        int j_end = 0;

        for (int j = 0; j < restart && total_iter < max_iter; ++j) {

            /* w = A * V[j], then apply left preconditioner */
            sparse_csr_matvec(A, V + j * n, w);
            if (diag_precond) {
                for (int k = 0; k < n; ++k) w[k] *= diag_precond[k];
            }

            /* Modified Gram-Schmidt orthogonalisation */
            for (int i = 0; i <= j; ++i) {
                double hij = dot_n(w, V + i * n, n);
                H[i * restart + j] = hij;
                for (int k = 0; k < n; ++k) {
                    w[k] -= hij * V[i * n + k];
                }
            }

            {
                double h_next = norm2_n(w, n);
                H[(j + 1) * restart + j] = h_next;

                if (h_next > eps) {
                    for (int k = 0; k < n; ++k) {
                        V[(j + 1) * n + k] = w[k] / h_next;
                    }
                }
            }

            /* Apply previous Givens rotations to column j of H */
            for (int i = 0; i < j; ++i) {
                double tmp =  cs[i] * H[i       * restart + j]
                            + sn[i] * H[(i + 1) * restart + j];
                H[(i + 1) * restart + j] = -sn[i] * H[i       * restart + j]
                                           +  cs[i] * H[(i + 1) * restart + j];
                H[i * restart + j] = tmp;
            }

            /* Compute new Givens rotation to zero H[j+1][j] */
            {
                double hj  = H[j       * restart + j];
                double hj1 = H[(j + 1) * restart + j];
                double rho = sqrt(hj * hj + hj1 * hj1);

                if (rho > eps) {
                    cs[j] = hj  / rho;
                    sn[j] = hj1 / rho;
                } else {
                    cs[j] = 1.0;
                    sn[j] = 0.0;
                }

                H[j       * restart + j] = rho;
                H[(j + 1) * restart + j] = 0.0;
            }

            /* Apply rotation to g */
            {
                double gj = g[j];
                g[j]     =  cs[j] * gj;
                g[j + 1] = -sn[j] * gj;
            }

            ++total_iter;
            j_end = j + 1;

            if (fabs(g[j + 1]) / beta0 < tol) {
                converged = 1;
                break;
            }
        }

        /* Back-substitute: solve upper triangular H[0..j_end-1] * y = g */
        for (int k = j_end - 1; k >= 0; --k) {
            double s = g[k];
            for (int l = k + 1; l < j_end; ++l) {
                s -= H[k * restart + l] * y[l];
            }
            {
                double diag = H[k * restart + k];
                y[k] = (fabs(diag) > eps) ? s / diag : 0.0;
            }
        }

        /* Update x = x + V * y */
        for (int i = 0; i < j_end; ++i) {
            for (int k = 0; k < n; ++k) {
                x[k] += y[i] * V[i * n + k];
            }
        }
    }

done:
    /* Compute final residual */
    sparse_csr_matvec(A, x, r);
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - r[i];
    }
    result->residual_final = norm2_n(r, n);
    result->converged      = converged;
    result->iterations     = total_iter;

    free(V); free(H); free(g);
    free(cs); free(sn); free(w); free(r); free(y);

    return GMRES_OK;
}
