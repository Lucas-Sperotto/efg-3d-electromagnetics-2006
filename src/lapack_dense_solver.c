#include "lapack_dense_solver.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef EFG_HAVE_LAPACK
extern void dgesv_(const int *n,
                   const int *nrhs,
                   double *a,
                   const int *lda,
                   int *ipiv,
                   double *b,
                   const int *ldb,
                   int *info);
#endif

static double elapsed_s(clock_t start)
{
    return (double)(clock() - start) / (double)CLOCKS_PER_SEC;
}

static void lapack_result_init(LapackDenseSolveResult *result)
{
    if (result == NULL) {
        return;
    }

#ifdef EFG_HAVE_LAPACK
    result->available = 1;
#else
    result->available = 0;
#endif
    result->lapack_info = 0;
    result->conversion_time_s = NAN;
    result->solve_time_s = NAN;
    result->residual_final = NAN;
}

static int checked_square_count(int n, size_t *matrix_count)
{
    const size_t n_size = (size_t)n;
    const size_t max_size = (size_t)-1;

    if (n <= 0 || matrix_count == NULL) {
        return 0;
    }
    if (n_size > max_size / n_size) {
        return 0;
    }

    *matrix_count = n_size * n_size;
    if (*matrix_count > max_size / sizeof(double)) {
        return 0;
    }

    return 1;
}

static double residual_norm_csr(const SparseCSR *A,
                                const double *b,
                                const double *x)
{
    double *Ax = NULL;
    double norm2 = 0.0;

    if (A == NULL || b == NULL || x == NULL || A->nrows <= 0) {
        return NAN;
    }

    Ax = malloc((size_t)A->nrows * sizeof(Ax[0]));
    if (Ax == NULL) {
        return NAN;
    }

    if (sparse_csr_matvec(A, x, Ax) != SPARSE_OK) {
        free(Ax);
        return NAN;
    }

    for (int i = 0; i < A->nrows; ++i) {
        const double r = b[i] - Ax[i];
        norm2 += r * r;
    }

    free(Ax);
    return sqrt(norm2);
}

int sparse_csr_to_dense_column_major(const SparseCSR *csr,
                                     double *dense,
                                     int leading_dim)
{
    if (csr == NULL || dense == NULL ||
        csr->row_ptr == NULL || csr->col_ind == NULL || csr->values == NULL ||
        csr->nrows <= 0 || csr->ncols <= 0 || leading_dim < csr->nrows) {
        return LAPACK_DENSE_SOLVER_INVALID_ARGUMENT;
    }

    memset(dense, 0,
           (size_t)leading_dim * (size_t)csr->ncols * sizeof(dense[0]));

    for (int row = 0; row < csr->nrows; ++row) {
        if (csr->row_ptr[row] > csr->row_ptr[row + 1]) {
            return LAPACK_DENSE_SOLVER_INVALID_ARGUMENT;
        }

        for (int p = csr->row_ptr[row]; p < csr->row_ptr[row + 1]; ++p) {
            const int col = csr->col_ind[p];
            if (col < 0 || col >= csr->ncols) {
                return LAPACK_DENSE_SOLVER_INVALID_ARGUMENT;
            }
            dense[row + col * leading_dim] += csr->values[p];
        }
    }

    return LAPACK_DENSE_SOLVER_OK;
}

int solve_augmented_lapack_dense(const SparseCSR *A,
                                 const double *b,
                                 double *x,
                                 LapackDenseSolveResult *result)
{
    lapack_result_init(result);

    if (A == NULL || b == NULL || x == NULL ||
        A->nrows <= 0 || A->ncols <= 0 || A->nrows != A->ncols) {
        return LAPACK_DENSE_SOLVER_INVALID_ARGUMENT;
    }

#ifndef EFG_HAVE_LAPACK
    return LAPACK_DENSE_SOLVER_UNAVAILABLE;
#else
    {
        const int n = A->nrows;
        const int nrhs = 1;
        const int lda = n;
        const int ldb = n;
        int info = 0;
        size_t matrix_count = 0;
        double *A_dense = NULL;
        double *rhs = NULL;
        int *ipiv = NULL;
        clock_t t_start;
        int status;

        if (!checked_square_count(n, &matrix_count)) {
            return LAPACK_DENSE_SOLVER_ALLOCATION_FAILED;
        }

        A_dense = calloc(matrix_count, sizeof(A_dense[0]));
        rhs = malloc((size_t)n * sizeof(rhs[0]));
        ipiv = malloc((size_t)n * sizeof(ipiv[0]));
        if (A_dense == NULL || rhs == NULL || ipiv == NULL) {
            free(A_dense);
            free(rhs);
            free(ipiv);
            return LAPACK_DENSE_SOLVER_ALLOCATION_FAILED;
        }

        t_start = clock();
        status = sparse_csr_to_dense_column_major(A, A_dense, lda);
        if (result != NULL) {
            result->conversion_time_s = elapsed_s(t_start);
        }
        if (status != LAPACK_DENSE_SOLVER_OK) {
            free(A_dense);
            free(rhs);
            free(ipiv);
            return status;
        }

        memcpy(rhs, b, (size_t)n * sizeof(rhs[0]));

        t_start = clock();
        dgesv_(&n, &nrhs, A_dense, &lda, ipiv, rhs, &ldb, &info);
        if (result != NULL) {
            result->solve_time_s = elapsed_s(t_start);
            result->lapack_info = info;
        }

        if (info != 0) {
            free(A_dense);
            free(rhs);
            free(ipiv);
            return (info > 0)
                   ? LAPACK_DENSE_SOLVER_SINGULAR
                   : LAPACK_DENSE_SOLVER_LAPACK_ERROR;
        }

        memcpy(x, rhs, (size_t)n * sizeof(x[0]));

        if (result != NULL) {
            result->residual_final = residual_norm_csr(A, b, x);
        }

        free(A_dense);
        free(rhs);
        free(ipiv);
        return LAPACK_DENSE_SOLVER_OK;
    }
#endif
}
