#include "dense_solver.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define DENSE_SOLVER_PIVOT_TOLERANCE 1.0e-14

static int checked_sizes(int n, size_t *matrix_count, size_t *vector_count)
{
    const size_t n_size = (size_t)n;
    const size_t max_size = (size_t)-1;

    if (n_size > max_size / n_size) {
        return 0;
    }

    *matrix_count = n_size * n_size;
    *vector_count = n_size;

    if (*matrix_count > max_size / sizeof(double)) {
        return 0;
    }

    if (*vector_count > max_size / sizeof(double)) {
        return 0;
    }

    return 1;
}

static void swap_rows(double *A, double *b, int n, int row_a, int row_b)
{
    if (row_a == row_b) {
        return;
    }

    for (int col = 0; col < n; ++col) {
        const double temporary = A[(row_a * n) + col];

        A[(row_a * n) + col] = A[(row_b * n) + col];
        A[(row_b * n) + col] = temporary;
    }

    {
        const double temporary = b[row_a];

        b[row_a] = b[row_b];
        b[row_b] = temporary;
    }
}

/*
 * Solve A x = b by Gaussian elimination with partial pivoting.
 *
 * A is stored in row-major order. The routine works on internal dense copies,
 * leaving the caller's A and b unchanged. This is a didactic direct solver for
 * small systems before introducing the GMRES implementation mentioned in the
 * article.
 */
int dense_solve(const double *A, const double *b, int n, double *x)
{
    size_t matrix_count = 0;
    size_t vector_count = 0;
    double *work_A = NULL;
    double *work_b = NULL;
    double *work_x = NULL;

    if (A == NULL || b == NULL || x == NULL || n <= 0) {
        return DENSE_SOLVER_INVALID_ARGUMENT;
    }

    if (!checked_sizes(n, &matrix_count, &vector_count)) {
        return DENSE_SOLVER_ALLOCATION_FAILED;
    }

    work_A = malloc(matrix_count * sizeof(work_A[0]));
    work_b = malloc(vector_count * sizeof(work_b[0]));
    work_x = malloc(vector_count * sizeof(work_x[0]));

    if (work_A == NULL || work_b == NULL || work_x == NULL) {
        free(work_A);
        free(work_b);
        free(work_x);
        return DENSE_SOLVER_ALLOCATION_FAILED;
    }

    for (size_t i = 0; i < matrix_count; ++i) {
        work_A[i] = A[i];
    }

    for (size_t i = 0; i < vector_count; ++i) {
        work_b[i] = b[i];
        work_x[i] = 0.0;
    }

    for (int pivot_col = 0; pivot_col < n; ++pivot_col) {
        int pivot_row = pivot_col;
        double pivot_abs = fabs(work_A[(pivot_col * n) + pivot_col]);

        for (int row = pivot_col + 1; row < n; ++row) {
            const double candidate_abs = fabs(work_A[(row * n) + pivot_col]);

            if (candidate_abs > pivot_abs) {
                pivot_abs = candidate_abs;
                pivot_row = row;
            }
        }

        if (pivot_abs == 0.0) {
            free(work_A);
            free(work_b);
            free(work_x);
            return DENSE_SOLVER_SINGULAR;
        }

        if (pivot_abs < DENSE_SOLVER_PIVOT_TOLERANCE) {
            free(work_A);
            free(work_b);
            free(work_x);
            return DENSE_SOLVER_NEAR_SINGULAR;
        }

        swap_rows(work_A, work_b, n, pivot_col, pivot_row);

        for (int row = pivot_col + 1; row < n; ++row) {
            const double factor =
                work_A[(row * n) + pivot_col] / work_A[(pivot_col * n) + pivot_col];

            work_A[(row * n) + pivot_col] = 0.0;
            for (int col = pivot_col + 1; col < n; ++col) {
                work_A[(row * n) + col] -= factor * work_A[(pivot_col * n) + col];
            }
            work_b[row] -= factor * work_b[pivot_col];
        }
    }

    for (int row = n - 1; row >= 0; --row) {
        double sum = work_b[row];
        const double diagonal = work_A[(row * n) + row];

        if (diagonal == 0.0) {
            free(work_A);
            free(work_b);
            free(work_x);
            return DENSE_SOLVER_SINGULAR;
        }

        if (fabs(diagonal) < DENSE_SOLVER_PIVOT_TOLERANCE) {
            free(work_A);
            free(work_b);
            free(work_x);
            return DENSE_SOLVER_NEAR_SINGULAR;
        }

        for (int col = row + 1; col < n; ++col) {
            sum -= work_A[(row * n) + col] * work_x[col];
        }

        work_x[row] = sum / diagonal;
    }

    for (int i = 0; i < n; ++i) {
        x[i] = work_x[i];
    }

    free(work_A);
    free(work_b);
    free(work_x);
    return DENSE_SOLVER_OK;
}
